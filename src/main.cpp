// Isosurface Singularity Lab
//
// Type a scalar field F(x,y,z,k), drag k in [-1,1], and watch the
// isosurface {F=0} and its gradient field react in real time.
//
// Designed for exploring "unfoldings" of a singular isosurface, e.g.
//   x^3 - k*y                 (cuspidal edge, singular only at k=0)
//   x^2 + y^2 - z^2 - k       (cone at k=0 -> hyperboloid for k != 0)
//   x^2*z - y^2 - k*z         (Whitney-umbrella-like family)
//
// At F=0 with grad F = 0 simultaneously, the surface has a singular
// point (the implicit function theorem fails there). The app's job is
// to make that fact visible and quantitative:
//   - the isosurface is colored by |grad F| (dark/blue near-zero = near
//     a singularity; bright/red = steep, well-behaved region)
//   - optional gradient glyphs show the exact analytic gradient
//     direction at sampled surface points (via forward-mode AD, see
//     Dual3.h -- not a finite-difference estimate)
//   - a small live plot sweeps k from -1 to 1 and shows |grad F| (and
//     each component) at a chosen "point of interest" (default: the
//     origin), with a cursor at the current k, so you can watch the
//     gradient collapse to zero and recover as k crosses the
//     bifurcation value.
#include <GLFW/glfw3.h>
#include "imgui.h"
#include "backends/imgui_impl_glfw.h"
#include "backends/imgui_impl_opengl2.h"

#include "Expr.h"
#include "Isosurface.h"
#include "Camera.h"

#include <vector>
#include <string>
#include <cstdio>
#include <cmath>
#include <algorithm>

// ---------------------------------------------------------------------
// Cool-to-warm colormap: blue/purple (|grad F| near 0, i.e. near a
// singularity) -> cyan -> green -> yellow -> red (steep / regular).
static void colormap(float t, float& r, float& g, float& b) {
    t = std::clamp(t, 0.0f, 1.0f);
    static const float stops[5][3] = {
        {0.15f, 0.05f, 0.45f}, // near-singular
        {0.10f, 0.45f, 0.75f},
        {0.25f, 0.75f, 0.35f},
        {0.95f, 0.80f, 0.15f},
        {0.90f, 0.15f, 0.15f}, // steep / regular
    };
    float scaled = t * 4.0f;
    int i0 = std::min(3, (int)scaled);
    float f = scaled - i0;
    r = stops[i0][0] + f * (stops[i0+1][0] - stops[i0][0]);
    g = stops[i0][1] + f * (stops[i0+1][1] - stops[i0][1]);
    b = stops[i0][2] + f * (stops[i0+1][2] - stops[i0][2]);
}

struct AppMesh {
    std::vector<float> pos, col, norm;
    size_t triCount = 0;
    double maxGradMag = 1e-9;
    double buildMillis = 0;
};

static void rebuildMesh(const Node* ast, double k, int res, double halfExtent,
                         float gradGamma, AppMesh& out) {
    IsosurfaceResult r = Isosurface::build(ast, k, res, halfExtent);
    out.pos.clear(); out.col.clear(); out.norm.clear();
    out.pos.reserve(r.vertices.size()*3);
    out.col.reserve(r.vertices.size()*3);
    out.norm.reserve(r.vertices.size()*3);
    out.maxGradMag = std::max(r.maxGradMag, 1e-9);
    out.buildMillis = r.buildMillis;
    out.triCount = r.vertices.size() / 3;
    for (auto& v : r.vertices) {
        out.pos.push_back(v.px); out.pos.push_back(v.py); out.pos.push_back(v.pz);
        out.norm.push_back(v.nx); out.norm.push_back(v.ny); out.norm.push_back(v.nz);
        float t = (float)(std::log(1.0 + v.gradMag) / std::log(1.0 + out.maxGradMag));
        t = std::pow(std::clamp(t,0.0f,1.0f), gradGamma);
        float cr,cg,cb; colormap(t, cr,cg,cb);
        out.col.push_back(cr); out.col.push_back(cg); out.col.push_back(cb);
    }
}

// k-sweep diagnostic data: |grad F| and components at a fixed spatial
// point, evaluated for a dense sweep of k in [-1,1].
struct KSweep {
    static constexpr int N = 201;
    float kAxis[N];
    float gradMagArr[N];
    float gx[N], gy[N], gz[N];
    float fVal[N];
};
static void rebuildKSweep(const Node* ast, double px, double py, double pz, KSweep& out) {
    for (int i = 0; i < KSweep::N; i++) {
        double k = -1.0 + 2.0 * i / (KSweep::N - 1);
        Dual3 xd(px,1,0,0), yd(py,0,1,0), zd(pz,0,0,1), kd(k,0,0,0);
        Dual3 r = evalExpr<Dual3>(ast, xd, yd, zd, kd);
        out.kAxis[i] = (float)k;
        out.gx[i] = (float)r.dx; out.gy[i] = (float)r.dy; out.gz[i] = (float)r.dz;
        out.gradMagArr[i] = (float)std::sqrt(r.dx*r.dx+r.dy*r.dy+r.dz*r.dz);
        out.fVal[i] = (float)r.v;
    }
}

// Draws a PlotLines-style widget with an extra vertical cursor line at
// the current k value, so you can see exactly where the live 3D view
// sits within the k in [-1,1] sweep.
static void scrollCallback(GLFWwindow* window, double /*xoff*/, double yoff) {
    ImGui_ImplGlfw_ScrollCallback(window, 0, yoff);
    if (ImGui::GetIO().WantCaptureMouse) return;
    auto* cam = static_cast<OrbitCamera*>(glfwGetWindowUserPointer(window));
    if (!cam) return;
    cam->distance *= std::pow(0.9f, (float)yoff);
    cam->distance = std::clamp(cam->distance, 0.1f, 50.0f);
}

static void plotWithCursor(const char* label, const float* ys, int count,
                            float curK, float yMin, float yMax, float height) {
    ImGui::PlotLines(label, ys, count, 0, nullptr, yMin, yMax, ImVec2(0, height));
    ImVec2 rmin = ImGui::GetItemRectMin();
    ImVec2 rmax = ImGui::GetItemRectMax();
    float t = (curK + 1.0f) * 0.5f; // k in [-1,1] -> [0,1]
    float cx = rmin.x + t * (rmax.x - rmin.x);
    ImDrawList* dl = ImGui::GetWindowDrawList();
    dl->AddLine(ImVec2(cx, rmin.y), ImVec2(cx, rmax.y), IM_COL32(255,255,255,180), 1.5f);
}

int main() {
    if (!glfwInit()) { std::fprintf(stderr, "glfwInit failed\n"); return 1; }
    glfwWindowHint(GLFW_SAMPLES, 4);
    GLFWwindow* window = glfwCreateWindow(1400, 900, "Isosurface Singularity Lab", nullptr, nullptr);
    if (!window) { std::fprintf(stderr, "glfwCreateWindow failed\n"); glfwTerminate(); return 1; }
    glfwMakeContextCurrent(window);
    glfwSwapInterval(1);

    OrbitCamera cam;
    glfwSetWindowUserPointer(window, &cam);
    glfwSetScrollCallback(window, scrollCallback);

    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    ImGui::StyleColorsDark();
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    ImGui_ImplOpenGL2_Init();

    // ---- app state ----
    char exprBuf[512] = "x*x*x - k*y";
    std::string lastGoodExpr;
    std::string parseError;
    NodePtr ast = parseExpression(exprBuf);
    lastGoodExpr = exprBuf;

    float k = 0.0f;
    int gridRes = 48;
    float halfExtent = 1.0f;
    float gradGamma = 0.6f; // perceptual curve for the gradient-magnitude colormap
    float poi[3] = { 0.0f, 0.0f, 0.0f }; // "point of interest" for the k-sweep plot
    bool showGradientGlyphs = true;
    bool showWireframeBox = true;
    float glyphScale = 0.15f;
    int glyphStride = 23; // draw every Nth triangle-vertex as a glyph (kept prime-ish to avoid grid aliasing patterns)

    AppMesh mesh;
    KSweep sweep;
    bool dirtyMesh = true;
    bool dirtySweep = true;

    cam.distance = 3.2f * halfExtent;
    bool draggingCam = false;
    double lastMouseX = 0, lastMouseY = 0;

    const char* presets[] = {
        "x*x*x - k*y",
        "x*x + y*y - z*z - k",
        "x*x*z - y*y - k*z",
        "x*x*x - 3*x*y*y - k",
        "x*x*x + y*y*y + z*z*z - 3*k*x*y*z",
    };
    const char* presetNames[] = {
        "Cuspidal: x^3 - k y",
        "Cone unfolding: x^2+y^2-z^2 - k",
        "Whitney-umbrella-like: x^2 z - y^2 - k z",
        "Monkey saddle: x^3 - 3xy^2 - k",
        "Cubic surface family",
    };
    int presetIdx = 0;

    while (!glfwWindowShouldClose(window)) {
        glfwPollEvents();

        // ---- mouse orbit controls (skip while ImGui wants the mouse) ----
        ImGuiIO& io2 = ImGui::GetIO();
        if (!io2.WantCaptureMouse) {
            double mx, my; glfwGetCursorPos(window, &mx, &my);
            bool leftDown = glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_PRESS;
            if (leftDown && !draggingCam) { draggingCam = true; lastMouseX = mx; lastMouseY = my; }
            if (!leftDown) draggingCam = false;
            if (draggingCam) {
                double dx = mx - lastMouseX, dy = my - lastMouseY;
                cam.yaw   += (float)dx * 0.01f;
                cam.pitch += (float)dy * 0.01f;
                cam.pitch = std::clamp(cam.pitch, -1.5f, 1.5f);
                lastMouseX = mx; lastMouseY = my;
            }
        }

        ImGui_ImplOpenGL2_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // ---- control panel ----
        ImGui::SetNextWindowPos(ImVec2(12, 12), ImGuiCond_FirstUseEver);
        ImGui::SetNextWindowSize(ImVec2(420, 820), ImGuiCond_FirstUseEver);
        ImGui::Begin("Isosurface Singularity Lab");

        ImGui::TextWrapped("F(x, y, z, k) = 0 — variables: x y z k. "
                            "Functions: sin cos tan exp log sqrt abs pow min max, ^ for power.");
        bool exprChanged = ImGui::InputText("F(x,y,z,k)", exprBuf, sizeof(exprBuf), ImGuiInputTextFlags_EnterReturnsTrue);
        ImGui::SameLine();
        bool applyClicked = ImGui::Button("Apply");

        if (ImGui::Combo("Preset", &presetIdx, presetNames, IM_ARRAYSIZE(presetNames))) {
            std::snprintf(exprBuf, sizeof(exprBuf), "%s", presets[presetIdx]);
            exprChanged = true;
        }

        if (exprChanged || applyClicked) {
            try {
                NodePtr newAst = parseExpression(exprBuf);
                ast = std::move(newAst);
                lastGoodExpr = exprBuf;
                parseError.clear();
                dirtyMesh = true;
                dirtySweep = true;
            } catch (const ExprError& e) {
                parseError = e.what();
            }
        }
        if (!parseError.empty()) {
            ImGui::TextColored(ImVec4(1,0.4f,0.4f,1), "Parse error: %s", parseError.c_str());
            ImGui::TextDisabled("(still showing the last valid surface for \"%s\")", lastGoodExpr.c_str());
        }

        ImGui::Separator();
        dirtyMesh |= ImGui::SliderFloat("k", &k, -1.0f, 1.0f, "%.3f");
        dirtyMesh |= ImGui::SliderInt("grid resolution", &gridRes, 8, 110);
        dirtyMesh |= ImGui::SliderFloat("domain half-extent", &halfExtent, 0.2f, 3.0f, "%.2f");
        dirtyMesh |= ImGui::SliderFloat("color gamma", &gradGamma, 0.2f, 2.0f, "%.2f");
        ImGui::Checkbox("show gradient glyphs", &showGradientGlyphs);
        ImGui::SameLine();
        ImGui::Checkbox("show domain box", &showWireframeBox);
        if (showGradientGlyphs) {
            ImGui::SliderFloat("glyph length", &glyphScale, 0.02f, 0.6f, "%.2f");
            ImGui::SliderInt("glyph stride", &glyphStride, 1, 80);
        }

        ImGui::Separator();
        ImGui::Text("Point of interest (for the k-sweep plot below)");
        dirtySweep |= ImGui::InputFloat3("P", poi, "%.3f");
        if (ImGui::Button("Snap P to origin")) { poi[0]=poi[1]=poi[2]=0; dirtySweep = true; }

        if (dirtyMesh) {
            try {
                rebuildMesh(ast.get(), k, gridRes, halfExtent, gradGamma, mesh);
            } catch (...) { /* evaluation never throws; parse already validated */ }
            dirtyMesh = false;
        }
        if (dirtySweep) {
            rebuildKSweep(ast.get(), poi[0], poi[1], poi[2], sweep);
            dirtySweep = false;
        }

        float maxG = 1e-6f;
        for (float v : sweep.gradMagArr) maxG = std::max(maxG, v);

        // current readout at P for the live k
        {
            Dual3 xd(poi[0],1,0,0), yd(poi[1],0,1,0), zd(poi[2],0,0,1), kd(k,0,0,0);
            Dual3 r = evalExpr<Dual3>(ast.get(), xd, yd, zd, kd);
            double gm = std::sqrt(r.dx*r.dx+r.dy*r.dy+r.dz*r.dz);
            ImGui::Text("At P, k=%.3f:  F = % .5f    |grad F| = %.5f", k, r.v, gm);
            ImGui::Text("grad F = (% .5f, % .5f, % .5f)", r.dx, r.dy, r.dz);
            if (gm < 0.02 * std::max(1.0, (double)maxG))
                ImGui::TextColored(ImVec4(1,0.85f,0.2f,1), "  -> near-singular at this (P, k): the gradient is ~vanishing.");
        }

        ImGui::Separator();
        ImGui::Text("|grad F(P,k)| across k in [-1,1]  (white line = current k)");
        plotWithCursor("##gradmag", sweep.gradMagArr, KSweep::N, k, 0.0f, maxG * 1.05f, 70.0f);

        if (ImGui::CollapsingHeader("Gradient components & F(P,k) vs k")) {
            float lo = 1e9f, hi = -1e9f;
            for (int i=0;i<KSweep::N;i++) {
                lo = std::min({lo, sweep.gx[i], sweep.gy[i], sweep.gz[i]});
                hi = std::max({hi, sweep.gx[i], sweep.gy[i], sweep.gz[i]});
            }
            ImGui::TextUnformatted("dF/dx"); plotWithCursor("##gx", sweep.gx, KSweep::N, k, lo, hi, 45.0f);
            ImGui::TextUnformatted("dF/dy"); plotWithCursor("##gy", sweep.gy, KSweep::N, k, lo, hi, 45.0f);
            ImGui::TextUnformatted("dF/dz"); plotWithCursor("##gz", sweep.gz, KSweep::N, k, lo, hi, 45.0f);
            float flo=1e9f, fhi=-1e9f;
            for (float v: sweep.fVal) { flo=std::min(flo,v); fhi=std::max(fhi,v); }
            ImGui::TextUnformatted("F(P,k)"); plotWithCursor("##fv", sweep.fVal, KSweep::N, k, flo, fhi, 45.0f);
        }

        ImGui::Separator();
        ImGui::Text("Surface: %zu triangles, built in %.1f ms", mesh.triCount, mesh.buildMillis);
        ImGui::TextDisabled("Color: blue/purple = |grad F| near 0 (singular), red = steep.");
        ImGui::TextDisabled("Drag with left mouse to orbit. Scroll to zoom.");
        if (ImGui::Button("Reset camera")) { cam.yaw=0.9f; cam.pitch=0.5f; cam.distance = 3.2f*halfExtent; }

        ImGui::End();

        // ---- 3D render ----
        int fbw, fbh; glfwGetFramebufferSize(window, &fbw, &fbh);
        glViewport(0, 0, fbw, fbh);
        glClearColor(0.07f, 0.08f, 0.10f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_MULTISAMPLE);
        glDisable(GL_CULL_FACE);

        glMatrixMode(GL_PROJECTION);
        glLoadIdentity();
        float aspect = fbh > 0 ? (float)fbw/(float)fbh : 1.0f;
        cam.applyProjection(aspect);

        glMatrixMode(GL_MODELVIEW);
        glLoadIdentity();
        cam.applyView();

        // simple directional "lighting" baked per-vertex on the CPU side is
        // already folded into vertex color via the colormap; add a mild
        // geometric shading term here using the fixed-function light so
        // curvature still reads clearly.
        glEnable(GL_LIGHTING);
        glEnable(GL_LIGHT0);
        glEnable(GL_COLOR_MATERIAL);
        glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);
        GLfloat lpos[4] = { 2.0f, 3.0f, 2.5f, 0.0f };
        glLightfv(GL_LIGHT0, GL_POSITION, lpos);
        GLfloat amb[4] = { 0.55f,0.55f,0.55f,1 };
        GLfloat diff[4] = { 0.75f,0.75f,0.75f,1 };
        glLightfv(GL_LIGHT0, GL_AMBIENT, amb);
        glLightfv(GL_LIGHT0, GL_DIFFUSE, diff);

        if (!mesh.pos.empty()) {
            glEnableClientState(GL_VERTEX_ARRAY);
            glEnableClientState(GL_COLOR_ARRAY);
            glEnableClientState(GL_NORMAL_ARRAY);
            glVertexPointer(3, GL_FLOAT, 0, mesh.pos.data());
            glColorPointer(3, GL_FLOAT, 0, mesh.col.data());
            glNormalPointer(GL_FLOAT, 0, mesh.norm.data());
            glDrawArrays(GL_TRIANGLES, 0, (GLsizei)(mesh.pos.size()/3));
            glDisableClientState(GL_VERTEX_ARRAY);
            glDisableClientState(GL_COLOR_ARRAY);
            glDisableClientState(GL_NORMAL_ARRAY);
        }
        glDisable(GL_LIGHTING);

        // gradient glyphs
        if (showGradientGlyphs && !mesh.pos.empty()) {
            glLineWidth(1.5f);
            glBegin(GL_LINES);
            size_t nv = mesh.pos.size()/3;
            for (size_t vi = 0; vi < nv; vi += (size_t)std::max(1, glyphStride)) {
                float px = mesh.pos[vi*3+0], py = mesh.pos[vi*3+1], pz = mesh.pos[vi*3+2];
                float nx = mesh.norm[vi*3+0], ny = mesh.norm[vi*3+1], nz = mesh.norm[vi*3+2];
                float r = mesh.col[vi*3+0], g = mesh.col[vi*3+1], b = mesh.col[vi*3+2];
                glColor3f(r,g,b);
                glVertex3f(px,py,pz);
                glVertex3f(px+nx*glyphScale, py+ny*glyphScale, pz+nz*glyphScale);
            }
            glEnd();
        }

        // domain wireframe box + axes + point-of-interest marker
        if (showWireframeBox) {
            float h = halfExtent;
            glColor4f(0.5f,0.5f,0.55f,0.5f);
            glLineWidth(1.0f);
            glBegin(GL_LINES);
            float c[8][3] = {{-h,-h,-h},{h,-h,-h},{h,h,-h},{-h,h,-h},{-h,-h,h},{h,-h,h},{h,h,h},{-h,h,h}};
            int edges[12][2] = {{0,1},{1,2},{2,3},{3,0},{4,5},{5,6},{6,7},{7,4},{0,4},{1,5},{2,6},{3,7}};
            for (auto& e : edges) {
                glVertex3fv(c[e[0]]); glVertex3fv(c[e[1]]);
            }
            glEnd();

            float axisLen = h * 1.15f;
            glLineWidth(2.0f);
            glBegin(GL_LINES);
            glColor3f(1,0.3f,0.3f); glVertex3f(0,0,0); glVertex3f(axisLen,0,0);
            glColor3f(0.3f,1,0.3f); glVertex3f(0,0,0); glVertex3f(0,axisLen,0);
            glColor3f(0.3f,0.5f,1); glVertex3f(0,0,0); glVertex3f(0,0,axisLen);
            glEnd();
        }
        {
            glDisable(GL_DEPTH_TEST);
            glPointSize(9.0f);
            glColor3f(1.0f, 1.0f, 1.0f);
            glBegin(GL_POINTS);
            glVertex3f(poi[0], poi[1], poi[2]);
            glEnd();
            glEnable(GL_DEPTH_TEST);
        }

        ImGui::Render();
        ImGui_ImplOpenGL2_RenderDrawData(ImGui::GetDrawData());
        glfwSwapBuffers(window);
    }

    ImGui_ImplOpenGL2_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();
    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}
