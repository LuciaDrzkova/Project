// IsosurfaceViewer core: construction, the per-frame loop (idle throttling,
// surface/sweep rebuilds, overlays), mouse picking and overlay geometry.

#include "IsosurfaceViewer.h"

#include "ColorMap.h"
#include "core/Grid.h"

#include <GLFW/glfw3.h>
#include <imgui.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <limits>
#include <string>
#include <utility>

namespace iso_app {

IsosurfaceViewer::IsosurfaceViewer(const char* title, int width, int height)
    : pmp::MeshViewer(title, width, height)
{
    std::snprintf(expression_buffer_.data(), expression_buffer_.size(),
                  "%s", "x*x*x - k*y");
    std::snprintf(slice_output_directory_.data(), slice_output_directory_.size(),
                  "%s", "isosurface_slices");
    std::snprintf(batch_f_expression_buffer_.data(), batch_f_expression_buffer_.size(),
                  "%s", "x^3+y^3+z^3");
    std::snprintf(batch_g_expressions_buffer_.data(), batch_g_expressions_buffer_.size(),
                  "%s", "1\nx\ny\nz\nx+y\nx+z\ny+z\nx-y\nx-z\ny-z\n"
                         "x^2\ny^2\nz^2\nx*y\nx*z\ny*z\n"
                         "x^2+y^2+z^2\nx*y+y*z+z*x\nx^3\ny^3\nz^3\n"
                         "x*y*z\nx^2*y\nx^2*z\ny^2*z\n");
    ast_ = parseExpression(expression_buffer_.data());
    add_help_item("Left drag", "Orbit the isosurface");
    add_help_item("Click grid vertex", "Select a sample vertex and inspect its F value");
    add_help_item("Scroll", "Zoom the viewer");
    add_help_item("?", "Show key bindings");
    setup_open_ = true;
}

IsosurfaceViewer::~IsosurfaceViewer()
{
    // Background jobs touch no viewer state; their destructors wait for them.
    // The derived object is destroyed before the base Window destroys the
    // OpenGL context, so GPU resources can still be released here.
    overlay_.release();
}

// PMP's Window::run() redraws continuously at vsync rate (glfwPollEvents),
// which keeps the CPU/GPU busy even when nothing changes. Sleep here until
// input arrives, then render a short burst of frames so ImGui can settle
// (hover states, popups, layout). Keep rendering while work is pending.
void IsosurfaceViewer::throttle_when_idle()
{
    const ImGuiIO& io = ImGui::GetIO();
    bool mouse_held = false;
    for (bool down : io.MouseDown)
        mouse_held = mouse_held || down;

    // Background jobs need frames so their completion is noticed and shown.
    const bool busy = mesh_dirty_ || sweep_dirty_ || batch_job_.running() ||
                      convergence_job_.running() || mouse_held;
    if (busy)
    {
        active_frames_left_ = kActiveFramesAfterEvent;
        return;
    }
    if (active_frames_left_ > 0)
    {
        --active_frames_left_;
        return;
    }

    constexpr double timeout_seconds = 0.5;
    const double before = glfwGetTime();
    glfwWaitEventsTimeout(timeout_seconds);
    if (glfwGetTime() - before < timeout_seconds * 0.95)
        active_frames_left_ = kActiveFramesAfterEvent; // woken by an event
}

void IsosurfaceViewer::do_processing()
{
    throttle_when_idle();
    if (mesh_dirty_)
        rebuild_surface();
    if (sweep_dirty_)
        rebuild_k_sweep();
}

// Select a visible sampling-grid vertex with a click (not a drag). PMP's
// Window filters clicks over ImGui controls before dispatching this callback.
void IsosurfaceViewer::mouse(int button, int action, int mods)
{
    pmp::MeshViewer::mouse(button, action, mods);

    if (button != GLFW_MOUSE_BUTTON_LEFT)
        return;

    double x = 0.0, y = 0.0;
    cursor_pos(x, y); // PMP returns framebuffer-scaled coordinates.

    if (action == GLFW_PRESS)
    {
        selection_click_candidate_ = true;
        selection_click_start_ = {x, y};
        return;
    }

    if (action == GLFW_RELEASE && selection_click_candidate_)
    {
        selection_click_candidate_ = false;
        const double dx = x - selection_click_start_[0];
        const double dy = y - selection_click_start_[1];
        const double click_tolerance = 5.0 * std::max(1.0f, high_dpi_scaling());
        if (dx * dx + dy * dy <= click_tolerance * click_tolerance)
            select_grid_vertex_near(x, y);
    }
}

void IsosurfaceViewer::draw(const std::string& draw_mode)
{
    pmp::MeshViewer::draw(draw_mode);
    if (!overlay_.initialized())
    {
        overlay_.initialize();
        grid_geometry_dirty_ = cell_geometry_dirty_ = glyph_geometry_dirty_ = true;
        point_of_interest_geometry_dirty_ = selected_vertex_geometry_dirty_ = true;
        slice_reference_geometry_dirty_ = true;
    }
    refresh_overlays();

    using Layer = OverlayRenderer;
    using Depth = OverlayRenderer::Depth;
    const bool slices_current = slice_result_available_ && !slice_data_stale_;
    overlay_.begin(projection_matrix_, modelview_matrix_);
    if (show_grid_)
        overlay_.draw(Layer::Grid, {GL_LINES, 0.62f, 1.0f});
    // Contours lie exactly on F=0 and may share depth with the surface.
    if (show_slice_contours_ && slices_current)
        overlay_.draw(Layer::SliceContours, {GL_LINES, 1.0f, 2.0f, Depth::LessEqual});
    if (show_slice_grid_in_3d_ && slices_current)
        overlay_.draw(Layer::SlicePlane, {GL_LINES, 0.82f, 1.2f, Depth::LessEqual});
    // Keep the selected slice node visible while rotating.
    if (show_slice_selected_point_in_3d_ && slices_current)
        overlay_.draw(Layer::SliceSelectedPoint, {GL_LINES, 1.0f, 3.5f, Depth::Off});
    if (show_gradient_glyphs_)
        overlay_.draw(Layer::GradientGlyphs, {GL_LINES, 0.95f, 1.3f});
    if (show_selected_cell_)
    {
        overlay_.draw(Layer::CellEdges, {GL_LINES, 1.0f, 2.5f});
        overlay_.draw(Layer::CellCorners, {GL_POINTS, 1.0f, 7.0f});
    }
    overlay_.draw(Layer::PointOfInterest, {GL_POINTS, 1.0f, 9.0f});
    // The selected grid node is drawn on top of the surface so it stays easy
    // to locate after rotating, even when occluded. It is a real triangle
    // sphere, not a GL_POINTS sprite (some drivers only support 1px points).
    if (has_selected_grid_vertex_)
        overlay_.draw(Layer::SelectedGridVertex, {GL_TRIANGLES, 1.0f, 1.0f, Depth::Off});
    overlay_.end();
}

void IsosurfaceViewer::refresh_overlays()
{
    if (grid_geometry_dirty_)
    {
        overlay_.set(OverlayRenderer::Grid, build_grid_geometry());
        grid_geometry_dirty_ = false;
    }
    if (cell_geometry_dirty_)
    {
        std::vector<OverlayVertex> edges, corners;
        build_cell_geometry(edges, corners);
        overlay_.set(OverlayRenderer::CellEdges, std::move(edges));
        overlay_.set(OverlayRenderer::CellCorners, std::move(corners));
        cell_geometry_dirty_ = false;
    }
    if (glyph_geometry_dirty_)
    {
        overlay_.set(OverlayRenderer::GradientGlyphs, build_glyph_geometry());
        glyph_geometry_dirty_ = false;
    }
    if (point_of_interest_geometry_dirty_)
    {
        overlay_.set(OverlayRenderer::PointOfInterest, build_point_of_interest_geometry());
        point_of_interest_geometry_dirty_ = false;
    }
    if (selected_vertex_geometry_dirty_)
    {
        overlay_.set(OverlayRenderer::SelectedGridVertex, build_selected_grid_vertex_geometry());
        selected_vertex_geometry_dirty_ = false;
    }
    if (slice_reference_geometry_dirty_)
    {
        std::vector<OverlayVertex> plane, point;
        build_slice_reference_geometry(plane, point);
        overlay_.set(OverlayRenderer::SlicePlane, std::move(plane));
        overlay_.set(OverlayRenderer::SliceSelectedPoint, std::move(point));
        slice_reference_geometry_dirty_ = false;
    }
}

void IsosurfaceViewer::sanitize_parameters()
{
    if (!std::isfinite(k_)) k_ = 0.0f;
    k_ = std::clamp(k_, -1.0f, 1.0f);
    grid_resolution_ = std::clamp(grid_resolution_, 1, kMaxGridResolution);
    if (!std::isfinite(half_extent_)) half_extent_ = 1.0f;
    half_extent_ = std::clamp(half_extent_, 0.2f, 3.0f);
    if (!std::isfinite(gradient_gamma_)) gradient_gamma_ = 0.6f;
    gradient_gamma_ = std::clamp(gradient_gamma_, 0.2f, 2.0f);
    grid_stride_ = std::clamp(grid_stride_, 1, grid_resolution_);
    for (int axis = 0; axis < 3; ++axis)
    {
        selected_cell_[axis] = std::clamp(selected_cell_[axis], 0, grid_resolution_ - 1);
        const int clamped_vertex = std::clamp(selected_grid_vertex_[axis], 0, grid_resolution_);
        if (clamped_vertex != selected_grid_vertex_[axis])
        {
            selected_grid_vertex_[axis] = clamped_vertex;
            selected_vertex_geometry_dirty_ = true;
        }
    }
}

void IsosurfaceViewer::apply_expression()
{
    try
    {
        NodePtr next = parseExpression(expression_buffer_.data());
        ast_ = std::move(next);
        last_good_expression_ = expression_buffer_.data();
        parse_error_.clear();
        runtime_message_.clear();
        invalidate_slice_overlay();
        mesh_dirty_ = true;
        sweep_dirty_ = true;
        glyph_geometry_dirty_ = true;
    }
    catch (const ExprError& error)
    {
        parse_error_ = error.what();
        runtime_message_ = "Invalid expression; keeping the last valid surface.";
    }
}

void IsosurfaceViewer::rebuild_surface()
{
    mesh_dirty_ = false;
    runtime_message_.clear();
    const auto start = std::chrono::steady_clock::now();
    try
    {
        IsosurfaceResult result = Isosurface::build(
            ast_.get(), static_cast<double>(k_), grid_resolution_, half_extent_);
        if (result.vertices.empty())
        {
            runtime_message_ = "No zero-level surface was found in this domain; showing the last non-empty surface.";
            last_build_ms_ = result.buildMillis;
            return;
        }

        // The extractor emits triangle soup (three corners per triangle)
        // plus, per corner, the id of the grid-edge vertex it came from.
        // Reconstruct edge-connected vertex fans before handing it to PMP:
        // smooth normals/colors are shared on regular parts, while
        // disconnected fans at singularities remain separate.
        mesh_.clear();
        auto vertex_colors = mesh_.vertex_property<pmp::Color>(
            "v:color", pmp::Color(0.78f, 0.80f, 0.84f));
        surface_color_denominator_ = std::log(1.0 + std::max(result.maxGradMag, 1e-9));
        surface_vertex_grad_mag_.clear();

        const std::size_t triangle_count = result.vertices.size() / 3;
        const std::size_t corner_count = triangle_count * 3;
        const bool indexed = result.vertexIds.size() == corner_count;
        auto geometric_id = [&](std::size_t corner) -> std::uint32_t {
            return indexed ? result.vertexIds[corner] : static_cast<std::uint32_t>(corner);
        };

        // Every undirected edge with the face using it; sorting groups the
        // uses of each edge together (replaces a map of per-edge vectors).
        struct EdgeUse
        {
            std::uint64_t key;
            std::uint32_t face;
            bool operator<(const EdgeUse& other) const
            {
                return key != other.key ? key < other.key : face < other.face;
            }
        };
        std::vector<EdgeUse> edge_uses;
        edge_uses.reserve(corner_count);
        for (std::size_t face = 0; face < triangle_count; ++face)
        {
            for (std::size_t local = 0; local < 3; ++local)
            {
                std::uint64_t a = geometric_id(face * 3 + local);
                std::uint64_t b = geometric_id(face * 3 + (local + 1) % 3);
                if (a > b) std::swap(a, b);
                edge_uses.push_back({(a << 32) | b, static_cast<std::uint32_t>(face)});
            }
        }
        std::sort(edge_uses.begin(), edge_uses.end());

        std::vector<std::size_t> parent(corner_count);
        std::vector<unsigned int> tree_rank(corner_count, 0);
        for (std::size_t i = 0; i < corner_count; ++i) parent[i] = i;
        auto find_root = [&](std::size_t x)
        {
            while (parent[x] != x)
            {
                parent[x] = parent[parent[x]];
                x = parent[x];
            }
            return x;
        };
        auto unite = [&](std::size_t a, std::size_t b)
        {
            a = find_root(a);
            b = find_root(b);
            if (a == b) return;
            if (tree_rank[a] < tree_rank[b]) std::swap(a, b);
            parent[b] = a;
            if (tree_rank[a] == tree_rank[b]) ++tree_rank[a];
        };
        auto corner_for = [&](std::size_t face, std::uint64_t id)
        {
            const std::size_t first = face * 3;
            for (std::size_t local = 0; local < 3; ++local)
                if (geometric_id(first + local) == id)
                    return first + local;
            return first; // an edge endpoint always belongs to its face
        };
        for (std::size_t begin = 0; begin < edge_uses.size();)
        {
            std::size_t end = begin + 1;
            while (end < edge_uses.size() && edge_uses[end].key == edge_uses[begin].key) ++end;
            // Only two-face edges define a shared manifold fan. Boundary
            // edges remain open; edges with >2 faces are deliberately split.
            if (end - begin == 2)
            {
                const std::uint64_t key = edge_uses[begin].key;
                const std::size_t f0 = edge_uses[begin].face, f1 = edge_uses[begin + 1].face;
                for (const std::uint64_t id : {key >> 32, key & 0xffffffffu})
                    unite(corner_for(f0, id), corner_for(f1, id));
            }
            begin = end;
        }

        std::vector<pmp::Vertex> fan_vertex(corner_count); // indexed by union-find root
        std::array<pmp::Vertex, 3> triangle;
        for (std::size_t face = 0; face < triangle_count; ++face)
        {
            for (std::size_t local = 0; local < 3; ++local)
            {
                const std::size_t corner = face * 3 + local;
                pmp::Vertex& pv = fan_vertex[find_root(corner)];
                if (!pv.is_valid())
                {
                    const SurfaceVertex& source = result.vertices[corner];
                    pv = mesh_.add_vertex(pmp::Point(source.px, source.py, source.pz));
                    vertex_colors[pv] = surface_color(source.gradMag);
                    surface_vertex_grad_mag_.push_back(source.gradMag);
                }
                triangle[local] = pv;
            }
            mesh_.add_triangle(triangle[0], triangle[1], triangle[2]);
        }
        renderer_.set_use_colors(true);
        update_mesh();
        last_surface_result_ = std::move(result);
        last_triangle_count_ = last_surface_result_.vertices.size() / 3;
        last_build_ms_ = std::chrono::duration<double, std::milli>(
            std::chrono::steady_clock::now() - start).count();
        glyph_geometry_dirty_ = true;
        cell_geometry_dirty_ = true;
        if (!first_mesh_built_)
        {
            set_scene(pmp::vec3(0.0f, 0.0f, 0.0f), std::max(0.5f, 1.8f * half_extent_));
            first_mesh_built_ = true;
        }
    }
    catch (const std::exception& error)
    {
        runtime_message_ = std::string("Surface rebuild failed: ") + error.what();
    }
}

// Recolor the displayed surface from its stored |grad F| (gamma changed).
void IsosurfaceViewer::recolor_surface()
{
    if (surface_vertex_grad_mag_.size() != mesh_.n_vertices())
        return;
    auto colors = mesh_.vertex_property<pmp::Color>("v:color");
    for (const pmp::Vertex v : mesh_.vertices())
        colors[v] = surface_color(surface_vertex_grad_mag_[v.idx()]);
    update_mesh();
}

// Surface/glyph color for a gradient magnitude: log-scaled against the
// surface maximum, then gamma-corrected (blue/purple = flat, red = steep).
pmp::Color IsosurfaceViewer::surface_color(float grad_mag) const
{
    float t = static_cast<float>(std::log(1.0 + std::max(0.0f, grad_mag)) / surface_color_denominator_);
    t = std::pow(std::clamp(t, 0.0f, 1.0f), gradient_gamma_);
    return gradient_color(t);
}

void IsosurfaceViewer::rebuild_k_sweep()
{
    sweep_dirty_ = false;
    for (int i = 0; i < KSweep::N; ++i)
    {
        const double sweep_k = -1.0 + 2.0 * static_cast<double>(i) / (KSweep::N - 1);
        Dual3 x(point_of_interest_[0], 1, 0, 0);
        Dual3 y(point_of_interest_[1], 0, 1, 0);
        Dual3 z(point_of_interest_[2], 0, 0, 1);
        Dual3 k(sweep_k, 0, 0, 0);
        const Dual3 value = evalExpr<Dual3>(ast_.get(), x, y, z, k);
        sweep_.gx[i] = static_cast<float>(value.dx);
        sweep_.gy[i] = static_cast<float>(value.dy);
        sweep_.gz[i] = static_cast<float>(value.dz);
        sweep_.grad_magnitude[i] = static_cast<float>(
            std::sqrt(value.dx * value.dx + value.dy * value.dy + value.dz * value.dz));
        sweep_.f_value[i] = static_cast<float>(value.v);
    }
}

void IsosurfaceViewer::get_cell_corners(GridCorner corners[8]) const
{
    static constexpr const char* names[8] = {"C000", "C100", "C110", "C010", "C001", "C101", "C111", "C011"};
    static constexpr int offsets[8][3] = {
        {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
        {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
    for (int i = 0; i < 8; ++i)
    {
        const int ix = selected_cell_[0] + offsets[i][0];
        const int iy = selected_cell_[1] + offsets[i][1];
        const int iz = selected_cell_[2] + offsets[i][2];
        corners[i].name = names[i];
        corners[i].ix = ix;
        corners[i].iy = iy;
        corners[i].iz = iz;
        corners[i].x = gridCoordinate(ix, grid_resolution_, half_extent_);
        corners[i].y = gridCoordinate(iy, grid_resolution_, half_extent_);
        corners[i].z = gridCoordinate(iz, grid_resolution_, half_extent_);
        corners[i].value = evalExpr<double>(ast_.get(), corners[i].x, corners[i].y, corners[i].z, k_);
    }
}

std::vector<int> IsosurfaceViewer::visible_grid_samples() const
{
    const int stride = std::max(1, grid_stride_);
    std::vector<int> samples;
    for (int i = 0; i <= grid_resolution_; i += stride)
        samples.push_back(i);
    if (samples.empty() || samples.back() != grid_resolution_)
        samples.push_back(grid_resolution_);
    return samples;
}

void IsosurfaceViewer::select_grid_vertex_near(double mouse_x, double mouse_y)
{
    if (!show_grid_)
        return;

    const std::vector<int> samples = visible_grid_samples();
    if (samples.empty() || width() <= 0 || height() <= 0)
        return;

    const pmp::mat4 model_projection = projection_matrix_ * modelview_matrix_;
    const double pick_radius = 18.0 * std::max(1.0f, high_dpi_scaling());
    double best_distance_squared = pick_radius * pick_radius;
    double best_depth = std::numeric_limits<double>::max();
    std::array<int, 3> best_index{-1, -1, -1};

    // Only sample nodes at line intersections currently shown by Grid line stride.
    // If two nodes overlap in screen space, prefer the one nearer the camera.
    for (int ix : samples)
    {
        const double x = gridCoordinate(ix, grid_resolution_, half_extent_);
        for (int iy : samples)
        {
            const double y = gridCoordinate(iy, grid_resolution_, half_extent_);
            for (int iz : samples)
            {
                const double z = gridCoordinate(iz, grid_resolution_, half_extent_);
                const pmp::vec4 clip = model_projection *
                    pmp::vec4(static_cast<float>(x), static_cast<float>(y),
                              static_cast<float>(z), 1.0f);
                const double w = static_cast<double>(clip[3]);
                if (!std::isfinite(w) || w <= 1e-6)
                    continue;

                const double ndc_x = static_cast<double>(clip[0]) / w;
                const double ndc_y = static_cast<double>(clip[1]) / w;
                const double ndc_z = static_cast<double>(clip[2]) / w;
                if (!std::isfinite(ndc_x) || !std::isfinite(ndc_y) ||
                    !std::isfinite(ndc_z) || ndc_x < -1.0 || ndc_x > 1.0 ||
                    ndc_y < -1.0 || ndc_y > 1.0 || ndc_z < -1.0 || ndc_z > 1.0)
                    continue;

                const double screen_x = (ndc_x * 0.5 + 0.5) * width();
                const double screen_y = (0.5 - ndc_y * 0.5) * height();
                const double dx = screen_x - mouse_x;
                const double dy = screen_y - mouse_y;
                const double distance_squared = dx * dx + dy * dy;
                if (distance_squared < best_distance_squared - 0.25 ||
                    (std::abs(distance_squared - best_distance_squared) <= 0.25 &&
                     ndc_z < best_depth))
                {
                    best_distance_squared = distance_squared;
                    best_depth = ndc_z;
                    best_index = {ix, iy, iz};
                }
            }
        }
    }

    if (best_index[0] < 0)
        return;

    selected_grid_vertex_ = best_index;
    has_selected_grid_vertex_ = true;
    for (int axis = 0; axis < 3; ++axis)
        selected_cell_[axis] = std::min(selected_grid_vertex_[axis], grid_resolution_ - 1);

    cell_geometry_dirty_ = true;
    selected_vertex_geometry_dirty_ = true;
    runtime_message_ = "Selected grid vertex [" + std::to_string(best_index[0]) + ", " +
        std::to_string(best_index[1]) + ", " + std::to_string(best_index[2]) + "].";
}

std::vector<OverlayVertex> IsosurfaceViewer::build_grid_geometry() const
{
    std::vector<OverlayVertex> vertices;
    std::vector<double> c; // coordinates of the visible grid lines
    for (int index : visible_grid_samples())
        c.push_back(gridCoordinate(index, grid_resolution_, half_extent_));
    const double h = half_extent_;
    const float gray_r = 0.24f, gray_g = 0.27f, gray_b = 0.32f;
    auto add_line = [&](double ax, double ay, double az, double bx, double by, double bz)
    {
        vertices.push_back(overlay_vertex(ax, ay, az, gray_r, gray_g, gray_b));
        vertices.push_back(overlay_vertex(bx, by, bz, gray_r, gray_g, gray_b));
    };
    vertices.reserve(6 * c.size() * c.size());
    for (double a : c)
        for (double b : c)
        {
            add_line(a, b, -h, a, b, h); // z-parallel
            add_line(a, -h, b, a, h, b); // y-parallel
            add_line(-h, a, b, h, a, b); // x-parallel
        }
    return vertices;
}

void IsosurfaceViewer::build_cell_geometry(std::vector<OverlayVertex>& edge_lines,
                                           std::vector<OverlayVertex>& corner_points) const
{
    edge_lines.clear();
    corner_points.clear();
    GridCorner corners[8];
    get_cell_corners(corners);
    static constexpr int edges[12][2] = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
        {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
    for (const auto& edge : edges)
    {
        const auto& a = corners[edge[0]];
        const auto& b = corners[edge[1]];
        edge_lines.push_back(overlay_vertex(a.x, a.y, a.z, 0.05f, 0.90f, 1.0f));
        edge_lines.push_back(overlay_vertex(b.x, b.y, b.z, 0.05f, 0.90f, 1.0f));
    }
    for (const auto& corner : corners)
        corner_points.push_back(overlay_vertex(corner.x, corner.y, corner.z, 1.0f, 1.0f, 1.0f));
}

std::vector<OverlayVertex> IsosurfaceViewer::build_glyph_geometry() const
{
    std::vector<OverlayVertex> vertices;
    if (!show_gradient_glyphs_ || last_surface_result_.vertices.empty())
        return vertices;
    const std::size_t stride = static_cast<std::size_t>(std::max(1, glyph_stride_));
    for (std::size_t i = 0; i < last_surface_result_.vertices.size(); i += stride)
    {
        const SurfaceVertex& v = last_surface_result_.vertices[i];
        const pmp::Color color = surface_color(v.gradMag);
        vertices.push_back(overlay_vertex(v.px, v.py, v.pz, color[0], color[1], color[2]));
        vertices.push_back(overlay_vertex(
            v.px + v.nx * glyph_scale_, v.py + v.ny * glyph_scale_, v.pz + v.nz * glyph_scale_,
            color[0], color[1], color[2]));
    }
    return vertices;
}

std::vector<OverlayVertex> IsosurfaceViewer::build_point_of_interest_geometry() const
{
    std::vector<OverlayVertex> vertices;
    vertices.push_back(overlay_vertex(
        point_of_interest_[0], point_of_interest_[1], point_of_interest_[2],
        1.0f, 0.72f, 0.15f));
    return vertices;
}

std::vector<OverlayVertex> IsosurfaceViewer::build_selected_grid_vertex_geometry() const
{
    std::vector<OverlayVertex> vertices;
    if (!has_selected_grid_vertex_)
        return vertices;

    // Center the marker exactly on the selected sample. Use a triangle sphere
    // instead of GL_POINTS: some macOS/Metal OpenGL implementations clamp
    // point sprites to one pixel, making a point-size-based marker invisible.
    const double cx = gridCoordinate(selected_grid_vertex_[0], grid_resolution_, half_extent_);
    const double cy = gridCoordinate(selected_grid_vertex_[1], grid_resolution_, half_extent_);
    const double cz = gridCoordinate(selected_grid_vertex_[2], grid_resolution_, half_extent_);
    const double spacing = (2.0 * static_cast<double>(half_extent_)) /
                           static_cast<double>(std::max(1, grid_resolution_));
    //const double radius = std::max(spacing * 0.18, static_cast<double>(half_extent_) * 0.012);
    const double radius = spacing * 0.06;
    constexpr int rings = 10;
    constexpr int sectors = 16;
    constexpr float red = 1.0f, green = 0.12f, blue = 0.025f;

    auto point = [&](int ring, int sector)
    {
        const double theta = 3.14159265358979323846 * static_cast<double>(ring) / rings;
        const double phi = 2.0 * 3.14159265358979323846 * static_cast<double>(sector) / sectors;
        const double st = std::sin(theta);
        return std::array<double, 3>{
            cx + radius * st * std::cos(phi),
            cy + radius * std::cos(theta),
            cz + radius * st * std::sin(phi)};
    };
    auto add_triangle = [&](const std::array<double, 3>& a,
                            const std::array<double, 3>& b,
                            const std::array<double, 3>& c)
    {
        vertices.push_back(overlay_vertex(a[0], a[1], a[2], red, green, blue));
        vertices.push_back(overlay_vertex(b[0], b[1], b[2], red, green, blue));
        vertices.push_back(overlay_vertex(c[0], c[1], c[2], red, green, blue));
    };

    for (int r = 0; r < rings; ++r)
    {
        for (int sec = 0; sec < sectors; ++sec)
        {
            const auto p00 = point(r, sec);
            const auto p01 = point(r, sec + 1);
            const auto p10 = point(r + 1, sec);
            const auto p11 = point(r + 1, sec + 1);
            if (r != 0)
                add_triangle(p00, p10, p01);
            if (r != rings - 1)
                add_triangle(p01, p10, p11);
        }
    }
    return vertices;
}

} // namespace iso_app
