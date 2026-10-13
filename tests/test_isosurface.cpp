// Standalone correctness check for the marching-tetrahedra isosurface
#include <map>
#include <algorithm>
// extractor + AD gradients, independent of any windowing/GL code so it
// can run in a headless sandbox. Not part of the shipped app.
#include "Expr.h"
#include "Isosurface.h"
#include <cstdio>
#include <cmath>

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::printf("FAIL: %s\n", msg); failures++; } } while(0)

int main() {
    // --- Test 1: sphere x^2+y^2+z^2 - 0.25 = 0  (radius 0.5) ---
    // Every generated vertex should lie very close to radius 0.5, and its
    // analytic gradient should point radially outward (dot(normal, pos_hat) ~ 1).
    {
        auto ast = parseExpression("x*x + y*y + z*z - 0.25");
        auto res = Isosurface::build(ast.get(), 0.0, 48, 1.0);
        std::printf("[sphere] %zu verts (%zu tris), build %.2fms, maxGrad=%.4f\n",
                    res.vertices.size(), res.vertices.size()/3, res.buildMillis, res.maxGradMag);
        CHECK(res.vertices.size() > 1000, "sphere produced too few vertices");
        double maxRadiusErr = 0, minDot = 1.0;
        for (auto& v : res.vertices) {
            double r = std::sqrt(v.px*v.px + v.py*v.py + v.pz*v.pz);
            maxRadiusErr = std::max(maxRadiusErr, std::fabs(r - 0.5));
            double rl = r > 1e-9 ? 1.0/r : 0.0;
            double dot = (v.px*rl)*v.nx + (v.py*rl)*v.ny + (v.pz*rl)*v.nz;
            minDot = std::min(minDot, dot);
        }
        std::printf("[sphere] max radius error = %.5f, min(normal . radial_hat) = %.5f\n", maxRadiusErr, minDot);
        CHECK(maxRadiusErr < 0.03, "sphere vertices deviate too far from r=0.5");
        CHECK(minDot > 0.99, "sphere analytic gradient not radial / not outward-oriented");
        // expected analytic |grad| on this sphere is 2r = 1.0 everywhere
        CHECK(std::fabs(res.maxGradMag - 1.0) < 0.05, "sphere |grad F| should be ~1.0 (2r at r=0.5) everywhere");
    }

    // --- Test 2: cuspidal family x^3 - k*y. At k=0, every point of the
    // surface (the plane x=0) must have |grad F| == 0. For k != 0 it must not.
    {
        auto ast = parseExpression("x*x*x - k*y");
        auto res0 = Isosurface::build(ast.get(), 0.0, 24, 1.0);
        std::printf("[cusp k=0] %zu verts, maxGrad=%.6f\n", res0.vertices.size(), res0.maxGradMag);
        CHECK(res0.maxGradMag < 1e-6, "cusp at k=0 should have |grad F|=0 everywhere on the surface"); // note: IsosurfaceResult::maxGradMag floors at 1e-9 to avoid a div-by-zero in color normalization

        auto res1 = Isosurface::build(ast.get(), 0.4, 24, 1.0);
        std::printf("[cusp k=0.4] %zu verts, maxGrad=%.6f\n", res1.vertices.size(), res1.maxGradMag);
        CHECK(res1.maxGradMag > 0.3, "cusp at k=0.4 should have nonzero gradient (>=|k|)");
    }

    // --- Test 3: watertightness smoke test -- every edge of the triangle
    // soup should be shared by exactly 2 triangles (closed orientable
    // surface inside the domain), for a simple closed surface fully
    // contained in the box.
    {
        auto ast = parseExpression("x*x + y*y + z*z - 0.16"); // r=0.4, well inside [-1,1]^3
        auto res = Isosurface::build(ast.get(), 0.0, 32, 1.0);
        struct EdgeKey { long a, b; bool operator<(const EdgeKey& o) const { return a<o.a || (a==o.a&&b<o.b); } };
        auto quant = [](float v){ return (long)std::lround(v * 100000.0f); };
        std::map<EdgeKey,int> edgeCount;
        size_t ntri = res.vertices.size()/3;
        for (size_t t = 0; t < ntri; t++) {
            long id[3];
            for (int c = 0; c < 3; c++) {
                auto& v = res.vertices[t*3+c];
                // hash position to an integer id by coordinate quantization+lookup is
                // overkill here; instead just compare edges by quantized endpoint pairs.
                id[c] = 0; (void)v;
            }
            for (int e = 0; e < 3; e++) {
                auto& va = res.vertices[t*3+e];
                auto& vb = res.vertices[t*3+((e+1)%3)];
                long ax=quant(va.px), ay=quant(va.py), az=quant(va.pz);
                long bx=quant(vb.px), by=quant(vb.py), bz=quant(vb.pz);
                long A = ax*1000003 + ay*1009 + az;
                long B = bx*1000003 + by*1009 + bz;
                EdgeKey k{std::min(A,B), std::max(A,B)};
                edgeCount[k]++;
            }
        }
        int bad = 0;
        for (auto& [k,c] : edgeCount) if (c != 2) bad++;
        std::printf("[watertight] %zu unique edges, %d not shared by exactly 2 triangles (of %zu total)\n",
                    edgeCount.size(), bad, edgeCount.size());
        CHECK(bad == 0, "closed sphere mesh is not watertight (some edge not shared by exactly 2 triangles)");
    }

    // --- Test 4: parser error handling ---
    {
        bool threw = false;
        try { parseExpression("x + + *"); } catch (const ExprError&) { threw = true; }
        CHECK(threw, "malformed expression should throw ExprError");

        bool threw2 = false;
        try { parseExpression("x + q"); } catch (const ExprError&) { threw2 = true; }
        CHECK(threw2, "unknown identifier should throw ExprError");

        auto ast = parseExpression("sin(x) + cos(y)*k - sqrt(abs(z)) + pow(x,2) - x^2");
        Dual3 r = evalExpr<Dual3>(ast.get(), Dual3(0.3,1,0,0), Dual3(0.4,0,1,0), Dual3(-0.2,0,0,1), Dual3(0.1,0,0,0));
        std::printf("[parser] mixed-function expr evaluated: v=%.6f\n", r.v);
        CHECK(std::isfinite(r.v), "mixed-function expression should evaluate to a finite value");
    }

    if (failures == 0) { std::printf("\nALL TESTS PASSED\n"); return 0; }
    std::printf("\n%d TEST(S) FAILED\n", failures);
    return 1;
}
