// Standalone correctness check for the (cube-based) Marching Cubes
// isosurface extractor + AD gradients, independent of any windowing/GL
// code so it can run in a headless sandbox. Not part of the shipped app.
#include "Expr.h"
#include "Isosurface.h"
#include <cstdio>
#include <cmath>
#include <map>
#include <array>
#include <algorithm>

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::printf("FAIL: %s\n", msg); failures++; } } while(0)

// ---------------------------------------------------------------------
// Shared mesh-topology helper, used by several tests below. It welds the
// (non-indexed) triangle soup's vertices by quantized position, then
// reports vertex/edge/triangle counts and, for a fully closed surface,
// the Euler characteristic V - E + F (2 for a sphere, 0 for a torus) --
// a strong, independent check that the 256-case triangle table is wired
// up correctly, since a wrong table entry anywhere typically shows up as
// either a hole (chi too high) or a non-manifold/extra edge (chi wrong
// and nonManifoldEdges > 0).
// ---------------------------------------------------------------------
struct MeshTopology {
    size_t numTriangles = 0;
    size_t numUniqueVertices = 0;
    size_t numUniqueEdges = 0;
    int boundaryEdges = 0;  // shared by exactly 1 triangle -- expected where the surface exits the domain box
    int badEdges = 0;       // shared by neither 1 nor 2 triangles -- a real crack or degenerate sliver
    long eulerCharacteristic() const {
        return (long)numUniqueVertices - (long)numUniqueEdges + (long)numTriangles;
    }
};

static MeshTopology analyzeTopology(const std::vector<SurfaceVertex>& verts) {
    auto quant = [](float v) { return (long long)std::llround((double)v * 100000.0); };
    std::map<std::array<long long,3>, int> vertexId;
    std::map<std::pair<int,int>, int> edgeCount;
    size_t ntri = verts.size() / 3;
    auto getId = [&](const SurfaceVertex& v) -> int {
        std::array<long long,3> key{ quant(v.px), quant(v.py), quant(v.pz) };
        auto it = vertexId.find(key);
        if (it != vertexId.end()) return it->second;
        int id = (int)vertexId.size();
        vertexId[key] = id;
        return id;
    };
    for (size_t t = 0; t < ntri; t++) {
        int ids[3];
        for (int c = 0; c < 3; c++) ids[c] = getId(verts[t*3 + c]);
        for (int e = 0; e < 3; e++) {
            int a = ids[e], b = ids[(e + 1) % 3];
            edgeCount[{std::min(a,b), std::max(a,b)}]++;
        }
    }
    MeshTopology topo;
    topo.numTriangles = ntri;
    topo.numUniqueVertices = vertexId.size();
    topo.numUniqueEdges = edgeCount.size();
    for (auto& kv : edgeCount) {
        if (kv.second == 1) topo.boundaryEdges++;
        else if (kv.second != 2) topo.badEdges++;
    }
    return topo;
}

int main() {
    // --- Test 1: sphere x^2+y^2+z^2 - 0.25 = 0  (radius 0.5) ---
    // Every generated vertex should lie very close to radius 0.5, and its
    // analytic gradient should point radially outward (dot(normal, pos_hat) ~ 1).
    {
        auto ast = parseExpression("x*x + y*y + z*z - 0.2209");
        auto res = Isosurface::build(ast.get(), 0.0, 48, 1.0);
        std::printf("[sphere] %zu verts (%zu tris), build %.2fms, maxGrad=%.4f\n",
                    res.vertices.size(), res.vertices.size()/3, res.buildMillis, res.maxGradMag);
        CHECK(res.vertices.size() > 1000, "sphere produced too few vertices");
        double maxRadiusErr = 0, minDot = 1.0;
        for (auto& v : res.vertices) {
            double r = std::sqrt(v.px*v.px + v.py*v.py + v.pz*v.pz);
            maxRadiusErr = std::max(maxRadiusErr, std::fabs(r - 0.47));
            double rl = r > 1e-9 ? 1.0/r : 0.0;
            double dot = (v.px*rl)*v.nx + (v.py*rl)*v.ny + (v.pz*rl)*v.nz;
            minDot = std::min(minDot, dot);
        }
        std::printf("[sphere] max radius error = %.5f, min(normal . radial_hat) = %.5f\n", maxRadiusErr, minDot);
        CHECK(maxRadiusErr < 0.01, "sphere vertices deviate too far from r=0.47");
        CHECK(minDot > 0.99, "sphere analytic gradient not radial / not outward-oriented");
        // For F = x*x + y*y + z*z - r*r, |grad F| = 2r.
        CHECK(std::fabs(res.maxGradMag - 0.94) < 0.05,"sphere |grad F| should be approximately 0.94 (2r at r=0.47)");

        // Topology: a sphere is a closed genus-0 surface, Euler characteristic 2.
        auto topo = analyzeTopology(res.vertices);
        std::printf("[sphere] V=%zu E=%zu F=%zu boundaryEdges=%d badEdges=%d chi=%ld (expect 2)\n",
                    topo.numUniqueVertices, topo.numUniqueEdges, topo.numTriangles,
                    topo.boundaryEdges, topo.badEdges, topo.eulerCharacteristic());
        CHECK(res.skippedDegenerateTriangles == 0,
              "grid-aligned sphere emitted skipped degenerate triangles");
        CHECK(topo.boundaryEdges == 0 && topo.badEdges == 0, "sphere mesh is not a closed manifold (marching cubes crack/degenerate triangle/hole)");
        CHECK(topo.eulerCharacteristic() == 2, "sphere mesh Euler characteristic should be 2 (genus 0)");
        std::printf("[sphere] skipped degenerate triangles = %zu\n", res.skippedDegenerateTriangles);

    }

    // --- Exact-zero grid-alignment regression (sphere radius 0.5) ---
    // At resolution 48 over [-1,1], the grid contains exact samples on this
    // sphere. This explicitly exercises the exact-zero tie-break, rather than
    // relying only on a sphere whose radius misses the grid points.
    {
        auto exactZeroAst = parseExpression("x*x + y*y + z*z - 0.25");
        auto exactZeroRes = Isosurface::build(exactZeroAst.get(), 0.0, 48, 1.0);
        const auto exactZeroTopo = analyzeTopology(exactZeroRes.vertices);
        std::printf(
            "[exact-zero sphere] exactZeroSamples=%zu skippedDegenerate=%zu "
            "V=%zu E=%zu F=%zu boundaryEdges=%d badEdges=%d chi=%ld\n",
            exactZeroRes.exactZeroGridSamples,
            exactZeroRes.skippedDegenerateTriangles,
            exactZeroTopo.numUniqueVertices, exactZeroTopo.numUniqueEdges,
            exactZeroTopo.numTriangles, exactZeroTopo.boundaryEdges,
            exactZeroTopo.badEdges, exactZeroTopo.eulerCharacteristic());
        CHECK(exactZeroRes.exactZeroGridSamples > 0,
              "exact-zero sphere regression did not hit any exact-zero grid samples");
        CHECK(exactZeroRes.skippedDegenerateTriangles == 0,
              "exact-zero sphere still emitted skipped degenerate triangles");
        CHECK(exactZeroTopo.boundaryEdges == 0 && exactZeroTopo.badEdges == 0,
              "exact-zero sphere mesh is not a closed manifold");
        CHECK(exactZeroTopo.eulerCharacteristic() == 2,
              "exact-zero sphere Euler characteristic should be 2 (genus 0)");
    }

    // --- Test 2: torus (R=0.35 tube radius r=0.15), a genus-1 surface. ---
    // This exercises a much richer set of the 256 marching-cubes cases
    // than a sphere does (saddle regions on the inner/outer equator), so
    // it's a good stress test of the triangle table specifically -- a
    // wrong or transposed table entry anywhere tends to show up here as
    // either a hole (nonzero boundary edges / chi too high) or a crack
    // (badEdges > 0), even if a sphere alone still looked fine.
    {
        auto ast = parseExpression("pow(sqrt(x*x + y*y) - 0.35, 2) + z*z - 0.0225");
        for (int res : {32, 48, 64}) {
            auto result = Isosurface::build(ast.get(), 0.0, res, 1.0);
            auto topo = analyzeTopology(result.vertices);
            std::printf("[torus res=%d] V=%zu E=%zu F=%zu boundaryEdges=%d badEdges=%d chi=%ld (expect 0, genus 1)\n",
                        res, topo.numUniqueVertices, topo.numUniqueEdges, topo.numTriangles,
                        topo.boundaryEdges, topo.badEdges, topo.eulerCharacteristic());
            CHECK(topo.boundaryEdges == 0 && topo.badEdges == 0, "torus mesh is not a closed manifold (marching cubes crack/degenerate triangle/hole)");
            CHECK(topo.eulerCharacteristic() == 0, "torus mesh Euler characteristic should be 0 (genus 1)");
        }
    }

    // --- Test 3: cuspidal family x^3 - k*y -- the whole point of this app.
    // At k=0 the family is singular (F = x^3, zero set is the plane x=0,
    // and grad F = (3x^2,0,0) = 0 identically on that plane: the surface
    // "folds" with no well-defined normal anywhere). For k != 0 the k*y
    // term unfolds it: the zero set becomes x = cbrt(k*y), a smooth graph
    // over (y,z), and grad F = (3x^2, -k, 0) is nonzero everywhere on it
    // (in fact |grad F| >= |k| everywhere on the surface, with equality
    // exactly at the single point x=y=0, since |grad F|^2 = 9x^4 + k^2).
    {
        auto ast = parseExpression("x*x*x - k*y");

        auto res0 = Isosurface::build(ast.get(), 0.0, 24, 1.0);
        std::printf("[cusp k=0] %zu verts, maxGrad=%.6f\n", res0.vertices.size(), res0.maxGradMag);
        CHECK(res0.maxGradMag < 1e-6, "cusp at k=0 should have |grad F|=0 everywhere on the surface"); // note: IsosurfaceResult::maxGradMag floors at 1e-9 to avoid a div-by-zero in color normalization

        auto res1 = Isosurface::build(ast.get(), 0.4, 24, 1.0);
        std::printf("[cusp k=0.4] %zu verts, maxGrad=%.6f\n", res1.vertices.size(), res1.maxGradMag);
        CHECK(res1.maxGradMag > 0.3, "cusp at k=0.4 should have nonzero gradient (>=|k|)");

        // --- 3a: min|grad F| on the surface tracks |k| almost exactly.
        // This is the quantitative version of "k=0.01 has (almost) no
        // singularity, it's just very close to one": the worst point on
        // the whole surface still has gradient magnitude ~|k|, not ~0.
        for (double k : { 0.1, 0.01, -0.1, -0.01 }) {
            auto res = Isosurface::build(ast.get(), k, 48, 1.0);
            double minGrad = 1e18;
            for (auto& v : res.vertices) minGrad = std::min(minGrad, (double)v.gradMag);
            std::printf("[cusp k=%+6.3f] min|grad F| on surface = %.6f (|k|=%.6f)\n", k, minGrad, std::fabs(k));
            // allow some slack for grid interpolation error, but it should
            // land very close to |k| at this resolution (see dev notes:
            // observed ratio is ~1.000 for k=0.1/0.01 at res=48).
            CHECK(minGrad > std::fabs(k) * 0.9 && minGrad < std::fabs(k) * 1.2,
                  "min|grad F| on the unfolded cusp surface should track |k| closely");
        }

        // --- 3b: how far the surface strays from the singular plane x=0
        // shrinks like |k|^(1/3), not linearly in k. This is the
        // characteristic "unfolding" scaling for x^3=k*y: solving for x
        // gives x=cbrt(k*y), so over a fixed y-range the maximum |x|
        // reached scales as |k|^(1/3). It's *this* sublinear scaling --
        // not a linear one -- that makes small-|k| surfaces look like
        // they're still "almost through the singularity" even though
        // they're technically smooth: halving k only shrinks the surface's
        // excursion from x=0 by a factor of 2^(1/3)=~0.79, not by half.
        {
            const double h = 1.0;
            double prevMaxAbsX = -1.0, prevK = 0.0;
            for (double k : { 0.1, 0.0125 }) { // ratio 8 -> cube-root ratio should be ~2
                auto res = Isosurface::build(ast.get(), k, 64, h);
                double maxAbsX = 0;
                for (auto& v : res.vertices) maxAbsX = std::max(maxAbsX, (double)std::fabs(v.px));
                double predicted = std::cbrt(std::fabs(k) * h);
                std::printf("[cusp k=%.4f] max|x| on surface = %.6f (predicted cbrt(|k|*h) = %.6f)\n",
                            k, maxAbsX, predicted);
                CHECK(std::fabs(maxAbsX - predicted) < 0.02,
                      "max|x| excursion of the unfolded cusp surface should follow the cbrt(|k|) scaling law");
                if (prevMaxAbsX > 0) {
                    double ratio = prevMaxAbsX / maxAbsX;
                    std::printf("   excursion ratio for 8x change in |k| = %.4f (expect cbrt(8)=2.0)\n", ratio);
                    CHECK(ratio > 1.8 && ratio < 2.2, "excursion ratio across an 8x change in |k| should be close to cbrt(8)=2");
                }
                prevMaxAbsX = maxAbsX; prevK = k;
            }
            (void)prevK;
        }

        // --- 3c: mesh integrity holds up even very close to the
        // singularity (k small but nonzero). This is exactly the regime
        // where naive marching cubes is most at risk of cracking, because
        // the surface becomes nearly tangent to many grid cell faces near
        // x=y=0 -- so this specifically targets the degenerate/near-degenerate
        // triangle case. An open surface clipped by the domain box will
        // have some boundary edges (shared by 1 triangle); the check is
        // that no edge is EVER shared by anything other than 1 or 2.
        for (double k : { 0.0, 0.1, 0.01, 0.001, -0.01 }) {
            auto res = Isosurface::build(ast.get(), k, 48, 1.0);
            auto topo = analyzeTopology(res.vertices);
            std::printf("[cusp k=%+6.3f] V=%zu E=%zu F=%zu boundaryEdges=%d badEdges=%d (boundary edges expected: surface exits the domain box)\n",
                        k, topo.numUniqueVertices, topo.numUniqueEdges, topo.numTriangles, topo.boundaryEdges, topo.badEdges);
            CHECK(topo.badEdges == 0,
                  "cusp surface mesh has edges shared by neither 1 nor 2 triangles (crack or degenerate sliver near the singularity)");
        }
    }

    // --- Test 4: watertightness smoke test -- every edge of the triangle
    // soup should be shared by exactly 2 triangles (closed orientable
    // surface inside the domain), for a simple closed surface fully
    // contained in the box.
    {
        auto ast = parseExpression("x*x + y*y + z*z - 0.16"); // r=0.4, well inside [-1,1]^3
        auto res = Isosurface::build(ast.get(), 0.0, 32, 1.0);
        auto topo = analyzeTopology(res.vertices);
        std::printf("[watertight] %zu unique edges, %d not shared by exactly 2 triangles (of %zu total)\n",
                    topo.numUniqueEdges, topo.boundaryEdges + topo.badEdges, topo.numUniqueEdges);
        CHECK(topo.boundaryEdges == 0 && topo.badEdges == 0, "closed sphere mesh is not watertight (some edge not shared by exactly 2 triangles)");
    }

    // --- Test 5: parser error handling ---
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
