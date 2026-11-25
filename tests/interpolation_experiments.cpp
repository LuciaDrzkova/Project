// Can the singular surface S_0 = {F = 0} be reconstructed by linear
// interpolation between the perturbed surfaces S_-eps and S_+eps?
// (Thesis report, section "Linear interpolation of perturbed surfaces".)
//
// For every reference surface F of the report (all of degree 3, isolated
// singular point at the origin P), a few perturbations G, several eps and the
// grid resolutions 32, 64, 128 on [-0.4, 0.4]^3, the following candidate
// reconstructions of S_0 are compared with the EXACT surface S_0 inside the
// local box Q = [-0.18, 0.18]^3 (distances as in ExactDistance.h):
//
//   endpoint_minus / endpoint_plus  the mesh of S_-eps or S_+eps itself
//   direct_S0                       Marching Cubes of F (the reference mesh)
//   vertex_mid_from_minus           vertices v of mesh(S_-eps) moved to
//                                   (v + c(v)) / 2, where c(v) is the closest
//                                   point of mesh(S_+eps); connectivity of
//                                   mesh(S_-eps) (and symmetrically _from_plus)
//   sdf_average                     zero level of (phi_- + phi_+) / 2, phi_+- the
//                                   signed distance to mesh(S_+-eps) (sign from
//                                   the grid samples of H_+-eps), extracted by
//                                   Marching Cubes on the same grid
//   exact_vertex_mid                as vertex_mid_from_minus, but on the exact
//                                   surfaces: v is first moved to the closest
//                                   point of S_-eps, c(v) is the closest point of
//                                   S_+eps; separates the interpolation error
//                                   from the extraction error of the meshes
//   normal_mid                      as exact_vertex_mid, but the correspondence
//                                   follows the gradient line of F
//   hybrid_cone                     vertex_mid_from_minus away from P, the
//                                   weighted-cone patch over its boundary loop
//                                   near P (no access to F)
//   singular_cone                   Marching Cubes of F away from P, the
//                                   weighted-cone patch projected onto F = 0
//                                   near P (the proposed extraction)
// All reconstructions are implemented in src/analysis/Reconstruction.h.
//
// Interpolating the implicit functions themselves is exact, since H_k = F + kG
// is affine in k: (H_-eps + H_+eps) / 2 = F. The program checks this on the
// grid (column field_average_max_error).
//
//   isosurface_interpolation_experiments [output-directory] [--quick | --slices]
//
// --slices writes, for four illustrative cases (N = 128, eps = 0.01), the
// intersections of the reconstructions with a coordinate plane through P
// (<dir>/slices/*.csv) for the figures of the report.
//
// Output: <dir>/interpolation_all.csv, one row per (F, G, eps, N, method).
#include "analysis/ExactDistance.h"
#include "analysis/LocalAlgebra.h"
#include "analysis/Reconstruction.h"
#include "analysis/StationaryPoints.h"
#include "core/CsvFormat.h"
#include "core/Grid.h"
#include "core/Parallel.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <map>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using iso_geometry::Point3;

namespace {

struct SurfaceSpec {
    std::string folder, type, f;
};

const std::vector<SurfaceSpec> kSurfaces = {
    {"A1_x2+y2-z2+xyz", "A1", "x^2+y^2-z^2+x*y*z"},
    {"A2_x2+y2+z3", "A2", "x^2+y^2+z^3"},
    {"A2_x2-y2+z3", "A2", "x^2-y^2+z^3"},
    {"A2_x3+y3+z2+x2", "A2", "x^3+y^3+z^2+x^2"},
    {"D4_x2+y3+z3", "D4", "x^2+y^3+z^3"},
    {"D4_x2+y3-3yz2", "D4", "x^2+y^3-3*y*z^2"},
    {"E6_x3+y3+z3", "E6~", "x^3+y^3+z^3"},
    {"E6_x3+y3-z3", "E6~", "x^3+y^3-z^3"},
    {"E6_x3+y3+z3+xyz", "E6~", "x^3+y^3+z^3+x*y*z"},
    {"E6_x3+2y3+3z3+2xyz", "E6~", "x^3+2*y^3+3*z^3+2*x*y*z"},
};
const std::vector<std::string> kG = {"1", "x", "x+y", "x^2+y^2-z^2"};

constexpr double kHalfExtent = 0.4;   // meshing domain [-0.4, 0.4]^3
constexpr double kBox = 0.18;         // local box Q
constexpr double kNear = 0.06;        // |p - P|_inf < kNear: "near P", otherwise "away from P"

using iso_reconstruct::Mesh;
using iso_reconstruct::bvh;
using iso_reconstruct::fromResult;
using iso_reconstruct::Method;

Point3 centroidOf(const Mesh& m, const std::array<std::uint32_t, 3>& t) {
    return (m.v[t[0]] + m.v[t[1]] + m.v[t[2]]) * (1.0 / 3.0);
}
double inf(const Point3& p) { return std::max({std::fabs(p.x), std::fabs(p.y), std::fabs(p.z)}); }

// ----- exact closest points ---------------------------------------------------

// The closest point of the exact surface s to a (same candidates as
// iso_exact::distanceTo, but returning the point).
std::optional<Point3> closestOn(const iso_exact::Surface& s, const Point3& a) {
    double best = std::numeric_limits<double>::infinity();
    Point3 bestPoint{};
    auto consider = [&](const Point3& q) {
        const double d = iso_geometry::distance(a, q);
        if (d < best) { best = d; bestPoint = q; }
    };
    for (const Point3& q : s.singularPoints) consider(q);
    auto refine = [&](const Point3& start) {
        consider(start);
        if (const auto p = iso_exact::detail::lagrange(*s.field, a, start, s.step)) consider(*p);
    };
    Point3 start;
    if (s.mesh.nearestPoint(a, start) && iso_exact::detail::project(*s.field, start, s.step)) refine(start);
    if (s.exactPoints.nearestPoint(a, start)) refine(start);
    if (!std::isfinite(best)) return std::nullopt;
    return bestPoint;
}

// ----- one exact surface ------------------------------------------------------

struct Exact {
    std::unique_ptr<iso_exact::Field> field;
    IsosurfaceResult mesh;
    iso_exact::Surface surface;
};

std::unique_ptr<Exact> makeExact(const Node* ast, double k, int resolution, int mu) {
    auto e = std::make_unique<Exact>();
    const double h = 2.0 * kHalfExtent / resolution;
    e->field = std::make_unique<iso_exact::Field>(ast, k);
    e->mesh = Isosurface::build(ast, k, resolution, kHalfExtent);
    std::vector<Point3> singular, real;
    if (const auto cubic = iso_poly::Cubic::fromExpression(ast, k)) {
        iso_stationary::Options so;
        so.regionRadius = std::sqrt(3.0) * kBox;
        const auto r = iso_stationary::findStationaryPoints(*cubic, mu, so);
        for (const auto& p : r.points) {
            if (p.onZeroLevel) singular.push_back(p.point);
            if (p.real) real.push_back(p.point);
        }
    }
    std::vector<Point3> seeds = iso_exact::gradedSeeds(ast, k, {0, 0, 0}, kBox);
    const auto more = iso_exact::stationarySeeds(real, {0, 0, 0}, kBox);
    seeds.insert(seeds.end(), more.begin(), more.end());
    e->surface = iso_exact::prepare(*e->field, e->mesh, {0, 0, 0}, kBox, h, std::move(singular), seeds);
    return e;
}

// ----- measurements -----------------------------------------------------------

struct Errors {
    double all = 0.0, near = 0.0, far = 0.0;   // Hausdorff in Q, split by |p|_inf < kNear
    double meshToS0 = 0.0, s0ToMesh = 0.0;
    double distP = iso_exact::kNaN;             // dist(P, mesh)
    std::size_t failures = 0;
};

using iso_reconstruct::Topology;

// Mesh points in Q: vertices, centroids, and a barycentric grid near P.
std::vector<Point3> meshPoints(const Mesh& m, double h) {
    std::vector<Point3> out;
    for (const auto& t : m.t) {
        const Point3 c = centroidOf(m, t);
        if (inf(c) > kBox) continue;
        if (inf(c) <= 3.0 * h) {
            constexpr int n = 6;
            for (int u = 0; u <= n; ++u)
                for (int w = 0; w <= n - u; ++w)
                    out.push_back(m.v[t[0]] * (double(u) / n) + m.v[t[1]] * (double(w) / n) + m.v[t[2]] * (double(n - u - w) / n));
        } else {
            out.push_back(c);
            if (inf(m.v[t[0]]) <= kBox) out.push_back(m.v[t[0]]);
        }
    }
    return out;
}

Errors measure(const Mesh& m, const iso_exact::Surface& s0, double h) {
    Errors e;
    const auto tree = bvh(m);
    const std::vector<Point3> pts = meshPoints(m, h);
    std::vector<double> d1(pts.size()), d2(s0.samples.size());
    iso_parallel::parallelFor(pts.size(), [&](std::size_t i) { d1[i] = iso_exact::distanceTo(s0, pts[i]); });
    iso_parallel::parallelFor(s0.samples.size(), [&](std::size_t i) { d2[i] = tree.nearestDistance(s0.samples[i]); });
    auto add = [&](const Point3& p, double d, double& directed) {
        if (!std::isfinite(d)) { ++e.failures; return; }
        directed = std::max(directed, d);
        e.all = std::max(e.all, d);
        (inf(p) < kNear ? e.near : e.far) = std::max(inf(p) < kNear ? e.near : e.far, d);
    };
    for (std::size_t i = 0; i < pts.size(); ++i) add(pts[i], d1[i], e.meshToS0);
    for (std::size_t i = 0; i < s0.samples.size(); ++i) add(s0.samples[i], d2[i], e.s0ToMesh);
    if (pts.empty()) e.all = e.near = e.far = e.meshToS0 = iso_exact::kNaN;
    e.distP = tree.nearestDistance({0, 0, 0});
    return e;
}

Topology topology(const Mesh& m, double h, const Mesh* source = nullptr) {
    return iso_reconstruct::topology(m, {0, 0, 0}, kBox, h, source);
}

// ----- reconstructions --------------------------------------------------------

iso_reconstruct::Region nearQ(double h) { return {{0, 0, 0}, kBox + 2.0 * h}; }

iso_reconstruct::ClosestFn closestFn(const iso_exact::Surface& s) {
    return [&s](const Point3& a) { return closestOn(s, a); };
}

Mesh vertexMid(const Mesh& from, const iso_geometry::TriangleBvh& to, double h) {
    return iso_reconstruct::vertexMid(from, to, nearQ(h));
}

Mesh exactVertexMid(const Mesh& from, const iso_exact::Surface& minus, const iso_exact::Surface& plus, double h,
                    std::size_t& failures) {
    return iso_reconstruct::exactVertexMid(from, closestFn(minus), closestFn(plus), nearQ(h), failures);
}

Mesh normalMid(const Mesh& from, const Exact& minus, const Exact& f, const Exact& plus, double h, std::size_t& failures) {
    std::size_t fallbacks = 0;
    Mesh m = iso_reconstruct::normalMid(from, *minus.field, *f.field, *plus.field, closestFn(plus.surface), h, nearQ(h),
                                        fallbacks, failures);
    failures += fallbacks; // reported together: vertices without a gradient-line partner
    return m;
}

// Zero level of the average of the two signed distance fields on the main grid.
Mesh sdfAverage(const Node* familyAst, double eps, int resolution, const iso_geometry::TriangleBvh& minus,
                const iso_geometry::TriangleBvh& plus) {
    const auto hm = Isosurface::sampleField(familyAst, -eps, resolution, kHalfExtent);
    const auto hp = Isosurface::sampleField(familyAst, +eps, resolution, kHalfExtent);
    return iso_reconstruct::sdfAverage(familyAst, resolution, kHalfExtent, hm, hp, minus, plus);
}

constexpr double kPatchCells = 4.0;

Mesh hybridCone(const Mesh& interpolated, const iso_reconstruct::LocalStructure& local, double eps, double h,
                iso_reconstruct::PatchInfo& info) {
    Mesh m = interpolated;
    const double r = iso_reconstruct::hybridRadius(local, eps, h, kPatchCells, 0.6 * kHalfExtent);
    info = iso_reconstruct::conePatch(m, {0, 0, 0}, local.weights, r, h, nullptr);
    return m;
}

Mesh singularCone(const Mesh& directF, const iso_reconstruct::LocalStructure& local, const iso_exact::Field& f, double h,
                  iso_reconstruct::PatchInfo& info) {
    Mesh m = directF;
    info = iso_reconstruct::conePatch(m, {0, 0, 0}, local.weights,
                                      iso_reconstruct::patchRadius(local, h, kPatchCells, 0.6 * kHalfExtent), h, &f);
    return m;
}

std::string fmt(double v) { return iso_csv::number(v); }

} // namespace

// ----- slices for the figures --------------------------------------------------

// Intersection of a mesh with the coordinate plane {p_axis = 0}: segments in the
// two remaining coordinates (in increasing axis order).
void writeSlice(const fs::path& path, const Mesh& m, int axis) {
    std::ofstream o(path);
    o << "u0,v0,u1,v1\n";
    auto comp = [&](const Point3& p, int a) { return a == 0 ? p.x : (a == 1 ? p.y : p.z); };
    const int ua = axis == 0 ? 1 : 0, va = axis == 2 ? 1 : 2;
    for (const auto& t : m.t) {
        std::vector<Point3> cuts;
        for (int c = 0; c < 3; ++c) {
            const Point3& a = m.v[t[static_cast<std::size_t>(c)]];
            const Point3& b = m.v[t[static_cast<std::size_t>((c + 1) % 3)]];
            const double fa = comp(a, axis), fb = comp(b, axis);
            if ((fa < 0.0) != (fb < 0.0)) cuts.push_back(a + (b - a) * (fa / (fa - fb)));
        }
        if (cuts.size() == 2)
            o << fmt(comp(cuts[0], ua)) << ',' << fmt(comp(cuts[0], va)) << ',' << fmt(comp(cuts[1], ua)) << ','
              << fmt(comp(cuts[1], va)) << '\n';
    }
}

struct SliceCase { std::string tag, f, g; double eps; int axis; };

int runSlices(const fs::path& root) {
    const std::vector<SliceCase> cases = {
        {"A1_1", "x^2+y^2-z^2+x*y*z", "1", 1e-2, 0},
        {"A1_x", "x^2+y^2-z^2+x*y*z", "x", 1e-2, 0},
        {"A2a_x", "x^2+y^2+z^3", "x", 1e-2, 1},
        {"E6a_x", "x^3+y^3+z^3", "x", 1e-2, 2},
    };
    const int resolution = 128;
    const double h = 2.0 * kHalfExtent / resolution;
    fs::create_directories(root / "slices");
    for (const SliceCase& c : cases) {
        const NodePtr fAst = parseExpression(c.f);
        const int mu = iso_algebra::localAlgebra(*iso_poly::Cubic::fromExpression(fAst.get()), iso_algebra::IdealKind::Milnor).dimension;
        const NodePtr ast = parseExpression("(" + c.f + ")+k*(" + c.g + ")");
        const auto s0 = makeExact(fAst.get(), 0.0, resolution, mu);
        const auto sm = makeExact(ast.get(), -c.eps, resolution, mu);
        const auto sp = makeExact(ast.get(), +c.eps, resolution, mu);
        const Mesh m0 = fromResult(s0->mesh, h), mm = fromResult(sm->mesh, h), mp = fromResult(sp->mesh, h);
        const auto bm = bvh(mm), bp = bvh(mp);
        const auto local = iso_reconstruct::localStructure(ast.get(), {0, 0, 0});
        std::size_t fails = 0;
        iso_reconstruct::PatchInfo info;
        const Mesh interpolated = vertexMid(mm, bp, h);
        const std::vector<std::pair<std::string, Mesh>> meshes = {
            {"direct_S0", m0},
            {"normal_mid", normalMid(mm, *sm, *s0, *sp, h, fails)},
            {"hybrid_cone", hybridCone(interpolated, local, c.eps, h, info)},
            {"singular_cone", singularCone(m0, local, *s0->field, h, info)},
            {"vertex_mid_from_minus", vertexMid(mm, bp, h)},
            {"vertex_mid_from_plus", vertexMid(mp, bm, h)},
            {"sdf_average", sdfAverage(ast.get(), c.eps, resolution, bm, bp)},
            {"exact_vertex_mid", exactVertexMid(mm, sm->surface, sp->surface, h, fails)},
        };
        for (const auto& [name, mesh] : meshes) writeSlice(root / "slices" / (c.tag + "_" + name + ".csv"), mesh, c.axis);
        std::printf("slices %s written\n", c.tag.c_str());
    }
    return 0;
}

int main(int argc, char** argv) {
    try {
        fs::path root = "interpolation_results";
        bool quick = false, slices = false;
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--quick") quick = true;
            else if (a == "--slices") slices = true;
            else root = a;
        }
        fs::create_directories(root);
        if (slices) return runSlices(root);
        const std::vector<int> resolutions = quick ? std::vector<int>{32} : std::vector<int>{32, 64, 128};
        const std::vector<double> epsilons = quick ? std::vector<double>{1e-2, 1e-3}
                                                   : std::vector<double>{1e-2, std::pow(10.0, -2.5), 1e-3, std::pow(10.0, -3.5), 1e-4};
        std::ofstream out(root / "interpolation_all.csv");
        out << "surface,type,F,G,eps,resolution,grid_spacing,method,hausdorff_Q,hausdorff_near_P,hausdorff_away_from_P,"
               "mesh_to_S0_max,S0_to_mesh_max,dist_P_to_mesh,distance_failures,triangles_in_Q,components,euler,"
               "hole_edges,nonmanifold_edges,degenerate_triangles,flipped_triangles,field_average_max_error,"
               "correspondence_failures,patch_radius,patch_loops,patch_open_chains,patch_projection_failures\n";
        const auto t0 = std::chrono::steady_clock::now();

        for (const SurfaceSpec& s : kSurfaces) {
            const NodePtr fAst = parseExpression(s.f);
            const int mu = iso_algebra::localAlgebra(*iso_poly::Cubic::fromExpression(fAst.get()), iso_algebra::IdealKind::Milnor).dimension;
            for (int resolution : resolutions) {
                const double h = 2.0 * kHalfExtent / resolution;
                const auto s0 = makeExact(fAst.get(), 0.0, resolution, mu);
                const Mesh m0 = fromResult(s0->mesh, h);
                const Errors e0 = measure(m0, s0->surface, h);
                const Topology tp0 = topology(m0, h);
                for (const std::string& g : kG) {
                    const NodePtr ast = parseExpression("(" + s.f + ")+k*(" + g + ")");
                    const auto local = iso_reconstruct::localStructure(ast.get(), {0, 0, 0});
                    for (double eps : epsilons) {
                        const auto sm = makeExact(ast.get(), -eps, resolution, mu);
                        const auto sp = makeExact(ast.get(), +eps, resolution, mu);
                        const Mesh mm = fromResult(sm->mesh, h), mp = fromResult(sp->mesh, h);
                        const auto bm = bvh(mm), bp = bvh(mp);

                        // (H_-eps + H_+eps)/2 - F on the grid: zero up to rounding.
                        double fieldError = 0.0;
                        {
                            const auto a = Isosurface::sampleField(ast.get(), -eps, resolution / 4, kHalfExtent);
                            const auto b = Isosurface::sampleField(ast.get(), +eps, resolution / 4, kHalfExtent);
                            const auto f = Isosurface::sampleField(fAst.get(), 0.0, resolution / 4, kHalfExtent);
                            for (std::size_t i = 0; i < f.size(); ++i) fieldError = std::max(fieldError, std::fabs(0.5 * (a[i] + b[i]) - f[i]));
                        }

                        struct Row {
                            std::string method;
                            Mesh mesh;
                            const Mesh* source;
                            std::size_t corrFailures;
                            const iso_reconstruct::PatchInfo* patch = nullptr;
                        };
                        std::vector<Row> rows;
                        rows.reserve(9); // rows[2] is referenced below
                        rows.push_back({"endpoint_minus", mm, nullptr, 0});
                        rows.push_back({"endpoint_plus", mp, nullptr, 0});
                        rows.push_back({"vertex_mid_from_minus", vertexMid(mm, bp, h), &mm, 0});
                        rows.push_back({"vertex_mid_from_plus", vertexMid(mp, bm, h), &mp, 0});
                        rows.push_back({"sdf_average", sdfAverage(ast.get(), eps, resolution, bm, bp), nullptr, 0});
                        std::size_t fails = 0;
                        Mesh ex = exactVertexMid(mm, sm->surface, sp->surface, h, fails);
                        rows.push_back({"exact_vertex_mid", std::move(ex), &mm, fails});
                        std::size_t normalFails = 0;
                        Mesh nm = normalMid(mm, *sm, *s0, *sp, h, normalFails);
                        rows.push_back({"normal_mid", std::move(nm), &mm, normalFails});
                        iso_reconstruct::PatchInfo hybridInfo, coneInfo;
                        rows.push_back({"hybrid_cone", hybridCone(rows[2].mesh, local, eps, h, hybridInfo), nullptr, 0,
                                        &hybridInfo});
                        rows.push_back({"singular_cone", singularCone(m0, local, *s0->field, h, coneInfo), nullptr, 0,
                                        &coneInfo});

                        auto write = [&](const std::string& method, const Errors& e, const Topology& tp, std::size_t corr,
                                         const iso_reconstruct::PatchInfo* patch) {
                            out << iso_csv::quote(s.folder) << ',' << iso_csv::quote(s.type) << ',' << iso_csv::quote(s.f) << ','
                                << iso_csv::quote(g) << ',' << fmt(eps) << ',' << resolution << ',' << fmt(h) << ',' << method << ','
                                << fmt(e.all) << ',' << fmt(e.near) << ',' << fmt(e.far) << ',' << fmt(e.meshToS0) << ','
                                << fmt(e.s0ToMesh) << ',' << fmt(e.distP) << ',' << e.failures << ',' << tp.triangles << ','
                                << tp.components << ',' << tp.euler << ',' << tp.holeEdges << ',' << tp.nonManifoldEdges << ','
                                << tp.degenerate << ',';
                            if (tp.flipped == static_cast<std::size_t>(-1)) out << ""; else out << tp.flipped;
                            out << ',' << fmt(fieldError) << ',' << corr << ',';
                            if (patch && patch->applied)
                                out << fmt(patch->radius) << ',' << patch->loops << ',' << patch->openChains << ','
                                    << patch->projectionFailures;
                            else
                                out << ",,,";
                            out << '\n';
                        };
                        write("direct_S0", e0, tp0, 0, nullptr);
                        for (const Row& r : rows)
                            write(r.method, measure(r.mesh, s0->surface, h), topology(r.mesh, h, r.source), r.corrFailures, r.patch);
                        out.flush();
                    }
                    const double sec = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                    std::printf("%s N=%d G=%s done (%.0f s)\n", s.f.c_str(), resolution, g.c_str(), sec);
                    std::fflush(stdout);
                }
            }
        }
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "interpolation experiments failed: " << e.what() << '\n';
        return 1;
    }
}
