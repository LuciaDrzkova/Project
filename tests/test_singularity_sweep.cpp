// Sweep all built-in surface families over the requested k values.
// Headless test: depends only on Expr.h and Isosurface.h, not on GLFW/ImGui/PMP.
// It reports (rather than hiding) cases where no regular, topology-clean pair exists.

#include "Expr.h"
#include "Isosurface.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace {

constexpr int kDefaultResolution = 48;
int gResolution = kDefaultResolution;
constexpr double kHalfExtent = 1.0;
constexpr double kPreviewRadius = 0.45;
constexpr double kGradientEpsilon = 1e-9;
constexpr double kZeroTolerance = 1e-9;
constexpr double kWeldScale = 1'000'000.0;

// Preserve the order requested in the issue so the CSV is easy to compare
// against the user's experiment notes.
constexpr std::array<double, 9> kValues = {
    -0.0001, -0.001, -0.01, -0.1, 0.0, 0.1, 0.01, 0.001, 0.0001};

enum class Behavior {
    RegularizedByNonzeroK, // the intended isolated singularity disappears for k != 0
    WhitneyUmbrella,       // k > 0 moves singularities to (+/-sqrt(k), 0, 0)
    CubicFamily            // the origin remains singular for all k
};

struct SurfaceSpec {
    const char* name;
    const char* slug;
    const char* expression;
    const char* singularity_description;
    Behavior behavior;
};

// Keep in sync with the 17 selectable formulas in src/main.cpp.
const std::array<SurfaceSpec, 17> kSurfaces = {{
    {"Plane example: x^3-k*y", "plane_x3_ky", "x*x*x - k*y",
     "At k=0 the singular set is the plane x=0 (not an isolated point).", Behavior::RegularizedByNonzeroK},
    {"Isolated cusp X", "isolated_cusp_x", "x*x*x + y*y + z*z + k*x",
     "At k=0 the isolated real singular point is (0,0,0).", Behavior::RegularizedByNonzeroK},
    {"Isolated cusp Y", "isolated_cusp_y", "x*x + y*y*y + z*z + k*y",
     "At k=0 the isolated real singular point is (0,0,0).", Behavior::RegularizedByNonzeroK},
    {"Isolated cusp Z", "isolated_cusp_z", "x*x + y*y + z*z*z + k*z",
     "At k=0 the isolated real singular point is (0,0,0).", Behavior::RegularizedByNonzeroK},
    {"Quartic cusp X", "quartic_cusp_x", "x*x*x*x + y*y + z*z + k*x",
     "At k=0 the isolated real singular point is (0,0,0).", Behavior::RegularizedByNonzeroK},
    {"Quartic cusp Y", "quartic_cusp_y", "x*x + y*y*y*y + z*z + k*y",
     "At k=0 the isolated real singular point is (0,0,0).", Behavior::RegularizedByNonzeroK},
    {"Quartic cusp Z", "quartic_cusp_z", "x*x + y*y + z*z*z*z + k*z",
     "At k=0 the isolated real singular point is (0,0,0).", Behavior::RegularizedByNonzeroK},
    {"Higher-order cusp", "higher_order_cusp", "x*x*x*x*x + y*y + z*z + k*x",
     "At k=0 the isolated real singular point is (0,0,0).", Behavior::RegularizedByNonzeroK},
    {"Mixed cusp Y", "mixed_cusp_y", "x*x*x + y*y*y*y + z*z + k*y",
     "At k=0 the isolated real singular point is (0,0,0).", Behavior::RegularizedByNonzeroK},
    {"Mixed cusp Z", "mixed_cusp_z", "x*x*x + y*y + z*z*z*z + k*z",
     "At k=0 the isolated real singular point is (0,0,0).", Behavior::RegularizedByNonzeroK},
    {"Mixed cusp X", "mixed_cusp_x", "x*x + y*y*y + z*z*z*z + k*x",
     "At k=0 the isolated real singular point is (0,0,0).", Behavior::RegularizedByNonzeroK},
    {"Mixed-order Z", "mixed_order_z", "x*x*x*x + y*y*y + z*z + k*z",
     "At k=0 the isolated real singular point is (0,0,0).", Behavior::RegularizedByNonzeroK},
    {"Mixed cubic", "mixed_cubic", "x*x*x + y*y*y + z*z + k*z",
     "At k=0 the isolated real singular point is (0,0,0).", Behavior::RegularizedByNonzeroK},
    {"Level-shift cone", "level_shift_cone", "x*x + y*y - z*z - k",
     "The cone point is the origin at k=0; for k!=0 the origin is not on F=0.", Behavior::RegularizedByNonzeroK},
    {"Whitney umbrella", "whitney_umbrella", "x*x*z - y*y - k*z",
     "At k=0 the singular set is the z-axis; for k>0 singular points are (+/-sqrt(k),0,0).", Behavior::WhitneyUmbrella},
    {"Monkey saddle", "monkey_saddle", "x*x*x - 3*x*y*y - k",
     "At k=0 the singular set is the z-axis; for k!=0 that line is not on F=0.", Behavior::RegularizedByNonzeroK},
    {"Cubic family", "cubic_family", "x*x*x + y*y*y + z*z*z - 3*k*x*y*z",
     "The origin remains a singular point on F=0 for every k; k=1 has additional special behavior.", Behavior::CubicFamily}
}};

struct Point3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct Eval {
    double f = std::numeric_limits<double>::quiet_NaN();
    double gx = std::numeric_limits<double>::quiet_NaN();
    double gy = std::numeric_limits<double>::quiet_NaN();
    double gz = std::numeric_limits<double>::quiet_NaN();
    double grad = std::numeric_limits<double>::quiet_NaN();
};

Eval evaluate(const Node* ast, const Point3& p, double k)
{
    const Dual3 x(p.x, 1.0, 0.0, 0.0);
    const Dual3 y(p.y, 0.0, 1.0, 0.0);
    const Dual3 z(p.z, 0.0, 0.0, 1.0);
    const Dual3 kd(k, 0.0, 0.0, 0.0);
    const Dual3 result = evalExpr<Dual3>(ast, x, y, z, kd);
    Eval out;
    out.f = result.v;
    out.gx = result.dx;
    out.gy = result.dy;
    out.gz = result.dz;
    out.grad = std::sqrt(out.gx * out.gx + out.gy * out.gy + out.gz * out.gz);
    return out;
}

bool finiteEval(const Eval& e)
{
    return std::isfinite(e.f) && std::isfinite(e.gx) && std::isfinite(e.gy) &&
           std::isfinite(e.gz) && std::isfinite(e.grad);
}

bool analyticSurfaceRegular(const SurfaceSpec& surface, double k)
{
    if (std::fabs(k) <= 1e-15)
        return false; // k=0 is the distinguished singular member for all presets.
    if (surface.behavior == Behavior::CubicFamily)
        return false; // grad F(0,0,0)=0 and F(0,0,0)=0 for every k.
    if (surface.behavior == Behavior::WhitneyUmbrella && k > 0.0)
        return false; // two singular points remain on F=0 at +/-sqrt(k).
    return true;
}

std::vector<std::pair<std::string, Point3>> singularityProbePoints(
    const SurfaceSpec& surface, double k)
{
    std::vector<std::pair<std::string, Point3>> points;
    points.push_back({"origin_reference", {0.0, 0.0, 0.0}});
    if (surface.behavior == Behavior::WhitneyUmbrella && k > 0.0) {
        const double sx = std::sqrt(k);
        points.push_back({"whitney_singularity_plus", {sx, 0.0, 0.0}});
        points.push_back({"whitney_singularity_minus", {-sx, 0.0, 0.0}});
    }
    return points;
}

using PointKey = std::array<long long, 3>;
using EdgeKey = std::pair<PointKey, PointKey>;

PointKey pointKey(const SurfaceVertex& v)
{
    return {{std::llround(static_cast<double>(v.px) * kWeldScale),
             std::llround(static_cast<double>(v.py) * kWeldScale),
             std::llround(static_cast<double>(v.pz) * kWeldScale)}};
}

EdgeKey edgeKey(PointKey a, PointKey b)
{
    if (b < a) std::swap(a, b);
    return {a, b};
}

struct TopologyMetrics {
    std::size_t triangles = 0;
    std::size_t unique_vertices = 0;
    std::size_t unique_edges = 0;
    std::size_t connected_components = 0;
    long long euler_characteristic = 0;
    std::size_t degenerate_triangles = 0;
    std::size_t domain_boundary_edges = 0;
    std::size_t interior_boundary_edges = 0;
    std::size_t nonmanifold_edges = 0;
    std::size_t orientation_conflict_edges = 0;
    bool clean = false;
};

bool onSameDomainFace(const SurfaceVertex& a, const SurfaceVertex& b)
{
    const double tol = 4e-5 * std::max(1.0, kHalfExtent);
    const double av[3] = {a.px, a.py, a.pz};
    const double bv[3] = {b.px, b.py, b.pz};
    for (int axis = 0; axis < 3; ++axis) {
        const bool aMinus = std::fabs(av[axis] + kHalfExtent) <= tol;
        const bool bMinus = std::fabs(bv[axis] + kHalfExtent) <= tol;
        const bool aPlus = std::fabs(av[axis] - kHalfExtent) <= tol;
        const bool bPlus = std::fabs(bv[axis] - kHalfExtent) <= tol;
        if ((aMinus && bMinus) || (aPlus && bPlus)) return true;
    }
    return false;
}

TopologyMetrics analyzeTopology(const IsosurfaceResult& mesh)
{
    TopologyMetrics out;
    out.triangles = mesh.vertices.size() / 3;
    std::map<PointKey, SurfaceVertex> pointReps;
    std::map<EdgeKey, std::size_t> edgeUseCount;
    std::map<EdgeKey, int> orientationBalance;
    const double step = (2.0 * kHalfExtent) / gResolution;
    const double areaTolerance = step * step * 1e-8;

    for (std::size_t t = 0; t + 2 < mesh.vertices.size(); t += 3) {
        const SurfaceVertex& a = mesh.vertices[t];
        const SurfaceVertex& b = mesh.vertices[t + 1];
        const SurfaceVertex& c = mesh.vertices[t + 2];
        const double e1x = static_cast<double>(b.px) - a.px;
        const double e1y = static_cast<double>(b.py) - a.py;
        const double e1z = static_cast<double>(b.pz) - a.pz;
        const double e2x = static_cast<double>(c.px) - a.px;
        const double e2y = static_cast<double>(c.py) - a.py;
        const double e2z = static_cast<double>(c.pz) - a.pz;
        const double cx = e1y * e2z - e1z * e2y;
        const double cy = e1z * e2x - e1x * e2z;
        const double cz = e1x * e2y - e1y * e2x;
        const double area = 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
        if (!std::isfinite(area) || area <= areaTolerance)
            ++out.degenerate_triangles;

        const SurfaceVertex tri[3] = {a, b, c};
        PointKey keys[3];
        for (int j = 0; j < 3; ++j) {
            keys[j] = pointKey(tri[j]);
            pointReps.emplace(keys[j], tri[j]);
        }
        for (int edge = 0; edge < 3; ++edge) {
            const PointKey& from = keys[edge];
            const PointKey& to = keys[(edge + 1) % 3];
            const EdgeKey key = edgeKey(from, to);
            ++edgeUseCount[key];
            orientationBalance[key] += (from < to ? 1 : -1);
        }
    }

    out.unique_vertices = pointReps.size();
    out.unique_edges = edgeUseCount.size();
    out.euler_characteristic = static_cast<long long>(out.unique_vertices) -
        static_cast<long long>(out.unique_edges) + static_cast<long long>(out.triangles);

    // Count edge-connected components and orientation conflicts using welded
    // geometric coordinates, rather than the extractor's duplicated triangle soup.
    std::map<PointKey, std::size_t> vertexId;
    std::size_t nextId = 0;
    for (const auto& entry : pointReps) vertexId.emplace(entry.first, nextId++);
    std::vector<std::size_t> parent(nextId), rank(nextId, 0);
    for (std::size_t i = 0; i < nextId; ++i) parent[i] = i;
    auto findRoot = [&](auto&& self, std::size_t x) -> std::size_t {
        if (parent[x] != x) parent[x] = self(self, parent[x]);
        return parent[x];
    };
    auto unite = [&](std::size_t a, std::size_t b) {
        a = findRoot(findRoot, a);
        b = findRoot(findRoot, b);
        if (a == b) return;
        if (rank[a] < rank[b]) std::swap(a, b);
        parent[b] = a;
        if (rank[a] == rank[b]) ++rank[a];
    };
    for (const auto& entry : edgeUseCount)
        unite(vertexId.at(entry.first.first), vertexId.at(entry.first.second));
    std::map<std::size_t, bool> roots;
    for (std::size_t i = 0; i < nextId; ++i) roots[findRoot(findRoot, i)] = true;
    out.connected_components = roots.size();

    for (const auto& entry : edgeUseCount) {
        const std::size_t uses = entry.second;
        if (uses == 1) {
            const auto ia = pointReps.find(entry.first.first);
            const auto ib = pointReps.find(entry.first.second);
            if (ia != pointReps.end() && ib != pointReps.end() &&
                onSameDomainFace(ia->second, ib->second))
                ++out.domain_boundary_edges; // expected where the crop box cuts a surface
            else
                ++out.interior_boundary_edges; // unexpected hole/crack inside the domain
        } else if (uses > 2) {
            ++out.nonmanifold_edges;
        } else if (uses == 2 && orientationBalance[entry.first] != 0) {
            ++out.orientation_conflict_edges;
        }
    }

    out.clean = out.triangles > 0 && out.degenerate_triangles == 0 &&
                out.interior_boundary_edges == 0 && out.nonmanifold_edges == 0 &&
                out.orientation_conflict_edges == 0;
    return out;
}

struct CaseMetrics {
    const SurfaceSpec* surface = nullptr;
    double k = 0.0;
    std::size_t vertex_count = 0;
    double build_millis = 0.0;
    TopologyMetrics topology;
    Eval origin;
    std::size_t near_vertex_count = 0;
    double nearest_vertex_radius = std::numeric_limits<double>::infinity();
    double near_grad_min = std::numeric_limits<double>::quiet_NaN();
    double near_grad_p10 = std::numeric_limits<double>::quiet_NaN();
    double near_grad_median = std::numeric_limits<double>::quiet_NaN();
    bool finite_vertices = true;
    bool analytic_regular = false;
    bool candidate_ok = false;
    std::string reason;
};

CaseMetrics analyzeCase(const SurfaceSpec& surface, const Node* ast, double k)
{
    CaseMetrics row;
    row.surface = &surface;
    row.k = k;
    IsosurfaceResult mesh = Isosurface::build(ast, k, gResolution, kHalfExtent);
    row.vertex_count = mesh.vertices.size();
    row.build_millis = mesh.buildMillis;
    for (const SurfaceVertex& v : mesh.vertices) {
        if (!std::isfinite(v.px) || !std::isfinite(v.py) || !std::isfinite(v.pz) ||
            !std::isfinite(v.nx) || !std::isfinite(v.ny) || !std::isfinite(v.nz) ||
            !std::isfinite(v.gradMag) || v.gradMag < 0.0f) {
            row.finite_vertices = false;
            break;
        }
    }
    row.topology = row.finite_vertices ? analyzeTopology(mesh) : TopologyMetrics{};
    row.origin = evaluate(ast, {0.0, 0.0, 0.0}, k);
    row.analytic_regular = analyticSurfaceRegular(surface, k);

    std::vector<double> localGradients;
    localGradients.reserve(mesh.vertices.size());
    for (const SurfaceVertex& v : mesh.vertices) {
        const double d = std::sqrt(static_cast<double>(v.px) * v.px +
                                   static_cast<double>(v.py) * v.py +
                                   static_cast<double>(v.pz) * v.pz);
        row.nearest_vertex_radius = std::min(row.nearest_vertex_radius, d);
        if (d <= kPreviewRadius) {
            ++row.near_vertex_count;
            if (std::isfinite(v.gradMag)) localGradients.push_back(v.gradMag);
        }
    }
    if (!localGradients.empty()) {
        std::sort(localGradients.begin(), localGradients.end());
        row.near_grad_min = localGradients.front();
        row.near_grad_p10 = localGradients[static_cast<std::size_t>(
            0.10 * static_cast<double>(localGradients.size() - 1))];
        const std::size_t mid = localGradients.size() / 2;
        row.near_grad_median = localGradients.size() % 2 == 0
            ? 0.5 * (localGradients[mid - 1] + localGradients[mid])
            : localGradients[mid];
    }

    std::vector<std::string> reasons;
    if (std::fabs(k) <= 1e-15) reasons.push_back("k=0 singular member, not a regularization candidate");
    if (!row.analytic_regular) reasons.push_back("an analytic singularity remains on F=0");
    if (!row.finite_vertices) reasons.push_back("non-finite/invalid mesh vertex data");
    if (mesh.vertices.empty()) reasons.push_back("no surface triangles in the sampled domain");
    if (row.near_vertex_count == 0) reasons.push_back("no mesh vertices within local preview radius");
    if (!localGradients.empty() && !(row.near_grad_min > kGradientEpsilon))
        reasons.push_back("a local surface vertex has zero/nearly-zero gradient");
    if (!row.topology.clean) {
        std::ostringstream topo;
        topo << "topology defects (degenerate=" << row.topology.degenerate_triangles
             << ", interior-boundary=" << row.topology.interior_boundary_edges
             << ", nonmanifold=" << row.topology.nonmanifold_edges
             << ", orientation-conflicts=" << row.topology.orientation_conflict_edges << ')';
        reasons.push_back(topo.str());
    }
    row.candidate_ok = std::fabs(k) > 1e-15 && row.analytic_regular &&
        row.finite_vertices && row.vertex_count > 0 && row.near_vertex_count > 0 &&
        std::isfinite(row.near_grad_min) && row.near_grad_min > kGradientEpsilon &&
        row.topology.clean;
    if (reasons.empty()) row.reason = "eligible";
    else {
        for (std::size_t i = 0; i < reasons.size(); ++i) {
            if (i) row.reason += "; ";
            row.reason += reasons[i];
        }
    }
    return row;
}

std::string csvQuote(const std::string& s)
{
    std::string out = "\"";
    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out += c;
    }
    out += '"';
    return out;
}

void writeCaseCsv(const std::filesystem::path& path, const std::vector<CaseMetrics>& rows)
{
    std::ofstream out(path);
    out << "surface,expression,k,triangles,triangle_soup_vertices,unique_vertices,unique_edges,connected_components,euler_characteristic,build_ms,origin_F,origin_gx,origin_gy,origin_gz,origin_grad_mag,origin_on_F0,analytic_regular_on_surface,near_vertex_count,nearest_vertex_radius,near_grad_min,near_grad_p10,near_grad_median,degenerate_triangles,domain_boundary_edges,interior_boundary_edges,nonmanifold_edges,orientation_conflict_edges,topology_clean,candidate_ok,decision_reason\n";
    out << std::setprecision(12);
    for (const auto& r : rows) {
        out << csvQuote(r.surface->name) << ',' << csvQuote(r.surface->expression) << ',' << r.k << ','
            << r.topology.triangles << ',' << r.vertex_count << ','
            << r.topology.unique_vertices << ',' << r.topology.unique_edges << ','
            << r.topology.connected_components << ',' << r.topology.euler_characteristic << ','
            << r.build_millis << ','
            << r.origin.f << ',' << r.origin.gx << ',' << r.origin.gy << ',' << r.origin.gz << ','
            << r.origin.grad << ',' << (std::fabs(r.origin.f) <= kZeroTolerance ? 1 : 0) << ','
            << (r.analytic_regular ? 1 : 0) << ',' << r.near_vertex_count << ',';
        if (std::isfinite(r.nearest_vertex_radius)) out << r.nearest_vertex_radius;
        out << ',';
        if (std::isfinite(r.near_grad_min)) out << r.near_grad_min;
        out << ',';
        if (std::isfinite(r.near_grad_p10)) out << r.near_grad_p10;
        out << ',';
        if (std::isfinite(r.near_grad_median)) out << r.near_grad_median;
        out << ',' << r.topology.degenerate_triangles << ',' << r.topology.domain_boundary_edges << ','
            << r.topology.interior_boundary_edges << ',' << r.topology.nonmanifold_edges << ','
            << r.topology.orientation_conflict_edges << ','
            << (r.topology.clean ? 1 : 0) << ',' << (r.candidate_ok ? 1 : 0) << ','
            << csvQuote(r.reason) << '\n';
    }
}

void writeSingularityPointCsv(const std::filesystem::path& path)
{
    std::ofstream out(path);
    out << "surface,k,probe_name,x,y,z,F,gx,gy,gz,grad_mag,on_F0,stationary_on_F0\n";
    out << std::setprecision(12);
    for (const auto& s : kSurfaces) {
        NodePtr ast = parseExpression(s.expression);
        for (double k : kValues) {
            for (const auto& candidate : singularityProbePoints(s, k)) {
                const Eval e = evaluate(ast.get(), candidate.second, k);
                const bool onSurface = std::fabs(e.f) <= kZeroTolerance;
                const bool stationary = onSurface && e.grad <= kGradientEpsilon;
                out << csvQuote(s.name) << ',' << k << ',' << csvQuote(candidate.first) << ','
                    << candidate.second.x << ',' << candidate.second.y << ',' << candidate.second.z << ','
                    << e.f << ',' << e.gx << ',' << e.gy << ',' << e.gz << ',' << e.grad << ','
                    << (onSurface ? 1 : 0) << ',' << (stationary ? 1 : 0) << '\n';
            }
        }
    }
}

const CaseMetrics* chooseBest(const std::vector<CaseMetrics>& rows,
                              const SurfaceSpec& surface, bool wantPositive)
{
    const CaseMetrics* best = nullptr;
    for (const CaseMetrics& row : rows) {
        if (row.surface != &surface || !row.candidate_ok) continue;
        if (wantPositive != (row.k > 0.0)) continue;
        if (!best || std::fabs(row.k) < std::fabs(best->k)) best = &row;
    }
    return best;
}

bool writeLocalPairPly(const std::filesystem::path& path,
                       const SurfaceSpec& surface,
                       const CaseMetrics* negative,
                       const CaseMetrics* positive)
{
    struct PlyVertex { float x, y, z, nx, ny, nz, grad; unsigned char r, g, b; };
    std::vector<PlyVertex> vertices;
    std::vector<std::array<std::uint32_t, 3>> faces;

    auto appendMesh = [&](double k, unsigned char r, unsigned char g, unsigned char b) {
        NodePtr ast = parseExpression(surface.expression);
        IsosurfaceResult mesh = Isosurface::build(ast.get(), k, gResolution, kHalfExtent);
        for (std::size_t i = 0; i + 2 < mesh.vertices.size(); i += 3) {
            const SurfaceVertex tri[3] = {mesh.vertices[i], mesh.vertices[i + 1], mesh.vertices[i + 2]};
            bool local = true;
            for (const SurfaceVertex& v : tri) {
                const double d = std::sqrt(static_cast<double>(v.px) * v.px +
                                           static_cast<double>(v.py) * v.py +
                                           static_cast<double>(v.pz) * v.pz);
                if (d > kPreviewRadius) { local = false; break; }
            }
            if (!local) continue;
            const std::uint32_t base = static_cast<std::uint32_t>(vertices.size());
            for (const SurfaceVertex& v : tri)
                vertices.push_back({v.px, v.py, v.pz, v.nx, v.ny, v.nz, v.gradMag, r, g, b});
            faces.push_back({base, static_cast<std::uint32_t>(base + 1),
                             static_cast<std::uint32_t>(base + 2)});
        }
    };

    if (negative) appendMesh(negative->k, 55, 150, 255);  // negative k: blue
    if (positive) appendMesh(positive->k, 255, 105, 55); // positive k: orange-red
    if (faces.empty()) return false;

    std::ofstream out(path);
    if (!out) return false;
    out << "ply\nformat ascii 1.0\ncomment " << surface.name
        << " near-singularity preview; blue=negative k, orange=positive k\n"
        << "element vertex " << vertices.size() << "\n"
        << "property float x\nproperty float y\nproperty float z\n"
        << "property float nx\nproperty float ny\nproperty float nz\n"
        << "property float grad_mag\n"
        << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
        << "element face " << faces.size() << "\n"
        << "property list uchar uint vertex_indices\nend_header\n";
    out << std::setprecision(8);
    for (const auto& v : vertices)
        out << v.x << ' ' << v.y << ' ' << v.z << ' '
            << v.nx << ' ' << v.ny << ' ' << v.nz << ' ' << v.grad << ' '
            << static_cast<int>(v.r) << ' ' << static_cast<int>(v.g) << ' '
            << static_cast<int>(v.b) << '\n';
    for (const auto& f : faces)
        out << "3 " << f[0] << ' ' << f[1] << ' ' << f[2] << '\n';
    return true;
}

// Local common-grid samples near F=0 for the selected negative/positive pair.
// The scalar field and its AD gradient are evaluated at identical (x,y,z)
// coordinates for k-, k=0 and k+. This avoids falsely assuming that vertices
// from different meshes have a one-to-one correspondence.
struct LocalInterpolationStats {
    std::size_t sample_count = 0;
    double max_abs_f_error = 0.0;
    double max_gradient_error = 0.0;
};

struct LocalGridNode {
    Point3 p;
    std::array<Eval, 3> values; // negative winner, k=0 reference, positive winner
    std::array<bool, 3> near_surface{{false, false, false}};
};

using GridIndex = std::array<int, 3>;

LocalInterpolationStats writeNearSurfaceSamples(
    const std::filesystem::path& path, const SurfaceSpec& surface,
    const Node* ast, double kNegative, double kPositive)
{
    LocalInterpolationStats stats;
    const double step = (2.0 * kHalfExtent) / static_cast<double>(gResolution);
    const double roi = kPreviewRadius + 2.0 * step;
    const int lower = std::max(0, static_cast<int>(std::floor((kHalfExtent - roi) / step)));
    const int upper = std::min(gResolution, static_cast<int>(std::ceil((kHalfExtent + roi) / step)));
    const std::array<double, 3> ks{{kNegative, 0.0, kPositive}};
    const double interpDenom = kPositive - kNegative;
    if (!(interpDenom > 0.0)) return stats;
    const double weightNegative = kPositive / interpDenom;
    const double weightPositive = -kNegative / interpDenom;

    std::map<GridIndex, LocalGridNode> nodes;
    for (int i = lower; i <= upper; ++i) {
        const double x = -kHalfExtent + step * i;
        for (int j = lower; j <= upper; ++j) {
            const double y = -kHalfExtent + step * j;
            for (int l = lower; l <= upper; ++l) {
                const double z = -kHalfExtent + step * l;
                LocalGridNode node;
                node.p = {x, y, z};
                for (int ki = 0; ki < 3; ++ki)
                    node.values[static_cast<std::size_t>(ki)] = evaluate(ast, node.p, ks[static_cast<std::size_t>(ki)]);
                nodes.emplace(GridIndex{{i, j, l}}, node);
            }
        }
    }

    const int cornerOffsets[8][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}};
    auto minAbsInterval = [](double a, double b) {
        if (a <= 0.0 && b >= 0.0) return 0.0;
        return std::min(std::fabs(a), std::fabs(b));
    };

    for (int i = lower; i < upper; ++i) {
        const double x0 = -kHalfExtent + step * i, x1 = x0 + step;
        for (int j = lower; j < upper; ++j) {
            const double y0 = -kHalfExtent + step * j, y1 = y0 + step;
            for (int l = lower; l < upper; ++l) {
                const double z0 = -kHalfExtent + step * l, z1 = z0 + step;
                const double dx = minAbsInterval(x0, x1);
                const double dy = minAbsInterval(y0, y1);
                const double dz = minAbsInterval(z0, z1);
                if (dx * dx + dy * dy + dz * dz >
                    (kPreviewRadius + 0.25 * step) * (kPreviewRadius + 0.25 * step))
                    continue;

                std::array<GridIndex, 8> keys;
                std::array<const LocalGridNode*, 8> corners{};
                for (int c = 0; c < 8; ++c) {
                    keys[static_cast<std::size_t>(c)] = GridIndex{{
                        i + cornerOffsets[c][0], j + cornerOffsets[c][1], l + cornerOffsets[c][2]}};
                    corners[static_cast<std::size_t>(c)] = &nodes.at(keys[static_cast<std::size_t>(c)]);
                }
                for (int ki = 0; ki < 3; ++ki) {
                    double lo = std::numeric_limits<double>::infinity();
                    double hi = -std::numeric_limits<double>::infinity();
                    bool finite = true;
                    for (const LocalGridNode* corner : corners) {
                        const double f = corner->values[static_cast<std::size_t>(ki)].f;
                        if (!std::isfinite(f)) { finite = false; break; }
                        lo = std::min(lo, f);
                        hi = std::max(hi, f);
                    }
                    // Exactly the cells the marching-cubes sign test can cut.
                    if (finite && lo <= 0.0 && hi >= 0.0) {
                        for (const GridIndex& key : keys)
                            nodes.at(key).near_surface[static_cast<std::size_t>(ki)] = true;
                    }
                }
            }
        }
    }

    std::ofstream out(path);
    if (!out) return stats;
    out << "surface,expression,k_negative,k_zero,k_positive,sample_kind,near_negative_surface,near_zero_surface,near_positive_surface,x,y,z,F_negative,F_zero_direct,F_positive,F_zero_interpolated,F_interpolation_error,gx_negative,gy_negative,gz_negative,gradmag_negative,gx_zero_direct,gy_zero_direct,gz_zero_direct,gradmag_zero_direct,gx_positive,gy_positive,gz_positive,gradmag_positive,gx_zero_interpolated,gy_zero_interpolated,gz_zero_interpolated,gradmag_zero_interpolated,gradient_interpolation_error,min_abs_F_at_three_k\n";
    out << std::setprecision(12);

    bool wroteOrigin = false;
    auto writeSample = [&](const Point3& p, const std::array<bool, 3>& nearFlags,
                           const char* kind, const std::array<Eval, 3>& e) {
        const double fInterp = weightNegative * e[0].f + weightPositive * e[2].f;
        const double gxInterp = weightNegative * e[0].gx + weightPositive * e[2].gx;
        const double gyInterp = weightNegative * e[0].gy + weightPositive * e[2].gy;
        const double gzInterp = weightNegative * e[0].gz + weightPositive * e[2].gz;
        const double gradInterp = std::sqrt(gxInterp * gxInterp + gyInterp * gyInterp + gzInterp * gzInterp);
        const double fError = fInterp - e[1].f;
        const double gxErr = gxInterp - e[1].gx, gyErr = gyInterp - e[1].gy, gzErr = gzInterp - e[1].gz;
        const double gradError = std::sqrt(gxErr * gxErr + gyErr * gyErr + gzErr * gzErr);
        const double minAbsF = std::min({std::fabs(e[0].f), std::fabs(e[1].f), std::fabs(e[2].f)});
        out << csvQuote(surface.name) << ',' << csvQuote(surface.expression) << ','
            << kNegative << ",0," << kPositive << ',' << kind << ','
            << (nearFlags[0] ? 1 : 0) << ',' << (nearFlags[1] ? 1 : 0) << ','
            << (nearFlags[2] ? 1 : 0) << ',' << p.x << ',' << p.y << ',' << p.z << ','
            << e[0].f << ',' << e[1].f << ',' << e[2].f << ',' << fInterp << ',' << fError << ','
            << e[0].gx << ',' << e[0].gy << ',' << e[0].gz << ',' << e[0].grad << ','
            << e[1].gx << ',' << e[1].gy << ',' << e[1].gz << ',' << e[1].grad << ','
            << e[2].gx << ',' << e[2].gy << ',' << e[2].gz << ',' << e[2].grad << ','
            << gxInterp << ',' << gyInterp << ',' << gzInterp << ',' << gradInterp << ','
            << gradError << ',' << minAbsF << '\n';
        ++stats.sample_count;
        stats.max_abs_f_error = std::max(stats.max_abs_f_error, std::fabs(fError));
        stats.max_gradient_error = std::max(stats.max_gradient_error, gradError);
    };

    for (const auto& entry : nodes) {
        const LocalGridNode& node = entry.second;
        const double d = std::sqrt(node.p.x * node.p.x + node.p.y * node.p.y + node.p.z * node.p.z);
        if (d > kPreviewRadius + 1.8 * step) continue;
        if (!(node.near_surface[0] || node.near_surface[1] || node.near_surface[2])) continue;
        const bool isOrigin = d <= 1e-12;
        writeSample(node.p, node.near_surface, isOrigin ? "origin_reference" : "grid_vertex", node.values);
        if (isOrigin) wroteOrigin = true;
    }
    // Always record the analytic singular point, even if the grid resolution
    // is odd and the origin is not itself a grid vertex or neighboring cell.
    if (!wroteOrigin) {
        const Point3 origin{0.0, 0.0, 0.0};
        std::array<Eval, 3> e;
        for (int ki = 0; ki < 3; ++ki) e[static_cast<std::size_t>(ki)] = evaluate(ast, origin, ks[static_cast<std::size_t>(ki)]);
        writeSample(origin, {{false, false, false}}, "origin_reference", e);
    }
    return stats;
}

bool writeLocalInterpolationTripletPly(const std::filesystem::path& path,
                                       const SurfaceSpec& surface,
                                       const CaseMetrics* negative,
                                       const CaseMetrics* positive)
{
    if (!negative || !positive) return false;
    struct PlyVertex { float x, y, z, nx, ny, nz, grad; unsigned char r, g, b; };
    std::vector<PlyVertex> vertices;
    std::vector<std::array<std::uint32_t, 3>> faces;
    auto appendMesh = [&](double k, unsigned char r, unsigned char g, unsigned char b) {
        NodePtr ast = parseExpression(surface.expression);
        IsosurfaceResult mesh = Isosurface::build(ast.get(), k, gResolution, kHalfExtent);
        for (std::size_t i = 0; i + 2 < mesh.vertices.size(); i += 3) {
            const SurfaceVertex tri[3] = {mesh.vertices[i], mesh.vertices[i + 1], mesh.vertices[i + 2]};
            bool local = true;
            for (const SurfaceVertex& v : tri) {
                const double d = std::sqrt(static_cast<double>(v.px) * v.px +
                                           static_cast<double>(v.py) * v.py +
                                           static_cast<double>(v.pz) * v.pz);
                if (d > kPreviewRadius) { local = false; break; }
            }
            if (!local) continue;
            const std::uint32_t base = static_cast<std::uint32_t>(vertices.size());
            for (const SurfaceVertex& v : tri)
                vertices.push_back({v.px, v.py, v.pz, v.nx, v.ny, v.nz, v.gradMag, r, g, b});
            faces.push_back({base, static_cast<std::uint32_t>(base + 1),
                             static_cast<std::uint32_t>(base + 2)});
        }
    };
    appendMesh(negative->k, 55, 150, 255);  // negative k: blue
    appendMesh(0.0, 205, 185, 245);         // singular target member: pale violet
    appendMesh(positive->k, 255, 105, 55);  // positive k: orange-red
    if (faces.empty()) return false;
    std::ofstream out(path);
    if (!out) return false;
    out << "ply\nformat ascii 1.0\ncomment " << surface.name
        << " local interpolation triplet; blue=k-, violet=k=0, orange=k+\n"
        << "element vertex " << vertices.size() << "\n"
        << "property float x\nproperty float y\nproperty float z\n"
        << "property float nx\nproperty float ny\nproperty float nz\n"
        << "property float grad_mag\n"
        << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
        << "element face " << faces.size() << "\n"
        << "property list uchar uint vertex_indices\nend_header\n";
    out << std::setprecision(8);
    for (const auto& v : vertices)
        out << v.x << ' ' << v.y << ' ' << v.z << ' '
            << v.nx << ' ' << v.ny << ' ' << v.nz << ' ' << v.grad << ' '
            << static_cast<int>(v.r) << ' ' << static_cast<int>(v.g) << ' '
            << static_cast<int>(v.b) << '\n';
    for (const auto& f : faces) out << "3 " << f[0] << ' ' << f[1] << ' ' << f[2] << '\n';
    return true;
}

void writeInterpolationValidation(const std::filesystem::path& path,
                                 const std::map<std::string, LocalInterpolationStats>& statsBySlug)
{
    std::ofstream out(path);
    out << "Pointwise interpolation check using common spatial sample points\n"
        << "F_interp(x,0) = w_minus*F(x,k_minus) + w_plus*F(x,k_plus), where weights place the interpolation at k=0.\n"
        << "The same weighted interpolation is applied to the analytic gradient components.\n"
        << "This is a scalar-field interpolation check; it does not claim that independently meshed vertices correspond one-to-one.\n\n";
    out << std::setprecision(12);
    for (const auto& entry : statsBySlug) {
        out << entry.first << ": samples=" << entry.second.sample_count
            << ", max |F_interpolated(0)-F_direct(0)|=" << entry.second.max_abs_f_error
            << ", max gradient-vector interpolation error=" << entry.second.max_gradient_error << '\n';
    }
}

void writeSummary(const std::filesystem::path& path,
                  const std::vector<CaseMetrics>& rows)
{
    std::ofstream out(path);
    out << "k sweep for the 17 built-in surfaces\n"
        << "k values: -0.0001, -0.001, -0.01, -0.1, 0, 0.1, 0.01, 0.001, 0.0001\n"
        << "A candidate must: have a nonzero k, have no known analytic singular point on F=0, "
        << "have sampled nearby vertices with nonzero gradient, and have no detected local/global "
        << "topology defects (domain-crop boundary edges are allowed).\n\n";
    for (const auto& s : kSurfaces) {
        const CaseMetrics* neg = chooseBest(rows, s, false);
        const CaseMetrics* pos = chooseBest(rows, s, true);
        out << s.name << "\n  expression: " << s.expression
            << "\n  singularity: " << s.singularity_description << '\n';
        out << "  best negative k: ";
        if (neg) out << std::setprecision(8) << neg->k << " (min |grad F| near origin="
                     << neg->near_grad_min << ", triangles=" << neg->topology.triangles << ')';
        else out << "NONE (no negative candidate met all criteria)";
        out << "\n  best positive k: ";
        if (pos) out << std::setprecision(8) << pos->k << " (min |grad F| near origin="
                     << pos->near_grad_min << ", triangles=" << pos->topology.triangles << ')';
        else out << "NONE (no positive candidate met all criteria)";
        out << '\n';
        if (neg && pos) out << "  pair preview: " << s.slug << "_winning_pair.ply\n"
                    << "  interpolation triplet: " << s.slug << "_interpolation_triplet.ply\n"
                    << "  near-surface samples: " << s.slug << "_near_surface_samples.csv\n";
        else out << "  preview: no two-sided pair; see CSV for the exact rejection reason.\n";
        out << '\n';
    }
}

} // namespace

int main(int argc, char** argv)
{
    namespace fs = std::filesystem;
    if (argc > 1) {
        try {
            gResolution = std::stoi(argv[1]);
        } catch (...) {
            std::cerr << "Usage: " << argv[0] << " [grid-resolution 8..128]\n";
            return 2;
        }
        if (gResolution < 8 || gResolution > 128) {
            std::cerr << "Grid resolution must be between 8 and 128.\n";
            return 2;
        }
    }
    const fs::path outputDir("singularity_sweep_artifacts");
    std::error_code ec;
    fs::create_directories(outputDir, ec);
    if (ec) {
        std::cerr << "ERROR: cannot create " << outputDir << ": " << ec.message() << '\n';
        return 2;
    }

    std::vector<CaseMetrics> rows;
    rows.reserve(kSurfaces.size() * kValues.size());
    int failures = 0;

    std::cout << "Singularity k sweep: " << kSurfaces.size() << " surface families x "
              << kValues.size() << " k values; grid resolution=" << gResolution
              << ", domain=[-" << kHalfExtent << ',' << kHalfExtent << "]^3\n";
    for (const SurfaceSpec& surface : kSurfaces) {
        NodePtr ast;
        try {
            ast = parseExpression(surface.expression);
        } catch (const std::exception& e) {
            std::cerr << "FAIL: cannot parse " << surface.name << ": " << e.what() << '\n';
            ++failures;
            continue;
        }
        std::cout << "\n[" << surface.name << "] " << surface.expression << '\n';
        for (double k : kValues) {
            try {
                CaseMetrics row = analyzeCase(surface, ast.get(), k);
                if (row.vertex_count % 3 != 0) {
                    std::cerr << "FAIL: " << surface.name << " k=" << k
                              << " emitted a vertex count not divisible by 3\n";
                    ++failures;
                }
                if (!row.finite_vertices) {
                    std::cerr << "FAIL: non-finite/invalid mesh vertex for "
                              << surface.name << " k=" << k << '\n';
                    ++failures;
                }
                if (!finiteEval(row.origin)) {
                    std::cerr << "FAIL: non-finite value/gradient at origin for "
                              << surface.name << " k=" << k << '\n';
                    ++failures;
                }
                if (std::fabs(k) <= 1e-15 &&
                    (std::fabs(row.origin.f) > 1e-9 || row.origin.grad > 1e-8)) {
                    std::cerr << "FAIL: the k=0 preset is expected to have F(0)=0 and grad F(0)=0: "
                              << surface.name << '\n';
                    ++failures;
                }
                if (surface.behavior == Behavior::WhitneyUmbrella && k > 0.0) {
                    const double sx = std::sqrt(k);
                    for (double sign : {-1.0, 1.0}) {
                        const Eval shifted = evaluate(ast.get(), {sign * sx, 0.0, 0.0}, k);
                        if (std::fabs(shifted.f) > 1e-8 || shifted.grad > 1e-8) {
                            std::cerr << "FAIL: expected Whitney singular point at x=" << sign * sx
                                      << " for k=" << k << '\n';
                            ++failures;
                        }
                    }
                }
                std::cout << "  k=" << std::setw(9) << std::setprecision(6) << k
                          << "  tris=" << std::setw(6) << row.topology.triangles
                          << "  |grad(0)|=" << std::setw(10) << row.origin.grad
                          << "  nearMin|grad|=";
                if (std::isfinite(row.near_grad_min)) std::cout << row.near_grad_min;
                else std::cout << "n/a";
                std::cout << "  topology=" << (row.topology.clean ? "clean" : "DEFECT")
                          << "  candidate=" << (row.candidate_ok ? "yes" : "no") << '\n';
                rows.push_back(std::move(row));
            } catch (const std::exception& e) {
                std::cerr << "FAIL: sweep crashed for " << surface.name << " k=" << k
                          << ": " << e.what() << '\n';
                ++failures;
            }
        }
    }

    writeCaseCsv(outputDir / "singularity_sweep.csv", rows);
    writeSingularityPointCsv(outputDir / "singularity_points.csv");
    writeSummary(outputDir / "singularity_sweep_summary.txt", rows);

    std::cout << "\n=== Winning +/- k pairs ===\n";
    std::map<std::string, LocalInterpolationStats> interpolationStats;
    for (const SurfaceSpec& surface : kSurfaces) {
        const CaseMetrics* neg = chooseBest(rows, surface, false);
        const CaseMetrics* pos = chooseBest(rows, surface, true);
        std::cout << surface.name << ": negative=";
        if (neg) std::cout << neg->k; else std::cout << "NONE";
        std::cout << ", positive=";
        if (pos) std::cout << pos->k; else std::cout << "NONE";
        std::cout << '\n';
        if (neg || pos) {
            const fs::path preview = outputDir / (std::string(surface.slug) + "_winning_pair.ply");
            if (writeLocalPairPly(preview, surface, neg, pos))
                std::cout << "  local colored PLY: " << preview.string() << '\n';
            else
                std::cout << "  no local triangles written to pair preview.\n";
        }
        if (neg && pos) {
            const fs::path triplet = outputDir / (std::string(surface.slug) + "_interpolation_triplet.ply");
            if (writeLocalInterpolationTripletPly(triplet, surface, neg, pos))
                std::cout << "  negative / singular k=0 / positive triplet: " << triplet.string() << '\n';
            NodePtr ast = parseExpression(surface.expression);
            const fs::path samples = outputDir / (std::string(surface.slug) + "_near_surface_samples.csv");
            const LocalInterpolationStats stats = writeNearSurfaceSamples(samples, surface, ast.get(), neg->k, pos->k);
            interpolationStats.emplace(surface.slug, stats);
            std::cout << "  nearby shared-grid samples: " << samples.string()
                      << " (" << stats.sample_count << " points; max |F interpolation error|="
                      << std::scientific << stats.max_abs_f_error << std::defaultfloat << ")\n";
        }
    }
    writeInterpolationValidation(outputDir / "interpolation_validation.txt", interpolationStats);

    std::cout << "\nArtifacts written to: " << outputDir.string() << '\n'
              << "  singularity_sweep.csv: all 153 mesh/k cases and topology/gradient metrics\n"
              << "  singularity_points.csv: exact origin and moving Whitney singularity probes\n"
              << "  singularity_sweep_summary.txt: selected smallest eligible negative/positive k\n"
              << "  *_winning_pair.ply: two regularized surfaces, blue=-k and orange=+k\n"
              << "  *_interpolation_triplet.ply: k-, actual singular target k=0, and k+ surfaces in one local view\n"
              << "  *_near_surface_samples.csv: shared grid vertices in cells intersected by any of the three isosurfaces; F and gradients at all three k values\n"
              << "  interpolation_validation.txt: pointwise field / gradient interpolation errors\n";

    if (failures == 0) {
        std::cout << "\nSWEEP TEST COMPLETED (results are diagnostics; 'no winning pair' is a valid finding).\n";
        return 0;
    }
    std::cout << "\n" << failures << " TEST FAILURE(S)\n";
    return 1;
}
