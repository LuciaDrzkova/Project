#pragma once

// Reconstructions of the singular surface S_0 = {F = 0} near an isolated
// singular point P, from the perturbed surfaces S_-eps and S_+eps of the family
// H_k = F + kG, or from F itself. Shared by the viewer ("Reconstruct S_0") and
// by tests/interpolation_experiments.cpp.
//
// Linear interpolation (the thesis report, section "Reconstructing S_0 by
// linear interpolation"):
//   FieldAverage    Marching Cubes of (H_-eps + H_+eps)/2 on the grid. H_k is
//                   affine in k, so this IS the grid of F (up to rounding) and
//                   the mesh is the direct Marching Cubes mesh of F.
//   VertexMid*      every vertex v of mesh(S_-eps) (or S_+eps) moved to the
//                   midpoint of v and its closest point on the other mesh.
//   ExactVertexMid  the same with closest points on the exact surfaces.
//   SdfAverage      zero level of the average of the two signed distance fields.
//
// Beyond linear interpolation:
//   NormalMid       the correspondence follows the gradient line of F instead
//                   of the closest point: from a on S_-eps along
//                   n = grad F(a)/|grad F(a)| to S_+eps. This is the expansion
//                   s(k) = -k G/|grad F| + c k^2 of the report, so the midpoint
//                   cancels the first-order term by construction.
//   SingularCone    the proposed singularity-aware extraction. Outside a small
//                   box B around P (patchRadius): Marching Cubes of F. Inside B: the
//                   surface is rebuilt from its own boundary loop by the
//                   weighted scaling Phi_t(p) = (t^w1 x, t^w2 y, t^w3 z) that
//                   maps the principal part F_0 of F into itself
//                   (F_0(Phi_t p) = t F_0(p)), each new vertex projected onto
//                   {F = 0} by Newton's method, and the last ring is closed at
//                   P. For a nondegenerate weighted-homogeneous principal part
//                   every orbit of Phi_t enters B exactly once, so S_0 n B is
//                   the union of the orbits through its boundary curve: the
//                   patch has the topology of the cone over the link (the conic
//                   structure theorem) and contains P as a vertex.
//   HybridCone      the same patch built from the vertex interpolation of
//                   S_-eps and S_+eps instead of from F, without projection:
//                   the interpolation is used only away from P, where it is
//                   second order in eps, and the part within the inner scale
//                   lambda(eps) is replaced by the scaled boundary loop.

#include "ConvergenceClassifier.h"
#include "ExactDistance.h"
#include "SingularityResearch.h"
#include "core/Geometry.h"
#include "core/Grid.h"
#include "core/Isosurface.h"
#include "core/Parallel.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <functional>
#include <limits>
#include <map>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace iso_reconstruct {

using iso_geometry::Point3;
constexpr double kInf = std::numeric_limits<double>::infinity();

// ----- indexed meshes in double precision -------------------------------------

struct Mesh {
    std::vector<Point3> v;
    std::vector<std::array<std::uint32_t, 3>> t;
};

inline double infDistance(const Point3& p, const Point3& c) {
    return std::max({std::fabs(p.x - c.x), std::fabs(p.y - c.y), std::fabs(p.z - c.z)});
}

// Welds corners at identical positions. The extractor's own vertex ids are not
// used: when the field is (numerically) zero at a grid node, several edge
// crossings fall onto that node, the extractor drops the resulting zero-area
// triangles, and with per-edge ids this would leave artificial cracks.
// Triangles that collapse to repeated corners after welding are dropped.
inline Mesh fromResult(const IsosurfaceResult& r, double cell) {
    Mesh m;
    const double tol = 2e-4 * cell; // the extractor drops triangles whose corners are this close
    std::map<std::array<long long, 3>, std::uint32_t> ids;
    for (std::size_t i = 0; i + 2 < r.vertices.size(); i += 3) {
        std::array<std::uint32_t, 3> tri{};
        for (int c = 0; c < 3; ++c) {
            const Point3 p = iso_research::surfacePoint(r.vertices[i + static_cast<std::size_t>(c)]);
            const std::array<long long, 3> key{std::llround(p.x / tol), std::llround(p.y / tol), std::llround(p.z / tol)};
            const auto [it, inserted] = ids.emplace(key, static_cast<std::uint32_t>(m.v.size()));
            if (inserted) m.v.push_back(p);
            tri[static_cast<std::size_t>(c)] = it->second;
        }
        if (tri[0] != tri[1] && tri[1] != tri[2] && tri[0] != tri[2]) m.t.push_back(tri);
    }
    return m;
}

inline iso_geometry::TriangleBvh bvh(const Mesh& m) {
    std::vector<iso_geometry::Triangle> tris;
    tris.reserve(m.t.size());
    for (const auto& t : m.t) tris.push_back(iso_geometry::makeTriangle(m.v[t[0]], m.v[t[1]], m.v[t[2]]));
    return iso_geometry::TriangleBvh(std::move(tris));
}

// Triangle soup with shared vertex ids, normals and |grad F| from `field`, for
// the viewer's display code.
inline IsosurfaceResult toResult(const Mesh& m, const iso_exact::Field& field) {
    IsosurfaceResult r;
    std::vector<SurfaceVertex> vertex(m.v.size());
    iso_parallel::parallelFor(m.v.size(), [&](std::size_t i) {
        const auto d = field.evaluate(m.v[i]);
        const double g = std::sqrt(d.gradient[0] * d.gradient[0] + d.gradient[1] * d.gradient[1] + d.gradient[2] * d.gradient[2]);
        SurfaceVertex& s = vertex[i];
        s.px = static_cast<float>(m.v[i].x);
        s.py = static_cast<float>(m.v[i].y);
        s.pz = static_cast<float>(m.v[i].z);
        s.gradMag = static_cast<float>(g);
        const double inv = g > 0.0 ? 1.0 / g : 0.0;
        s.nx = static_cast<float>(d.gradient[0] * inv);
        s.ny = static_cast<float>(d.gradient[1] * inv);
        s.nz = static_cast<float>(d.gradient[2] * inv);
    });
    r.vertices.reserve(m.t.size() * 3);
    r.vertexIds.reserve(m.t.size() * 3);
    for (const auto& t : m.t)
        for (std::uint32_t id : t) {
            r.vertices.push_back(vertex[id]);
            r.vertexIds.push_back(id);
            r.maxGradMag = std::max(r.maxGradMag, static_cast<double>(vertex[id].gradMag));
        }
    r.uniqueVertexCount = m.v.size();
    return r;
}

// ----- topology inside a box ---------------------------------------------------

struct Topology {
    std::size_t triangles = 0, vertices = 0, edges = 0, components = 0;
    long long euler = 0;
    std::size_t holeEdges = 0;        // boundary edges well inside the box (not created by clipping)
    std::size_t nonManifoldEdges = 0; // edges with more than two triangles
    std::size_t degenerate = 0;       // triangles of (almost) zero area
    std::size_t flipped = static_cast<std::size_t>(-1); // orientation reversed w.r.t. `source` (same triangles)

    // Near an isolated singular point S_0 is the cone over its link
    // (conic structure theorem), which is contractible: one component, Euler
    // characteristic 1, no holes.
    bool contractible() const { return components == 1 && euler == 1 && holeEdges == 0 && nonManifoldEdges == 0; }
};

// Triangles with centroid in center + [-box, box]^3.
inline Topology topology(const Mesh& m, const Point3& center, double box, double h, const Mesh* source = nullptr) {
    auto cross = [](const Point3& a, const Point3& b) {
        return Point3{a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x};
    };
    auto normal = [&](const Mesh& mesh, const std::array<std::uint32_t, 3>& t) {
        return cross(mesh.v[t[1]] - mesh.v[t[0]], mesh.v[t[2]] - mesh.v[t[0]]);
    };
    Topology tp;
    std::vector<std::size_t> kept;
    for (std::size_t i = 0; i < m.t.size(); ++i) {
        const auto& t = m.t[i];
        if (infDistance((m.v[t[0]] + m.v[t[1]] + m.v[t[2]]) * (1.0 / 3.0), center) <= box) kept.push_back(i);
    }
    tp.triangles = kept.size();
    std::vector<std::uint32_t> parent(m.v.size());
    for (std::uint32_t i = 0; i < parent.size(); ++i) parent[i] = i;
    auto find = [&](std::uint32_t a) {
        while (parent[a] != a) a = parent[a] = parent[parent[a]];
        return a;
    };
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> edges;
    std::vector<char> used(m.v.size(), 0);
    if (source) tp.flipped = 0;
    for (std::size_t i : kept) {
        const auto& t = m.t[i];
        for (int c = 0; c < 3; ++c) {
            const std::uint32_t a = t[static_cast<std::size_t>(c)], b = t[static_cast<std::size_t>((c + 1) % 3)];
            used[a] = 1;
            ++edges[{std::min(a, b), std::max(a, b)}];
            parent[find(a)] = find(b);
        }
        const Point3 n = normal(m, t);
        if (iso_geometry::norm(n) < 1e-12 * h * h) ++tp.degenerate;
        if (source && i < source->t.size() && iso_geometry::dot(n, normal(*source, source->t[i])) < 0.0) ++tp.flipped;
    }
    std::vector<char> root(m.v.size(), 0);
    for (std::size_t v = 0; v < m.v.size(); ++v)
        if (used[v]) {
            ++tp.vertices;
            const std::uint32_t r = find(static_cast<std::uint32_t>(v));
            if (!root[r]) { root[r] = 1; ++tp.components; }
        }
    tp.edges = edges.size();
    for (const auto& [e, count] : edges) {
        if (count > 2) ++tp.nonManifoldEdges;
        if (count == 1 && infDistance(m.v[e.first], center) < box - 2.0 * h && infDistance(m.v[e.second], center) < box - 2.0 * h)
            ++tp.holeEdges;
    }
    tp.euler = static_cast<long long>(tp.vertices) - static_cast<long long>(tp.edges) + static_cast<long long>(tp.triangles);
    return tp;
}

// ----- methods ----------------------------------------------------------------

enum class Method : int {
    FieldAverage,
    VertexMidMinus,
    VertexMidPlus,
    ExactVertexMid,
    NormalMid,
    SdfAverage,
    HybridCone,
    SingularCone,
    Count
};

inline const char* methodId(Method m) {
    switch (m) {
    case Method::FieldAverage: return "field_average";
    case Method::VertexMidMinus: return "vertex_mid_from_minus";
    case Method::VertexMidPlus: return "vertex_mid_from_plus";
    case Method::ExactVertexMid: return "exact_vertex_mid";
    case Method::NormalMid: return "normal_mid";
    case Method::SdfAverage: return "sdf_average";
    case Method::HybridCone: return "hybrid_cone";
    case Method::SingularCone: return "singular_cone";
    default: return "unknown";
    }
}

inline const char* methodLabel(Method m) {
    switch (m) {
    case Method::FieldAverage: return "Average of implicit values (= MC of F)";
    case Method::VertexMidMinus: return "Vertex interpolation, mesh (from S-eps)";
    case Method::VertexMidPlus: return "Vertex interpolation, mesh (from S+eps)";
    case Method::ExactVertexMid: return "Vertex interpolation, exact surfaces";
    case Method::NormalMid: return "Vertex interpolation along grad F";
    case Method::SdfAverage: return "Average of signed distance fields";
    case Method::HybridCone: return "Hybrid: interpolation + cone patch at P";
    case Method::SingularCone: return "Singular cone patch (MC of F + cone at P)";
    default: return "unknown";
    }
}

// Closest point on an exact surface (or nothing).
using ClosestFn = std::function<std::optional<Point3>(const Point3&)>;

// Only vertices within `halfSize` (max norm) of `center` are moved; the others
// are irrelevant for local measurements and are left in place.
struct Region {
    Point3 center{0, 0, 0};
    double halfSize = kInf;
    bool contains(const Point3& p) const { return infDistance(p, center) <= halfSize; }
};

// Closest point of {H = 0} to a, started from the closest point of the mesh
// of that surface: moved onto the surface, then refined by the Lagrange
// conditions of min |p - a|^2.
inline ClosestFn closestOnField(const iso_exact::Field& field, const iso_geometry::TriangleBvh& mesh, double maxStep) {
    return [&field, &mesh, maxStep](const Point3& a) -> std::optional<Point3> {
        Point3 start;
        if (!mesh.nearestPoint(a, start) || !iso_exact::detail::project(field, start, maxStep)) return std::nullopt;
        Point3 best = start;
        if (const auto p = iso_exact::detail::lagrange(field, a, start, maxStep))
            if (iso_geometry::distance(a, *p) < iso_geometry::distance(a, best)) best = *p;
        return best;
    };
}

// Vertices of `from` moved halfway to their closest point on the mesh `to`.
inline Mesh vertexMid(const Mesh& from, const iso_geometry::TriangleBvh& to, const Region& region = {}) {
    Mesh m = from;
    iso_parallel::parallelFor(m.v.size(), [&](std::size_t i) {
        if (!region.contains(from.v[i])) return;
        Point3 c;
        if (to.nearestPoint(from.v[i], c)) m.v[i] = (from.v[i] + c) * 0.5;
    });
    return m;
}

// The same on the exact surfaces: v is first moved to its closest point a on
// S_-eps, then to the midpoint of a and the closest point of a on S_+eps.
inline Mesh exactVertexMid(const Mesh& from, const ClosestFn& closestMinus, const ClosestFn& closestPlus,
                           const Region& region, std::size_t& failures) {
    Mesh m = from;
    std::vector<char> failed(m.v.size(), 0);
    iso_parallel::parallelFor(m.v.size(), [&](std::size_t i) {
        if (!region.contains(from.v[i])) return;
        const auto a = closestMinus(from.v[i]);
        if (!a) { failed[i] = 1; return; }
        const auto b = closestPlus(*a);
        if (!b) { failed[i] = 1; return; }
        m.v[i] = (*a + *b) * 0.5;
    });
    failures = static_cast<std::size_t>(std::count(failed.begin(), failed.end(), 1));
    return m;
}

// Correspondence along the gradient line of F. Each vertex is moved onto
// S_-eps (point a), then the line a + s n, n = grad F(a)/|grad F(a)|, is
// followed to S_+eps by Newton's method in s, and the vertex goes to
// a + (s/2) n. Where grad F vanishes or the line does not reach S_+eps near a
// (|s| more than three times the closest-point distance plus two cells), the
// closest point is used instead; `fallbacks` counts these vertices.
inline Mesh normalMid(const Mesh& from, const iso_exact::Field& fieldMinus, const iso_exact::Field& fieldF,
                      const iso_exact::Field& fieldPlus, const ClosestFn& closestPlus, double h, const Region& region,
                      std::size_t& fallbacks, std::size_t& failures) {
    Mesh m = from;
    std::vector<char> fellBack(m.v.size(), 0), failed(m.v.size(), 0);
    iso_parallel::parallelFor(m.v.size(), [&](std::size_t i) {
        if (!region.contains(from.v[i])) return;
        Point3 a = from.v[i];
        if (!iso_exact::detail::project(fieldMinus, a, h)) { failed[i] = 1; return; }
        const auto c = closestPlus(a);
        const double closestDistance = c ? iso_geometry::distance(a, *c) : kInf;
        const auto df = fieldF.evaluate(a);
        Point3 n{df.gradient[0], df.gradient[1], df.gradient[2]};
        const double gn = iso_geometry::norm(n);
        bool ok = gn > 0.0 && std::isfinite(gn);
        double s = 0.0;
        if (ok) {
            n = n * (1.0 / gn);
            ok = false;
            for (int it = 0; it < 60; ++it) {
                const auto d = fieldPlus.evaluate(a + n * s);
                const double slope = d.gradient[0] * n.x + d.gradient[1] * n.y + d.gradient[2] * n.z;
                if (!(std::fabs(slope) > 0.0) || !std::isfinite(d.value)) break;
                double ds = d.value / slope;
                if (std::fabs(ds) > h) ds = std::copysign(h, ds);
                s -= ds;
                if (std::fabs(ds) <= 1e-14 * (1.0 + std::fabs(s))) { ok = true; break; }
            }
            ok = ok && std::fabs(s) <= 3.0 * closestDistance + 2.0 * h;
        }
        if (ok) {
            m.v[i] = a + n * (0.5 * s);
        } else if (c) {
            m.v[i] = (a + *c) * 0.5;
            fellBack[i] = 1;
        } else {
            failed[i] = 1;
        }
    });
    fallbacks = static_cast<std::size_t>(std::count(fellBack.begin(), fellBack.end(), 1));
    failures = static_cast<std::size_t>(std::count(failed.begin(), failed.end(), 1));
    return m;
}

// Zero level of (phi_- + phi_+)/2, phi_+- the signed distance to mesh(S_+-eps)
// with the sign of the grid samples of H_+-eps, extracted by Marching Cubes on
// the same grid. Distances are computed only at nodes that touch a grid edge on
// which the average can change sign; at every other node both signs agree with
// all six neighbours, so only the sign matters there.
inline Mesh sdfAverage(const Node* ast, int res, double halfExtent, const std::vector<double>& hMinus,
                       const std::vector<double>& hPlus, const iso_geometry::TriangleBvh& minus,
                       const iso_geometry::TriangleBvh& plus) {
    const int n = res + 1;
    const std::size_t nn = static_cast<std::size_t>(n);
    auto index = [nn](int i, int j, int l) {
        return (static_cast<std::size_t>(i) * nn + static_cast<std::size_t>(j)) * nn + static_cast<std::size_t>(l);
    };
    auto sign = [](double v) { return v < 0.0 ? -1 : 1; };
    std::vector<char> need(hMinus.size(), 0);
    iso_parallel::parallelFor(nn, [&](std::size_t si) {
        const int i = static_cast<int>(si);
        for (int j = 0; j < n; ++j)
            for (int l = 0; l < n; ++l) {
                const std::size_t id = index(i, j, l);
                const int sm = sign(hMinus[id]), sp = sign(hPlus[id]);
                bool flag = sm != sp;
                const int nb[6][3] = {{i - 1, j, l}, {i + 1, j, l}, {i, j - 1, l}, {i, j + 1, l}, {i, j, l - 1}, {i, j, l + 1}};
                for (const auto& q : nb) {
                    if (flag) break;
                    if (q[0] < 0 || q[1] < 0 || q[2] < 0 || q[0] >= n || q[1] >= n || q[2] >= n) continue;
                    const std::size_t qd = index(q[0], q[1], q[2]);
                    flag = sign(hMinus[qd]) != sm || sign(hPlus[qd]) != sp || sign(hMinus[qd]) != sign(hPlus[qd]);
                }
                need[id] = flag;
            }
    });
    std::vector<double> phi(hMinus.size());
    iso_parallel::parallelFor(nn, [&](std::size_t si) {
        const int i = static_cast<int>(si);
        for (int j = 0; j < n; ++j)
            for (int l = 0; l < n; ++l) {
                const std::size_t id = index(i, j, l);
                if (!need[id]) {
                    phi[id] = hMinus[id] < 0.0 ? -1.0 : 1.0;
                    continue;
                }
                const Point3 p{gridCoordinate(i, res, halfExtent), gridCoordinate(j, res, halfExtent),
                               gridCoordinate(l, res, halfExtent)};
                const double dm = (hMinus[id] < 0.0 ? -1.0 : 1.0) * minus.nearestDistance(p);
                const double dp = (hPlus[id] < 0.0 ? -1.0 : 1.0) * plus.nearestDistance(p);
                phi[id] = 0.5 * (dm + dp);
            }
    });
    return fromResult(Isosurface::extract(ast, 0.0, res, halfExtent, std::move(phi)), 2.0 * halfExtent / res);
}

// ----- local structure at P ---------------------------------------------------

struct LocalStructure {
    bool singular = false;          // F(P) = 0 and grad F(P) = 0 (relative to the 3-jet)
    bool weightsFound = false;
    bool nondegenerate = false;     // the principal part has an isolated critical point
    std::array<double, 3> weights{iso_exact::kNaN, iso_exact::kNaN, iso_exact::kNaN};
    std::string principalPart;
    double weightedOrderG = kInf;   // weighted order of G at P
    double innerExponent = 1.0;     // the deformation has size ~ eps^innerExponent near P
    std::string note;
};

// For the family H_k = F + kG (affine in k): type data of F at P and the
// weighted order of G = H_1 - H_0.
inline LocalStructure localStructure(const Node* familyAst, const Point3& P) {
    LocalStructure s;
    const auto f = iso_research::evaluateJet(familyAst, P, 0.0).taylor;
    const auto f1 = iso_research::evaluateJet(familyAst, P, 1.0).taylor;
    std::array<double, 20> g{};
    for (std::size_t i = 0; i < 20; ++i) g[i] = f1[i] - f[i];
    const double scale = std::max(iso_classify::detail::jetScale(f), 1e-300);
    s.singular = std::fabs(f[0]) <= 1e-9 * scale && std::fabs(f[1]) <= 1e-9 * scale && std::fabs(f[2]) <= 1e-9 * scale &&
                 std::fabs(f[3]) <= 1e-9 * scale;
    if (!s.singular) {
        s.note = "P is not a singular point of S_0 (F or grad F is nonzero at P)";
        return s;
    }
    const auto pp = iso_classify::detail::principalPart(f);
    s.weightsFound = pp.found;
    s.nondegenerate = pp.nondegenerate;
    s.weights = pp.weights;
    s.principalPart = pp.text;
    if (!pp.found) {
        s.note = "no weighted principal part within the 3-jet: the cone patch is not applicable";
        return s;
    }
    if (!pp.nondegenerate) s.note = "principal part numerically degenerate: the cone patch may be wrong";
    const auto& e = iso_research::exponents();
    const double sg = std::max(iso_classify::detail::jetScale(g), 1e-300);
    for (std::size_t i = 0; i < 20; ++i)
        if (std::fabs(g[i]) > iso_classify::detail::kZero * sg)
            s.weightedOrderG = std::min(s.weightedOrderG, iso_classify::detail::weightedDegree(e[i], s.weights));
    const double wMin = std::min({s.weights[0], s.weights[1], s.weights[2]});
    if (s.weightedOrderG < 1.0 - 1e-9) s.innerExponent = std::min(1.0, wMin / (1.0 - s.weightedOrderG));
    return s;
}

// ----- the weighted-cone patch ------------------------------------------------

struct PatchInfo {
    bool applied = false;
    double radius = 0.0;            // half-size of the box B around P
    int loops = 0;                  // closed boundary loops on the box
    int rings = 0;                  // scaled copies per loop
    std::size_t loopVertices = 0;
    std::size_t openChains = 0;     // boundary chains near B that did not close (holes)
    std::size_t projectionFailures = 0;
    std::string note;
};

// Replaces the part of `m` inside the box P + [-radius, radius]^3 by the cone
// over its boundary loops under the weighted scaling Phi_t. With `field`, the
// loop and every new vertex are projected onto {field = 0}; without it the
// scaled loops are used as they are.
inline PatchInfo conePatch(Mesh& m, const Point3& P, const std::array<double, 3>& w, double radius, double h,
                           const iso_exact::Field* field) {
    PatchInfo info;
    info.radius = radius;
    const double wMin = std::min({w[0], w[1], w[2]});
    if (!(wMin > 0.0) || !(radius > h)) {
        info.note = "invalid weights or box";
        return info;
    }
    const auto original = m.t;
    std::vector<std::array<std::uint32_t, 3>> kept;
    kept.reserve(m.t.size());
    for (const auto& t : m.t)
        if (infDistance(m.v[t[0]], P) >= radius && infDistance(m.v[t[1]], P) >= radius && infDistance(m.v[t[2]], P) >= radius)
            kept.push_back(t);
    if (kept.size() == m.t.size()) {
        info.note = "the mesh does not enter the box around P";
        return info;
    }
    m.t = std::move(kept);

    // Directed boundary edges (a -> b as in their triangle) next to the box.
    std::map<std::pair<std::uint32_t, std::uint32_t>, int> uses;
    for (const auto& t : m.t)
        for (int c = 0; c < 3; ++c) {
            const std::uint32_t a = t[static_cast<std::size_t>(c)], b = t[static_cast<std::size_t>((c + 1) % 3)];
            ++uses[{std::min(a, b), std::max(a, b)}];
        }
    std::map<std::uint32_t, std::vector<std::uint32_t>> next;
    std::size_t boundaryEdges = 0;
    for (const auto& t : m.t)
        for (int c = 0; c < 3; ++c) {
            const std::uint32_t a = t[static_cast<std::size_t>(c)], b = t[static_cast<std::size_t>((c + 1) % 3)];
            if (uses[{std::min(a, b), std::max(a, b)}] != 1) continue;
            if (infDistance(m.v[a], P) > radius + 2.0 * h || infDistance(m.v[b], P) > radius + 2.0 * h) continue;
            next[a].push_back(b);
            ++boundaryEdges;
        }

    std::vector<std::vector<std::uint32_t>> loops;
    while (true) {
        auto it = std::find_if(next.begin(), next.end(), [](const auto& e) { return !e.second.empty(); });
        if (it == next.end()) break;
        const std::uint32_t start = it->first;
        std::vector<std::uint32_t> loop{start};
        std::uint32_t cur = start;
        bool closed = false;
        while (true) {
            auto& out = next[cur];
            if (out.empty()) break;
            const std::uint32_t b = out.back();
            out.pop_back();
            if (b == start) { closed = true; break; }
            loop.push_back(b);
            cur = b;
            if (loop.size() > boundaryEdges + 1) break;
        }
        if (closed && loop.size() >= 3) loops.push_back(std::move(loop));
        else ++info.openChains;
    }
    if (loops.empty()) {
        m.t = original;
        info.note = "no closed boundary loop around P: patch not applied";
        return info;
    }

    // Loop vertices onto the surface.
    if (field) {
        std::vector<std::uint32_t> ids;
        for (const auto& l : loops) ids.insert(ids.end(), l.begin(), l.end());
        std::sort(ids.begin(), ids.end());
        ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
        for (std::uint32_t id : ids) {
            Point3 p = m.v[id];
            if (iso_exact::detail::project(*field, p, 0.5 * h)) m.v[id] = p;
            else ++info.projectionFailures;
        }
    }

    // Rings: the slowest axis shrinks by q = 1 - h/radius per ring, so the
    // first ring lies inside the box and the rings are similar copies, with the
    // spacing of the loop's edges. They stop at a twentieth of a cell.
    const double q = std::max(0.25, 1.0 - h / radius);
    const double rho = std::pow(q, 1.0 / wMin);
    const std::uint32_t center = static_cast<std::uint32_t>(m.v.size());
    m.v.push_back(P);
    for (const auto& loop : loops) {
        const std::size_t n = loop.size();
        std::vector<Point3> base(n);
        double maxInf = 0.0;
        for (std::size_t i = 0; i < n; ++i) {
            base[i] = m.v[loop[i]] - P;
            maxInf = std::max(maxInf, infDistance(m.v[loop[i]], P));
        }
        const int rings = std::max(1, static_cast<int>(std::ceil(std::log(0.05 * h / maxInf) / std::log(q))));
        info.rings = std::max(info.rings, rings);
        std::vector<std::uint32_t> prev(loop.begin(), loop.end()), cur(n);
        for (int j = 1; j <= rings; ++j) {
            const double t = std::pow(rho, j);
            const double s[3] = {std::pow(t, w[0]), std::pow(t, w[1]), std::pow(t, w[2])};
            std::vector<Point3> ring(n);
            std::vector<char> failed(n, 0);
            iso_parallel::parallelFor(n, [&](std::size_t i) {
                Point3 p = P + Point3{s[0] * base[i].x, s[1] * base[i].y, s[2] * base[i].z};
                if (field) {
                    Point3 projected = p;
                    if (iso_exact::detail::project(*field, projected, 0.25 * infDistance(p, P) + 1e-300)) p = projected;
                    else failed[i] = 1;
                }
                ring[i] = p;
            });
            info.projectionFailures += static_cast<std::size_t>(std::count(failed.begin(), failed.end(), 1));
            for (std::size_t i = 0; i < n; ++i) {
                cur[i] = static_cast<std::uint32_t>(m.v.size());
                m.v.push_back(ring[i]);
            }
            for (std::size_t i = 0; i < n; ++i) {
                const std::size_t k = (i + 1) % n;
                m.t.push_back({prev[k], prev[i], cur[i]});
                m.t.push_back({prev[k], cur[i], cur[k]});
            }
            prev = cur;
        }
        for (std::size_t i = 0; i < n; ++i) m.t.push_back({prev[(i + 1) % n], prev[i], center});
        info.loopVertices += n;
    }
    info.loops = static_cast<int>(loops.size());
    info.applied = true;
    if (info.openChains > 0) info.note = std::to_string(info.openChains) + " boundary chain(s) did not close";
    return info;
}

// Box half-size of the cone patch: at least `cells` cells, and large enough
// that the surface is at least two cells thick on the boundary of the box. In
// the box of size t^{w_min} (the slowest coordinate) the thinnest direction of
// S_0 has the size t^{w_max}; it is two cells when t^{w_max} = 2h, i.e. for the
// half-size (2h)^{w_min/w_max}. Inside, Marching Cubes would miss parts of the
// surface (the tip of a horn cusp is thinner than a cell). For equal weights
// (the cubic cones, A_1) this is 2h and the minimum of `cells` cells decides.
inline double patchRadius(const LocalStructure& s, double h, double cells, double maxRadius) {
    const double wMin = std::min({s.weights[0], s.weights[1], s.weights[2]});
    const double wMax = std::max({s.weights[0], s.weights[1], s.weights[2]});
    double r = cells * h;
    if (wMin > 0.0 && wMax > 0.0) r = std::max(r, std::pow(2.0 * h, wMin / wMax));
    return std::min(maxRadius, r);
}

// Box half-size of the hybrid patch: as patchRadius, and at least twice the
// size eps^innerExponent of the deformation near P (inside it the two sides
// differ and the interpolation does not converge), at most `maxRadius`.
inline double hybridRadius(const LocalStructure& s, double eps, double h, double cells, double maxRadius) {
    return std::min(maxRadius, std::max(patchRadius(s, h, cells, maxRadius), 2.0 * std::pow(eps, s.innerExponent)));
}

// ----- one reconstruction for the viewer --------------------------------------

struct Settings {
    double eps = 0.01;
    int resolution = 64;
    double halfExtent = 1.0;
    Point3 center{0, 0, 0}; // the singular point P
    double patchCells = 4.0;
};

struct Output {
    Mesh mesh, minus, plus;
    LocalStructure local;
    PatchInfo patch;
    std::size_t correspondenceFallbacks = 0, correspondenceFailures = 0;
    std::string note;
    double millis = 0.0;
};

// Builds S_-eps, S_+eps and the reconstruction of S_0 with `method`.
// `familyAst` is H(x, y, z, k), assumed affine in k.
inline Output reconstruct(const Node* familyAst, Method method, const Settings& s) {
    const auto t0 = std::chrono::steady_clock::now();
    Output out;
    const int res = std::max(s.resolution, 2);
    const double h = 2.0 * s.halfExtent / res;
    std::vector<double> hm = Isosurface::sampleField(familyAst, -s.eps, res, s.halfExtent);
    std::vector<double> hp = Isosurface::sampleField(familyAst, +s.eps, res, s.halfExtent);
    out.minus = fromResult(Isosurface::extract(familyAst, -s.eps, res, s.halfExtent, hm), h);
    out.plus = fromResult(Isosurface::extract(familyAst, +s.eps, res, s.halfExtent, hp), h);
    const auto bm = bvh(out.minus), bp = bvh(out.plus);
    const iso_exact::Field fieldMinus(familyAst, -s.eps), fieldPlus(familyAst, +s.eps), fieldF(familyAst, 0.0);
    out.local = localStructure(familyAst, s.center);

    auto fieldAverage = [&] {
        std::vector<double> avg(hm.size());
        for (std::size_t i = 0; i < avg.size(); ++i) avg[i] = 0.5 * (hm[i] + hp[i]);
        return fromResult(Isosurface::extract(familyAst, 0.0, res, s.halfExtent, std::move(avg)), h);
    };
    auto needCone = [&]() -> bool {
        if (out.local.singular && out.local.weightsFound) return true;
        out.note = out.local.note.empty() ? "cone patch not applicable at P" : out.local.note;
        return false;
    };

    switch (method) {
    case Method::FieldAverage: out.mesh = fieldAverage(); break;
    case Method::VertexMidMinus: out.mesh = vertexMid(out.minus, bp); break;
    case Method::VertexMidPlus: out.mesh = vertexMid(out.plus, bm); break;
    case Method::ExactVertexMid:
        out.mesh = exactVertexMid(out.minus, closestOnField(fieldMinus, bm, h), closestOnField(fieldPlus, bp, h), {},
                                  out.correspondenceFailures);
        break;
    case Method::NormalMid:
        out.mesh = normalMid(out.minus, fieldMinus, fieldF, fieldPlus, closestOnField(fieldPlus, bp, h), h, {},
                             out.correspondenceFallbacks, out.correspondenceFailures);
        break;
    case Method::SdfAverage: out.mesh = sdfAverage(familyAst, res, s.halfExtent, hm, hp, bm, bp); break;
    case Method::HybridCone:
        out.mesh = vertexMid(out.minus, bp);
        if (needCone()) {
            const double r = hybridRadius(out.local, s.eps, h, s.patchCells, 0.6 * s.halfExtent);
            out.patch = conePatch(out.mesh, s.center, out.local.weights, r, h, nullptr);
        }
        break;
    case Method::SingularCone:
        out.mesh = fieldAverage();
        if (needCone())
            out.patch = conePatch(out.mesh, s.center, out.local.weights,
                                  patchRadius(out.local, h, s.patchCells, 0.6 * s.halfExtent), h, &fieldF);
        break;
    default: break;
    }
    if (!out.patch.note.empty()) out.note += (out.note.empty() ? "" : "; ") + out.patch.note;
    out.millis = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - t0).count();
    return out;
}

} // namespace iso_reconstruct
