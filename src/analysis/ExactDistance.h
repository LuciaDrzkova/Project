#pragma once

// Distances measured on the exact surfaces S_k = {H_k = 0}, not on meshes.
//
// Comparing two extracted meshes mixes the change of the surface with the
// extraction error of both meshes. Here the meshes are used only to find
// starting points:
//   * Exact samples of S_k inside the local box Q: mesh points of S_k moved
//     onto {H_k = 0} by Newton's method along the gradient, plus the singular
//     points of S_k found by the certified stationary-point search.
//   * Distance from a point a to an exact surface S: start at the closest
//     point of S's mesh, move it onto S, then solve the Lagrange conditions of
//     min |p - a|^2 subject to H(p) = 0 by Newton's method. Every candidate is
//     an exact point of S (and the singular points of S are candidates too),
//     so the result is an upper bound of dist(a, S) that is attained up to the
//     choice of the starting point.
// The same machinery measures the geometric accuracy of the extraction: the
// distance from mesh points to the exact surface and back.

#include "CubicPolynomial.h"
#include "core/Geometry.h"
#include "core/Isosurface.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <optional>
#include <vector>

namespace iso_exact {

using Point3 = iso_geometry::Point3;
using iso_geometry::distance;
using iso_geometry::dot;
using iso_geometry::norm;
constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// H_k as an exact cubic when possible, otherwise through the expression tree.
class Field {
public:
    Field(const Node* ast, double k) : ast_(ast), k_(k), cubic_(iso_poly::Cubic::fromExpression(ast, k)) {}

    iso_poly::Derivatives<double> evaluate(const Point3& p) const {
        if (cubic_) return cubic_->evaluate(p);
        const iso_research::DifferentialData d = iso_research::evaluateJet(ast_, p, k_);
        iso_poly::Derivatives<double> out;
        out.value = d.value;
        out.gradient = {d.gradient.x, d.gradient.y, d.gradient.z};
        out.hessian = d.hessian;
        return out;
    }
    double termScale(const Point3& p) const {
        return cubic_ ? cubic_->termScale(p) : iso_research::evalAbsScale(ast_, p, k_);
    }
    bool exactPolynomial() const { return cubic_.has_value(); }

private:
    const Node* ast_;
    double k_;
    std::optional<iso_poly::Cubic> cubic_;
};

namespace detail {

inline Point3 toPoint(const std::array<double, 3>& a) { return {a[0], a[1], a[2]}; }

// Newton's method along the gradient onto {H = 0}; steps are capped at maxStep.
// Most points converge in a few steps; across a thin feature (a spike whose
// width is far below the step cap) the iteration zig-zags and needs more.
inline bool project(const Field& f, Point3& p, double maxStep) {
    for (int it = 0; it < 150; ++it) {
        const auto d = f.evaluate(p);
        const Point3 g = toPoint(d.gradient);
        const double g2 = dot(g, g);
        if (!(g2 > 0.0) || !std::isfinite(d.value)) return false;
        Point3 step = g * (d.value / g2);
        const double length = norm(step);
        if (length > maxStep) step = step * (maxStep / length);
        p = p - step;
        if (length <= 1e-15 * (1.0 + norm(p))) break;
    }
    const double value = f.evaluate(p).value;
    return std::fabs(value) <= 1e-10 * std::max(f.termScale(p), 1e-300);
}

// 4x4 Gaussian elimination with partial pivoting.
inline bool solve4(std::array<std::array<double, 5>, 4> m, std::array<double, 4>& x) {
    for (int c = 0; c < 4; ++c) {
        int p = c;
        for (int r = c + 1; r < 4; ++r) if (std::fabs(m[r][c]) > std::fabs(m[p][c])) p = r;
        if (std::fabs(m[p][c]) < 1e-300) return false;
        std::swap(m[p], m[c]);
        for (int r = c + 1; r < 4; ++r) {
            const double f = m[r][c] / m[c][c];
            for (int j = c; j < 5; ++j) m[r][j] -= f * m[c][j];
        }
    }
    for (int r = 3; r >= 0; --r) {
        double s = m[r][4];
        for (int j = r + 1; j < 4; ++j) s -= m[r][j] * x[j];
        x[r] = s / m[r][r];
    }
    return std::isfinite(x[0]) && std::isfinite(x[1]) && std::isfinite(x[2]) && std::isfinite(x[3]);
}

// Newton on the Lagrange system  p - a + lambda grad H(p) = 0,  H(p) = 0,
// started from a point p of the surface. Returns the stationary point of the
// distance (normally the closest point), or nothing.
inline std::optional<Point3> lagrange(const Field& f, const Point3& a, Point3 p, double maxStep) {
    auto d = f.evaluate(p);
    Point3 g = toPoint(d.gradient);
    if (!(dot(g, g) > 0.0)) return std::nullopt;
    double lambda = -dot(p - a, g) / dot(g, g);
    for (int it = 0; it < 30; ++it) {
        d = f.evaluate(p);
        g = toPoint(d.gradient);
        const Point3 r = p - a + g * lambda;
        std::array<std::array<double, 5>, 4> m{};
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) m[i][j] = (i == j ? 1.0 : 0.0) + lambda * d.hessian[i][j];
            m[i][3] = d.gradient[i];
            m[3][i] = d.gradient[i];
            m[i][4] = iso_geometry::component(r, i);
        }
        m[3][3] = 0.0;
        m[3][4] = d.value;
        std::array<double, 4> delta{};
        if (!solve4(m, delta)) return std::nullopt;
        Point3 step{delta[0], delta[1], delta[2]};
        const double length = norm(step);
        if (length > maxStep) { step = step * (maxStep / length); delta[3] *= maxStep / length; }
        p = p - step;
        lambda -= delta[3];
        if (length <= 1e-14 * (1.0 + norm(p))) break;
    }
    if (!project(f, p, maxStep)) return std::nullopt;
    return p;
}

} // namespace detail

// Everything needed to measure distances to one exact surface.
struct Surface {
    const Field* field = nullptr;
    std::vector<Point3> samples;          // exact points of the surface inside the box
    std::vector<Point3> singularPoints;   // exact singular points (all of them, not only in the box)
    iso_geometry::TriangleBvh mesh{{}};   // the whole extracted mesh, for starting points
    iso_geometry::KdTree exactPoints{{}}; // every exact point found (box plus one cell), second starting points
    std::size_t failedProjections = 0;
    double step = 0.0;                    // grid spacing: cap for Newton steps
};

inline std::vector<iso_geometry::Triangle> triangles(const IsosurfaceResult& mesh) {
    std::vector<iso_geometry::Triangle> out;
    out.reserve(mesh.vertices.size() / 3);
    for (std::size_t i = 0; i + 2 < mesh.vertices.size(); i += 3)
        out.push_back(iso_geometry::makeTriangle(iso_research::surfacePoint(mesh.vertices[i]),
                                                 iso_research::surfacePoint(mesh.vertices[i + 1]),
                                                 iso_research::surfacePoint(mesh.vertices[i + 2])));
    return out;
}

// Mesh points used as seeds/samples: the vertices and the centroid of every
// triangle in the box, plus a finer barycentric grid (6 subdivisions) on the
// triangles within three cells of the centre, where the deformation is
// concentrated.
inline std::vector<Point3> meshSamples(const IsosurfaceResult& mesh, const Point3& center, double halfSize, double h) {
    std::vector<Point3> out;
    const double core = 3.0 * h;
    for (std::size_t i = 0; i + 2 < mesh.vertices.size(); i += 3) {
        const Point3 a = iso_research::surfacePoint(mesh.vertices[i]);
        const Point3 b = iso_research::surfacePoint(mesh.vertices[i + 1]);
        const Point3 c = iso_research::surfacePoint(mesh.vertices[i + 2]);
        const Point3 centroid = (a + b + c) * (1.0 / 3.0);
        if (iso_geometry::maxNormDistance(centroid, center) > halfSize + h) continue;
        if (iso_geometry::maxNormDistance(centroid, center) <= core) {
            constexpr int n = 6;
            for (int u = 0; u <= n; ++u)
                for (int v = 0; v <= n - u; ++v)
                    out.push_back(a * (double(u) / n) + b * (double(v) / n) + c * (double(n - u - v) / n));
        } else {
            out.push_back(a);
            out.push_back(centroid);
        }
    }
    return out;
}

// `extraSeeds`: further points near the surface (from finer local meshes) that
// are moved onto it and used like the mesh samples.
// Graded refinement towards the centre (the singular point). Features of S_k
// near P can be much smaller than the grid spacing h (their size is the inner
// scale lambda(k) of the thesis); a mesh with spacing h misses them, and so
// would the samples. The surface is therefore extracted again on nested boxes
// centre + [-a/4^l, a/4^l]^3, l = 1..levels, with `resolution` cells each
// (for a = 0.18: spacings down to 4.4e-5),
// and the vertices of these local meshes serve as extra seeds. Only the
// samples are refined; the extraction accuracy is still measured on the
// original mesh. Requires the centre to be the origin of the coordinates
// (the extractor's grid is centred there); otherwise no seeds are produced.
inline std::vector<Point3> gradedSeeds(const Node* ast, double k, const Point3& center, double halfSize,
                                       int levels = 4, int resolution = 32) {
    std::vector<Point3> out;
    if (iso_geometry::norm(center) != 0.0) return out;
    double size = halfSize;
    for (int l = 1; l <= levels; ++l) {
        size *= 0.25;
        const IsosurfaceResult local = Isosurface::build(ast, k, resolution, size);
        // The mesh is a triangle soup: take every welded vertex once.
        const bool indexed = local.vertexIds.size() == local.vertices.size();
        std::vector<char> seen(indexed ? local.uniqueVertexCount : 0, 0);
        for (std::size_t i = 0; i < local.vertices.size(); ++i) {
            if (indexed) {
                char& flag = seen[local.vertexIds[i]];
                if (flag) continue;
                flag = 1;
            }
            out.push_back(iso_research::surfacePoint(local.vertices[i]));
        }
    }
    return out;
}

// Seeds at the real stationary points of H. By Morse theory the level sets
// change their topology only at stationary points, so thin features that a
// mesh misses (a spike, a small bubble, a narrow neck) start there. From the
// point itself and from points 10^-6 * a away along each axis, the projection
// (Newton along the gradient, steps capped) runs along the feature onto it.
inline std::vector<Point3> stationarySeeds(const std::vector<Point3>& stationary, const Point3& center, double halfSize) {
    std::vector<Point3> out;
    const double delta = 1e-6 * halfSize;
    for (const Point3& q : stationary) {
        if (iso_geometry::maxNormDistance(q, center) > halfSize) continue;
        out.push_back(q);
        for (int axis = 0; axis < 3; ++axis)
            for (double sign : {-1.0, 1.0}) {
                Point3 p = q;
                (axis == 0 ? p.x : (axis == 1 ? p.y : p.z)) += sign * delta;
                out.push_back(p);
            }
    }
    return out;
}

inline Surface prepare(const Field& field, const IsosurfaceResult& mesh, const Point3& center, double halfSize,
                       double h, std::vector<Point3> singularPoints, const std::vector<Point3>& extraSeeds = {}) {
    Surface s;
    s.field = &field;
    s.step = h;
    s.mesh = iso_geometry::TriangleBvh(triangles(mesh));
    s.singularPoints = std::move(singularPoints);
    std::vector<Point3> all;
    std::vector<Point3> seeds = meshSamples(mesh, center, halfSize, h);
    seeds.insert(seeds.end(), extraSeeds.begin(), extraSeeds.end());
    for (Point3 p : seeds) {
        if (!detail::project(field, p, h)) { ++s.failedProjections; continue; }
        all.push_back(p);
        if (iso_geometry::maxNormDistance(p, center) <= halfSize) s.samples.push_back(p);
    }
    s.exactPoints = iso_geometry::KdTree(std::move(all));
    for (const Point3& q : s.singularPoints)
        if (iso_geometry::maxNormDistance(q, center) <= halfSize) s.samples.push_back(q);
    return s;
}

// Upper bound of dist(a, S), attained up to the starting points (see top).
// Two starts: the nearest point of the mesh (moved onto S) and the nearest
// exact sample of S; the singular points of S are candidates as well.
inline double distanceTo(const Surface& s, const Point3& a) {
    double best = std::numeric_limits<double>::infinity();
    for (const Point3& q : s.singularPoints) best = std::min(best, distance(a, q));
    auto refine = [&](const Point3& start) {
        best = std::min(best, distance(a, start));
        if (const auto p = detail::lagrange(*s.field, a, start, s.step)) best = std::min(best, distance(a, *p));
    };
    Point3 start;
    if (s.mesh.nearestPoint(a, start) && detail::project(*s.field, start, s.step)) refine(start);
    if (s.exactPoints.nearestPoint(a, start)) refine(start);
    return best;
}

struct Directed { double max = kNaN, rms = kNaN; std::size_t samples = 0, failures = 0; };

inline Directed directed(const std::vector<Point3>& source, const Surface& target) {
    Directed out;
    double maxd = 0.0, sum = 0.0;
    for (const Point3& a : source) {
        const double d = distanceTo(target, a);
        if (!std::isfinite(d)) { ++out.failures; continue; }
        maxd = std::max(maxd, d);
        sum += d * d;
        ++out.samples;
    }
    if (out.samples) { out.max = maxd; out.rms = std::sqrt(sum / static_cast<double>(out.samples)); }
    return out;
}

inline Directed directedToMesh(const std::vector<Point3>& source, const iso_geometry::TriangleBvh& mesh) {
    Directed out;
    double maxd = 0.0, sum = 0.0;
    for (const Point3& a : source) {
        const double d = mesh.nearestDistance(a);
        if (!std::isfinite(d)) { ++out.failures; continue; }
        maxd = std::max(maxd, d);
        sum += d * d;
        ++out.samples;
    }
    if (out.samples) { out.max = maxd; out.rms = std::sqrt(sum / static_cast<double>(out.samples)); }
    return out;
}

// All exact measurements for one k.
struct Measurement {
    double k = 0.0;
    std::size_t perturbedSamples = 0, referenceSamples = 0, failedProjections = 0;
    Directed perturbedToReference, referenceToPerturbed;   // S_k -> S_0 and S_0 -> S_k
    double symmetric = kNaN;                              // local Hausdorff distance d_Q(S_k, S_0)
    double probeDistance = kNaN;                          // dist(P, S_k)
    Directed meshToSurface, surfaceToMesh;                // extraction accuracy of the S_k mesh
    double meshHausdorff = kNaN;                          // max of the two
    double singularPointToMesh = kNaN;                    // max over singular points of S_k in the box
    bool comparable = false;
};

inline Measurement measure(double k, const Surface& perturbed, const Surface& reference,
                           const IsosurfaceResult& perturbedMesh, const Point3& center, double halfSize, double h) {
    Measurement m;
    m.k = k;
    m.perturbedSamples = perturbed.samples.size();
    m.referenceSamples = reference.samples.size();
    m.failedProjections = perturbed.failedProjections;
    m.probeDistance = distanceTo(perturbed, center);
    if (!perturbed.samples.empty() && !reference.samples.empty()) {
        m.perturbedToReference = directed(perturbed.samples, reference);
        m.referenceToPerturbed = directed(reference.samples, perturbed);
        if (m.perturbedToReference.samples && m.referenceToPerturbed.samples) {
            m.symmetric = std::max(m.perturbedToReference.max, m.referenceToPerturbed.max);
            m.comparable = true;
        }
    }
    // Extraction accuracy: mesh points (vertices and centroids in the box) to S_k, and back.
    std::vector<Point3> meshPoints;
    for (const Point3& p : meshSamples(perturbedMesh, center, halfSize, h))
        if (iso_geometry::maxNormDistance(p, center) <= halfSize) meshPoints.push_back(p);
    m.meshToSurface = directed(meshPoints, perturbed);
    m.surfaceToMesh = directedToMesh(perturbed.samples, perturbed.mesh);
    if (m.meshToSurface.samples && m.surfaceToMesh.samples)
        m.meshHausdorff = std::max(m.meshToSurface.max, m.surfaceToMesh.max);
    double worst = kNaN;
    for (const Point3& q : perturbed.singularPoints) {
        if (iso_geometry::maxNormDistance(q, center) > halfSize) continue;
        const double d = perturbed.mesh.nearestDistance(q);
        worst = std::isfinite(worst) ? std::max(worst, d) : d;
    }
    m.singularPointToMesh = worst;
    return m;
}

} // namespace iso_exact
