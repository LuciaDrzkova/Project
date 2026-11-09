#pragma once

// Small double-precision geometry kit shared by the analysis headers:
// points, closest-point distances, a point KD-tree and a triangle BVH.

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <numeric>
#include <vector>

namespace iso_geometry {

struct Point3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

inline Point3 operator+(const Point3& a, const Point3& b) { return {a.x + b.x, a.y + b.y, a.z + b.z}; }
inline Point3 operator-(const Point3& a, const Point3& b) { return {a.x - b.x, a.y - b.y, a.z - b.z}; }
inline Point3 operator*(const Point3& a, double s) { return {a.x * s, a.y * s, a.z * s}; }
inline double dot(const Point3& a, const Point3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double normSquared(const Point3& p) { return dot(p, p); }
inline double norm(const Point3& p) { return std::sqrt(normSquared(p)); }
inline double distance(const Point3& a, const Point3& b) { return norm(a - b); }
inline double component(const Point3& p, int axis) { return axis == 0 ? p.x : (axis == 1 ? p.y : p.z); }
inline Point3 minPoint(const Point3& a, const Point3& b) {
    return {std::min(a.x, b.x), std::min(a.y, b.y), std::min(a.z, b.z)};
}
inline Point3 maxPoint(const Point3& a, const Point3& b) {
    return {std::max(a.x, b.x), std::max(a.y, b.y), std::max(a.z, b.z)};
}

inline Point3 closestPointOnSegment(const Point3& p, const Point3& a, const Point3& b) {
    const Point3 ab = b - a;
    const double den = normSquared(ab);
    if (den <= 1e-30) return a;
    return a + ab * std::clamp(dot(p - a, ab) / den, 0.0, 1.0);
}

// Closest point to p on triangle abc (Ericson, "Real-Time Collision
// Detection", closest point by Voronoi region), with a segment fallback for
// degenerate triangles.
inline Point3 closestPointOnTriangle(const Point3& p, const Point3& a, const Point3& b, const Point3& c) {
    const Point3 ab = b - a;
    const Point3 ac = c - a;
    const Point3 ap = p - a;
    const double d1 = dot(ab, ap);
    const double d2 = dot(ac, ap);
    if (d1 <= 0.0 && d2 <= 0.0) return a;

    const Point3 bp = p - b;
    const double d3 = dot(ab, bp);
    const double d4 = dot(ac, bp);
    if (d3 >= 0.0 && d4 <= d3) return b;

    const double vc = d1 * d4 - d3 * d2;
    if (vc <= 0.0 && d1 >= 0.0 && d3 <= 0.0) return a + ab * (d1 / (d1 - d3));

    const Point3 cp = p - c;
    const double d5 = dot(ab, cp);
    const double d6 = dot(ac, cp);
    if (d6 >= 0.0 && d5 <= d6) return c;

    const double vb = d5 * d2 - d1 * d6;
    if (vb <= 0.0 && d2 >= 0.0 && d6 <= 0.0) return a + ac * (d2 / (d2 - d6));

    const double va = d3 * d6 - d5 * d4;
    if (va <= 0.0 && (d4 - d3) >= 0.0 && (d5 - d6) >= 0.0)
        return b + (c - b) * ((d4 - d3) / ((d4 - d3) + (d5 - d6)));

    const double denom = va + vb + vc;
    if (std::fabs(denom) <= 1e-30) {
        const Point3 candidates[3] = {closestPointOnSegment(p, a, b), closestPointOnSegment(p, b, c),
                                      closestPointOnSegment(p, c, a)};
        Point3 best = candidates[0];
        for (const Point3& q : candidates)
            if (normSquared(p - q) < normSquared(p - best)) best = q;
        return best;
    }
    const double inv = 1.0 / denom;
    return a + ab * (vb * inv) + ac * (vc * inv);
}

inline double pointTriangleDistanceSquared(const Point3& p, const Point3& a, const Point3& b, const Point3& c) {
    return normSquared(p - closestPointOnTriangle(p, a, b, c));
}

// Chebyshev (maximum-norm) distance: |p - q|_inf <= r  <=>  q lies in the
// axis-aligned box of half-size r centred at p.
inline double maxNormDistance(const Point3& p, const Point3& q) {
    return std::max({std::fabs(p.x - q.x), std::fabs(p.y - q.y), std::fabs(p.z - q.z)});
}

inline double pointBoxDistanceSquared(const Point3& p, const Point3& lo, const Point3& hi) {
    double s = 0.0;
    for (int axis = 0; axis < 3; ++axis) {
        const double v = component(p, axis), l = component(lo, axis), h = component(hi, axis);
        const double d = v < l ? l - v : (v > h ? v - h : 0.0);
        s += d * d;
    }
    return s;
}

// Nearest-neighbour queries over a fixed point set.
class KdTree {
public:
    explicit KdTree(std::vector<Point3> points) : points_(std::move(points)) {
        nodes_.reserve(points_.size());
        root_ = points_.empty() ? -1 : build(0, points_.size(), 0);
    }

    double nearestDistance(const Point3& query) const {
        if (root_ < 0) return std::numeric_limits<double>::infinity();
        double best2 = std::numeric_limits<double>::infinity();
        Point3 point{};
        nearest(root_, query, best2, point);
        return std::sqrt(best2);
    }

    // Nearest stored point to the query; false for an empty set.
    bool nearestPoint(const Point3& query, Point3& out) const {
        if (root_ < 0) return false;
        double best2 = std::numeric_limits<double>::infinity();
        nearest(root_, query, best2, out);
        return true;
    }

private:
    struct NodeEntry { Point3 p; int left = -1; int right = -1; int axis = 0; };
    std::vector<Point3> points_;
    std::vector<NodeEntry> nodes_;
    int root_ = -1;

    int build(std::size_t first, std::size_t last, int axis) {
        const std::size_t mid = first + (last - first) / 2;
        std::nth_element(points_.begin() + static_cast<std::ptrdiff_t>(first),
                         points_.begin() + static_cast<std::ptrdiff_t>(mid),
                         points_.begin() + static_cast<std::ptrdiff_t>(last),
                         [axis](const Point3& a, const Point3& b) {
                             return component(a, axis) < component(b, axis);
                         });
        const int id = static_cast<int>(nodes_.size());
        nodes_.push_back({points_[mid], -1, -1, axis});
        const int nextAxis = (axis + 1) % 3;
        if (first < mid) {
            const int child = build(first, mid, nextAxis);
            nodes_[static_cast<std::size_t>(id)].left = child;
        }
        if (mid + 1 < last) {
            const int child = build(mid + 1, last, nextAxis);
            nodes_[static_cast<std::size_t>(id)].right = child;
        }
        return id;
    }

    void nearest(int nodeId, const Point3& query, double& best2, Point3& bestPoint) const {
        if (nodeId < 0) return;
        const NodeEntry& node = nodes_[static_cast<std::size_t>(nodeId)];
        const double d2 = normSquared(query - node.p);
        if (d2 < best2) { best2 = d2; bestPoint = node.p; }
        const double delta = component(query, node.axis) - component(node.p, node.axis);
        nearest(delta < 0.0 ? node.left : node.right, query, best2, bestPoint);
        if (delta * delta < best2) nearest(delta < 0.0 ? node.right : node.left, query, best2, bestPoint);
    }
};

struct Triangle { Point3 a, b, c, centroid, lo, hi; };

inline Triangle makeTriangle(const Point3& a, const Point3& b, const Point3& c) {
    Triangle t;
    t.a = a; t.b = b; t.c = c;
    t.centroid = (a + b + c) * (1.0 / 3.0);
    t.lo = minPoint(a, minPoint(b, c));
    t.hi = maxPoint(a, maxPoint(b, c));
    return t;
}

// Exact point-to-triangle-set distance, accelerated by a median-split BVH.
class TriangleBvh {
public:
    explicit TriangleBvh(std::vector<Triangle> triangles) : triangles_(std::move(triangles)) {
        order_.resize(triangles_.size());
        std::iota(order_.begin(), order_.end(), std::size_t{0});
        if (!order_.empty()) {
            nodes_.reserve(order_.size() * 2);
            root_ = build(0, order_.size());
        }
    }

    bool empty() const { return triangles_.empty(); }

    double nearestDistance(const Point3& p) const {
        if (root_ < 0) return std::numeric_limits<double>::infinity();
        double best = std::numeric_limits<double>::infinity();
        Point3 point{};
        query(root_, p, best, point);
        return std::sqrt(best);
    }

    // Closest point of the triangle set to p; false for an empty set.
    bool nearestPoint(const Point3& p, Point3& out) const {
        if (root_ < 0) return false;
        double best = std::numeric_limits<double>::infinity();
        query(root_, p, best, out);
        return std::isfinite(best);
    }

private:
    struct BvhNode { Point3 lo, hi; int left = -1, right = -1; std::size_t begin = 0, count = 0; };
    std::vector<Triangle> triangles_;
    std::vector<std::size_t> order_;
    std::vector<BvhNode> nodes_;
    int root_ = -1;

    int build(std::size_t begin, std::size_t end) {
        constexpr double inf = std::numeric_limits<double>::infinity();
        BvhNode n;
        n.begin = begin;
        n.count = end - begin;
        n.lo = {inf, inf, inf};
        n.hi = {-inf, -inf, -inf};
        Point3 clo = n.lo, chi = n.hi;
        for (std::size_t i = begin; i < end; ++i) {
            const Triangle& t = triangles_[order_[i]];
            n.lo = minPoint(n.lo, t.lo);
            n.hi = maxPoint(n.hi, t.hi);
            clo = minPoint(clo, t.centroid);
            chi = maxPoint(chi, t.centroid);
        }
        const int id = static_cast<int>(nodes_.size());
        nodes_.push_back(n);
        if (end - begin <= 8) return id;

        const Point3 span = chi - clo;
        int axis = span.y > span.x ? 1 : 0;
        if (span.z > component(span, axis)) axis = 2;
        const std::size_t mid = begin + (end - begin) / 2;
        std::nth_element(order_.begin() + static_cast<std::ptrdiff_t>(begin),
                         order_.begin() + static_cast<std::ptrdiff_t>(mid),
                         order_.begin() + static_cast<std::ptrdiff_t>(end),
                         [&](std::size_t a, std::size_t b) {
                             return component(triangles_[a].centroid, axis) <
                                    component(triangles_[b].centroid, axis);
                         });
        const int left = build(begin, mid);
        const int right = build(mid, end);
        BvhNode& self = nodes_[static_cast<std::size_t>(id)];
        self.left = left;
        self.right = right;
        self.count = 0;
        return id;
    }

    void query(int id, const Point3& p, double& best, Point3& bestPoint) const {
        const BvhNode& n = nodes_[static_cast<std::size_t>(id)];
        if (pointBoxDistanceSquared(p, n.lo, n.hi) >= best) return;
        if (n.left < 0) {
            for (std::size_t i = n.begin; i < n.begin + n.count; ++i) {
                const Triangle& t = triangles_[order_[i]];
                const Point3 q = closestPointOnTriangle(p, t.a, t.b, t.c);
                const double d = normSquared(p - q);
                if (d < best) { best = d; bestPoint = q; }
            }
            return;
        }
        const BvhNode& a = nodes_[static_cast<std::size_t>(n.left)];
        const BvhNode& b = nodes_[static_cast<std::size_t>(n.right)];
        const bool leftFirst = pointBoxDistanceSquared(p, a.lo, a.hi) < pointBoxDistanceSquared(p, b.lo, b.hi);
        query(leftFirst ? n.left : n.right, p, best, bestPoint);
        query(leftFirst ? n.right : n.left, p, best, bestPoint);
    }
};

} // namespace iso_geometry
