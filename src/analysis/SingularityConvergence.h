#pragma once

// Numerical diagnostics for the family S_k = {F(x,y,z,k)=0} near a point p.
//
// IMPORTANT: These measurements are experimental mesh diagnostics, not proofs
// of local Hausdorff convergence or of smoothness/singularity. The input AST
// may describe any parameter dependence supported by Expr.h; for the thesis
// model F+kG, the expression should be affine in k.

#include "core/CsvFormat.h"
#include "core/Expr.h"
#include "core/Geometry.h"
#include "core/Isosurface.h"
#include "core/Parallel.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <ctime>
#include <fstream>
#include <functional>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace iso_convergence {

using Point3 = iso_geometry::Point3;
using iso_geometry::distance;
using iso_geometry::dot;
using iso_geometry::norm;
using iso_geometry::normSquared;

struct Options {
    Point3 probe{0.0, 0.0, 0.0};
    int resolution = 80;
    double halfExtent = 1.0;
    // Half-size of the local box Q = P + [-a, a]^3 in which all local
    // quantities are measured.
    double localHalfSize = 0.4;
    // Backwards-compatible mode: analyse both signs of each positive magnitude.
    // Used when kValues is empty.
    std::vector<double> magnitudes{1e-1, 3e-2, 1e-2, 3e-3, 1e-3, 3e-4};
    // Optional explicit k samples. When nonempty, these values are analysed as-is
    // (rather than generating +/- pairs from magnitudes). The GUI fills this from
    // the user's selected [kMinimum, kMaximum] and sample count.
    std::vector<double> kValues{};
    // Exact expression text entered in the viewer, e.g. "x*x+y*y-z*z+k*y".
    // It is recorded verbatim and decomposed into F(x,y,z)+kG(x,y,z) when
    // the AST can be shown to be structurally affine in k.
    std::string familyExpression = "UNKNOWN_NOT_SUPPLIED";
};

// Generate exactly sampleCount linearly spaced samples, including both endpoints.
// This makes the UI's range/count controls predictable; samples are not log-spaced.
inline std::vector<double> makeLinearKValues(double kMinimum, double kMaximum, int sampleCount) {
    if (!std::isfinite(kMinimum) || !std::isfinite(kMaximum))
        throw std::invalid_argument("The k range endpoints must be finite.");
    if (!(kMinimum < kMaximum))
        throw std::invalid_argument("k minimum must be smaller than k maximum.");
    if (sampleCount < 3 || sampleCount > 101)
        throw std::invalid_argument("Number of k samples must be between 3 and 101.");
    std::vector<double> values;
    values.reserve(static_cast<std::size_t>(sampleCount));
    const long double lo = static_cast<long double>(kMinimum);
    const long double hi = static_cast<long double>(kMaximum);
    for (int i = 0; i < sampleCount; ++i) {
        if (i == 0) values.push_back(kMinimum);
        else if (i == sampleCount - 1) values.push_back(kMaximum);
        else {
            const long double t = static_cast<long double>(i) / static_cast<long double>(sampleCount - 1);
            values.push_back(static_cast<double>(lo + (hi - lo) * t));
        }
    }
    for (std::size_t i = 1; i < values.size(); ++i) {
        if (!(values[i] > values[i - 1]))
            throw std::invalid_argument("The selected k range is too narrow for distinct double-precision samples.");
    }
    return values;
}

struct PointEvaluation {
    double value = std::numeric_limits<double>::quiet_NaN();
    Point3 gradient{};
    double gradientNorm = std::numeric_limits<double>::quiet_NaN();
};

struct Sample {
    double k = 0.0;
    std::size_t totalVertices = 0;
    std::size_t localVertices = 0;
    std::size_t triangles = 0;
    std::size_t exactZeroGridSamples = 0;
    std::size_t skippedDegenerateTriangles = 0;
    double meshBuildMilliseconds = std::numeric_limits<double>::quiet_NaN();
    double meshMaximumGradientMagnitude = std::numeric_limits<double>::quiet_NaN();
    bool hasSurfaceInLocalBox = false;
    double distanceProbeToMesh = std::numeric_limits<double>::infinity();
    double minimumGradientLocal = std::numeric_limits<double>::quiet_NaN();
    double maxAbsFieldResidualLocal = std::numeric_limits<double>::quiet_NaN();
    double meshToReference = std::numeric_limits<double>::quiet_NaN();
    double referenceToMesh = std::numeric_limits<double>::quiet_NaN();
    double symmetricHausdorff = std::numeric_limits<double>::quiet_NaN();
    bool geometryComparable = false;
};

struct MetricTrend {
    std::size_t usableSamples = 0;
    double logLogSlope = std::numeric_limits<double>::quiet_NaN();
    double coarseValue = std::numeric_limits<double>::quiet_NaN();
    double fineValue = std::numeric_limits<double>::quiet_NaN();
    std::string classification;
};

struct SideSummary {
    int sign = 0; // -1 or +1
    MetricTrend pointDistance;
    MetricTrend localHausdorff;
};

struct RunMetadata {
    std::string timestampUtc = "UNKNOWN";
    std::string familyExpression = "UNKNOWN";
    std::string fExpression = "UNKNOWN";
    std::string gExpression = "UNKNOWN";
    std::string kAffineStatus = "UNKNOWN";
    std::string degreeF = "UNKNOWN";
    std::string degreeG = "UNKNOWN";
    std::string degreeStatus = "UNKNOWN";
    std::string polynomialDegreeConvention = "total degree in x,y,z; k treated as parameter";

    Point3 probe{};
    double probeFAtZero = std::numeric_limits<double>::quiet_NaN();
    double probeGradientNormAtZero = std::numeric_limits<double>::quiet_NaN();

    double domainXMin = 0.0, domainXMax = 0.0;
    double domainYMin = 0.0, domainYMax = 0.0;
    double domainZMin = 0.0, domainZMax = 0.0;
    int gridResolutionCellsPerAxis = 0;
    int gridSamplesPerAxis = 0;
    double gridSpacingX = 0.0, gridSpacingY = 0.0, gridSpacingZ = 0.0;
    double localHalfSize = 0.0;

    double kMinimum = std::numeric_limits<double>::quiet_NaN();
    double kMaximum = std::numeric_limits<double>::quiet_NaN();
    std::size_t kSampleCount = 0;
    std::size_t negativeKSampleCount = 0;
    std::size_t zeroKSampleCount = 0;
    std::size_t positiveKSampleCount = 0;
    bool kValuesUniformlySpaced = false;
    double kSpacingIfUniform = std::numeric_limits<double>::quiet_NaN();

    double referenceK = 0.0;
    std::size_t referenceTotalVertices = 0;
    std::size_t referenceLocalVertices = 0;
    std::size_t referenceTriangles = 0;
    std::size_t referenceExactZeroGridSamples = 0;
    std::size_t referenceSkippedDegenerateTriangles = 0;

    std::string scalarFieldGrid = "regular Cartesian grid, nodes at -h+2h*i/N; scalar values evaluated in double precision";
    std::string isosurfaceExtractor = "Marching Cubes, classic 256-case lookup table";
    std::string edgeInterpolation = "linear interpolation along grid edges, each edge interpolated once from its lower to its upper node";
    std::string ambiguityHandling = "standard lookup-table resolution; no asymptotic decider";
    std::string exactZeroGridHandling = "grid nodes with F exactly 0, or |F| <= 1e-10 times the largest finite face-neighbour |F|, count as on-surface; an on-surface node with two or more face neighbours of opposite sign receives a deterministic positive tie-break bias of 5% of the smallest finite nonzero off-surface face-neighbour magnitude; original analytic field is retained for gradient evaluation and residual computation";
    double exactZeroBiasFraction = Isosurface::kZeroBiasFraction;
    std::string degenerateTriangleHandling = "reject a candidate triangle when squared cross-product norm < (step^2 * 1e-8)^2; rejected candidate count is recorded";
    std::string gradientEvaluation = "analytic forward-mode automatic differentiation (Dual3) at extracted vertices";
    std::string meshRepresentation = "non-indexed triangle soup; three SurfaceVertex records per triangle; counts are triangle-corner occurrences, not welded unique vertices";
    std::string pointToMeshDistance = "minimum Euclidean point-to-triangle distance over all mesh triangles";
    std::string localHausdorffDistance = "sampled-vertex proxy: directed nearest-neighbour distances between local mesh vertices using a KD-tree; symmetric value is the maximum of the two directed values; not exact triangle-surface Hausdorff distance";
    std::string localVertexSelection = "mesh vertices inside the probe-centred box P + [-a,a]^3, a = local_box_half_size (maximum-norm distance <= a)";
    std::string domainBoundaryTreatment = "mesh clipped to the axis-aligned cube; local comparison uses vertices inside the probe-centred box";

    std::string referenceMeshSettings = "k=0 reference mesh uses exactly the same grid resolution, domain, extractor, and local box as perturbed meshes";
    std::string numericTypes = "field evaluation and distances: double; stored mesh position and gradMag fields: float";
    double doubleEpsilon = std::numeric_limits<double>::epsilon();
    double floatEpsilon = std::numeric_limits<float>::epsilon();
    int doubleDigits10 = std::numeric_limits<double>::digits10;
    int floatDigits10 = std::numeric_limits<float>::digits10;
    std::string csvDoublePrecision = "17 significant digits (round-trip precision for IEEE-754 double in typical implementations)";
    std::string compiler = "UNKNOWN";
    std::string platform = "UNKNOWN";
    std::string cxxStandard = "UNKNOWN";
    std::string buildMode = "UNKNOWN";
    std::string sourceRevision = "not embedded by build system";
    std::string applicationVersion = "not embedded by build system";
};

struct Analysis {
    Options options;
    PointEvaluation probeAtZero;
    std::vector<Sample> samples;
    SideSummary negativeSide;
    SideSummary positiveSide;
    double gridSpacing = 0.0;
    RunMetadata metadata;
};

inline PointEvaluation evaluate(const Node* ast, const Point3& p, double k) {
    const Dual3 q = evalExpr<Dual3>(ast, Dual3(p.x, 1.0, 0.0, 0.0), Dual3(p.y, 0.0, 1.0, 0.0),
                                    Dual3(p.z, 0.0, 0.0, 1.0), Dual3(k, 0.0, 0.0, 0.0));
    PointEvaluation out;
    out.value = q.v;
    out.gradient = {q.dx, q.dy, q.dz};
    out.gradientNorm = std::sqrt(q.dx * q.dx + q.dy * q.dy + q.dz * q.dz);
    return out;
}

// ----- Symbolic helpers: expression text, affine-in-k split, degrees --------

inline std::string numberText(double value) {
    std::ostringstream out;
    out << std::setprecision(17) << value;
    return out.str();
}

namespace detail {
inline std::string parenthesize(const std::string& value) { return "(" + value + ")"; }
inline bool isZeroExpr(const std::string& value) { return value == "0" || value == "0.0"; }
inline bool isOneExpr(const std::string& value) { return value == "1" || value == "1.0"; }
inline std::string negExpr(const std::string& value) {
    if (isZeroExpr(value)) return "0";
    if (value.size() > 1 && value.front() == '-') return value.substr(1);
    return "-" + parenthesize(value);
}
inline std::string addExpr(const std::string& a, const std::string& b) {
    if (isZeroExpr(a)) return b;
    if (isZeroExpr(b)) return a;
    if (a == b) return "2*" + parenthesize(a);
    return parenthesize(a) + "+" + parenthesize(b);
}
inline std::string subExpr(const std::string& a, const std::string& b) {
    if (isZeroExpr(b)) return a;
    if (a == b) return "0";
    if (isZeroExpr(a)) return "-" + parenthesize(b);
    return parenthesize(a) + "-" + parenthesize(b);
}
inline std::string mulExpr(const std::string& a, const std::string& b) {
    if (isZeroExpr(a) || isZeroExpr(b)) return "0";
    if (isOneExpr(a)) return b;
    if (isOneExpr(b)) return a;
    return parenthesize(a) + "*" + parenthesize(b);
}
inline std::string divExpr(const std::string& a, const std::string& b) {
    if (isZeroExpr(a)) return "0";
    if (isOneExpr(b)) return a;
    return parenthesize(a) + "/" + parenthesize(b);
}
inline std::string powExpr(const std::string& a, const std::string& b) {
    if (isZeroExpr(b)) return "1";
    if (isOneExpr(b)) return a;
    return parenthesize(a) + "^" + parenthesize(b);
}
} // namespace detail

inline std::string expressionText(const Node* node, bool substituteK = false, double kValue = 0.0) {
    if (!node) return "UNKNOWN";
    auto kid = [&](std::size_t i) {
        return detail::parenthesize(expressionText(node->kids[i].get(), substituteK, kValue));
    };
    switch (node->type) {
        case Node::NUM: return numberText(node->num);
        case Node::VAR_X: return "x";
        case Node::VAR_Y: return "y";
        case Node::VAR_Z: return "z";
        case Node::VAR_K: return substituteK ? numberText(kValue) : "k";
        case Node::NEG: return "-" + kid(0);
        case Node::ADD: return kid(0) + "+" + kid(1);
        case Node::SUB: return kid(0) + "-" + kid(1);
        case Node::MUL: return kid(0) + "*" + kid(1);
        case Node::DIV: return kid(0) + "/" + kid(1);
        case Node::POW: return kid(0) + "^" + kid(1);
        case Node::CALL: {
            std::string result = node->fname + "(";
            for (std::size_t i = 0; i < node->kids.size(); ++i) {
                if (i) result += ",";
                result += expressionText(node->kids[i].get(), substituteK, kValue);
            }
            return result + ")";
        }
    }
    return "UNKNOWN";
}

struct AffineParts {
    bool valid = true;
    std::string f = "0";
    std::string g = "0";
};

// Syntactic affine-in-k decomposition. It rejects k in denominators, nonlinear
// products, nontrivial k-powers, and function arguments containing k. This is
// deliberately conservative: it may reject expressions that simplify to an
// affine expression after algebraic cancellation.
inline AffineParts decomposeAffineK(const Node* node) {
    using namespace detail;
    const AffineParts invalid{false, "UNKNOWN", "UNKNOWN"};
    if (!node) return invalid;
    switch (node->type) {
        case Node::NUM: return {true, numberText(node->num), "0"};
        case Node::VAR_X: return {true, "x", "0"};
        case Node::VAR_Y: return {true, "y", "0"};
        case Node::VAR_Z: return {true, "z", "0"};
        case Node::VAR_K: return {true, "0", "1"};
        case Node::NEG: {
            const auto a = decomposeAffineK(node->kids[0].get());
            return {a.valid, negExpr(a.f), negExpr(a.g)};
        }
        case Node::ADD:
        case Node::SUB: {
            const auto a = decomposeAffineK(node->kids[0].get());
            const auto b = decomposeAffineK(node->kids[1].get());
            if (!a.valid || !b.valid) return invalid;
            if (node->type == Node::ADD) return {true, addExpr(a.f, b.f), addExpr(a.g, b.g)};
            return {true, subExpr(a.f, b.f), subExpr(a.g, b.g)};
        }
        case Node::MUL: {
            const auto a = decomposeAffineK(node->kids[0].get());
            const auto b = decomposeAffineK(node->kids[1].get());
            if (!a.valid || !b.valid || !isZeroExpr(mulExpr(a.g, b.g))) return invalid;
            return {true, mulExpr(a.f, b.f), addExpr(mulExpr(a.f, b.g), mulExpr(a.g, b.f))};
        }
        case Node::DIV: {
            const auto a = decomposeAffineK(node->kids[0].get());
            const auto b = decomposeAffineK(node->kids[1].get());
            if (!a.valid || !b.valid || !isZeroExpr(b.g)) return invalid;
            return {true, divExpr(a.f, b.f), divExpr(a.g, b.f)};
        }
        case Node::POW: {
            const auto a = decomposeAffineK(node->kids[0].get());
            const auto b = decomposeAffineK(node->kids[1].get());
            if (!a.valid || !b.valid || !isZeroExpr(b.g)) return invalid;
            // If the base has a k coefficient, only powers 0 and 1 preserve affinity.
            if (!isZeroExpr(a.g)) {
                try {
                    const double exponent = std::stod(b.f);
                    if (exponent == 0.0) return {true, "1", "0"};
                    if (exponent == 1.0) return a;
                } catch (...) {}
                return invalid;
            }
            return {true, powExpr(a.f, b.f), "0"};
        }
        case Node::CALL: {
            std::string f = node->fname + "(";
            for (std::size_t i = 0; i < node->kids.size(); ++i) {
                const auto part = decomposeAffineK(node->kids[i].get());
                if (!part.valid || !isZeroExpr(part.g)) return invalid;
                if (i) f += ",";
                f += part.f;
            }
            return {true, f + ")", "0"};
        }
    }
    return invalid;
}

inline bool dependsOnXYZ(const Node* node) {
    return exprContains(node, Node::VAR_X) || exprContains(node, Node::VAR_Y) ||
           exprContains(node, Node::VAR_Z);
}

struct DegreeInfo { bool polynomial = false; int degree = -1; };

inline DegreeInfo spatialPolynomialDegree(const Node* node) {
    if (!node) return {};
    auto constantZero = [](const Node* n) {
        if (dependsOnXYZ(n)) return false;
        try { return evalExpr<double>(n, 0.0, 0.0, 0.0, 0.0) == 0.0; }
        catch (...) { return false; }
    };
    switch (node->type) {
        case Node::NUM: return {true, 0};
        case Node::VAR_X: case Node::VAR_Y: case Node::VAR_Z: return {true, 1};
        case Node::VAR_K: return {true, 0};
        case Node::NEG: return spatialPolynomialDegree(node->kids[0].get());
        case Node::ADD: case Node::SUB: {
            if (constantZero(node->kids[0].get())) return spatialPolynomialDegree(node->kids[1].get());
            if (constantZero(node->kids[1].get())) return spatialPolynomialDegree(node->kids[0].get());
            const auto a = spatialPolynomialDegree(node->kids[0].get());
            const auto b = spatialPolynomialDegree(node->kids[1].get());
            if (!a.polynomial || !b.polynomial) return {};
            return {true, std::max(a.degree, b.degree)};
        }
        case Node::MUL: {
            if (constantZero(node->kids[0].get()) || constantZero(node->kids[1].get())) return {true, 0};
            const auto a = spatialPolynomialDegree(node->kids[0].get());
            const auto b = spatialPolynomialDegree(node->kids[1].get());
            if (!a.polynomial || !b.polynomial) return {};
            return {true, a.degree + b.degree};
        }
        case Node::DIV: {
            const auto denominator = spatialPolynomialDegree(node->kids[1].get());
            if (!denominator.polynomial || denominator.degree != 0 || constantZero(node->kids[1].get())) return {};
            return spatialPolynomialDegree(node->kids[0].get());
        }
        case Node::POW: {
            const auto base = spatialPolynomialDegree(node->kids[0].get());
            if (!base.polynomial || dependsOnXYZ(node->kids[1].get())) return {};
            double exponent = 0.0;
            try { exponent = evalExpr<double>(node->kids[1].get(), 0.0, 0.0, 0.0, 0.0); }
            catch (...) { return {}; }
            if (!std::isfinite(exponent) || exponent < 0.0 || std::floor(exponent) != exponent || exponent > 1000.0)
                return {};
            return {true, base.degree * static_cast<int>(exponent)};
        }
        case Node::CALL:
            if (!dependsOnXYZ(node)) return {true, 0};
            return {};
    }
    return {};
}

inline std::string degreeText(const Node* expression) {
    if (expression && !dependsOnXYZ(expression)) {
        try {
            if (evalExpr<double>(expression, 0.0, 0.0, 0.0, 0.0) == 0.0)
                return "ZERO_POLYNOMIAL (degree undefined)";
        } catch (...) {}
    }
    const DegreeInfo info = spatialPolynomialDegree(expression);
    return info.polynomial ? std::to_string(info.degree) : "UNKNOWN_NONPOLYNOMIAL_OR_RATIONAL";
}

namespace detail {

struct MeshPoint {
    Point3 p;
    double gradient = std::numeric_limits<double>::quiet_NaN();
    double absFieldResidual = std::numeric_limits<double>::quiet_NaN();
};

// Mesh vertices inside the probe box. `points` holds each welded vertex once
// (duplicated soup corners cannot change a min/max or a nearest-neighbour
// distance); `cornerCount` keeps the documented corner-occurrence count.
struct LocalPatch {
    std::size_t cornerCount = 0;
    std::vector<MeshPoint> points;
};

inline Point3 vertexPoint(const SurfaceVertex& v) {
    return {static_cast<double>(v.px), static_cast<double>(v.py), static_cast<double>(v.pz)};
}

inline LocalPatch localPatch(const Node* ast, const IsosurfaceResult& mesh, double k,
                             const Point3& center, double halfSize) {
    LocalPatch out;
    const bool indexed = mesh.vertexIds.size() == mesh.vertices.size();
    std::vector<char> seen(indexed ? mesh.uniqueVertexCount : 0, 0);
    for (std::size_t corner = 0; corner < mesh.vertices.size(); ++corner) {
        const SurfaceVertex& v = mesh.vertices[corner];
        const Point3 p = vertexPoint(v);
        if (iso_geometry::maxNormDistance(p, center) > halfSize) continue;
        ++out.cornerCount;
        if (indexed) {
            char& flag = seen[mesh.vertexIds[corner]];
            if (flag) continue;
            flag = 1;
        }
        const double fValue = evalExpr<double>(ast, p.x, p.y, p.z, k);
        out.points.push_back({p, static_cast<double>(v.gradMag), std::fabs(fValue)});
    }
    return out;
}

inline iso_geometry::KdTree makeTree(const LocalPatch& patch) {
    std::vector<Point3> points;
    points.reserve(patch.points.size());
    for (const MeshPoint& p : patch.points) points.push_back(p.p);
    return iso_geometry::KdTree(std::move(points));
}

inline double pointMeshDistance(const IsosurfaceResult& mesh, const Point3& p) {
    if (mesh.vertices.size() < 3) return std::numeric_limits<double>::infinity();
    double best = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i + 2 < mesh.vertices.size(); i += 3) {
        best = std::min(best, iso_geometry::pointTriangleDistanceSquared(
            p, vertexPoint(mesh.vertices[i]), vertexPoint(mesh.vertices[i + 1]),
            vertexPoint(mesh.vertices[i + 2])));
        if (best <= 1e-28) return 0.0;
    }
    return std::sqrt(best);
}

inline double directedHausdorff(const std::vector<MeshPoint>& source, const iso_geometry::KdTree& target) {
    if (source.empty()) return std::numeric_limits<double>::quiet_NaN();
    double largest = 0.0;
    for (const MeshPoint& p : source) {
        const double d = target.nearestDistance(p.p);
        if (!std::isfinite(d)) return std::numeric_limits<double>::quiet_NaN();
        largest = std::max(largest, d);
    }
    return largest;
}

inline double linearLogSlope(const std::vector<std::pair<double, double>>& data) {
    if (data.size() < 3) return std::numeric_limits<double>::quiet_NaN();
    double mx = 0.0, my = 0.0;
    for (const auto& xy : data) { mx += xy.first; my += xy.second; }
    mx /= static_cast<double>(data.size());
    my /= static_cast<double>(data.size());
    double numerator = 0.0, denominator = 0.0;
    for (const auto& xy : data) {
        numerator += (xy.first - mx) * (xy.second - my);
        denominator += (xy.first - mx) * (xy.first - mx);
    }
    if (denominator <= 1e-30) return std::numeric_limits<double>::quiet_NaN();
    return numerator / denominator;
}

inline MetricTrend trendFor(const std::vector<Sample>& samples, int sign, bool useHausdorff) {
    std::vector<const Sample*> sideSamples;
    for (const Sample& s : samples) {
        if (s.k == 0.0) continue; // k=0 is the reference fibre, not a one-sided sample.
        if ((s.k > 0.0 ? 1 : -1) != sign) continue;
        // A finite distance to a mesh elsewhere in the plotting domain is not
        // evidence of local approach if the surface does not meet the local box.
        if (!useHausdorff && !s.hasSurfaceInLocalBox) continue;
        sideSamples.push_back(&s);
    }
    // Fit against decreasing |k| irrespective of how the user ordered the range.
    std::sort(sideSamples.begin(), sideSamples.end(), [](const Sample* a, const Sample* b) {
        return std::fabs(a->k) > std::fabs(b->k);
    });
    std::vector<std::pair<double, double>> logData;
    std::vector<double> values;
    for (const Sample* s : sideSamples) {
        const double value = useHausdorff ? s->symmetricHausdorff : s->distanceProbeToMesh;
        if (!std::isfinite(value) || value <= 1e-12) continue;
        logData.emplace_back(std::log(std::fabs(s->k)), std::log(value));
        values.push_back(value);
    }

    MetricTrend out;
    out.usableSamples = values.size();
    if (!values.empty()) {
        out.coarseValue = values.front();
        out.fineValue = values.back();
    }
    out.logLogSlope = linearLogSlope(logData);
    if (values.size() < 3) {
        out.classification = "insufficient positive, resolved samples";
        return out;
    }

    // Samples are sorted by decreasing |k|, so a positive slope means the
    // distance decreases as |k| decreases. This is an empirical heuristic only.
    const double ratio = out.fineValue / out.coarseValue;
    if (std::isfinite(out.logLogSlope) && out.logLogSlope > 0.15 && ratio < 0.8)
        out.classification = "distance decreases over sampled range (heuristic)";
    else if (std::isfinite(out.logLogSlope) && out.logLogSlope >= -0.05 && ratio < 1.25)
        out.classification = "plateau / likely resolution-limited";
    else
        out.classification = "no stable decrease observed";
    return out;
}

// ----- Reproducibility metadata ---------------------------------------------

inline std::string utcTimestamp() {
    const std::time_t now = std::time(nullptr);
    std::tm tmValue{};
#if defined(_WIN32)
    gmtime_s(&tmValue, &now);
#else
    gmtime_r(&now, &tmValue);
#endif
    char buffer[32]{};
    std::strftime(buffer, sizeof(buffer), "%Y-%m-%dT%H:%M:%SZ", &tmValue);
    return buffer;
}

inline std::string compilerText() {
#if defined(__clang__)
    return std::string("Clang ") + __clang_version__;
#elif defined(__GNUC__)
    return std::string("GCC ") + __VERSION__;
#elif defined(_MSC_VER)
    return std::string("MSVC ") + std::to_string(_MSC_VER);
#else
    return "UNKNOWN_COMPILER";
#endif
}

inline std::string platformText() {
#if defined(_WIN32)
    return "Windows";
#elif defined(__APPLE__)
    return "macOS";
#elif defined(__linux__)
    return "Linux";
#else
    return "unknown_platform";
#endif
}

inline std::string buildModeText() {
#ifdef NDEBUG
    return "NDEBUG defined (typically Release; exact build flags not embedded)";
#else
    return "NDEBUG not defined (Debug or unspecified; exact build flags not embedded)";
#endif
}

inline std::string sourceRevisionText() {
#ifdef ISO_PROJECT_GIT_COMMIT
    return ISO_PROJECT_GIT_COMMIT;
#else
    return "not embedded by build system";
#endif
}

inline std::string applicationVersionText() {
#ifdef ISO_PROJECT_VERSION
    return ISO_PROJECT_VERSION;
#else
    return "not embedded by build system";
#endif
}

inline RunMetadata makeMetadata(const Node* ast, const Options& options,
                                const PointEvaluation& probeAtZero,
                                const IsosurfaceResult& reference,
                                std::size_t referenceLocalCorners,
                                double gridSpacing) {
    RunMetadata m;
    m.timestampUtc = utcTimestamp();
    m.familyExpression = options.familyExpression.empty() ? "UNKNOWN_NOT_SUPPLIED" : options.familyExpression;
    m.probe = options.probe;
    m.probeFAtZero = probeAtZero.value;
    m.probeGradientNormAtZero = probeAtZero.gradientNorm;
    m.domainXMin = m.domainYMin = m.domainZMin = -options.halfExtent;
    m.domainXMax = m.domainYMax = m.domainZMax = options.halfExtent;
    m.gridResolutionCellsPerAxis = options.resolution;
    m.gridSamplesPerAxis = options.resolution + 1;
    m.gridSpacingX = m.gridSpacingY = m.gridSpacingZ = gridSpacing;
    m.localHalfSize = options.localHalfSize;
    m.kSampleCount = options.kValues.size();
    if (!options.kValues.empty()) {
        const auto bounds = std::minmax_element(options.kValues.begin(), options.kValues.end());
        m.kMinimum = *bounds.first;
        m.kMaximum = *bounds.second;
        for (double k : options.kValues) {
            if (k < 0.0) ++m.negativeKSampleCount;
            else if (k > 0.0) ++m.positiveKSampleCount;
            else ++m.zeroKSampleCount;
        }
        if (options.kValues.size() >= 2) {
            const double step = options.kValues[1] - options.kValues[0];
            const double tolerance = std::max(1.0, std::fabs(step)) * 1e-10;
            bool uniform = true;
            for (std::size_t i = 2; i < options.kValues.size(); ++i)
                if (std::fabs((options.kValues[i] - options.kValues[i - 1]) - step) > tolerance) uniform = false;
            m.kValuesUniformlySpaced = uniform;
            if (uniform) m.kSpacingIfUniform = step;
        }
    }
    m.referenceK = 0.0;
    m.referenceTotalVertices = reference.vertices.size();
    m.referenceLocalVertices = referenceLocalCorners;
    m.referenceTriangles = reference.vertices.size() / 3;
    m.referenceExactZeroGridSamples = reference.exactZeroGridSamples;
    m.referenceSkippedDegenerateTriangles = reference.skippedDegenerateTriangles;

    const AffineParts affine = decomposeAffineK(ast);
    if (affine.valid) {
        m.kAffineStatus = "STRUCTURALLY_AFFINE_IN_K";
        m.fExpression = affine.f;
        m.gExpression = affine.g;
        try {
            const NodePtr fAst = parseExpression(m.fExpression);
            const NodePtr gAst = parseExpression(m.gExpression);
            m.degreeF = degreeText(fAst.get());
            m.degreeG = degreeText(gAst.get());
            m.degreeStatus = "total spatial polynomial degrees computed from the affine AST coefficients; degrees are UNKNOWN when a coefficient is not polynomial";
        } catch (...) {
            m.degreeStatus = "could not parse derived affine coefficient expressions";
        }
    } else {
        m.kAffineStatus = "NOT_VERIFIED_AS_AFFINE_IN_K (conservative structural check failed)";
        m.fExpression = expressionText(ast, true, 0.0);
        m.gExpression = "NOT_AVAILABLE: family is not structurally affine in k";
        try {
            const NodePtr fAst = parseExpression(m.fExpression);
            m.degreeF = degreeText(fAst.get());
        } catch (...) {}
        m.degreeG = "UNKNOWN";
        m.degreeStatus = "G coefficient is not reported because F_k=F+kG could not be verified structurally";
    }

    m.compiler = compilerText();
    m.platform = platformText();
    m.cxxStandard = std::to_string(static_cast<long long>(__cplusplus));
    m.buildMode = buildModeText();
    m.sourceRevision = sourceRevisionText();
    m.applicationVersion = applicationVersionText();
    return m;
}

// Validates the options and returns the distinct k samples to analyse.
inline std::vector<double> resolveKValues(Options& options) {
    if (!(options.halfExtent > 0.0) || !std::isfinite(options.halfExtent))
        throw std::invalid_argument("halfExtent must be finite and positive.");
    if (!(options.localHalfSize > 0.0) || !std::isfinite(options.localHalfSize))
        throw std::invalid_argument("localHalfSize must be finite and positive.");
    std::vector<double> kValues;
    if (!options.kValues.empty()) {
        // Preserve the range order for the results table and CSV, removing only
        // exact duplicate values (possible when endpoints are nearly identical).
        for (double kValue : options.kValues) {
            if (!std::isfinite(kValue))
                throw std::invalid_argument("Every selected k value must be finite.");
            if (std::find(kValues.begin(), kValues.end(), kValue) == kValues.end())
                kValues.push_back(kValue);
        }
    } else {
        if (options.magnitudes.empty())
            throw std::invalid_argument("Provide explicit k values or at least one |k| magnitude.");
        for (double a : options.magnitudes) {
            if (!(a > 0.0) || !std::isfinite(a))
                throw std::invalid_argument("Every k magnitude must be finite and positive.");
        }
        std::sort(options.magnitudes.begin(), options.magnitudes.end(), std::greater<double>());
        options.magnitudes.erase(std::unique(options.magnitudes.begin(), options.magnitudes.end()),
                                 options.magnitudes.end());
        for (double magnitude : options.magnitudes) {
            kValues.push_back(-magnitude);
            kValues.push_back(+magnitude);
        }
    }
    if (kValues.empty())
        throw std::invalid_argument("No distinct k samples were provided.");
    return kValues;
}

} // namespace detail

// Called once per analysed k (from worker threads, possibly concurrently) with
// the sample index, k, that k's mesh and the k=0 reference mesh. It lets other
// analyses reuse each mesh instead of extracting it a second time; it must only
// write to per-index storage.
using MeshVisitor = std::function<void(std::size_t sampleIndex, double k, const IsosurfaceResult& mesh,
                                       const IsosurfaceResult& reference)>;

inline Analysis analyze(const Node* ast, const Options& inputOptions, const MeshVisitor& visitMesh = {}) {
    if (!ast) throw std::invalid_argument("Convergence analysis needs a parsed expression.");
    Options options = inputOptions;
    options.resolution = std::clamp(options.resolution, 8, Isosurface::kMaxAnalysisResolution);
    options.kValues = detail::resolveKValues(options);

    Analysis out;
    out.options = options;
    out.gridSpacing = 2.0 * options.halfExtent / options.resolution;
    out.probeAtZero = evaluate(ast, options.probe, 0.0);

    const IsosurfaceResult reference = Isosurface::build(ast, 0.0, options.resolution, options.halfExtent);
    const detail::LocalPatch referenceLocal =
        detail::localPatch(ast, reference, 0.0, options.probe, options.localHalfSize);
    const iso_geometry::KdTree referenceTree = detail::makeTree(referenceLocal);
    out.metadata = detail::makeMetadata(ast, options, out.probeAtZero, reference,
                                        referenceLocal.cornerCount, out.gridSpacing);

    // Every k is independent: extract and measure them in parallel. Each
    // worker owns one mesh at a time, and results land in per-index slots so
    // the table/CSV order is the requested k order regardless of scheduling.
    out.samples.resize(options.kValues.size());
    iso_parallel::parallelFor(options.kValues.size(), [&](std::size_t index) {
        const double k = options.kValues[index];
        const IsosurfaceResult mesh = Isosurface::build(ast, k, options.resolution, options.halfExtent);
        const detail::LocalPatch local = detail::localPatch(ast, mesh, k, options.probe, options.localHalfSize);

        Sample& row = out.samples[index];
        row.k = k;
        row.totalVertices = mesh.vertices.size();
        row.localVertices = local.cornerCount;
        row.triangles = mesh.vertices.size() / 3;
        row.exactZeroGridSamples = mesh.exactZeroGridSamples;
        row.skippedDegenerateTriangles = mesh.skippedDegenerateTriangles;
        row.meshBuildMilliseconds = mesh.buildMillis;
        row.meshMaximumGradientMagnitude = mesh.maxGradMag;
        row.distanceProbeToMesh = detail::pointMeshDistance(mesh, options.probe);
        row.hasSurfaceInLocalBox = !local.points.empty();
        double maxResidual = 0.0;
        bool sawFiniteResidual = false;
        for (const detail::MeshPoint& p : local.points) {
            if (std::isfinite(p.gradient) &&
                (!std::isfinite(row.minimumGradientLocal) || p.gradient < row.minimumGradientLocal))
                row.minimumGradientLocal = p.gradient;
            if (std::isfinite(p.absFieldResidual)) {
                maxResidual = std::max(maxResidual, p.absFieldResidual);
                sawFiniteResidual = true;
            }
        }
        if (sawFiniteResidual) row.maxAbsFieldResidualLocal = maxResidual;

        if (!local.points.empty() && !referenceLocal.points.empty()) {
            const iso_geometry::KdTree localTree = detail::makeTree(local);
            row.meshToReference = detail::directedHausdorff(local.points, referenceTree);
            row.referenceToMesh = detail::directedHausdorff(referenceLocal.points, localTree);
            if (std::isfinite(row.meshToReference) && std::isfinite(row.referenceToMesh)) {
                row.symmetricHausdorff = std::max(row.meshToReference, row.referenceToMesh);
                row.geometryComparable = true;
            }
        }
        if (visitMesh) visitMesh(index, k, mesh, reference);
    });

    out.negativeSide.sign = -1;
    out.negativeSide.pointDistance = detail::trendFor(out.samples, -1, false);
    out.negativeSide.localHausdorff = detail::trendFor(out.samples, -1, true);
    out.positiveSide.sign = +1;
    out.positiveSide.pointDistance = detail::trendFor(out.samples, +1, false);
    out.positiveSide.localHausdorff = detail::trendFor(out.samples, +1, true);
    return out;
}

inline bool writeCsv(const std::string& path, const Analysis& a, std::string* error = nullptr) {
    using iso_csv::number;
    using iso_csv::quote;
    std::ofstream out(path);
    if (!out) {
        if (error) *error = "Could not open output file: " + path;
        return false;
    }
    const RunMetadata& m = a.metadata;
    out << "sample_index,k,abs_k,sign,total_vertices,local_vertices,triangles,exact_zero_grid_samples,"
           "skipped_degenerate_triangles,mesh_build_milliseconds,mesh_maximum_gradient_magnitude,"
           "has_surface_in_local_box,distance_probe_to_mesh,min_gradient_local,max_abs_field_residual_local,mesh_to_reference,"
           "reference_to_mesh,symmetric_local_hausdorff,geometry_comparable,"
           "metadata_timestamp_utc,metadata_family_expression,metadata_F_expression,metadata_G_expression,"
           "metadata_k_affine_status,metadata_degree_F,metadata_degree_G,metadata_degree_status,"
           "metadata_degree_convention,metadata_probe_x,metadata_probe_y,metadata_probe_z,"
           "metadata_probe_F_at_k0,metadata_probe_gradient_norm_at_k0,"
           "metadata_domain_x_min,metadata_domain_x_max,metadata_domain_y_min,metadata_domain_y_max,"
           "metadata_domain_z_min,metadata_domain_z_max,metadata_grid_resolution_cells_per_axis,"
           "metadata_grid_samples_per_axis,metadata_grid_spacing_x,metadata_grid_spacing_y,metadata_grid_spacing_z,"
           "metadata_local_box_half_size,metadata_k_min,metadata_k_max,metadata_k_sample_count,"
           "metadata_negative_k_sample_count,metadata_zero_k_sample_count,metadata_positive_k_sample_count,"
           "metadata_k_values_uniformly_spaced,metadata_k_spacing_if_uniform,metadata_reference_k,metadata_reference_total_vertices,metadata_reference_local_vertices,"
           "metadata_reference_triangles,metadata_reference_exact_zero_grid_samples,"
           "metadata_reference_skipped_degenerate_triangles,metadata_scalar_field_grid,"
           "metadata_isosurface_extractor,metadata_edge_interpolation,metadata_ambiguity_handling,"
           "metadata_exact_zero_grid_handling,metadata_exact_zero_bias_fraction,metadata_degenerate_triangle_handling,"
           "metadata_gradient_evaluation,metadata_mesh_representation,metadata_point_to_mesh_distance,"
           "metadata_local_hausdorff_distance,metadata_local_vertex_selection,metadata_domain_boundary_treatment,"
           "metadata_reference_mesh_settings,metadata_numeric_types,metadata_double_epsilon,"
           "metadata_float_epsilon,metadata_double_digits10,metadata_float_digits10,metadata_csv_double_precision,"
           "metadata_compiler,metadata_platform,metadata_cxx_standard,metadata_build_mode,"
           "metadata_source_revision,metadata_application_version\n";

    // The metadata block is identical on every row: format it once.
    std::ostringstream meta;
    meta << quote(m.timestampUtc) << ',' << quote(m.familyExpression) << ','
         << quote(m.fExpression) << ',' << quote(m.gExpression) << ','
         << quote(m.kAffineStatus) << ',' << quote(m.degreeF) << ','
         << quote(m.degreeG) << ',' << quote(m.degreeStatus) << ','
         << quote(m.polynomialDegreeConvention) << ','
         << number(m.probe.x) << ',' << number(m.probe.y) << ',' << number(m.probe.z) << ','
         << number(m.probeFAtZero) << ',' << number(m.probeGradientNormAtZero) << ','
         << number(m.domainXMin) << ',' << number(m.domainXMax) << ','
         << number(m.domainYMin) << ',' << number(m.domainYMax) << ','
         << number(m.domainZMin) << ',' << number(m.domainZMax) << ','
         << m.gridResolutionCellsPerAxis << ',' << m.gridSamplesPerAxis << ','
         << number(m.gridSpacingX) << ',' << number(m.gridSpacingY) << ','
         << number(m.gridSpacingZ) << ',' << number(m.localHalfSize) << ','
         << number(m.kMinimum) << ',' << number(m.kMaximum) << ',' << m.kSampleCount << ','
         << m.negativeKSampleCount << ',' << m.zeroKSampleCount << ',' << m.positiveKSampleCount << ','
         << (m.kValuesUniformlySpaced ? 1 : 0) << ',' << number(m.kSpacingIfUniform) << ','
         << number(m.referenceK) << ',' << m.referenceTotalVertices << ','
         << m.referenceLocalVertices << ',' << m.referenceTriangles << ','
         << m.referenceExactZeroGridSamples << ',' << m.referenceSkippedDegenerateTriangles << ','
         << quote(m.scalarFieldGrid) << ',' << quote(m.isosurfaceExtractor) << ','
         << quote(m.edgeInterpolation) << ',' << quote(m.ambiguityHandling) << ','
         << quote(m.exactZeroGridHandling) << ',' << number(m.exactZeroBiasFraction) << ','
         << quote(m.degenerateTriangleHandling) << ','
         << quote(m.gradientEvaluation) << ',' << quote(m.meshRepresentation) << ','
         << quote(m.pointToMeshDistance) << ',' << quote(m.localHausdorffDistance) << ','
         << quote(m.localVertexSelection) << ',' << quote(m.domainBoundaryTreatment) << ','
         << quote(m.referenceMeshSettings) << ',' << quote(m.numericTypes) << ','
         << number(m.doubleEpsilon) << ',' << number(m.floatEpsilon) << ','
         << m.doubleDigits10 << ',' << m.floatDigits10 << ',' << quote(m.csvDoublePrecision) << ','
         << quote(m.compiler) << ',' << quote(m.platform) << ','
         << quote(m.cxxStandard) << ',' << quote(m.buildMode) << ','
         << quote(m.sourceRevision) << ',' << quote(m.applicationVersion);
    const std::string metadataColumns = meta.str();

    for (std::size_t index = 0; index < a.samples.size(); ++index) {
        const Sample& s = a.samples[index];
        out << index << ',' << number(s.k) << ',' << number(std::fabs(s.k)) << ','
            << (s.k < 0.0 ? -1 : (s.k > 0.0 ? 1 : 0)) << ','
            << s.totalVertices << ',' << s.localVertices << ',' << s.triangles << ','
            << s.exactZeroGridSamples << ',' << s.skippedDegenerateTriangles << ','
            << number(s.meshBuildMilliseconds) << ',' << number(s.meshMaximumGradientMagnitude) << ','
            << (s.hasSurfaceInLocalBox ? 1 : 0) << ','
            << number(s.distanceProbeToMesh) << ',' << number(s.minimumGradientLocal) << ','
            << number(s.maxAbsFieldResidualLocal) << ',' << number(s.meshToReference) << ','
            << number(s.referenceToMesh) << ',' << number(s.symmetricHausdorff) << ','
            << (s.geometryComparable ? 1 : 0) << ',' << metadataColumns << '\n';
    }
    if (!out) {
        if (error) *error = "An error occurred while writing: " + path;
        return false;
    }
    return true;
}

} // namespace iso_convergence
