#pragma once

// Research-grade local diagnostics for F_k(x,y,z)=F(x,y,z)+kG(x,y,z).
//
// Uses the project's Expr.h, Isosurface.h and SingularityConvergence.h APIs.
// All numerical searches are explicitly labelled as non-certified.
//
// analyzeFamily() is the entry point used by the viewer and the batch runner:
// it runs the convergence analysis and this research analysis on one shared
// set of meshes (each k is extracted once, in parallel).

#include "core/CsvFormat.h"
#include "core/Geometry.h"
#include "core/Parallel.h"
#include "SingularityConvergence.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <fstream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace iso_research {

using Point3 = iso_geometry::Point3;
using iso_geometry::distance;
using iso_geometry::norm;

struct Exponent3 { int x, y, z; };

inline const std::array<Exponent3, 20>& exponents() {
    static const std::array<Exponent3, 20> e{{
        {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1},
        {2,0,0}, {1,1,0}, {1,0,1}, {0,2,0}, {0,1,1}, {0,0,2},
        {3,0,0}, {2,1,0}, {2,0,1}, {1,2,0}, {1,1,1},
        {1,0,2}, {0,3,0}, {0,2,1}, {0,1,2}, {0,0,3}
    }};
    return e;
}

inline int exponentIndex(int x, int y, int z) {
    const auto& e = exponents();
    for (int i = 0; i < static_cast<int>(e.size()); ++i) {
        const Exponent3& m = e[static_cast<std::size_t>(i)];
        if (m.x == x && m.y == y && m.z == z) return i;
    }
    return -1;
}

inline int factorial(int n) { return n < 2 ? 1 : (n == 2 ? 2 : 6); }

// productIndex[i][j]: coefficient slot of monomial_i * monomial_j, or -1 when
// the product's total degree exceeds 3 (truncated). Built once.
inline const std::array<std::array<std::int8_t, 20>, 20>& jetProductTable() {
    static const auto table = [] {
        std::array<std::array<std::int8_t, 20>, 20> t{};
        const auto& e = exponents();
        for (std::size_t i = 0; i < 20; ++i)
            for (std::size_t j = 0; j < 20; ++j) {
                const int x = e[i].x + e[j].x, y = e[i].y + e[j].y, z = e[i].z + e[j].z;
                t[i][j] = static_cast<std::int8_t>(x + y + z > 3 ? -1 : exponentIndex(x, y, z));
            }
        return t;
    }();
    return table;
}

// A normalized multivariate Taylor jet through total degree 3.
// c[i,j,k] stores (partial^(i+j+k) f)/(i!j!k!).
struct Jet3 {
    std::array<double, 20> c{};
    bool smooth = true;
    Jet3() = default;
    explicit Jet3(double value) { c[0] = value; }
    static Jet3 variable(double value, int axis) {
        Jet3 q(value);
        if (axis >= 0 && axis < 3) q.c[static_cast<std::size_t>(1 + axis)] = 1.0;
        return q;
    }
    double value() const { return c[0]; }
};

inline Jet3 operator+(const Jet3& a, const Jet3& b) {
    Jet3 r; r.smooth = a.smooth && b.smooth;
    for (std::size_t i = 0; i < r.c.size(); ++i) r.c[i] = a.c[i] + b.c[i];
    return r;
}
inline Jet3 operator-(const Jet3& a, const Jet3& b) {
    Jet3 r; r.smooth = a.smooth && b.smooth;
    for (std::size_t i = 0; i < r.c.size(); ++i) r.c[i] = a.c[i] - b.c[i];
    return r;
}
inline Jet3 operator-(const Jet3& a) {
    Jet3 r; r.smooth = a.smooth;
    for (std::size_t i = 0; i < r.c.size(); ++i) r.c[i] = -a.c[i];
    return r;
}
inline Jet3 operator*(const Jet3& a, const Jet3& b) {
    Jet3 r; r.smooth = a.smooth && b.smooth;
    const auto& product = jetProductTable();
    for (std::size_t i = 0; i < 20; ++i) {
        if (a.c[i] == 0.0) continue;
        const auto& row = product[i];
        for (std::size_t j = 0; j < 20; ++j) {
            if (b.c[j] == 0.0 || row[j] < 0) continue;
            r.c[static_cast<std::size_t>(row[j])] += a.c[i] * b.c[j];
        }
    }
    return r;
}
// Scalar multiple with the same zero-skipping as a product with a constant jet.
inline Jet3 scaled(const Jet3& a, double s) {
    Jet3 r; r.smooth = a.smooth;
    if (s == 0.0) return r;
    for (std::size_t i = 0; i < r.c.size(); ++i) r.c[i] = a.c[i] == 0.0 ? 0.0 : a.c[i] * s;
    return r;
}
// f(a) for a scalar function with derivatives f0..f3 at a.value().
inline Jet3 composeUnary(const Jet3& a, double f0, double f1, double f2, double f3) {
    Jet3 u = a; u.c[0] = 0.0;
    const Jet3 u2 = u * u;
    const Jet3 u3 = u2 * u;
    Jet3 r(f0);
    r = r + scaled(u, f1);
    r = r + scaled(u2, 0.5 * f2);
    r = r + scaled(u3, f3 / 6.0);
    r.smooth = a.smooth;
    return r;
}
inline Jet3 nonSmoothNaN() {
    Jet3 r(std::numeric_limits<double>::quiet_NaN());
    r.smooth = false;
    return r;
}
inline Jet3 reciprocal(const Jet3& a) {
    const double x = a.value();
    if (x == 0.0 || !std::isfinite(x)) return nonSmoothNaN();
    return composeUnary(a, 1.0 / x, -1.0 / (x * x), 2.0 / (x * x * x), -6.0 / (x * x * x * x));
}
inline Jet3 operator/(const Jet3& a, const Jet3& b) { return a * reciprocal(b); }
inline Jet3 mc_sin(const Jet3& a) {
    const double x = a.value();
    return composeUnary(a, std::sin(x), std::cos(x), -std::sin(x), -std::cos(x));
}
inline Jet3 mc_cos(const Jet3& a) {
    const double x = a.value();
    return composeUnary(a, std::cos(x), -std::sin(x), -std::cos(x), std::sin(x));
}
inline Jet3 mc_tan(const Jet3& a) {
    const double x = a.value(), t = std::tan(x), sec2 = 1.0 / (std::cos(x) * std::cos(x));
    return composeUnary(a, t, sec2, 2.0 * t * sec2, 2.0 * sec2 * (sec2 + 2.0 * t * t));
}
inline Jet3 mc_exp(const Jet3& a) {
    const double e = std::exp(a.value());
    return composeUnary(a, e, e, e, e);
}
inline Jet3 mc_log(const Jet3& a) {
    const double x = a.value();
    if (!(x > 0.0)) return nonSmoothNaN();
    return composeUnary(a, std::log(x), 1.0 / x, -1.0 / (x * x), 2.0 / (x * x * x));
}
inline Jet3 mc_sqrt(const Jet3& a) {
    const double x = a.value();
    if (!(x > 0.0)) {
        if (x == 0.0 && std::all_of(a.c.begin() + 1, a.c.end(), [](double v) { return v == 0.0; }))
            return Jet3(0.0);
        return nonSmoothNaN();
    }
    const double r = std::sqrt(x);
    return composeUnary(a, r, 0.5 / r, -0.25 / (x * r), 0.375 / (x * x * r));
}
inline Jet3 mc_abs(const Jet3& a) {
    if (a.value() == 0.0) { Jet3 r(0.0); r.smooth = false; return r; }
    return a.value() > 0.0 ? a : -a;
}
inline Jet3 mc_pow(const Jet3& a, const Jet3& b) {
    const bool exponentConstant = std::all_of(b.c.begin() + 1, b.c.end(), [](double v) { return v == 0.0; });
    if (!exponentConstant) {
        if (!(a.value() > 0.0)) return nonSmoothNaN();
        return mc_exp(b * mc_log(a));
    }
    const double p = b.value(), x = a.value();
    if (x == 0.0 && p < 0.0) return nonSmoothNaN();
    if (x < 0.0 && std::fabs(p - std::round(p)) > 1e-12) return nonSmoothNaN();
    const double f0 = std::pow(x, p);
    const double f1 = (p == 0.0) ? 0.0 : p * std::pow(x, p - 1.0);
    const double f2 = (p == 0.0 || p == 1.0) ? 0.0 : p * (p - 1.0) * std::pow(x, p - 2.0);
    const double f3 = (p == 0.0 || p == 1.0 || p == 2.0) ? 0.0 : p * (p - 1.0) * (p - 2.0) * std::pow(x, p - 3.0);
    return composeUnary(a, f0, f1, f2, f3);
}
inline Jet3 mc_min(const Jet3& a, const Jet3& b) {
    if (a.value() == b.value()) { Jet3 r = a; r.smooth = false; return r; }
    return a.value() < b.value() ? a : b;
}
inline Jet3 mc_max(const Jet3& a, const Jet3& b) {
    if (a.value() == b.value()) { Jet3 r = a; r.smooth = false; return r; }
    return a.value() > b.value() ? a : b;
}

struct DifferentialData {
    double value = std::numeric_limits<double>::quiet_NaN();
    Point3 gradient{};
    std::array<std::array<double, 3>, 3> hessian{};
    std::array<double, 10> third{}; // xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz
    std::array<double, 20> taylor{};
    bool smooth = false;
};

inline DifferentialData extractDifferential(const Jet3& j) {
    DifferentialData d;
    d.value = j.c[0];
    d.smooth = j.smooth;
    d.taylor = j.c;
    d.gradient = {j.c[1], j.c[2], j.c[3]};
    d.hessian[0][0] = 2 * j.c[4];
    d.hessian[0][1] = d.hessian[1][0] = j.c[5];
    d.hessian[0][2] = d.hessian[2][0] = j.c[6];
    d.hessian[1][1] = 2 * j.c[7];
    d.hessian[1][2] = d.hessian[2][1] = j.c[8];
    d.hessian[2][2] = 2 * j.c[9];
    const auto& e = exponents();
    for (std::size_t i = 10; i < 20; ++i)
        d.third[i - 10] = j.c[i] * factorial(e[i].x) * factorial(e[i].y) * factorial(e[i].z);
    return d;
}

inline DifferentialData evaluateJet(const Node* ast, const Point3& p, double k) {
    const Jet3 x = Jet3::variable(p.x, 0), y = Jet3::variable(p.y, 1), z = Jet3::variable(p.z, 2), kk(k);
    return extractDifferential(evalExpr<Jet3>(ast, x, y, z, kk));
}

inline double gradientNorm(const Point3& g) { return std::sqrt(g.x * g.x + g.y * g.y + g.z * g.z); }

inline double determinant3(const std::array<std::array<double, 3>, 3>& a) {
    return a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
         - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
         + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
}

// Jacobi eigenvalue iteration for a symmetric 3x3 matrix, sorted ascending.
inline std::array<double, 3> symmetricEigenvalues(std::array<std::array<double, 3>, 3> a) {
    for (int sweep = 0; sweep < 24; ++sweep) {
        int p = 0, q = 1;
        if (std::fabs(a[0][2]) > std::fabs(a[p][q])) { p = 0; q = 2; }
        if (std::fabs(a[1][2]) > std::fabs(a[p][q])) { p = 1; q = 2; }
        if (std::fabs(a[p][q]) < 1e-15) break;
        const double tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
        const double t = std::copysign(1.0, tau) / (std::fabs(tau) + std::sqrt(1.0 + tau * tau));
        const double c = 1.0 / std::sqrt(1.0 + t * t), s = t * c;
        const double app = a[p][p], aqq = a[q][q], apq = a[p][q];
        a[p][p] = app - t * apq;
        a[q][q] = aqq + t * apq;
        a[p][q] = a[q][p] = 0.0;
        for (int r = 0; r < 3; ++r) {
            if (r == p || r == q) continue;
            const double arp = a[r][p], arq = a[r][q];
            a[r][p] = a[p][r] = c * arp - s * arq;
            a[r][q] = a[q][r] = s * arp + c * arq;
        }
    }
    std::array<double, 3> e{{a[0][0], a[1][1], a[2][2]}};
    std::sort(e.begin(), e.end());
    return e;
}

inline int firstNonzeroJetOrder(const std::array<double, 20>& c) {
    double scale = 0.0;
    for (double v : c) scale = std::max(scale, std::fabs(v));
    if (scale == 0.0) return -1;
    const double tol = 512.0 * std::numeric_limits<double>::epsilon() * scale;
    for (std::size_t i = 0; i < c.size(); ++i) {
        if (std::fabs(c[i]) > tol) {
            const Exponent3 e = exponents()[i];
            return e.x + e.y + e.z;
        }
    }
    return -1;
}

// A conservative magnitude scale of the terms of F at p, used to make
// floating-point zero tests relative to the size of the cancelling terms.
inline double evalAbsScale(const Node* n, const Point3& p, double k) {
    if (!n) return 0.0;
    constexpr double kHuge = std::numeric_limits<double>::max() / 16.0;
    auto kid = [&](std::size_t i) { return evalAbsScale(n->kids[i].get(), p, k); };
    switch (n->type) {
        case Node::NUM: return std::fabs(n->num);
        case Node::VAR_X: return std::fabs(p.x);
        case Node::VAR_Y: return std::fabs(p.y);
        case Node::VAR_Z: return std::fabs(p.z);
        case Node::VAR_K: return std::fabs(k);
        case Node::NEG: return kid(0);
        case Node::ADD: case Node::SUB: return kid(0) + kid(1);
        case Node::MUL: return kid(0) * kid(1);
        case Node::DIV: return kid(0) / std::max(kid(1), 1e-300);
        case Node::POW: {
            const double b = evalExpr<double>(n->kids[1].get(), p.x, p.y, p.z, k);
            const double v = std::pow(std::max(kid(0), 1e-300), b);
            return std::isfinite(v) ? v : kHuge;
        }
        case Node::CALL: {
            const double a = kid(0);
            switch (n->func) {
                case Node::SIN: case Node::COS: return 1.0;
                case Node::EXP: { const double v = std::exp(a); return std::isfinite(v) ? v : kHuge; }
                case Node::SQRT: return std::sqrt(a);
                case Node::ABS: return a;
                case Node::LOG: return std::max(1.0, std::fabs(std::log(std::max(a, 1e-300))));
                default: break;
            }
            if (n->kids.size() == 2) return std::max(a, kid(1));
            return std::max(1.0, a);
        }
    }
    return 1.0;
}

struct IsolationAssessment {
    std::string status = "UNKNOWN";
    double fieldValue = std::numeric_limits<double>::quiet_NaN();
    double gradientMagnitude = std::numeric_limits<double>::quiet_NaN();
    double hessianDeterminant = std::numeric_limits<double>::quiet_NaN();
    std::array<double, 3> hessianEigenvalues{};
    std::string explanation;
};

inline double hessianMaxAbs(const std::array<std::array<double, 3>, 3>& h) {
    double m = 0.0;
    for (const auto& row : h) for (double v : row) m = std::max(m, std::fabs(v));
    return m;
}

inline IsolationAssessment assessReferenceSingularity(const Node* fAst, const Point3& p, double radius) {
    constexpr double eps = std::numeric_limits<double>::epsilon();
    IsolationAssessment out;
    const DifferentialData d = evaluateJet(fAst, p, 0.0);
    out.fieldValue = d.value;
    out.gradientMagnitude = gradientNorm(d.gradient);
    out.hessianDeterminant = determinant3(d.hessian);
    out.hessianEigenvalues = symmetricEigenvalues(d.hessian);
    const double fTol = 512.0 * eps * evalAbsScale(fAst, p, 0.0) + 1e-30;
    const double hScale = hessianMaxAbs(d.hessian);
    const double gTol = 1e-10 * std::max(hScale * std::max(radius, 1e-6), 1e-12);
    if (!d.smooth) {
        out.status = "NOT_CERTIFIED_NONSMOOTH_OR_UNSUPPORTED";
        out.explanation = "The expression is non-smooth or its third-order jet is unsupported at the probe point.";
        return out;
    }
    if (std::fabs(d.value) > fTol || out.gradientMagnitude > gTol) {
        out.status = "PROBE_NOT_CONFIRMED_SINGULAR";
        out.explanation = "The probe does not satisfy F(p)=0 and grad F(p)=0 within scale-aware floating-point tolerances.";
        return out;
    }
    const double detTol = 512.0 * eps * std::max(hScale * hScale * hScale, 1e-300);
    if (std::fabs(out.hessianDeterminant) > detTol) {
        out.status = "NONDEGENERATE_SINGULAR_POINT_LOCAL_ISOLATION_FOLLOWS";
        out.explanation = "The probe is numerically a singular point and det(H_F(p)) is nonzero. For exact input data, the inverse function theorem makes p an isolated zero of grad F, hence an isolated hypersurface singularity locally.";
        return out;
    }
    // A coarse witness scan may detect another singular point, but cannot prove
    // non-isolation. Degenerate Hessian cases remain explicitly inconclusive.
    int extraWitnesses = 0;
    for (int i = -2; i <= 2; ++i) for (int j = -2; j <= 2; ++j) for (int l = -2; l <= 2; ++l) {
        const Point3 q{p.x + radius * 0.5 * i, p.y + radius * 0.5 * j, p.z + radius * 0.5 * l};
        if (distance(q, p) < 1e-7 * std::max(1.0, radius)) continue;
        const DifferentialData v = evaluateJet(fAst, q, 0.0);
        const double ft = 512.0 * eps * evalAbsScale(fAst, q, 0.0) + 1e-30;
        if (std::fabs(v.value) <= ft && gradientNorm(v.gradient) <= 1e-10 * std::max(1.0, hScale)) ++extraWitnesses;
    }
    if (extraWitnesses > 0) {
        out.status = "NONISOLATED_SINGULAR_WITNESSES_DETECTED_NOT_CERTIFIED";
        out.explanation = "Additional grid witnesses approximately satisfy F=0 and grad F=0. This is evidence of a non-isolated singular locus, not a proof.";
    } else {
        out.status = "DEGENERATE_HESSIAN_ISOLATION_INCONCLUSIVE";
        out.explanation = "The Hessian is singular. A singular Hessian neither proves nor disproves isolatedness; higher-order/algebraic analysis is required.";
    }
    return out;
}

struct CriticalPoint {
    double k = 0.0;
    Point3 point{};
    double fieldValue = std::numeric_limits<double>::quiet_NaN();
    double fieldScale = std::numeric_limits<double>::quiet_NaN();
    double fieldTolerance = std::numeric_limits<double>::quiet_NaN();
    double gradientMagnitude = std::numeric_limits<double>::quiet_NaN();
    double gradientTolerance = std::numeric_limits<double>::quiet_NaN();
    double newtonStepTolerance = std::numeric_limits<double>::quiet_NaN();
    double deduplicationTolerance = std::numeric_limits<double>::quiet_NaN();
    double hessianDeterminant = std::numeric_limits<double>::quiet_NaN();
    int iterations = 0;
    bool probePointCandidate = false;
    std::string searchSource = "MULTISTART_DAMPED_NEWTON";
    std::string status;
};

struct CriticalPointTolerances {
    double field = std::numeric_limits<double>::denorm_min() * 64.0;
    double gradient = std::numeric_limits<double>::denorm_min() * 64.0;
    double newtonStep = 0.0;
    double deduplication = 0.0;
};

// These tolerances have different meanings and are intentionally not shared:
// field/gradient residuals test equations, step tolerance tests Newton progress,
// and deduplication tolerance decides whether two numerically found roots are
// the same candidate. No arbitrary 1e-30 field floor is added: that floor can
// incorrectly turn a small but nonzero critical value into a zero-level root.
inline CriticalPointTolerances criticalPointTolerances(const Node* ast, double k, const Point3& q,
                                                       const Point3& center, double radius,
                                                       const DifferentialData& d) {
    constexpr double eps = std::numeric_limits<double>::epsilon();
    const double denormFloor = 64.0 * std::numeric_limits<double>::denorm_min();
    const double safeRadius = std::max(std::fabs(radius), 1e-12);
    const double localLength = std::max(safeRadius, distance(q, center));
    const double coordinateScale = std::max({std::fabs(q.x), std::fabs(q.y), std::fabs(q.z), safeRadius});
    CriticalPointTolerances t;
    t.field = std::max(denormFloor, 512.0 * eps * evalAbsScale(ast, q, k));
    t.gradient = std::max(denormFloor, 512.0 * eps * hessianMaxAbs(d.hessian) * localLength);
    t.newtonStep = std::max(128.0 * eps * std::max(coordinateScale, 1e-12), 1e-12 * safeRadius);
    t.deduplication = std::max(512.0 * eps * std::max(coordinateScale, 1e-12), 1e-10 * safeRadius);
    return t;
}

inline double criticalPointQuality(const CriticalPoint& c) {
    const double floor = 64.0 * std::numeric_limits<double>::denorm_min();
    return std::hypot(c.gradientMagnitude / std::max(c.gradientTolerance, floor),
                      std::fabs(c.fieldValue) / std::max(c.fieldTolerance, floor));
}

inline void addOrImproveCriticalPoint(std::vector<CriticalPoint>& candidates, CriticalPoint candidate) {
    auto it = std::find_if(candidates.begin(), candidates.end(), [&](const CriticalPoint& existing) {
        const double tol = std::max(existing.deduplicationTolerance, candidate.deduplicationTolerance);
        return distance(existing.point, candidate.point) <= tol;
    });
    if (it == candidates.end()) {
        candidates.push_back(std::move(candidate));
        return;
    }
    // Keep the candidate with the smallest normalized pair of residuals. If
    // quality ties, prefer direct evaluation of the selected probe point.
    const double oldQuality = criticalPointQuality(*it);
    const double newQuality = criticalPointQuality(candidate);
    if (newQuality < oldQuality ||
        (newQuality == oldQuality && candidate.probePointCandidate && !it->probePointCandidate))
        *it = std::move(candidate);
}

inline CriticalPoint makeCriticalPoint(const Node* ast, double k, const Point3& q, const Point3& probe,
                                       double radius, int iterations, bool isProbePoint,
                                       const std::string& source) {
    const DifferentialData d = evaluateJet(ast, q, k);
    const CriticalPointTolerances tolerances = criticalPointTolerances(ast, k, q, probe, radius, d);
    CriticalPoint row;
    row.k = k;
    row.point = q;
    row.fieldValue = d.value;
    row.fieldScale = evalAbsScale(ast, q, k);
    row.fieldTolerance = tolerances.field;
    row.gradientMagnitude = gradientNorm(d.gradient);
    row.gradientTolerance = tolerances.gradient;
    row.newtonStepTolerance = tolerances.newtonStep;
    row.deduplicationTolerance = tolerances.deduplication;
    row.hessianDeterminant = determinant3(d.hessian);
    row.iterations = iterations;
    row.probePointCandidate = isProbePoint;
    row.searchSource = source;
    const bool stationary = row.gradientMagnitude <= row.gradientTolerance;
    const bool onLevel = std::fabs(row.fieldValue) <= row.fieldTolerance;
    if (stationary && onLevel) row.status = "NUMERICAL_SINGULAR_CANDIDATE_ON_ZERO_LEVEL_NOT_CERTIFIED";
    else if (stationary) row.status = "STATIONARY_POINT_OFF_ZERO_LEVEL";
    else row.status = "STATIONARY_RESIDUAL_ABOVE_SCALE_AWARE_TOLERANCE";
    return row;
}

// Gauss-Jordan with partial pivoting; false for a numerically singular matrix.
inline bool solve3x3(const std::array<std::array<double, 3>, 3>& a, const Point3& b, Point3& x) {
    double m[3][4] = {{a[0][0], a[0][1], a[0][2], b.x},
                      {a[1][0], a[1][1], a[1][2], b.y},
                      {a[2][0], a[2][1], a[2][2], b.z}};
    double matrixScale = 0.0;
    for (int r = 0; r < 3; ++r) for (int c = 0; c < 3; ++c) matrixScale = std::max(matrixScale, std::fabs(m[r][c]));
    if (!(matrixScale > 0.0) || !std::isfinite(matrixScale)) return false;
    const double pivotTolerance = 64.0 * std::numeric_limits<double>::epsilon() * matrixScale;
    for (int c = 0; c < 3; ++c) {
        int pivot = c;
        for (int r = c + 1; r < 3; ++r) if (std::fabs(m[r][c]) > std::fabs(m[pivot][c])) pivot = r;
        if (std::fabs(m[pivot][c]) <= pivotTolerance) return false;
        if (pivot != c) for (int j = c; j < 4; ++j) std::swap(m[pivot][j], m[c][j]);
        const double div = m[c][c];
        for (int j = c; j < 4; ++j) m[c][j] /= div;
        for (int r = 0; r < 3; ++r) {
            if (r == c) continue;
            const double f = m[r][c];
            for (int j = c; j < 4; ++j) m[r][j] -= f * m[c][j];
        }
    }
    x = {m[0][3], m[1][3], m[2][3]};
    return std::isfinite(x.x) && std::isfinite(x.y) && std::isfinite(x.z);
}

// Damped Newton on grad F_k = 0 inside the box |x - center|_inf <= radius.
inline bool newtonStationary(const Node* ast, double k, const Point3& center, double radius,
                             Point3 seed, Point3& root, int& iterations, double& gradResidual) {
    Point3 x = seed;
    for (int it = 0; it < 60; ++it) {
        const DifferentialData d = evaluateJet(ast, x, k);
        const double gn = gradientNorm(d.gradient);
        iterations = it;
        if (!std::isfinite(gn) || !d.smooth) return false;
        if (gn <= criticalPointTolerances(ast, k, x, center, radius, d).gradient) {
            root = x;
            gradResidual = gn;
            return true;
        }
        Point3 delta;
        if (!solve3x3(d.hessian, d.gradient, delta)) return false;
        const double dn = norm(delta);
        if (dn > radius * 0.5 && dn > 0.0) delta = delta * (radius * 0.5 / dn);
        bool accepted = false;
        double alpha = 1.0;
        for (int ls = 0; ls < 16; ++ls, alpha *= 0.5) {
            const Point3 candidate = x - delta * alpha;
            if (iso_geometry::maxNormDistance(candidate, center) > radius) continue;
            const DifferentialData trial = evaluateJet(ast, candidate, k);
            const double trialNorm = gradientNorm(trial.gradient);
            if (std::isfinite(trialNorm) && trial.smooth && trialNorm < gn) {
                x = candidate;
                accepted = true;
                break;
            }
        }
        // A small Newton step alone is not accepted as a stationary point: the
        // gradient equation must independently meet its own residual test.
        if (!accepted) return false;
    }
    const DifferentialData d = evaluateJet(ast, x, k);
    gradResidual = gradientNorm(d.gradient);
    root = x;
    return std::isfinite(gradResidual) &&
           gradResidual <= criticalPointTolerances(ast, k, x, center, radius, d).gradient;
}

inline std::vector<CriticalPoint> criticalPointsForK(const Node* ast, double k, const Point3& probe, double radius) {
    std::vector<CriticalPoint> candidates;

    // Always evaluate the selected probe exactly before Newton starts. This
    // prevents an approximate root from a neighbouring seed from displacing
    // an exact stationary probe during root deduplication.
    const DifferentialData probeJet = evaluateJet(ast, probe, k);
    const CriticalPointTolerances probeTol = criticalPointTolerances(ast, k, probe, probe, radius, probeJet);
    if (probeJet.smooth && gradientNorm(probeJet.gradient) <= probeTol.gradient)
        addOrImproveCriticalPoint(candidates,
            makeCriticalPoint(ast, k, probe, probe, radius, 0, true, "DIRECT_PROBE_EVALUATION"));

    // Local 3x3x3 multi-start search. The exact center is intentionally
    // omitted because it was already evaluated above.
    for (int ix = -1; ix <= 1; ++ix) for (int iy = -1; iy <= 1; ++iy) for (int iz = -1; iz <= 1; ++iz) {
        if (ix == 0 && iy == 0 && iz == 0) continue;
        const Point3 seed{probe.x + 0.45 * radius * ix, probe.y + 0.45 * radius * iy, probe.z + 0.45 * radius * iz};
        Point3 root;
        int iterations = 0;
        double residual = std::numeric_limits<double>::infinity();
        if (!newtonStationary(ast, k, probe, radius, seed, root, iterations, residual)) continue;
        addOrImproveCriticalPoint(candidates,
            makeCriticalPoint(ast, k, root, probe, radius, iterations, false, "MULTISTART_DAMPED_NEWTON"));
    }

    if (candidates.empty()) {
        CriticalPoint row;
        row.k = k;
        row.status = "NO_STATIONARY_POINT_FOUND_NOT_CERTIFIED";
        return {row};
    }
    std::sort(candidates.begin(), candidates.end(), [&](const CriticalPoint& a, const CriticalPoint& b) {
        return distance(a.point, probe) < distance(b.point, probe);
    });
    return candidates;
}

// Runs the k values in parallel; rows stay grouped in the given k order.
inline std::vector<CriticalPoint> findCriticalPoints(const Node* ast, const std::vector<double>& kValues,
                                                     const Point3& probe, double radius) {
    std::vector<std::vector<CriticalPoint>> perK(kValues.size());
    iso_parallel::parallelFor(kValues.size(), [&](std::size_t i) {
        perK[i] = criticalPointsForK(ast, kValues[i], probe, radius);
    });
    std::vector<CriticalPoint> out;
    for (auto& rows : perK) out.insert(out.end(), rows.begin(), rows.end());
    return out;
}

inline Point3 surfacePoint(const SurfaceVertex& v) {
    return {static_cast<double>(v.px), static_cast<double>(v.py), static_cast<double>(v.pz)};
}

// Triangles whose centroid lies in the box of half-size `radius` around `center`.
inline std::vector<iso_geometry::Triangle> meshTriangles(const IsosurfaceResult& mesh, const Point3& center,
                                                         double radius) {
    std::vector<iso_geometry::Triangle> out;
    for (std::size_t i = 0; i + 2 < mesh.vertices.size(); i += 3) {
        iso_geometry::Triangle t = iso_geometry::makeTriangle(
            surfacePoint(mesh.vertices[i]), surfacePoint(mesh.vertices[i + 1]), surfacePoint(mesh.vertices[i + 2]));
        if (iso_geometry::maxNormDistance(t.centroid, center) <= radius) out.push_back(t);
    }
    return out;
}

// Barycentric samples (n subdivisions per edge) on triangles near `center`.
inline std::vector<Point3> triangleInteriorSamples(const IsosurfaceResult& mesh, const Point3& center,
                                                   double sampleHalfSize, int subdivisions = 3) {
    std::vector<Point3> out;
    const int n = std::max(1, subdivisions);
    for (std::size_t i = 0; i + 2 < mesh.vertices.size(); i += 3) {
        const Point3 a = surfacePoint(mesh.vertices[i]);
        const Point3 b = surfacePoint(mesh.vertices[i + 1]);
        const Point3 c = surfacePoint(mesh.vertices[i + 2]);
        if (iso_geometry::maxNormDistance((a + b + c) * (1.0 / 3.0), center) > sampleHalfSize + 1e-12) continue;
        for (int u = 0; u <= n; ++u) for (int v = 0; v <= n - u; ++v) {
            const double wa = static_cast<double>(u) / n, wb = static_cast<double>(v) / n, wc = 1.0 - wa - wb;
            const Point3 q = a * wa + b * wb + c * wc;
            if (iso_geometry::maxNormDistance(q, center) <= sampleHalfSize) out.push_back(q);
        }
    }
    return out;
}

struct GeometrySample {
    double k = 0.0;
    double meshToReference = std::numeric_limits<double>::quiet_NaN();
    double referenceToMesh = std::numeric_limits<double>::quiet_NaN();
    double symmetricEstimate = std::numeric_limits<double>::quiet_NaN();
    double meshToReferenceRms = std::numeric_limits<double>::quiet_NaN();
    double referenceToMeshRms = std::numeric_limits<double>::quiet_NaN();
    std::size_t sampleCountPerturbed = 0, sampleCountReference = 0;
    double interiorSamplingHalfSize = 0.0;
    bool comparable = false;
};

// (max, rms) of source-point distances to the target triangle set.
inline std::pair<double, double> directedTriangleEstimate(const std::vector<Point3>& source,
                                                          const iso_geometry::TriangleBvh& target) {
    constexpr double nan = std::numeric_limits<double>::quiet_NaN();
    if (source.empty() || target.empty()) return {nan, nan};
    double maxd = 0.0, sum = 0.0;
    for (const Point3& p : source) {
        const double d = target.nearestDistance(p);
        if (!std::isfinite(d)) return {nan, nan};
        maxd = std::max(maxd, d);
        sum += d * d;
    }
    return {maxd, std::sqrt(sum / static_cast<double>(source.size()))};
}

inline GeometrySample compareMeshes(double k, const IsosurfaceResult& mesh, const IsosurfaceResult& reference,
                                    const Point3& center, double radius, double gridSpacing, int subdivisions = 3) {
    GeometrySample out;
    out.k = k;
    // Erode the requested box by two grid widths. This avoids most of the
    // artificial mismatch caused by abruptly clipping triangle patches at the
    // box boundary. The requested and interior half-sizes are both recorded.
    out.interiorSamplingHalfSize = std::max(0.0, radius - 2.0 * gridSpacing);
    const iso_geometry::TriangleBvh meshBvh(meshTriangles(mesh, center, radius + 2.0 * gridSpacing));
    const iso_geometry::TriangleBvh refBvh(meshTriangles(reference, center, radius + 2.0 * gridSpacing));
    const auto meshSamples = triangleInteriorSamples(mesh, center, out.interiorSamplingHalfSize, subdivisions);
    const auto refSamples = triangleInteriorSamples(reference, center, out.interiorSamplingHalfSize, subdivisions);
    out.sampleCountPerturbed = meshSamples.size();
    out.sampleCountReference = refSamples.size();
    if (meshSamples.empty() || refSamples.empty() || meshBvh.empty() || refBvh.empty()) return out;
    const auto a = directedTriangleEstimate(meshSamples, refBvh);
    const auto b = directedTriangleEstimate(refSamples, meshBvh);
    out.meshToReference = a.first;
    out.meshToReferenceRms = a.second;
    out.referenceToMesh = b.first;
    out.referenceToMeshRms = b.second;
    if (std::isfinite(a.first) && std::isfinite(b.first)) {
        out.symmetricEstimate = std::max(a.first, b.first);
        out.comparable = true;
    }
    return out;
}

struct TaylorRow { std::string name; DifferentialData data; int firstNonzeroOrder = -1; };

struct ResearchAnalysis {
    Point3 probe{};
    double requestedLocalHalfSize = 0.0;
    std::string affineStatus = "NOT_VERIFIED_AFFINE_IN_K";
    int gridResolution = 0;
    double halfExtent = 0.0, gridSpacing = 0.0;
    int triangleSubdivisions = 3;
    std::string fExpression = "UNKNOWN", gExpression = "UNKNOWN";
    TaylorRow fTaylor{"F", {}, -1}, gTaylor{"G", {}, -1};
    IsolationAssessment isolation;
    std::vector<CriticalPoint> criticalPoints;
    std::vector<GeometrySample> geometry;
    std::string summary;
};

// `precomputedGeometry`, when given, must hold one GeometrySample per entry
// of options.kValues (in that order); otherwise the meshes are extracted here.
inline ResearchAnalysis analyze(const Node* familyAst, const iso_convergence::Options& options,
                                int triangleSubdivisions = 3,
                                const std::vector<GeometrySample>* precomputedGeometry = nullptr) {
    if (!familyAst) throw std::invalid_argument("Research analysis needs a parsed family expression.");
    if (options.kValues.empty()) throw std::invalid_argument("Research analysis needs explicit k values.");
    if (precomputedGeometry && precomputedGeometry->size() != options.kValues.size())
        throw std::invalid_argument("Precomputed geometry does not match the k samples.");
    ResearchAnalysis out;
    out.probe = options.probe;
    out.requestedLocalHalfSize = options.localHalfSize;
    out.gridResolution = options.resolution;
    out.halfExtent = options.halfExtent;
    out.gridSpacing = 2.0 * options.halfExtent / options.resolution;
    out.triangleSubdivisions = triangleSubdivisions;

    const iso_convergence::AffineParts affine = iso_convergence::decomposeAffineK(familyAst);
    if (!affine.valid) {
        out.summary = "Taylor F/G analysis unavailable: the expression could not be structurally decomposed as F(x,y,z)+kG(x,y,z). The full k-dependent stationary search remains available.";
    } else {
        out.affineStatus = "STRUCTURALLY_AFFINE_IN_K";
        out.fExpression = affine.f;
        out.gExpression = affine.g;
        const NodePtr fAst = parseExpression(affine.f);
        const NodePtr gAst = parseExpression(affine.g);
        out.fTaylor.data = evaluateJet(fAst.get(), options.probe, 0.0);
        out.fTaylor.firstNonzeroOrder = firstNonzeroJetOrder(out.fTaylor.data.taylor);
        out.gTaylor.data = evaluateJet(gAst.get(), options.probe, 0.0);
        out.gTaylor.firstNonzeroOrder = firstNonzeroJetOrder(out.gTaylor.data.taylor);
        out.isolation = assessReferenceSingularity(fAst.get(), options.probe, options.localHalfSize);
    }
    out.criticalPoints = findCriticalPoints(familyAst, options.kValues, options.probe, options.localHalfSize);

    if (precomputedGeometry) {
        out.geometry = *precomputedGeometry;
    } else {
        const IsosurfaceResult reference = Isosurface::build(familyAst, 0.0, options.resolution, options.halfExtent);
        out.geometry.resize(options.kValues.size());
        iso_parallel::parallelFor(options.kValues.size(), [&](std::size_t i) {
            const double k = options.kValues[i];
            if (k == 0.0) {
                out.geometry[i] = compareMeshes(k, reference, reference, options.probe, options.localHalfSize,
                                                out.gridSpacing, triangleSubdivisions);
                return;
            }
            const IsosurfaceResult mesh = Isosurface::build(familyAst, k, options.resolution, options.halfExtent);
            out.geometry[i] = compareMeshes(k, mesh, reference, options.probe, options.localHalfSize,
                                            out.gridSpacing, triangleSubdivisions);
        });
    }

    std::size_t singularCandidates = 0, stationary = 0;
    for (const CriticalPoint& c : out.criticalPoints) {
        if (c.status.find("NUMERICAL_SINGULAR_CANDIDATE") != std::string::npos) ++singularCandidates;
        if (c.status.find("NO_STATIONARY") != 0) ++stationary;
    }
    std::ostringstream report;
    report << "Local isolation: " << out.isolation.status << ". ";
    if (out.affineStatus == "STRUCTURALLY_AFFINE_IN_K")
        report << "Taylor jet orders at probe: first nonzero F order=" << out.fTaylor.firstNonzeroOrder
               << ", first nonzero G order=" << out.gTaylor.firstNonzeroOrder << ". ";
    report << "Stationary-point search: " << stationary << " candidate rows, " << singularCandidates
           << " numerical on-zero-level singular candidates. Triangle-sampled geometry comparisons: "
           << out.geometry.size() << " k samples. All root-search labels and mesh distances are numerical diagnostics, not certified completeness or exact Hausdorff distances.";
    out.summary = report.str();
    return out;
}

// ----- CSV export ------------------------------------------------------------

inline bool writeTaylorCsv(const std::string& path, const ResearchAnalysis& a, std::string* error = nullptr) {
    using iso_csv::number;
    using iso_csv::quote;
    std::ofstream o(path);
    if (!o) { if (error) *error = "Cannot open " + path; return false; }
    o << "component,expression,first_nonzero_jet_order,smooth_jet,value,gx,gy,gz,hxx,hxy,hxz,hyy,hyz,hzz,"
         "hessian_determinant,eigenvalue_1,eigenvalue_2,eigenvalue_3,";
    static const char* names[20] = {"c000","c100","c010","c001","c200","c110","c101","c020","c011","c002",
                                    "c300","c210","c201","c120","c111","c102","c030","c021","c012","c003"};
    for (const char* n : names) o << n << ',';
    o << "affine_status,reference_isolation_status,probe_x,probe_y,probe_z\n";
    auto row = [&](const TaylorRow& t, const std::string& expr) {
        const DifferentialData& d = t.data;
        const auto eig = symmetricEigenvalues(d.hessian);
        o << quote(t.name) << ',' << quote(expr) << ',' << t.firstNonzeroOrder << ',' << (d.smooth ? 1 : 0) << ','
          << number(d.value) << ',' << number(d.gradient.x) << ',' << number(d.gradient.y) << ',' << number(d.gradient.z) << ','
          << number(d.hessian[0][0]) << ',' << number(d.hessian[0][1]) << ',' << number(d.hessian[0][2]) << ','
          << number(d.hessian[1][1]) << ',' << number(d.hessian[1][2]) << ',' << number(d.hessian[2][2]) << ','
          << number(determinant3(d.hessian)) << ',' << number(eig[0]) << ',' << number(eig[1]) << ',' << number(eig[2]) << ',';
        for (double c : d.taylor) o << number(c) << ',';
        o << quote(a.affineStatus) << ',' << quote(a.isolation.status) << ','
          << number(a.probe.x) << ',' << number(a.probe.y) << ',' << number(a.probe.z) << "\n";
    };
    row(a.fTaylor, a.fExpression);
    row(a.gTaylor, a.gExpression);
    return static_cast<bool>(o);
}

inline bool writeCriticalPointsCsv(const std::string& path, const ResearchAnalysis& a, std::string* error = nullptr) {
    using iso_csv::number;
    using iso_csv::quote;
    std::ofstream o(path);
    if (!o) { if (error) *error = "Cannot open " + path; return false; }
    o << "k,critical_index,x,y,z,Fk_value,absolute_field_scale,field_zero_tolerance,gradient_norm,gradient_tolerance,"
         "newton_step_tolerance,deduplication_tolerance,hessian_determinant,newton_iterations,probe_point_candidate,"
         "search_source,status,search_method,reference_isolation_status,family_F,family_G\n";
    const std::string method = quote("direct probe evaluation plus 26-start damped Newton on grad(F_k)=0; local and uncertified");
    std::map<double, int> index;
    for (const CriticalPoint& c : a.criticalPoints) {
        const int id = index[c.k]++;
        o << number(c.k) << ',' << id << ',' << number(c.point.x) << ',' << number(c.point.y) << ',' << number(c.point.z) << ','
          << number(c.fieldValue) << ',' << number(c.fieldScale) << ',' << number(c.fieldTolerance) << ','
          << number(c.gradientMagnitude) << ',' << number(c.gradientTolerance) << ','
          << number(c.newtonStepTolerance) << ',' << number(c.deduplicationTolerance) << ','
          << number(c.hessianDeterminant) << ',' << c.iterations << ',' << (c.probePointCandidate ? 1 : 0) << ','
          << quote(c.searchSource) << ',' << quote(c.status) << ',' << method << ','
          << quote(a.isolation.status) << ',' << quote(a.fExpression) << ',' << quote(a.gExpression) << "\n";
    }
    return static_cast<bool>(o);
}

inline bool writeTriangleGeometryCsv(const std::string& path, const ResearchAnalysis& a, std::string* error = nullptr) {
    using iso_csv::number;
    using iso_csv::quote;
    std::ofstream o(path);
    if (!o) { if (error) *error = "Cannot open " + path; return false; }
    o << "k,mesh_to_reference_triangle_sampled_max,reference_to_mesh_triangle_sampled_max,symmetric_triangle_sampled_estimate,"
         "mesh_to_reference_rms,reference_to_mesh_rms,perturbed_source_samples,reference_source_samples,interior_sampling_half_size,comparable,"
         "distance_method,local_box_half_size,grid_resolution_cells_per_axis,grid_spacing,domain_half_extent,triangle_subdivisions\n";
    const std::string method = quote("barycentric samples on source triangles; exact point-to-target-triangle distance accelerated by BVH; sampled estimate, not continuous Hausdorff");
    for (const GeometrySample& g : a.geometry)
        o << number(g.k) << ',' << number(g.meshToReference) << ',' << number(g.referenceToMesh) << ','
          << number(g.symmetricEstimate) << ',' << number(g.meshToReferenceRms) << ',' << number(g.referenceToMeshRms) << ','
          << g.sampleCountPerturbed << ',' << g.sampleCountReference << ',' << number(g.interiorSamplingHalfSize) << ','
          << (g.comparable ? 1 : 0) << ',' << method << ',' << number(a.requestedLocalHalfSize) << ','
          << a.gridResolution << ',' << number(a.gridSpacing) << ',' << number(a.halfExtent) << ','
          << a.triangleSubdivisions << "\n";
    return static_cast<bool>(o);
}

} // namespace iso_research
