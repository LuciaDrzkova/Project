#pragma once

// Classifies HOW the surfaces S_k = {F + kG = 0} approach the singular surface
// S_0 = {F = 0} near a singular point P as k -> 0, separately for k < 0 and
// k > 0. Input is one iso_research::analyzeFamily() run (all k samples).
//
// Four questions are answered per side:
//   1. Surface:  does d_H(S_k, S_0) -> 0 near P, and at which rate |k|^alpha?
//                d_H is the local Hausdorff distance in the box Q = P + [-a,a]^3,
//                measured on the exact surfaces (ExactDistance.h) when F and G
//                are cubic polynomials, otherwise between the meshes (then
//                fitted only inside the resolvable window h/2 <= d_H <= 0.6a).
//                Every fitted exponent comes with its standard error and a 95%
//                confidence interval.
//   2. Point P:  does P stay on S_k, and is it regular or singular there?
//   3. Singular points of S_k near P: removed, persistent at P, or moving
//                into P like |k|^gamma? With cubic F and G these come from the
//                certified stationary-point search (StationaryPoints.h), which
//                also reports their multiplicities (Milnor numbers).
//   4. Gradient: how fast min |grad F_k| on S_k near P grows with |k|.
//
// The exponent is predicted from the Taylor jets at P (scaling / weighted
// blow-up argument): if the weighted-homogeneous principal part of F has
// weights w (P(t^w x) = t P(x)) and is nondegenerate, and G has weighted order
// d_w(G) < 1, then the inner scale is lambda = |k|^(1/(1 - d_w(G))): the core
// of S_k is S_0's core moved by ~lambda^(w_i) along axis i. Hence
//     alpha in [ w_min, w_max ] / (1 - d_w(G))      (each end capped at 1);
// which end is attained depends on the real geometry of the principal part
// (an indefinite principal part lets S_0 shadow the w_min direction). For a
// nondegenerate homogeneous F of degree d and G of order m < d both ends equal
// the familiar alpha = 1/(d - m). Everything here is a numerical diagnostic.

#include "FamilyAnalysis.h"
#include "core/CsvFormat.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace iso_classify {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// ----- k sampling -------------------------------------------------------------

// Symmetric logarithmic sampling for rate fits: `perSide` values of |k| from
// kMaxAbs down to kMaxAbs * 10^-decades on each side (no k = 0 sample).
inline std::vector<double> makeSymmetricLogKValues(double kMaxAbs, int perSide, double decades) {
    if (!(kMaxAbs > 0.0) || !std::isfinite(kMaxAbs))
        throw std::invalid_argument("|k| maximum must be finite and positive.");
    if (perSide < 3 || perSide > 50)
        throw std::invalid_argument("Logarithmic sampling needs 3 to 50 values per side.");
    if (!(decades > 0.0) || decades > 8.0)
        throw std::invalid_argument("The |k| range must span between 0 and 8 decades.");
    std::vector<double> out;
    for (int sign : {-1, 1})
        for (int i = 0; i < perSide; ++i) {
            const double t = static_cast<double>(i) / static_cast<double>(perSide - 1);
            out.push_back(sign * kMaxAbs * std::pow(10.0, -decades * (sign < 0 ? t : 1.0 - t)));
        }
    return out; // ascending: -kMax ... -kMin, +kMin ... +kMax
}

// ----- Results ------------------------------------------------------------------

// Least-squares line  log(value) = log(coefficient) + exponent * log|k|.
// standardError is the usual standard error of the slope,
//     SE = sqrt( sum r_i^2 / (n - 2) / sum (x_i - mean x)^2 ),
// and [ciLow, ciHigh] = exponent -/+ t_{0.975, n-2} * SE is the 95% confidence
// interval under the regression model (independent, normally distributed
// scatter around the line). It measures how well the samples determine the
// slope of the sampled range, not how far that range is from the asymptotic
// regime k -> 0.
struct RateFit {
    double exponent = kNaN;    // slope of log(value) against log|k|
    double coefficient = kNaN; // value ~ coefficient * |k|^exponent
    double rSquared = kNaN;
    double standardError = kNaN;
    double ciLow = kNaN, ciHigh = kNaN;
    std::size_t samples = 0;
    double kLow = kNaN, kHigh = kNaN; // |k| range used
    bool valid() const { return samples >= 3 && std::isfinite(exponent); }
};

enum class SurfaceClass {
    Identical,          // S_k = S_0 near P (e.g. G proportional to F)
    Linear,             // d_H ~ |k|            (Lipschitz)
    Sublinear,          // d_H ~ |k|^alpha, 0 < alpha < 1   (Hoelder; singular core)
    Superlinear,        // d_H ~ |k|^alpha, alpha > 1
    NotConverging,      // d_H resolved but does not decrease
    NoLocalSurface,     // S_k has no surface near P on this side (one-sided)
    ResolutionLimited,  // whole deformation below the grid spacing
    ExceedsBall,        // deformation larger than the analysis ball for all k
    Insufficient        // too few samples inside the resolvable window
};

enum class PointRelation {
    OffSurface,         // G(P) != 0: P leaves S_k
    OnSurfaceRegular,   // G(P) = 0, grad G(P) != 0: P stays on S_k, becomes regular
    OnSurfaceSingular,  // G(P) = 0, grad G(P) = 0: P stays a singular point of S_k
    NotSingularAtZero   // P is not a singular point of S_0: nothing to unfold
};

enum class SingularBehavior {
    Regularized,        // no singular point of S_k near P (k != 0)
    PersistentAtProbe,  // P itself stays singular
    Moving,             // singular points away from P that approach it as k -> 0
    PersistentAndMoving,
    Intermittent        // found for some k only (inspect the CSV)
};

struct Prediction {
    int orderF = -1;            // first nonzero Taylor order of F at P
    int orderG = -1;            // first nonzero Taylor order of G at P (-1: zero 3-jet)
    bool weightsFound = false;
    bool principalNondegenerate = false;
    std::array<double, 3> weights{kNaN, kNaN, kNaN};
    double weightedOrderG = kNaN;
    double exponent = kNaN;     // predicted alpha for d_H(S_k, S_0): lower end of the range
    double exponentMax = kNaN;  // upper end (equal to `exponent` for homogeneous F)
    std::string regime = "UNDETERMINED";
    std::string principalPart;  // text of the weighted principal part of F
    std::string note;
};

struct SideResult {
    int sign = 0;
    std::size_t samples = 0;
    SurfaceClass surface = SurfaceClass::Insufficient;
    RateFit hausdorff;          // local Hausdorff distance (exact surfaces when available)
    bool exactDistances = false; // hausdorff was measured on the exact surfaces
    RateFit meshHausdorff;      // the same, measured between the two meshes (for comparison)
    double extractionErrorMax = kNaN; // largest mesh-to-exact-surface Hausdorff distance on this side
    double hausdorffAtSmallestK = kNaN;
    RateFit probeDistance;      // d(P, S_k), when P leaves the surface
    std::size_t samplesWithoutLocalSurface = 0;
    SingularBehavior singular = SingularBehavior::Regularized;
    std::size_t kWithSingularPoints = 0;
    std::size_t maxSingularPoints = 0;
    RateFit singularDistance;   // distance of moving singular points to P
    int movingSingularCorank = -1; // Hessian corank at a moving singular point (0: A1, 1: A_k k>=2, 2: D/E, 3: cone)
    int movingSingularMultiplicity = -1; // Milnor number of a moving singular point (certified search)
    std::size_t certifiedK = 0;  // k values whose stationary-point count is certified complete
    RateFit minGradient;        // min |grad F_k| on S_k near P
    int hessianCorankAtP = -1;  // corank of Hess F_k(P), smallest |k| on this side
    bool agreesWithPrediction = false;
    std::string verdict;
};

struct Classification {
    std::string familyExpression, fExpression, gExpression;
    iso_geometry::Point3 probe{};
    double gridSpacing = kNaN, localHalfSize = kNaN;
    double windowLow = kNaN, windowHigh = kNaN; // window of d_H values used in the rate fit
    int milnorNumber = -1, tjurinaNumber = -1;
    std::string tjurinaBasis;    // monomial basis of the Tjurina algebra of F at P
    std::string essentialCoordinates; // class of G in that basis: "coefficient*monomial" terms
    double essentialNorm = kNaN; // 0: G is, to first order, a coordinate change (inessential)
    double gAtProbe = kNaN, gradGAtProbe = kNaN;
    int hessianCorankF = -1;
    PointRelation point = PointRelation::NotSingularAtZero;
    Prediction prediction;
    SideResult negative, positive;
    std::string summary;
};

inline const char* toString(SurfaceClass c) {
    switch (c) {
        case SurfaceClass::Identical: return "IDENTICAL";
        case SurfaceClass::Linear: return "LINEAR";
        case SurfaceClass::Sublinear: return "SUBLINEAR_HOELDER";
        case SurfaceClass::Superlinear: return "SUPERLINEAR";
        case SurfaceClass::NotConverging: return "NOT_CONVERGING";
        case SurfaceClass::NoLocalSurface: return "NO_LOCAL_SURFACE";
        case SurfaceClass::ResolutionLimited: return "RESOLUTION_LIMITED";
        case SurfaceClass::ExceedsBall: return "EXCEEDS_ANALYSIS_BALL";
        case SurfaceClass::Insufficient: return "INSUFFICIENT_WINDOW";
    }
    return "UNKNOWN";
}
inline const char* toString(PointRelation p) {
    switch (p) {
        case PointRelation::OffSurface: return "P_LEAVES_SURFACE";
        case PointRelation::OnSurfaceRegular: return "P_STAYS_REGULAR";
        case PointRelation::OnSurfaceSingular: return "P_STAYS_SINGULAR";
        case PointRelation::NotSingularAtZero: return "P_NOT_SINGULAR_AT_K0";
    }
    return "UNKNOWN";
}
inline const char* toString(SingularBehavior s) {
    switch (s) {
        case SingularBehavior::Regularized: return "REGULARIZED";
        case SingularBehavior::PersistentAtProbe: return "PERSISTENT_AT_P";
        case SingularBehavior::Moving: return "MOVING_INTO_P";
        case SingularBehavior::PersistentAndMoving: return "PERSISTENT_AND_MOVING";
        case SingularBehavior::Intermittent: return "INTERMITTENT";
    }
    return "UNKNOWN";
}

namespace detail {

constexpr double kZero = 1e-12; // relative tolerance for "zero" Taylor data

// 0.975 quantile of Student's t distribution with `dof` degrees of freedom.
inline double studentT975(std::size_t dof) {
    static const double table[] = {0.0, 12.706, 4.303, 3.182, 2.776, 2.571, 2.447, 2.365, 2.306, 2.262, 2.228,
                                   2.201, 2.179, 2.160, 2.145, 2.131, 2.120, 2.110, 2.101, 2.093, 2.086,
                                   2.080, 2.074, 2.069, 2.064, 2.060, 2.056, 2.052, 2.048, 2.045, 2.042};
    if (dof == 0) return kNaN;
    if (dof <= 30) return table[dof];
    return 1.960 + 2.4 / static_cast<double>(dof); // close to the exact value for dof > 30
}

// Least-squares fit of log(value) = log(c) + alpha log|k|.
inline RateFit fitRate(const std::vector<std::pair<double, double>>& kv) {
    RateFit fit;
    std::vector<std::pair<double, double>> logs;
    for (const auto& [k, v] : kv)
        if (std::fabs(k) > 0.0 && std::isfinite(v) && v > 0.0) logs.emplace_back(std::log(std::fabs(k)), std::log(v));
    fit.samples = logs.size();
    if (logs.size() < 3) return fit;
    double mx = 0, my = 0;
    for (const auto& [x, y] : logs) { mx += x; my += y; }
    mx /= logs.size(); my /= logs.size();
    double sxx = 0, sxy = 0, syy = 0;
    for (const auto& [x, y] : logs) { sxx += (x - mx) * (x - mx); sxy += (x - mx) * (y - my); syy += (y - my) * (y - my); }
    // Require at least a factor 2 in |k| so the slope is meaningful.
    double lo = logs.front().first, hi = lo;
    for (const auto& p : logs) { lo = std::min(lo, p.first); hi = std::max(hi, p.first); }
    fit.kLow = std::exp(lo);
    fit.kHigh = std::exp(hi);
    if (sxx <= 0.0 || hi - lo < std::log(2.0)) { fit.samples = 0; return fit; }
    fit.exponent = sxy / sxx;
    fit.coefficient = std::exp(my - fit.exponent * mx);
    fit.rSquared = syy > 0.0 ? (sxy * sxy) / (sxx * syy) : 1.0;
    const std::size_t n = logs.size();
    if (n > 2) {
        double rss = 0.0;
        for (const auto& [x, y] : logs) {
            const double r = y - (my + fit.exponent * (x - mx));
            rss += r * r;
        }
        fit.standardError = std::sqrt(rss / static_cast<double>(n - 2) / sxx);
        const double t = studentT975(n - 2);
        fit.ciLow = fit.exponent - t * fit.standardError;
        fit.ciHigh = fit.exponent + t * fit.standardError;
    }
    return fit;
}

// Taylor coefficient arrays c[20] (see iso_research::exponents()).
inline double jetScale(const std::array<double, 20>& c) {
    double s = 0.0;
    for (double v : c) s = std::max(s, std::fabs(v));
    return s;
}

inline std::string monomialText(const iso_research::Exponent3& e) {
    std::string s;
    auto add = [&](const char* name, int p) {
        if (p == 0) return;
        if (!s.empty()) s += '*';
        s += name;
        if (p > 1) s += '^' + std::to_string(p);
    };
    add("x", e.x); add("y", e.y); add("z", e.z);
    return s.empty() ? "1" : s;
}

inline double weightedDegree(const iso_research::Exponent3& e, const std::array<double, 3>& w) {
    return e.x * w[0] + e.y * w[1] + e.z * w[2];
}

// Gradient of the polynomial sum_{i in terms} c_i x^{e_i}.
inline std::array<double, 3> gradient(const std::vector<std::pair<iso_research::Exponent3, double>>& terms,
                                      const std::array<double, 3>& x) {
    std::array<double, 3> g{0, 0, 0};
    for (const auto& [e, c] : terms) {
        const int p[3] = {e.x, e.y, e.z};
        for (int i = 0; i < 3; ++i) {
            if (p[i] == 0) continue;
            double v = c * p[i];
            for (int j = 0; j < 3; ++j) v *= std::pow(x[j], j == i ? p[j] - 1 : p[j]);
            g[i] += v;
        }
    }
    return g;
}

// Numerical check that a weighted-homogeneous polynomial has no real critical
// point other than 0: min |grad| on the weighted unit sphere (sampled, then
// locally refined) must stay a clear fraction of its maximum.
inline bool weightedNondegenerate(const std::vector<std::pair<iso_research::Exponent3, double>>& terms,
                                  const std::array<double, 3>& w) {
    auto onSphere = [&](std::array<double, 3> u) {
        const double n = std::sqrt(u[0] * u[0] + u[1] * u[1] + u[2] * u[2]);
        for (double& v : u) v /= n;
        double s = 0.0;
        for (int i = 0; i < 3; ++i) s += std::pow(std::fabs(u[i]), 2.0 / w[i]);
        const double t = 1.0 / std::sqrt(s);
        return std::array<double, 3>{std::pow(t, w[0]) * u[0], std::pow(t, w[1]) * u[1], std::pow(t, w[2]) * u[2]};
    };
    auto gradNorm = [&](const std::array<double, 3>& u) {
        const auto g = gradient(terms, onSphere(u));
        return std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
    };
    constexpr int kSamples = 4000;
    std::vector<std::pair<double, std::array<double, 3>>> values;
    values.reserve(kSamples);
    double maxNorm = 0.0;
    const double golden = 3.14159265358979323846 * (3.0 - std::sqrt(5.0));
    for (int i = 0; i < kSamples; ++i) {
        const double z = 1.0 - 2.0 * (i + 0.5) / kSamples, r = std::sqrt(1.0 - z * z);
        const std::array<double, 3> u{r * std::cos(golden * i), r * std::sin(golden * i), z};
        const double g = gradNorm(u);
        maxNorm = std::max(maxNorm, g);
        values.emplace_back(g, u);
    }
    if (!(maxNorm > 0.0)) return false;
    std::partial_sort(values.begin(), values.begin() + 24, values.end(),
                      [](const auto& a, const auto& b) { return a.first < b.first; });
    double best = values.front().first;
    for (int s = 0; s < 24; ++s) { // coordinate-wise pattern search from the 24 smallest samples
        auto [g, u] = values[static_cast<std::size_t>(s)];
        for (double step = 0.05; step > 1e-7; step *= 0.5)
            for (int improved = 1; improved;) {
                improved = 0;
                for (int i = 0; i < 3; ++i)
                    for (double d : {step, -step}) {
                        auto v = u;
                        v[i] += d;
                        const double gv = gradNorm(v);
                        if (gv < g) { g = gv; u = v; improved = 1; }
                    }
            }
        best = std::min(best, g);
    }
    return best > 1e-4 * maxNorm;
}

// Weighted principal part of F from its 3-jet: weights w > 0 such that every
// monomial has weighted degree >= 1 and the degree-1 monomials (the principal
// part) involve all three variables. Chooses a nondegenerate one if possible.
inline Prediction predict(const iso_research::ResearchAnalysis& r) {
    Prediction p;
    const auto& e = iso_research::exponents();
    const auto& cf = r.fTaylor.data.taylor;
    const auto& cg = r.gTaylor.data.taylor;
    p.orderF = r.fTaylor.firstNonzeroOrder;
    p.orderG = r.gTaylor.firstNonzeroOrder;
    if (r.affineStatus != "STRUCTURALLY_AFFINE_IN_K") {
        p.note = "family is not verified to be F + kG";
        return p;
    }
    if (p.orderF < 2) {
        p.note = "P is not a singular point of S_0 (F has a nonzero value or gradient at P)";
        return p;
    }
    const double sf = jetScale(cf);
    std::vector<std::pair<iso_research::Exponent3, double>> terms; // order >= 2 terms of F
    for (std::size_t i = 4; i < 20; ++i)
        if (std::fabs(cf[i]) > kZero * sf) terms.emplace_back(e[i], cf[i]);

    std::vector<std::array<double, 3>> candidates;
    for (std::size_t a = 0; a < terms.size(); ++a)
        for (std::size_t b = a + 1; b < terms.size(); ++b)
            for (std::size_t c = b + 1; c < terms.size(); ++c) {
                const std::array<std::array<double, 3>, 3> m{{
                    {double(terms[a].first.x), double(terms[a].first.y), double(terms[a].first.z)},
                    {double(terms[b].first.x), double(terms[b].first.y), double(terms[b].first.z)},
                    {double(terms[c].first.x), double(terms[c].first.y), double(terms[c].first.z)}}};
                iso_geometry::Point3 w;
                if (!iso_research::solve3x3(m, {1, 1, 1}, w)) continue;
                const std::array<double, 3> ww{w.x, w.y, w.z};
                if (ww[0] <= 1e-9 || ww[1] <= 1e-9 || ww[2] <= 1e-9) continue;
                bool ok = true;
                for (const auto& t : terms) ok = ok && weightedDegree(t.first, ww) >= 1.0 - 1e-9;
                if (!ok) continue;
                bool duplicate = false;
                for (const auto& q : candidates)
                    duplicate = duplicate || (std::fabs(q[0] - ww[0]) + std::fabs(q[1] - ww[1]) + std::fabs(q[2] - ww[2]) < 1e-9);
                if (!duplicate) candidates.push_back(ww);
            }

    for (const auto& w : candidates) {
        std::vector<std::pair<iso_research::Exponent3, double>> principal;
        for (const auto& t : terms)
            if (std::fabs(weightedDegree(t.first, w) - 1.0) < 1e-9) principal.push_back(t);
        const bool nondegenerate = weightedNondegenerate(principal, w);
        if (!p.weightsFound || (nondegenerate && !p.principalNondegenerate)) {
            p.weightsFound = true;
            p.weights = w;
            p.principalNondegenerate = nondegenerate;
            std::ostringstream text;
            for (const auto& [mono, coeff] : principal)
                text << (coeff < 0 ? " - " : (text.tellp() > 0 ? " + " : "")) << iso_csv::shortNumber(std::fabs(coeff))
                     << '*' << monomialText(mono);
            p.principalPart = text.str();
        }
        if (p.principalNondegenerate) break;
    }
    if (!p.weightsFound) {
        p.note = "no weighted principal part within the 3-jet (degenerate or higher-order singularity)";
        return p;
    }

    const double sg = std::max(jetScale(cg), 1e-300);
    p.weightedOrderG = std::numeric_limits<double>::infinity();
    for (std::size_t i = 0; i < 20; ++i)
        if (std::fabs(cg[i]) > kZero * sg) p.weightedOrderG = std::min(p.weightedOrderG, weightedDegree(e[i], p.weights));
    const double wMin = std::min({p.weights[0], p.weights[1], p.weights[2]});
    const double wMax = std::max({p.weights[0], p.weights[1], p.weights[2]});
    p.exponentMax = 1.0;
    if (p.weightedOrderG < 1.0 - 1e-9) {
        p.exponent = std::min(1.0, wMin / (1.0 - p.weightedOrderG));
        p.exponentMax = std::min(1.0, wMax / (1.0 - p.weightedOrderG));
        p.regime = "SCALING";
        p.note = "inner scale lambda = |k|^(1/(1-d_w(G))); rescaled S_k -> model {P(u) + sign(k) G_w(u) = 0}";
    } else if (p.weightedOrderG < 1.0 + 1e-9) {
        p.exponent = 1.0;
        p.regime = "SAME_WEIGHT";
        p.note = "G deforms the principal part itself (equisingular-type): d_H ~ |k| on a fixed ball";
    } else {
        p.exponent = 1.0;
        p.regime = "HIGHER_WEIGHT";
        p.note = "principal part unchanged; only the outer surface moves: d_H ~ |k| on a fixed ball";
    }
    if (!p.principalNondegenerate)
        p.note += " [principal part numerically degenerate: prediction unreliable]";
    return p;
}

inline int corank(const std::array<std::array<double, 3>, 3>& h, double tolerance) {
    const auto ev = iso_research::symmetricEigenvalues(h);
    int n = 0;
    for (double v : ev) n += std::fabs(v) <= tolerance ? 1 : 0;
    return n;
}

inline std::string predictionText(const Prediction& p) {
    if (!std::isfinite(p.exponent)) return "n/a";
    auto r = [](double v) { return iso_csv::shortNumber(std::round(v * 1000.0) / 1000.0); };
    if (std::fabs(p.exponentMax - p.exponent) < 1e-9) return "|k|^" + r(p.exponent);
    return "|k|^alpha, alpha in [" + r(p.exponent) + ", " + r(p.exponentMax) + "]";
}

inline std::string rateText(const RateFit& f) {
    if (!f.valid()) return "n/a";
    std::ostringstream s;
    s << "|k|^" << iso_csv::shortNumber(std::round(f.exponent * 100.0) / 100.0);
    return s.str();
}

inline SideResult classifySide(int sign, const Classification& c, const iso_research::FamilyAnalysis& a,
                               const Node* familyAst) {
    SideResult side;
    side.sign = sign;
    const auto& samples = a.convergence.samples;
    const auto& requested = a.research.geometry; // one row per requested k (same values)
    auto onSide = [sign](double k) { return k != 0.0 && (k > 0.0) == (sign > 0); };
    auto exactFor = [&](double k) -> const iso_exact::Measurement* {
        for (const auto& m : a.exact) if (m.k == k) return &m;
        return nullptr;
    };
    auto stationaryFor = [&](double k) -> const iso_stationary::Result* {
        for (std::size_t i = 0; i < a.exact.size() && i < a.stationary.size(); ++i)
            if (a.exact[i].k == k) return &a.stationary[i];
        return nullptr;
    };
    const bool exact = a.local.available; // exact surfaces and certified stationary points

    std::vector<std::pair<double, double>> hausdorff, inWindow, meshInWindow, probe, gradient;
    std::vector<double> sideK;
    double smallest = std::numeric_limits<double>::infinity();
    const double meshWindowLow = 0.5 * c.gridSpacing;
    for (const auto& s : samples) {
        if (!onSide(s.k)) continue;
        ++side.samples;
        sideK.push_back(s.k);
        if (!s.hasSurfaceInLocalBox) ++side.samplesWithoutLocalSurface;
        double meshValue = kNaN;
        for (const auto& g : requested)
            if (g.k == s.k && g.comparable) meshValue = g.symmetricEstimate;
        if (!std::isfinite(meshValue) && s.geometryComparable) meshValue = s.symmetricHausdorff;
        if (std::isfinite(meshValue) && meshValue >= meshWindowLow && meshValue <= c.windowHigh)
            meshInWindow.emplace_back(s.k, meshValue);
        const iso_exact::Measurement* m = exactFor(s.k);
        if (m && std::isfinite(m->meshHausdorff))
            side.extractionErrorMax = std::isfinite(side.extractionErrorMax) ? std::max(side.extractionErrorMax, m->meshHausdorff)
                                                                             : m->meshHausdorff;
        const double h = exact ? (m && m->comparable ? m->symmetric : kNaN) : meshValue;
        if (std::isfinite(h)) {
            hausdorff.emplace_back(s.k, h);
            if (h >= c.windowLow && h <= c.windowHigh) inWindow.emplace_back(s.k, h);
            if (std::fabs(s.k) < smallest) { smallest = std::fabs(s.k); side.hausdorffAtSmallestK = h; }
        }
        const double pd = exact && m ? m->probeDistance : s.distanceProbeToMesh;
        if (s.hasSurfaceInLocalBox && std::isfinite(pd) && pd >= c.windowLow && pd <= c.windowHigh)
            probe.emplace_back(s.k, pd);
        if (std::isfinite(s.minimumGradientLocal)) gradient.emplace_back(s.k, s.minimumGradientLocal);
    }
    side.exactDistances = exact;
    side.hausdorff = fitRate(inWindow);
    side.meshHausdorff = fitRate(meshInWindow);
    side.probeDistance = fitRate(probe);
    side.minGradient = fitRate(gradient);

    // Surface class.
    const bool allZero = !hausdorff.empty() &&
        std::all_of(hausdorff.begin(), hausdorff.end(), [](const auto& kv) { return kv.second <= 1e-9; });
    if (side.samples == 0) side.surface = SurfaceClass::Insufficient;
    else if (side.samplesWithoutLocalSurface * 2 > side.samples && hausdorff.size() * 2 < side.samples)
        side.surface = SurfaceClass::NoLocalSurface;
    else if (allZero) side.surface = SurfaceClass::Identical;
    else if (!side.hausdorff.valid()) {
        const bool allBelow = std::all_of(hausdorff.begin(), hausdorff.end(), [&](const auto& kv) { return kv.second < c.windowLow; });
        const bool allAbove = std::all_of(hausdorff.begin(), hausdorff.end(), [&](const auto& kv) { return kv.second > c.windowHigh; });
        side.surface = hausdorff.empty() ? SurfaceClass::NoLocalSurface
                     : allBelow ? SurfaceClass::ResolutionLimited
                     : allAbove ? SurfaceClass::ExceedsBall : SurfaceClass::Insufficient;
    } else {
        const double alpha = side.hausdorff.exponent;
        side.surface = alpha <= 0.1 ? SurfaceClass::NotConverging
                     : alpha < 0.9 ? SurfaceClass::Sublinear
                     : alpha <= 1.15 ? SurfaceClass::Linear : SurfaceClass::Superlinear;
        side.agreesWithPrediction = std::isfinite(c.prediction.exponent) &&
                                    alpha >= c.prediction.exponent - 0.1 && alpha <= c.prediction.exponentMax + 0.1;
    }

    // Singular points of S_k near P. Certified search: exact points with
    // multiplicities. Otherwise the real Newton search of the research
    // analysis, where one degenerate point can be reported several times
    // ~1e-8 apart: candidates closer than 1e-3 grid cells are merged.
    const double cluster = exact ? 1e-7 : 1e-3 * c.gridSpacing;
    std::size_t atProbe = 0, moving = 0;
    std::vector<std::pair<double, double>> movingDistance;
    iso_geometry::Point3 movingExample{};
    double movingExampleK = 0.0;
    for (double k : sideK) {
        std::vector<iso_geometry::Point3> distinct;
        std::vector<int> multiplicity;
        std::vector<int> corankOf;
        if (exact) {
            if (const iso_stationary::Result* r = stationaryFor(k)) {
                if (r->certified()) ++side.certifiedK;
                for (const auto& p : r->points)
                    if (p.onZeroLevel && iso_geometry::maxNormDistance(p.point, c.probe) <= c.localHalfSize) {
                        distinct.push_back(p.point);
                        multiplicity.push_back(p.multiplicity);
                        corankOf.push_back(p.hessianCorank);
                    }
            }
        } else {
            for (const auto& cp : a.research.criticalPoints) {
                if (cp.k != k || cp.status.find("NUMERICAL_SINGULAR_CANDIDATE") == std::string::npos) continue;
                if (std::any_of(distinct.begin(), distinct.end(),
                                [&](const iso_geometry::Point3& q) { return iso_geometry::distance(q, cp.point) <= cluster; }))
                    continue;
                distinct.push_back(cp.point);
                multiplicity.push_back(-1);
                corankOf.push_back(-1);
            }
        }
        bool probeHere = false;
        double nearestOther = std::numeric_limits<double>::infinity();
        std::size_t nearestIndex = 0;
        for (std::size_t i = 0; i < distinct.size(); ++i) {
            const double d = iso_geometry::distance(distinct[i], c.probe);
            if (d <= cluster) probeHere = true;
            else if (d < nearestOther) { nearestOther = d; nearestIndex = i; }
        }
        if (!distinct.empty()) ++side.kWithSingularPoints;
        side.maxSingularPoints = std::max(side.maxSingularPoints, distinct.size());
        if (probeHere) ++atProbe;
        if (std::isfinite(nearestOther)) {
            ++moving;
            movingDistance.emplace_back(k, nearestOther);
            if (movingExampleK == 0.0 || std::fabs(k) < std::fabs(movingExampleK)) {
                movingExampleK = k;
                movingExample = distinct[nearestIndex];
                side.movingSingularMultiplicity = multiplicity[nearestIndex];
                side.movingSingularCorank = corankOf[nearestIndex];
            }
        }
    }
    if (!exact && movingExampleK != 0.0 && familyAst) {
        const auto d = iso_research::evaluateJet(familyAst, movingExample, movingExampleK);
        // Newton resolves a flat direction only to ~sqrt(tolerance): judge "zero"
        // eigenvalues relative to the largest one with a matching tolerance.
        side.movingSingularCorank = corank(d.hessian, 1e-4 * std::max(iso_research::hessianMaxAbs(d.hessian), 1e-300));
    }
    side.singularDistance = fitRate(movingDistance);
    const std::size_t n = sideK.size();
    if (side.kWithSingularPoints == 0) side.singular = SingularBehavior::Regularized;
    else if (atProbe == n && moving == 0) side.singular = SingularBehavior::PersistentAtProbe;
    else if (atProbe == n && moving > 0) side.singular = SingularBehavior::PersistentAndMoving;
    else if (moving * 2 >= n && atProbe == 0) side.singular = SingularBehavior::Moving;
    else side.singular = SingularBehavior::Intermittent;

    // Hessian corank of F_k at P for the smallest |k| of this side.
    if (!sideK.empty() && familyAst) {
        const double k = *std::min_element(sideK.begin(), sideK.end(),
                                           [](double x, double y) { return std::fabs(x) < std::fabs(y); });
        const auto d = iso_research::evaluateJet(familyAst, c.probe, k);
        const double scale = std::max(iso_research::hessianMaxAbs(d.hessian), 1e-300);
        side.hessianCorankAtP = corank(d.hessian, 1e-9 * std::max(scale, std::fabs(k)));
    }

    // One-sentence verdict.
    std::ostringstream v;
    v << (sign < 0 ? "k<0: " : "k>0: ");
    switch (side.surface) {
        case SurfaceClass::Linear: case SurfaceClass::Sublinear: case SurfaceClass::Superlinear:
            v << "S_k -> S_0 near P like " << rateText(side.hausdorff);
            if (std::isfinite(side.hausdorff.ciLow))
                v << " (95% CI " << iso_csv::shortNumber(std::round(side.hausdorff.ciLow * 1000.0) / 1000.0) << ".."
                  << iso_csv::shortNumber(std::round(side.hausdorff.ciHigh * 1000.0) / 1000.0) << ")";
            v << " (" << toString(side.surface) << (side.exactDistances ? ", exact surfaces" : ", meshes") << ")";
            if (std::isfinite(c.prediction.exponent))
                v << (side.agreesWithPrediction ? ", within" : ", outside") << " the predicted "
                  << predictionText(c.prediction);
            break;
        case SurfaceClass::Identical: v << "S_k coincides with S_0 near P (trivial deformation)"; break;
        case SurfaceClass::NotConverging: v << "S_k does NOT approach S_0 near P over the resolved range"; break;
        case SurfaceClass::NoLocalSurface: v << "S_k has no surface near P (one-sided: convergence fails from this side)"; break;
        case SurfaceClass::ResolutionLimited: v << "deformation is below the grid spacing: refine the grid or use larger |k|"; break;
        case SurfaceClass::ExceedsBall: v << "deformation exceeds the analysis ball: use smaller |k| or a larger domain"; break;
        case SurfaceClass::Insufficient: v << "too few k inside the measurement window to fit a rate"; break;
    }
    v << "; singular points: " << toString(side.singular);
    if (side.singular == SingularBehavior::Moving || side.singular == SingularBehavior::PersistentAndMoving)
        v << " (" << side.maxSingularPoints << " max, distance " << rateText(side.singularDistance)
          << ", Hessian corank " << side.movingSingularCorank
          << (side.movingSingularMultiplicity > 0 ? ", Milnor number " + std::to_string(side.movingSingularMultiplicity) : std::string())
          << ")";
    if (side.exactDistances)
        v << " Stationary points certified complete for " << side.certifiedK << " of " << side.samples << " k.";
    v << '.';
    side.verdict = v.str();
    return side;
}

} // namespace detail

// Classifies one finished analyzeFamily() run. `familyAst` is the parsed F + kG.
inline Classification classify(const Node* familyAst, const iso_research::FamilyAnalysis& a) {
    Classification c;
    const auto& m = a.convergence.metadata;
    c.familyExpression = m.familyExpression;
    c.fExpression = m.fExpression;
    c.gExpression = m.gExpression;
    c.probe = a.convergence.options.probe;
    c.gridSpacing = a.convergence.gridSpacing;
    c.localHalfSize = a.convergence.options.localHalfSize;
    // Exact distances are meaningful down to rounding level; mesh distances
    // only above half a grid cell.
    c.windowLow = a.local.available ? 1e-12 : 0.5 * c.gridSpacing;
    c.windowHigh = 0.6 * c.localHalfSize;
    if (a.local.available) {
        c.milnorNumber = a.local.milnor.dimension;
        c.tjurinaNumber = a.local.tjurina.dimension;
        c.tjurinaBasis = a.local.tjurina.basisText();
        c.essentialNorm = a.local.essentialNorm;
        std::ostringstream e;
        for (std::size_t i = 0; i < a.local.essential.size(); ++i) {
            const double v = a.local.essential[i];
            if (std::fabs(v) <= 1e-12 * std::max(1.0, a.local.essentialNorm)) continue;
            e << (e.tellp() > 0 ? (v < 0 ? " - " : " + ") : (v < 0 ? "-" : ""))
              << (std::fabs(std::fabs(v) - 1.0) < 1e-12 ? std::string() : iso_csv::shortNumber(std::fabs(v)) + "*")
              << a.local.tjurina.basis[i].text();
        }
        c.essentialCoordinates = e.tellp() > 0 ? e.str() : "0";
    }

    const auto& r = a.research;
    c.prediction = detail::predict(r);
    c.gAtProbe = r.gTaylor.data.value;
    c.gradGAtProbe = iso_research::gradientNorm(r.gTaylor.data.gradient);
    const double fScale = std::max(detail::jetScale(r.fTaylor.data.taylor), 1e-300);
    const double gScale = std::max(detail::jetScale(r.gTaylor.data.taylor), 1e-300);
    c.hessianCorankF = detail::corank(r.fTaylor.data.hessian, 1e-9 * fScale);
    if (c.prediction.orderF < 2) c.point = PointRelation::NotSingularAtZero;
    else if (std::fabs(c.gAtProbe) > detail::kZero * gScale) c.point = PointRelation::OffSurface;
    else if (c.gradGAtProbe > detail::kZero * gScale) c.point = PointRelation::OnSurfaceRegular;
    else c.point = PointRelation::OnSurfaceSingular;

    c.negative = detail::classifySide(-1, c, a, familyAst);
    c.positive = detail::classifySide(+1, c, a, familyAst);

    std::ostringstream s;
    s << "F order " << c.prediction.orderF << " (Hessian corank " << c.hessianCorankF << "), G order "
      << c.prediction.orderG << "; " << toString(c.point);
    if (c.point == PointRelation::OnSurfaceRegular)
        s << " with |grad F_k(P)| = |k|*" << iso_csv::shortNumber(c.gradGAtProbe) << " (linear gradient shift)";
    s << ". Prediction: " << c.prediction.regime;
    if (std::isfinite(c.prediction.exponent))
        s << ", d_H ~ " << detail::predictionText(c.prediction);
    if (c.milnorNumber >= 0)
        s << ". Milnor number " << c.milnorNumber << ", Tjurina number " << c.tjurinaNumber
          << "; essential part of G: " << c.essentialCoordinates;
    s << ".\n" << c.negative.verdict << '\n' << c.positive.verdict;
    c.summary = s.str();
    return c;
}

// ----- Export -------------------------------------------------------------------

inline std::string csvHeader() {
    return "family,F,G,probe_x,probe_y,probe_z,grid_spacing,local_box_half_size,window_low,window_high,"
           "order_F,order_G,hessian_corank_F,G_at_P,grad_G_at_P,point_relation,"
           "weights_x,weights_y,weights_z,principal_part,principal_nondegenerate,weighted_order_G,"
           "predicted_regime,predicted_exponent_min,predicted_exponent_max,prediction_note,"
           "side,samples,surface_class,hausdorff_exponent,hausdorff_coefficient,hausdorff_r2,"
           "hausdorff_fit_samples,fit_k_low,fit_k_high,hausdorff_at_smallest_k,agrees_with_prediction,"
           "probe_distance_exponent,samples_without_local_surface,singular_behavior,k_with_singular_points,"
           "max_singular_points,singular_distance_exponent,moving_singular_hessian_corank,min_gradient_exponent,"
           "hessian_corank_Fk_at_P,distances_on_exact_surfaces,hausdorff_exponent_std_error,hausdorff_exponent_ci95_low,"
           "hausdorff_exponent_ci95_high,mesh_hausdorff_exponent,mesh_hausdorff_exponent_ci95_low,"
           "mesh_hausdorff_exponent_ci95_high,extraction_error_max,singular_distance_exponent_ci95_low,"
           "singular_distance_exponent_ci95_high,moving_singular_milnor_number,certified_k_count,milnor_number,"
           "tjurina_number,tjurina_basis,G_essential_coordinates,G_essential_norm,verdict";
}

inline void writeCsvRows(std::ostream& out, const Classification& c) {
    using iso_csv::number;
    using iso_csv::quote;
    for (const SideResult* s : {&c.negative, &c.positive}) {
        out << quote(c.familyExpression) << ',' << quote(c.fExpression) << ',' << quote(c.gExpression) << ','
            << number(c.probe.x) << ',' << number(c.probe.y) << ',' << number(c.probe.z) << ','
            << number(c.gridSpacing) << ',' << number(c.localHalfSize) << ',' << number(c.windowLow) << ','
            << number(c.windowHigh) << ',' << c.prediction.orderF << ',' << c.prediction.orderG << ','
            << c.hessianCorankF << ',' << number(c.gAtProbe) << ',' << number(c.gradGAtProbe) << ','
            << toString(c.point) << ',' << number(c.prediction.weights[0]) << ',' << number(c.prediction.weights[1])
            << ',' << number(c.prediction.weights[2]) << ',' << quote(c.prediction.principalPart) << ','
            << (c.prediction.principalNondegenerate ? 1 : 0) << ',' << number(c.prediction.weightedOrderG) << ','
            << c.prediction.regime << ',' << number(c.prediction.exponent) << ','
            << number(c.prediction.exponentMax) << ',' << quote(c.prediction.note) << ','
            << (s->sign < 0 ? "k<0" : "k>0") << ',' << s->samples << ',' << toString(s->surface) << ','
            << number(s->hausdorff.exponent) << ',' << number(s->hausdorff.coefficient) << ','
            << number(s->hausdorff.rSquared) << ',' << s->hausdorff.samples << ',' << number(s->hausdorff.kLow) << ','
            << number(s->hausdorff.kHigh) << ',' << number(s->hausdorffAtSmallestK) << ','
            << (s->agreesWithPrediction ? 1 : 0) << ',' << number(s->probeDistance.exponent) << ','
            << s->samplesWithoutLocalSurface << ',' << toString(s->singular) << ',' << s->kWithSingularPoints << ','
            << s->maxSingularPoints << ',' << number(s->singularDistance.exponent) << ','
            << s->movingSingularCorank << ',' << number(s->minGradient.exponent) << ',' << s->hessianCorankAtP << ','
            << (s->exactDistances ? 1 : 0) << ',' << number(s->hausdorff.standardError) << ','
            << number(s->hausdorff.ciLow) << ',' << number(s->hausdorff.ciHigh) << ','
            << number(s->meshHausdorff.exponent) << ',' << number(s->meshHausdorff.ciLow) << ','
            << number(s->meshHausdorff.ciHigh) << ',' << number(s->extractionErrorMax) << ','
            << number(s->singularDistance.ciLow) << ',' << number(s->singularDistance.ciHigh) << ','
            << s->movingSingularMultiplicity << ',' << s->certifiedK << ',' << c.milnorNumber << ','
            << c.tjurinaNumber << ',' << quote(c.tjurinaBasis) << ',' << quote(c.essentialCoordinates) << ','
            << number(c.essentialNorm) << ',' << quote(s->verdict) << '\n';
    }
}

inline bool writeCsv(const std::string& path, const Classification& c, std::string* error = nullptr) {
    std::ofstream out(path);
    if (!out) { if (error) *error = "Cannot open " + path; return false; }
    out << csvHeader() << '\n';
    writeCsvRows(out, c);
    return static_cast<bool>(out);
}

} // namespace iso_classify
