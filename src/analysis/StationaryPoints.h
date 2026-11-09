#pragma once

// Certified count of the stationary points of H = F + kG near a singular
// point of F (the origin), for polynomials of degree <= 3.
//
// Why a count can be certified. Let F have an isolated singular point at the
// origin with Milnor number mu. For every small k, the stationary points of
// H_k (the solutions of grad H_k = 0, real AND complex) inside a small ball
// around the origin, counted with multiplicity, add up to exactly mu
// (conservation of the Milnor number; Milnor 1968, Greuel-Lossen-Shustin
// "Introduction to Singularities and Deformations", Ch. I). So if a search finds stationary points whose multiplicities
// add up to mu, it has found all of them, and nothing real can be missing.
//
// Method.
//  1. Perturb the gradient equations by a tiny random complex vector,
//     grad H(q) + eps*v = 0. For generic v every solution is simple, so
//     Newton's method converges to each of them quadratically; a solution of
//     multiplicity m of the unperturbed system splits into m simple ones at
//     distance ~eps^(1/2) from it.
//  2. Find the simple solutions by complex Newton from random starts at many
//     scales (from the size of the region down to 1e-6), stopping as soon as
//     mu distinct solutions have been found.
//  3. Group solutions closer than `clusterRadius`: each group is one
//     stationary point of H, its size is the multiplicity. Simple points are
//     polished by Newton on the unperturbed system.
//  4. A point is real when its imaginary parts vanish (below clusterRadius),
//     and singular (a singular point of the surface H = 0) when H also
//     vanishes relative to the size of its terms.
// The certificate is "the multiplicities add up to mu". It is numerical, not
// a proof in exact arithmetic: it relies on the tolerances below being
// meaningful, which holds when distinct stationary points are farther apart
// than clusterRadius.

#include "CubicPolynomial.h"
#include "LocalAlgebra.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <random>
#include <string>
#include <vector>

namespace iso_stationary {

using Complex = std::complex<double>;
using CPoint = std::array<Complex, 3>;
using Point3 = iso_geometry::Point3;

struct Options {
    double regionRadius = 0.3;     // complex norm |q| = sqrt(sum |q_i|^2) of the counted region
    double perturbation = 1e-14;   // eps, relative to the largest coefficient of grad H
    double clusterRadius = 1e-5;   // solutions closer than this are one stationary point
    double onLevelTolerance = 1e-10; // |H| relative to the size of its terms below which H counts as zero
    int startsPerScale = 24;
    int scales = 40;
    std::uint64_t seed = 0x5eed1234abcdULL;
};

struct StationaryPoint {
    CPoint q{};
    int multiplicity = 1;
    bool real = false;
    Point3 point{};                 // real part
    double value = std::numeric_limits<double>::quiet_NaN();      // H at the real point
    double valueScale = std::numeric_limits<double>::quiet_NaN(); // sum of |terms| of H there
    bool onZeroLevel = false;       // real and H = 0: a singular point of the surface
    int hessianCorank = -1;         // real points only
    int morseIndex = -1;            // real points only: number of negative Hessian eigenvalues
};

struct Result {
    int milnorNumber = -1;          // expected total multiplicity
    int totalMultiplicity = 0;      // found
    std::vector<StationaryPoint> points;
    std::string status = "NOT_RUN";
    bool certified() const { return status == "CERTIFIED_COMPLETE"; }
    std::size_t realCount() const { return std::count_if(points.begin(), points.end(), [](const auto& p) { return p.real; }); }
    std::size_t singularCount() const { return std::count_if(points.begin(), points.end(), [](const auto& p) { return p.onZeroLevel; }); }
    int singularMultiplicity() const {
        int m = 0;
        for (const auto& p : points) if (p.onZeroLevel) m += p.multiplicity;
        return m;
    }
};

namespace detail {

inline double cnorm(const CPoint& q) { return std::sqrt(std::norm(q[0]) + std::norm(q[1]) + std::norm(q[2])); }
inline double cdistance(const CPoint& a, const CPoint& b) {
    return std::sqrt(std::norm(a[0] - b[0]) + std::norm(a[1] - b[1]) + std::norm(a[2] - b[2]));
}

// Solves A x = b for a complex 3x3 matrix (partial pivoting).
inline bool solve(std::array<std::array<Complex, 3>, 3> a, std::array<Complex, 3> b, std::array<Complex, 3>& x) {
    for (int c = 0; c < 3; ++c) {
        int p = c;
        for (int r = c + 1; r < 3; ++r) if (std::abs(a[r][c]) > std::abs(a[p][c])) p = r;
        if (std::abs(a[p][c]) < 1e-300) return false;
        std::swap(a[p], a[c]);
        std::swap(b[p], b[c]);
        for (int r = c + 1; r < 3; ++r) {
            const Complex f = a[r][c] / a[c][c];
            for (int j = c; j < 3; ++j) a[r][j] -= f * a[c][j];
            b[r] -= f * b[c];
        }
    }
    for (int r = 2; r >= 0; --r) {
        Complex s = b[r];
        for (int j = r + 1; j < 3; ++j) s -= a[r][j] * x[j];
        x[r] = s / a[r][r];
    }
    return std::isfinite(std::abs(x[0])) && std::isfinite(std::abs(x[1])) && std::isfinite(std::abs(x[2]));
}

inline double residual(const iso_poly::Cubic& h, const CPoint& q, const CPoint& shift) {
    const auto d = h.evaluate<Complex>(q);
    return std::sqrt(std::norm(d.gradient[0] + shift[0]) + std::norm(d.gradient[1] + shift[1]) +
                     std::norm(d.gradient[2] + shift[2]));
}

// Damped Newton for grad H(q) + shift = 0. Returns true on convergence. Near
// a cluster of simple roots the Jacobian is nearly singular, so a small
// residual alone does not mean an accurate root: the iteration continues
// until the Newton step itself is at rounding level.
inline bool newton(const iso_poly::Cubic& h, const CPoint& shift, CPoint& q, double scale, double limit) {
    double r = residual(h, q, shift);
    for (int it = 0; it < 200; ++it) {
        const auto d = h.evaluate<Complex>(q);
        std::array<Complex, 3> step{};
        if (!solve(d.hessian, {d.gradient[0] + shift[0], d.gradient[1] + shift[1], d.gradient[2] + shift[2]}, step))
            return r <= 1e-14 * scale;
        const double stepNorm = cnorm(step);
        if (stepNorm <= 1e-15 * std::max(cnorm(q), 1e-12)) return r <= 1e-14 * scale;
        double t = 1.0;
        bool accepted = false;
        for (int ls = 0; ls < 30; ++ls, t *= 0.5) {
            const CPoint trial{q[0] - t * step[0], q[1] - t * step[1], q[2] - t * step[2]};
            if (cnorm(trial) > limit) continue;
            const double rt = residual(h, trial, shift);
            if (rt < r || (t == 1.0 && rt <= 1e-14 * scale)) { q = trial; r = rt; accepted = true; break; }
        }
        if (!accepted) return r <= 1e-14 * scale;
    }
    return r <= 1e-14 * scale;
}

// Number of zero (|ev| <= 1e-6 max|ev|) and of negative Hessian eigenvalues.
inline std::pair<int, int> corankAndIndex(const std::array<std::array<double, 3>, 3>& m) {
    const auto ev = iso_research::symmetricEigenvalues(m);
    double big = 0.0;
    for (double v : ev) big = std::max(big, std::fabs(v));
    if (big == 0.0) return {3, 0};
    int zero = 0, negative = 0;
    for (double v : ev) {
        if (std::fabs(v) <= 1e-6 * big) ++zero;
        else if (v < 0.0) ++negative;
    }
    return {zero, negative};
}

} // namespace detail

// Stationary points of h near the origin. `milnorNumber` is mu(F, 0) of the
// unperturbed F (iso_algebra::localAlgebra(F, Milnor).dimension).
inline Result findStationaryPoints(const iso_poly::Cubic& h, int milnorNumber, const Options& options = {}) {
    Result out;
    out.milnorNumber = milnorNumber;
    // Scale of grad H: its largest coefficient (degree <= 2 polynomial).
    double scale = 0.0;
    for (std::size_t i = 1; i < 20; ++i) scale = std::max(scale, std::fabs(h.c[i]));
    if (!(scale > 0.0)) { out.status = "ZERO_GRADIENT_POLYNOMIAL"; return out; }

    std::mt19937_64 rng(options.seed);
    std::uniform_real_distribution<double> unit(-1.0, 1.0);
    auto randomComplex = [&] { return Complex(unit(rng), unit(rng)); };
    const double eps = options.perturbation * scale;
    const CPoint shift{eps * randomComplex(), eps * randomComplex(), eps * randomComplex()};
    const double limit = 4.0 * options.regionRadius;

    // Step 1-2: simple solutions of the perturbed system.
    std::vector<CPoint> simple;
    const int wanted = milnorNumber > 0 ? milnorNumber : std::numeric_limits<int>::max();
    auto inRegion = [&](const CPoint& q) { return detail::cnorm(q) <= options.regionRadius; };
    for (int s = 0; s < options.scales && static_cast<int>(simple.size()) < wanted; ++s) {
        const double radius = options.regionRadius * std::pow(1e-6 / options.regionRadius,
                                                              static_cast<double>(s) / std::max(1, options.scales - 1));
        for (int i = 0; i < options.startsPerScale && static_cast<int>(simple.size()) < wanted; ++i) {
            CPoint q{radius * randomComplex(), radius * randomComplex(), radius * randomComplex()};
            if (!detail::newton(h, shift, q, scale, limit) || !inRegion(q)) continue;
            const bool known = std::any_of(simple.begin(), simple.end(), [&](const CPoint& p) {
                return detail::cdistance(p, q) <= 1e-3 * std::sqrt(eps / scale) + 1e-13;
            });
            if (!known) simple.push_back(q);
        }
    }

    // Step 3: group into stationary points of the unperturbed system.
    std::vector<int> group(simple.size(), -1);
    int groups = 0;
    for (std::size_t i = 0; i < simple.size(); ++i) {
        if (group[i] >= 0) continue;
        group[i] = groups;
        for (bool grew = true; grew;) { // single linkage
            grew = false;
            for (std::size_t j = 0; j < simple.size(); ++j) {
                if (group[j] >= 0) continue;
                for (std::size_t l = 0; l < simple.size(); ++l)
                    if (group[l] == groups && detail::cdistance(simple[j], simple[l]) <= options.clusterRadius) {
                        group[j] = groups;
                        grew = true;
                        break;
                    }
            }
        }
        ++groups;
    }
    const CPoint zero{};
    for (int g = 0; g < groups; ++g) {
        StationaryPoint p;
        p.multiplicity = 0;
        for (std::size_t i = 0; i < simple.size(); ++i) {
            if (group[i] != g) continue;
            for (int a = 0; a < 3; ++a) p.q[a] += simple[i][a];
            ++p.multiplicity;
        }
        for (int a = 0; a < 3; ++a) p.q[a] /= static_cast<double>(p.multiplicity);
        // Spread of the group: the averaged position of a degenerate point is
        // only accurate to about this size.
        double spread = 0.0;
        for (std::size_t i = 0; i < simple.size(); ++i)
            if (group[i] == g) spread = std::max(spread, detail::cdistance(simple[i], p.q));
        if (p.multiplicity == 1) { // polish on the unperturbed system
            CPoint q = p.q;
            if (detail::newton(h, zero, q, scale, limit) && detail::cdistance(q, p.q) <= options.clusterRadius) p.q = q;
        }
        p.real = std::fabs(p.q[0].imag()) <= options.clusterRadius && std::fabs(p.q[1].imag()) <= options.clusterRadius &&
                 std::fabs(p.q[2].imag()) <= options.clusterRadius;
        if (p.real) {
            p.point = {p.q[0].real(), p.q[1].real(), p.q[2].real()};
            const auto d = h.evaluate(p.point);
            p.value = d.value;
            // Size of the cancelling terms: at the point itself, and at the
            // resolution limit clusterRadius around it (at a point where every
            // term vanishes, such as a persistent singular point at the
            // origin, only the second is meaningful).
            const double r = options.clusterRadius;
            p.valueScale = std::max(h.termScale(p.point), h.termScale(p.point + Point3{r, r, r}));
            // For a degenerate point the position error (up to 10 x spread)
            // changes H by up to about |Hess| (10 spread)^2 / 2.
            double hessian = 0.0;
            for (const auto& row : d.hessian) for (double v : row) hessian = std::max(hessian, std::fabs(v));
            const double positionError = p.multiplicity > 1 ? 0.5 * hessian * (10.0 * spread) * (10.0 * spread) : 0.0;
            p.onZeroLevel = std::fabs(p.value) <= options.onLevelTolerance * std::max(p.valueScale, 1e-300) + positionError;
            const auto [zero, negative] = detail::corankAndIndex(d.hessian);
            p.hessianCorank = zero;
            p.morseIndex = negative;
        }
        out.totalMultiplicity += p.multiplicity;
        out.points.push_back(p);
    }
    std::sort(out.points.begin(), out.points.end(), [](const StationaryPoint& a, const StationaryPoint& b) {
        return detail::cnorm(a.q) < detail::cnorm(b.q);
    });
    if (milnorNumber < 0) out.status = "MILNOR_NUMBER_UNKNOWN";
    else if (out.totalMultiplicity == milnorNumber) out.status = "CERTIFIED_COMPLETE";
    else if (out.totalMultiplicity < milnorNumber) out.status = "INCOMPLETE_SEARCH";
    else out.status = "EXCESS_IN_REGION";
    return out;
}

} // namespace iso_stationary
