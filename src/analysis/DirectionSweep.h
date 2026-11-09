#pragma once

// Sweep over the directions of a linear perturbation G = c . (p - P), |c| = 1.
//
// For a fixed k, the family F + kG depends on the unit vector c only. For
// every c on a grid of the sphere the certified stationary-point search gives
// the real stationary points of H_k near P, and for each of them its Morse
// index (number of negative Hessian eigenvalues) and the sign of H_k there.
// By Morse theory (Milnor, "Morse Theory", 1963) these data determine how the
// level sets of H_k change near P, so within a region of the sphere where the
// counts by (index, sign) - the "signature" - are constant, the surfaces S_k
// have the same local structure. (Counting signs alone is not enough: for a
// homogeneous cubic F and linear G, H(-p) = -H(p), so the stationary points
// come in pairs with opposite values and opposite indices i and 3 - i, and at
// a crossing both values change sign together.) Two kinds of boundaries
// separate the regions:
//   * the number of real stationary points changes: two of them meet and
//     become complex (a property of the function H_k, not of the surface);
//   * the number of real stationary points is the same but the sign of H_k
//     at one of them changes: there H_k = 0 and grad H_k = 0, so S_k has a
//     singular point. These directions form the discriminant: the directions
//     c for which the line F + kG crosses the set of singular surfaces.
// Boundaries of the second kind are located to machine precision by bisection
// along the great circle between two neighbouring grid directions.
//
// The grid is the Lambert cylindrical equal-area grid: longitude phi in
// [-180, 180) and height z = c_z uniform in [-1, 1], so every cell has the
// same area on the sphere.

#include "LocalAlgebra.h"
#include "StationaryPoints.h"
#include "core/CsvFormat.h"
#include "core/Parallel.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>

namespace iso_sweep {

using Direction = std::array<double, 3>;

struct Signature {
    int real = 0, positive = 0, negative = 0, zero = 0;
    std::array<int, 8> byIndexAndSign{}; // [2*index + (value > 0)], nondegenerate off-level points
    int degenerate = 0;                  // real off-level points with a singular Hessian
    bool certified = false;
    bool sameRegion(const Signature& o) const {
        return real == o.real && positive == o.positive && negative == o.negative &&
               byIndexAndSign == o.byIndexAndSign && degenerate == o.degenerate;
    }
};

struct Cell {
    double longitudeDegrees = 0.0, height = 0.0;
    Direction c{};
    Signature signature;
};

struct BoundaryPoint {
    Direction c{};
    double longitudeDegrees = 0.0, height = 0.0;
    bool singularFound = false;     // a real stationary point with H = 0 at the located direction
    int singularMultiplicity = 0;
    int singularHessianCorank = -1;
    double smallestRelativeValue = 0.0; // min |H| / term scale over the real stationary points
};

struct Sweep {
    double k = 0.0;
    int longitudes = 0, heights = 0, milnorNumber = -1;
    std::vector<Cell> cells;        // row-major: height index, then longitude index
    std::vector<BoundaryPoint> boundary;
};

inline Direction direction(double longitudeDegrees, double height) {
    const double phi = longitudeDegrees * 3.14159265358979323846 / 180.0;
    const double r = std::sqrt(std::max(0.0, 1.0 - height * height));
    return {r * std::cos(phi), r * std::sin(phi), height};
}

inline iso_poly::Cubic perturbed(const iso_poly::Cubic& f, double k, const Direction& c) {
    iso_poly::Cubic h = f;
    h.c[1] += k * c[0];
    h.c[2] += k * c[1];
    h.c[3] += k * c[2];
    return h;
}

inline iso_stationary::Options sweepSolverOptions(double regionRadius) {
    iso_stationary::Options o;
    o.regionRadius = regionRadius;
    o.startsPerScale = 12;
    o.scales = 30;
    // A narrow on-level band: bisection locates a discriminant direction only
    // up to the direction at which the value enters the band.
    o.onLevelTolerance = 1e-14;
    return o;
}

inline Signature signature(const iso_stationary::Result& r) {
    Signature s;
    s.certified = r.certified();
    for (const auto& p : r.points) {
        if (!p.real) continue;
        ++s.real;
        if (p.onZeroLevel) { ++s.zero; continue; }
        if (p.value > 0.0) ++s.positive;
        else ++s.negative;
        if (p.hessianCorank > 0) ++s.degenerate;
        else ++s.byIndexAndSign[static_cast<std::size_t>(2 * p.morseIndex + (p.value > 0.0 ? 1 : 0))];
    }
    return s;
}

namespace detail {

inline Direction slerp(const Direction& a, const Direction& b, double t) {
    const double cosine = std::clamp(a[0] * b[0] + a[1] * b[1] + a[2] * b[2], -1.0, 1.0);
    const double angle = std::acos(cosine);
    if (angle < 1e-15) return a;
    const double wa = std::sin((1.0 - t) * angle) / std::sin(angle), wb = std::sin(t * angle) / std::sin(angle);
    return {wa * a[0] + wb * b[0], wa * a[1] + wb * b[1], wa * a[2] + wb * b[2]};
}

} // namespace detail

// f: Taylor polynomial of F at P (P = origin of its coordinates), mu its Milnor number.
inline Sweep sweep(const iso_poly::Cubic& f, int mu, double k, int longitudes, int heights, double regionRadius) {
    Sweep out;
    out.k = k;
    out.longitudes = longitudes;
    out.heights = heights;
    out.milnorNumber = mu;
    const iso_stationary::Options options = sweepSolverOptions(regionRadius);
    auto solve = [&](const Direction& c) { return iso_stationary::findStationaryPoints(perturbed(f, k, c), mu, options); };

    out.cells.resize(static_cast<std::size_t>(longitudes * heights));
    iso_parallel::parallelFor(out.cells.size(), [&](std::size_t index) {
        Cell& cell = out.cells[index];
        const int i = static_cast<int>(index) % longitudes, j = static_cast<int>(index) / longitudes;
        cell.longitudeDegrees = -180.0 + 360.0 * (i + 0.5) / longitudes;
        cell.height = -1.0 + 2.0 * (j + 0.5) / heights;
        cell.c = direction(cell.longitudeDegrees, cell.height);
        cell.signature = signature(solve(cell.c));
    });

    // Neighbouring cells with the same number of real stationary points but a
    // different sign pattern enclose a discriminant direction.
    std::vector<std::pair<std::size_t, std::size_t>> pairs;
    for (int j = 0; j < heights; ++j)
        for (int i = 0; i < longitudes; ++i) {
            const std::size_t a = static_cast<std::size_t>(j * longitudes + i);
            const std::size_t east = static_cast<std::size_t>(j * longitudes + (i + 1) % longitudes);
            const Signature& sa = out.cells[a].signature;
            auto consider = [&](std::size_t b) {
                const Signature& sb = out.cells[b].signature;
                if (sa.real == sb.real && !sa.sameRegion(sb) && sa.zero == 0 && sb.zero == 0) pairs.emplace_back(a, b);
            };
            consider(east);
            if (j + 1 < heights) consider(static_cast<std::size_t>((j + 1) * longitudes + i));
        }
    out.boundary.resize(pairs.size());
    iso_parallel::parallelFor(pairs.size(), [&](std::size_t n) {
        Direction lo = out.cells[pairs[n].first].c, hi = out.cells[pairs[n].second].c;
        const Signature target = out.cells[pairs[n].first].signature;
        for (int it = 0; it < 48; ++it) {
            const Direction mid = detail::slerp(lo, hi, 0.5);
            if (signature(solve(mid)).sameRegion(target)) lo = mid;
            else hi = mid;
        }
        BoundaryPoint& b = out.boundary[n];
        b.c = detail::slerp(lo, hi, 0.5);
        b.height = b.c[2];
        b.longitudeDegrees = std::atan2(b.c[1], b.c[0]) * 180.0 / 3.14159265358979323846;
        const iso_stationary::Result r = solve(b.c);
        b.smallestRelativeValue = std::numeric_limits<double>::infinity();
        for (const auto& p : r.points) {
            if (!p.real) continue;
            const double rel = std::fabs(p.value) / std::max(p.valueScale, 1e-300);
            if (rel < b.smallestRelativeValue) {
                b.smallestRelativeValue = rel;
                b.singularMultiplicity = p.multiplicity;
                b.singularHessianCorank = p.hessianCorank;
            }
        }
        // After 48 halvings the direction is within ~1e-14 rad of the boundary;
        // the value at the stationary point is then zero up to rounding.
        b.singularFound = b.smallestRelativeValue <= 1e-12;
    });
    return out;
}

// Discriminant directions on the meridian through `longitudeDegrees`:
// c(theta) = (cos theta cos phi, cos theta sin phi, sin theta), theta in
// (-90, 90) degrees, scanned in `steps` steps and refined by bisection.
inline std::vector<BoundaryPoint> meridianBoundary(const iso_poly::Cubic& f, int mu, double k, double longitudeDegrees,
                                                   double regionRadius, int steps = 720) {
    const iso_stationary::Options options = sweepSolverOptions(regionRadius);
    auto solve = [&](const Direction& c) { return iso_stationary::findStationaryPoints(perturbed(f, k, c), mu, options); };
    auto at = [&](double t) { return direction(longitudeDegrees, std::sin(t)); };
    const double pi = 3.14159265358979323846;
    std::vector<double> theta(static_cast<std::size_t>(steps) + 1);
    std::vector<Signature> sig(theta.size());
    iso_parallel::parallelFor(theta.size(), [&](std::size_t i) {
        theta[i] = -0.5 * pi + pi * (static_cast<double>(i) + 0.5) / (steps + 1.0);
        sig[i] = signature(solve(at(theta[i])));
    });
    std::vector<BoundaryPoint> out;
    for (std::size_t i = 0; i + 1 < theta.size(); ++i) {
        if (sig[i].real != sig[i + 1].real || sig[i].sameRegion(sig[i + 1]) || sig[i].zero || sig[i + 1].zero) continue;
        double lo = theta[i], hi = theta[i + 1];
        for (int it = 0; it < 60; ++it) {
            const double mid = 0.5 * (lo + hi);
            if (signature(solve(at(mid))).sameRegion(sig[i])) lo = mid;
            else hi = mid;
        }
        BoundaryPoint b;
        b.c = at(0.5 * (lo + hi));
        b.height = b.c[2];
        b.longitudeDegrees = longitudeDegrees;
        b.smallestRelativeValue = std::numeric_limits<double>::infinity();
        for (const auto& p : solve(b.c).points) {
            if (!p.real) continue;
            const double rel = std::fabs(p.value) / std::max(p.valueScale, 1e-300);
            if (rel < b.smallestRelativeValue) {
                b.smallestRelativeValue = rel;
                b.singularMultiplicity = p.multiplicity;
                b.singularHessianCorank = p.hessianCorank;
            }
        }
        b.singularFound = b.smallestRelativeValue <= 1e-12;
        out.push_back(b);
    }
    return out;
}

inline bool writeCsv(const std::string& mapPath, const std::string& boundaryPath, const Sweep& s) {
    using iso_csv::number;
    std::ofstream m(mapPath);
    m << "longitude_deg,height,cx,cy,cz,real_stationary,positive,negative,zero,certified,k,milnor_number\n";
    for (const Cell& c : s.cells)
        m << number(c.longitudeDegrees) << ',' << number(c.height) << ',' << number(c.c[0]) << ',' << number(c.c[1])
          << ',' << number(c.c[2]) << ',' << c.signature.real << ',' << c.signature.positive << ','
          << c.signature.negative << ',' << c.signature.zero << ',' << (c.signature.certified ? 1 : 0) << ','
          << number(s.k) << ',' << s.milnorNumber << '\n';
    std::ofstream b(boundaryPath);
    b << "longitude_deg,height,cx,cy,cz,singular_found,singular_multiplicity,singular_hessian_corank,"
         "smallest_relative_value,k\n";
    for (const BoundaryPoint& p : s.boundary)
        b << number(p.longitudeDegrees) << ',' << number(p.height) << ',' << number(p.c[0]) << ',' << number(p.c[1])
          << ',' << number(p.c[2]) << ',' << (p.singularFound ? 1 : 0) << ',' << p.singularMultiplicity << ','
          << p.singularHessianCorank << ',' << number(p.smallestRelativeValue) << ',' << number(s.k) << '\n';
    return static_cast<bool>(m) && static_cast<bool>(b);
}

} // namespace iso_sweep
