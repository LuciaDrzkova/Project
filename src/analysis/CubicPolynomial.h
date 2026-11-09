#pragma once

// Exact representation of a polynomial of total degree <= 3 in x, y, z, with
// fast evaluation of value, gradient and Hessian in real or complex
// arithmetic. Every surface studied in the thesis is of this kind; the exact
// form is what makes the certified stationary-point count
// (StationaryPoints.h), the local algebra (LocalAlgebra.h) and the exact
// surface distances (ExactDistance.h) possible.
//
// Coefficients use the monomial order of iso_research::exponents():
//   1, x, y, z, x^2, xy, xz, y^2, yz, z^2, x^3, x^2y, x^2z, xy^2, xyz, xz^2, y^3, y^2z, yz^2, z^3.

#include "SingularityResearch.h"

#include <array>
#include <cmath>
#include <complex>
#include <optional>
#include <string>

namespace iso_poly {

using Point3 = iso_geometry::Point3;

template <class T>
struct Derivatives {
    T value{};
    std::array<T, 3> gradient{};
    std::array<std::array<T, 3>, 3> hessian{};
};

struct Cubic {
    std::array<double, 20> c{};

    // Exact coefficients of `ast` at parameter k, when the expression is a
    // polynomial of total degree <= 3 in x, y, z. The Taylor jet at the origin
    // of such a polynomial is the polynomial itself.
    static std::optional<Cubic> fromExpression(const Node* ast, double k = 0.0) {
        if (!ast) return std::nullopt;
        const iso_convergence::DegreeInfo degree = iso_convergence::spatialPolynomialDegree(ast);
        if (!degree.polynomial || degree.degree > 3) return std::nullopt;
        const iso_research::DifferentialData d = iso_research::evaluateJet(ast, {0.0, 0.0, 0.0}, k);
        if (!d.smooth) return std::nullopt;
        for (double v : d.taylor)
            if (!std::isfinite(v)) return std::nullopt;
        Cubic p;
        p.c = d.taylor;
        return p;
    }

    Cubic operator+(const Cubic& o) const { Cubic r; for (int i = 0; i < 20; ++i) r.c[i] = c[i] + o.c[i]; return r; }
    Cubic operator*(double s) const { Cubic r; for (int i = 0; i < 20; ++i) r.c[i] = c[i] * s; return r; }

    double maxAbsCoefficient() const {
        double m = 0.0;
        for (double v : c) m = std::max(m, std::fabs(v));
        return m;
    }

    // Value, gradient and Hessian at q (T = double or std::complex<double>).
    template <class T>
    Derivatives<T> evaluate(const std::array<T, 3>& q) const {
        const T x = q[0], y = q[1], z = q[2];
        Derivatives<T> d;
        // value
        d.value = T(c[0]) + c[1] * x + c[2] * y + c[3] * z
                + c[4] * x * x + c[5] * x * y + c[6] * x * z + c[7] * y * y + c[8] * y * z + c[9] * z * z
                + c[10] * x * x * x + c[11] * x * x * y + c[12] * x * x * z + c[13] * x * y * y
                + c[14] * x * y * z + c[15] * x * z * z + c[16] * y * y * y + c[17] * y * y * z
                + c[18] * y * z * z + c[19] * z * z * z;
        // gradient
        d.gradient[0] = T(c[1]) + 2.0 * c[4] * x + c[5] * y + c[6] * z
                      + 3.0 * c[10] * x * x + 2.0 * c[11] * x * y + 2.0 * c[12] * x * z + c[13] * y * y
                      + c[14] * y * z + c[15] * z * z;
        d.gradient[1] = T(c[2]) + c[5] * x + 2.0 * c[7] * y + c[8] * z
                      + c[11] * x * x + 2.0 * c[13] * x * y + c[14] * x * z + 3.0 * c[16] * y * y
                      + 2.0 * c[17] * y * z + c[18] * z * z;
        d.gradient[2] = T(c[3]) + c[6] * x + c[8] * y + 2.0 * c[9] * z
                      + c[12] * x * x + c[14] * x * y + 2.0 * c[15] * x * z + c[17] * y * y
                      + 2.0 * c[18] * y * z + 3.0 * c[19] * z * z;
        // Hessian
        const T hxx = 2.0 * c[4] + 6.0 * c[10] * x + 2.0 * c[11] * y + 2.0 * c[12] * z;
        const T hxy = T(c[5]) + 2.0 * c[11] * x + 2.0 * c[13] * y + c[14] * z;
        const T hxz = T(c[6]) + 2.0 * c[12] * x + c[14] * y + 2.0 * c[15] * z;
        const T hyy = 2.0 * c[7] + 2.0 * c[13] * x + 6.0 * c[16] * y + 2.0 * c[17] * z;
        const T hyz = T(c[8]) + c[14] * x + 2.0 * c[17] * y + 2.0 * c[18] * z;
        const T hzz = 2.0 * c[9] + 2.0 * c[15] * x + 2.0 * c[18] * y + 6.0 * c[19] * z;
        d.hessian = {{{hxx, hxy, hxz}, {hxy, hyy, hyz}, {hxz, hyz, hzz}}};
        return d;
    }

    Derivatives<double> evaluate(const Point3& p) const { return evaluate<double>({p.x, p.y, p.z}); }

    // Sum of the absolute values of the terms at p: the size of the numbers
    // that cancel when the value is (close to) zero.
    double termScale(const Point3& p) const {
        const auto& e = iso_research::exponents();
        double s = 0.0;
        for (std::size_t i = 0; i < 20; ++i)
            if (c[i] != 0.0)
                s += std::fabs(c[i]) * std::pow(std::fabs(p.x), e[i].x) * std::pow(std::fabs(p.y), e[i].y) *
                     std::pow(std::fabs(p.z), e[i].z);
        return s;
    }

    std::string text() const {
        const auto& e = iso_research::exponents();
        std::string out;
        for (std::size_t i = 0; i < 20; ++i) {
            if (c[i] == 0.0) continue;
            const double a = std::fabs(c[i]);
            out += out.empty() ? (c[i] < 0 ? "-" : "") : (c[i] < 0 ? " - " : " + ");
            std::string mono;
            auto add = [&](const char* v, int p) {
                if (!p) return;
                if (!mono.empty()) mono += '*';
                mono += v;
                if (p > 1) mono += '^' + std::to_string(p);
            };
            add("x", e[i].x); add("y", e[i].y); add("z", e[i].z);
            if (mono.empty()) out += iso_csv::shortNumber(a);
            else out += (a == 1.0 ? "" : iso_csv::shortNumber(a) + "*") + mono;
        }
        return out.empty() ? "0" : out;
    }
};

} // namespace iso_poly
