// Dual3.h
//
// A tiny forward-mode automatic-differentiation number.
//
// Dual3 carries a value together with its three partial derivatives
// (d/dx, d/dy, d/dz). If you seed x = Dual3(x0, 1,0,0), y = Dual3(y0, 0,1,0),
// z = Dual3(z0, 0,0,1) and then evaluate an arithmetic expression using the
// overloads below, the result's .dx/.dy/.dz fields are the EXACT partial
// derivatives of the expression at (x0,y0,z0) -- no finite differences,
// no step-size tuning, no truncation error. This is what lets the app show
// the true analytic gradient of F at every point on the isosurface.
//
// k (the unfolding parameter) is passed in as a Dual3 with zero derivative
// components, since we only differentiate with respect to space here.
#pragma once
#include <cmath>
#include <algorithm>

struct Dual3 {
    double v  = 0.0;   // value
    double dx = 0.0;
    double dy = 0.0;
    double dz = 0.0;

    Dual3() = default;
    explicit Dual3(double val) : v(val), dx(0), dy(0), dz(0) {}
    Dual3(double val, double ddx, double ddy, double ddz)
        : v(val), dx(ddx), dy(ddy), dz(ddz) {}
};

// ---- arithmetic -----------------------------------------------------

inline Dual3 operator+(const Dual3& a, const Dual3& b) {
    return Dual3(a.v + b.v, a.dx + b.dx, a.dy + b.dy, a.dz + b.dz);
}
inline Dual3 operator-(const Dual3& a, const Dual3& b) {
    return Dual3(a.v - b.v, a.dx - b.dx, a.dy - b.dy, a.dz - b.dz);
}
inline Dual3 operator-(const Dual3& a) {
    return Dual3(-a.v, -a.dx, -a.dy, -a.dz);
}
inline Dual3 operator*(const Dual3& a, const Dual3& b) {
    // product rule: d(ab) = a'b + ab'
    return Dual3(a.v * b.v,
                 a.dx * b.v + a.v * b.dx,
                 a.dy * b.v + a.v * b.dy,
                 a.dz * b.v + a.v * b.dz);
}
inline Dual3 operator/(const Dual3& a, const Dual3& b) {
    // quotient rule: d(a/b) = (a'b - ab') / b^2
    double inv = 1.0 / b.v;
    double inv2 = inv * inv;
    return Dual3(a.v * inv,
                 (a.dx * b.v - a.v * b.dx) * inv2,
                 (a.dy * b.v - a.v * b.dy) * inv2,
                 (a.dz * b.v - a.v * b.dz) * inv2);
}

// ---- elementary functions, overloaded for plain double and for Dual3.
// The expression evaluator (Expr.h) is a template that calls these by
// unqualified name, so the same AST-walking code works whether it is
// producing a plain value (T = double) or a value-with-gradient (T = Dual3).

inline double mc_sin(double a)  { return std::sin(a); }
inline double mc_cos(double a)  { return std::cos(a); }
inline double mc_tan(double a)  { return std::tan(a); }
inline double mc_exp(double a)  { return std::exp(a); }
inline double mc_log(double a)  { return std::log(a); }
inline double mc_sqrt(double a) { return std::sqrt(a); }
inline double mc_abs(double a)  { return std::fabs(a); }
inline double mc_pow(double a, double b) { return std::pow(a, b); }
inline double mc_min(double a, double b) { return std::min(a, b); }
inline double mc_max(double a, double b) { return std::max(a, b); }

inline Dual3 mc_sin(const Dual3& a) {
    double c = std::cos(a.v);
    return Dual3(std::sin(a.v), c * a.dx, c * a.dy, c * a.dz);
}
inline Dual3 mc_cos(const Dual3& a) {
    double s = -std::sin(a.v);
    return Dual3(std::cos(a.v), s * a.dx, s * a.dy, s * a.dz);
}
inline Dual3 mc_tan(const Dual3& a) {
    double c = std::cos(a.v);
    double sec2 = 1.0 / (c * c);
    return Dual3(std::tan(a.v), sec2 * a.dx, sec2 * a.dy, sec2 * a.dz);
}
inline Dual3 mc_exp(const Dual3& a) {
    double e = std::exp(a.v);
    return Dual3(e, e * a.dx, e * a.dy, e * a.dz);
}
inline Dual3 mc_log(const Dual3& a) {
    double inv = 1.0 / a.v;
    return Dual3(std::log(a.v), inv * a.dx, inv * a.dy, inv * a.dz);
}
inline Dual3 mc_sqrt(const Dual3& a) {
    double r = std::sqrt(a.v);
    double coeff = (r > 1e-12) ? 0.5 / r : 0.0; // guard the singular case
    return Dual3(r, coeff * a.dx, coeff * a.dy, coeff * a.dz);
}
inline Dual3 mc_abs(const Dual3& a) {
    double s = (a.v > 0) ? 1.0 : (a.v < 0 ? -1.0 : 0.0);
    return Dual3(std::fabs(a.v), s * a.dx, s * a.dy, s * a.dz);
}
inline Dual3 mc_pow(const Dual3& a, const Dual3& b) {
    // General d/dt[a^b] = a^b * ( b' ln(a) + b a'/a ).
    // Most expressions use a constant exponent (b.dx==b.dy==b.dz==0), in
    // which case this correctly reduces to the ordinary power rule
    // b * a^(b-1) * a' for any sign of a.
    double val = std::pow(a.v, b.v);
    double coeffA = b.v * std::pow(a.v, b.v - 1.0);
    if (b.dx == 0.0 && b.dy == 0.0 && b.dz == 0.0) {
        return Dual3(val, coeffA * a.dx, coeffA * a.dy, coeffA * a.dz);
    }
    double lnA = std::log(a.v); // only reached with a variable exponent; needs a.v>0
    double coeffB = val * lnA;
    return Dual3(val,
                 coeffA * a.dx + coeffB * b.dx,
                 coeffA * a.dy + coeffB * b.dy,
                 coeffA * a.dz + coeffB * b.dz);
}
inline Dual3 mc_min(const Dual3& a, const Dual3& b) { return (a.v <= b.v) ? a : b; }
inline Dual3 mc_max(const Dual3& a, const Dual3& b) { return (a.v >= b.v) ? a : b; }
