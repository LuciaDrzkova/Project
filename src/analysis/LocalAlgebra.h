#pragma once

// Local algebra of an isolated hypersurface singularity at the origin.
//
// For a polynomial F with F(0) = 0 the Milnor algebra and the Tjurina algebra
// at the origin are the finite-dimensional vector spaces
//     M(F) = O / (dF/dx, dF/dy, dF/dz),        mu  = dim M(F)   (Milnor number)
//     T(F) = O / (F, dF/dx, dF/dy, dF/dz),     tau = dim T(F)   (Tjurina number)
// where O is the ring of power series (germs) at the origin. A monomial basis
// of T(F) is the list of "essential directions" of a perturbation: the class
// of G in T(F) is the first-order part of G that cannot be removed by a change
// of coordinates and a multiplication of F by a unit (Greuel, Lossen, Shustin,
// "Introduction to Singularities and Deformations", Sec. II.1).
//
// Computation. Working modulo m^N (m = the maximal ideal), the quotient
// O/(I + m^N) is a finite linear-algebra problem on the monomials of degree
// < N. Its dimension d_N is nondecreasing in N; if d_N = d_{N+1}, then
// I + m^N = I + m^{N+1}, and Nakayama's lemma gives m^N in I, hence
// dim O/I = d_N exactly. Gaussian elimination with the columns ordered by
// increasing degree takes the lowest-degree term of every element as its
// leading term (a local ordering), so the non-pivot monomials form a basis of
// O/I made of monomials of the smallest possible degrees.

#include "CubicPolynomial.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <string>
#include <vector>

namespace iso_algebra {

struct Monomial {
    int x = 0, y = 0, z = 0;
    int degree() const { return x + y + z; }
    std::string text() const {
        std::string s;
        auto add = [&](const char* v, int p) {
            if (!p) return;
            if (!s.empty()) s += '*';
            s += v;
            if (p > 1) s += '^' + std::to_string(p);
        };
        add("x", x); add("y", y); add("z", z);
        return s.empty() ? "1" : s;
    }
};

using Polynomial = std::vector<std::pair<Monomial, double>>;

enum class IdealKind { Milnor, Tjurina };

namespace detail {

// All monomials of degree < n, ordered by degree, then x^a y^b z^c with a, then b, descending.
inline std::vector<Monomial> monomialsBelow(int n) {
    std::vector<Monomial> out;
    for (int d = 0; d < n; ++d)
        for (int a = d; a >= 0; --a)
            for (int b = d - a; b >= 0; --b) out.push_back({a, b, d - a - b});
    return out;
}

inline Polynomial fromCubic(const iso_poly::Cubic& p) {
    Polynomial out;
    const auto& e = iso_research::exponents();
    for (std::size_t i = 0; i < 20; ++i)
        if (p.c[i] != 0.0) out.push_back({{e[i].x, e[i].y, e[i].z}, p.c[i]});
    return out;
}

inline Polynomial partial(const Polynomial& p, int axis) {
    Polynomial out;
    for (const auto& [m, c] : p) {
        const int power = axis == 0 ? m.x : (axis == 1 ? m.y : m.z);
        if (power == 0) continue;
        Monomial d = m;
        (axis == 0 ? d.x : (axis == 1 ? d.y : d.z)) -= 1;
        out.push_back({d, c * power});
    }
    return out;
}

// Row-reduced echelon form over the monomials of degree < n, built row by row.
class Echelon {
public:
    explicit Echelon(int n) : n_(n), monomials_(monomialsBelow(n)) {
        index_.assign(static_cast<std::size_t>(n * n * n), -1);
        for (std::size_t i = 0; i < monomials_.size(); ++i) {
            const Monomial& m = monomials_[i];
            index_[static_cast<std::size_t>((m.x * n + m.y) * n + m.z)] = static_cast<int>(i);
        }
        pivotRow_.assign(monomials_.size(), -1);
    }

    std::size_t columns() const { return monomials_.size(); }
    std::size_t rank() const { return rows_.size(); }
    const std::vector<Monomial>& monomials() const { return monomials_; }
    bool isPivot(std::size_t column) const { return pivotRow_[column] >= 0; }

    // Dense vector of p truncated to degree < n.
    std::vector<double> vector(const Polynomial& p, const Monomial& shift = {}) const {
        std::vector<double> v(monomials_.size(), 0.0);
        for (const auto& [m, c] : p) {
            const Monomial s{m.x + shift.x, m.y + shift.y, m.z + shift.z};
            if (s.degree() >= n_) continue;
            v[static_cast<std::size_t>(index_[static_cast<std::size_t>((s.x * n_ + s.y) * n_ + s.z)])] += c;
        }
        return v;
    }

    // Subtracts the stored rows; what remains lives on the non-pivot columns.
    void reduce(std::vector<double>& v) const {
        for (std::size_t col = 0; col < v.size(); ++col) {
            if (v[col] == 0.0 || pivotRow_[col] < 0) continue;
            const std::vector<double>& row = rows_[static_cast<std::size_t>(pivotRow_[col])];
            const double f = v[col];
            for (std::size_t j = col; j < v.size(); ++j) v[j] -= f * row[j];
            v[col] = 0.0;
        }
    }

    void add(std::vector<double> v) {
        double scale = 0.0;
        for (double a : v) scale = std::max(scale, std::fabs(a));
        if (scale == 0.0) return;
        reduce(v);
        const double tolerance = 1e-10 * scale;
        std::size_t pivot = v.size();
        for (std::size_t j = 0; j < v.size(); ++j) {
            if (std::fabs(v[j]) > tolerance) { pivot = j; break; }
            v[j] = 0.0;
        }
        if (pivot == v.size()) return;
        const double inv = 1.0 / v[pivot];
        for (std::size_t j = pivot; j < v.size(); ++j) v[j] *= inv;
        v[pivot] = 1.0;
        for (std::vector<double>& row : rows_) { // keep the form fully reduced
            const double f = row[pivot];
            if (f == 0.0) continue;
            for (std::size_t j = pivot; j < row.size(); ++j) row[j] -= f * v[j];
            row[pivot] = 0.0;
        }
        pivotRow_[pivot] = static_cast<int>(rows_.size());
        rows_.push_back(std::move(v));
    }

private:
    int n_;
    std::vector<Monomial> monomials_;
    std::vector<int> index_;
    std::vector<int> pivotRow_;
    std::vector<std::vector<double>> rows_;
};

inline Echelon idealModulo(const std::vector<Polynomial>& generators, int n) {
    Echelon e(n);
    const std::vector<Monomial> shifts = monomialsBelow(n);
    for (const Polynomial& g : generators) {
        int low = n;
        for (const auto& [m, c] : g) low = std::min(low, m.degree());
        for (const Monomial& s : shifts)
            if (s.degree() + low < n) e.add(e.vector(g, s));
    }
    return e;
}

} // namespace detail

struct LocalAlgebra {
    IdealKind kind = IdealKind::Tjurina;
    int dimension = -1;              // mu or tau; -1: not finite within the degree limit
    int truncationDegree = 0;        // N with m^N contained in the ideal
    std::vector<Monomial> basis;     // monomial basis of O/I, lowest degrees first
    bool valid() const { return dimension >= 0; }

    std::string basisText() const {
        std::string s;
        for (const Monomial& m : basis) s += (s.empty() ? "" : " ") + m.text();
        return s;
    }

    // Coordinates of the class of g in O/I with respect to `basis`.
    std::vector<double> coordinates(const iso_poly::Cubic& g) const {
        if (!valid()) return {};
        std::vector<double> v = echelon_.vector(detail::fromCubic(g));
        echelon_.reduce(v);
        std::vector<double> out;
        for (std::size_t col = 0; col < v.size(); ++col)
            if (!echelon_.isPivot(col)) out.push_back(v[col]);
        return out;
    }

    detail::Echelon echelon_{1};
};

// Milnor (gradient ideal) or Tjurina (F and its gradient) algebra of F at the
// origin. F(0) must be 0 for the Tjurina algebra to be the local one.
inline LocalAlgebra localAlgebra(const iso_poly::Cubic& f, IdealKind kind, int maxDegree = 12) {
    const Polynomial p = detail::fromCubic(f);
    std::vector<Polynomial> generators = {detail::partial(p, 0), detail::partial(p, 1), detail::partial(p, 2)};
    if (kind == IdealKind::Tjurina) generators.push_back(p);
    generators.erase(std::remove_if(generators.begin(), generators.end(),
                                    [](const Polynomial& g) { return g.empty(); }),
                     generators.end());
    LocalAlgebra out;
    out.kind = kind;
    std::size_t previous = static_cast<std::size_t>(-1);
    for (int n = 1; n <= maxDegree; ++n) {
        detail::Echelon e = detail::idealModulo(generators, n);
        const std::size_t dim = e.columns() - e.rank();
        if (dim == previous) {
            out.dimension = static_cast<int>(dim);
            out.truncationDegree = n;
            for (std::size_t col = 0; col < e.columns(); ++col)
                if (!e.isPivot(col)) out.basis.push_back(e.monomials()[col]);
            out.echelon_ = std::move(e);
            return out;
        }
        previous = dim;
    }
    return out; // not stabilised: the singularity is not isolated (or mu is very large)
}

} // namespace iso_algebra
