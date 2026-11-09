// Local algebra, certified stationary points, exact surface distances and the
// direction sweep, checked against results known in closed form.
#include "analysis/ConvergenceClassifier.h"
#include "analysis/DirectionSweep.h"

#include <cmath>
#include <iostream>
#include <string>

namespace {
int failures = 0;

void check(bool condition, const std::string& message)
{
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

bool near(double a, double b, double tolerance) { return std::isfinite(a) && std::fabs(a - b) <= tolerance; }

iso_poly::Cubic cubic(const std::string& text, double k = 0.0)
{
    const NodePtr ast = parseExpression(text);
    return *iso_poly::Cubic::fromExpression(ast.get(), k);
}

int milnor(const std::string& f) { return iso_algebra::localAlgebra(cubic(f), iso_algebra::IdealKind::Milnor).dimension; }
int tjurina(const std::string& f) { return iso_algebra::localAlgebra(cubic(f), iso_algebra::IdealKind::Tjurina).dimension; }

void testLocalAlgebra()
{
    // Milnor numbers of the ten reference surfaces: A1 = 1, A2 = 2, D4 = 4, E6~ = 8.
    check(milnor("x^2+y^2-z^2+x*y*z") == 1, "A1: mu = 1");
    check(milnor("x^2+y^2+z^3") == 2 && milnor("x^2-y^2+z^3") == 2, "A2: mu = 2");
    check(milnor("x^3+y^3+z^2+x^2") == 2, "x^3+y^3+z^2+x^2 is A2 (y^3 is the cubic direction)");
    check(milnor("x^2+y^3+z^3") == 4 && milnor("x^2+y^3-3*y*z^2") == 4, "D4: mu = 4");
    check(milnor("x^3+y^3+z^3") == 8 && milnor("x^3+y^3-z^3") == 8, "Fermat cones: mu = 8");
    check(milnor("x^3+y^3+z^3+x*y*z") == 8 && milnor("x^3+2*y^3+3*z^3+2*x*y*z") == 8, "other cubic cones: mu = 8");
    // Quasi-homogeneous singularities have tau = mu (Saito).
    check(tjurina("x^3+y^3+z^3") == 8 && tjurina("x^2+y^3+z^3") == 4, "tau = mu for quasi-homogeneous F");
    check(milnor("x^2+y^2") == -1, "non-isolated singularity (a line) is detected");

    // Essential coordinates for the Fermat cone: basis 1,x,y,z,xy,xz,yz,xyz.
    const auto t = iso_algebra::localAlgebra(cubic("x^3+y^3+z^3"), iso_algebra::IdealKind::Tjurina);
    check(t.basisText() == "1 x y z x*y x*z y*z x*y*z", "Fermat cone Tjurina basis, got " + t.basisText());
    auto coords = [&](const std::string& g) { return t.coordinates(cubic(g)); };
    const auto cx2 = coords("x^2");
    double n2 = 0.0;
    for (double v : cx2) n2 += v * v;
    check(n2 < 1e-20, "x^2 is inessential for the Fermat cone (it lies in the gradient ideal)");
    const auto cxy = coords("2*x*y+3");
    check(near(cxy[0], 3.0, 1e-12) && near(cxy[4], 2.0, 1e-12), "2xy+3 has coordinates 3 (on 1) and 2 (on xy)");
    const auto cF = coords("x^3+y^3+z^3");
    double nF = 0.0;
    for (double v : cF) nF += v * v;
    check(nF < 1e-20, "G = F is inessential (F lies in its Tjurina ideal)");
}

void testCertifiedStationaryPoints()
{
    // F_+ + k(x+y), k < 0: two cusps (A2, multiplicity 2 each) on the zero
    // level at (+-s, -+s, 0), s = sqrt(-k/3); all 8 = mu points accounted for.
    const double k = -0.01, s = std::sqrt(-k / 3.0);
    const auto r = iso_stationary::findStationaryPoints(cubic("x^3+y^3+z^3+k*(x+y)", k), 8);
    check(r.certified(), "x+y: certified complete, status " + r.status);
    check(r.singularCount() == 2 && r.singularMultiplicity() == 4, "x+y: two cusps of multiplicity 2");
    for (const auto& p : r.points)
        if (p.onZeroLevel)
            check(near(std::fabs(p.point.x), s, 1e-9) && near(p.point.x + p.point.y, 0.0, 1e-9) &&
                      near(p.point.z, 0.0, 1e-9) && p.hessianCorank == 1,
                  "x+y: cusp position and Hessian corank 1");
    // k > 0: the stationary points are complex; nothing real.
    const auto rp = iso_stationary::findStationaryPoints(cubic("x^3+y^3+z^3+k*(x+y)", 0.01), 8);
    check(rp.certified() && rp.realCount() == 0, "x+y, k>0: certified, no real stationary point");
    // G = x^2: two D4 points (multiplicity 4 each), only the origin on the zero level.
    const auto rq = iso_stationary::findStationaryPoints(cubic("x^3+y^3+z^3+k*x^2", -0.01), 8);
    check(rq.certified() && rq.points.size() == 2 && rq.singularCount() == 1, "x^2: two D4 points, origin singular");
    // G = 1: the origin stays the only stationary point (multiplicity 8) and leaves the surface.
    const auto r1 = iso_stationary::findStationaryPoints(cubic("x^3+y^3+z^3+k", 0.01), 8);
    check(r1.certified() && r1.points.size() == 1 && r1.points[0].multiplicity == 8 && r1.singularCount() == 0,
          "G = 1: one stationary point of multiplicity 8, off the zero level");
    // A1 + k: one nondegenerate stationary point.
    const auto ra = iso_stationary::findStationaryPoints(cubic("x^2+y^2-z^2+x*y*z+k", 0.01), 1);
    check(ra.certified() && ra.points.size() == 1 && ra.points[0].multiplicity == 1, "A1 + k: one simple point");
}

void testExactDistances()
{
    // dist(0, {x^3+y^3+z^3 = -k}) = |k|^(1/3) exactly.
    const std::string family = "x^3+y^3+z^3+k";
    const NodePtr ast = parseExpression(family);
    iso_convergence::Options o;
    o.resolution = 48;
    o.halfExtent = 0.4;
    o.localHalfSize = 0.18;
    o.kValues = {-0.001, 0.001, 0.008};
    o.familyExpression = family;
    const auto a = iso_research::analyzeFamily(ast.get(), o, 2);
    for (const auto& m : a.exact) {
        check(near(m.probeDistance, std::cbrt(std::fabs(m.k)), 1e-10),
              "exact dist(0, S_k) = |k|^(1/3), got " + std::to_string(m.probeDistance));
        check(m.comparable && m.meshHausdorff > 0.0 && m.meshHausdorff < 0.5 * 0.8 / 48,
              "extraction error below half a cell");
    }
    // S_k = |k|^(1/3) S_{sign k} and S_0 is a cone, so d_Q scales exactly like |k|^(1/3)
    // (while S_k stays inside the box), and S_{-k} = -S_k gives the same value for -k.
    check(near(a.exact[0].symmetric, a.exact[1].symmetric, 1e-12), "G = 1: d_Q(-k) = d_Q(k)");
    check(near(a.exact[2].symmetric / a.exact[1].symmetric, 2.0, 1e-6), "G = 1: d_Q(8k) = 2 d_Q(k)");
    // G = F: identical surfaces, distance 0.
    const NodePtr same = parseExpression("(x^3+y^3+z^3)+k*(x^3+y^3+z^3)");
    o.kValues = {0.01};
    const auto b = iso_research::analyzeFamily(same.get(), o, 2);
    check(b.exact[0].comparable && b.exact[0].symmetric < 1e-9, "G = F: exact local Hausdorff distance 0");
}

void testClassification()
{
    const std::string family = "(x^3+y^3+z^3)+k*(x+y)";
    const NodePtr ast = parseExpression(family);
    iso_convergence::Options o;
    o.resolution = 48;
    o.halfExtent = 0.4;
    o.localHalfSize = 0.18;
    o.kValues = iso_classify::makeSymmetricLogKValues(0.01, 6, 2.0);
    o.familyExpression = family;
    const auto c = iso_classify::classify(ast.get(), iso_research::analyzeFamily(ast.get(), o, 2));
    check(c.milnorNumber == 8 && c.tjurinaNumber == 8, "classification reports mu = tau = 8");
    check(c.essentialCoordinates == "x + y", "essential part of x+y is x+y, got " + c.essentialCoordinates);
    check(c.negative.exactDistances && near(c.negative.hausdorff.exponent, 0.5, 0.05),
          "x+y, k<0: exact rate 1/2, got " + std::to_string(c.negative.hausdorff.exponent));
    check(c.negative.hausdorff.ciLow <= c.negative.hausdorff.exponent &&
              c.negative.hausdorff.exponent <= c.negative.hausdorff.ciHigh,
          "confidence interval contains the estimate");
    check(c.negative.certifiedK == 6 && c.negative.maxSingularPoints == 2 && c.negative.movingSingularMultiplicity == 2,
          "x+y, k<0: certified, two moving cusps of Milnor number 2");
    check(near(c.negative.singularDistance.exponent, 0.5, 1e-6), "cusps approach P like |k|^(1/2)");
    check(c.positive.singular == iso_classify::SingularBehavior::Regularized, "x+y, k>0: no singular points");
}

void testSweep()
{
    // F_+ at k < 0: the discriminant on the meridian through (1,1,0)
    // contains the balanced direction (1, 1, 2^(2/3)) / |.|.
    const auto boundary = iso_sweep::meridianBoundary(cubic("x^3+y^3+z^3"), 8, -0.01, 45.0, 0.31, 360);
    bool found = false;
    for (const auto& b : boundary)
        if (b.singularFound && b.c[2] > 0.0)
            found = found || near(b.c[2] / b.c[0], std::cbrt(4.0), 1e-9);
    check(found, "meridian sweep finds x + y + 2^(2/3) z");
}

} // namespace

int main()
{
    testLocalAlgebra();
    testCertifiedStationaryPoints();
    testExactDistances();
    testClassification();
    testSweep();
    if (failures) {
        std::cerr << failures << " check(s) failed\n";
        return 1;
    }
    std::cout << "exact analysis tests passed\n";
    return 0;
}
