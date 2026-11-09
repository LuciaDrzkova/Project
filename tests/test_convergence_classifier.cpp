// Convergence classifier: predicted exponents for known singularity types and
// end-to-end classification of families whose behaviour is known analytically.
#include "analysis/ConvergenceClassifier.h"

#include <cmath>
#include <iostream>
#include <string>

namespace {
int failures = 0;

void check(bool condition, const std::string& message)
{
    if (!condition)
    {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

bool near(double a, double b, double tolerance) { return std::isfinite(a) && std::fabs(a - b) <= tolerance; }

// Prediction only: Taylor jets of F and G at the origin, no meshes.
iso_classify::Prediction predictFor(const std::string& f, const std::string& g)
{
    const NodePtr fAst = parseExpression(f);
    const NodePtr gAst = parseExpression(g);
    iso_research::ResearchAnalysis r;
    r.affineStatus = "STRUCTURALLY_AFFINE_IN_K";
    r.fTaylor.data = iso_research::evaluateJet(fAst.get(), {0, 0, 0}, 0.0);
    r.fTaylor.firstNonzeroOrder = iso_research::firstNonzeroJetOrder(r.fTaylor.data.taylor);
    r.gTaylor.data = iso_research::evaluateJet(gAst.get(), {0, 0, 0}, 0.0);
    r.gTaylor.firstNonzeroOrder = iso_research::firstNonzeroJetOrder(r.gTaylor.data.taylor);
    return iso_classify::detail::predict(r);
}

// Full pipeline: analyzeFamily on log-spaced k, then classify.
iso_classify::Classification classifyFamily(const std::string& f, const std::string& g, double kMax)
{
    const std::string family = "(" + f + ")+k*(" + g + ")";
    const NodePtr ast = parseExpression(family);
    iso_convergence::Options o;
    o.probe = {0, 0, 0};
    o.resolution = 64;
    o.halfExtent = 0.5;
    o.localHalfSize = 0.45 * o.halfExtent;
    o.kValues = iso_classify::makeSymmetricLogKValues(kMax, 9, 2.0);
    o.familyExpression = family;
    return iso_classify::classify(ast.get(), iso_research::analyzeFamily(ast.get(), o, 2));
}

void testPredictions()
{
    using iso_classify::Prediction;
    // E6-tilde (cone over a smooth plane cubic): homogeneous weights 1/3, alpha = 1/(3-m).
    Prediction p = predictFor("x^3+y^3+z^3", "1");
    check(p.principalNondegenerate && near(p.weights[0], 1.0 / 3.0, 1e-9), "Fermat cubic has weights 1/3");
    check(near(p.exponent, 1.0 / 3.0, 1e-9), "cubic + k: alpha = 1/3");
    check(near(predictFor("x^3+y^3+z^3", "x").exponent, 0.5, 1e-9), "cubic + kx: alpha = 1/2");
    check(near(predictFor("x^3+y^3+z^3", "2*x+4*y+6*z").exponent, 0.5, 1e-9), "cubic + k(linear): alpha = 1/2");
    check(near(predictFor("x^3+y^3+z^3", "x^2").exponent, 1.0, 1e-9), "cubic + kx^2: alpha = 1");
    check(predictFor("x^3+y^3+z^3", "x*y*z").regime == "SAME_WEIGHT", "cubic + k*cubic is a same-weight deformation");
    check(predictFor("x^3+y^3+z^3+x*y*z", "1").principalNondegenerate, "Hesse cubic x^3+y^3+z^3+xyz is nondegenerate");

    // A1 (Morse cone): weights 1/2.
    check(near(predictFor("x^2+y^2-z^2+x*y*z", "1").exponent, 0.5, 1e-9), "A1 + k: alpha = 1/2");
    check(near(predictFor("x^2+y^2-z^2+x*y*z", "x").exponent, 1.0, 1e-9), "A1 + kx: alpha = 1");

    // A2 (x^2+y^2+z^3): weights (1/2, 1/2, 1/3).
    p = predictFor("x^2+y^2+z^3", "x");
    check(near(p.weights[2], 1.0 / 3.0, 1e-9) && near(p.weights[0], 0.5, 1e-9), "A2 weights (1/2,1/2,1/3)");
    check(near(p.exponent, 2.0 / 3.0, 1e-9), "A2 + kx: alpha = 2/3 (lambda=|k|^2, z-scale lambda^(1/3))");
    check(near(predictFor("x^2+y^2+z^3", "z").exponent, 0.5, 1e-9), "A2 + kz: alpha = 1/2");
    p = predictFor("x^2+y^2+z^3", "1");
    check(near(p.exponent, 1.0 / 3.0, 1e-9) && near(p.exponentMax, 0.5, 1e-9), "A2 + k: alpha in [1/3, 1/2]");
    check(near(predictFor("x^3+y^3+z^3", "x").exponentMax, 0.5, 1e-9), "homogeneous F: a single exponent");

    // D4 (x^2 + y^3 + z^3 and its real form x^2 + y^3 - 3yz^2): weights (1/2, 1/3, 1/3).
    check(near(predictFor("x^2+y^3+z^3", "1").exponent, 1.0 / 3.0, 1e-9), "D4 + k: alpha = 1/3");
    check(near(predictFor("x^2+y^3-3*y*z^2", "y").exponent, 0.5, 1e-9), "D4- + ky: alpha = 1/2");

    // Degenerate principal part: x^2 + y^2 is singular along the whole z-axis.
    check(!std::isfinite(predictFor("x^2+y^2", "1").exponent), "non-isolated singularity gets no prediction");
}

void testRegularizingLinearPerturbation()
{
    // P stays on S_k, becomes regular with |grad F_k(P)| = |k|; d_H ~ |k|^(1/2).
    const auto c = classifyFamily("x^3+y^3+z^3", "x", 0.01);
    check(c.point == iso_classify::PointRelation::OnSurfaceRegular, "x: P stays on S_k as a regular point");
    check(near(c.gradGAtProbe, 1.0, 1e-12), "x: gradient shift coefficient |grad G(P)| = 1");
    for (const auto* s : {&c.negative, &c.positive})
    {
        check(s->surface == iso_classify::SurfaceClass::Sublinear, "x: Hoelder convergence on both sides: " + s->verdict);
        check(near(s->hausdorff.exponent, 0.5, 0.15), "x: measured alpha near 1/2: " + s->verdict);
        check(s->singular == iso_classify::SingularBehavior::Regularized, "x: singularity removed: " + s->verdict);
        check(s->hessianCorankAtP == 3, "x: Hessian of F_k at P stays zero (corank 3)");
    }
}

void testMovingSingularPoints()
{
    // k<0: two singular points (±sqrt(-k/3), ∓sqrt(-k/3), 0) on S_k move into P.
    const auto c = classifyFamily("x^3+y^3+z^3", "x+y", 0.01);
    check(c.negative.singular == iso_classify::SingularBehavior::Moving, "x+y, k<0: moving singular points: " + c.negative.verdict);
    check(c.negative.maxSingularPoints == 2, "x+y, k<0: exactly two singular points");
    check(near(c.negative.singularDistance.exponent, 0.5, 0.05), "x+y, k<0: they approach P like |k|^(1/2)");
    check(c.negative.movingSingularCorank == 1, "x+y, k<0: the moving singular points are cusps (Hessian corank 1)");
    check(c.positive.singular == iso_classify::SingularBehavior::Regularized, "x+y, k>0: no singular points: " + c.positive.verdict);
}

void testPersistentSingularity()
{
    // G = x^2 vanishes to second order: P stays singular for every k.
    const auto c = classifyFamily("x^3+y^3+z^3", "x^2", 0.05);
    check(c.point == iso_classify::PointRelation::OnSurfaceSingular, "x^2: P stays singular");
    check(c.negative.singular == iso_classify::SingularBehavior::PersistentAtProbe ||
          c.negative.singular == iso_classify::SingularBehavior::PersistentAndMoving, "x^2, k<0: persistent at P: " + c.negative.verdict);
    check(c.positive.hessianCorankAtP == 2, "x^2: Hessian corank drops 3 -> 2 (E6-tilde simplifies to D4)");
}

void testTrivialDeformation()
{
    const auto c = classifyFamily("x^3+y^3+z^3", "x^3+y^3+z^3", 0.01);
    check(c.negative.surface == iso_classify::SurfaceClass::Identical &&
          c.positive.surface == iso_classify::SurfaceClass::Identical, "G = F: S_k = S_0 (identical)");
}

void testOneSidedFamily()
{
    // x^2+y^2+z^2 + k: S_k is empty for k > 0 (no real points).
    const auto c = classifyFamily("x^2+y^2+z^2", "1", 0.01);
    check(c.point == iso_classify::PointRelation::OffSurface, "definite A1 + k: P leaves S_k");
    check(c.positive.surface == iso_classify::SurfaceClass::NoLocalSurface, "definite A1, k>0: no local surface: " + c.positive.verdict);
}

void testLogSampling()
{
    const auto k = iso_classify::makeSymmetricLogKValues(0.01, 5, 3.0);
    check(k.size() == 10, "5 values per side");
    check(near(k.front(), -0.01, 1e-15) && near(k.back(), 0.01, 1e-15), "range ends at +-kMax");
    check(near(k[4], -1e-5, 1e-18) && near(k[5], 1e-5, 1e-18), "3 decades below kMax next to zero");
}
} // namespace

int main()
{
    testPredictions();
    testLogSampling();
    testRegularizingLinearPerturbation();
    testMovingSingularPoints();
    testPersistentSingularity();
    testTrivialDeformation();
    testOneSidedFamily();
    if (failures)
    {
        std::cerr << failures << " convergence classifier test(s) failed.\n";
        return 1;
    }
    std::cout << "Convergence classifier tests passed.\n";
    return 0;
}
