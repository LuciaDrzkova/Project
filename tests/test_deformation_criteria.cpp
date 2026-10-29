#include "Expr.h"

#include <cmath>
#include <cstdio>

static int failures = 0;

#define CHECK(condition, message)                                    \
    do {                                                             \
        if (!(condition)) {                                          \
            std::fprintf(stderr, "FAIL: %s\n", message);             \
            ++failures;                                              \
        }                                                            \
    } while (false)

static Dual3 evaluate(
    const Node* ast,
    double x, double y, double z,
    double k)
{
    return evalExpr<Dual3>(
        ast,
        Dual3(x, 1.0, 0.0, 0.0),
        Dual3(y, 0.0, 1.0, 0.0),
        Dual3(z, 0.0, 0.0, 1.0),
        Dual3(k, 0.0, 0.0, 0.0));
}

static double gradNorm(const Dual3& q)
{
    return std::sqrt(
        q.dx * q.dx +
        q.dy * q.dy +
        q.dz * q.dz);
}

// Singularity requires BOTH the zero-level equation and
// a vanishing spatial gradient.
static bool isSingular(
    const Dual3& q,
    double tolerance = 1e-10)
{
    return std::abs(q.v) <= tolerance &&
           gradNorm(q) <= tolerance;
}

int main()
{
    // TEST 1:
    // Cubic ordinary double point:
    //
    // F = x^2 + y^2 - z^2 + x^3
    // f_k = F + k.
    //
    // The origin is singular at k=0. Nearby surface points
    // approach the origin from BOTH signs of k, although
    // the nearby critical point is no longer on the zero level.

    {
        auto ast = parseExpression(
            "x*x + y*y - z*z + x*x*x + k");

        CHECK(
            isSingular(evaluate(ast.get(), 0, 0, 0, 0)),
            "The origin must be singular at k=0.");

        const double tValues[] = {1e-2, 1e-3, 1e-4};

        for (double t : tValues) {
            // Exact zero-level points for k<0:
            // p(t) = (t,0,0), k = -(t^2+t^3).
            const double kMinus = -(t*t + t*t*t);
            const Dual3 minusPoint =
                evaluate(ast.get(), t, 0, 0, kMinus);

            CHECK(
                std::abs(minusPoint.v) < 1e-12,
                "Negative-k point must lie on the surface.");

            CHECK(
                gradNorm(minusPoint) > 1e-10,
                "Negative-k point must be regular.");

            // Exact zero-level points for k>0:
            // p(t) = (0,0,t), k = t^2.
            const double kPlus = t*t;
            const Dual3 plusPoint =
                evaluate(ast.get(), 0, 0, t, kPlus);

            CHECK(
                std::abs(plusPoint.v) < 1e-12,
                "Positive-k point must lie on the surface.");

            CHECK(
                gradNorm(plusPoint) > 1e-10,
                "Positive-k point must be regular.");
        }

        const Dual3 originPlus =
            evaluate(ast.get(), 0, 0, 0, 0.1);

        CHECK(
            std::abs(originPlus.v) > 1e-3,
            "For k!=0, the critical origin is off the zero level.");

        CHECK(
            !isSingular(originPlus),
            "An off-level critical point is not a surface singularity.");
    }

    // TEST 2:
    // An isolated cubic cusp whose singularity is removed:
    //
    // f_k = x^3 + y^2 + z^2 + k*x.
    //
    // The origin remains on every surface, but is regular
    // whenever k!=0.

    {
        auto ast = parseExpression(
            "x*x*x + y*y + z*z + k*x");

        CHECK(
            isSingular(evaluate(ast.get(), 0, 0, 0, 0)),
            "The cusp must be singular at k=0.");

        const double kValues[] = {-0.04, 0.04};

        for (double k : kValues) {
            const Dual3 origin =
                evaluate(ast.get(), 0, 0, 0, k);

            CHECK(
                std::abs(origin.v) < 1e-12,
                "The origin remains on the cusp family.");

            CHECK(
                std::abs(gradNorm(origin) - std::abs(k)) < 1e-12,
                "The gradient norm at the origin must equal |k|.");

            CHECK(
                !isSingular(origin),
                "The origin must be regular for nonzero k.");

            if (k < 0.0) {
                // Critical points exist at x=+/-sqrt(-k/3).
                // Their critical values must be nonzero.
                const double a = std::sqrt(-k / 3.0);

                for (double x : {a, -a}) {
                    const Dual3 critical =
                        evaluate(ast.get(), x, 0, 0, k);

                    CHECK(
                        gradNorm(critical) < 1e-10,
                        "The candidate must be a critical point.");

                    CHECK(
                        std::abs(critical.v) > 1e-8,
                        "The critical point must be off the zero level.");
                }
            }

            // For k>0, 3*x^2+k cannot be zero over R.
            // Thus no real critical points exist for this family.
        }
    }

    // TEST 3:
    // Persistent singularity in a cubic family:
    //
    // f_k = x^3 + y^3 + z^3 - 3*k*x*y*z.
    //
    // The origin satisfies both singularity equations for
    // every k. This is a regression test for persistence,
    // not a general singularity solver.

    {
        auto ast = parseExpression(
            "x*x*x + y*y*y + z*z*z - 3*k*x*y*z");

        const double kValues[] = {-0.5, 0.0, 0.5};

        for (double k : kValues) {
            CHECK(
                isSingular(evaluate(ast.get(), 0, 0, 0, k)),
                "The origin must remain singular for this family.");
        }
    }

    // TEST 4:
    // Intersection of the +k and -k surfaces.
    //
    // F = x^2 + y^2 - z^2 + x^3
    // G = x.
    //
    // Every point (0,t,t) belongs to both surfaces because
    // F=G=0 there, independently of the nonzero value of k.

    {
        auto ast = parseExpression(
            "x*x + y*y - z*z + x*x*x + k*x");

        const double t = 0.1;
        const double kValues[] = {-0.2, 0.2};

        for (double k : kValues) {
            const Dual3 plus =
                evaluate(ast.get(), 0, t, t, k);

            const Dual3 minus =
                evaluate(ast.get(), 0, t, t, -k);

            CHECK(
                std::abs(plus.v) < 1e-12,
                "The common point must belong to X_k.");

            CHECK(
                std::abs(minus.v) < 1e-12,
                "The common point must belong to X_-k.");

            CHECK(
                gradNorm(plus) > 1e-10 &&
                gradNorm(minus) > 1e-10,
                "The chosen intersection point must be regular.");
        }
    }

    if (failures != 0) {
        std::fprintf(stderr, "%d test(s) failed.\n", failures);
        return 1;
    }

    std::printf("All deformation-criteria tests passed.\n");
    return 0;
}