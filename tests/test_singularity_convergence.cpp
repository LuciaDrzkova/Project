// Headless regression tests for the k -> 0 convergence diagnostics.
// These tests verify numerical trends for selected cubic examples; they do not
// prove a general convergence theorem or certify isolation of arbitrary germs.
#include "SingularityConvergence.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cstdio>
#include <string>

namespace {
int failures = 0;

void check(bool condition, const std::string& message) {
    if (!condition) {
        std::cerr << "FAIL: " << message << '\n';
        ++failures;
    }
}

bool near(double a, double b, double tolerance) {
    return std::isfinite(a) && std::isfinite(b) && std::fabs(a - b) <= tolerance;
}

void testIndefiniteCubicApproachesFromBothSides() {
    // Degree-three family F_k = x^2 + y^2 - z^2 + x^3 + k.
    // At k=0 the origin is a singular point. Along the z-axis the positive
    // side has points at z= +/-sqrt(k); on the negative side the x-axis gives
    // roots approaching 0. The local geometry is a cone-like singular fibre.
    NodePtr ast = parseExpression("x*x + y*y - z*z + x*x*x + k");
    const auto p0 = iso_convergence::evaluate(ast.get(), {0,0,0}, 0.0);
    check(near(p0.value, 0.0, 1e-12), "indefinite cubic: F(0,0,0,0)=0");
    check(p0.gradientNorm < 1e-12, "indefinite cubic: origin gradient vanishes at k=0");

    iso_convergence::Options options;
    options.probe = {0,0,0};
    options.resolution = 64;
    options.halfExtent = 0.6;
    options.localRadius = 0.42;
    options.magnitudes = {0.09, 0.04, 0.01, 0.0025};
    const auto a = iso_convergence::analyze(ast.get(), options);
    check(a.samples.size() == 8, "both signs are analysed for every magnitude");

    check(std::isfinite(a.negativeSide.pointDistance.logLogSlope) &&
          a.negativeSide.pointDistance.logLogSlope > 0.1,
          "negative-k point-distance samples should show a positive empirical log-log slope");
    check(std::isfinite(a.positiveSide.pointDistance.logLogSlope) &&
          a.positiveSide.pointDistance.logLogSlope > 0.1,
          "positive-k point-distance samples should show a positive empirical log-log slope");

    for (int sign : {-1, +1}) {
        double previous = std::numeric_limits<double>::infinity();
        int count = 0;
        for (const auto& s : a.samples) {
            if ((s.k < 0.0 ? -1 : 1) != sign) continue;
            check(s.hasSurfaceWithinLocalRadius,
                  "indefinite cubic should have mesh near the origin at k=" + std::to_string(s.k));
            if (std::isfinite(s.distanceProbeToMesh)) {
                // Values are coarse-grid estimates; only assert a broad trend.
                if (count > 0)
                    check(s.distanceProbeToMesh <= previous + 0.08,
                          "probe-to-mesh distance should not jump dramatically as |k| shrinks");
                previous = s.distanceProbeToMesh;
            }
            ++count;
        }
        check(count == 4, "each sign has four samples");
    }
}

void testDefiniteCubicIsOneSidedInSmallBox() {
    // F_k = x^2+y^2+z^2+x^3+k. For |x|<0.6, x^2+x^3 >= 0.4 x^2;
    // hence no nearby real zero exists for positive k, while negative k gives
    // a small real surface surrounding the origin.
    NodePtr ast = parseExpression("x*x + y*y + z*z + x*x*x + k");
    const auto p0 = iso_convergence::evaluate(ast.get(), {0,0,0}, 0.0);
    check(near(p0.value, 0.0, 1e-12), "definite cubic: origin is on S_0");
    check(p0.gradientNorm < 1e-12, "definite cubic: origin is singular at k=0");

    iso_convergence::Options options;
    options.probe = {0,0,0};
    options.resolution = 64;
    options.halfExtent = 0.6;
    options.localRadius = 0.42;
    options.magnitudes = {0.09, 0.04, 0.01};
    const auto a = iso_convergence::analyze(ast.get(), options);
    for (const auto& s : a.samples) {
        if (s.k > 0.0) {
            check(!s.hasSurfaceWithinLocalRadius,
                  "positive k should have no local real surface for the definite cubic");
        } else {
            check(s.hasSurfaceWithinLocalRadius,
                  "negative k should create a local real surface for the definite cubic");
        }
    }
}


void testSelectedRangeAndCount() {
    const auto values = iso_convergence::makeLinearKValues(-0.2, 0.2, 5);
    check(values.size() == 5, "selected range generates exactly the requested number of k samples");
    check(near(values.front(), -0.2, 1e-14), "selected range includes its lower endpoint");
    check(near(values.back(), 0.2, 1e-14), "selected range includes its upper endpoint");
    check(near(values[2], 0.0, 1e-14), "odd symmetric sample count includes k=0");

    NodePtr ast = parseExpression("x*x + y*y + z*z + k");
    iso_convergence::Options options;
    options.probe = {0,0,0};
    options.resolution = 24;
    options.halfExtent = 0.6;
    options.localRadius = 0.4;
    options.kValues = values;
    const auto a = iso_convergence::analyze(ast.get(), options);
    check(a.samples.size() == 5, "analysis extracts exactly the selected number of k samples");
    int zeroCount = 0;
    for (const auto& sample : a.samples)
        if (sample.k == 0.0) ++zeroCount;
    check(zeroCount == 1, "k=0 is represented once when it is in the selected range");
}

void testKRangeValidation() {
    bool threw = false;
    try { (void)iso_convergence::makeLinearKValues(1.0, -1.0, 5); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "reversed k range is rejected");

    threw = false;
    try { (void)iso_convergence::makeLinearKValues(-1.0, 1.0, 2); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "sample count below 3 is rejected");
}


void testCsvContainsReproducibilityMetadata() {
    const std::string family = "x*x + y*y - z*z + x*x*x + k*(x+y)";
    NodePtr ast = parseExpression(family);
    iso_convergence::Options options;
    options.probe = {0.125, -0.25, 0.0};
    options.resolution = 16;
    options.halfExtent = 0.8;
    options.localRadius = 0.3;
    options.kValues = {-0.02, 0.0, 0.02};
    options.familyExpression = family;

    const auto analysis = iso_convergence::analyze(ast.get(), options);
    const auto& meta = analysis.metadata;
    check(meta.familyExpression == family, "metadata preserves the exact family expression text");
    check(meta.kAffineStatus == "STRUCTURALLY_AFFINE_IN_K", "family expression is identified as structurally affine in k");
    check(meta.degreeF == "3", "degree of F(x,y,z) is recorded as total degree three");
    check(meta.degreeG == "1", "degree of G(x,y,z) is recorded as total degree one");
    check(meta.gridResolutionCellsPerAxis == 16 && meta.gridSamplesPerAxis == 17,
          "grid cell resolution and sample count are recorded");
    check(near(meta.domainXMin, -0.8, 1e-12) && near(meta.domainXMax, 0.8, 1e-12),
          "domain bounds are recorded from the half extent");
    check(near(meta.gridSpacingX, 0.1, 1e-12), "grid spacing is recorded");
    check(near(meta.localRadius, 0.3, 1e-12), "local comparison radius is recorded");
    check(near(meta.probe.x, 0.125, 1e-12) && near(meta.probe.y, -0.25, 1e-12),
          "probe coordinates are recorded");
    check(meta.kSampleCount == 3 && meta.negativeKSampleCount == 1 &&
          meta.zeroKSampleCount == 1 && meta.positiveKSampleCount == 1,
          "sample count and positive/zero/negative counts are recorded");
    check(meta.referenceK == 0.0 && meta.referenceTriangles == meta.referenceTotalVertices / 3,
          "reference fibre and reference mesh sizes are recorded");
    check(meta.isosurfaceExtractor.find("Marching Cubes") != std::string::npos,
          "isosurface extraction algorithm is documented");
    check(meta.localHausdorffDistance.find("not exact") != std::string::npos,
          "the sampled-vertex Hausdorff limitation is documented");
    check(meta.kValuesUniformlySpaced && near(meta.kSpacingIfUniform, 0.02, 1e-12),
          "uniform k sample spacing is recorded");
    check(near(meta.exactZeroBiasFraction, 0.05, 1e-12) &&
          meta.degenerateTriangleHandling.find("step^2") != std::string::npos,
          "exact-zero tie-break and degenerate-triangle policy are documented");

    const std::string path = "singularity_convergence_metadata_test.csv";
    std::string error;
    check(iso_convergence::writeCsv(path, analysis, &error), "metadata CSV is written successfully: " + error);
    std::ifstream input(path);
    std::stringstream buffer;
    buffer << input.rdbuf();
    const std::string csv = buffer.str();
    check(csv.find("metadata_family_expression") != std::string::npos,
          "CSV header contains family-expression metadata");
    check(csv.find("metadata_probe_x") != std::string::npos &&
          csv.find("metadata_domain_z_max") != std::string::npos,
          "CSV header contains probe and domain metadata");
    check(csv.find("metadata_reference_triangles") != std::string::npos,
          "CSV header contains reference-mesh metadata");
    check(csv.find("exact_zero_grid_samples") != std::string::npos &&
          csv.find("skipped_degenerate_triangles") != std::string::npos &&
          csv.find("mesh_build_milliseconds") != std::string::npos,
          "CSV header contains per-mesh extraction diagnostics");
    check(csv.find("metadata_exact_zero_grid_handling") != std::string::npos &&
          csv.find("metadata_degenerate_triangle_handling") != std::string::npos,
          "CSV header contains marching-cubes tie-break and filtering metadata");
    check(csv.find("metadata_compiler") != std::string::npos &&
          csv.find("metadata_double_epsilon") != std::string::npos,
          "CSV header contains software and precision metadata");
    check(csv.find("STRUCTURALLY_AFFINE_IN_K") != std::string::npos,
          "CSV rows contain affine-decomposition status");
    check(csv.find("x*x + y*y - z*z + x*x*x + k*(x+y)") != std::string::npos,
          "CSV rows contain the original polynomial family text");
    std::remove(path.c_str());
}


void testNonAffineFamilyIsNotMisrepresentedAsFPlusKG() {
    const std::string family = "x*x + k*k";
    NodePtr ast = parseExpression(family);
    iso_convergence::Options options;
    options.resolution = 12;
    options.halfExtent = 0.5;
    options.localRadius = 0.25;
    options.kValues = {-0.1, 0.1};
    options.familyExpression = family;
    const auto analysis = iso_convergence::analyze(ast.get(), options);
    check(analysis.metadata.kAffineStatus.find("NOT_VERIFIED") != std::string::npos,
          "a nonlinear k dependence is not silently classified as an affine unfolding");
    check(analysis.metadata.gExpression.find("NOT_AVAILABLE") != std::string::npos,
          "G is not reported as an exact coefficient for a non-affine family");
}

void testInputValidation() {
    NodePtr ast = parseExpression("x*x+y*y-z*z+k");
    iso_convergence::Options options;
    options.resolution = 4;
    options.magnitudes = {0.1, 0.0};
    bool threw = false;
    try { (void)iso_convergence::analyze(ast.get(), options); }
    catch (const std::invalid_argument&) { threw = true; }
    check(threw, "zero k magnitude is rejected rather than entering log-log fits");
}
} // namespace

int main() {
    testIndefiniteCubicApproachesFromBothSides();
    testDefiniteCubicIsOneSidedInSmallBox();
    testSelectedRangeAndCount();
    testKRangeValidation();
    testCsvContainsReproducibilityMetadata();
    testNonAffineFamilyIsNotMisrepresentedAsFPlusKG();
    testInputValidation();
    if (failures != 0) {
        std::cerr << failures << " convergence test(s) failed.\n";
        return EXIT_FAILURE;
    }
    std::cout << "Singularity convergence regression tests passed.\n";
    return EXIT_SUCCESS;
}
