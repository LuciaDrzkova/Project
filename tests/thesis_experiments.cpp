// Reproduces the complete experiment set of the thesis (chapter 5):
//
//   * ten reference surfaces F with an isolated singular point at the origin
//     (types A1, A2, D4 and the cubic cones E6~);
//   * for each F the common perturbations G, the balanced linear
//     perturbation (cubic cones only, located on the discriminant by the
//     direction sweep) and G = F;
//   * grid resolutions 32, 64 and 128 on the domain [-0.4, 0.4]^3, local box
//     [-0.18, 0.18]^3, |k| log-spaced in [1e-4, 1e-2] (9 values per side);
//   * for each F the direction sweep of linear perturbations at k = -0.01.
//
//   isosurface_thesis_experiments [output-directory] [--quick]
//
// Output: <dir>/<surface>/res<N>/<G>/  (all CSV reports of the run)
//         <dir>/<surface>/direction_sweep_{map,boundary}.csv
//         <dir>/classification_all.csv  (one table, two rows per run)
//         <dir>/surfaces.csv            (Milnor/Tjurina numbers, balanced G)
#include "analysis/ConvergenceClassifier.h"
#include "analysis/DirectionSweep.h"
#include "core/FileNames.h"

#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

struct Surface {
    std::string folder; // <type>_<F>
    std::string type;
    std::string f;
};

const std::vector<Surface>& surfaces()
{
    static const std::vector<Surface> list = {
        {"A1_x2+y2-z2+xyz", "A1", "x^2+y^2-z^2+x*y*z"},
        {"A2_x2+y2+z3", "A2", "x^2+y^2+z^3"},
        {"A2_x2-y2+z3", "A2", "x^2-y^2+z^3"},
        {"A2_x3+y3+z2+x2", "A2", "x^3+y^3+z^2+x^2"},
        {"D4_x2+y3+z3", "D4", "x^2+y^3+z^3"},
        {"D4_x2+y3-3yz2", "D4", "x^2+y^3-3*y*z^2"},
        {"E6_x3+y3+z3", "E6~", "x^3+y^3+z^3"},
        {"E6_x3+y3-z3", "E6~", "x^3+y^3-z^3"},
        {"E6_x3+y3+z3+xyz", "E6~", "x^3+y^3+z^3+x*y*z"},
        {"E6_x3+2y3+3z3+2xyz", "E6~", "x^3+2*y^3+3*z^3+2*x*y*z"},
    };
    return list;
}

const std::vector<std::string> kCommonG = {
    "1", "x", "y", "z", "x+y", "x+z", "x+y+z", "2*x-y+3*z",
    "x^2", "x^2+y^2-z^2", "x*y", "x^3", "x^2*y+y^2*z+z^2*x"};

constexpr double kHalfExtent = 0.4;
constexpr double kLocalHalfSize = 0.18;

iso_convergence::Options options(const std::string& family, int resolution, const std::vector<double>& k)
{
    iso_convergence::Options o;
    o.probe = {0, 0, 0};
    o.resolution = resolution;
    o.halfExtent = kHalfExtent;
    o.localHalfSize = kLocalHalfSize;
    o.kValues = k;
    o.familyExpression = family;
    return o;
}

// Balanced linear perturbation of a cubic cone: a discriminant direction on
// the meridian through (1,1,0), written as x + b*y + c*z.
std::string balancedG(const iso_poly::Cubic& f, int mu, std::ostream& log)
{
    const auto boundary = iso_sweep::meridianBoundary(f, mu, -0.01, 45.0, std::sqrt(3.0) * kLocalHalfSize);
    const iso_sweep::BoundaryPoint* best = nullptr;
    for (const auto& b : boundary) {
        if (!b.singularFound) continue;
        const double score = b.c[0] + b.c[1] + b.c[2];
        if (!best || score > best->c[0] + best->c[1] + best->c[2] + 1e-12) best = &b;
    }
    if (!best) {
        log << "  no discriminant direction on the meridian through (1,1,0)\n";
        return {};
    }
    // On this meridian c_x = c_y, so G = x + y + t z with t = c_z / c_x
    // (located to double precision by bisection and secant refinement).
    const double t = best->c[2] / best->c[0];
    std::ostringstream g;
    g << "x+y" << (t < 0 ? "-" : "+") << std::setprecision(17) << std::fabs(t) << "*z";
    return g.str();
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        fs::path root = "experiments_thesis";
        bool quick = false;
        for (int i = 1; i < argc; ++i) {
            const std::string a = argv[i];
            if (a == "--quick") quick = true;
            else root = a;
        }
        fs::create_directories(root);
        const std::vector<int> resolutions = quick ? std::vector<int>{24} : std::vector<int>{32, 64, 128};
        const std::vector<double> kValues = iso_classify::makeSymmetricLogKValues(0.01, quick ? 4 : 9, 2.0);

        std::ofstream all(root / "classification_all.csv");
        all << "surface,type,resolution," << iso_classify::csvHeader() << '\n';
        std::ofstream summary(root / "surfaces.csv");
        summary << "surface,type,F,milnor_number,tjurina_number,tjurina_basis,balanced_G,sweep_boundary_points,"
                   "sweep_singular_boundary_points\n";

        for (const Surface& surface : surfaces()) {
            const fs::path dir = root / surface.folder;
            fs::create_directories(dir);
            std::printf("== %s (%s)\n", surface.f.c_str(), surface.type.c_str());
            std::fflush(stdout);

            const NodePtr fAst = parseExpression(surface.f);
            const iso_poly::Cubic f = *iso_poly::Cubic::fromExpression(fAst.get());
            const auto milnor = iso_algebra::localAlgebra(f, iso_algebra::IdealKind::Milnor);
            const auto tjurina = iso_algebra::localAlgebra(f, iso_algebra::IdealKind::Tjurina);

            // Direction sweep of linear perturbations (k < 0; k > 0 is the antipodal map).
            const iso_sweep::Sweep sweep = iso_sweep::sweep(f, milnor.dimension, -0.01, quick ? 24 : 72,
                                                            quick ? 12 : 36, std::sqrt(3.0) * kLocalHalfSize);
            iso_sweep::writeCsv((dir / "direction_sweep_map.csv").string(),
                                (dir / "direction_sweep_boundary.csv").string(), sweep);
            std::size_t singularBoundary = 0;
            for (const auto& b : sweep.boundary) singularBoundary += b.singularFound ? 1 : 0;

            std::vector<std::string> gList = kCommonG;
            std::string balanced;
            if (surface.type == "E6~") {
                balanced = balancedG(f, milnor.dimension, std::cout);
                if (!balanced.empty()) gList.push_back(balanced);
            }
            gList.push_back(surface.f); // G = F: the trivial deformation (1+k)F
            summary << iso_csv::quote(surface.folder) << ',' << iso_csv::quote(surface.type) << ','
                    << iso_csv::quote(surface.f) << ',' << milnor.dimension << ',' << tjurina.dimension << ','
                    << iso_csv::quote(tjurina.basisText()) << ',' << iso_csv::quote(balanced) << ','
                    << sweep.boundary.size() << ',' << singularBoundary << '\n';
            summary.flush();

            for (const std::string& g : gList) {
                const std::string family = "(" + surface.f + ")+k*(" + g + ")";
                const NodePtr ast = parseExpression(family);
                const std::string gFolder = g == balanced ? "balanced" : (g == surface.f ? "G_equals_F"
                                                                                         : iso_files::safeFileComponent(g, 60));
                for (int resolution : resolutions) {
                    const auto a = iso_research::analyzeFamily(ast.get(), options(family, resolution, kValues), 3);
                    const auto c = iso_classify::classify(ast.get(), a);
                    const fs::path runDir = dir / ("res" + std::to_string(resolution)) / gFolder;
                    fs::create_directories(runDir);
                    std::string error;
                    if (!iso_research::writeFamilyCsvs(runDir.string(), a, &error) ||
                        !iso_classify::writeCsv((runDir / "convergence_classification.csv").string(), c, &error))
                        throw std::runtime_error(error);
                    std::ostringstream rows;
                    iso_classify::writeCsvRows(rows, c);
                    std::istringstream lines(rows.str());
                    for (std::string line; std::getline(lines, line);)
                        all << iso_csv::quote(surface.folder) << ',' << iso_csv::quote(surface.type) << ','
                            << resolution << ',' << line << '\n';
                    all.flush();
                }
                std::printf("  G = %s done\n", g.c_str());
                std::fflush(stdout);
            }
        }
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "thesis experiments failed: " << e.what() << '\n';
        return 1;
    }
}
