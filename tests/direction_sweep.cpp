// Direction sweep of linear perturbations G = c . (p - P) for one F:
//
//   isosurface_direction_sweep --F "x^3+y^3+z^3" [--k -0.01] [--longitudes 72]
//       [--heights 36] [--box 0.18] [--out sweep]
//
// Writes <out>_map.csv (signature of every grid direction) and
// <out>_boundary.csv (discriminant directions located by bisection), and
// prints the Milnor/Tjurina numbers and the number of boundary points.
#include "analysis/DirectionSweep.h"

#include <cmath>
#include <cstdio>
#include <iostream>
#include <string>

int main(int argc, char** argv)
{
    try
    {
        std::string f = "x^3+y^3+z^3", out = "sweep";
        double k = -0.01, box = 0.18;
        int longitudes = 72, heights = 36;
        for (int i = 1; i + 1 < argc; i += 2) {
            const std::string key = argv[i], value = argv[i + 1];
            if (key == "--F") f = value;
            else if (key == "--k") k = std::stod(value);
            else if (key == "--longitudes") longitudes = std::stoi(value);
            else if (key == "--heights") heights = std::stoi(value);
            else if (key == "--box") box = std::stod(value);
            else if (key == "--out") out = value;
            else throw std::invalid_argument("unknown option " + key);
        }
        const NodePtr ast = parseExpression(f);
        const auto cubic = iso_poly::Cubic::fromExpression(ast.get());
        if (!cubic) throw std::invalid_argument("F must be a polynomial of degree <= 3");
        const auto milnor = iso_algebra::localAlgebra(*cubic, iso_algebra::IdealKind::Milnor);
        const auto tjurina = iso_algebra::localAlgebra(*cubic, iso_algebra::IdealKind::Tjurina);
        if (!milnor.valid()) throw std::invalid_argument("the origin is not an isolated singular point of F");
        const iso_sweep::Sweep s = iso_sweep::sweep(*cubic, milnor.dimension, k, longitudes, heights, std::sqrt(3.0) * box);
        iso_sweep::writeCsv(out + "_map.csv", out + "_boundary.csv", s);
        std::size_t singular = 0, certified = 0;
        for (const auto& b : s.boundary) singular += b.singularFound ? 1 : 0;
        for (const auto& c : s.cells) certified += c.signature.certified ? 1 : 0;
        std::printf("F = %s: mu = %d, tau = %d (basis %s)\n", f.c_str(), milnor.dimension, tjurina.dimension,
                    tjurina.basisText().c_str());
        std::printf("%zu directions (%zu certified), %zu boundary points, %zu with a singular point of S_k\n",
                    s.cells.size(), certified, s.boundary.size(), singular);
        return 0;
    }
    catch (const std::exception& e)
    {
        std::cerr << "direction sweep failed: " << e.what() << '\n';
        return 1;
    }
}
