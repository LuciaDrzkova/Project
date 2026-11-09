// Headless k -> 0 convergence classification for one F and a list of G:
//
//   isosurface_convergence_classification --F "x^3+y^3+z^3" --G "1;x;x+y;x^2;x*y*z"
//       [--kmax 0.01] [--per-side 9] [--decades 2] [--resolution 96]
//       [--half-extent 0.4] [--out classification.csv]
//
// k is sampled log-spaced on both sides; every family is analysed exactly as by
// the viewer's "Classify convergence" button. Prints one line per G and side.
#include "analysis/ConvergenceClassifier.h"

#include <cstdio>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct Arguments {
    std::string f = "x^3+y^3+z^3";
    std::vector<std::string> g{"1", "x", "x+y", "x+y+z", "x^2", "x*y", "x^3", "x*y*z"};
    double kMax = 0.01;
    int perSide = 9;
    double decades = 2.0;
    int resolution = 96;
    double halfExtent = 0.4;
    std::string out;
};

std::vector<std::string> split(const std::string& text, char separator)
{
    std::vector<std::string> parts;
    std::stringstream in(text);
    for (std::string part; std::getline(in, part, separator);)
        if (!part.empty()) parts.push_back(part);
    return parts;
}

Arguments parse(int argc, char** argv)
{
    Arguments a;
    for (int i = 1; i + 1 < argc; i += 2)
    {
        const std::string key = argv[i], value = argv[i + 1];
        if (key == "--F") a.f = value;
        else if (key == "--G") a.g = split(value, ';');
        else if (key == "--kmax") a.kMax = std::stod(value);
        else if (key == "--per-side") a.perSide = std::stoi(value);
        else if (key == "--decades") a.decades = std::stod(value);
        else if (key == "--resolution") a.resolution = std::stoi(value);
        else if (key == "--half-extent") a.halfExtent = std::stod(value);
        else if (key == "--out") a.out = value;
        else throw std::invalid_argument("unknown option " + key);
    }
    return a;
}

// "0.5" for a single predicted exponent, "0.33-0.5" for a range.
std::string predicted(const iso_classify::Prediction& p)
{
    auto r = [](double v) { return iso_csv::shortNumber(std::round(v * 100.0) / 100.0); };
    if (!std::isfinite(p.exponent)) return "-";
    return std::fabs(p.exponentMax - p.exponent) < 1e-9 ? r(p.exponent) : r(p.exponent) + "-" + r(p.exponentMax);
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        const Arguments args = parse(argc, argv);
        std::ofstream csv;
        if (!args.out.empty())
        {
            csv.open(args.out);
            csv << iso_classify::csvHeader() << '\n';
        }
        std::printf("F = %s, |k| in [%g, %g] (%d per side, log), resolution %d, domain [-%g, %g]^3\n",
                    args.f.c_str(), args.kMax * std::pow(10.0, -args.decades), args.kMax, args.perSide,
                    args.resolution, args.halfExtent, args.halfExtent);
        std::printf("%-14s %-18s %-12s %-9s | %-20s %6s | %-20s %6s | %-22s %-4s %6s %s\n", "G", "point", "regime",
                    "pred", "surface k<0", "alpha", "surface k>0", "alpha", "singular k<0 / k>0", "n", "gamma",
                    "corank");
        for (const std::string& g : args.g)
        {
            const std::string family = "(" + args.f + ")+k*(" + g + ")";
            const NodePtr ast = parseExpression(family);
            iso_convergence::Options o;
            o.resolution = args.resolution;
            o.halfExtent = args.halfExtent;
            o.localHalfSize = 0.45 * args.halfExtent;
            o.kValues = iso_classify::makeSymmetricLogKValues(args.kMax, args.perSide, args.decades);
            o.familyExpression = family;
            const iso_classify::Classification c =
                iso_classify::classify(ast.get(), iso_research::analyzeFamily(ast.get(), o, 3));
            const auto& n = c.negative;
            const auto& p = c.positive;
            const iso_classify::SideResult& moving = n.maxSingularPoints >= p.maxSingularPoints ? n : p;
            std::printf("%-14s %-18s %-12s %-9s | %-20s %6s | %-20s %6s | %-11s/%-10s %-4zu %6s %d\n", g.c_str(),
                        iso_classify::toString(c.point), c.prediction.regime.c_str(),
                        predicted(c.prediction).c_str(),
                        iso_classify::toString(n.surface),
                        iso_csv::shortNumber(std::round(n.hausdorff.exponent * 100) / 100).c_str(),
                        iso_classify::toString(p.surface),
                        iso_csv::shortNumber(std::round(p.hausdorff.exponent * 100) / 100).c_str(),
                        iso_classify::toString(n.singular), iso_classify::toString(p.singular),
                        moving.maxSingularPoints,
                        iso_csv::shortNumber(std::round(moving.singularDistance.exponent * 100) / 100).c_str(),
                        moving.movingSingularCorank);
            std::fflush(stdout);
            if (csv) iso_classify::writeCsvRows(csv, c);
        }
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "classification failed: " << error.what() << '\n';
        return 1;
    }
}
