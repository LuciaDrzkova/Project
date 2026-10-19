#include "Expr.h"
#include "IsosurfaceSlices.h"

#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv)
{
    try
    {
        std::string out = "slice_extraction_test_output";
        if (argc == 3 && std::string(argv[1]) == "--out")
            out = argv[2];

        namespace fs = std::filesystem;
        fs::remove_all(out);
        const NodePtr expression = parseExpression("x");
        // Four cells mean five grid coordinates {-1,-0.5,0,0.5,1}.
        // Cuts are exactly on those Z grid planes, not cell centers.
        const iso_slices::Result result = iso_slices::extract_and_export(
            expression.get(), "x", 0.0, 4, 1.0, 2, out);

        if (result.plane_count != 5)
            throw std::runtime_error("Expected resolution+1 grid-aligned slice planes.");
        if (result.sample_count != 125 || result.finite_sample_count != 125 ||
            result.samples.size() != 125)
            throw std::runtime_error("Expected 125 finite samples for 5 slices of 5x5 nodes.");
        if (result.contour_segment_count != 20)
            throw std::runtime_error("Expected 4 contour segments on each of 5 grid-aligned planes.");

        const auto at = [&](int slice, int iu, int iv) -> const iso_slices::GridSample& {
            const std::size_t side = static_cast<std::size_t>(result.resolution + 1);
            return result.samples[static_cast<std::size_t>(slice) * side * side +
                                  static_cast<std::size_t>(iv) * side +
                                  static_cast<std::size_t>(iu)];
        };
        if (std::fabs(at(0, 2, 1).point.z + 1.0) > 1e-12 ||
            std::fabs(at(2, 2, 1).point.z) > 1e-12 ||
            std::fabs(at(4, 2, 1).point.z - 1.0) > 1e-12)
            throw std::runtime_error("Slice coordinates do not land on the original grid planes.");
        if (at(2, 2, 1).ix != 2 || at(2, 2, 1).iy != 1 || at(2, 2, 1).iz != 2)
            throw std::runtime_error("A sampled 2D grid node did not preserve its 3D lattice indices.");
        if (std::fabs(at(2, 2, 1).value) > 1e-12)
            throw std::runtime_error("Expected F=0 at the selected x=0 grid node.");

        if (!fs::exists(fs::path(out) / "slice_grid_samples.csv") ||
            !fs::exists(fs::path(out) / "slice_contour_segments.csv") ||
            !fs::exists(fs::path(out) / "slice_summary.csv") ||
            !fs::exists(fs::path(out) / "slice_metadata.txt"))
            throw std::runtime_error("One or more output files are missing.");

        std::ifstream grid_csv(fs::path(out) / "slice_grid_samples.csv");
        std::string header;
        std::getline(grid_csv, header);
        if (header.find("grid_ix,grid_iy,grid_iz") == std::string::npos ||
            header.find("x,y,z,F") == std::string::npos)
            throw std::runtime_error("CSV does not include both 3D lattice indices and coordinates/isovalues.");

        std::cout << "2D grid-aligned slice extraction passed: " << result.plane_count
                  << " planes on axis coordinates, " << result.sample_count
                  << " 3D-mapped grid samples, " << result.contour_segment_count
                  << " contour segments.\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << "2D slice extraction test failed: " << error.what() << '\n';
        return 1;
    }
}
