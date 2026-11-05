// Slice snapshot: validation, naming, and the SVG/CSV files it writes.
#include "core/Expr.h"
#include "core/IsosurfaceSlices.h"
#include "core/SliceSnapshot.h"

#include <filesystem>
#include <fstream>
#include <iostream>
#include <sstream>
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

std::string read_file(const std::filesystem::path& path)
{
    std::ifstream in(path);
    std::stringstream buffer;
    buffer << in.rdbuf();
    return buffer.str();
}
} // namespace

int main(int argc, char** argv)
{
    namespace fs = std::filesystem;
    const fs::path out = (argc == 3 && std::string(argv[1]) == "--out") ? fs::path(argv[2])
                                                                       : fs::path("slice_snapshot_test_output");
    fs::remove_all(out);

    // Four cells -> grid coordinates {-1, -0.5, 0, 0.5, 1}; Z-normal slices.
    const std::string expression = "x*x+y*y-0.3+k*z";
    const NodePtr ast = parseExpression(expression);
    const iso_slices::Result result =
        iso_slices::extract_and_export(ast.get(), expression, 0.1, 4, 1.0, 2, (out / "slices").string());
    check(result.expression == expression, "slice result records its expression");

    iso_slices::SnapshotRequest request;
    request.slices = {3, 1, 2}; // order does not matter
    request.iu_min = 1;
    request.iu_max = 3;
    request.iv_min = 0;
    request.iv_max = 2;
    check(iso_slices::validate_snapshot(result, request).empty(), "valid request is accepted");
    check(iso_slices::slice_list_label({1, 2, 3, 6}) == "2-4_7", "slice label collapses consecutive runs");
    check(iso_slices::snapshot_base_name(result, request) == "xxx+yxy-0.3+kxz_slices_2-4_u1-3_v0-2",
          "base name encodes expression, slices and grid window: " + iso_slices::snapshot_base_name(result, request));

    auto rejected = [&](iso_slices::SnapshotRequest bad) { return !iso_slices::validate_snapshot(result, bad).empty(); };
    iso_slices::SnapshotRequest bad = request;
    bad.slices = {0, 1};
    check(rejected(bad), "fewer than three slices are rejected");
    bad = request;
    bad.slices = {0, 1, 1};
    check(rejected(bad), "duplicate slices are rejected");
    bad = request;
    bad.slices = {0, 1, 5};
    check(rejected(bad), "slices outside the result are rejected");
    bad = request;
    bad.iu_max = bad.iu_min;
    check(rejected(bad), "a degenerate (non-rectangular) window is rejected");
    bad = request;
    bad.iv_max = 5;
    check(rejected(bad), "a window outside the grid is rejected");

    const iso_slices::SnapshotFiles files = iso_slices::write_snapshot(result, request, out / "snapshots");
    check(fs::exists(files.svg) && fs::exists(files.csv), "snapshot writes an SVG and a CSV");

    const std::string csv = read_file(files.csv);
    std::size_t lines = 0;
    for (char c : csv) lines += (c == '\n') ? 1 : 0;
    check(lines == 1 + 3 * 3 * 3, "CSV has a header plus 3 slices x 3x3 nodes");
    check(csv.find("slice_number,slice_index,axis,plane_coordinate,iu,iv,grid_ix,grid_iy,grid_iz,x,y,z,F") == 0,
          "CSV header lists indices, coordinates and F");

    const std::string svg = read_file(files.svg);
    check(svg.find("<svg") != std::string::npos && svg.find("</svg>") != std::string::npos, "SVG is complete");
    check(svg.find("Slice 2 ") != std::string::npos && svg.find("Slice 4 ") != std::string::npos,
          "SVG has one panel per selected slice (1-based numbers)");
    // Slice 2 (index 1) is the plane z = -0.5; corners are nodes (iu,iv) = (1,0) and (3,2).
    check(svg.find("bottom-left [1,0,1] (-0.5, -1, -0.5)") != std::string::npos, "bottom-left corner is labelled");
    check(svg.find("top-right [3,2,1] (0.5, 0, -0.5)") != std::string::npos, "top-right corner is labelled");

    const iso_slices::SnapshotFiles second = iso_slices::write_snapshot(result, request, out / "snapshots");
    check(second.svg != files.svg && second.svg.filename().string().find("_02.svg") != std::string::npos,
          "a repeated snapshot never overwrites an earlier one");

    if (failures)
    {
        std::cerr << failures << " slice snapshot test(s) failed.\n";
        return 1;
    }
    std::cout << "Slice snapshot tests passed.\n";
    return 0;
}
