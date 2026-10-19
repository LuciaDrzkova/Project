#pragma once

// Grid-aligned 2D cross-sections of the current implicit field F(x,y,z,k)=0.
// Each selected plane is sampled on the same grid as the 3D extractor, and
// marching squares extracts contour line segments from those samples.
#include "Expr.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace iso_slices
{
struct Point3
{
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

struct ContourSegment
{
    int slice_index = 0;
    int cell_u = 0;
    int cell_v = 0;
    double plane_coordinate = 0.0;
    // Corner order: (u,v), (u+1,v), (u+1,v+1), (u,v+1).
    std::array<double, 4> corner_values{};
    Point3 a;
    Point3 b;
};

struct GridSample
{
    int slice_index = 0;
    int iu = 0;
    int iv = 0;
    int ix = 0;
    int iy = 0;
    int iz = 0;
    Point3 point;
    double value = std::numeric_limits<double>::quiet_NaN();
};

struct Result
{
    int axis = 2; // 0=X (YZ planes), 1=Y (XZ planes), 2=Z (XY planes)
    int resolution = 0;
    double k = 0.0;
    double half_extent = 1.0;
    double spacing = 0.0;
    std::size_t plane_count = 0;
    std::size_t sample_count = 0;
    std::size_t finite_sample_count = 0;
    std::size_t invalid_sample_count = 0;
    std::size_t contour_segment_count = 0;
    std::size_t skipped_nonfinite_cells = 0;
    std::size_t unresolved_degenerate_cells = 0;
    // Flat order: slice, then iv, then iu. Every row maps to an exact 3D grid node.
    std::vector<GridSample> samples;
    std::vector<ContourSegment> segments;
};

inline const char* axis_name(int axis)
{
    switch (axis)
    {
    case 0: return "X";
    case 1: return "Y";
    default: return "Z";
    }
}

inline double grid_coord(int index, int resolution, double extent)
{
    return -extent + (2.0 * extent * static_cast<double>(index) /
                      static_cast<double>(std::max(1, resolution)));
}

inline Point3 point_from_indices(const std::array<int, 3>& index,
                                 int resolution, double extent)
{
    return {grid_coord(index[0], resolution, extent),
            grid_coord(index[1], resolution, extent),
            grid_coord(index[2], resolution, extent)};
}

inline std::array<double, 3> as_array(const Point3& p)
{
    return {p.x, p.y, p.z};
}

inline double evaluate(const Node* ast, const Point3& p, double k)
{
    return evalExpr<double>(ast, p.x, p.y, p.z, k);
}

inline Point3 interpolate_edge(const Point3& a, const Point3& b,
                               double fa, double fb)
{
    const double denom = fa - fb;
    double t = std::abs(denom) > std::numeric_limits<double>::min()
        ? fa / denom : 0.5;
    t = std::clamp(t, 0.0, 1.0);
    return {a.x + t * (b.x - a.x),
            a.y + t * (b.y - a.y),
            a.z + t * (b.z - a.z)};
}

inline double squared_distance(const Point3& a, const Point3& b)
{
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    const double dz = a.z - b.z;
    return dx * dx + dy * dy + dz * dz;
}

inline void write_point(std::ostream& out, const Point3& p)
{
    out << p.x << ',' << p.y << ',' << p.z;
}

// Extracts N+1 parallel slices exactly on the original grid planes along the
// selected axis, including both domain boundaries. Each plane is sampled using
// the existing (N+1)x(N+1) in-plane grid, then marching squares extracts F=0.
// No viewer mesh or existing viewer functionality is modified.
inline Result extract_and_export(const Node* ast,
                                 const std::string& expression,
                                 double k,
                                 int resolution,
                                 double half_extent,
                                 int axis,
                                 const std::string& output_directory)
{
    if (!ast)
        throw std::invalid_argument("Cannot slice an empty expression.");
    if (resolution < 1)
        throw std::invalid_argument("Grid resolution must be at least 1.");
    if (!(half_extent > 0.0) || !std::isfinite(half_extent))
        throw std::invalid_argument("Domain half-extent must be finite and positive.");
    if (axis < 0 || axis > 2)
        throw std::invalid_argument("Slice axis must be X (0), Y (1), or Z (2).");
    if (output_directory.empty())
        throw std::invalid_argument("Output folder must not be empty.");

    namespace fs = std::filesystem;
    const fs::path out_dir(output_directory);
    fs::create_directories(out_dir);
    std::ofstream grid_csv(out_dir / "slice_grid_samples.csv");
    std::ofstream segments_csv(out_dir / "slice_contour_segments.csv");
    std::ofstream summary_csv(out_dir / "slice_summary.csv");
    std::ofstream metadata(out_dir / "slice_metadata.txt");
    if (!grid_csv || !segments_csv || !summary_csv || !metadata)
        throw std::runtime_error("Could not create one or more slice output files in: " + output_directory);

    Result result;
    result.axis = axis;
    result.resolution = resolution;
    result.k = k;
    result.half_extent = half_extent;
    result.spacing = (2.0 * half_extent) / resolution;
    result.plane_count = static_cast<std::size_t>(resolution + 1);
    const int side = resolution + 1;
    result.sample_count = result.plane_count * static_cast<std::size_t>(side) * side;
    result.samples.reserve(result.sample_count);
    result.segments.reserve(static_cast<std::size_t>(2) * resolution * resolution * result.plane_count);

    // The slice normal is the selected axis. The other two axes are sampled
    // at their original grid vertices, so each slice is a true lattice plane.
    const int actual_u_axis = (axis == 0) ? 1 : 0;
    const int actual_v_axis = (axis == 2) ? 1 : 2;

    grid_csv << std::setprecision(17);
    segments_csv << std::setprecision(17);
    summary_csv << std::setprecision(17);
    grid_csv << "slice_index,axis,plane_coordinate,iu,iv,grid_ix,grid_iy,grid_iz,x,y,z,F,is_negative\n";
    segments_csv << "slice_index,axis,plane_coordinate,cell_u,cell_v,"
                    "F00,F10,F11,F01,x0,y0,z0,x1,y1,z1\n";
    summary_csv << "slice_index,axis,plane_coordinate,grid_samples,finite_samples,"
                   "nonfinite_samples,contour_segments,skipped_nonfinite_cells,"
                   "unresolved_degenerate_cells\n";

    const std::array<std::array<int, 2>, 4> edge_corners{{
        {{0, 1}}, {{1, 2}}, {{2, 3}}, {{3, 0}}
    }};
    const double point_tolerance = std::max(1e-24, result.spacing * result.spacing * 1e-20);

    for (int slice = 0; slice <= resolution; ++slice)
    {
        // Each cut coincides with a grid plane: -extent + slice * spacing.
        const double plane_coordinate = grid_coord(slice, resolution, half_extent);
        std::vector<double> values(static_cast<std::size_t>(side) * side,
                                   std::numeric_limits<double>::quiet_NaN());
        std::size_t finite_on_slice = 0;
        std::size_t invalid_on_slice = 0;
        std::size_t segments_on_slice = 0;
        std::size_t skipped_nonfinite_on_slice = 0;
        std::size_t unresolved_degenerate_on_slice = 0;

        auto sample_point = [&](int iu, int iv) {
            std::array<int, 3> idx{0, 0, 0};
            idx[axis] = slice; // exact grid-plane index along the selected axis
            idx[actual_u_axis] = iu;
            idx[actual_v_axis] = iv;
            std::array<double, 3> xyz{
                grid_coord(idx[0], resolution, half_extent),
                grid_coord(idx[1], resolution, half_extent),
                grid_coord(idx[2], resolution, half_extent)};
            // Do not offset the selected coordinate: it is the exact grid node.
            return Point3{xyz[0], xyz[1], xyz[2]};
        };

        for (int iv = 0; iv <= resolution; ++iv)
        {
            for (int iu = 0; iu <= resolution; ++iu)
            {
                const Point3 p = sample_point(iu, iv);
                const double f = evaluate(ast, p, k);
                values[static_cast<std::size_t>(iv) * side + iu] = f;
                std::array<int, 3> grid_index{0, 0, 0};
                grid_index[axis] = slice;
                grid_index[actual_u_axis] = iu;
                grid_index[actual_v_axis] = iv;
                result.samples.push_back(GridSample{
                    slice, iu, iv, grid_index[0], grid_index[1], grid_index[2], p, f});
                grid_csv << slice << ',' << axis_name(axis) << ',' << plane_coordinate
                         << ',' << iu << ',' << iv << ',' << grid_index[0] << ','
                         << grid_index[1] << ',' << grid_index[2] << ','
                         << p.x << ',' << p.y << ',' << p.z << ',' << f << ',';
                if (std::isfinite(f))
                {
                    grid_csv << (f < 0.0 ? 1 : 0) << '\n';
                    ++finite_on_slice;
                    ++result.finite_sample_count;
                }
                else
                {
                    grid_csv << "undefined\n";
                    ++invalid_on_slice;
                    ++result.invalid_sample_count;
                }
            }
        }

        auto value_at = [&](int iu, int iv) {
            return values[static_cast<std::size_t>(iv) * side + iu];
        };
        auto emit_segment = [&](int cell_u, int cell_v,
                                const std::array<double, 4>& f,
                                const Point3& a, const Point3& b) {
            if (squared_distance(a, b) <= point_tolerance)
            {
                ++unresolved_degenerate_on_slice;
                ++result.unresolved_degenerate_cells;
                return;
            }
            ContourSegment seg;
            seg.slice_index = slice;
            seg.cell_u = cell_u;
            seg.cell_v = cell_v;
            seg.plane_coordinate = plane_coordinate;
            seg.corner_values = f;
            seg.a = a;
            seg.b = b;
            result.segments.push_back(seg);
            segments_csv << slice << ',' << axis_name(axis) << ',' << plane_coordinate
                         << ',' << cell_u << ',' << cell_v << ','
                         << f[0] << ',' << f[1] << ',' << f[2] << ',' << f[3] << ',';
            write_point(segments_csv, a);
            segments_csv << ',';
            write_point(segments_csv, b);
            segments_csv << '\n';
            ++segments_on_slice;
            ++result.contour_segment_count;
        };

        for (int iv = 0; iv < resolution; ++iv)
        {
            for (int iu = 0; iu < resolution; ++iu)
            {
                const std::array<double, 4> f{{
                    value_at(iu, iv), value_at(iu + 1, iv),
                    value_at(iu + 1, iv + 1), value_at(iu, iv + 1)
                }};
                if (!std::all_of(f.begin(), f.end(), [](double v) { return std::isfinite(v); }))
                {
                    ++skipped_nonfinite_on_slice;
                    ++result.skipped_nonfinite_cells;
                    continue;
                }

                const std::array<Point3, 4> p{{
                    sample_point(iu, iv), sample_point(iu + 1, iv),
                    sample_point(iu + 1, iv + 1), sample_point(iu, iv + 1)
                }};
                std::array<Point3, 4> edge_points{};
                std::array<bool, 4> crossed{{false, false, false, false}};
                std::vector<Point3> unique_points;
                for (std::size_t edge = 0; edge < edge_corners.size(); ++edge)
                {
                    const int a = edge_corners[edge][0];
                    const int b = edge_corners[edge][1];
                    // Half-open sign convention avoids double-counting exact-zero vertices.
                    if ((f[a] < 0.0) == (f[b] < 0.0))
                        continue;
                    crossed[edge] = true;
                    edge_points[edge] = interpolate_edge(p[a], p[b], f[a], f[b]);
                    bool duplicate = false;
                    for (const Point3& existing : unique_points)
                    {
                        if (squared_distance(existing, edge_points[edge]) <= point_tolerance)
                        {
                            duplicate = true;
                            break;
                        }
                    }
                    if (!duplicate)
                        unique_points.push_back(edge_points[edge]);
                }

                if (unique_points.size() == 2)
                {
                    emit_segment(iu, iv, f, unique_points[0], unique_points[1]);
                    continue;
                }

                int crossing_count = 0;
                for (bool was_crossed : crossed)
                    crossing_count += was_crossed ? 1 : 0;
                if (crossing_count == 4 && unique_points.size() == 4)
                {
                    // Resolve the checkerboard case with the cell-center sign.
                    const Point3 center{
                        (p[0].x + p[2].x) * 0.5,
                        (p[0].y + p[2].y) * 0.5,
                        (p[0].z + p[2].z) * 0.5};
                    const double fc = evaluate(ast, center, k);
                    if (!std::isfinite(fc))
                    {
                        ++skipped_nonfinite_on_slice;
                        ++result.skipped_nonfinite_cells;
                        continue;
                    }
                    if ((fc < 0.0) == (f[0] < 0.0))
                    {
                        emit_segment(iu, iv, f, edge_points[0], edge_points[1]);
                        emit_segment(iu, iv, f, edge_points[2], edge_points[3]);
                    }
                    else
                    {
                        emit_segment(iu, iv, f, edge_points[3], edge_points[0]);
                        emit_segment(iu, iv, f, edge_points[1], edge_points[2]);
                    }
                }
                else if (!unique_points.empty())
                {
                    // A contour passing exactly through a grid corner can create
                    // one-point or three-point patterns. Do not fabricate a line
                    // through them; record them as unresolved so the CSV remains honest.
                    ++unresolved_degenerate_on_slice;
                    ++result.unresolved_degenerate_cells;
                }
            }
        }

        summary_csv << slice << ',' << axis_name(axis) << ',' << plane_coordinate
                    << ',' << static_cast<std::size_t>(side) * side << ','
                    << finite_on_slice << ',' << invalid_on_slice << ','
                    << segments_on_slice << ',' << skipped_nonfinite_on_slice << ','
                    << unresolved_degenerate_on_slice << '\n';
    }

    if (!grid_csv || !segments_csv || !summary_csv)
        throw std::runtime_error("A slice CSV could not be written completely.");

    metadata << "2D isosurface slice extraction\n"
             << "Expression: " << expression << '\n'
             << "Isovalue: F(x,y,z,k) = 0\n"
             << "k: " << std::setprecision(17) << k << '\n'
             << "Axis: " << axis_name(axis) << " (planes normal to this axis)\n"
             << "Domain: [-" << half_extent << ", " << half_extent << "] on each axis\n"
             << "Grid resolution: " << resolution << " cells per axis\n"
             << "Grid spacing: " << result.spacing << '\n'
             << "Parallel slice planes: " << result.plane_count
             << " (exactly on grid coordinates, including both domain boundaries)\n"
             << "Samples per slice: " << (resolution + 1) * (resolution + 1) << '\n'
             << "Total slice-plane sample rows: " << result.sample_count << '\n'
             << "Finite sample values: " << result.finite_sample_count << '\n'
             << "Non-finite sample values: " << result.invalid_sample_count << '\n'
             << "Contour segments: " << result.contour_segment_count << '\n'
             << "Skipped non-finite cells: " << result.skipped_nonfinite_cells << '\n'
             << "Unresolved exact-zero/degenerate cells: " << result.unresolved_degenerate_cells << '\n'
             << "Files: slice_grid_samples.csv, slice_contour_segments.csv, slice_summary.csv\n";
    if (!metadata)
        throw std::runtime_error("Could not write slice metadata.");
    return result;
}
} // namespace iso_slices
