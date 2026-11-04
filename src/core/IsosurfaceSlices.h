#pragma once

// Grid-aligned 2D cross-sections of the current implicit field F(x,y,z,k)=0.
// Each selected plane is sampled on the same lattice as the 3D extractor, and
// marching squares extracts contour line segments from those samples.
#include "CsvFormat.h"
#include "Expr.h"
#include "ExprProgram.h"
#include "Geometry.h"
#include "Grid.h"
#include "Parallel.h"

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
using Point3 = iso_geometry::Point3;

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

// The two in-plane axes for slices normal to `axis` (u first, then v).
inline std::array<int, 2> in_plane_axes(int axis)
{
    return {axis == 0 ? 1 : 0, axis == 2 ? 1 : 2};
}

// Linear zero crossing between a and b (t clamped to the segment).
inline Point3 interpolate_edge(const Point3& a, const Point3& b, double fa, double fb)
{
    const double denom = fa - fb;
    double t = std::abs(denom) > std::numeric_limits<double>::min() ? fa / denom : 0.5;
    t = std::clamp(t, 0.0, 1.0);
    return {a.x + t * (b.x - a.x), a.y + t * (b.y - a.y), a.z + t * (b.z - a.z)};
}

namespace detail
{
// Buffers CSV text and writes it to the stream in large chunks.
class BufferedCsv
{
public:
    explicit BufferedCsv(std::ofstream& out) : out_(out) { line_.reserve(kFlushBytes + 4096); }
    ~BufferedCsv() { flush(); }
    BufferedCsv(const BufferedCsv&) = delete;
    BufferedCsv& operator=(const BufferedCsv&) = delete;

    template <typename T>
    BufferedCsv& operator<<(const T& value)
    {
        iso_csv::append(line_, value);
        return *this;
    }
    BufferedCsv& operator<<(char c)
    {
        line_ += c;
        if (c == '\n' && line_.size() >= kFlushBytes) flush();
        return *this;
    }
    BufferedCsv& operator<<(const Point3& p)
    {
        return *this << p.x << ',' << p.y << ',' << p.z;
    }
    void flush()
    {
        out_.write(line_.data(), static_cast<std::streamsize>(line_.size()));
        line_.clear();
    }

private:
    static constexpr std::size_t kFlushBytes = 1 << 20;
    std::ofstream& out_;
    std::string line_;
};
} // namespace detail

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
    std::ofstream grid_file(out_dir / "slice_grid_samples.csv", std::ios::binary);
    std::ofstream segments_file(out_dir / "slice_contour_segments.csv", std::ios::binary);
    std::ofstream summary_file(out_dir / "slice_summary.csv", std::ios::binary);
    std::ofstream metadata(out_dir / "slice_metadata.txt");
    if (!grid_file || !segments_file || !summary_file || !metadata)
        throw std::runtime_error("Could not create one or more slice output files in: " + output_directory);

    Result result;
    result.axis = axis;
    result.resolution = resolution;
    result.k = k;
    result.half_extent = half_extent;
    result.spacing = (2.0 * half_extent) / resolution;
    result.plane_count = static_cast<std::size_t>(resolution + 1);
    const int side = resolution + 1;
    const std::size_t per_slice = static_cast<std::size_t>(side) * side;
    result.sample_count = result.plane_count * per_slice;

    // The slice normal is the selected axis. The other two axes are sampled
    // at their original grid vertices, so each slice is a true lattice plane.
    const std::array<int, 2> uv = in_plane_axes(axis);
    std::vector<double> coords(static_cast<std::size_t>(side));
    for (int i = 0; i < side; ++i) coords[static_cast<std::size_t>(i)] = gridCoordinate(i, resolution, half_extent);

    auto lattice_index = [&](int slice, int iu, int iv) {
        std::array<int, 3> idx{0, 0, 0};
        idx[axis] = slice; // exact grid-plane index along the selected axis
        idx[uv[0]] = iu;
        idx[uv[1]] = iv;
        return idx;
    };
    auto lattice_point = [&](const std::array<int, 3>& idx) {
        return Point3{coords[static_cast<std::size_t>(idx[0])], coords[static_cast<std::size_t>(idx[1])],
                      coords[static_cast<std::size_t>(idx[2])]};
    };

    // Evaluate every slice in parallel; samples are stored in their final order.
    // Each row of constant iv is one batch for the compiled program (only the
    // u coordinate varies along it).
    const ExprProgram program(ast);
    result.samples.resize(result.sample_count);
    iso_parallel::parallelFor(result.plane_count, [&](std::size_t slice_index) {
        const int slice = static_cast<int>(slice_index);
        GridSample* out = result.samples.data() + slice_index * per_slice;
        std::vector<double> row(static_cast<std::size_t>(side));
        std::vector<double> scratch;
        for (int iv = 0; iv <= resolution; ++iv)
        {
            ExprProgram::Lane lane[3];
            const std::array<int, 3> first = lattice_index(slice, 0, iv);
            for (int a = 0; a < 3; ++a)
                lane[a] = a == uv[0] ? ExprProgram::Lane::array(coords.data())
                                     : ExprProgram::Lane::shared(coords[static_cast<std::size_t>(first[a])]);
            program.evaluate(lane[0], lane[1], lane[2], k, row.size(), row.data(), scratch);
            for (int iu = 0; iu <= resolution; ++iu)
            {
                const std::array<int, 3> idx = lattice_index(slice, iu, iv);
                *out++ = GridSample{slice, iu, iv, idx[0], idx[1], idx[2], lattice_point(idx),
                                    row[static_cast<std::size_t>(iu)]};
            }
        }
    });

    detail::BufferedCsv grid_csv(grid_file);
    detail::BufferedCsv segments_csv(segments_file);
    detail::BufferedCsv summary_csv(summary_file);
    grid_csv << "slice_index,axis,plane_coordinate,iu,iv,grid_ix,grid_iy,grid_iz,x,y,z,F,is_negative\n";
    segments_csv << "slice_index,axis,plane_coordinate,cell_u,cell_v,"
                    "F00,F10,F11,F01,x0,y0,z0,x1,y1,z1\n";
    summary_csv << "slice_index,axis,plane_coordinate,grid_samples,finite_samples,"
                   "nonfinite_samples,contour_segments,skipped_nonfinite_cells,"
                   "unresolved_degenerate_cells\n";

    constexpr int edge_corners[4][2] = {{0, 1}, {1, 2}, {2, 3}, {3, 0}};
    const double point_tolerance = std::max(1e-24, result.spacing * result.spacing * 1e-20);
    const char* axis_label = axis_name(axis);

    for (int slice = 0; slice <= resolution; ++slice)
    {
        // Each cut coincides with a grid plane: -extent + slice * spacing.
        const double plane_coordinate = coords[static_cast<std::size_t>(slice)];
        const GridSample* slice_samples = result.samples.data() + static_cast<std::size_t>(slice) * per_slice;
        std::size_t finite_on_slice = 0;
        std::size_t invalid_on_slice = 0;
        std::size_t segments_on_slice = 0;
        std::size_t skipped_nonfinite_on_slice = 0;
        std::size_t unresolved_degenerate_on_slice = 0;

        for (std::size_t s = 0; s < per_slice; ++s)
        {
            const GridSample& g = slice_samples[s];
            grid_csv << slice << ',' << axis_label << ',' << plane_coordinate << ',' << g.iu << ','
                     << g.iv << ',' << g.ix << ',' << g.iy << ',' << g.iz << ',' << g.point << ','
                     << g.value << ',';
            if (std::isfinite(g.value))
            {
                grid_csv << (g.value < 0.0 ? "1" : "0") << '\n';
                ++finite_on_slice;
            }
            else
            {
                grid_csv << "undefined" << '\n';
                ++invalid_on_slice;
            }
        }
        result.finite_sample_count += finite_on_slice;
        result.invalid_sample_count += invalid_on_slice;

        auto sample_at = [&](int iu, int iv) -> const GridSample& {
            return slice_samples[static_cast<std::size_t>(iv) * side + static_cast<std::size_t>(iu)];
        };
        auto emit_segment = [&](int cell_u, int cell_v, const std::array<double, 4>& f,
                                const Point3& a, const Point3& b) {
            if (iso_geometry::normSquared(a - b) <= point_tolerance)
            {
                ++unresolved_degenerate_on_slice;
                return;
            }
            result.segments.push_back(ContourSegment{slice, cell_u, cell_v, plane_coordinate, f, a, b});
            segments_csv << slice << ',' << axis_label << ',' << plane_coordinate << ',' << cell_u << ','
                         << cell_v << ',' << f[0] << ',' << f[1] << ',' << f[2] << ',' << f[3] << ','
                         << a << ',' << b << '\n';
            ++segments_on_slice;
        };

        for (int iv = 0; iv < resolution; ++iv)
        {
            for (int iu = 0; iu < resolution; ++iu)
            {
                const GridSample* corner[4] = {&sample_at(iu, iv), &sample_at(iu + 1, iv),
                                               &sample_at(iu + 1, iv + 1), &sample_at(iu, iv + 1)};
                const std::array<double, 4> f{{corner[0]->value, corner[1]->value,
                                               corner[2]->value, corner[3]->value}};
                if (!std::all_of(f.begin(), f.end(), [](double v) { return std::isfinite(v); }))
                {
                    ++skipped_nonfinite_on_slice;
                    continue;
                }

                std::array<Point3, 4> edge_points{};
                std::array<Point3, 4> unique_points{};
                int crossing_count = 0;
                std::size_t unique_count = 0;
                for (int edge = 0; edge < 4; ++edge)
                {
                    const int a = edge_corners[edge][0];
                    const int b = edge_corners[edge][1];
                    // Half-open sign convention avoids double-counting exact-zero vertices.
                    if ((f[a] < 0.0) == (f[b] < 0.0))
                        continue;
                    ++crossing_count;
                    edge_points[edge] = interpolate_edge(corner[a]->point, corner[b]->point, f[a], f[b]);
                    const bool duplicate = std::any_of(
                        unique_points.begin(), unique_points.begin() + static_cast<std::ptrdiff_t>(unique_count),
                        [&](const Point3& existing) {
                            return iso_geometry::normSquared(existing - edge_points[edge]) <= point_tolerance;
                        });
                    if (!duplicate)
                        unique_points[unique_count++] = edge_points[edge];
                }

                if (unique_count == 2)
                {
                    emit_segment(iu, iv, f, unique_points[0], unique_points[1]);
                }
                else if (crossing_count == 4 && unique_count == 4)
                {
                    // Resolve the checkerboard case with the cell-center sign.
                    const Point3 center = (corner[0]->point + corner[2]->point) * 0.5;
                    const double fc = evalExpr<double>(ast, center.x, center.y, center.z, k);
                    if (!std::isfinite(fc))
                    {
                        ++skipped_nonfinite_on_slice;
                    }
                    else if ((fc < 0.0) == (f[0] < 0.0))
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
                else if (unique_count != 0)
                {
                    // A contour passing exactly through a grid corner can create
                    // one-point or three-point patterns. Do not fabricate a line
                    // through them; record them as unresolved so the CSV remains honest.
                    ++unresolved_degenerate_on_slice;
                }
            }
        }

        result.contour_segment_count += segments_on_slice;
        result.skipped_nonfinite_cells += skipped_nonfinite_on_slice;
        result.unresolved_degenerate_cells += unresolved_degenerate_on_slice;
        summary_csv << slice << ',' << axis_label << ',' << plane_coordinate << ',' << per_slice << ','
                    << finite_on_slice << ',' << invalid_on_slice << ',' << segments_on_slice << ','
                    << skipped_nonfinite_on_slice << ',' << unresolved_degenerate_on_slice << '\n';
    }

    grid_csv.flush();
    segments_csv.flush();
    summary_csv.flush();
    if (!grid_file || !segments_file || !summary_file)
        throw std::runtime_error("A slice CSV could not be written completely.");

    metadata << "2D isosurface slice extraction\n"
             << "Expression: " << expression << '\n'
             << "Isovalue: F(x,y,z,k) = 0\n"
             << "k: " << std::setprecision(17) << k << '\n'
             << "Axis: " << axis_label << " (planes normal to this axis)\n"
             << "Domain: [-" << half_extent << ", " << half_extent << "] on each axis\n"
             << "Grid resolution: " << resolution << " cells per axis\n"
             << "Grid spacing: " << result.spacing << '\n'
             << "Parallel slice planes: " << result.plane_count
             << " (exactly on grid coordinates, including both domain boundaries)\n"
             << "Samples per slice: " << per_slice << '\n'
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
