#pragma once

// Snapshot of a rectangular window of grid nodes across several extracted 2D
// slices: an SVG figure (one panel per slice, every node labelled with its F
// value, the F=0 contour, and the 3D coordinates of the bottom-left and
// top-right corners) plus a CSV with the same nodes at full precision.
//
// Works only on an existing iso_slices::Result, so the snapshot shows exactly
// the data the slice gallery displays (same k, same expression).

#include "CsvFormat.h"
#include "FileNames.h"
#include "IsosurfaceSlices.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace iso_slices
{

inline constexpr std::size_t kMinSnapshotSlices = 3;

struct SnapshotRequest
{
    std::vector<int> slices; // 0-based slice indices (the gallery shows them 1-based)
    int iu_min = 0;          // inclusive grid-index range along the in-plane u axis
    int iu_max = 0;
    int iv_min = 0;          // inclusive grid-index range along the in-plane v axis
    int iv_max = 0;
};

struct SnapshotFiles
{
    std::filesystem::path svg;
    std::filesystem::path csv;
};

// Empty when `request` can be taken from `result`, otherwise the reason it cannot.
inline std::string validate_snapshot(const Result& result, const SnapshotRequest& request)
{
    if (result.samples.empty() || result.plane_count == 0)
        return "Extract slices first.";
    std::vector<int> slices = request.slices;
    std::sort(slices.begin(), slices.end());
    if (std::adjacent_find(slices.begin(), slices.end()) != slices.end())
        return "Each slice can be selected only once.";
    if (slices.size() < kMinSnapshotSlices)
        return "Select at least " + std::to_string(kMinSnapshotSlices) + " slices.";
    if (slices.front() < 0 || slices.back() >= static_cast<int>(result.plane_count))
        return "A selected slice does not exist in the current result.";
    const int n = result.resolution;
    if (request.iu_min < 0 || request.iv_min < 0 || request.iu_max > n || request.iv_max > n)
        return "Grid indices must lie between 0 and " + std::to_string(n) + ".";
    if (request.iu_min >= request.iu_max || request.iv_min >= request.iv_max)
        return "The grid window must be a rectangle: min index < max index on both axes.";
    return {};
}

// 1-based slice numbers with consecutive runs collapsed: {1,2,3,6} -> "2-4_7".
inline std::string slice_list_label(std::vector<int> slices)
{
    std::sort(slices.begin(), slices.end());
    std::string label;
    for (std::size_t i = 0; i < slices.size();)
    {
        std::size_t j = i;
        while (j + 1 < slices.size() && slices[j + 1] == slices[j] + 1) ++j;
        if (!label.empty()) label += '_';
        label += std::to_string(slices[i] + 1);
        if (j > i) label += '-' + std::to_string(slices[j] + 1);
        i = j + 1;
    }
    return label;
}

// e.g. "x^3-kxy_slices_2-4_u3-5_v4-6" (expression, slices, grid window).
inline std::string snapshot_base_name(const Result& result, const SnapshotRequest& request)
{
    return iso_files::safeFileComponent(result.expression.empty() ? "F" : result.expression, 48) +
           "_slices_" + slice_list_label(request.slices) +
           "_u" + std::to_string(request.iu_min) + '-' + std::to_string(request.iu_max) +
           "_v" + std::to_string(request.iv_min) + '-' + std::to_string(request.iv_max);
}

namespace detail
{
inline const GridSample& sample_at(const Result& r, int slice, int iu, int iv)
{
    const std::size_t side = static_cast<std::size_t>(r.resolution) + 1;
    return r.samples[static_cast<std::size_t>(slice) * side * side + static_cast<std::size_t>(iv) * side +
                     static_cast<std::size_t>(iu)];
}

inline std::string short_number(double value)
{
    if (!std::isfinite(value)) return "undefined";
    char buffer[32];
    std::snprintf(buffer, sizeof(buffer), "%.6g", value);
    return buffer;
}

inline std::string xml_escape(const std::string& text)
{
    std::string out;
    for (char c : text)
    {
        switch (c)
        {
        case '&': out += "&amp;"; break;
        case '<': out += "&lt;"; break;
        case '>': out += "&gt;"; break;
        case '"': out += "&quot;"; break;
        default: out += c;
        }
    }
    return out;
}

inline std::string node_label(const GridSample& s)
{
    return "[" + std::to_string(s.ix) + "," + std::to_string(s.iy) + "," + std::to_string(s.iz) + "] (" +
           short_number(s.point.x) + ", " + short_number(s.point.y) + ", " + short_number(s.point.z) + ")";
}

// Fill colour by sign, matching the slice gallery (blue < 0, red > 0).
inline const char* value_color(double f)
{
    if (!std::isfinite(f)) return "#d83a3a";
    if (std::fabs(f) <= 1e-8) return "#c99a06";
    return f < 0.0 ? "#2f6fd6" : "#d4512e";
}

// First unused "<directory>/<base><suffix>.svg" / ".csv" pair.
inline SnapshotFiles unique_files(const std::filesystem::path& directory, const std::string& base)
{
    for (int n = 1;; ++n)
    {
        std::string name = base;
        if (n > 1)
        {
            char suffix[16];
            std::snprintf(suffix, sizeof(suffix), "_%02d", n);
            name += suffix;
        }
        SnapshotFiles files{directory / (name + ".svg"), directory / (name + ".csv")};
        if (!std::filesystem::exists(files.svg) && !std::filesystem::exists(files.csv)) return files;
    }
}
} // namespace detail

inline void write_snapshot_csv(std::ostream& out, const Result& r, const SnapshotRequest& request,
                               const std::vector<int>& slices)
{
    out << "slice_number,slice_index,axis,plane_coordinate,iu,iv,grid_ix,grid_iy,grid_iz,x,y,z,F,"
           "expression,k\n";
    const std::string expression = iso_csv::quote(r.expression);
    for (int slice : slices)
        for (int iv = request.iv_min; iv <= request.iv_max; ++iv)
            for (int iu = request.iu_min; iu <= request.iu_max; ++iu)
            {
                const GridSample& s = detail::sample_at(r, slice, iu, iv);
                out << slice + 1 << ',' << slice << ',' << axis_name(r.axis) << ','
                    << iso_csv::number(gridCoordinate(slice, r.resolution, r.half_extent)) << ','
                    << iu << ',' << iv << ',' << s.ix << ',' << s.iy << ',' << s.iz << ','
                    << iso_csv::number(s.point.x) << ',' << iso_csv::number(s.point.y) << ','
                    << iso_csv::number(s.point.z) << ',' << iso_csv::number(s.value) << ','
                    << expression << ',' << iso_csv::number(r.k) << '\n';
            }
}

inline void write_snapshot_svg(std::ostream& out, const Result& r, const SnapshotRequest& request,
                               const std::vector<int>& slices)
{
    using detail::short_number;
    using detail::xml_escape;
    // Layout (px): panels in rows of up to three; each node gets a cell.
    constexpr double cell_w = 96.0, cell_h = 58.0, pad_x = 52.0, pad_top = 64.0, pad_bottom = 52.0;
    constexpr double margin = 24.0, header_h = 92.0, gap = 28.0;
    const int cols = request.iu_max - request.iu_min;
    const int rows = request.iv_max - request.iv_min;
    const double plot_w = cols * cell_w, plot_h = rows * cell_h;
    const double panel_w = std::max(plot_w + 2.0 * pad_x, 380.0);
    const double panel_h = plot_h + pad_top + pad_bottom;
    const int per_row = static_cast<int>(std::min<std::size_t>(slices.size(), 3));
    const int panel_rows = static_cast<int>((slices.size() + per_row - 1) / per_row);
    const double width = 2.0 * margin + per_row * panel_w + (per_row - 1) * gap;
    const double height = 2.0 * margin + header_h + panel_rows * panel_h + (panel_rows - 1) * gap;
    const std::array<int, 2> uv = in_plane_axes(r.axis);
    const char* u_name = axis_name(uv[0]);
    const char* v_name = axis_name(uv[1]);

    out << std::setprecision(6) << std::fixed;
    out << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"" << width << "\" height=\"" << height
        << "\" viewBox=\"0 0 " << width << ' ' << height << "\" font-family=\"Helvetica, Arial, sans-serif\">\n"
        << "<rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n";
    auto text = [&](double x, double y, const std::string& s, int size, const char* fill = "#1f2933",
                    const char* anchor = "start", const char* weight = "normal") {
        out << "<text x=\"" << x << "\" y=\"" << y << "\" font-size=\"" << size << "\" fill=\"" << fill
            << "\" text-anchor=\"" << anchor << "\" font-weight=\"" << weight << "\">" << xml_escape(s)
            << "</text>\n";
    };

    double y = margin + 20.0;
    text(margin, y, "Slice snapshot: F(x,y,z,k) = " + r.expression, 18, "#1f2933", "start", "bold");
    text(margin, y += 24.0, "k = " + short_number(r.k) + "   |   " + std::string(axis_name(r.axis)) +
         "-normal slices   |   grid resolution " + std::to_string(r.resolution) + ", domain [-" +
         short_number(r.half_extent) + ", " + short_number(r.half_extent) + "]^3", 13);
    text(margin, y += 20.0, "Window: iu " + std::to_string(request.iu_min) + ".." + std::to_string(request.iu_max) +
         " (+u = " + u_name + "), iv " + std::to_string(request.iv_min) + ".." + std::to_string(request.iv_max) +
         " (+v = " + v_name + ").  Blue: F < 0, red: F > 0, teal: F = 0 contour.", 13, "#52606d");

    for (std::size_t p = 0; p < slices.size(); ++p)
    {
        const int slice = slices[p];
        const double px = margin + static_cast<double>(p % per_row) * (panel_w + gap);
        const double py = margin + header_h + static_cast<double>(p / per_row) * (panel_h + gap);
        const double left = px + (panel_w - plot_w) / 2.0, top = py + pad_top;
        auto node_x = [&](int iu) { return left + (iu - request.iu_min) * cell_w; };
        auto node_y = [&](int iv) { return top + plot_h - (iv - request.iv_min) * cell_h; };

        out << "<rect x=\"" << px << "\" y=\"" << py << "\" width=\"" << panel_w << "\" height=\"" << panel_h
            << "\" rx=\"6\" fill=\"#f7f9fc\" stroke=\"#c9d2dc\"/>\n";
        text(px + 12.0, py + 22.0, "Slice " + std::to_string(slice + 1) + "  |  " + axis_name(r.axis) + " = " +
             short_number(gridCoordinate(slice, r.resolution, r.half_extent)), 14, "#1f2933", "start", "bold");

        // Grid lines through the window's nodes.
        out << "<g stroke=\"#d5dde6\" stroke-width=\"1\">\n";
        for (int iu = request.iu_min; iu <= request.iu_max; ++iu)
            out << "<line x1=\"" << node_x(iu) << "\" y1=\"" << top << "\" x2=\"" << node_x(iu)
                << "\" y2=\"" << top + plot_h << "\"/>\n";
        for (int iv = request.iv_min; iv <= request.iv_max; ++iv)
            out << "<line x1=\"" << left << "\" y1=\"" << node_y(iv) << "\" x2=\"" << left + plot_w
                << "\" y2=\"" << node_y(iv) << "\"/>\n";
        out << "</g>\n";

        // F = 0 contour segments of the cells inside the window.
        const double u0 = gridCoordinate(request.iu_min, r.resolution, r.half_extent);
        const double v0 = gridCoordinate(request.iv_min, r.resolution, r.half_extent);
        const double du = r.spacing > 0.0 ? cell_w / r.spacing : 0.0;
        const double dv = r.spacing > 0.0 ? cell_h / r.spacing : 0.0;
        auto map_x = [&](const Point3& q) { return left + (iso_geometry::component(q, uv[0]) - u0) * du; };
        auto map_y = [&](const Point3& q) { return top + plot_h - (iso_geometry::component(q, uv[1]) - v0) * dv; };
        out << "<g stroke=\"#0f9e92\" stroke-width=\"2.5\" stroke-linecap=\"round\">\n";
        for (const ContourSegment& s : r.segments)
        {
            if (s.slice_index != slice || s.cell_u < request.iu_min || s.cell_u >= request.iu_max ||
                s.cell_v < request.iv_min || s.cell_v >= request.iv_max)
                continue;
            out << "<line x1=\"" << map_x(s.a) << "\" y1=\"" << map_y(s.a) << "\" x2=\"" << map_x(s.b)
                << "\" y2=\"" << map_y(s.b) << "\"/>\n";
        }
        out << "</g>\n";

        // Nodes with their isovalues.
        for (int iv = request.iv_min; iv <= request.iv_max; ++iv)
            for (int iu = request.iu_min; iu <= request.iu_max; ++iu)
            {
                const GridSample& s = detail::sample_at(r, slice, iu, iv);
                const char* color = detail::value_color(s.value);
                out << "<circle cx=\"" << node_x(iu) << "\" cy=\"" << node_y(iv) << "\" r=\"4\" fill=\""
                    << color << "\"/>\n";
                text(node_x(iu), node_y(iv) - 8.0, short_number(s.value), 11, color, "middle");
            }

        // Corner coordinates: 3D grid index and point.
        const GridSample& bottom_left = detail::sample_at(r, slice, request.iu_min, request.iv_min);
        const GridSample& top_right = detail::sample_at(r, slice, request.iu_max, request.iv_max);
        text(left + plot_w, top - 26.0, "top-right " + detail::node_label(top_right), 11, "#52606d", "end");
        text(left, top + plot_h + 22.0, "bottom-left " + detail::node_label(bottom_left), 11, "#52606d");
        text(left + plot_w, top + plot_h + 40.0, "+u = " + std::string(u_name) + " →, +v = " + v_name +
             " ↑", 10, "#9aa5b1", "end");
    }
    out << "</svg>\n";
}

// The selected slices (1-based, as the gallery shows them) on which no part of
// the F = 0 contour lies inside the window. A snapshot of such a slice shows
// only node values and no surface; the viewer warns about it before saving.
inline std::vector<int> slices_without_contour(const Result& result, const SnapshotRequest& request)
{
    std::vector<int> empty;
    for (int slice : request.slices)
    {
        const bool found = std::any_of(result.segments.begin(), result.segments.end(), [&](const ContourSegment& s) {
            return s.slice_index == slice && s.cell_u >= request.iu_min && s.cell_u < request.iu_max &&
                   s.cell_v >= request.iv_min && s.cell_v < request.iv_max;
        });
        if (!found) empty.push_back(slice + 1);
    }
    std::sort(empty.begin(), empty.end());
    return empty;
}

// Writes "<base>.svg" and "<base>.csv" into `directory` (created if needed;
// an existing snapshot of the same name gets a _02, _03, ... suffix).
inline SnapshotFiles write_snapshot(const Result& result, const SnapshotRequest& request,
                                    const std::filesystem::path& directory)
{
    const std::string error = validate_snapshot(result, request);
    if (!error.empty())
        throw std::invalid_argument(error);
    std::vector<int> slices = request.slices;
    std::sort(slices.begin(), slices.end());

    std::filesystem::create_directories(directory);
    const SnapshotFiles files = detail::unique_files(directory, snapshot_base_name(result, request));
    std::ofstream svg(files.svg);
    std::ofstream csv(files.csv);
    if (!svg || !csv)
        throw std::runtime_error("Could not create snapshot files in " + directory.string());
    write_snapshot_svg(svg, result, request, slices);
    write_snapshot_csv(csv, result, request, slices);
    if (!svg || !csv)
        throw std::runtime_error("Could not write snapshot files in " + directory.string());
    return files;
}

} // namespace iso_slices
