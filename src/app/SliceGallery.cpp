// 2D slice extraction: state handling, the 3D slice overlays and the
// side-by-side slice gallery window.

#include "IsosurfaceViewer.h"

#include "ColorMap.h"
#include "core/Grid.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <utility>

namespace iso_app {

void IsosurfaceViewer::invalidate_slice_overlay()
{
    // Parameter changes retain the last extracted slice snapshot so it can be compared
    // or re-extracted. The gallery clearly marks it stale until refreshed.
    if (slice_result_available_) slice_data_stale_ = true;
}

void IsosurfaceViewer::clear_slice_overlay()
{
    slice_result_data_ = iso_slices::Result{};
    slice_segment_ranges_.clear();
    slice_selected_grid_point_valid_ = false;
    slice_result_available_ = false;
    slice_data_stale_ = false;
    slice_reference_geometry_dirty_ = true;
    overlay_.set(OverlayRenderer::SliceContours, {});
}

iso_slices::Result IsosurfaceViewer::extract_slices_at_current_parameters()
{
    sanitize_parameters();
    return iso_slices::extract_and_export(
        ast_.get(), last_good_expression_, static_cast<double>(k_),
        grid_resolution_, static_cast<double>(half_extent_), slice_axis_,
        slice_output_directory_.data());
}

// Installs a fresh slice result: 3D contour overlay (colored by slice
// index), gallery state and counters. Selection indices are clamped.
void IsosurfaceViewer::apply_slice_result(iso_slices::Result&& result)
{
    slice_result_data_ = std::move(result);
    slice_data_stale_ = false;
    std::vector<OverlayVertex> contours;
    contours.reserve(slice_result_data_.segments.size() * 2);
    const double denom = static_cast<double>(std::max(1, slice_result_data_.resolution));
    for (const iso_slices::ContourSegment& segment : slice_result_data_.segments)
    {
        const pmp::Color color = gradient_color(static_cast<float>(segment.slice_index / denom));
        contours.push_back(overlay_vertex(segment.a.x, segment.a.y, segment.a.z, color));
        contours.push_back(overlay_vertex(segment.b.x, segment.b.y, segment.b.z, color));
    }
    overlay_.set(OverlayRenderer::SliceContours, std::move(contours));
    slice_result_available_ = true;
    show_slice_contours_ = true;
    const int resolution = std::max(1, slice_result_data_.resolution);
    slice_selected_slice_ = std::clamp(slice_selected_slice_, 0, resolution);
    slice_selected_iu_ = std::clamp(slice_selected_iu_, 0, resolution);
    slice_selected_iv_ = std::clamp(slice_selected_iv_, 0, resolution);
    slice_selected_grid_point_valid_ = !slice_result_data_.samples.empty();
    rebuild_slice_segment_ranges();
    slice_reference_geometry_dirty_ = true;
    show_slice_gallery_ = true;
}

void IsosurfaceViewer::rebuild_slice_segment_ranges()
{
    slice_segment_ranges_.assign(slice_result_data_.plane_count, {0u, 0u});
    std::size_t cursor = 0;
    for (std::size_t slice = 0; slice < slice_result_data_.plane_count; ++slice)
    {
        const std::size_t begin = cursor;
        while (cursor < slice_result_data_.segments.size() &&
               slice_result_data_.segments[cursor].slice_index == static_cast<int>(slice))
            ++cursor;
        slice_segment_ranges_[slice] = {begin, cursor};
    }
}

// Draw the grid plane currently represented by the selected 2D panel and a
// 3D cross at the exact grid node selected in the gallery.
void IsosurfaceViewer::build_slice_reference_geometry(std::vector<OverlayVertex>& plane,
                                                      std::vector<OverlayVertex>& selected_point) const
{
    plane.clear();
    selected_point.clear();
    if (!slice_result_available_ || slice_result_data_.plane_count == 0 ||
        slice_result_data_.resolution < 1)
        return;

    const int resolution = slice_result_data_.resolution;
    const int side = resolution + 1;
    const int slice = std::clamp(slice_selected_slice_, 0, resolution);
    const double extent = slice_result_data_.half_extent;
    const int axis = slice_result_data_.axis;
    const int u_axis = (axis == 0) ? 1 : 0;
    const int v_axis = (axis == 2) ? 1 : 2;

    auto point_at = [&](int iu, int iv) -> iso_slices::Point3 {
        std::array<int, 3> idx{0, 0, 0};
        idx[axis] = slice;
        idx[u_axis] = iu;
        idx[v_axis] = iv;
        return {gridCoordinate(idx[0], resolution, extent),
                gridCoordinate(idx[1], resolution, extent),
                gridCoordinate(idx[2], resolution, extent)};
    };
    auto add_line = [&](const iso_slices::Point3& a, const iso_slices::Point3& b) {
        const float r = 0.22f, g = 0.55f, bl = 0.95f;
        plane.push_back(overlay_vertex(a.x, a.y, a.z, r, g, bl));
        plane.push_back(overlay_vertex(b.x, b.y, b.z, r, g, bl));
    };
    // Keep the 3D guide readable at high resolutions while always including
    // the box edges. The slice data and 2D panel still retain every node.
    const int stride = std::max(1, resolution / 24);
    for (int i = 0; i <= resolution; i += stride)
    {
        add_line(point_at(i, 0), point_at(i, resolution));
        add_line(point_at(0, i), point_at(resolution, i));
    }
    if (resolution % stride != 0)
    {
        add_line(point_at(resolution, 0), point_at(resolution, resolution));
        add_line(point_at(0, resolution), point_at(resolution, resolution));
    }

    if (slice_selected_grid_point_valid_)
    {
        const int iu = std::clamp(slice_selected_iu_, 0, resolution);
        const int iv = std::clamp(slice_selected_iv_, 0, resolution);
        const std::size_t sample_index = static_cast<std::size_t>(slice) * side * side +
            static_cast<std::size_t>(iv) * side + static_cast<std::size_t>(iu);
        if (sample_index < slice_result_data_.samples.size())
        {
            const iso_slices::Point3 p = slice_result_data_.samples[sample_index].point;
            const double radius = std::max(slice_result_data_.spacing * 0.22,
                                           extent * 0.012);
            const float r = 1.0f, g = 0.24f, bl = 0.08f;
            selected_point.push_back(overlay_vertex(p.x-radius, p.y, p.z, r,g,bl));
            selected_point.push_back(overlay_vertex(p.x+radius, p.y, p.z, r,g,bl));
            selected_point.push_back(overlay_vertex(p.x, p.y-radius, p.z, r,g,bl));
            selected_point.push_back(overlay_vertex(p.x, p.y+radius, p.z, r,g,bl));
            selected_point.push_back(overlay_vertex(p.x, p.y, p.z-radius, r,g,bl));
            selected_point.push_back(overlay_vertex(p.x, p.y, p.z+radius, r,g,bl));
        }
    }
}

void IsosurfaceViewer::draw_slice_gallery()
{
    if (!show_slice_gallery_ || !slice_result_available_ ||
        slice_result_data_.plane_count == 0 || slice_result_data_.samples.empty())
        return;

    const ImVec2 display = ImGui::GetIO().DisplaySize;
    const float gallery_x = std::min(395.0f, std::max(0.0f, display.x - 500.0f));
    const float window_width = std::max(360.0f, display.x - gallery_x - 12.0f);
    const float window_height = std::max(420.0f, display.y - 54.0f);
    ImGui::SetNextWindowPos(ImVec2(gallery_x, 30.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(window_width, window_height), ImGuiCond_FirstUseEver);
    if (!ImGui::Begin("2D Slice Gallery##readable_k_refresh", &show_slice_gallery_, ImGuiWindowFlags_NoCollapse))
    {
        ImGui::End();
        return;
    }

    const int resolution = std::max(1, slice_result_data_.resolution);
    const int side = resolution + 1;
    ImGui::Text("%zu grid-aligned %s-normal planes | k = %.17g | resolution = %d",
                slice_result_data_.plane_count,
                iso_slices::axis_name(slice_result_data_.axis),
                slice_result_data_.k, resolution);
    ImGui::TextDisabled("Cuts are exactly on grid coordinates -L + i*h (including both domain ends), not at cell centers.");
    ImGui::Separator();
    if (draw_k_text_editor("Viewer k", gallery_k_error_))
    {
        sanitize_parameters();
        mesh_dirty_ = true;
        sweep_dirty_ = true;
        invalidate_slice_overlay();
        if (slice_result_available_) slice_data_stale_ = true;
    }
    ImGui::SameLine();
    if (ImGui::Button("Re-extract slices at current k"))
    {
        try
        {
            apply_slice_result(extract_slices_at_current_parameters());
            runtime_message_ = "2D slices refreshed at current k and CSV files updated in '" +
                std::string(slice_output_directory_.data()) + "'.";
        }
        catch (const std::exception& error)
        {
            runtime_message_ = std::string("2D slice refresh failed: ") + error.what();
            if (slice_result_available_) slice_data_stale_ = true;
        }
    }
    ImGui::SameLine();
    if (ImGui::Button("Create snapshot..."))
        open_snapshot_dialog();
    if (slice_data_stale_)
        ImGui::TextColored(ImVec4(1.0f, 0.65f, 0.18f, 1.0f),
            "STALE SNAPSHOT: slices are at k = %.17g; viewer k = %.9g. Re-extract to update contours and CSV.",
            slice_result_data_.k, static_cast<double>(k_));
    else
        ImGui::TextColored(ImVec4(0.25f, 0.85f, 0.45f, 1.0f),
            "Slices are current at k = %.17g", slice_result_data_.k);

    ImGui::Checkbox("Show grid-node isovalues", &show_slice_grid_values_);
    ImGui::SameLine();
    ImGui::Checkbox("Full precision labels (may overlap)", &show_slice_full_precision_labels_);
    if (slice_label_resolution_ != resolution)
    {
        slice_label_resolution_ = resolution;
        slice_value_label_stride_ = std::max(1, (resolution + 6) / 7);
    }
    ImGui::SetNextItemWidth(210.0f);
    ImGui::SliderInt("Label stride (1 = every node)", &slice_value_label_stride_, 1, resolution + 1);
    ImGui::SameLine();
    ImGui::TextDisabled("Hover/click for the exact 17-digit value.");
    ImGui::SameLine();
    ImGui::Checkbox("Show selected slice grid in 3D", &show_slice_grid_in_3d_);
    ImGui::SameLine();
    ImGui::Checkbox("Show selected point in 3D", &show_slice_selected_point_in_3d_);

    if (slice_selected_grid_point_valid_)
    {
        const int selected_slice = std::clamp(slice_selected_slice_, 0, resolution);
        const int selected_u = std::clamp(slice_selected_iu_, 0, resolution);
        const int selected_v = std::clamp(slice_selected_iv_, 0, resolution);
        const std::size_t idx = static_cast<std::size_t>(selected_slice) * side * side +
            static_cast<std::size_t>(selected_v) * side + static_cast<std::size_t>(selected_u);
        if (idx < slice_result_data_.samples.size())
        {
            const auto& sample = slice_result_data_.samples[idx];
            ImGui::TextColored(ImVec4(1.0f, 0.55f, 0.22f, 1.0f),
                "Selected 3D grid node [%d,%d,%d] = (% .17g, % .17g, % .17g), F = % .17g",
                sample.ix, sample.iy, sample.iz,
                sample.point.x, sample.point.y, sample.point.z, sample.value);
        }
    }
    ImGui::TextDisabled("Click any grid node in a panel. The orange 3D cross marks that exact (x,y,z); the blue grid is its slice plane.");

    const float canvas_size = std::clamp(display.y - 420.0f, 300.0f, 480.0f);
    const float panel_height = canvas_size + 44.0f;
    ImGui::BeginChild("SliceGalleryHorizontalScroll",
                      ImVec2(0.0f, panel_height + 12.0f), true,
                      ImGuiWindowFlags_HorizontalScrollbar);

    if (slice_segment_ranges_.size() == slice_result_data_.plane_count)
    {
        for (std::size_t slice_index = 0; slice_index < slice_result_data_.plane_count; ++slice_index)
        {
            if (slice_index != 0)
                ImGui::SameLine(0.0f, 10.0f);
            ImGui::BeginGroup();
            const double plane_coordinate = gridCoordinate(
                static_cast<int>(slice_index), resolution, slice_result_data_.half_extent);
            const char* normal_axis = iso_slices::axis_name(slice_result_data_.axis);
            ImGui::Text("Slice %zu | %s = %+.17g | grid %zu/%d",
                slice_index + 1, normal_axis, plane_coordinate, slice_index, resolution);
            ImGui::PushID(static_cast<int>(slice_index));
            const ImVec2 canvas_origin = ImGui::GetCursorScreenPos();
            ImGui::InvisibleButton("SliceCanvas", ImVec2(canvas_size, canvas_size));
            if (!ImGui::IsItemVisible())
            {
                // Scrolled out of view: keep the layout slot, skip the
                // (N+1)^2 nodes, labels and contour lines of this panel.
                ImGui::PopID();
                ImGui::EndGroup();
                continue;
            }
            const bool canvas_hovered = ImGui::IsItemHovered();
            const bool canvas_clicked = canvas_hovered && ImGui::IsMouseClicked(ImGuiMouseButton_Left);
            const ImVec2 mouse_position = ImGui::GetIO().MousePos;

            ImDrawList* draw = ImGui::GetWindowDrawList();
            const ImVec2 lo(canvas_origin.x + 5.0f, canvas_origin.y + 5.0f);
            const ImVec2 hi(canvas_origin.x + canvas_size - 5.0f,
                            canvas_origin.y + canvas_size - 5.0f);
            draw->AddRectFilled(lo, hi, IM_COL32(14, 19, 29, 255), 2.0f);
            const float plot_w = hi.x - lo.x;
            const float plot_h = hi.y - lo.y;
            const int guide_step = std::max(1, resolution / 12);
            const ImU32 grid_color = IM_COL32(75, 89, 110, 150);
            for (int g = 0; g <= resolution; g += guide_step)
            {
                const float t = static_cast<float>(g) / static_cast<float>(resolution);
                const float x = lo.x + t * plot_w;
                const float y = hi.y - t * plot_h;
                draw->AddLine(ImVec2(x, lo.y), ImVec2(x, hi.y), grid_color, 1.0f);
                draw->AddLine(ImVec2(lo.x, y), ImVec2(hi.x, y), grid_color, 1.0f);
            }
            if (resolution % guide_step != 0)
            {
                draw->AddLine(ImVec2(hi.x, lo.y), hi, grid_color, 1.0f);
                draw->AddLine(ImVec2(lo.x, lo.y), ImVec2(hi.x, lo.y), grid_color, 1.0f);
            }
            draw->AddRect(lo, hi, IM_COL32(140, 155, 174, 220), 2.0f, 0, 1.0f);
            draw->AddText(ImVec2(lo.x + 3.0f, lo.y + 2.0f), IM_COL32(180,195,211,230), "+v");
            draw->AddText(ImVec2(hi.x - 18.0f, hi.y - 15.0f), IM_COL32(180,195,211,230), "+u");

            const double extent = slice_result_data_.half_extent;
            const double diameter = std::max(2.0 * extent, 1e-12);
            auto project_u = [&](const iso_slices::Point3& p) {
                const double u = slice_result_data_.axis == 0 ? p.y : p.x;
                return lo.x + static_cast<float>((u + extent) / diameter) * plot_w;
            };
            auto project_v = [&](const iso_slices::Point3& p) {
                const double v = slice_result_data_.axis == 2 ? p.y : p.z;
                return hi.y - static_cast<float>((v + extent) / diameter) * plot_h;
            };
            const std::size_t sample_base = slice_index * static_cast<std::size_t>(side) * side;
            auto sample_at = [&](int iu, int iv) -> const iso_slices::GridSample& {
                return slice_result_data_.samples[sample_base +
                    static_cast<std::size_t>(iv) * side + static_cast<std::size_t>(iu)];
            };

            // Numeric isovalues are placed at the actual sampled grid nodes.
            // At high resolutions the canvas is too dense for readable labels;
            // hover any node for its full value and exact 3D coordinate.
            const bool draw_value_labels = show_slice_grid_values_;
            const int label_stride = std::clamp(slice_value_label_stride_, 1, resolution + 1);
            const int point_stride = std::max(1, resolution / 24);
            for (int iv = 0; iv <= resolution; ++iv)
            {
                for (int iu = 0; iu <= resolution; ++iu)
                {
                    const auto& sample = sample_at(iu, iv);
                    const ImVec2 pos(lo.x + static_cast<float>(iu) / resolution * plot_w,
                                     hi.y - static_cast<float>(iv) / resolution * plot_h);
                    const bool selected = slice_selected_grid_point_valid_ &&
                        slice_selected_slice_ == static_cast<int>(slice_index) &&
                        slice_selected_iu_ == iu && slice_selected_iv_ == iv;
                    const bool point_this_node =
                        (iu % point_stride == 0 || iu == resolution) &&
                        (iv % point_stride == 0 || iv == resolution);
                    const bool label_this_node = draw_value_labels &&
                        (iu % label_stride == 0 || iu == resolution) &&
                        (iv % label_stride == 0 || iv == resolution);
                    const ImU32 value_color = !std::isfinite(sample.value) ? IM_COL32(240,80,80,255) :
                        (std::fabs(sample.value) <= 1e-8 ? IM_COL32(255,220,85,255) :
                         (sample.value < 0.0 ? IM_COL32(135,195,255,255) : IM_COL32(255,165,145,255)));
                    if (point_this_node || selected || label_this_node)
                        draw->AddCircleFilled(pos, selected ? 5.0f : 2.4f,
                            selected ? IM_COL32(255,105,35,255) : IM_COL32(225,235,248,235));
                    if (label_this_node)
                    {
                        char label[96];
                        if (std::isfinite(sample.value))
                            std::snprintf(label, sizeof(label),
                                show_slice_full_precision_labels_ ? "%.17g" : "%.9g", sample.value);
                        else
                            std::snprintf(label, sizeof(label), "undefined");
                        const ImVec2 text_size = ImGui::CalcTextSize(label);
                        float text_x = pos.x - text_size.x * 0.5f;
                        text_x = std::clamp(text_x, lo.x + 2.0f, hi.x - text_size.x - 2.0f);
                        float text_y = pos.y - text_size.y - 7.0f;
                        if (text_y < lo.y + 2.0f) text_y = pos.y + 7.0f;
                        text_y = std::clamp(text_y, lo.y + 2.0f, hi.y - text_size.y - 2.0f);
                        const ImVec2 text_min(text_x - 3.0f, text_y - 2.0f);
                        const ImVec2 text_max(text_x + text_size.x + 3.0f, text_y + text_size.y + 2.0f);
                        draw->AddRectFilled(text_min, text_max, IM_COL32(3, 8, 15, 238), 3.0f);
                        draw->AddRect(text_min, text_max, IM_COL32(75, 95, 122, 240), 3.0f, 0, 1.0f);
                        draw->AddText(ImVec2(text_x, text_y), value_color, label);
                    }
                }
            }

            // A click maps to the nearest exact lattice node on this slice.
            if (canvas_clicked && mouse_position.x >= lo.x && mouse_position.x <= hi.x &&
                mouse_position.y >= lo.y && mouse_position.y <= hi.y)
            {
                const double fu = (mouse_position.x - lo.x) / plot_w * resolution;
                const double fv = (hi.y - mouse_position.y) / plot_h * resolution;
                slice_selected_slice_ = static_cast<int>(slice_index);
                slice_selected_iu_ = std::clamp(static_cast<int>(std::lround(fu)), 0, resolution);
                slice_selected_iv_ = std::clamp(static_cast<int>(std::lround(fv)), 0, resolution);
                slice_selected_grid_point_valid_ = true;
                slice_reference_geometry_dirty_ = true;
            }
            if (canvas_hovered && mouse_position.x >= lo.x && mouse_position.x <= hi.x &&
                mouse_position.y >= lo.y && mouse_position.y <= hi.y)
            {
                const int iu = std::clamp(static_cast<int>(std::lround(
                    (mouse_position.x - lo.x) / plot_w * resolution)), 0, resolution);
                const int iv = std::clamp(static_cast<int>(std::lround(
                    (hi.y - mouse_position.y) / plot_h * resolution)), 0, resolution);
                const auto& sample = sample_at(iu, iv);
                ImGui::BeginTooltip();
                ImGui::Text("Grid index: (%d, %d, %d)", sample.ix, sample.iy, sample.iz);
                ImGui::Text("Point: (% .17g, % .17g, % .17g)", sample.point.x, sample.point.y, sample.point.z);
                if (std::isfinite(sample.value)) ImGui::Text("F(x,y,z,k) = % .17g", sample.value);
                else ImGui::TextUnformatted("F(x,y,z,k) = undefined");
                ImGui::EndTooltip();
            }

            const auto range = slice_segment_ranges_[slice_index];
            for (std::size_t seg_index = range.first;
                 seg_index < range.second && seg_index < slice_result_data_.segments.size(); ++seg_index)
            {
                const iso_slices::ContourSegment& segment = slice_result_data_.segments[seg_index];
                if (!std::isfinite(segment.a.x) || !std::isfinite(segment.a.y) ||
                    !std::isfinite(segment.a.z) || !std::isfinite(segment.b.x) ||
                    !std::isfinite(segment.b.y) || !std::isfinite(segment.b.z))
                    continue;
                draw->AddLine(ImVec2(project_u(segment.a), project_v(segment.a)),
                              ImVec2(project_u(segment.b), project_v(segment.b)),
                              IM_COL32(70,232,221,255), 2.0f);
            }
            ImGui::PopID();
            ImGui::EndGroup();
        }
    }
    ImGui::EndChild();
    ImGui::TextDisabled("Label stride controls clutter on dense grids. Hover any node or click it for exact coordinates and full-precision F; CSV preserves all values.");
    draw_snapshot_popup();
    ImGui::End();
}

namespace {
constexpr const char* kSnapshotPopup = "Slice snapshot";
} // namespace

// Prepares the dialog: keeps a previous selection if it still fits the current
// result, otherwise preselects the gallery's selected slice with its two
// neighbours and a 3x3 node window around the selected grid node.
void IsosurfaceViewer::open_snapshot_dialog()
{
    const int resolution = std::max(1, slice_result_data_.resolution);
    if (snapshot_slice_selected_.size() != slice_result_data_.plane_count)
    {
        snapshot_slice_selected_.assign(slice_result_data_.plane_count, 0);
        const int first = std::clamp(slice_selected_slice_ - 1, 0, std::max(0, resolution - 2));
        for (int s = first; s < first + 3 && s < static_cast<int>(snapshot_slice_selected_.size()); ++s)
            snapshot_slice_selected_[static_cast<std::size_t>(s)] = 1;
    }
    auto window_around = [resolution](int center) {
        const int lo = std::clamp(center - 1, 0, std::max(0, resolution - 2));
        return std::array<int, 2>{lo, std::min(resolution, lo + 2)};
    };
    auto fits = [resolution](const std::array<int, 2>& range) {
        return range[0] >= 0 && range[0] < range[1] && range[1] <= resolution;
    };
    if (!fits(snapshot_u_range_)) snapshot_u_range_ = window_around(slice_selected_iu_);
    if (!fits(snapshot_v_range_)) snapshot_v_range_ = window_around(slice_selected_iv_);
    snapshot_status_.clear();
    ImGui::OpenPopup(kSnapshotPopup);
}

iso_slices::SnapshotRequest IsosurfaceViewer::snapshot_request() const
{
    iso_slices::SnapshotRequest request;
    for (std::size_t s = 0; s < snapshot_slice_selected_.size(); ++s)
        if (snapshot_slice_selected_[s]) request.slices.push_back(static_cast<int>(s));
    request.iu_min = snapshot_u_range_[0];
    request.iu_max = snapshot_u_range_[1];
    request.iv_min = snapshot_v_range_[0];
    request.iv_max = snapshot_v_range_[1];
    return request;
}

void IsosurfaceViewer::draw_snapshot_popup()
{
    if (!ImGui::BeginPopupModal(kSnapshotPopup, nullptr, ImGuiWindowFlags_AlwaysAutoResize))
        return;

    const iso_slices::Result& result = slice_result_data_;
    const int resolution = result.resolution;
    const std::array<int, 2> uv = iso_slices::in_plane_axes(result.axis);
    const char* u_name = iso_slices::axis_name(uv[0]);
    const char* v_name = iso_slices::axis_name(uv[1]);

    ImGui::TextWrapped("Choose at least %zu slices and a rectangle of grid nodes. The snapshot shows the F value "
                       "of every node in the rectangle on each selected slice, with the corner coordinates.",
                       iso_slices::kMinSnapshotSlices);
    ImGui::Text("Data: k = %.9g, %s-normal slices, resolution %d", result.k,
                iso_slices::axis_name(result.axis), resolution);
    if (slice_data_stale_)
        ImGui::TextColored(ImVec4(1.0f, 0.65f, 0.18f, 1.0f),
                           "The viewer k changed since extraction; the snapshot uses the extracted slices.");
    ImGui::Separator();

    // Slice selection: a 4-column grid of checkboxes.
    const std::size_t selected = static_cast<std::size_t>(
        std::count(snapshot_slice_selected_.begin(), snapshot_slice_selected_.end(), 1));
    ImGui::Text("Slices: %zu selected", selected);
    ImGui::SameLine();
    if (ImGui::SmallButton("Clear"))
        std::fill(snapshot_slice_selected_.begin(), snapshot_slice_selected_.end(), 0);
    ImGui::BeginChild("SnapshotSlices", ImVec2(560.0f, 160.0f), true);
    if (ImGui::BeginTable("SnapshotSliceTable", 4))
    {
        for (std::size_t s = 0; s < snapshot_slice_selected_.size(); ++s)
        {
            ImGui::TableNextColumn();
            bool checked = snapshot_slice_selected_[s] != 0;
            char label[64];
            std::snprintf(label, sizeof(label), "Slice %zu (%s=%.4g)", s + 1, iso_slices::axis_name(result.axis),
                          gridCoordinate(static_cast<int>(s), resolution, result.half_extent));
            if (ImGui::Checkbox(label, &checked))
                snapshot_slice_selected_[s] = checked ? 1 : 0;
        }
        ImGui::EndTable();
    }
    ImGui::EndChild();

    // Grid window: inclusive index ranges along the in-plane axes.
    char u_label[48], v_label[48];
    std::snprintf(u_label, sizeof(u_label), "u = %s grid index", u_name);
    std::snprintf(v_label, sizeof(v_label), "v = %s grid index", v_name);
    ImGui::SetNextItemWidth(300.0f);
    ImGui::DragIntRange2(u_label, &snapshot_u_range_[0], &snapshot_u_range_[1], 0.1f, 0, resolution,
                         "from %d", "to %d");
    ImGui::SetNextItemWidth(300.0f);
    ImGui::DragIntRange2(v_label, &snapshot_v_range_[0], &snapshot_v_range_[1], 0.1f, 0, resolution,
                         "from %d", "to %d");
    auto coordinate = [&](int index) { return gridCoordinate(index, resolution, result.half_extent); };
    ImGui::Text("Bottom-left: %s = %.6g, %s = %.6g    Top-right: %s = %.6g, %s = %.6g",
                u_name, coordinate(snapshot_u_range_[0]), v_name, coordinate(snapshot_v_range_[0]),
                u_name, coordinate(snapshot_u_range_[1]), v_name, coordinate(snapshot_v_range_[1]));
    ImGui::Text("%d x %d nodes per slice", snapshot_u_range_[1] - snapshot_u_range_[0] + 1,
                snapshot_v_range_[1] - snapshot_v_range_[0] + 1);

    const iso_slices::SnapshotRequest request = snapshot_request();
    const std::string problem = iso_slices::validate_snapshot(result, request);
    const std::filesystem::path directory = project_output_directory("snapshots");
    ImGui::Separator();
    if (problem.empty())
        ImGui::TextWrapped("Saves %s.svg and .csv to %s", iso_slices::snapshot_base_name(result, request).c_str(),
                           directory.string().c_str());
    else
        ImGui::TextColored(ImVec4(1.0f, 0.45f, 0.35f, 1.0f), "%s", problem.c_str());

    ImGui::BeginDisabled(!problem.empty());
    if (ImGui::Button("Save snapshot", ImVec2(160.0f, 0.0f)))
    {
        try
        {
            const iso_slices::SnapshotFiles files = iso_slices::write_snapshot(result, request, directory);
            snapshot_status_ = "Saved " + files.svg.filename().string() + " and " +
                files.csv.filename().string() + " in " + directory.string();
        }
        catch (const std::exception& error)
        {
            snapshot_status_ = std::string("Snapshot failed: ") + error.what();
        }
        runtime_message_ = snapshot_status_;
    }
    ImGui::EndDisabled();
    ImGui::SameLine();
    if (ImGui::Button("Close", ImVec2(100.0f, 0.0f)))
        ImGui::CloseCurrentPopup();
    if (!snapshot_status_.empty())
        ImGui::TextWrapped("%s", snapshot_status_.c_str());
    ImGui::EndPopup();
}

} // namespace iso_app
