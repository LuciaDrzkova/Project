// ImGui control panel of the viewer: the PMP side panel sections and the
// startup setup dialog.

#include "IsosurfaceViewer.h"

#include "core/SurfaceCatalogue.h"

#include <imgui.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <string>
#include <vector>

namespace iso_app {

namespace {

template <std::size_t N>
float array_max(const std::array<float, N>& values, float floor)
{
    float result = floor;
    for (float value : values)
        if (std::isfinite(value)) result = std::max(result, value);
    return result;
}

} // namespace

void IsosurfaceViewer::process_imgui()
{
    poll_batch_job();
    poll_convergence_job();
    // Open the adapted startup setup as a modal in the PMP viewer's
    // standard ImGui panel. This program generates the isosurface
    // directly, so there is no input mesh to load.
    if (setup_open_ && !setup_popup_requested_)
    {
        ImGui::OpenPopup("Isosurface Setup");
        setup_popup_requested_ = true;
    }

    pmp::MeshViewer::process_imgui();

    // The inherited PMP control window can remember an overly wide size
    // from an earlier version. Clamp it to a practical width, while
    // preserving its current height and allowing the user to make it
    // narrower. The fixed-width controls below ensure content no longer
    // forces the window to expand again.
    const ImVec2 panel_size = ImGui::GetWindowSize();
    constexpr float kMaxPanelWidth = 375.0f;
    if (panel_size.x > kMaxPanelWidth)
        ImGui::SetWindowSize(ImVec2(kMaxPanelWidth, panel_size.y), ImGuiCond_Always);

    ImGui::Separator();
    ImGui::TextUnformatted("Isosurface Singularity Lab");
    ImGui::TextWrapped("Explore F(x,y,z,k) = 0. The zero isovalue is fixed; k is the unfolding parameter.");

    if (ImGui::Button("Setup..."))
    {
        setup_open_ = true;
        setup_popup_requested_ = false;
    }

    draw_function_panel();
    draw_sampling_panel();
    draw_slice_panel();
    draw_cell_panel();
    draw_gradient_panel();
    draw_rendering_panel();

    ImGui::Separator();
    ImGui::Text("Surface: %zu triangles", last_triangle_count_);
    if (last_build_ms_ >= 0.0)
        ImGui::Text("Last rebuild: %.2f ms", last_build_ms_);
    if (!runtime_message_.empty())
        ImGui::TextWrapped("%s", runtime_message_.c_str());

    draw_slice_gallery();
    draw_setup_popup();
    draw_batch_experiment_popup();
}

// Expression editor, presets and k.
void IsosurfaceViewer::draw_function_panel()
{
    if (ImGui::CollapsingHeader("Function and k", ImGuiTreeNodeFlags_DefaultOpen))
    {
        draw_formula_editor("##FunctionEditor");
        ImGui::Separator();
        bool k_changed = false;
        k_changed |= draw_k_text_editor("k value", main_k_error_);
        // Keep this main-panel slider compact. A width of -1.0f fills
        // the entire available panel width, which is too wide here.
        ImGui::SetNextItemWidth(210.0f);
        k_changed |= ImGui::SliderFloat("k slider", &k_, -1.0f, 1.0f, "%.4f");
        if (k_changed)
        {
            sanitize_parameters();
            mesh_dirty_ = true;
            invalidate_slice_overlay();
        }
        ImGui::Text("Target isovalue: F = 0");
    }
}

// Grid resolution, domain and color gamma.
void IsosurfaceViewer::draw_sampling_panel()
{
    if (ImGui::CollapsingHeader("Sampling grid", ImGuiTreeNodeFlags_DefaultOpen))
    {
        bool resolution_changed = ImGui::SliderInt(
            "Grid resolution (cells/axis)", &grid_resolution_, 1, 110);
        bool extent_changed = ImGui::SliderFloat(
            "Domain half-extent", &half_extent_, 0.2f, 3.0f, "%.2f");
        bool gamma_changed = ImGui::SliderFloat(
            "Gradient color gamma", &gradient_gamma_, 0.2f, 2.0f, "%.2f");
        if (resolution_changed || extent_changed)
        {
            sanitize_parameters();
            mesh_dirty_ = true;
            grid_geometry_dirty_ = true;
            cell_geometry_dirty_ = true;
            selected_vertex_geometry_dirty_ = true;
            invalidate_slice_overlay();
        }
        if (gamma_changed)
        {
            // Gamma only affects colors: recolor instead of re-extracting.
            sanitize_parameters();
            recolor_surface();
            glyph_geometry_dirty_ = true;
        }

        ImGui::Checkbox("Show sampling grid", &show_grid_);
        ImGui::SameLine();
        ImGui::Checkbox("Highlight selected cell", &show_selected_cell_);
        if (show_grid_)
        {
            const bool stride_changed = ImGui::SliderInt(
                "Grid line stride (1 = every line)", &grid_stride_, 1,
                std::max(1, grid_resolution_));
            if (stride_changed)
                grid_geometry_dirty_ = true;
        }
        ImGui::Text("Grid samples per axis: %d", grid_resolution_ + 1);
        ImGui::Text("Voxel cells: %d x %d x %d",
                    grid_resolution_, grid_resolution_, grid_resolution_);
        ImGui::Text("Voxel spacing: %.5g", (2.0 * half_extent_) / grid_resolution_);
    }
}

// Grid-aligned 2D slice export.
void IsosurfaceViewer::draw_slice_panel()
{
    if (ImGui::CollapsingHeader("2D slice extraction", ImGuiTreeNodeFlags_DefaultOpen))
    {
        static const char* slice_axis_labels[] = {
            "X-normal cuts (YZ planes)",
            "Y-normal cuts (XZ planes)",
            "Z-normal cuts (XY planes)"};
        if (ImGui::Combo("Slice direction", &slice_axis_, slice_axis_labels, 3))
            invalidate_slice_overlay();
        ImGui::TextWrapped("Resolution N creates N+1 cuts exactly on grid coordinates -L+i*h along the selected axis, including both domain endpoints. Each cut uses the original (N+1)x(N+1) in-plane grid, labels F at grid nodes, and extracts the F=0 contour.");
        ImGui::SetNextItemWidth(225.0f);
        ImGui::InputText("Output folder", slice_output_directory_.data(), slice_output_directory_.size());
        if (ImGui::Button("Extract 2D slices + export CSV"))
        {
            try
            {
                apply_slice_result(extract_slices_at_current_parameters());
                slice_selected_slice_ = static_cast<int>(slice_result_data_.plane_count / 2);
                slice_selected_iu_ = slice_result_data_.resolution / 2;
                slice_selected_iv_ = slice_result_data_.resolution / 2;
                runtime_message_ = "2D slices exported to '" + std::string(slice_output_directory_.data()) +
                    "' (" + std::to_string(slice_result_data_.plane_count) + " planes, " +
                    std::to_string(slice_result_data_.contour_segment_count) + " contour segments).";
            }
            catch (const std::exception& error)
            {
                invalidate_slice_overlay();
                runtime_message_ = std::string("2D slice extraction failed: ") + error.what();
            }
        }
        if (ImGui::Button("Clear 2D contours"))
            clear_slice_overlay();
        ImGui::Checkbox("Show extracted contours in 3D", &show_slice_contours_);
        ImGui::Checkbox("Show side-by-side 2D slice gallery", &show_slice_gallery_);
        if (slice_result_available_)
            ImGui::Text("Last result: %zu planes, %zu grid values, %zu contour segments",
                        slice_result_data_.plane_count, slice_result_data_.sample_count, slice_result_data_.contour_segment_count);
        ImGui::TextDisabled("CSV output includes every sampled grid value, each contour segment endpoint, and the four corner isovalues of its 2D cell.");
    }
}

// Corner values of the selected grid cell.
void IsosurfaceViewer::draw_cell_panel()
{
    if (ImGui::CollapsingHeader("Selected cell corner isovalues", ImGuiTreeNodeFlags_DefaultOpen))
    {
        ImGui::TextWrapped("Choose a cell to inspect F at its eight grid corners. Coordinates and values are evaluated at the current k.");
        bool cell_changed = false;
        cell_changed |= ImGui::SliderInt("Cell i (x)", &selected_cell_[0], 0, grid_resolution_ - 1);
        cell_changed |= ImGui::SliderInt("Cell j (y)", &selected_cell_[1], 0, grid_resolution_ - 1);
        cell_changed |= ImGui::SliderInt("Cell l (z)", &selected_cell_[2], 0, grid_resolution_ - 1);
        if (cell_changed)
        {
            has_selected_grid_vertex_ = false;
            cell_geometry_dirty_ = true;
            selected_vertex_geometry_dirty_ = true;
        }

        GridCorner corners[8];
        get_cell_corners(corners);
        if (ImGui::BeginTable("GridCellCornerValues", 4,
                              ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                              ImGuiTableFlags_SizingStretchProp))
        {
            ImGui::TableSetupColumn("Corner");
            ImGui::TableSetupColumn("Position (x,y,z)");
            ImGui::TableSetupColumn("F value");
            ImGui::TableSetupColumn("Sign vs 0");
            ImGui::TableHeadersRow();
            for (const GridCorner& corner : corners)
            {
                ImGui::TableNextRow();
                bool selected_corner = has_selected_grid_vertex_ &&
                    corner.ix == selected_grid_vertex_[0] &&
                    corner.iy == selected_grid_vertex_[1] &&
                    corner.iz == selected_grid_vertex_[2];
                if (selected_corner)
                    ImGui::TableSetBgColor(ImGuiTableBgTarget_RowBg0,
                                           IM_COL32(190, 70, 30, 135));
                ImGui::TableSetColumnIndex(0);
                // Make the whole corner row clickable. Selecting a row in this
                // table also selects that exact grid vertex in the 3D viewer.
                if (ImGui::Selectable(corner.name, selected_corner,
                                      ImGuiSelectableFlags_SpanAllColumns |
                                      ImGuiSelectableFlags_AllowOverlap))
                {
                    selected_grid_vertex_ = {corner.ix, corner.iy, corner.iz};
                    has_selected_grid_vertex_ = true;
                    selected_vertex_geometry_dirty_ = true;
                    selected_corner = true;
                }
                ImGui::TableNextColumn();
                ImGui::Text("(% .3f, % .3f, % .3f)", corner.x, corner.y, corner.z);
                ImGui::TableNextColumn();
                if (std::isfinite(corner.value))
                    ImGui::Text("% .6g", corner.value);
                else
                    ImGui::TextUnformatted("undefined");
                ImGui::TableNextColumn();
                if (!std::isfinite(corner.value))
                    ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "undefined");
                else if (std::fabs(corner.value) < 1e-6)
                    ImGui::TextColored(ImVec4(1.0f, 0.78f, 0.2f, 1.0f), "~= 0");
                else if (corner.value < 0.0)
                    ImGui::TextColored(ImVec4(0.35f, 0.65f, 1.0f, 1.0f), "negative");
                else
                    ImGui::TextColored(ImVec4(1.0f, 0.45f, 0.35f, 1.0f), "positive");
            }
            ImGui::EndTable();
        }
        ImGui::TextDisabled("The cyan outline identifies this cell; its corners are white, and the selected corner is red-orange in the 3D viewer.");
        ImGui::TextWrapped("Click any corner row in this table to highlight that vertex in the 3D view. You can also click a visible grid intersection in the viewer.");
        if (has_selected_grid_vertex_)
        {
            const double sx = gridCoordinate(selected_grid_vertex_[0], grid_resolution_, half_extent_);
            const double sy = gridCoordinate(selected_grid_vertex_[1], grid_resolution_, half_extent_);
            const double sz = gridCoordinate(selected_grid_vertex_[2], grid_resolution_, half_extent_);
            const double f = evalExpr<double>(ast_.get(), sx, sy, sz, k_);
            ImGui::TextColored(ImVec4(1.0f, 0.62f, 0.15f, 1.0f),
                "Selected vertex [%d, %d, %d]: (% .4f, % .4f, % .4f)",
                selected_grid_vertex_[0], selected_grid_vertex_[1], selected_grid_vertex_[2],
                sx, sy, sz);
            if (std::isfinite(f))
                ImGui::Text("F(x,y,z,k) = % .8g  (k = %.4f)", f, k_);
            else
                ImGui::TextUnformatted("F(x,y,z,k) = undefined");
        }
        else
            ImGui::TextDisabled("Selected vertex: none. Click a corner row here or a grid intersection in the 3D view.");
    }
}

// Gradient at P, the k sweep plots and the convergence/batch tools.
void IsosurfaceViewer::draw_gradient_panel()
{
    if (ImGui::CollapsingHeader("Gradient diagnostics", false))
    {
        ImGui::TextUnformatted("Point of interest P: used for the k-sweep plots");
        if (ImGui::InputFloat3("P", point_of_interest_.data(), "%.3f"))
        {
            sweep_dirty_ = true;
            point_of_interest_geometry_dirty_ = true;
        }
        if (ImGui::Button("Snap P to origin"))
        {
            point_of_interest_ = {0.0f, 0.0f, 0.0f};
            sweep_dirty_ = true;
            point_of_interest_geometry_dirty_ = true;
        }
        Dual3 x(point_of_interest_[0], 1, 0, 0);
        Dual3 y(point_of_interest_[1], 0, 1, 0);
        Dual3 z(point_of_interest_[2], 0, 0, 1);
        Dual3 kk(k_, 0, 0, 0);
        const Dual3 value = evalExpr<Dual3>(ast_.get(), x, y, z, kk);
        const double grad_norm = std::sqrt(value.dx * value.dx + value.dy * value.dy + value.dz * value.dz);
        ImGui::Text("At P, k=%.4f: F = % .5g, |grad F| = %.5g", k_, value.v, grad_norm);
        ImGui::Text("grad F = (% .5g, % .5g, % .5g)", value.dx, value.dy, value.dz);
        if (!parse_error_.empty())
            ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "Parse error: %s", parse_error_.c_str());

        ImGui::Separator();
        // Convergence analysis and batch tools.
        const float max_grad = array_max(sweep_.grad_magnitude, 1e-5f);
        ImGui::Text("|grad F(P,k)| over k in [-1,1]");
        ImGui::SameLine();
        if (ImGui::SmallButton("Analyze surface convergence"))
            start_convergence_job();
        ImGui::SameLine();
        if (ImGui::SmallButton("Batch 25 F + kG tests..."))
        {
            batch_setup_open_ = true;
            batch_popup_requested_ = false;
        }

        // Select the actual k samples for the geometric convergence experiment.
        ImGui::SetNextItemWidth(92.0f);
        ImGui::InputDouble("k minimum", &convergence_k_min_, 0.0, 0.0, "%.8g");
        ImGui::SameLine();
        ImGui::SetNextItemWidth(92.0f);
        ImGui::InputDouble("k maximum", &convergence_k_max_, 0.0, 0.0, "%.8g");
        ImGui::SetNextItemWidth(145.0f);
        ImGui::InputInt("Number of k samples", &convergence_sample_count_);
        convergence_sample_count_ = std::clamp(convergence_sample_count_, 3, 101);
        ImGui::TextDisabled("Uniformly spaced; endpoints included. Max 101 samples.");
        ImGui::SetNextItemWidth(210.0f);
        ImGui::SliderInt("Convergence mesh resolution", &convergence_resolution_, 16, 110);
        draw_singularity_convergence_report();

        plot_with_cursor("##gradmag", sweep_.grad_magnitude.data(), KSweep::N, k_, 0.0f, max_grad * 1.05f, 75.0f);
        if (ImGui::CollapsingHeader("Gradient components and F(P,k)", false))
        {
            float lo = std::numeric_limits<float>::max();
            float hi = std::numeric_limits<float>::lowest();
            for (int i = 0; i < KSweep::N; ++i)
            {
                lo = std::min({lo, sweep_.gx[i], sweep_.gy[i], sweep_.gz[i]});
                hi = std::max({hi, sweep_.gx[i], sweep_.gy[i], sweep_.gz[i]});
            }
            if (!(hi > lo)) { lo -= 1.0f; hi += 1.0f; }
            ImGui::TextUnformatted("dF/dx"); plot_with_cursor("##gx", sweep_.gx.data(), KSweep::N, k_, lo, hi, 50.0f);
            ImGui::TextUnformatted("dF/dy"); plot_with_cursor("##gy", sweep_.gy.data(), KSweep::N, k_, lo, hi, 50.0f);
            ImGui::TextUnformatted("dF/dz"); plot_with_cursor("##gz", sweep_.gz.data(), KSweep::N, k_, lo, hi, 50.0f);
            float flo = std::numeric_limits<float>::max();
            float fhi = std::numeric_limits<float>::lowest();
            for (float v : sweep_.f_value) { flo = std::min(flo, v); fhi = std::max(fhi, v); }
            if (!(fhi > flo)) { flo -= 1.0f; fhi += 1.0f; }
            ImGui::TextUnformatted("F(P,k)"); plot_with_cursor("##fvalue", sweep_.f_value.data(), KSweep::N, k_, flo, fhi, 50.0f);
        }
    }
}

// Gradient glyph options.
void IsosurfaceViewer::draw_rendering_panel()
{
    if (ImGui::CollapsingHeader("Rendering", false))
    {
        if (ImGui::Checkbox("Show gradient glyphs", &show_gradient_glyphs_))
            glyph_geometry_dirty_ = true;
        bool glyph_changed = false;
        glyph_changed |= ImGui::SliderFloat("Glyph length", &glyph_scale_, 0.02f, 0.6f, "%.2f");
        glyph_changed |= ImGui::SliderInt("Glyph stride", &glyph_stride_, 1, 80);
        if (glyph_changed)
            glyph_geometry_dirty_ = true;
        ImGui::TextDisabled("Surface color: blue/purple = small |grad F|, red = steep.");
    }
}

void IsosurfaceViewer::draw_formula_editor(const char* id_suffix)
{
    ImGui::TextWrapped("Expression variables: x, y, z, k. Operators: + - * / ^. Functions: sin cos tan exp log sqrt abs pow min max.");

    // Use explicit widths here. In particular, an implicit-width InputText
    // beside an Apply button and a long preset label made ImGui grow the
    // entire PMP information panel when this section was expanded.
    ImGui::TextUnformatted("F(x,y,z,k)");
    ImGui::SetNextItemWidth(218.0f);
    const std::string expression_id = std::string("##Expression") + id_suffix;
    const bool entered = ImGui::InputText(expression_id.c_str(), expression_buffer_.data(),
                                          expression_buffer_.size(), ImGuiInputTextFlags_EnterReturnsTrue);
    ImGui::SameLine(0.0f, 6.0f);
    const std::string apply_id = std::string("Apply") + id_suffix;
    const bool apply_clicked = ImGui::Button(apply_id.c_str(), ImVec2(58.0f, 0.0f));

            // The viewer and diagnostics tests consume the same parameterized catalogue.
    const auto& surface_presets = iso_surfaces::catalogue();
    static std::vector<const char*> preset_names;
    if (preset_names.empty()) {
        preset_names.reserve(surface_presets.size());
        for (const auto& preset : surface_presets)
            preset_names.push_back(preset.label.c_str());
    }
    if (preset_index_ < 0 || preset_index_ >= static_cast<int>(surface_presets.size()))
        preset_index_ = 0;
    ImGui::TextUnformatted("Preset");
    ImGui::SetNextItemWidth(282.0f);
    const std::string preset_id = std::string("##Preset") + id_suffix;
    if (ImGui::Combo(preset_id.c_str(), &preset_index_, preset_names.data(),
                     static_cast<int>(preset_names.size())))
    {
        const auto& selected = surface_presets[static_cast<std::size_t>(preset_index_)];
        std::snprintf(expression_buffer_.data(), expression_buffer_.size(), "%s",
                      selected.expression.c_str());
        // Some positive-even leading forms collapse to a point at k=0;
        // their preset-specific nonzero default makes the surface visible.
        k_ = static_cast<float>(selected.defaultK);
        apply_expression();
    }
    ImGui::TextWrapped("%s", surface_presets[static_cast<std::size_t>(preset_index_)].note.c_str());

    if (entered || apply_clicked)
        apply_expression();
    if (!parse_error_.empty())
    {
        ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "Parse error: %s", parse_error_.c_str());
        ImGui::TextDisabled("Showing the last valid expression: %s", last_good_expression_.c_str());
    }
}

// Edits k directly as a scientific-notation-friendly number field.
bool IsosurfaceViewer::draw_k_text_editor(const char* label, std::string& error_message)
{
    ImGui::SetNextItemWidth(130.0f);
    const bool changed = ImGui::InputScalar(
        label, ImGuiDataType_Float, &k_, nullptr, nullptr, "%.9g",
        ImGuiInputTextFlags_CharsScientific);

    if (changed)
    {
        error_message.clear();
        if (!std::isfinite(k_))
        {
            k_ = 0.0f;
            error_message = "Invalid k. Enter a finite number, e.g. -0.001 or 1e-4.";
        }
    }

    if (!error_message.empty())
        ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f),
                           "%s", error_message.c_str());
    return changed;
}

void IsosurfaceViewer::draw_setup_popup()
{
    if (!setup_open_)
        return;
    if (ImGui::BeginPopupModal("Isosurface Setup", &setup_open_, ImGuiWindowFlags_AlwaysAutoResize))
    {
        ImGui::TextWrapped("Configure the implicit function and sampling domain, then open the interactive viewer. This app generates F(x,y,z,k)=0 directly; no input mesh is needed.");
        ImGui::Separator();
        draw_formula_editor("##SetupEditor");
        ImGui::SetNextItemWidth(110.0f);
        bool changed = ImGui::InputFloat("k value", &k_, 0.01f, 0.1f, "%.4f");
        ImGui::SetNextItemWidth(210.0f);
        changed |= ImGui::SliderFloat("k slider", &k_, -1.0f, 1.0f, "%.4f");
        changed |= ImGui::SliderInt("Grid resolution", &grid_resolution_, 1, kMaxGridResolution);
        changed |= ImGui::SliderFloat("Domain half-extent", &half_extent_, 0.2f, 3.0f, "%.2f");
        ImGui::Checkbox("Show sampling grid", &show_grid_);
        if (changed)
        {
            sanitize_parameters();
            mesh_dirty_ = true;
            sweep_dirty_ = true;
            grid_geometry_dirty_ = true;
            cell_geometry_dirty_ = true;
            selected_vertex_geometry_dirty_ = true;
            invalidate_slice_overlay();
        }
        if (ImGui::Button("Open Viewer", ImVec2(160.0f, 0.0f)))
        {
            sanitize_parameters();
            mesh_dirty_ = true;
            grid_geometry_dirty_ = true;
            cell_geometry_dirty_ = true;
            setup_open_ = false;
            ImGui::CloseCurrentPopup();
        }
        ImGui::SameLine();
        if (ImGui::Button("Use defaults"))
        {
            std::snprintf(expression_buffer_.data(), expression_buffer_.size(), "%s", "x*x*x - k*y");
            preset_index_ = 0;
            k_ = 0.0f;
            grid_resolution_ = 48;
            half_extent_ = 1.0f;
            grid_stride_ = 4;
            show_grid_ = false;
            apply_expression();
            sanitize_parameters();
            grid_geometry_dirty_ = true;
            cell_geometry_dirty_ = true;
            selected_vertex_geometry_dirty_ = true;
            mesh_dirty_ = true;
        }
        ImGui::EndPopup();
    }
}

void IsosurfaceViewer::plot_with_cursor(const char* label, const float* ys, int count,
                                        float current_k, float y_min, float y_max, float height)
{
    ImGui::PlotLines(label, ys, count, 0, nullptr, y_min, y_max, ImVec2(0, height));
    const ImVec2 min = ImGui::GetItemRectMin();
    const ImVec2 max = ImGui::GetItemRectMax();
    const float t = std::clamp((current_k + 1.0f) * 0.5f, 0.0f, 1.0f);
    const float x = min.x + t * (max.x - min.x);
    ImGui::GetWindowDrawList()->AddLine(
        ImVec2(x, min.y), ImVec2(x, max.y), IM_COL32(255, 255, 255, 180), 1.5f);
}

} // namespace iso_app
