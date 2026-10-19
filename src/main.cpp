// Isosurface Singularity Lab — complete PMP viewer application.
// The full viewer implementation and main() are kept in this file.

// PMP-based viewer layer for Isosurface Singularity Lab.
// The expression parser and marching-tetrahedra extractor remain the ones
// from the original Project; only the window/viewer layer uses pmp-library.

#include "Expr.h"
#include "Isosurface.h"
#include "SurfaceCatalogue.h"
#include "IsosurfaceSlices.h"

#include <GL/glew.h>
#include <GLFW/glfw3.h>
#include <pmp/visualization/gl.h>
#include <pmp/visualization/mesh_viewer.h>
#include <pmp/visualization/shader.h>
#include <imgui.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <chrono>
#include <exception>
#include <iostream>
#include <cstddef>
#include <cstdio>
#include <cstdlib>
#include <cctype>
#include <cstring>
#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

namespace iso
{

class MyViewer : public pmp::MeshViewer
{
public:
    MyViewer(const char* title, int width, int height)
        : pmp::MeshViewer(title, width, height)
    {
        std::snprintf(expression_buffer_.data(), expression_buffer_.size(),
                      "%s", "x*x*x - k*y");
        std::snprintf(slice_output_directory_.data(), slice_output_directory_.size(),
                      "%s", "isosurface_slices");
        ast_ = parseExpression(expression_buffer_.data());
        add_help_item("Left drag", "Orbit the isosurface");
        add_help_item("Click grid vertex", "Select a sample vertex and inspect its F value");
        add_help_item("Scroll", "Zoom the viewer");
        add_help_item("?", "Show key bindings");
        setup_open_ = true;
    }

    ~MyViewer() override
    {
        // The derived object is destroyed before the base Window destroys the
        // OpenGL context, so it is safe to release these GPU resources here.
        delete_overlay_buffer(grid_buffer_);
        delete_overlay_buffer(cell_line_buffer_);
        delete_overlay_buffer(cell_point_buffer_);
        delete_overlay_buffer(glyph_buffer_);
        delete_overlay_buffer(point_of_interest_buffer_);
        delete_overlay_buffer(selected_grid_vertex_buffer_);
        delete_overlay_buffer(slice_buffer_);
        delete_overlay_buffer(slice_plane_buffer_);
        delete_overlay_buffer(slice_selected_point_buffer_);
    }

protected:
    void do_processing() override
    {
        if (mesh_dirty_)
            rebuild_surface();
        if (sweep_dirty_)
            rebuild_k_sweep();
    }

    // Select a visible sampling-grid vertex with a click (not a drag). PMP's
    // Window filters clicks over ImGui controls before dispatching this callback.
    void mouse(int button, int action, int mods) override
    {
        pmp::MeshViewer::mouse(button, action, mods);

        if (button != GLFW_MOUSE_BUTTON_LEFT)
            return;

        double x = 0.0, y = 0.0;
        cursor_pos(x, y); // PMP returns framebuffer-scaled coordinates.

        if (action == GLFW_PRESS)
        {
            selection_click_candidate_ = true;
            selection_click_start_ = {x, y};
            return;
        }

        if (action == GLFW_RELEASE && selection_click_candidate_)
        {
            selection_click_candidate_ = false;
            const double dx = x - selection_click_start_[0];
            const double dy = y - selection_click_start_[1];
            const double click_tolerance = 5.0 * std::max(1.0f, high_dpi_scaling());
            if (dx * dx + dy * dy <= click_tolerance * click_tolerance)
                select_grid_vertex_near(x, y);
        }
    }

    void process_imgui() override
    {
        // Open the adapted startup setup as a modal in the PMP viewer's
        // standard ImGui panel. Unlike the remeshing project, this program
        // generates an isosurface directly and therefore has no input mesh.
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

        if (ImGui::CollapsingHeader("Function and k", ImGuiTreeNodeFlags_DefaultOpen))
        {
            draw_formula_editor("##FunctionEditor");
            ImGui::Separator();
            bool k_changed = false;
            ImGui::SetNextItemWidth(110.0f);
            k_changed |= draw_k_text_editor("k value", main_k_text_buffer_, main_k_text_editing_, "Apply k##main", main_k_text_error_);
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
                mesh_dirty_ = true;

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
                    sanitize_parameters();
                    const iso_slices::Result slice_result = iso_slices::extract_and_export(
                        ast_.get(), last_good_expression_, static_cast<double>(k_),
                        grid_resolution_, static_cast<double>(half_extent_), slice_axis_,
                        slice_output_directory_.data());
                    slice_result_data_ = slice_result;
                    slice_data_stale_ = false;
                    slice_selected_slice_ = static_cast<int>(slice_result_data_.plane_count / 2);
                    slice_selected_iu_ = slice_result_data_.resolution / 2;
                    slice_selected_iv_ = slice_result_data_.resolution / 2;
                    slice_selected_grid_point_valid_ = !slice_result_data_.samples.empty();
                    rebuild_slice_segment_ranges();
                    slice_reference_geometry_dirty_ = true;
                    show_slice_gallery_ = true;
                    slice_contour_vertices_.clear();
                    const double denom = static_cast<double>(std::max(1, grid_resolution_));
                    for (const iso_slices::ContourSegment& segment : slice_result.segments)
                    {
                        const float t = static_cast<float>(segment.slice_index / denom);
                        const pmp::Color color = gradient_color(t);
                        slice_contour_vertices_.push_back(overlay_vertex(
                            segment.a.x, segment.a.y, segment.a.z, color[0], color[1], color[2]));
                        slice_contour_vertices_.push_back(overlay_vertex(
                            segment.b.x, segment.b.y, segment.b.z, color[0], color[1], color[2]));
                    }
                    slice_plane_count_ = slice_result.plane_count;
                    slice_sample_count_ = slice_result.sample_count;
                    slice_segment_count_ = slice_result.contour_segment_count;
                    slice_result_available_ = true;
                    show_slice_contours_ = true;
                    slice_geometry_dirty_ = true;
                    runtime_message_ = "2D slices exported to '" + std::string(slice_output_directory_.data()) +
                        "' (" + std::to_string(slice_plane_count_) + " planes, " +
                        std::to_string(slice_segment_count_) + " contour segments).";
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
                            slice_plane_count_, slice_sample_count_, slice_segment_count_);
            ImGui::TextDisabled("CSV output includes every sampled grid value, each contour segment endpoint, and the four corner isovalues of its 2D cell.");
        }

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
                const double sx = grid_coord(selected_grid_vertex_[0], grid_resolution_, half_extent_);
                const double sy = grid_coord(selected_grid_vertex_[1], grid_resolution_, half_extent_);
                const double sz = grid_coord(selected_grid_vertex_[2], grid_resolution_, half_extent_);
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
            const float max_grad = array_max(sweep_.grad_magnitude, 1e-5f);
            ImGui::Text("|grad F(P,k)| over k in [-1,1]");
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

        ImGui::Separator();
        ImGui::Text("Surface: %zu triangles", last_triangle_count_);
        if (last_build_ms_ >= 0.0)
            ImGui::Text("Last rebuild: %.2f ms", last_build_ms_);
        if (!runtime_message_.empty())
            ImGui::TextWrapped("%s", runtime_message_.c_str());

        draw_slice_gallery();
        draw_setup_popup();
    }

    void draw(const std::string& draw_mode) override
    {
        pmp::MeshViewer::draw(draw_mode);
        if (!overlay_initialized_)
            initialize_overlay();

        if (grid_geometry_dirty_)
        {
            build_grid_geometry();
            upload_overlay_buffer(grid_buffer_, grid_vertices_);
            grid_geometry_dirty_ = false;
        }
        if (cell_geometry_dirty_)
        {
            build_cell_geometry();
            upload_overlay_buffer(cell_line_buffer_, cell_line_vertices_);
            upload_overlay_buffer(cell_point_buffer_, cell_point_vertices_);
            cell_geometry_dirty_ = false;
        }
        if (glyph_geometry_dirty_)
        {
            build_glyph_geometry();
            upload_overlay_buffer(glyph_buffer_, glyph_vertices_);
            glyph_geometry_dirty_ = false;
        }
        if (point_of_interest_geometry_dirty_)
        {
            build_point_of_interest_geometry();
            upload_overlay_buffer(point_of_interest_buffer_, point_of_interest_vertices_);
            point_of_interest_geometry_dirty_ = false;
        }
        if (selected_vertex_geometry_dirty_)
        {
            build_selected_grid_vertex_geometry();
            upload_overlay_buffer(selected_grid_vertex_buffer_, selected_grid_vertex_vertices_);
            selected_vertex_geometry_dirty_ = false;
        }

        if (slice_geometry_dirty_)
        {
            upload_overlay_buffer(slice_buffer_, slice_contour_vertices_);
            slice_geometry_dirty_ = false;
        }
        if (slice_reference_geometry_dirty_)
        {
            build_slice_reference_geometry();
            upload_overlay_buffer(slice_plane_buffer_, slice_plane_vertices_);
            upload_overlay_buffer(slice_selected_point_buffer_, slice_selected_point_vertices_);
            slice_reference_geometry_dirty_ = false;
        }
        overlay_shader_.use();
        overlay_shader_.set_uniform("projection_matrix", projection_matrix_);
        overlay_shader_.set_uniform("modelview_matrix", modelview_matrix_);
        overlay_shader_.set_uniform("point_size", supported_point_size(1.0f));
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        glDepthMask(GL_FALSE);
        glEnable(GL_PROGRAM_POINT_SIZE);

        if (show_grid_ && !grid_vertices_.empty())
        {
            overlay_shader_.set_uniform("alpha", 0.62f);
            glLineWidth(supported_line_width(1.0f));
            glBindVertexArray(grid_buffer_.vao);
            glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(grid_vertices_.size()));
        }
        if (show_slice_contours_ && slice_result_available_ && !slice_data_stale_ && !slice_contour_vertices_.empty())
        {
            overlay_shader_.set_uniform("alpha", 1.0f);
            glLineWidth(supported_line_width(2.0f));
            GLint previous_depth_func = GL_LESS;
            glGetIntegerv(GL_DEPTH_FUNC, &previous_depth_func);
            glDepthFunc(GL_LEQUAL); // contours lie directly on F=0 and may share depth with the surface
            glBindVertexArray(slice_buffer_.vao);
            glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(slice_contour_vertices_.size()));
            glDepthFunc(static_cast<GLenum>(previous_depth_func));
        }
        if (show_slice_grid_in_3d_ && slice_result_available_ && !slice_data_stale_ && !slice_plane_vertices_.empty())
        {
            overlay_shader_.set_uniform("alpha", 0.82f);
            glLineWidth(supported_line_width(1.2f));
            GLint previous_depth_func = GL_LESS;
            glGetIntegerv(GL_DEPTH_FUNC, &previous_depth_func);
            glDepthFunc(GL_LEQUAL);
            glBindVertexArray(slice_plane_buffer_.vao);
            glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(slice_plane_vertices_.size()));
            glDepthFunc(static_cast<GLenum>(previous_depth_func));
        }
        if (show_slice_selected_point_in_3d_ && slice_result_available_ && !slice_data_stale_ && !slice_selected_point_vertices_.empty())
        {
            const GLboolean depth_was_enabled = glIsEnabled(GL_DEPTH_TEST);
            glDisable(GL_DEPTH_TEST); // make the matching 3D grid node visible while rotating
            overlay_shader_.set_uniform("alpha", 1.0f);
            glLineWidth(supported_line_width(3.5f));
            glBindVertexArray(slice_selected_point_buffer_.vao);
            glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(slice_selected_point_vertices_.size()));
            if (depth_was_enabled) glEnable(GL_DEPTH_TEST);
        }
        if (show_gradient_glyphs_ && !glyph_vertices_.empty())
        {
            overlay_shader_.set_uniform("alpha", 0.95f);
            glLineWidth(supported_line_width(1.3f));
            glBindVertexArray(glyph_buffer_.vao);
            glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(glyph_vertices_.size()));
        }
        if (show_selected_cell_ && !cell_line_vertices_.empty())
        {
            overlay_shader_.set_uniform("alpha", 1.0f);
            glLineWidth(supported_line_width(2.5f));
            glBindVertexArray(cell_line_buffer_.vao);
            glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(cell_line_vertices_.size()));
            overlay_shader_.set_uniform("point_size", supported_point_size(7.0f));
            glBindVertexArray(cell_point_buffer_.vao);
            glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(cell_point_vertices_.size()));
        }
        if (!point_of_interest_vertices_.empty())
        {
            overlay_shader_.set_uniform("alpha", 1.0f);
            overlay_shader_.set_uniform("point_size", supported_point_size(9.0f));
            glBindVertexArray(point_of_interest_buffer_.vao);
            glDrawArrays(GL_POINTS, 0, static_cast<GLsizei>(point_of_interest_vertices_.size()));
        }

        // Draw the selected grid node on top of the surface so it remains easy
        // to locate after rotating the model, even when that node is occluded.
        if (has_selected_grid_vertex_ && !selected_grid_vertex_vertices_.empty())
        {
            glDisable(GL_DEPTH_TEST);
            overlay_shader_.set_uniform("alpha", 1.0f);
            // This marker is a real triangle-mesh sphere, not a GL_POINTS sprite.
            // It stays visible even on drivers that support only 1-pixel points.
            glBindVertexArray(selected_grid_vertex_buffer_.vao);
            glDrawArrays(GL_TRIANGLES, 0, static_cast<GLsizei>(selected_grid_vertex_vertices_.size()));
            glEnable(GL_DEPTH_TEST);
        }

        glBindVertexArray(0);
        glDepthMask(GL_TRUE);
        glDisable(GL_BLEND);
        glDisable(GL_PROGRAM_POINT_SIZE);
        overlay_shader_.disable();
    }

private:
    struct GridCorner
    {
        const char* name = "";
        int ix = 0, iy = 0, iz = 0;
        double x = 0.0, y = 0.0, z = 0.0, value = 0.0;
    };

    struct OverlayVertex
    {
        float x, y, z;
        float r, g, b;
    };

    struct OverlayBuffer
    {
        GLuint vao = 0;
        GLuint vbo = 0;
    };

    struct KSweep
    {
        static constexpr int N = 201;
        std::array<float, N> grad_magnitude{};
        std::array<float, N> gx{};
        std::array<float, N> gy{};
        std::array<float, N> gz{};
        std::array<float, N> f_value{};
    };

    static constexpr int kMaxGridResolution = 110;
    std::array<char, 512> expression_buffer_{};
    std::string last_good_expression_ = "x*x*x - k*y";
    std::string parse_error_;
    std::string runtime_message_;
    NodePtr ast_;

    float k_ = 0.0f;
    int grid_resolution_ = 48; // number of cells; there are resolution+1 samples
    float half_extent_ = 1.0f;
    float gradient_gamma_ = 0.6f;
    std::array<float, 3> point_of_interest_{0.0f, 0.0f, 0.0f};
    int selected_cell_[3] = {24, 24, 24};
    std::array<int, 3> selected_grid_vertex_{0, 0, 0};
    std::array<double, 2> selection_click_start_{0.0, 0.0};
    int grid_stride_ = 4;
    float glyph_scale_ = 0.15f;
    int glyph_stride_ = 23;

    bool show_grid_ = false;
    bool show_selected_cell_ = true;
    bool show_gradient_glyphs_ = false;
    bool setup_open_ = false;
    bool setup_popup_requested_ = false;
    bool mesh_dirty_ = true;
    bool sweep_dirty_ = true;
    bool grid_geometry_dirty_ = true;
    bool cell_geometry_dirty_ = true;
    bool glyph_geometry_dirty_ = true;
    bool point_of_interest_geometry_dirty_ = true;
    bool selected_vertex_geometry_dirty_ = true;
    bool slice_geometry_dirty_ = true;
    bool slice_result_available_ = false;
    bool show_slice_contours_ = true;
    bool show_slice_gallery_ = false;
    bool show_slice_grid_values_ = true;
    bool show_slice_full_precision_labels_ = false;
    int slice_value_label_stride_ = 1;
    int slice_label_resolution_ = -1;
    bool slice_data_stale_ = false;
    bool show_slice_grid_in_3d_ = true;
    bool show_slice_selected_point_in_3d_ = true;
    bool slice_selected_grid_point_valid_ = false;
    bool slice_reference_geometry_dirty_ = true;
    int slice_selected_slice_ = 0;
    int slice_selected_iu_ = 0;
    int slice_selected_iv_ = 0;
    bool has_selected_grid_vertex_ = false;
    int slice_axis_ = 2;
    std::array<char, 512> slice_output_directory_{};
    std::size_t slice_plane_count_ = 0;
    std::size_t slice_sample_count_ = 0;
    std::size_t slice_segment_count_ = 0;
    bool selection_click_candidate_ = false;
    bool overlay_initialized_ = false;
    bool first_mesh_built_ = false;
    double last_build_ms_ = -1.0;
    std::size_t last_triangle_count_ = 0;

    IsosurfaceResult last_surface_result_;
    KSweep sweep_;

    pmp::Shader overlay_shader_;
    // Many Core Profile drivers (not just macOS: this also reproduces on
    // Mesa/Linux) report a generous GL_ALIASED_LINE_WIDTH_RANGE (e.g. up to
    // 255) purely for legacy/compatibility-profile compatibility, but then
    // raise GL_INVALID_VALUE from glLineWidth() for any width other than
    // 1.0 once the context is actually Core Profile. Trusting that range
    // query alone (as a prior version of this code did) silently breaks:
    // every glLineWidth(>1) call fails, leaving the line width at its
    // last-good value and spamming "GL error: invalid value" every frame
    // forever. Instead, actually try setting a wide line once at startup
    // and see whether the driver accepts it; fall back to 1px if not. This
    // mirrors the existing fallback for the selected-vertex marker, which
    // already avoids GL_POINTS for the same class of reason (point sprite
    // size restrictions) by drawing a real triangle-mesh sphere instead.
    bool wide_lines_supported_ = false;
    std::array<GLfloat, 2> line_width_range_{1.0f, 1.0f};
    std::array<GLfloat, 2> point_size_range_{1.0f, 1.0f};
    OverlayBuffer grid_buffer_;
    OverlayBuffer cell_line_buffer_;
    OverlayBuffer cell_point_buffer_;
    OverlayBuffer glyph_buffer_;
    OverlayBuffer point_of_interest_buffer_;
    OverlayBuffer selected_grid_vertex_buffer_;
    OverlayBuffer slice_buffer_;
    OverlayBuffer slice_plane_buffer_;
    OverlayBuffer slice_selected_point_buffer_;
    std::vector<OverlayVertex> slice_contour_vertices_;
    std::vector<OverlayVertex> slice_plane_vertices_;
    std::vector<OverlayVertex> slice_selected_point_vertices_;
    iso_slices::Result slice_result_data_;
    std::vector<std::pair<std::size_t, std::size_t>> slice_segment_ranges_;
    std::vector<OverlayVertex> grid_vertices_;
    std::vector<OverlayVertex> cell_line_vertices_;
    std::vector<OverlayVertex> cell_point_vertices_;
    std::vector<OverlayVertex> glyph_vertices_;
    std::vector<OverlayVertex> point_of_interest_vertices_;
    std::vector<OverlayVertex> selected_grid_vertex_vertices_;

    void invalidate_slice_overlay()
    {
        // Parameter changes retain the last extracted slice snapshot so it can be compared
        // or re-extracted. The gallery clearly marks it stale until refreshed.
        if (slice_result_available_) slice_data_stale_ = true;
    }

    void clear_slice_overlay()
    {
        slice_contour_vertices_.clear();
        slice_result_data_.segments.clear();
        slice_result_data_.samples.clear();
        slice_plane_vertices_.clear();
        slice_selected_point_vertices_.clear();
        slice_selected_grid_point_valid_ = false;
        slice_reference_geometry_dirty_ = true;
        slice_result_data_.plane_count = 0;
        slice_segment_ranges_.clear();
        slice_plane_count_ = 0;
        slice_sample_count_ = 0;
        slice_segment_count_ = 0;
        slice_result_available_ = false;
        slice_geometry_dirty_ = true;
        slice_data_stale_ = false;
    }

    static pmp::Color gradient_color(float t)
    {
        t = std::clamp(t, 0.0f, 1.0f);
        static constexpr float stops[5][3] = {
            {0.15f, 0.05f, 0.45f},
            {0.10f, 0.45f, 0.75f},
            {0.25f, 0.75f, 0.35f},
            {0.95f, 0.80f, 0.15f},
            {0.90f, 0.15f, 0.15f}};
        const float scaled = t * 4.0f;
        const int i0 = std::min(3, static_cast<int>(scaled));
        const float f = scaled - static_cast<float>(i0);
        return pmp::Color(
            stops[i0][0] + f * (stops[i0 + 1][0] - stops[i0][0]),
            stops[i0][1] + f * (stops[i0 + 1][1] - stops[i0][1]),
            stops[i0][2] + f * (stops[i0 + 1][2] - stops[i0][2]));
    }

    static OverlayVertex overlay_vertex(double x, double y, double z,
                                        float r, float g, float b)
    {
        return {static_cast<float>(x), static_cast<float>(y), static_cast<float>(z), r, g, b};
    }

    static double grid_coord(int sample_index, int resolution, double extent)
    {
        return -extent + (2.0 * extent * static_cast<double>(sample_index) /
                          static_cast<double>(std::max(1, resolution)));
    }

    static void plot_with_cursor(const char* label, const float* ys, int count,
                                 float current_k, float y_min, float y_max,
                                 float height)
    {
        ImGui::PlotLines(label, ys, count, 0, nullptr, y_min, y_max, ImVec2(0, height));
        const ImVec2 min = ImGui::GetItemRectMin();
        const ImVec2 max = ImGui::GetItemRectMax();
        const float t = std::clamp((current_k + 1.0f) * 0.5f, 0.0f, 1.0f);
        const float x = min.x + t * (max.x - min.x);
        ImGui::GetWindowDrawList()->AddLine(
            ImVec2(x, min.y), ImVec2(x, max.y), IM_COL32(255, 255, 255, 180), 1.5f);
    }

    template <std::size_t N>
    static float array_max(const std::array<float, N>& values, float floor)
    {
        float result = floor;
        for (float value : values)
            if (std::isfinite(value)) result = std::max(result, value);
        return result;
    }

    void sanitize_parameters()
    {
        if (!std::isfinite(k_)) k_ = 0.0f;
        k_ = std::clamp(k_, -1.0f, 1.0f);
        grid_resolution_ = std::clamp(grid_resolution_, 1, kMaxGridResolution);
        if (!std::isfinite(half_extent_)) half_extent_ = 1.0f;
        half_extent_ = std::clamp(half_extent_, 0.2f, 3.0f);
        if (!std::isfinite(gradient_gamma_)) gradient_gamma_ = 0.6f;
        gradient_gamma_ = std::clamp(gradient_gamma_, 0.2f, 2.0f);
        grid_stride_ = std::clamp(grid_stride_, 1, grid_resolution_);
        for (int axis = 0; axis < 3; ++axis)
        {
            selected_cell_[axis] = std::clamp(selected_cell_[axis], 0, grid_resolution_ - 1);
            const int clamped_vertex = std::clamp(selected_grid_vertex_[axis], 0, grid_resolution_);
            if (clamped_vertex != selected_grid_vertex_[axis])
            {
                selected_grid_vertex_[axis] = clamped_vertex;
                selected_vertex_geometry_dirty_ = true;
            }
        }
    }

    void apply_expression()
    {
        try
        {
            NodePtr next = parseExpression(expression_buffer_.data());
            ast_ = std::move(next);
            last_good_expression_ = expression_buffer_.data();
            parse_error_.clear();
            runtime_message_.clear();
            invalidate_slice_overlay();
            mesh_dirty_ = true;
            sweep_dirty_ = true;
            glyph_geometry_dirty_ = true;
        }
        catch (const ExprError& error)
        {
            parse_error_ = error.what();
            runtime_message_ = "Invalid expression; keeping the last valid surface.";
        }
    }

    void draw_formula_editor(const char* id_suffix)
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

    int preset_index_ = 0;

    void draw_setup_popup()
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

    void rebuild_surface()
    {
        mesh_dirty_ = false;
        runtime_message_.clear();
        const auto start = std::chrono::steady_clock::now();
        try
        {
            IsosurfaceResult result = Isosurface::build(
                ast_.get(), static_cast<double>(k_), grid_resolution_, half_extent_);
            if (result.vertices.empty())
            {
                runtime_message_ = "No zero-level surface was found in this domain; showing the last non-empty surface.";
                last_build_ms_ = result.buildMillis;
                return;
            }

            // Project's extractor emits triangle soup (three corners per
            // triangle). Reconstruct edge-connected vertex fans before handing
            // it to PMP: smooth normals/colors are shared on regular parts,
            // while disconnected fans at singularities remain separate.
            mesh_.clear();
            auto vertex_colors = mesh_.vertex_property<pmp::Color>(
                "v:color", pmp::Color(0.78f, 0.80f, 0.84f));
            const double denominator = std::log(1.0 + std::max(result.maxGradMag, 1e-9));
            auto color_for = [&](const SurfaceVertex& v)
            {
                float t = static_cast<float>(std::log(1.0 + std::max(0.0f, v.gradMag)) / denominator);
                t = std::pow(std::clamp(t, 0.0f, 1.0f), gradient_gamma_);
                return gradient_color(t);
            };

            using QuantizedPoint = std::array<long long, 3>;
            using GeometricEdge = std::pair<std::size_t, std::size_t>;
            const std::size_t triangle_count = result.vertices.size() / 3;
            const std::size_t corner_count = triangle_count * 3;
            std::map<QuantizedPoint, std::size_t> geometric_ids;
            std::vector<std::size_t> corner_geometric_id(corner_count);
            std::size_t next_geometric_id = 0;
            constexpr double weld_scale = 1000000.0;
            for (std::size_t corner = 0; corner < corner_count; ++corner)
            {
                const SurfaceVertex& v = result.vertices[corner];
                const QuantizedPoint key{
                    std::llround(static_cast<double>(v.px) * weld_scale),
                    std::llround(static_cast<double>(v.py) * weld_scale),
                    std::llround(static_cast<double>(v.pz) * weld_scale)};
                auto insertion = geometric_ids.emplace(key, next_geometric_id);
                if (insertion.second)
                    ++next_geometric_id;
                corner_geometric_id[corner] = insertion.first->second;
            }

            std::map<GeometricEdge, std::vector<std::size_t>> edge_faces;
            auto make_edge = [](std::size_t a, std::size_t b)
            {
                if (a > b) std::swap(a, b);
                return GeometricEdge{a, b};
            };
            for (std::size_t face = 0; face < triangle_count; ++face)
            {
                const std::size_t c = face * 3;
                const std::size_t a = corner_geometric_id[c];
                const std::size_t b = corner_geometric_id[c + 1];
                const std::size_t d = corner_geometric_id[c + 2];
                edge_faces[make_edge(a, b)].push_back(face);
                edge_faces[make_edge(b, d)].push_back(face);
                edge_faces[make_edge(d, a)].push_back(face);
            }

            std::vector<std::size_t> parent(corner_count);
            std::vector<unsigned int> tree_rank(corner_count, 0);
            for (std::size_t i = 0; i < corner_count; ++i) parent[i] = i;
            auto find_root = [&](std::size_t x)
            {
                while (parent[x] != x)
                {
                    parent[x] = parent[parent[x]];
                    x = parent[x];
                }
                return x;
            };
            auto unite = [&](std::size_t a, std::size_t b)
            {
                a = find_root(a);
                b = find_root(b);
                if (a == b) return;
                if (tree_rank[a] < tree_rank[b]) std::swap(a, b);
                parent[b] = a;
                if (tree_rank[a] == tree_rank[b]) ++tree_rank[a];
            };
            auto corner_for = [&](std::size_t face, std::size_t geometric_id)
            {
                const std::size_t first = face * 3;
                for (std::size_t local = 0; local < 3; ++local)
                    if (corner_geometric_id[first + local] == geometric_id)
                        return first + local;
                return first; // an edge endpoint always belongs to its face
            };
            for (const auto& entry : edge_faces)
            {
                const auto& edge = entry.first;
                const auto& faces = entry.second;
                // Only two-face edges define a shared manifold fan. Boundary
                // edges remain open; edges with >2 faces are deliberately split.
                if (faces.size() != 2) continue;
                unite(corner_for(faces[0], edge.first), corner_for(faces[1], edge.first));
                unite(corner_for(faces[0], edge.second), corner_for(faces[1], edge.second));
            }

            std::map<std::size_t, pmp::Vertex> fan_vertices;
            std::vector<pmp::Vertex> corner_vertices(corner_count);
            for (std::size_t face = 0; face < triangle_count; ++face)
            {
                for (std::size_t local = 0; local < 3; ++local)
                {
                    const std::size_t corner = face * 3 + local;
                    const std::size_t root = find_root(corner);
                    auto it = fan_vertices.find(root);
                    if (it == fan_vertices.end())
                    {
                        const SurfaceVertex& source_vertex = result.vertices[corner];
                        const pmp::Vertex pv = mesh_.add_vertex(
                            pmp::Point(source_vertex.px, source_vertex.py, source_vertex.pz));
                        vertex_colors[pv] = color_for(source_vertex);
                        it = fan_vertices.emplace(root, pv).first;
                    }
                    corner_vertices[corner] = it->second;
                }
                const std::size_t c = face * 3;
                mesh_.add_triangle(corner_vertices[c], corner_vertices[c + 1], corner_vertices[c + 2]);
            }
            renderer_.set_use_colors(true);
            update_mesh();
            last_surface_result_ = std::move(result);
            last_triangle_count_ = last_surface_result_.vertices.size() / 3;
            last_build_ms_ = std::chrono::duration<double, std::milli>(
                std::chrono::steady_clock::now() - start).count();
            glyph_geometry_dirty_ = true;
            cell_geometry_dirty_ = true;
            if (!first_mesh_built_)
            {
                set_scene(pmp::vec3(0.0f, 0.0f, 0.0f), std::max(0.5f, 1.8f * half_extent_));
                first_mesh_built_ = true;
            }
        }
        catch (const std::exception& error)
        {
            runtime_message_ = std::string("Surface rebuild failed: ") + error.what();
        }
    }

    void rebuild_k_sweep()
    {
        sweep_dirty_ = false;
        for (int i = 0; i < KSweep::N; ++i)
        {
            const double sweep_k = -1.0 + 2.0 * static_cast<double>(i) / (KSweep::N - 1);
            Dual3 x(point_of_interest_[0], 1, 0, 0);
            Dual3 y(point_of_interest_[1], 0, 1, 0);
            Dual3 z(point_of_interest_[2], 0, 0, 1);
            Dual3 k(sweep_k, 0, 0, 0);
            const Dual3 value = evalExpr<Dual3>(ast_.get(), x, y, z, k);
            sweep_.gx[i] = static_cast<float>(value.dx);
            sweep_.gy[i] = static_cast<float>(value.dy);
            sweep_.gz[i] = static_cast<float>(value.dz);
            sweep_.grad_magnitude[i] = static_cast<float>(
                std::sqrt(value.dx * value.dx + value.dy * value.dy + value.dz * value.dz));
            sweep_.f_value[i] = static_cast<float>(value.v);
        }
    }

    void get_cell_corners(GridCorner corners[8]) const
    {
        static constexpr const char* names[8] = {"C000", "C100", "C110", "C010", "C001", "C101", "C111", "C011"};
        static constexpr int offsets[8][3] = {
            {0, 0, 0}, {1, 0, 0}, {1, 1, 0}, {0, 1, 0},
            {0, 0, 1}, {1, 0, 1}, {1, 1, 1}, {0, 1, 1}};
        for (int i = 0; i < 8; ++i)
        {
            const int ix = selected_cell_[0] + offsets[i][0];
            const int iy = selected_cell_[1] + offsets[i][1];
            const int iz = selected_cell_[2] + offsets[i][2];
            corners[i].name = names[i];
            corners[i].ix = ix;
            corners[i].iy = iy;
            corners[i].iz = iz;
            corners[i].x = grid_coord(ix, grid_resolution_, half_extent_);
            corners[i].y = grid_coord(iy, grid_resolution_, half_extent_);
            corners[i].z = grid_coord(iz, grid_resolution_, half_extent_);
            corners[i].value = evalExpr<double>(ast_.get(), corners[i].x, corners[i].y, corners[i].z, k_);
        }
    }

    std::vector<int> visible_grid_samples() const
    {
        const int stride = std::max(1, grid_stride_);
        std::vector<int> samples;
        for (int i = 0; i <= grid_resolution_; i += stride)
            samples.push_back(i);
        if (samples.empty() || samples.back() != grid_resolution_)
            samples.push_back(grid_resolution_);
        return samples;
    }

    void rebuild_slice_segment_ranges()
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
    void build_slice_reference_geometry()
    {
        slice_plane_vertices_.clear();
        slice_selected_point_vertices_.clear();
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
            return {iso_slices::grid_coord(idx[0], resolution, extent),
                    iso_slices::grid_coord(idx[1], resolution, extent),
                    iso_slices::grid_coord(idx[2], resolution, extent)};
        };
        auto add_line = [&](const iso_slices::Point3& a, const iso_slices::Point3& b) {
            const float r = 0.22f, g = 0.55f, bl = 0.95f;
            slice_plane_vertices_.push_back(overlay_vertex(a.x, a.y, a.z, r, g, bl));
            slice_plane_vertices_.push_back(overlay_vertex(b.x, b.y, b.z, r, g, bl));
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
                slice_selected_point_vertices_.push_back(overlay_vertex(p.x-radius, p.y, p.z, r,g,bl));
                slice_selected_point_vertices_.push_back(overlay_vertex(p.x+radius, p.y, p.z, r,g,bl));
                slice_selected_point_vertices_.push_back(overlay_vertex(p.x, p.y-radius, p.z, r,g,bl));
                slice_selected_point_vertices_.push_back(overlay_vertex(p.x, p.y+radius, p.z, r,g,bl));
                slice_selected_point_vertices_.push_back(overlay_vertex(p.x, p.y, p.z-radius, r,g,bl));
                slice_selected_point_vertices_.push_back(overlay_vertex(p.x, p.y, p.z+radius, r,g,bl));
            }
        }
    }

    // The gallery displays each axis-aligned grid plane side-by-side. Clicking a
    // lattice node selects the exact corresponding 3D point and reports F there.
    void refresh_slice_result_for_current_parameters()
    {
        sanitize_parameters();
        const iso_slices::Result refreshed = iso_slices::extract_and_export(
            ast_.get(), last_good_expression_, static_cast<double>(k_),
            grid_resolution_, static_cast<double>(half_extent_), slice_axis_,
            slice_output_directory_.data());
        slice_result_data_ = refreshed;
        slice_data_stale_ = false;
        slice_contour_vertices_.clear();
        const double denom = static_cast<double>(std::max(1, grid_resolution_));
        for (const iso_slices::ContourSegment& segment : slice_result_data_.segments)
        {
            const float t = static_cast<float>(segment.slice_index / denom);
            const pmp::Color color = gradient_color(t);
            slice_contour_vertices_.push_back(overlay_vertex(
                segment.a.x, segment.a.y, segment.a.z, color[0], color[1], color[2]));
            slice_contour_vertices_.push_back(overlay_vertex(
                segment.b.x, segment.b.y, segment.b.z, color[0], color[1], color[2]));
        }
        slice_plane_count_ = slice_result_data_.plane_count;
        slice_sample_count_ = slice_result_data_.sample_count;
        slice_segment_count_ = slice_result_data_.contour_segment_count;
        slice_result_available_ = true;
        show_slice_contours_ = true;
        slice_geometry_dirty_ = true;
        const int resolution = std::max(1, slice_result_data_.resolution);
        slice_selected_slice_ = std::clamp(slice_selected_slice_, 0, resolution);
        slice_selected_iu_ = std::clamp(slice_selected_iu_, 0, resolution);
        slice_selected_iv_ = std::clamp(slice_selected_iv_, 0, resolution);
        slice_selected_grid_point_valid_ = !slice_result_data_.samples.empty();
        rebuild_slice_segment_ranges();
        slice_reference_geometry_dirty_ = true;
        show_slice_gallery_ = true;
        runtime_message_ = "2D slices refreshed at current k and CSV files updated in '" +
            std::string(slice_output_directory_.data()) + "'.";
    }

    std::array<char, 64> main_k_text_buffer_{};
    bool main_k_text_editing_ = false;
    std::string main_k_text_error_;
    std::array<char, 64> gallery_k_text_buffer_{};
    bool gallery_k_text_editing_ = false;
    std::string gallery_k_text_error_;

    bool draw_k_text_editor(const char* label,
                            std::array<char, 64>& buffer,
                            bool& editing,
                            const char* apply_button_label,
                            std::string& error_message)
    {
        // Do not reformat the text while the user is actively editing it.
        if (!editing)
            std::snprintf(buffer.data(), buffer.size(), "%.9g", static_cast<double>(k_));

        ImGui::SetNextItemWidth(130.0f);
        const bool enter_submitted = ImGui::InputText(
            label, buffer.data(), buffer.size(),
            ImGuiInputTextFlags_EnterReturnsTrue | ImGuiInputTextFlags_CharsScientific);
        const bool field_active = ImGui::IsItemActive();
        const bool deactivated_after_edit = ImGui::IsItemDeactivatedAfterEdit();
        ImGui::SameLine();
        const bool apply_clicked = ImGui::Button(apply_button_label);
        editing = field_active;

        if (enter_submitted || deactivated_after_edit || apply_clicked)
        {
            char* end = nullptr;
            const double parsed = std::strtod(buffer.data(), &end);
            while (end != nullptr && *end != '\0' && std::isspace(static_cast<unsigned char>(*end)))
                ++end;
            const bool valid = end != buffer.data() && end != nullptr && *end == '\0' &&
                std::isfinite(parsed) &&
                parsed <= static_cast<double>(std::numeric_limits<float>::max()) &&
                parsed >= -static_cast<double>(std::numeric_limits<float>::max());

            if (!valid)
            {
                error_message = "Invalid k. Enter a finite number, e.g. -0.001 or 1e-4.";
                std::snprintf(buffer.data(), buffer.size(), "%.9g", static_cast<double>(k_));
                editing = false;
                return false;
            }

            const float new_k = static_cast<float>(parsed);
            const bool changed = new_k != k_;
            k_ = new_k;
            std::snprintf(buffer.data(), buffer.size(), "%.9g", static_cast<double>(k_));
            editing = false;
            error_message.clear();
            return changed;
        }

        if (!error_message.empty())
            ImGui::TextColored(ImVec4(1.0f, 0.35f, 0.35f, 1.0f), "%s", error_message.c_str());
        return false;
    }

    void draw_slice_gallery()
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
        ImGui::SetNextItemWidth(170.0f);
        if (draw_k_text_editor("Viewer k", gallery_k_text_buffer_, gallery_k_text_editing_, "Apply k##gallery", gallery_k_text_error_))
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
            try { refresh_slice_result_for_current_parameters(); }
            catch (const std::exception& error)
            {
                runtime_message_ = std::string("2D slice refresh failed: ") + error.what();
                if (slice_result_available_) slice_data_stale_ = true;
            }
        }
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
                const double plane_coordinate = iso_slices::grid_coord(
                    static_cast<int>(slice_index), resolution, slice_result_data_.half_extent);
                const char* normal_axis = iso_slices::axis_name(slice_result_data_.axis);
                ImGui::Text("Slice %zu | %s = %+.17g | grid %zu/%d",
                    slice_index + 1, normal_axis, plane_coordinate, slice_index, resolution);
                ImGui::PushID(static_cast<int>(slice_index));
                const ImVec2 canvas_origin = ImGui::GetCursorScreenPos();
                ImGui::InvisibleButton("SliceCanvas", ImVec2(canvas_size, canvas_size));
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
        ImGui::End();
    }

    void select_grid_vertex_near(double mouse_x, double mouse_y)
    {
        if (!show_grid_)
            return;

        const std::vector<int> samples = visible_grid_samples();
        if (samples.empty() || width() <= 0 || height() <= 0)
            return;

        const pmp::mat4 model_projection = projection_matrix_ * modelview_matrix_;
        const double pick_radius = 18.0 * std::max(1.0f, high_dpi_scaling());
        double best_distance_squared = pick_radius * pick_radius;
        double best_depth = std::numeric_limits<double>::max();
        std::array<int, 3> best_index{-1, -1, -1};

        // Only sample nodes at line intersections currently shown by Grid line stride.
        // If two nodes overlap in screen space, prefer the one nearer the camera.
        for (int ix : samples)
        {
            const double x = grid_coord(ix, grid_resolution_, half_extent_);
            for (int iy : samples)
            {
                const double y = grid_coord(iy, grid_resolution_, half_extent_);
                for (int iz : samples)
                {
                    const double z = grid_coord(iz, grid_resolution_, half_extent_);
                    const pmp::vec4 clip = model_projection *
                        pmp::vec4(static_cast<float>(x), static_cast<float>(y),
                                  static_cast<float>(z), 1.0f);
                    const double w = static_cast<double>(clip[3]);
                    if (!std::isfinite(w) || w <= 1e-6)
                        continue;

                    const double ndc_x = static_cast<double>(clip[0]) / w;
                    const double ndc_y = static_cast<double>(clip[1]) / w;
                    const double ndc_z = static_cast<double>(clip[2]) / w;
                    if (!std::isfinite(ndc_x) || !std::isfinite(ndc_y) ||
                        !std::isfinite(ndc_z) || ndc_x < -1.0 || ndc_x > 1.0 ||
                        ndc_y < -1.0 || ndc_y > 1.0 || ndc_z < -1.0 || ndc_z > 1.0)
                        continue;

                    const double screen_x = (ndc_x * 0.5 + 0.5) * width();
                    const double screen_y = (0.5 - ndc_y * 0.5) * height();
                    const double dx = screen_x - mouse_x;
                    const double dy = screen_y - mouse_y;
                    const double distance_squared = dx * dx + dy * dy;
                    if (distance_squared < best_distance_squared - 0.25 ||
                        (std::abs(distance_squared - best_distance_squared) <= 0.25 &&
                         ndc_z < best_depth))
                    {
                        best_distance_squared = distance_squared;
                        best_depth = ndc_z;
                        best_index = {ix, iy, iz};
                    }
                }
            }
        }

        if (best_index[0] < 0)
            return;

        selected_grid_vertex_ = best_index;
        has_selected_grid_vertex_ = true;
        for (int axis = 0; axis < 3; ++axis)
            selected_cell_[axis] = std::min(selected_grid_vertex_[axis], grid_resolution_ - 1);

        cell_geometry_dirty_ = true;
        selected_vertex_geometry_dirty_ = true;
        runtime_message_ = "Selected grid vertex [" + std::to_string(best_index[0]) + ", " +
            std::to_string(best_index[1]) + ", " + std::to_string(best_index[2]) + "].";
    }

    void build_grid_geometry()
    {
        grid_vertices_.clear();
        const std::vector<int> samples = visible_grid_samples();
        const float gray_r = 0.24f, gray_g = 0.27f, gray_b = 0.32f;
        auto add_line = [&](double ax, double ay, double az, double bx, double by, double bz)
        {
            grid_vertices_.push_back(overlay_vertex(ax, ay, az, gray_r, gray_g, gray_b));
            grid_vertices_.push_back(overlay_vertex(bx, by, bz, gray_r, gray_g, gray_b));
        };
        for (int ix : samples)
            for (int iy : samples)
                add_line(grid_coord(ix, grid_resolution_, half_extent_),
                         grid_coord(iy, grid_resolution_, half_extent_), -half_extent_,
                         grid_coord(ix, grid_resolution_, half_extent_),
                         grid_coord(iy, grid_resolution_, half_extent_), half_extent_);
        for (int ix : samples)
            for (int iz : samples)
                add_line(grid_coord(ix, grid_resolution_, half_extent_), -half_extent_,
                         grid_coord(iz, grid_resolution_, half_extent_),
                         grid_coord(ix, grid_resolution_, half_extent_), half_extent_,
                         grid_coord(iz, grid_resolution_, half_extent_));
        for (int iy : samples)
            for (int iz : samples)
                add_line(-half_extent_, grid_coord(iy, grid_resolution_, half_extent_),
                         grid_coord(iz, grid_resolution_, half_extent_),
                         half_extent_, grid_coord(iy, grid_resolution_, half_extent_),
                         grid_coord(iz, grid_resolution_, half_extent_));
    }

    void build_cell_geometry()
    {
        cell_line_vertices_.clear();
        cell_point_vertices_.clear();
        GridCorner corners[8];
        get_cell_corners(corners);
        static constexpr int edges[12][2] = {
            {0, 1}, {1, 2}, {2, 3}, {3, 0}, {4, 5}, {5, 6},
            {6, 7}, {7, 4}, {0, 4}, {1, 5}, {2, 6}, {3, 7}};
        for (const auto& edge : edges)
        {
            const auto& a = corners[edge[0]];
            const auto& b = corners[edge[1]];
            cell_line_vertices_.push_back(overlay_vertex(a.x, a.y, a.z, 0.05f, 0.90f, 1.0f));
            cell_line_vertices_.push_back(overlay_vertex(b.x, b.y, b.z, 0.05f, 0.90f, 1.0f));
        }
        for (const auto& corner : corners)
            cell_point_vertices_.push_back(overlay_vertex(corner.x, corner.y, corner.z, 1.0f, 1.0f, 1.0f));
    }

    void build_glyph_geometry()
    {
        glyph_vertices_.clear();
        if (!show_gradient_glyphs_ || last_surface_result_.vertices.empty())
            return;
        const double denominator = std::log(1.0 + std::max(last_surface_result_.maxGradMag, 1e-9));
        const std::size_t stride = static_cast<std::size_t>(std::max(1, glyph_stride_));
        for (std::size_t i = 0; i < last_surface_result_.vertices.size(); i += stride)
        {
            const SurfaceVertex& v = last_surface_result_.vertices[i];
            float t = static_cast<float>(std::log(1.0 + std::max(0.0f, v.gradMag)) / denominator);
            t = std::pow(std::clamp(t, 0.0f, 1.0f), gradient_gamma_);
            const pmp::Color color = gradient_color(t);
            glyph_vertices_.push_back(overlay_vertex(v.px, v.py, v.pz, color[0], color[1], color[2]));
            glyph_vertices_.push_back(overlay_vertex(
                v.px + v.nx * glyph_scale_, v.py + v.ny * glyph_scale_, v.pz + v.nz * glyph_scale_,
                color[0], color[1], color[2]));
        }
    }

    void build_point_of_interest_geometry()
    {
        point_of_interest_vertices_.clear();
        point_of_interest_vertices_.push_back(overlay_vertex(
            point_of_interest_[0], point_of_interest_[1], point_of_interest_[2],
            1.0f, 0.72f, 0.15f));
    }

    void build_selected_grid_vertex_geometry()
    {
        selected_grid_vertex_vertices_.clear();
        if (!has_selected_grid_vertex_)
            return;

        // Center the marker exactly on the selected sample. Use a triangle sphere
        // instead of GL_POINTS: some macOS/Metal OpenGL implementations clamp
        // point sprites to one pixel, making a point-size-based marker invisible.
        const double cx = grid_coord(selected_grid_vertex_[0], grid_resolution_, half_extent_);
        const double cy = grid_coord(selected_grid_vertex_[1], grid_resolution_, half_extent_);
        const double cz = grid_coord(selected_grid_vertex_[2], grid_resolution_, half_extent_);
        const double spacing = (2.0 * static_cast<double>(half_extent_)) /
                               static_cast<double>(std::max(1, grid_resolution_));
        //const double radius = std::max(spacing * 0.18, static_cast<double>(half_extent_) * 0.012);
        const double radius = spacing * 0.06;
        constexpr int rings = 10;
        constexpr int sectors = 16;
        constexpr float red = 1.0f, green = 0.12f, blue = 0.025f;

        auto point = [&](int ring, int sector)
        {
            const double theta = 3.14159265358979323846 * static_cast<double>(ring) / rings;
            const double phi = 2.0 * 3.14159265358979323846 * static_cast<double>(sector) / sectors;
            const double st = std::sin(theta);
            return std::array<double, 3>{
                cx + radius * st * std::cos(phi),
                cy + radius * std::cos(theta),
                cz + radius * st * std::sin(phi)};
        };
        auto add_triangle = [&](const std::array<double, 3>& a,
                                const std::array<double, 3>& b,
                                const std::array<double, 3>& c)
        {
            selected_grid_vertex_vertices_.push_back(overlay_vertex(a[0], a[1], a[2], red, green, blue));
            selected_grid_vertex_vertices_.push_back(overlay_vertex(b[0], b[1], b[2], red, green, blue));
            selected_grid_vertex_vertices_.push_back(overlay_vertex(c[0], c[1], c[2], red, green, blue));
        };

        for (int r = 0; r < rings; ++r)
        {
            for (int sec = 0; sec < sectors; ++sec)
            {
                const auto p00 = point(r, sec);
                const auto p01 = point(r, sec + 1);
                const auto p10 = point(r + 1, sec);
                const auto p11 = point(r + 1, sec + 1);
                if (r != 0)
                    add_triangle(p00, p10, p01);
                if (r != rings - 1)
                    add_triangle(p01, p10, p11);
            }
        }
    }

    static void create_overlay_buffer(OverlayBuffer& buffer)
    {
        glGenVertexArrays(1, &buffer.vao);
        glGenBuffers(1, &buffer.vbo);
        glBindVertexArray(buffer.vao);
        glBindBuffer(GL_ARRAY_BUFFER, buffer.vbo);
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(OverlayVertex),
                              reinterpret_cast<void*>(offsetof(OverlayVertex, x)));
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(OverlayVertex),
                              reinterpret_cast<void*>(offsetof(OverlayVertex, r)));
        glEnableVertexAttribArray(1);
        glBindVertexArray(0);
    }

    static void upload_overlay_buffer(const OverlayBuffer& buffer,
                                      const std::vector<OverlayVertex>& vertices)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buffer.vbo);
        glBufferData(GL_ARRAY_BUFFER,
                     static_cast<GLsizeiptr>(vertices.size() * sizeof(OverlayVertex)),
                     vertices.empty() ? nullptr : vertices.data(), GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
    }

    float supported_line_width(float requested) const
    {
        if (!wide_lines_supported_)
            return 1.0f;
        return std::clamp(requested, line_width_range_[0], line_width_range_[1]);
    }

    float supported_point_size(float requested) const
    {
        return std::clamp(requested, point_size_range_[0], point_size_range_[1]);
    }

    static void delete_overlay_buffer(OverlayBuffer& buffer)
    {
        if (buffer.vbo)
            glDeleteBuffers(1, &buffer.vbo);
        if (buffer.vao)
            glDeleteVertexArrays(1, &buffer.vao);
        buffer.vbo = 0;
        buffer.vao = 0;
    }

    void initialize_overlay()
    {
        glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, line_width_range_.data());
        glGetFloatv(GL_ALIASED_POINT_SIZE_RANGE, point_size_range_.data());
        auto sanitize_range = [](std::array<GLfloat, 2>& range)
        {
            if (!std::isfinite(range[0]) || !std::isfinite(range[1]) || range[0] < 1.0f || range[1] < range[0])
                range = {1.0f, 1.0f};
        };
        // The range query above is not trustworthy on its own: some Core
        // Profile drivers report it generously (we have observed [1, 255]
        // on Mesa) but still raise GL_INVALID_VALUE for any non-1.0 width
        // once actually requested. Clear any pending error, try the widest
        // width this code ever actually requests (2.5, for the selected-
        // cell outline), and check whether the driver truly accepted it
        // before using a non-1.0 width in a real draw call.
        while (glGetError() != GL_NO_ERROR) {} // drain any unrelated pending error first
        const float probe_width = std::clamp(2.5f, line_width_range_[0], line_width_range_[1]);
        glLineWidth(probe_width);
        wide_lines_supported_ = (glGetError() == GL_NO_ERROR);
        glLineWidth(1.0f);
        glGetError(); // discard; a reset to 1.0f cannot itself meaningfully fail
        sanitize_range(line_width_range_);
        sanitize_range(point_size_range_);

        static constexpr const char* vertex_shader = R"GLSL(
            #version 330 core
            layout(location = 0) in vec3 vertex;
            layout(location = 1) in vec3 vertex_color;
            uniform mat4 projection_matrix;
            uniform mat4 modelview_matrix;
            uniform float point_size;
            out vec3 color;
            void main()
            {
                gl_Position = projection_matrix * modelview_matrix * vec4(vertex, 1.0);
                gl_PointSize = point_size;
                color = vertex_color;
            }
        )GLSL";
        static constexpr const char* fragment_shader = R"GLSL(
            #version 330 core
            in vec3 color;
            uniform float alpha;
            out vec4 frag_color;
            void main() { frag_color = vec4(color, alpha); }
        )GLSL";
        overlay_shader_.source(vertex_shader, fragment_shader);
        create_overlay_buffer(grid_buffer_);
        create_overlay_buffer(cell_line_buffer_);
        create_overlay_buffer(cell_point_buffer_);
        create_overlay_buffer(glyph_buffer_);
        create_overlay_buffer(point_of_interest_buffer_);
        create_overlay_buffer(selected_grid_vertex_buffer_);
        create_overlay_buffer(slice_buffer_);
        create_overlay_buffer(slice_plane_buffer_);
        create_overlay_buffer(slice_selected_point_buffer_);
        overlay_initialized_ = true;
        grid_geometry_dirty_ = true;
        cell_geometry_dirty_ = true;
        glyph_geometry_dirty_ = true;
        point_of_interest_geometry_dirty_ = true;
        selected_vertex_geometry_dirty_ = true;
    }
};

} // namespace iso

int main()
{
    try
    {
        iso::MyViewer viewer("Isosurface Singularity Lab", 1280, 800);
        return viewer.run();
    }
    catch (const std::exception& error)
    {
        std::cerr << "Isosurface Singularity Lab failed: " << error.what() << '\n';
        return 1;
    }
}
