#pragma once

// The interactive viewer: a pmp::MeshViewer showing the zero level set of
// F(x,y,z,k) with ImGui controls, overlays and analysis tools.
//
// Implementation is split by concern:
//   IsosurfaceViewer.cpp  frame loop, surface extraction, picking, 3D overlays
//   ViewerPanels.cpp      main ImGui control panel and the setup dialog
//   SliceGallery.cpp      2D slice extraction state and the slice gallery
//   AnalysisJobs.cpp      background convergence analysis and the F+kG batch

#include "BackgroundJob.h"
#include "OverlayRenderer.h"

#include "analysis/SingularityBatch.h"
#include "analysis/SingularityConvergence.h"
#include "core/Expr.h"
#include "core/Isosurface.h"
#include "core/IsosurfaceSlices.h"

#include <pmp/visualization/mesh_viewer.h>

#include <array>
#include <atomic>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace iso_app {

class IsosurfaceViewer : public pmp::MeshViewer
{
public:
    IsosurfaceViewer(const char* title, int width, int height);
    ~IsosurfaceViewer() override;

protected:
    void do_processing() override;
    void process_imgui() override;
    void draw(const std::string& draw_mode) override;
    void mouse(int button, int action, int mods) override;

private:
    struct GridCorner
    {
        const char* name = "";
        int ix = 0, iy = 0, iz = 0;
        double x = 0.0, y = 0.0, z = 0.0, value = 0.0;
    };

    // F and grad F at the point of interest for k in [-1, 1].
    struct KSweep
    {
        static constexpr int N = 201;
        std::array<float, N> grad_magnitude{};
        std::array<float, N> gx{};
        std::array<float, N> gy{};
        std::array<float, N> gz{};
        std::array<float, N> f_value{};
    };

    struct ConvergenceJobResult
    {
        iso_convergence::Analysis analysis;
        std::string message;
    };

    static constexpr int kMaxGridResolution = 110;
    static constexpr int kActiveFramesAfterEvent = 10;

    // ---- IsosurfaceViewer.cpp: frame loop, surface, picking, overlays -------
    void throttle_when_idle();
    void sanitize_parameters();
    void apply_expression();
    void rebuild_surface();
    void recolor_surface();
    pmp::Color surface_color(float grad_mag) const;
    void rebuild_k_sweep();
    void get_cell_corners(GridCorner corners[8]) const;
    std::vector<int> visible_grid_samples() const;
    void select_grid_vertex_near(double mouse_x, double mouse_y);
    void refresh_overlays(); // rebuilds every dirty overlay layer
    std::vector<OverlayVertex> build_grid_geometry() const;
    void build_cell_geometry(std::vector<OverlayVertex>& edge_lines, std::vector<OverlayVertex>& corner_points) const;
    std::vector<OverlayVertex> build_glyph_geometry() const;
    std::vector<OverlayVertex> build_point_of_interest_geometry() const;
    std::vector<OverlayVertex> build_selected_grid_vertex_geometry() const;

    // ---- ViewerPanels.cpp: ImGui control panel ------------------------------
    void draw_function_panel();
    void draw_sampling_panel();
    void draw_slice_panel();
    void draw_cell_panel();
    void draw_gradient_panel();
    void draw_rendering_panel();
    void draw_formula_editor(const char* id_suffix);
    bool draw_k_text_editor(const char* label, std::string& error_message);
    void draw_setup_popup();
    static void plot_with_cursor(const char* label, const float* ys, int count, float current_k,
                                 float y_min, float y_max, float height);

    // ---- SliceGallery.cpp: 2D slices ----------------------------------------
    iso_slices::Result extract_slices_at_current_parameters();
    void apply_slice_result(iso_slices::Result&& result);
    void invalidate_slice_overlay();
    void clear_slice_overlay();
    void rebuild_slice_segment_ranges();
    void build_slice_reference_geometry(std::vector<OverlayVertex>& plane,
                                        std::vector<OverlayVertex>& selected_point) const;
    void draw_slice_gallery();

    // ---- AnalysisJobs.cpp: background analyses ------------------------------
    static std::filesystem::path experiments_directory();
    void start_convergence_job();
    static ConvergenceJobResult run_convergence_job(const iso_convergence::Options& options,
                                                    const std::filesystem::path& experiment_root);
    void poll_convergence_job();
    void draw_singularity_convergence_report();
    void start_batch_experiment();
    void poll_batch_job();
    void draw_batch_experiment_popup();

    // ---- Expression -----------------------------------------------------------
    std::array<char, 512> expression_buffer_{};
    std::string last_good_expression_ = "x*x*x - k*y";
    std::string parse_error_;
    std::string runtime_message_;
    NodePtr ast_;
    int preset_index_ = 0;
    std::string main_k_error_;
    std::string gallery_k_error_;

    // ---- Sampling and surface -------------------------------------------------
    float k_ = 0.0f;
    int grid_resolution_ = 48; // number of cells; there are resolution+1 samples
    float half_extent_ = 1.0f;
    float gradient_gamma_ = 0.6f;
    bool mesh_dirty_ = true;
    bool first_mesh_built_ = false;
    double last_build_ms_ = -1.0;
    std::size_t last_triangle_count_ = 0;
    IsosurfaceResult last_surface_result_;
    // |grad F| per PMP vertex of the displayed surface and the log-scale color
    // normalizer, so a gamma change can recolor without re-extracting.
    std::vector<float> surface_vertex_grad_mag_;
    double surface_color_denominator_ = 1.0;
    std::array<float, 3> point_of_interest_{0.0f, 0.0f, 0.0f};
    KSweep sweep_;
    bool sweep_dirty_ = true;

    // ---- Dialogs and idle throttling -------------------------------------------
    bool setup_open_ = false;
    bool setup_popup_requested_ = false;
    int active_frames_left_ = kActiveFramesAfterEvent;

    // ---- Selection and 3D overlays ---------------------------------------------
    OverlayRenderer overlay_;
    int selected_cell_[3] = {24, 24, 24};
    std::array<int, 3> selected_grid_vertex_{0, 0, 0};
    bool has_selected_grid_vertex_ = false;
    std::array<double, 2> selection_click_start_{0.0, 0.0};
    bool selection_click_candidate_ = false;
    int grid_stride_ = 4;
    float glyph_scale_ = 0.15f;
    int glyph_stride_ = 23;
    bool show_grid_ = false;
    bool show_selected_cell_ = true;
    bool show_gradient_glyphs_ = false;
    bool grid_geometry_dirty_ = true;
    bool cell_geometry_dirty_ = true;
    bool glyph_geometry_dirty_ = true;
    bool point_of_interest_geometry_dirty_ = true;
    bool selected_vertex_geometry_dirty_ = true;

    // ---- 2D slices -------------------------------------------------------------
    iso_slices::Result slice_result_data_;
    std::vector<std::pair<std::size_t, std::size_t>> slice_segment_ranges_; // [begin, end) per slice
    std::array<char, 512> slice_output_directory_{};
    int slice_axis_ = 2;
    bool slice_result_available_ = false;
    bool slice_data_stale_ = false;
    bool slice_reference_geometry_dirty_ = true;
    int slice_selected_slice_ = 0;
    int slice_selected_iu_ = 0;
    int slice_selected_iv_ = 0;
    bool slice_selected_grid_point_valid_ = false;
    bool show_slice_contours_ = true;
    bool show_slice_gallery_ = false;
    bool show_slice_grid_values_ = true;
    bool show_slice_full_precision_labels_ = false;
    bool show_slice_grid_in_3d_ = true;
    bool show_slice_selected_point_in_3d_ = true;
    int slice_value_label_stride_ = 1;
    int slice_label_resolution_ = -1;

    // ---- Convergence analysis --------------------------------------------------
    double convergence_k_min_ = -0.1;
    double convergence_k_max_ = 0.1;
    int convergence_sample_count_ = 25;
    int convergence_resolution_ = 48;
    bool convergence_analysis_available_ = false;
    iso_convergence::Analysis convergence_analysis_{};
    std::string convergence_message_;
    BackgroundJob<ConvergenceJobResult> convergence_job_;

    // ---- F + kG batch ----------------------------------------------------------
    std::array<char, 512> batch_f_expression_buffer_{};
    std::array<char, 16384> batch_g_expressions_buffer_{};
    double batch_k_min_ = -0.1;
    double batch_k_max_ = 0.1;
    int batch_k_sample_count_ = 25;
    int batch_mesh_resolution_ = 48;
    double batch_domain_half_extent_ = 1.0;
    bool batch_setup_open_ = false;
    bool batch_popup_requested_ = false;
    std::string batch_status_;
    std::shared_ptr<std::atomic<int>> batch_progress_;
    BackgroundJob<iso_batch::BatchResult> batch_job_;
};

} // namespace iso_app
