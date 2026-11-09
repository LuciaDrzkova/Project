// Long-running analyses executed in the background: the convergence +
// research analysis of the current family, its k -> 0 convergence
// classification, and the 25-run F + kG batch.

#include "IsosurfaceViewer.h"

#include "ProjectPaths.h"
#include "analysis/ConvergenceClassifier.h"
#include "core/CsvFormat.h"

#include <imgui.h>

#include <algorithm>
#include <chrono>
#include <cmath>
#include <stdexcept>
#include <string>
#include <utility>

namespace iso_app {

// Generated output lives next to the sources (experiments/, snapshots/),
// independent of the working directory the viewer was started from.
std::filesystem::path IsosurfaceViewer::project_output_directory(const char* name)
{
#ifdef ISOSURFACE_PROJECT_ROOT
    return std::filesystem::path(ISOSURFACE_PROJECT_ROOT) / name;
#else
    return std::filesystem::current_path() / name;
#endif
}

// Convergence settings shared by the analysis buttons: P, the viewer domain,
// the selected k samples and the convergence mesh resolution.
iso_convergence::Options IsosurfaceViewer::convergence_options_from_ui(bool log_spaced_k) const
{
    iso_convergence::Options options;
    options.probe = {static_cast<double>(point_of_interest_[0]),
                     static_cast<double>(point_of_interest_[1]),
                     static_cast<double>(point_of_interest_[2])};
    options.resolution = std::clamp(convergence_resolution_, 16, kMaxGridResolution);
    options.halfExtent = static_cast<double>(half_extent_);
    options.localHalfSize = 0.45 * options.halfExtent;
    options.kValues = log_spaced_k
        ? iso_classify::makeSymmetricLogKValues(std::max(std::fabs(convergence_k_min_), std::fabs(convergence_k_max_)),
                                                std::clamp((convergence_sample_count_ + 1) / 2, 3, 50), 3.0)
        : iso_convergence::makeLinearKValues(convergence_k_min_, convergence_k_max_, convergence_sample_count_);
    options.familyExpression = last_good_expression_;
    return options;
}

// <root>/<prefix>_<microseconds>[_n]: a fresh directory per run.
std::filesystem::path IsosurfaceViewer::create_run_directory(const std::filesystem::path& root, const std::string& prefix)
{
    std::filesystem::create_directories(root);
    const auto stamp = std::chrono::duration_cast<std::chrono::microseconds>(
        std::chrono::system_clock::now().time_since_epoch()).count();
    std::filesystem::path directory = root / (prefix + "_" + std::to_string(stamp));
    for (std::size_t suffix = 1; std::filesystem::exists(directory); ++suffix)
        directory = root / (prefix + "_" + std::to_string(stamp) + "_" + std::to_string(suffix));
    std::filesystem::create_directories(directory);
    return directory;
}

// All CSV reports of one analysis (see iso_research::writeFamilyCsvs).
void IsosurfaceViewer::write_family_csvs(const std::filesystem::path& directory, const iso_research::FamilyAnalysis& family)
{
    std::string write_error;
    if (!iso_research::writeFamilyCsvs(directory.string(), family, &write_error))
        throw std::runtime_error(write_error);
}

// Runs the convergence + research analysis in the background. The worker
// parses its own copy of the expression and receives every setting by
// value, so the viewer stays responsive and may keep editing meanwhile.
void IsosurfaceViewer::start_convergence_job()
{
    if (convergence_job_.running())
        return;
    try
    {
        iso_convergence::Options options = convergence_options_from_ui(false);
        convergence_message_ = "Running convergence analysis (" +
            std::to_string(options.kValues.size()) + " k samples at resolution " +
            std::to_string(options.resolution) + ")...";
        convergence_job_.start([options = std::move(options), root = project_output_directory("experiments")]() {
            return run_convergence_job(options, root);
        });
    }
    catch (const std::exception& error)
    {
        convergence_message_ = std::string("Convergence analysis failed: ") + error.what();
    }
}

IsosurfaceViewer::ConvergenceJobResult IsosurfaceViewer::run_convergence_job(
    const iso_convergence::Options& options, const std::filesystem::path& experiment_root)
{
    const NodePtr ast = parseExpression(options.familyExpression);
    iso_research::FamilyAnalysis family = iso_research::analyzeFamily(ast.get(), options, 3);
    const std::filesystem::path run_directory = create_run_directory(experiment_root, "run");
    write_family_csvs(run_directory, family);

    ConvergenceJobResult result;
    result.message = "Saved four CSV files in " + run_directory.string() + ".\n" + family.research.summary +
        "\nThe Newton search and triangle-sampled distances are numerical diagnostics, not proofs.";
    result.analysis = std::move(family.convergence);
    return result;
}

void IsosurfaceViewer::poll_convergence_job()
{
    convergence_job_.poll(
        [this](ConvergenceJobResult result) {
            convergence_analysis_ = std::move(result.analysis);
            convergence_analysis_available_ = true;
            convergence_message_ = std::move(result.message);
        },
        [this](const std::string& error) {
            convergence_analysis_available_ = false;
            convergence_message_ = "Convergence analysis failed: " + error;
        });
}

void IsosurfaceViewer::draw_singularity_convergence_report()
{
    if (!convergence_analysis_available_)
    {
        if (!convergence_message_.empty())
            ImGui::TextWrapped("%s", convergence_message_.c_str());
        return;
    }

    const auto& a = convergence_analysis_;
    ImGui::Text("Samples: %zu; range [% .6g, % .6g]; grid spacing %.5g",
                a.samples.size(), a.metadata.kMinimum, a.metadata.kMaximum, a.gridSpacing);
    ImGui::TextWrapped("At P, k=0: F=% .4g and |grad F|=% .4g. These should be near zero if P is the singular point being investigated; this check does not prove isolation.",
                       a.probeAtZero.value, a.probeAtZero.gradientNorm);
    ImGui::TextWrapped("k < 0: point-distance trend: %s (log-log slope %s). Local mesh-distance trend: %s.",
                       a.negativeSide.pointDistance.classification.c_str(),
                       iso_csv::shortNumber(a.negativeSide.pointDistance.logLogSlope).c_str(),
                       a.negativeSide.localHausdorff.classification.c_str());
    ImGui::TextWrapped("k > 0: point-distance trend: %s (log-log slope %s). Local mesh-distance trend: %s.",
                       a.positiveSide.pointDistance.classification.c_str(),
                       iso_csv::shortNumber(a.positiveSide.pointDistance.logLogSlope).c_str(),
                       a.positiveSide.localHausdorff.classification.c_str());
    ImGui::TextDisabled("Trend labels are heuristics over the selected range. Increase mesh resolution and inspect the CSV before interpreting them.");

    if (ImGui::CollapsingHeader("Convergence sample results", false))
    {
        if (ImGui::BeginTable("ConvergenceResults", 7,
                              ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                              ImGuiTableFlags_SizingStretchProp))
        {
            ImGui::TableSetupColumn("k");
            ImGui::TableSetupColumn("Vertices");
            ImGui::TableSetupColumn("Local vertices");
            ImGui::TableSetupColumn("d(P, mesh)");
            ImGui::TableSetupColumn("min |grad F|");
            ImGui::TableSetupColumn("max |F| local");
            ImGui::TableSetupColumn("legacy vertex d_H");
            ImGui::TableHeadersRow();
            for (const auto& sample : a.samples)
            {
                ImGui::TableNextRow();
                ImGui::TableNextColumn(); ImGui::Text("% .6g", sample.k);
                ImGui::TableNextColumn(); ImGui::Text("%zu", sample.totalVertices);
                ImGui::TableNextColumn(); ImGui::Text("%zu", sample.localVertices);
                ImGui::TableNextColumn(); ImGui::Text("%s", iso_csv::shortNumber(sample.distanceProbeToMesh).c_str());
                ImGui::TableNextColumn(); ImGui::Text("%s", iso_csv::shortNumber(sample.minimumGradientLocal).c_str());
                ImGui::TableNextColumn(); ImGui::Text("%s", iso_csv::shortNumber(sample.maxAbsFieldResidualLocal).c_str());
                ImGui::TableNextColumn(); ImGui::Text("%s", iso_csv::shortNumber(sample.symmetricHausdorff).c_str());
            }
            ImGui::EndTable();
        }
    }
    if (!convergence_message_.empty())
        ImGui::TextWrapped("%s", convergence_message_.c_str());
}

// ----- k -> 0 convergence classification ---------------------------------------

void IsosurfaceViewer::start_classification_job()
{
    if (classification_job_.running())
        return;
    try
    {
        iso_convergence::Options options = convergence_options_from_ui(classification_log_spaced_k_);
        classification_message_ = "Classifying convergence (" + std::to_string(options.kValues.size()) +
            " k samples at resolution " + std::to_string(options.resolution) + ")...";
        classification_job_.start([options = std::move(options), root = project_output_directory("experiments")]() {
            return run_classification_job(options, root);
        });
    }
    catch (const std::exception& error)
    {
        classification_message_ = std::string("Classification failed: ") + error.what();
    }
}

IsosurfaceViewer::ClassificationJobResult IsosurfaceViewer::run_classification_job(
    const iso_convergence::Options& options, const std::filesystem::path& experiment_root)
{
    const NodePtr ast = parseExpression(options.familyExpression);
    const iso_research::FamilyAnalysis family = iso_research::analyzeFamily(ast.get(), options, 3);
    ClassificationJobResult result;
    result.classification = iso_classify::classify(ast.get(), family);

    const std::filesystem::path run_directory = create_run_directory(experiment_root, "classify");
    write_family_csvs(run_directory, family);
    std::string write_error;
    if (!iso_classify::writeCsv((run_directory / "convergence_classification.csv").string(),
                                result.classification, &write_error))
        throw std::runtime_error(write_error);
    result.message = "Saved convergence_classification.csv and the four analysis CSVs in " +
        run_directory.string() + ".";
    return result;
}

void IsosurfaceViewer::poll_classification_job()
{
    classification_job_.poll(
        [this](ClassificationJobResult result) {
            classification_ = std::move(result.classification);
            classification_available_ = true;
            classification_message_ = std::move(result.message);
        },
        [this](const std::string& error) {
            classification_available_ = false;
            classification_message_ = "Classification failed: " + error;
        });
}

void IsosurfaceViewer::draw_classification_report()
{
    poll_classification_job();
    if (classification_job_.running() || !classification_available_)
    {
        if (!classification_message_.empty())
            ImGui::TextWrapped("%s", classification_message_.c_str());
        return;
    }
    const iso_classify::Classification& c = classification_;
    if (!ImGui::CollapsingHeader("Convergence classification (k -> 0)", ImGuiTreeNodeFlags_DefaultOpen))
        return;
    ImGui::TextWrapped("%s; prediction %s: d_H ~ %s", iso_classify::toString(c.point), c.prediction.regime.c_str(),
                       iso_classify::detail::predictionText(c.prediction).c_str());
    if (!c.prediction.principalPart.empty())
        ImGui::TextWrapped("Principal part of F at P: %s, weights (%.3g, %.3g, %.3g)", c.prediction.principalPart.c_str(),
                           c.prediction.weights[0], c.prediction.weights[1], c.prediction.weights[2]);
    if (c.milnorNumber >= 0)
        ImGui::TextWrapped("Milnor number %d, Tjurina number %d (basis %s); essential part of G: %s",
                           c.milnorNumber, c.tjurinaNumber, c.tjurinaBasis.c_str(), c.essentialCoordinates.c_str());
    ImGui::TextDisabled("Local box P + [-a, a]^3 with a = %.3g; grid spacing h = %.3g. %s", c.localHalfSize,
                        c.gridSpacing,
                        c.negative.exactDistances ? "Distances measured on the exact surfaces."
                                                  : "Distances measured between meshes, fitted where h/2 <= d_H <= 0.6 a.");
    if (ImGui::BeginTable("ClassificationTable", 5, ImGuiTableFlags_Borders | ImGuiTableFlags_RowBg |
                                                       ImGuiTableFlags_SizingStretchProp))
    {
        ImGui::TableSetupColumn("side");
        ImGui::TableSetupColumn("surface");
        ImGui::TableSetupColumn("alpha");
        ImGui::TableSetupColumn("singular points");
        ImGui::TableSetupColumn("count, rate");
        ImGui::TableHeadersRow();
        for (const iso_classify::SideResult* s : {&c.negative, &c.positive})
        {
            ImGui::TableNextRow();
            ImGui::TableNextColumn(); ImGui::TextUnformatted(s->sign < 0 ? "k < 0" : "k > 0");
            ImGui::TableNextColumn(); ImGui::TextUnformatted(iso_classify::toString(s->surface));
            ImGui::TableNextColumn();
            ImGui::Text("%s +/- %s", iso_csv::shortNumber(s->hausdorff.exponent).c_str(),
                        iso_csv::shortNumber(s->hausdorff.ciHigh - s->hausdorff.exponent).c_str());
            ImGui::TableNextColumn(); ImGui::TextUnformatted(iso_classify::toString(s->singular));
            ImGui::TableNextColumn();
            ImGui::Text("%zu, %s", s->maxSingularPoints, iso_csv::shortNumber(s->singularDistance.exponent).c_str());
        }
        ImGui::EndTable();
    }
    ImGui::TextWrapped("%s", c.negative.verdict.c_str());
    ImGui::TextWrapped("%s", c.positive.verdict.c_str());
    if (!classification_message_.empty())
        ImGui::TextWrapped("%s", classification_message_.c_str());
}

void IsosurfaceViewer::start_batch_experiment()
{
    if (batch_job_.running())
        return;
    try
    {
        iso_batch::Settings settings;
        settings.fExpression = iso_batch::trim(batch_f_expression_buffer_.data());
        settings.gExpressions = iso_batch::parseGExpressions(batch_g_expressions_buffer_.data());
        settings.kMinimum = batch_k_min_;
        settings.kMaximum = batch_k_max_;
        settings.kSampleCount = batch_k_sample_count_;
        settings.resolution = batch_mesh_resolution_;
        settings.halfExtent = batch_domain_half_extent_;
        settings.probe = {static_cast<double>(point_of_interest_[0]),
                          static_cast<double>(point_of_interest_[1]),
                          static_cast<double>(point_of_interest_[2])};
        iso_batch::validateSettings(settings);

        batch_progress_ = std::make_shared<std::atomic<int>>(0);
        const auto progress = batch_progress_;
        const auto outputRoot = project_output_directory("experiments");
        batch_status_ = "Running 25 F + kG analyses...";
        batch_job_.start([settings = std::move(settings), outputRoot, progress]() {
            return iso_batch::runBatch(settings, outputRoot, progress.get());
        });
    }
    catch (const std::exception& error)
    {
        batch_status_ = std::string("Batch setup error: ") + error.what();
    }
}

void IsosurfaceViewer::poll_batch_job()
{
    batch_job_.poll(
        [this](const iso_batch::BatchResult& result) { batch_status_ = iso_batch::summarize(result); },
        [this](const std::string& error) { batch_status_ = "Batch execution failed: " + error; });
}

void IsosurfaceViewer::draw_batch_experiment_popup()
{
    const char* popupName = "Batch F + kG experiments";
    if (batch_setup_open_ && !batch_popup_requested_)
    {
        ImGui::OpenPopup(popupName);
        batch_popup_requested_ = true;
    }
    if (!batch_setup_open_)
        batch_popup_requested_ = false;

    ImGui::SetNextWindowSize(ImVec2(650.0f, 760.0f), ImGuiCond_FirstUseEver);
    if (!ImGui::BeginPopupModal(popupName, &batch_setup_open_, ImGuiWindowFlags_None))
        return;

    ImGui::TextWrapped("Enter F(x,y,z) and exactly 25 perturbations G_i(x,y,z). "
                       "Each run analyses F + k G_i over the same k range and grid.");
    ImGui::Separator();
    ImGui::SetNextItemWidth(-1.0f);
    ImGui::InputText("F(x,y,z)", batch_f_expression_buffer_.data(),
                     batch_f_expression_buffer_.size());
    ImGui::TextDisabled("Do not put k in F or G. The program constructs F + k*G automatically.");
    ImGui::Text("G expressions: one expression per line (exactly 25 non-empty lines)");
    ImGui::InputTextMultiline("##BatchGExpressions", batch_g_expressions_buffer_.data(),
                              batch_g_expressions_buffer_.size(), ImVec2(-1.0f, 235.0f),
                              ImGuiInputTextFlags_AllowTabInput);

    ImGui::SetNextItemWidth(145.0f);
    ImGui::InputDouble("k minimum", &batch_k_min_, 0.0, 0.0, "%.10g");
    ImGui::SameLine();
    ImGui::SetNextItemWidth(145.0f);
    ImGui::InputDouble("k maximum", &batch_k_max_, 0.0, 0.0, "%.10g");
    ImGui::SetNextItemWidth(145.0f);
    ImGui::InputInt("Number of k samples", &batch_k_sample_count_);
    batch_k_sample_count_ = std::clamp(batch_k_sample_count_, 3, 101);
    ImGui::SetNextItemWidth(145.0f);
    ImGui::InputInt("Mesh resolution (cells/axis)", &batch_mesh_resolution_);
    batch_mesh_resolution_ = std::clamp(batch_mesh_resolution_, 16, Isosurface::kMaxAnalysisResolution);
    ImGui::SetNextItemWidth(145.0f);
    ImGui::InputDouble("Domain half-extent", &batch_domain_half_extent_, 0.0, 0.0, "%.8g");
    ImGui::Text("P used for all 25 runs: (% .6g, % .6g, % .6g)",
                point_of_interest_[0], point_of_interest_[1], point_of_interest_[2]);
    ImGui::TextDisabled("Resolution: 16-%d cells per axis; k samples: 3-101.", Isosurface::kMaxAnalysisResolution);
    ImGui::Separator();

    if (batch_job_.running())
    {
        const int completed = batch_progress_ ? batch_progress_->load(std::memory_order_relaxed) : 0;
        ImGui::Text("Progress: %d / 25", completed);
        ImGui::ProgressBar(static_cast<float>(completed) / 25.0f, ImVec2(-1.0f, 0.0f));
        ImGui::TextWrapped("The 25 runs are executing in the background. The viewer remains responsive.");
    }
    else
    {
        if (ImGui::Button("Run all 25 tests", ImVec2(180.0f, 0.0f)))
            start_batch_experiment();
        ImGui::SameLine();
        if (ImGui::Button("Close", ImVec2(100.0f, 0.0f)))
        {
            batch_setup_open_ = false;
            ImGui::CloseCurrentPopup();
        }
    }

    if (!batch_status_.empty())
    {
        ImGui::Separator();
        ImGui::TextWrapped("%s", batch_status_.c_str());
    }
    ImGui::EndPopup();
}

} // namespace iso_app
