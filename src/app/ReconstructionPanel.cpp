// Reconstruction of the singular surface S_0 = {F = 0} from the perturbed
// surfaces S_-eps and S_+eps of the current family H_k = F + kG, with any of
// the methods of analysis/Reconstruction.h. The result replaces the displayed
// surface until F, k or the grid change; S_-eps and S_+eps are drawn as
// wireframes around P.

#include "IsosurfaceViewer.h"

#include <imgui.h>

#include <cmath>
#include <cstdio>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>

namespace iso_app {

namespace {

// Wireframe of the mesh edges with both ends in P + [-box, box]^3.
std::vector<OverlayVertex> wireframe(const iso_reconstruct::Mesh& m, const iso_geometry::Point3& P, double box,
                                     float r, float g, float b)
{
    std::set<std::pair<std::uint32_t, std::uint32_t>> edges;
    for (const auto& t : m.t)
        for (int c = 0; c < 3; ++c)
        {
            const std::uint32_t a = t[static_cast<std::size_t>(c)], e = t[static_cast<std::size_t>((c + 1) % 3)];
            if (iso_reconstruct::infDistance(m.v[a], P) <= box && iso_reconstruct::infDistance(m.v[e], P) <= box)
                edges.insert({std::min(a, e), std::max(a, e)});
        }
    std::vector<OverlayVertex> lines;
    lines.reserve(edges.size() * 2);
    for (const auto& [a, e] : edges)
    {
        lines.push_back(overlay_vertex(m.v[a].x, m.v[a].y, m.v[a].z, r, g, b));
        lines.push_back(overlay_vertex(m.v[e].x, m.v[e].y, m.v[e].z, r, g, b));
    }
    return lines;
}

// Largest distance from a used mesh vertex within `box` of P to {F = 0}
// (Newton projection along grad F).
double vertex_error(const iso_reconstruct::Mesh& m, const iso_exact::Field& f, const iso_geometry::Point3& P,
                    double box, double h)
{
    std::vector<char> used(m.v.size(), 0);
    for (const auto& t : m.t) used[t[0]] = used[t[1]] = used[t[2]] = 1;
    double worst = 0.0;
    for (std::size_t i = 0; i < m.v.size(); ++i)
    {
        if (!used[i] || iso_reconstruct::infDistance(m.v[i], P) > box) continue;
        iso_geometry::Point3 p = m.v[i];
        if (iso_exact::detail::project(f, p, h)) worst = std::max(worst, iso_geometry::distance(p, m.v[i]));
    }
    return worst;
}

std::string topology_text(const iso_reconstruct::Topology& t)
{
    char buffer[160];
    std::snprintf(buffer, sizeof(buffer), "%zu component(s), Euler char. %lld, %zu hole edge(s), %zu non-manifold edge(s) -> %s",
                  t.components, t.euler, t.holeEdges, t.nonManifoldEdges,
                  t.contractible() ? "correct (a cone)" : "NOT the topology of a cone");
    return buffer;
}

bool uses_patch(iso_reconstruct::Method m)
{
    return m == iso_reconstruct::Method::SingularCone || m == iso_reconstruct::Method::HybridCone;
}

} // namespace

void IsosurfaceViewer::draw_reconstruction_panel()
{
    if (!ImGui::CollapsingHeader("Reconstruct S0 from S-eps and S+eps", false))
        return;
    ImGui::TextWrapped("Builds S-eps and S+eps of the current family F + kG (k is ignored) and reconstructs "
                       "S0 = {F = 0} near the singular point P.");
    ImGui::SetNextItemWidth(210.0f);
    if (ImGui::InputFloat3("P (singular point)", point_of_interest_.data(), "%.3f"))
    {
        sweep_dirty_ = true;
        point_of_interest_geometry_dirty_ = true;
    }
    ImGui::SetNextItemWidth(120.0f);
    if (ImGui::InputFloat("eps", &reconstruction_eps_, 0.0f, 0.0f, "%.5f"))
        reconstruction_eps_ = std::clamp(std::isfinite(reconstruction_eps_) ? std::fabs(reconstruction_eps_) : 0.01f,
                                         1e-6f, 0.5f);

    const int count = static_cast<int>(iso_reconstruct::Method::Count);
    const char* current = iso_reconstruct::methodLabel(static_cast<iso_reconstruct::Method>(reconstruction_method_));
    ImGui::SetNextItemWidth(-1.0f);
    if (ImGui::BeginCombo("##ReconstructionMethod", current))
    {
        for (int m = 0; m < count; ++m)
        {
            const auto method = static_cast<iso_reconstruct::Method>(m);
            if (ImGui::Selectable(iso_reconstruct::methodLabel(method), m == reconstruction_method_))
                reconstruction_method_ = m;
        }
        ImGui::EndCombo();
    }
    const auto method = static_cast<iso_reconstruct::Method>(reconstruction_method_);
    if (uses_patch(method))
    {
        ImGui::SetNextItemWidth(160.0f);
        ImGui::SliderFloat("Min. patch radius (cells)", &reconstruction_patch_cells_, 2.0f, 12.0f, "%.1f");
    }

    const bool running = reconstruction_job_.running();
    if (running) ImGui::BeginDisabled();
    if (ImGui::Button("Reconstruct S0"))
        start_reconstruction_job();
    if (running) ImGui::EndDisabled();
    if (reconstruction_shown_)
    {
        ImGui::SameLine();
        if (ImGui::Button("Show S_k again"))
        {
            clear_reconstruction();
            mesh_dirty_ = true;
        }
    }
    if (running)
        ImGui::TextUnformatted("Reconstructing...");
    if (reconstruction_shown_)
    {
        ImGui::Checkbox("S-eps (blue)", &show_reconstruction_minus_);
        ImGui::SameLine();
        ImGui::Checkbox("S+eps (orange)", &show_reconstruction_plus_);
    }
    if (!reconstruction_message_.empty())
        ImGui::TextWrapped("%s", reconstruction_message_.c_str());
    if (!reconstruction_summary_.empty())
        ImGui::TextWrapped("%s", reconstruction_summary_.c_str());
}

void IsosurfaceViewer::start_reconstruction_job()
{
    iso_reconstruct::Settings settings;
    settings.eps = static_cast<double>(reconstruction_eps_);
    settings.resolution = grid_resolution_;
    settings.halfExtent = static_cast<double>(half_extent_);
    settings.center = {static_cast<double>(point_of_interest_[0]), static_cast<double>(point_of_interest_[1]),
                       static_cast<double>(point_of_interest_[2])};
    settings.patchCells = static_cast<double>(reconstruction_patch_cells_);
    const auto method = static_cast<iso_reconstruct::Method>(reconstruction_method_);
    const std::string expression = last_good_expression_; // the job parses its own copy
    reconstruction_message_.clear();
    reconstruction_job_.start([expression, method, settings] {
        return run_reconstruction_job(expression, method, settings);
    });
}

IsosurfaceViewer::ReconstructionJobResult IsosurfaceViewer::run_reconstruction_job(
    const std::string& expression, iso_reconstruct::Method method, const iso_reconstruct::Settings& settings)
{
    const NodePtr ast = parseExpression(expression);
    const iso_reconstruct::Output out = iso_reconstruct::reconstruct(ast.get(), method, settings);
    if (out.mesh.t.empty())
        throw std::runtime_error("the reconstruction is empty (no surface in the domain)");

    const double h = 2.0 * settings.halfExtent / std::max(settings.resolution, 2);
    const iso_geometry::Point3& P = settings.center;
    const iso_exact::Field f(ast.get(), 0.0);
    ReconstructionJobResult r;
    r.surface = iso_reconstruct::toResult(out.mesh, f);
    const double overlayBox = std::max(12.0 * h, 0.3 * settings.halfExtent);
    r.minus_lines = wireframe(out.minus, P, overlayBox, 0.25f, 0.5f, 1.0f);
    r.plus_lines = wireframe(out.plus, P, overlayBox, 1.0f, 0.55f, 0.15f);

    // Quality near P, and the same for Marching Cubes of F for comparison.
    const iso_reconstruct::Mesh direct =
        iso_reconstruct::fromResult(Isosurface::build(ast.get(), 0.0, settings.resolution, settings.halfExtent), h);
    const double topoBox = std::min(10.0 * h, settings.halfExtent);
    const double errorBox = 4.0 * h;
    char buffer[512];
    std::string text;
    std::snprintf(buffer, sizeof(buffer), "%s, eps = %g, N = %d: %zu triangles, %.0f ms.\n",
                  iso_reconstruct::methodLabel(method), settings.eps, settings.resolution, out.mesh.t.size(), out.millis);
    text += buffer;
    const auto& l = out.local;
    if (l.weightsFound)
    {
        std::snprintf(buffer, sizeof(buffer), "Principal part of F at P: %s, weights (%.3g, %.3g, %.3g)%s.\n",
                      l.principalPart.c_str(), l.weights[0], l.weights[1], l.weights[2],
                      l.nondegenerate ? "" : " (degenerate)");
        text += buffer;
    }
    if (out.patch.applied)
    {
        std::snprintf(buffer, sizeof(buffer), "Cone patch: box radius %.4g, %d loop(s), %d rings.\n", out.patch.radius,
                      out.patch.loops, out.patch.rings);
        text += buffer;
    }
    if (!out.note.empty()) text += "Note: " + out.note + ".\n";
    if (out.correspondenceFallbacks > 0 || out.correspondenceFailures > 0)
    {
        std::snprintf(buffer, sizeof(buffer), "Correspondence: %zu closest-point fallback(s), %zu failure(s).\n",
                      out.correspondenceFallbacks, out.correspondenceFailures);
        text += buffer;
    }
    std::snprintf(buffer, sizeof(buffer), "\nNear P (box of %.0f cells):\n  this mesh: %s\n  MC of F:   %s\n", topoBox / h,
                  topology_text(iso_reconstruct::topology(out.mesh, P, topoBox, h)).c_str(),
                  topology_text(iso_reconstruct::topology(direct, P, topoBox, h)).c_str());
    text += buffer;
    std::snprintf(buffer, sizeof(buffer),
                  "Largest vertex distance to S0 within 4 cells of P:\n  this mesh: %.2e\n  MC of F:   %.2e\n"
                  "(vertices only; a missing piece of surface is not seen by this number.)",
                  vertex_error(out.mesh, f, P, errorBox, h), vertex_error(direct, f, P, errorBox, h));
    text += buffer;
    r.summary = std::move(text);
    return r;
}

void IsosurfaceViewer::poll_reconstruction_job()
{
    reconstruction_job_.poll(
        [this](ReconstructionJobResult&& r) {
            show_surface(std::move(r.surface));
            overlay_.set(OverlayRenderer::PerturbedMinus, std::move(r.minus_lines));
            overlay_.set(OverlayRenderer::PerturbedPlus, std::move(r.plus_lines));
            reconstruction_summary_ = std::move(r.summary);
            reconstruction_shown_ = true;
            runtime_message_ = "Showing the reconstruction of S0, not S_k.";
        },
        [this](const std::string& error) { reconstruction_message_ = "Reconstruction failed: " + error; });
}

void IsosurfaceViewer::clear_reconstruction()
{
    if (!reconstruction_shown_)
        return;
    reconstruction_shown_ = false;
    overlay_.set(OverlayRenderer::PerturbedMinus, {});
    overlay_.set(OverlayRenderer::PerturbedPlus, {});
    if (runtime_message_ == "Showing the reconstruction of S0, not S_k.")
        runtime_message_.clear();
}

} // namespace iso_app
