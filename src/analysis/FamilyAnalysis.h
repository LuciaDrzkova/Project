#pragma once

// One complete analysis of a family H_k = F + kG near the point P.
//
// For every sampled k, each mesh is extracted once (in parallel) and used for
//   * the mesh diagnostics of SingularityConvergence.h,
//   * the mesh-to-mesh comparison of SingularityResearch.h (kept for
//     comparison with earlier results), and
//   * the exact-surface measurements of ExactDistance.h.
// When F and G are polynomials of degree <= 3 it also computes the local
// algebra of F at P (Milnor and Tjurina numbers, essential coordinates of G)
// and, per k, the certified list of stationary points of H_k near P.

#include "ExactDistance.h"
#include "LocalAlgebra.h"
#include "SingularityResearch.h"
#include "StationaryPoints.h"

#include <algorithm>
#include <fstream>
#include <mutex>
#include <optional>
#include <string>
#include <vector>

namespace iso_research {

struct LocalStructure {
    bool available = false;            // F and G are cubic polynomials and P is a singular point of F
    iso_poly::Cubic f, g;              // Taylor polynomials at P (coordinates u = p - P)
    iso_algebra::LocalAlgebra milnor, tjurina;
    std::vector<double> essential;     // coordinates of the class of G in the Tjurina basis
    double essentialNorm = 0.0;        // Euclidean norm of `essential`
    std::string note;
};

struct FamilyAnalysis {
    iso_convergence::Analysis convergence;
    ResearchAnalysis research;
    LocalStructure local;
    std::vector<iso_stationary::Result> stationary;   // one per requested k; points in global coordinates
    std::vector<iso_exact::Measurement> exact;         // one per requested k
};

namespace detail {

// Taylor polynomial of `ast` at `center` (coordinates u = p - center).
inline std::optional<iso_poly::Cubic> cubicAt(const Node* ast, const Point3& center, double k) {
    const auto global = iso_poly::Cubic::fromExpression(ast, k);
    if (!global) return std::nullopt;
    const DifferentialData d = evaluateJet(ast, center, k);
    iso_poly::Cubic local;
    local.c = d.taylor;
    return local;
}

inline LocalStructure localStructure(const iso_convergence::AffineParts& affine, const Point3& probe) {
    LocalStructure s;
    if (!affine.valid) { s.note = "family is not verified to be F + kG"; return s; }
    const NodePtr fAst = parseExpression(affine.f);
    const NodePtr gAst = parseExpression(affine.g);
    const auto f = cubicAt(fAst.get(), probe, 0.0);
    const auto g = cubicAt(gAst.get(), probe, 0.0);
    if (!f || !g) { s.note = "F or G is not a polynomial of degree <= 3"; return s; }
    s.f = *f;
    s.g = *g;
    const double scale = std::max(s.f.maxAbsCoefficient(), 1e-300);
    if (std::fabs(s.f.c[0]) > 1e-12 * scale || std::fabs(s.f.c[1]) + std::fabs(s.f.c[2]) + std::fabs(s.f.c[3]) > 1e-12 * scale) {
        s.note = "P is not a singular point of F";
        return s;
    }
    s.f.c[0] = s.f.c[1] = s.f.c[2] = s.f.c[3] = 0.0; // remove rounding noise of the Taylor evaluation
    s.milnor = iso_algebra::localAlgebra(s.f, iso_algebra::IdealKind::Milnor);
    s.tjurina = iso_algebra::localAlgebra(s.f, iso_algebra::IdealKind::Tjurina);
    if (!s.milnor.valid() || !s.tjurina.valid()) {
        s.note = "the singular point is not isolated (local algebra infinite-dimensional)";
        return s;
    }
    s.essential = s.tjurina.coordinates(s.g);
    for (double v : s.essential) s.essentialNorm += v * v;
    s.essentialNorm = std::sqrt(s.essentialNorm);
    s.available = true;
    return s;
}

inline std::vector<Point3> singularPoints(const iso_stationary::Result& r) {
    std::vector<Point3> out;
    for (const auto& p : r.points)
        if (p.onZeroLevel) out.push_back(p.point);
    return out;
}

inline std::vector<Point3> realStationaryPoints(const iso_stationary::Result& r) {
    std::vector<Point3> out;
    for (const auto& p : r.points)
        if (p.real) out.push_back(p.point);
    return out;
}

// Extra seeds for the exact samples: nested local meshes around P and the
// neighbourhoods of the real stationary points (ExactDistance.h).
inline std::vector<Point3> extraSeeds(const Node* ast, double k, const Point3& probe, double a,
                                      const iso_stationary::Result* stationary) {
    std::vector<Point3> seeds = iso_exact::gradedSeeds(ast, k, probe, a);
    if (stationary) {
        const auto more = iso_exact::stationarySeeds(realStationaryPoints(*stationary), probe, a);
        seeds.insert(seeds.end(), more.begin(), more.end());
    }
    return seeds;
}

} // namespace detail

inline FamilyAnalysis analyzeFamily(const Node* familyAst, const iso_convergence::Options& inputOptions,
                                    int triangleSubdivisions = 3) {
    iso_convergence::Options options = inputOptions;
    options.resolution = std::clamp(options.resolution, 8, Isosurface::kMaxAnalysisResolution);
    const double h = 2.0 * options.halfExtent / options.resolution;
    const Point3 probe = options.probe;
    const double a = options.localHalfSize;

    FamilyAnalysis out;
    out.local = detail::localStructure(iso_convergence::decomposeAffineK(familyAst), probe);

    // The k values exactly as the convergence analysis will use them.
    iso_convergence::Options resolved = options;
    const std::vector<double> kValues = iso_convergence::detail::resolveKValues(resolved);

    // Certified stationary points per k (cheap; no meshes involved).
    std::vector<iso_stationary::Result> stationary(kValues.size());
    iso_stationary::Result referenceStationary;
    if (out.local.available) {
        iso_stationary::Options so;
        so.regionRadius = std::sqrt(3.0) * a; // the local box lies inside this ball
        const int mu = out.local.milnor.dimension;
        auto solve = [&](double k) {
            const auto hk = detail::cubicAt(familyAst, probe, k);
            iso_stationary::Result r = iso_stationary::findStationaryPoints(*hk, mu, so);
            for (auto& p : r.points) p.point = p.point + probe; // back to global coordinates
            return r;
        };
        referenceStationary = solve(0.0);
        iso_parallel::parallelFor(kValues.size(), [&](std::size_t i) { stationary[i] = solve(kValues[i]); });
    }

    // Meshes: one pass feeds all three geometric analyses.
    const iso_exact::Field referenceField(familyAst, 0.0);
    std::once_flag referenceOnce;
    std::optional<iso_exact::Surface> reference;
    std::vector<GeometrySample> geometry(kValues.size());
    std::vector<iso_exact::Measurement> exact(kValues.size());
    out.convergence = iso_convergence::analyze(familyAst, options,
        [&](std::size_t index, double k, const IsosurfaceResult& mesh, const IsosurfaceResult& referenceMesh) {
            geometry[index] = compareMeshes(k, mesh, referenceMesh, probe, a, h, triangleSubdivisions);
            std::call_once(referenceOnce, [&] {
                std::vector<Point3> singular = out.local.available ? detail::singularPoints(referenceStationary)
                                                                   : std::vector<Point3>{probe};
                reference.emplace(iso_exact::prepare(
                    referenceField, referenceMesh, probe, a, h, std::move(singular),
                    detail::extraSeeds(familyAst, 0.0, probe, a, out.local.available ? &referenceStationary : nullptr)));
            });
            const iso_exact::Field field(familyAst, k);
            std::vector<Point3> singular;
            if (out.local.available) singular = detail::singularPoints(stationary[index]);
            else
                for (const CriticalPoint& c : criticalPointsForK(familyAst, k, probe, a))
                    if (c.status.find("NUMERICAL_SINGULAR_CANDIDATE") != std::string::npos) singular.push_back(c.point);
            const iso_exact::Surface surface = iso_exact::prepare(
                field, mesh, probe, a, h, std::move(singular),
                detail::extraSeeds(familyAst, k, probe, a, out.local.available ? &stationary[index] : nullptr));
            exact[index] = iso_exact::measure(k, surface, *reference, mesh, probe, a, h);
        });

    // Research rows and the per-k results follow the requested k order
    // (duplicates included); the convergence run analysed each distinct k once.
    const std::vector<double>& analysed = out.convergence.options.kValues;
    if (options.kValues.empty()) options.kValues = analysed;
    std::vector<GeometrySample> perRequestedK;
    for (double k : options.kValues) {
        const std::size_t i = static_cast<std::size_t>(std::find(analysed.begin(), analysed.end(), k) - analysed.begin());
        perRequestedK.push_back(geometry[i]);
        out.exact.push_back(exact[i]);
        out.stationary.push_back(out.local.available ? stationary[i] : iso_stationary::Result{});
    }
    out.research = analyze(familyAst, options, triangleSubdivisions, &perRequestedK);
    return out;
}

// ----- CSV export ---------------------------------------------------------------

inline bool writeExactGeometryCsv(const std::string& path, const FamilyAnalysis& a, std::string* error);
inline bool writeStationaryCsv(const std::string& path, const FamilyAnalysis& a, std::string* error);
inline bool writeLocalAlgebraCsv(const std::string& path, const FamilyAnalysis& a, std::string* error);

// Every CSV report of one analysis, written into `directory`:
//   singularity_convergence.csv        mesh diagnostics per k
//   singularity_taylor.csv             Taylor data of F and G at P
//   singularity_critical_points.csv    real Newton search (uncertified)
//   singularity_triangle_geometry.csv  mesh-to-mesh distances
//   exact_geometry.csv                 distances on the exact surfaces, extraction error
//   stationary_points_certified.csv    certified stationary points per k
//   local_algebra.csv                  Milnor/Tjurina numbers, essential coordinates of G
inline bool writeFamilyCsvs(const std::string& directory, const FamilyAnalysis& a, std::string* error = nullptr) {
    auto path = [&](const char* name) { return directory + "/" + name; };
    return iso_convergence::writeCsv(path("singularity_convergence.csv"), a.convergence, error) &&
           writeTaylorCsv(path("singularity_taylor.csv"), a.research, error) &&
           writeCriticalPointsCsv(path("singularity_critical_points.csv"), a.research, error) &&
           writeTriangleGeometryCsv(path("singularity_triangle_geometry.csv"), a.research, error) &&
           writeExactGeometryCsv(path("exact_geometry.csv"), a, error) &&
           writeStationaryCsv(path("stationary_points_certified.csv"), a, error) &&
           writeLocalAlgebraCsv(path("local_algebra.csv"), a, error);
}

inline bool writeExactGeometryCsv(const std::string& path, const FamilyAnalysis& a, std::string* error) {
    using iso_csv::number;
    std::ofstream o(path);
    if (!o) { if (error) *error = "Cannot open " + path; return false; }
    o << "k,local_hausdorff_exact,perturbed_to_reference_max,perturbed_to_reference_rms,reference_to_perturbed_max,"
         "reference_to_perturbed_rms,probe_distance_exact,mesh_hausdorff_to_exact_surface,mesh_to_surface_max,"
         "mesh_to_surface_rms,surface_to_mesh_max,surface_to_mesh_rms,singular_point_to_mesh_max,"
         "perturbed_samples,reference_samples,failed_projections,distance_failures,comparable,local_box_half_size,"
         "grid_spacing,method\n";
    const double h = a.convergence.gridSpacing;
    const std::string method = iso_csv::quote(
        "exact surfaces: mesh points moved onto H_k=0 by gradient Newton; point-to-surface distance by Newton on the "
        "Lagrange conditions started at the nearest mesh point, minimum with singular points; maximum-norm local box");
    for (const iso_exact::Measurement& m : a.exact)
        o << number(m.k) << ',' << number(m.symmetric) << ',' << number(m.perturbedToReference.max) << ','
          << number(m.perturbedToReference.rms) << ',' << number(m.referenceToPerturbed.max) << ','
          << number(m.referenceToPerturbed.rms) << ',' << number(m.probeDistance) << ',' << number(m.meshHausdorff)
          << ',' << number(m.meshToSurface.max) << ',' << number(m.meshToSurface.rms) << ','
          << number(m.surfaceToMesh.max) << ',' << number(m.surfaceToMesh.rms) << ',' << number(m.singularPointToMesh)
          << ',' << m.perturbedSamples << ',' << m.referenceSamples << ',' << m.failedProjections << ','
          << (m.perturbedToReference.failures + m.referenceToPerturbed.failures) << ',' << (m.comparable ? 1 : 0)
          << ',' << number(a.convergence.options.localHalfSize) << ',' << number(h) << ',' << method << '\n';
    return static_cast<bool>(o);
}

inline bool writeStationaryCsv(const std::string& path, const FamilyAnalysis& a, std::string* error) {
    using iso_csv::number;
    using iso_csv::quote;
    std::ofstream o(path);
    if (!o) { if (error) *error = "Cannot open " + path; return false; }
    o << "k,status,milnor_number,total_multiplicity,point_index,multiplicity,real,on_zero_level,x,y,z,"
         "re_x,im_x,re_y,im_y,re_z,im_z,H_value,H_term_scale,hessian_corank\n";
    for (std::size_t i = 0; i < a.stationary.size(); ++i) {
        const double k = i < a.exact.size() ? a.exact[i].k : std::numeric_limits<double>::quiet_NaN();
        const iso_stationary::Result& r = a.stationary[i];
        if (r.points.empty())
            o << number(k) << ',' << quote(r.status) << ',' << r.milnorNumber << ',' << r.totalMultiplicity
              << ",,,,,,,,,,,,,,,,\n";
        for (std::size_t j = 0; j < r.points.size(); ++j) {
            const auto& p = r.points[j];
            const Point3 global = p.point;
            o << number(k) << ',' << quote(r.status) << ',' << r.milnorNumber << ',' << r.totalMultiplicity << ','
              << j << ',' << p.multiplicity << ',' << (p.real ? 1 : 0) << ',' << (p.onZeroLevel ? 1 : 0) << ','
              << (p.real ? number(global.x) : "") << ',' << (p.real ? number(global.y) : "") << ','
              << (p.real ? number(global.z) : "") << ','
              << number(p.q[0].real()) << ',' << number(p.q[0].imag()) << ',' << number(p.q[1].real()) << ','
              << number(p.q[1].imag()) << ',' << number(p.q[2].real()) << ',' << number(p.q[2].imag()) << ','
              << number(p.value) << ',' << number(p.valueScale) << ',' << p.hessianCorank << '\n';
        }
    }
    return static_cast<bool>(o);
}

inline bool writeLocalAlgebraCsv(const std::string& path, const FamilyAnalysis& a, std::string* error) {
    using iso_csv::number;
    using iso_csv::quote;
    std::ofstream o(path);
    if (!o) { if (error) *error = "Cannot open " + path; return false; }
    const LocalStructure& s = a.local;
    o << "F_at_P,G_at_P,milnor_number,tjurina_number,tjurina_basis_monomial,G_essential_coordinate,note\n";
    if (!s.available) {
        o << ",,,,,," << quote(s.note) << '\n';
        return static_cast<bool>(o);
    }
    for (std::size_t i = 0; i < s.tjurina.basis.size(); ++i)
        o << quote(s.f.text()) << ',' << quote(s.g.text()) << ',' << s.milnor.dimension << ','
          << s.tjurina.dimension << ',' << quote(s.tjurina.basis[i].text()) << ','
          << number(i < s.essential.size() ? s.essential[i] : 0.0) << ",\n";
    return static_cast<bool>(o);
}

} // namespace iso_research
