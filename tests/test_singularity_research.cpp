#include "core/Expr.h"
#include "core/Isosurface.h"
#include "analysis/SingularityConvergence.h"
#include "analysis/SingularityResearch.h"
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

static int failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++failures; } } while (false)

static iso_convergence::Options optionsFor(const std::vector<double>& ks) {
    iso_convergence::Options o;
    o.probe = {0,0,0}; o.resolution = 16; o.halfExtent = 0.5;
    o.localRadius = 0.22; o.kValues = ks; o.familyExpression = "test family";
    return o;
}
static const iso_research::CriticalPoint* nearestForK(const iso_research::ResearchAnalysis& a,double k) {
    const iso_research::CriticalPoint* best=nullptr;
    for(const auto& c:a.criticalPoints) if(c.k==k && c.status.find("NO_STATIONARY")!=0) {
        if(!best || iso_convergence::norm(c.point)<iso_convergence::norm(best->point)) best=&c;
    }
    return best;
}

int main() {
    const double k=1e-2;

    // Shared reference cubic has a nondegenerate isolated ordinary double point.
    {
        auto ast=parseExpression("x^2+y^2-z^2+x^3+k");
        const auto a=iso_research::analyze(ast.get(),optionsFor({-k,0.0,k}),2);
        CHECK(a.affineStatus=="STRUCTURALLY_AFFINE_IN_K", "A affine decomposition");
        CHECK(a.isolation.status=="NONDEGENERATE_SINGULAR_POINT_LOCAL_ISOLATION_FOLLOWS", "A reference singularity isolated by nonsingular Hessian");
        CHECK(a.fTaylor.firstNonzeroOrder==2, "reference cubic begins at quadratic order at p");
        CHECK(a.gTaylor.firstNonzeroOrder==0, "A perturbation has nonzero constant jet");
        const auto* cp=nearestForK(a,k);
        CHECK(cp!=nullptr, "A stationary point at the origin should be found");
        if(cp) {
            CHECK(std::fabs(cp->point.x)<1e-7 && std::fabs(cp->point.y)<1e-7 && std::fabs(cp->point.z)<1e-7, "A stationary point is p");
            CHECK(cp->status.find("OFF_ZERO_LEVEL")!=std::string::npos, "A critical point leaves the zero level");
        }
        CHECK(a.geometry.size()==3, "A triangle geometry samples written for all selected k values");
    }

    // B: p remains on S_k and becomes regular; nearby stationary point is not on S_k.
    {
        auto ast=parseExpression("x^2+y^2-z^2+x^3+k*y");
        // Include a small nonzero parameter: its critical value is tiny but
        // nonzero, so a fixed absolute field-tolerance floor must not classify
        // it as a singular point on the zero fibre.
        const double smallK=1e-8;
        const auto a=iso_research::analyze(ast.get(),optionsFor({-k,k,smallK}),2);
        CHECK(a.gTaylor.firstNonzeroOrder==1, "B G has first nonzero linear Taylor term");
        const auto* cp=nearestForK(a,k);
        CHECK(cp!=nullptr, "B nearby stationary point should be found");
        if(cp) {
            CHECK(std::fabs(cp->point.x)<1e-6, "B stationary x coordinate");
            CHECK(std::fabs(cp->point.y + k/2.0)<1e-6, "B stationary branch y=-k/2");
            CHECK(std::fabs(cp->point.z)<1e-6, "B stationary z coordinate");
            CHECK(std::fabs(cp->fieldValue + k*k/4.0)<1e-8, "B critical value is -k^2/4");
            CHECK(cp->status.find("OFF_ZERO_LEVEL")!=std::string::npos, "B stationary point is not a surface singularity");
            CHECK(cp->fieldTolerance>0.0 && cp->gradientTolerance>0.0, "B reports separate positive field and gradient tolerances");
        }
        const auto* smallCp=nearestForK(a,smallK);
        CHECK(smallCp!=nullptr, "B small-k stationary point should be found");
        if(smallCp) {
            CHECK(std::fabs(smallCp->fieldValue + smallK*smallK/4.0)<1e-28,
                  "B small-k critical value remains -k^2/4");
            CHECK(smallCp->status.find("NUMERICAL_SINGULAR_CANDIDATE")==std::string::npos,
                  "B small but nonzero critical value must not be misclassified as zero-level singularity");
        }
    }

    // C: G starts at quadratic order; the origin remains a numerical singular candidate.
    {
        auto ast=parseExpression("x^2+y^2-z^2+x^3+k*(x^2+y^2-z^2)");
        const auto a=iso_research::analyze(ast.get(),optionsFor({-k,0.0,k}),2);
        CHECK(a.gTaylor.firstNonzeroOrder==2, "C G has first nonzero quadratic Taylor term");
        const auto* cp=nearestForK(a,k);
        CHECK(cp!=nullptr, "C stationary point should be found");
        if(cp) {
            CHECK(cp->point.x==0.0 && cp->point.y==0.0 && cp->point.z==0.0,
                  "C exact probe point is retained instead of a nearby approximate Newton root");
            CHECK(cp->probePointCandidate, "C singular candidate is marked as directly evaluated at the probe");
            CHECK(cp->searchSource=="DIRECT_PROBE_EVALUATION", "C candidate records direct probe evaluation as its source");
            CHECK(cp->fieldValue==0.0 && cp->gradientMagnitude==0.0,
                  "C exact probe has exactly zero field and gradient residuals");
            CHECK(cp->status.find("NUMERICAL_SINGULAR_CANDIDATE")!=std::string::npos,
                  "C origin remains on the perturbed zero fibre and stationary");
            CHECK(cp->fieldTolerance>0.0 && cp->gradientTolerance>0.0 && cp->deduplicationTolerance>0.0,
                  "C records separate field, gradient and deduplication tolerances");
        }
    }

    // Triangle-to-triangle distance must be zero for identical triangle meshes.
    {
        IsosurfaceResult m;
        auto push=[&](float x,float y,float z){ SurfaceVertex v{}; v.px=x;v.py=y;v.pz=z;m.vertices.push_back(v); };
        push(-0.1f,-0.1f,0.0f); push(0.1f,-0.1f,0.0f); push(0.0f,0.1f,0.0f);
        const auto g=iso_research::compareMeshes(0.1,m,m,{0,0,0},0.5,0.01,3);
        CHECK(g.comparable, "identical triangle meshes comparable");
        CHECK(std::fabs(g.symmetricEstimate)<1e-12, "identical triangle meshes have zero sampled triangle distance");
        CHECK(g.sampleCountPerturbed>3, "triangle interior barycentric samples are used, not only the three vertices");
    }

    if(failures) { std::cerr << failures << " research test(s) failed.\n"; return 1; }
    std::cout << "All singularity research tests passed.\n";
    return 0;
}
