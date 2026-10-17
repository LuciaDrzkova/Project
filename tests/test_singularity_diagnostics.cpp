// Expanded singularity / isosurface diagnostics for Isosurface Singularity Lab.
//
// This is an additive headless test target. It does not alter the existing
// test_isosurface.cpp. It provides:
//   * a broad, explicitly documented catalogue of parameterized surfaces;
//   * pointwise analytic regression tests for several delicate families;
//   * mesh topology metrics and coordinates for anomalous edges;
//   * grid-phase sensitivity tests implemented by translating F(x,y,z) by a
//     fraction of one cell (equivalent to shifting the grid phase relative to
//     the surface, without changing the production extractor's API);
//   * colored PLY exports so suspicious geometry can be inspected in MeshLab.
//
// Topology findings are diagnostic, not automatically interpreted as proof
// that a mathematical surface is singular. See tests/SINGULARITY_DIAGNOSTICS.md.

#include "Expr.h"
#include "Isosurface.h"
#include "SurfaceCatalogue.h"

#include <algorithm>
#include <array>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

static int g_failures = 0;
#define CHECK(cond, msg) do { if (!(cond)) { std::cerr << "FAIL: " << msg << "\n"; ++g_failures; } } while (false)

constexpr double kDomainHalfExtent = 1.0;

using Behavior = iso_surfaces::Behavior;
using SurfaceCase = iso_surfaces::SurfaceCase;
static const std::vector<SurfaceCase>& catalogue() { return iso_surfaces::catalogue(); }
struct Vec3 { double x, y, z; };
struct Color { unsigned char r, g, b; };

static std::string csvQuote(const std::string& s) {
    std::string out = "\"";
    for (char c : s) { if (c == '"') out += "\"\""; else out += c; }
    out += '"';
    return out;
}
static std::string safeId(std::string s) {
    for (char& c : s) if (!(std::isalnum(static_cast<unsigned char>(c)) || c == '_' || c == '-')) c = '_';
    return s;
}
static std::string number(double x) {
    std::ostringstream s; s << std::setprecision(12) << x; return s.str();
}

static Dual3 evaluate(const Node* ast, double x, double y, double z, double k) {
    return evalExpr<Dual3>(ast,
        Dual3(x, 1, 0, 0), Dual3(y, 0, 1, 0), Dual3(z, 0, 0, 1), Dual3(k, 0, 0, 0));
}
static double gradLength(const Dual3& q) { return std::sqrt(q.dx*q.dx + q.dy*q.dy + q.dz*q.dz); }

static std::string behaviorFor(const SurfaceCase& s, double k) {
    const double eps = 1e-14;
    switch (s.behavior) {
        case Behavior::RegularForNonzeroK:
            return std::fabs(k) <= eps ? "singular/degenerate at k=0 (expected)" : "analytically regular on F=0 for k!=0";
        case Behavior::WhitneyUmbrella:
            if (k > eps) return "singular points (+/-sqrt(k),0,0) lie on F=0";
            if (k < -eps) return "no real Whitney singular points expected on F=0";
            return "singular z-axis at k=0";
        case Behavior::PersistentOrigin:
            return "origin remains singular for every k";
    }
    return "unknown";
}

struct CriticalWitness {
    Vec3 point{};
    std::string locusNote;
};

static double signedNthRoot(double value, int n) {
    return std::copysign(std::pow(std::fabs(value), 1.0/static_cast<double>(n)), value);
}

static bool criticalWitnesses(const SurfaceCase& s, double k,
                              std::vector<CriticalWitness>& out, std::string& rule) {
    const double eps = 1e-14;
    const bool zero = std::fabs(k) <= eps;
    const auto add = [&](double x, double y, double z, const std::string& note = "isolated critical point") {
        out.push_back({Vec3{x,y,z}, note});
    };
    const std::string& id = s.id;

    if (id == "plane_x3_ky") {
        rule = zero ? "origin represents the singular plane x=0" : "dF/dy=-k is nonzero everywhere";
        if (zero) add(0,0,0,"representative of singular plane x=0");
    } else if (id == "isolated_cusp_x" || id == "isolated_cusp_y" || id == "isolated_cusp_z") {
        const char axis = id.back();
        if (k < -eps) {
            const double r = std::sqrt(-k/3.0);
            rule = "two critical points on the cubic perturbed axis, with nonzero critical values";
            if (axis == 'x') { add(r,0,0); add(-r,0,0); }
            else if (axis == 'y') { add(0,r,0); add(0,-r,0); }
            else { add(0,0,r); add(0,0,-r); }
        } else if (zero) { rule = "origin is the singular critical point at k=0"; add(0,0,0); }
        else rule = "the axis derivative cannot vanish for k>0";
    } else if (id == "quartic_cusp_x" || id == "quartic_cusp_y" || id == "quartic_cusp_z") {
        const double r = std::cbrt(-k/4.0);
        rule = "one critical point on the quartic perturbed axis; its value is nonzero for k!=0";
        if (id == "quartic_cusp_x") add(r,0,0);
        else if (id == "quartic_cusp_y") add(0,r,0);
        else add(0,0,r);
    } else if (id == "higher_order_cusp") {
        if (k < -eps) {
            const double r = std::pow(-k/5.0, 0.25);
            rule = "two fourth-root critical points on X, with nonzero values";
            add(r,0,0); add(-r,0,0);
        } else if (zero) { rule = "origin is the singular critical point at k=0"; add(0,0,0); }
        else rule = "5*x^4+k cannot vanish for k>0";
    } else if (id == "mixed_cusp_x") {
        rule = "critical point x=-k/2, with F=-k^2/4"; add(-k/2,0,0);
    } else if (id == "mixed_cusp_y") {
        rule = "critical point y=cbrt(-k/4), with nonzero critical value for k!=0"; add(0,std::cbrt(-k/4.0),0);
    } else if (id == "mixed_cusp_z") {
        rule = "critical point z=cbrt(-k/4), with nonzero critical value for k!=0"; add(0,0,std::cbrt(-k/4.0));
    } else if (id == "mixed_order_z") {
        rule = "critical point z=-k/2, with F=-k^2/4"; add(0,0,-k/2);
    } else if (id == "mixed_cubic" || id == "monkey_saddle_z_unfolding") {
        rule = "critical point (0,0,-k/2), with critical value -k^2/4"; add(0,0,-k/2);
    } else if (id == "level_shift_cone") {
        rule = "origin is the only critical point, with level value -k"; add(0,0,0);
    } else if (id == "whitney_umbrella") {
        if (k > eps) {
            rule = "two real singular points x=+/-sqrt(k) on z=0";
            add(std::sqrt(k),0,0,"positive-k singular point");
            add(-std::sqrt(k),0,0,"negative-x positive-k singular point");
        } else if (zero) {
            rule = "critical z-axis x=y=0; origin is one representative";
            add(0,0,0,"representative of singular z-axis");
        } else rule = "x^2=k has no real solution for k<0";
    } else if (id == "monkey_saddle") {
        rule = "critical z-axis x=y=0; F=-k along it";
        add(0,0,0,"representative of critical z-axis");
    } else if (id == "cubic_family") {
        rule = "origin is singular for every k; at k=1 there is an additional diagonal singular locus";
        add(0,0,0,"persistent singular origin");
        if (std::fabs(k-1.0) <= eps) add(1,1,1,"witness on additional singular line x=y=z at k=1");
    } else if (id == "higher_anisotropic_cusp") {
        rule = "critical point x=seventh_root(k/8), with nonzero critical value for k!=0"; add(signedNthRoot(k/8.0,7),0,0);
    } else if (id == "mixed_even_cusp_x") {
        rule = "critical point x=cbrt(k/4)"; add(std::cbrt(k/4.0),0,0);
    } else if (id == "mixed_even_cusp_y") {
        rule = "critical point y=cbrt(k/4)"; add(0,std::cbrt(k/4.0),0);
    } else if (id == "mixed_even_cusp_z") {
        rule = "critical point z=cbrt(k/4)"; add(0,0,std::cbrt(k/4.0));
    } else if (id == "mixed_order_x") {
        rule = "critical point x=fifth_root(k/6)"; add(signedNthRoot(k/6.0,5),0,0);
    } else if (id == "mixed_order_y") {
        rule = "critical point y=fifth_root(k/6)"; add(0,signedNthRoot(k/6.0,5),0);
    } else if (id == "mixed_even_order_z") {
        rule = "critical point z=k/2, with F=-k^2/4"; add(0,0,k/2);
    } else if (id == "crossing_planes_xy" || id == "quadratic_crossing_unfolding") {
        if (zero) { rule = "origin witnesses the crossing singular locus at k=0"; add(0,0,0,"representative of crossing locus"); }
        else rule = "dF/dz=-k is nonzero everywhere";
    } else if (id == "double_cone_x_shift") {
        rule = "critical point x=k/2, with critical value -k^2/4"; add(k/2,0,0);
    } else if (id == "double_cone_z_shift") {
        rule = "critical point z=-k/2, with critical value +k^2/4"; add(0,0,-k/2);
    } else if (id == "triple_plane_unfolding") {
        if (k > eps) {
            const double r = std::sqrt(k);
            rule = "two diagonal critical points for k>0, with nonzero values";
            add(r,r,r,"positive diagonal critical point"); add(-r,-r,-r,"negative diagonal critical point");
        } else if (zero) { rule = "origin witnesses the coordinate-plane singular axes"; add(0,0,0,"representative of singular coordinate axes"); }
        else rule = "for k<0 the three pairwise products cannot all equal k";
    } else if (id == "parabolic_cusp") {
        if (k < -eps) {
            const double r = std::sqrt(-k/3.0);
            rule = "critical points z=+/-sqrt(-k/3), with nonzero critical values";
            add(0,0,r); add(0,0,-r);
        } else if (zero) { rule = "origin is the singular critical point at k=0"; add(0,0,0); }
        else rule = "-3*z^2-k cannot vanish for k>0";
    } else if (id == "radial_quartic_x") {
        rule = "critical point x=cbrt(k/4)"; add(std::cbrt(k/4.0),0,0);
    } else if (id == "degree_six_isolated") {
        rule = "critical point x=fifth_root(k/6)"; add(signedNthRoot(k/6.0,5),0,0);
    } else if (id == "degree_six_mixed") {
        rule = "critical point z=k/2, with nonzero critical value for k!=0"; add(0,0,k/2);
    } else {
        rule = "no catalog-specific critical-point rule is implemented";
        return false;
    }
    return true;
}

static bool expectedSingularityFor(const SurfaceCase& s, double k) {
    const double eps=1e-14;
    switch(s.behavior) {
        case Behavior::RegularForNonzeroK: return std::fabs(k)<=eps;
        case Behavior::WhitneyUmbrella: return k>=-eps;
        case Behavior::PersistentOrigin: return true;
    }
    return false;
}

static void writeAnalyticAssessment(std::ofstream& f, const SurfaceCase& s, const Node* ast, double k) {
    std::vector<CriticalWitness> witnesses;
    std::string rule;
    const bool ruleKnown=criticalWitnesses(s,k,witnesses,rule);
    const Dual3 origin=evaluate(ast,0,0,0,k);
    const double originGrad=gradLength(origin);
    const double valueTolerance=1e-12*std::max(std::fabs(k),1e-6);
    const double gradientTolerance=1e-8*std::max(std::fabs(k),1e-6);
    std::size_t singularWitnesses=0;
    double maxGradResidual=0.0, maxAbsCriticalValue=0.0;
    std::ostringstream points;
    bool checksPass=ruleKnown;
    for(std::size_t i=0;i<witnesses.size();++i) {
        const auto& w=witnesses[i];
        const Dual3 q=evaluate(ast,w.point.x,w.point.y,w.point.z,k);
        const double gm=gradLength(q);
        maxGradResidual=std::max(maxGradResidual,gm);
        maxAbsCriticalValue=std::max(maxAbsCriticalValue,std::fabs(q.v));
        const bool critical=std::isfinite(gm) && gm<=gradientTolerance;
        const bool onZero=std::isfinite(q.v) && std::fabs(q.v)<=valueTolerance;
        if(critical && onZero) ++singularWitnesses;
        if(!critical) {
            checksPass=false;
            CHECK(false,"critical-point formula failed for surface "+s.id+" k="+number(k)+" at ("+
                number(w.point.x)+","+number(w.point.y)+","+number(w.point.z)+")");
        }
        if(i) points << "; ";
        points << '(' << number(w.point.x) << ',' << number(w.point.y) << ',' << number(w.point.z)
               << " F=" << number(q.v) << " |grad|=" << number(gm) << ", " << w.locusNote << ')';
    }
    const bool expected=expectedSingularityFor(s,k);
    if(ruleKnown) {
        CHECK(expected ? singularWitnesses>0 : singularWitnesses==0,
              "analytic singularity expectation mismatch for surface "+s.id+" k="+number(k)+
              " (rule: "+rule+")");
    }
    std::string classification;
    if(!ruleKnown) classification="UNASSESSED_NO_CRITICAL_POINT_RULE";
    else if(singularWitnesses>0 && expected) classification="EXPECTED_SINGULARITY_WITNESS_FOUND";
    else if(singularWitnesses>0 && !expected) classification="UNEXPECTED_SINGULARITY_WITNESS_FOUND";
    else if(expected) classification="EXPECTED_SINGULARITY_WITNESS_NOT_FOUND";
    else classification="NO_SINGULAR_WITNESS_ON_ZERO_LEVEL";
    f << csvQuote(s.id) << ',' << csvQuote(s.label) << ',' << csvQuote(s.expression) << ',' << k << ','
      << csvQuote(classification) << ',' << csvQuote(rule) << ',' << witnesses.size() << ',' << singularWitnesses << ','
      << maxGradResidual << ',' << maxAbsCriticalValue << ',' << origin.v << ',' << originGrad << ','
      << (std::fabs(origin.v)<=valueTolerance?1:0) << ',' << ((std::fabs(origin.v)<=valueTolerance && originGrad<=gradientTolerance)?1:0) << ','
      << (checksPass?1:0) << ',' << csvQuote(points.str()) << '\n';
    CHECK(ruleKnown,"missing critical-point rule for catalog surface "+s.id);
}

struct EdgeKey {
    int a, b;
    bool operator<(const EdgeKey& other) const { return a < other.a || (a == other.a && b < other.b); }
};
struct EdgeUse { size_t triangle; int orientation; };
struct DiagnosticEdge {
    Vec3 a{}, b{};
    std::size_t incidence = 0;
    int orientationBalance = 0;
    std::string type;
    std::vector<std::size_t> triangles;
};
struct MeshMetrics {
    std::size_t triangles = 0;
    std::size_t weldedVertices = 0;
    std::size_t uniqueEdges = 0;
    std::size_t boundaryEdges = 0;
    std::size_t interiorBoundaryEdges = 0;
    std::size_t nonmanifoldEdges = 0;
    std::size_t orientationConflicts = 0;
    std::size_t nearDegenerateTriangles = 0;
    double minTriangleArea = std::numeric_limits<double>::infinity();
    double minGradient = std::numeric_limits<double>::infinity();
    double maxGradient = 0.0;
    std::vector<DiagnosticEdge> edges;
    std::vector<bool> flaggedTriangle;
};

static long long quantize(double x) { return std::llround(x * 1000000.0); }
using VertexKey = std::array<long long, 3>;
static VertexKey vertexKey(const SurfaceVertex& v) {
    return {quantize(v.px), quantize(v.py), quantize(v.pz)};
}
static double triangleArea(const SurfaceVertex& a, const SurfaceVertex& b, const SurfaceVertex& c) {
    const double ux = b.px-a.px, uy = b.py-a.py, uz = b.pz-a.pz;
    const double vx = c.px-a.px, vy = c.py-a.py, vz = c.pz-a.pz;
    const double cx = uy*vz-uz*vy, cy = uz*vx-ux*vz, cz = ux*vy-uy*vx;
    return 0.5 * std::sqrt(cx*cx+cy*cy+cz*cz);
}
static bool liesOnDomainFace(const Vec3& a, const Vec3& b, double h, double tol) {
    const auto sameFace = [tol](double p, double q, double limit) {
        return std::fabs(std::fabs(p)-limit) <= tol && std::fabs(std::fabs(q)-limit) <= tol && p*q > 0;
    };
    return sameFace(a.x,b.x,h) || sameFace(a.y,b.y,h) || sameFace(a.z,b.z,h);
}

static MeshMetrics analyzeMesh(const IsosurfaceResult& mesh, int resolution, double h) {
    MeshMetrics m;
    m.triangles = mesh.vertices.size()/3;
    m.flaggedTriangle.assign(m.triangles, false);
    const double step = 2.0*h/std::max(1,resolution);
    const double nearAreaThreshold = step*step*1e-5;
    std::map<VertexKey, int> idByPos;
    std::vector<Vec3> positions;
    std::map<EdgeKey, std::vector<EdgeUse>> uses;

    auto getId = [&](const SurfaceVertex& v) {
        const VertexKey key = vertexKey(v);
        auto it = idByPos.find(key);
        if (it != idByPos.end()) return it->second;
        const int id = static_cast<int>(positions.size());
        idByPos.emplace(key, id);
        positions.push_back({v.px, v.py, v.pz});
        return id;
    };

    for (std::size_t t=0; t<m.triangles; ++t) {
        const SurfaceVertex& a = mesh.vertices[t*3];
        const SurfaceVertex& b = mesh.vertices[t*3+1];
        const SurfaceVertex& c = mesh.vertices[t*3+2];
        const int ids[3] = {getId(a),getId(b),getId(c)};
        const double area = triangleArea(a,b,c);
        m.minTriangleArea = std::min(m.minTriangleArea, area);
        if (area < nearAreaThreshold || ids[0] == ids[1] || ids[1] == ids[2] || ids[2] == ids[0]) {
            ++m.nearDegenerateTriangles;
            m.flaggedTriangle[t]=true;
        }
        for (const SurfaceVertex* v : {&a,&b,&c}) {
            m.minGradient = std::min(m.minGradient, static_cast<double>(v->gradMag));
            m.maxGradient = std::max(m.maxGradient, static_cast<double>(v->gradMag));
        }
        for (int e=0; e<3; ++e) {
            const int from=ids[e], to=ids[(e+1)%3];
            if (from == to) continue;
            EdgeKey key{std::min(from,to), std::max(from,to)};
            uses[key].push_back({t, from < to ? +1 : -1});
        }
    }
    m.weldedVertices = positions.size();
    m.uniqueEdges = uses.size();
    for (const auto& kv : uses) {
        DiagnosticEdge e;
        e.a = positions[kv.first.a];
        e.b = positions[kv.first.b];
        e.incidence = kv.second.size();
        for (const auto& u : kv.second) {
            e.orientationBalance += u.orientation;
            e.triangles.push_back(u.triangle);
        }
        if (e.incidence == 1) {
            ++m.boundaryEdges;
            e.type = liesOnDomainFace(e.a,e.b,h,step*1.6) ? "domain_boundary_expected" : "interior_boundary";
            if (e.type == "interior_boundary") {
                ++m.interiorBoundaryEdges;
                for (auto t : e.triangles) m.flaggedTriangle[t]=true;
            }
        } else if (e.incidence > 2) {
            ++m.nonmanifoldEdges;
            e.type = "nonmanifold_edge";
            for (auto t : e.triangles) m.flaggedTriangle[t]=true;
        } else if (e.incidence == 2 && e.orientationBalance != 0) {
            ++m.orientationConflicts;
            e.type = "orientation_conflict";
            for (auto t : e.triangles) m.flaggedTriangle[t]=true;
        } else {
            e.type = "manifold_edge";
        }
        m.edges.push_back(std::move(e));
    }
    if (!std::isfinite(m.minTriangleArea)) m.minTriangleArea = 0.0;
    if (!std::isfinite(m.minGradient)) m.minGradient = 0.0;
    return m;
}

// Re-extract an isosurface from a scalar grid that was reconstructed from the
// two parameter samples. This deliberately mirrors the production edge/table
// path so the reconstructed field is not just compared to itself algebraically.
static IsosurfaceResult extractFromSampledField(const Node* ast, double k, int res, double h,
                                                const std::vector<double>& values) {
    IsosurfaceResult result;
    const auto t0=std::chrono::steady_clock::now();
    if(res<2) res=2;
    const int N=res+1;
    const double step=2.0*h/res;
    const auto idx=[N](int i,int j,int z){return (static_cast<std::size_t>(i)*N+j)*N+z;};
    const auto pos=[&](int i,int j,int z){return std::array<double,3>{-h+i*step,-h+j*step,-h+z*step};};
    const Dual3 kd(k,0,0,0);
    const float degenerateAreaEps=static_cast<float>(step*step)*1e-8f;

    auto cutEdge=[&](const std::array<double,3>& pa,double fa,
                     const std::array<double,3>& pb,double fb) {
        const double denom=fa-fb;
        const double t=(std::fabs(denom)>0.0)?fa/denom:0.5;
        const double px=pa[0]+t*(pb[0]-pa[0]);
        const double py=pa[1]+t*(pb[1]-pa[1]);
        const double pz=pa[2]+t*(pb[2]-pa[2]);
        const Dual3 r=evalExpr<Dual3>(ast,Dual3(px,1,0,0),Dual3(py,0,1,0),Dual3(pz,0,0,1),kd);
        SurfaceVertex v{};
        v.px=static_cast<float>(px); v.py=static_cast<float>(py); v.pz=static_cast<float>(pz);
        const double gm=gradLength(r);
        v.gradMag=static_cast<float>(gm);
        if(gm>1e-12) { v.nx=static_cast<float>(r.dx/gm); v.ny=static_cast<float>(r.dy/gm); v.nz=static_cast<float>(r.dz/gm); }
        else { v.nx=0; v.ny=0; v.nz=1; }
        return v;
    };
    auto emitTri=[&](const SurfaceVertex& a,const SurfaceVertex& b,const SurfaceVertex& c) {
        const float e1x=b.px-a.px,e1y=b.py-a.py,e1z=b.pz-a.pz;
        const float e2x=c.px-a.px,e2y=c.py-a.py,e2z=c.pz-a.pz;
        const float nx=e1y*e2z-e1z*e2y, ny=e1z*e2x-e1x*e2z, nz=e1x*e2y-e1y*e2x;
        const float crossSq=nx*nx+ny*ny+nz*nz;
        if(crossSq<degenerateAreaEps*degenerateAreaEps) {
            ++result.skippedDegenerateTriangles;
            result.skippedDegenerateGeometry.push_back({a.px,a.py,a.pz,b.px,b.py,b.pz,c.px,c.py,c.pz});
            return;
        }
        const float gnx=a.nx+b.nx+c.nx,gny=a.ny+b.ny+c.ny,gnz=a.nz+b.nz+c.nz;
        if(nx*gnx+ny*gny+nz*gnz<0) {
            result.vertices.push_back(a); result.vertices.push_back(c); result.vertices.push_back(b);
        } else {
            result.vertices.push_back(a); result.vertices.push_back(b); result.vertices.push_back(c);
        }
        result.maxGradMag=std::max(result.maxGradMag,static_cast<double>(std::max({a.gradMag,b.gradMag,c.gradMag})));
    };

    for(int i=0;i<res;++i) for(int j=0;j<res;++j) for(int z=0;z<res;++z) {
        std::array<double,3> cp[8];
        double cv[8];
        int cubeindex=0;
        for(int c=0;c<8;++c) {
            const int ci=i+kCubeCornerOffset[c][0], cj=j+kCubeCornerOffset[c][1], cz=z+kCubeCornerOffset[c][2];
            cp[c]=pos(ci,cj,cz);
            cv[c]=values[idx(ci,cj,cz)];
            if(cv[c]<0.0) cubeindex|=(1<<c);
        }
        if(cubeindex==0 || cubeindex==255) continue;
        SurfaceVertex edgeVerts[12]; bool done[12]={false};
        auto edgeVertex=[&](int e)->const SurfaceVertex& {
            if(!done[e]) {
                const int ca=kEdgeCorners[e][0], cb=kEdgeCorners[e][1];
                edgeVerts[e]=cutEdge(cp[ca],cv[ca],cp[cb],cv[cb]);
                done[e]=true;
            }
            return edgeVerts[e];
        };
        const int8_t* tri=kTriTable[cubeindex];
        for(int t=0;tri[t]!=-1;t+=3)
            emitTri(edgeVertex(tri[t]),edgeVertex(tri[t+1]),edgeVertex(tri[t+2]));
    }
    result.buildMillis=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-t0).count();
    return result;
}

static double pointTriangleDistanceSquared(const Vec3& p,const Vec3& a,const Vec3& b,const Vec3& c) {
    auto dot=[](const Vec3& u,const Vec3& v){return u.x*v.x+u.y*v.y+u.z*v.z;};
    auto sub=[](const Vec3& u,const Vec3& v){return Vec3{u.x-v.x,u.y-v.y,u.z-v.z};};
    const Vec3 ab=sub(b,a), ac=sub(c,a), ap=sub(p,a);
    const double d1=dot(ab,ap), d2=dot(ac,ap);
    if(d1<=0.0 && d2<=0.0) return dot(ap,ap);
    const Vec3 bp=sub(p,b); const double d3=dot(ab,bp), d4=dot(ac,bp);
    if(d3>=0.0 && d4<=d3) return dot(bp,bp);
    const double vc=d1*d4-d3*d2;
    if(vc<=0.0 && d1>=0.0 && d3<=0.0) {
        const double v=d1/(d1-d3); const Vec3 q{a.x+v*ab.x,a.y+v*ab.y,a.z+v*ab.z}; const Vec3 e=sub(p,q); return dot(e,e);
    }
    const Vec3 cp=sub(p,c); const double d5=dot(ab,cp), d6=dot(ac,cp);
    if(d6>=0.0 && d5<=d6) return dot(cp,cp);
    const double vb=d5*d2-d1*d6;
    if(vb<=0.0 && d2>=0.0 && d6<=0.0) {
        const double w=d2/(d2-d6); const Vec3 q{a.x+w*ac.x,a.y+w*ac.y,a.z+w*ac.z}; const Vec3 e=sub(p,q); return dot(e,e);
    }
    const double va=d3*d6-d5*d4;
    if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0) {
        const double w=(d4-d3)/((d4-d3)+(d5-d6)); const Vec3 q{b.x+w*(c.x-b.x),b.y+w*(c.y-b.y),b.z+w*(c.z-b.z)}; const Vec3 e=sub(p,q); return dot(e,e);
    }
    const double denom=1.0/(va+vb+vc), v=vb*denom, w=vc*denom;
    const Vec3 q{a.x+ab.x*v+ac.x*w,a.y+ab.y*v+ac.y*w,a.z+ab.z*v+ac.z*w};
    const Vec3 e=sub(p,q); return dot(e,e);
}

static double sampledVertexToMeshDistance(const IsosurfaceResult& source,const IsosurfaceResult& target) {
    // Approximate the directed Hausdorff distance by sampling at most 200
    // source vertices, but test each sampled point against EVERY destination
    // triangle. Subsampling destination triangles can make two identical meshes
    // appear far apart simply because the matching triangle was skipped.
    if(source.vertices.empty() && target.vertices.empty()) return 0.0;
    if(source.vertices.empty() || target.vertices.empty()) return std::numeric_limits<double>::infinity();
    const std::size_t nt=target.vertices.size()/3;
    constexpr std::size_t maxSamples=200;
    const std::size_t sampleCount=std::min(maxSamples,source.vertices.size());
    double maxD2=0.0;
    for(std::size_t sample=0;sample<sampleCount;++sample) {
        const std::size_t i = sampleCount == 1 ? 0 : sample*(source.vertices.size()-1)/(sampleCount-1);
        const auto& v=source.vertices[i]; const Vec3 p{v.px,v.py,v.pz};
        double best=std::numeric_limits<double>::infinity();
        for(std::size_t t=0;t<nt;++t) {
            const auto& va=target.vertices[t*3]; const auto& vb=target.vertices[t*3+1]; const auto& vc=target.vertices[t*3+2];
            best=std::min(best,pointTriangleDistanceSquared(p,{va.px,va.py,va.pz},{vb.px,vb.py,vb.pz},{vc.px,vc.py,vc.pz}));
            if(best<1e-24) break;
        }
        maxD2=std::max(maxD2,best);
    }
    return std::sqrt(maxD2);
}

static double maxNormalAngleDegrees(const IsosurfaceResult& a,const IsosurfaceResult& b,
                                    std::size_t& compared,std::size_t& excluded) {
    compared=0; excluded=0;
    if(a.vertices.size()!=b.vertices.size()) return std::numeric_limits<double>::quiet_NaN();
    double maxAngle=0.0;
    constexpr double gradientEps=1e-8;
    for(std::size_t i=0;i<a.vertices.size();++i) {
        const auto& va=a.vertices[i]; const auto& vb=b.vertices[i];
        const double dx=static_cast<double>(va.px)-vb.px;
        const double dy=static_cast<double>(va.py)-vb.py;
        const double dz=static_cast<double>(va.pz)-vb.pz;
        // Only compare normals when the same triangle-soup index still refers
        // to a geometrically corresponding vertex in both extractions.
        if(dx*dx+dy*dy+dz*dz>1e-10 || va.gradMag<=gradientEps || vb.gradMag<=gradientEps) { ++excluded; continue; }
        double d=static_cast<double>(va.nx)*vb.nx+static_cast<double>(va.ny)*vb.ny+static_cast<double>(va.nz)*vb.nz;
        d=std::max(-1.0,std::min(1.0,d));
        maxAngle=std::max(maxAngle,std::acos(d)*180.0/3.14159265358979323846);
        ++compared;
    }
    return maxAngle;
}

static void writeColoredPly(const fs::path& path, const IsosurfaceResult& mesh, const MeshMetrics& metrics) {
    std::ofstream f(path);
    if (!f) return;
    const std::size_t triCount = mesh.vertices.size()/3;
    f << "ply\nformat ascii 1.0\n"
      << "element vertex " << mesh.vertices.size() << "\n"
      << "property float x\nproperty float y\nproperty float z\n"
      << "property uchar red\nproperty uchar green\nproperty uchar blue\n"
      << "element face " << triCount << "\nproperty list uchar int vertex_indices\nend_header\n";
    for (std::size_t i=0; i<mesh.vertices.size(); ++i) {
        const auto& v = mesh.vertices[i];
        Color c{175,180,190};
        const std::size_t t=i/3;
        if (t < metrics.flaggedTriangle.size() && metrics.flaggedTriangle[t]) c={235,45,45};
        f << v.px << ' ' << v.py << ' ' << v.pz << ' '
          << static_cast<int>(c.r) << ' ' << static_cast<int>(c.g) << ' ' << static_cast<int>(c.b) << '\n';
    }
    for (std::size_t t=0; t<triCount; ++t)
        f << "3 " << t*3 << ' ' << t*3+1 << ' ' << t*3+2 << '\n';
}

static std::string translatedExpression(const std::string& expression, double dx, double dy, double dz) {
    auto fmt = [](double v) { std::ostringstream o; o << std::setprecision(17) << v; return o.str(); };
    std::string out;
    for (std::size_t i=0; i<expression.size();) {
        const unsigned char c = static_cast<unsigned char>(expression[i]);
        if (std::isalpha(c) || expression[i]=='_') {
            std::size_t j=i+1;
            while (j<expression.size()) {
                const unsigned char cj=static_cast<unsigned char>(expression[j]);
                if (!(std::isalnum(cj) || expression[j]=='_')) break;
                ++j;
            }
            const std::string token=expression.substr(i,j-i);
            if (token=="x") out += "(x-("+fmt(dx)+"))";
            else if (token=="y") out += "(y-("+fmt(dy)+"))";
            else if (token=="z") out += "(z-("+fmt(dz)+"))";
            else out += token;
            i=j;
        } else {
            out += expression[i++];
        }
    }
    return out;
}

static void writeSummaryHeader(std::ofstream& f) {
    f << "surface_id,label,expression,k,resolution,phase_x_cells,phase_y_cells,phase_z_cells,translation_dx,translation_dy,translation_dz,analytic_expectation,triangles,welded_vertices,unique_edges,boundary_edges,interior_boundary_edges,nonmanifold_edges,orientation_conflicts,near_degenerate_triangles,skipped_degenerate_triangles,min_triangle_area,min_gradient,max_gradient,build_ms,diagnostic_mesh\n";
}
static void writeMetricRow(std::ofstream& f, const SurfaceCase& s, const std::string& expr,
                           double k, int resolution, double fx, double fy, double fz,
                           double dx, double dy, double dz, const MeshMetrics& m,
                           const IsosurfaceResult& mesh, const std::string& meshPath) {
    f << csvQuote(s.id) << ',' << csvQuote(s.label) << ',' << csvQuote(expr) << ','
      << std::setprecision(12) << k << ',' << resolution << ',' << fx << ',' << fy << ',' << fz << ','
      << dx << ',' << dy << ',' << dz << ',' << csvQuote(behaviorFor(s,k)) << ','
      << m.triangles << ',' << m.weldedVertices << ',' << m.uniqueEdges << ',' << m.boundaryEdges << ','
      << m.interiorBoundaryEdges << ',' << m.nonmanifoldEdges << ',' << m.orientationConflicts << ','
      << m.nearDegenerateTriangles << ',' << mesh.skippedDegenerateTriangles << ','
      << m.minTriangleArea << ',' << m.minGradient << ',' << m.maxGradient << ','
      << mesh.buildMillis << ',' << csvQuote(meshPath) << '\n';
}

static void writeEdgeHeader(std::ofstream& f) {
    f << "surface_id,k,resolution,phase_x_cells,phase_y_cells,phase_z_cells,edge_type,incidence,orientation_balance,x1,y1,z1,x2,y2,z2,triangle_ids\n";
}
static void writeSkippedHeader(std::ofstream& f) {
    f << "surface_id,k,resolution,phase_x_cells,phase_y_cells,phase_z_cells,candidate_id,ax,ay,az,bx,by,bz,cx,cy,cz,area\n";
}
static void writeSkippedRows(std::ofstream& f, const SurfaceCase& s, double k, int res,
                             double fx, double fy, double fz, const IsosurfaceResult& mesh) {
    for (std::size_t i=0;i<mesh.skippedDegenerateGeometry.size();++i) {
        const auto& q=mesh.skippedDegenerateGeometry[i];
        const double ux=q.bx-q.ax, uy=q.by-q.ay, uz=q.bz-q.az;
        const double vx=q.cx-q.ax, vy=q.cy-q.ay, vz=q.cz-q.az;
        const double cx=uy*vz-uz*vy, cy=uz*vx-ux*vz, cz=ux*vy-uy*vx;
        const double area=0.5*std::sqrt(cx*cx+cy*cy+cz*cz);
        f << csvQuote(s.id) << ',' << k << ',' << res << ',' << fx << ',' << fy << ',' << fz << ',' << i << ','
          << q.ax << ',' << q.ay << ',' << q.az << ',' << q.bx << ',' << q.by << ',' << q.bz << ','
          << q.cx << ',' << q.cy << ',' << q.cz << ',' << area << '\n';
    }
}
static void writeEdgeRows(std::ofstream& f, const SurfaceCase& s, double k, int res,
                          double fx, double fy, double fz, const MeshMetrics& m) {
    for (const auto& e : m.edges) {
        if (e.type == "manifold_edge" || e.type == "domain_boundary_expected") continue;
        f << csvQuote(s.id) << ',' << k << ',' << res << ',' << fx << ',' << fy << ',' << fz << ','
          << e.type << ',' << e.incidence << ',' << e.orientationBalance << ','
          << e.a.x << "," << e.a.y << "," << e.a.z << "," << e.b.x << "," << e.b.y << "," << e.b.z << ",\"";
        for (std::size_t i=0;i<e.triangles.size();++i) { if (i) f << ';'; f << e.triangles[i]; }
        f << "\"\n";
    }
}

static const SurfaceCase* findSurface(const std::string& id) {
    for (const auto& s : catalogue()) if (s.id == id) return &s;
    return nullptr;
}

static void runAnalyticSelfTests() {
    // Parser/evaluator smoke test for every catalogue expression at several k values.
    for (const auto& s : catalogue()) {
        try {
            auto ast=parseExpression(s.expression);
            for (double k : {-0.1, 0.0, 0.1}) {
                const Dual3 q=evaluate(ast.get(),0.13,-0.21,0.17,k);
                CHECK(std::isfinite(q.v) && std::isfinite(q.dx) && std::isfinite(q.dy) && std::isfinite(q.dz),
                      "non-finite evaluation for surface " + s.id);
            }
        } catch (const std::exception& e) {
            CHECK(false, "could not parse surface " + s.id + ": " + e.what());
        }
    }
    // F=x^3-k*y: at the origin, grad=(0,-k,0), so only k=0 is singular there.
    {
        auto ast=parseExpression("x^3-k*y");
        for(double k : {-0.1,0.0,0.1}) {
            const auto q=evaluate(ast.get(),0,0,0,k);
            CHECK(std::fabs(q.v)<1e-12, "x^3-k*y should vanish at origin");
            CHECK(std::fabs(q.dx)<1e-12 && std::fabs(q.dz)<1e-12 && std::fabs(q.dy+k)<1e-12,
                  "x^3-k*y origin gradient should be (0,-k,0)");
        }
    }
    // Whitney: for k>0, both (±sqrt(k),0,0) are singular points on F=0.
    {
        auto ast=parseExpression("x^2*z-y^2-k*z");
        for(double k : {0.01,0.1}) for(double sign : {-1.0,1.0}) {
            const auto q=evaluate(ast.get(),sign*std::sqrt(k),0,0,k);
            CHECK(std::fabs(q.v)<1e-10 && gradLength(q)<1e-9,
                  "Whitney umbrella known positive-k singular point must satisfy F=0 and grad F=0");
        }
    }
    // The cubic family retains the singular origin for all k.
    {
        auto ast=parseExpression("x^3+y^3+z^3-3*k*x*y*z");
        for(double k : {-0.1,0.0,0.1}) {
            const auto q=evaluate(ast.get(),0,0,0,k);
            CHECK(std::fabs(q.v)<1e-12 && gradLength(q)<1e-12,
                  "cubic-family origin should remain singular for every k");
        }
    }
    // Mixed cubic: critical point (0,0,-k/2) has value -k^2/4, hence is off F=0 for k!=0.
    // Monkey-saddle unfolding uses the same critical-point/value relationship.
    for (const std::string expr : {"x^3+y^3+z^2+k*z", "x^3-3*x*y^2+z^2+k*z"}) {
        auto ast=parseExpression(expr);
        for(double k : {-0.1,-0.01,0.01,0.1}) {
            const auto q=evaluate(ast.get(),0,0,-k/2,k);
            CHECK(gradLength(q)<1e-9, "known mixed-cubic critical point should have zero gradient");
            CHECK(std::fabs(q.v + k*k/4.0)<1e-10, "known mixed-cubic critical value should be -k^2/4");
            CHECK(std::fabs(q.v)>1e-8, "nonzero-k critical point should not lie on zero level");
        }
    }
    // A simple guaranteed regularization: F=xy-k*z has partial derivative -k in z.
    {
        auto ast=parseExpression("x*y-k*z");
        for(double k : {-0.1,-0.01,0.01,0.1}) {
            const auto q=evaluate(ast.get(),0.23,-0.41,0.07,k);
            CHECK(std::fabs(q.dz+k)<1e-12, "xy-k*z has constant z derivative -k");
        }
    }
}

struct Options {
    int resolution = 48;
    bool quick = false;
    bool runPhase = true;
    bool runInterpolation = true;
    bool customKValues = false;
    std::vector<std::string> surfaces;
    std::vector<std::string> interpolationSurfaces;
    std::vector<double> ks;
    std::vector<double> interpolationMagnitudes;
    bool customInterpolationMagnitudes = false;
    fs::path out = "singularity_diagnostics_artifacts";
};
static Options parseOptions(int argc, char** argv) {
    Options o;
    for(int i=1;i<argc;++i) {
        const std::string a=argv[i];
        if(a=="--resolution" && i+1<argc) o.resolution=std::max(8,std::stoi(argv[++i]));
        else if(a=="--out" && i+1<argc) o.out=argv[++i];
        else if(a=="--surface" && i+1<argc) o.surfaces.push_back(argv[++i]);
        else if(a=="--interpolation-surface" && i+1<argc) o.interpolationSurfaces.push_back(argv[++i]);
        else if(a=="--interpolation-magnitudes" && i+1<argc) {
            o.customInterpolationMagnitudes=true;
            std::stringstream ss(argv[++i]);
            std::string item;
            while(std::getline(ss,item,',')) {
                if(!item.empty()) {
                    const double magnitude=std::stod(item);
                    if(magnitude<=0.0) { std::cerr << "--interpolation-magnitudes requires positive magnitudes\n"; std::exit(2); }
                    o.interpolationMagnitudes.push_back(magnitude);
                }
            }
            if(o.interpolationMagnitudes.empty()) { std::cerr << "--interpolation-magnitudes needs a comma-separated positive list\n"; std::exit(2); }
        }
        else if(a=="--k-values" && i+1<argc) {
            o.customKValues=true;
            std::stringstream ss(argv[++i]);
            std::string item;
            while(std::getline(ss,item,',')) {
                if(!item.empty()) o.ks.push_back(std::stod(item));
            }
            if(o.ks.empty()) { std::cerr << "--k-values needs a comma-separated list, e.g. -0.1,0,0.1\n"; std::exit(2); }
        }
        else if(a=="--quick") o.quick=true;
        else if(a=="--skip-phase") o.runPhase=false;
        else if(a=="--skip-interpolation") o.runInterpolation=false;
        else if(a=="--help" || a=="-h") {
            std::cout << "Usage: isosurface_singularity_diagnostics [--resolution N] [--out DIR] [--quick] [--skip-phase] [--skip-interpolation]\n"
                      << "       [--surface ID (repeatable)] [--k-values -0.1,0,0.1]\n"
                      << "       [--interpolation-surface ID (repeatable)] [--interpolation-magnitudes 0.1,0.01]\n"
                      << "Default: 34 parameterized surfaces x 9 k values, plus grid-phase tests.\n"
                      << "--quick: reduced 12-surface / 3-k smoke suite. Custom --surface and --k-values restrict a run.\n";
            std::exit(0);
        } else {
            std::cerr << "Unknown argument: " << a << "\n"; std::exit(2);
        }
    }
    return o;
}

static bool quickSurface(const std::string& id) {
    static const std::set<std::string> ids={"plane_x3_ky","isolated_cusp_x","quartic_cusp_x","mixed_cubic",
        "level_shift_cone","whitney_umbrella","monkey_saddle","cubic_family","crossing_planes_xy",
        "higher_anisotropic_cusp","radial_quartic_x","triple_plane_unfolding"};
    return ids.count(id)>0;
}
static std::vector<double> kValues(const Options& options) {
    if (options.customKValues) return options.ks;
    if (options.quick) return {-0.1,0.0,0.1};
    return {-0.1,-0.01,-0.001,-0.0001,0.0,0.0001,0.001,0.01,0.1};
}
static bool isFocusSurface(const std::string& id) {
    return id=="mixed_cubic" || id=="mixed_order_z" || id=="mixed_cusp_x" ||
           id=="whitney_umbrella" || id=="plane_x3_ky" || id=="crossing_planes_xy";
}

struct PhaseSnapshot {
    std::array<double,3> phase{};
    std::size_t triangles=0, boundary=0, interiorBoundary=0, nonmanifold=0, orientation=0;
    std::size_t nearDegenerate=0, skippedDegenerate=0;
    double minGradient=0.0, maxGradient=0.0;
};
static void writePhaseVariationHeader(std::ofstream& f) {
    f << "surface_id,k,resolution,phase_count,triangles_min,triangles_max,interior_boundary_min,interior_boundary_max,nonmanifold_min,nonmanifold_max,orientation_conflicts_min,orientation_conflicts_max,near_degenerate_min,near_degenerate_max,skipped_degenerate_min,skipped_degenerate_max,min_gradient_min,min_gradient_max,triangle_count_varies,topology_phase_sensitive\n";
}
static void writePhaseVariationRow(std::ofstream& f, const std::string& id, double k, int res,
                                   const std::vector<PhaseSnapshot>& phases) {
    if(phases.empty()) return;
    auto range=[&](auto member) {
        std::size_t lo=std::numeric_limits<std::size_t>::max(), hi=0;
        for(const auto& p:phases) { lo=std::min(lo,member(p)); hi=std::max(hi,member(p)); }
        return std::pair<std::size_t,std::size_t>{lo,hi};
    };
    const auto tri=range([](const PhaseSnapshot& p){return p.triangles;});
    const auto interior=range([](const PhaseSnapshot& p){return p.interiorBoundary;});
    const auto nm=range([](const PhaseSnapshot& p){return p.nonmanifold;});
    const auto orient=range([](const PhaseSnapshot& p){return p.orientation;});
    const auto near=range([](const PhaseSnapshot& p){return p.nearDegenerate;});
    const auto skipped=range([](const PhaseSnapshot& p){return p.skippedDegenerate;});
    double minGrad=std::numeric_limits<double>::infinity(), maxGrad=0.0;
    for(const auto& p:phases) { minGrad=std::min(minGrad,p.minGradient); maxGrad=std::max(maxGrad,p.minGradient); }
    const bool topologySensitive = interior.first!=interior.second || nm.first!=nm.second ||
        orient.first!=orient.second || skipped.first!=skipped.second;
    f << csvQuote(id) << ',' << k << ',' << res << ',' << phases.size() << ','
      << tri.first << ',' << tri.second << ',' << interior.first << ',' << interior.second << ','
      << nm.first << ',' << nm.second << ',' << orient.first << ',' << orient.second << ','
      << near.first << ',' << near.second << ',' << skipped.first << ',' << skipped.second << ','
      << minGrad << ',' << maxGrad << ',' << (tri.first!=tri.second?1:0) << ',' << (topologySensitive?1:0) << '\n';
}

int main(int argc, char** argv) {
    const Options options=parseOptions(argc,argv);
    runAnalyticSelfTests();
    if(g_failures) {
        std::cerr << g_failures << " analytic regression test(s) failed; diagnostics still run.\n";
    }
    fs::create_directories(options.out);
    fs::create_directories(options.out/"diagnostic_meshes");
    std::ofstream summary(options.out/"summary.csv");
    std::ofstream edgeFile(options.out/"edge_defects.csv");
    std::ofstream skippedFile(options.out/"skipped_degenerate_triangles.csv");
    std::ofstream phaseFile(options.out/"phase_summary.csv");
    std::ofstream phaseVariationFile(options.out/"phase_variation.csv");
    std::ofstream interpolationFile(options.out/"interpolation_validation.csv");
    std::ofstream analyticFile(options.out/"singularity_assessment.csv");
    writeSummaryHeader(summary); writeEdgeHeader(edgeFile); writeSkippedHeader(skippedFile); writeSummaryHeader(phaseFile);
    writePhaseVariationHeader(phaseVariationFile);
    analyticFile << "surface_id,label,expression,k,classification,critical_point_rule,critical_witnesses_tested,critical_witnesses_on_zero_level,max_critical_gradient_residual,max_abs_critical_value,origin_F,origin_gradient_magnitude,origin_on_zero_level,origin_singular_on_zero_level,pointwise_checks_pass,witness_details\n";
    interpolationFile << "surface_id,label,expression,k_negative,k_positive,resolution,grid_samples,max_abs_scalar_error,rms_scalar_error,direct_triangles,reconstructed_triangles,direct_skipped_degenerate,reconstructed_skipped_degenerate,direct_interior_boundary_edges,reconstructed_interior_boundary_edges,direct_nonmanifold_edges,reconstructed_nonmanifold_edges,direct_orientation_conflicts,reconstructed_orientation_conflicts,hausdorff_direct_to_reconstructed_approx,hausdorff_reconstructed_to_direct_approx,symmetric_hausdorff_approx,normal_pairs_compared,normal_pairs_excluded,max_normal_angle_degrees,topology_metrics_match,geometry_agrees\n";

    const int resolution=options.quick ? std::min(options.resolution,24) : options.resolution;
    const double step=2.0*kDomainHalfExtent/resolution;
    const auto ks=kValues(options);
    std::size_t exports=0;
    constexpr std::size_t kMaxMeshExports=180;

    std::cout << "Singularity diagnostics: " << (options.quick ? "quick" : "full")
              << ", resolution=" << resolution << ", catalogue=" << catalogue().size()
              << ", k values=" << ks.size() << "\n";
    for(const auto& s : catalogue()) {
        if(!options.surfaces.empty() && std::find(options.surfaces.begin(),options.surfaces.end(),s.id)==options.surfaces.end()) continue;
        if(options.surfaces.empty() && options.quick && !quickSurface(s.id)) continue;
        std::cout << "  " << s.label << "\n";
        std::unique_ptr<Node> ast;
        try { ast=parseExpression(s.expression); }
        catch(const std::exception& e) { CHECK(false,"catalogue parse failed: "+s.id+": "+e.what()); continue; }
        for(double k : ks) {
            writeAnalyticAssessment(analyticFile,s,ast.get(),k);
            const auto mesh=Isosurface::build(ast.get(),k,resolution,kDomainHalfExtent);
            const auto metrics=analyzeMesh(mesh,resolution,kDomainHalfExtent);
            const bool suspicious=mesh.skippedDegenerateTriangles>0 || metrics.nearDegenerateTriangles>0 ||
                                  metrics.interiorBoundaryEdges>0 || metrics.nonmanifoldEdges>0 || metrics.orientationConflicts>0;
            std::string meshRel;
            const bool exportThis = exports<kMaxMeshExports && isFocusSurface(s.id) &&
                (suspicious || std::fabs(k)<1e-15 || options.quick);
            if(exportThis) {
                const fs::path p=options.out/"diagnostic_meshes"/
                    (safeId(s.id)+"_k"+safeId(number(k))+"_base.ply");
                writeColoredPly(p,mesh,metrics); meshRel=fs::relative(p,options.out).generic_string(); ++exports;
                writeEdgeRows(edgeFile,s,k,resolution,0,0,0,metrics);
                writeSkippedRows(skippedFile,s,k,resolution,0,0,0,mesh);
            }
            writeMetricRow(summary,s,s.expression,k,resolution,0,0,0,0,0,0,metrics,mesh,meshRel);
            if(s.behavior==Behavior::RegularForNonzeroK && std::fabs(k)>1e-14 && mesh.vertices.empty()) {
                // Report, but don't fail: some high-order, tiny-|k| surfaces can lie below this grid's resolution.
                std::cout << "    NOTE: empty extracted mesh for k=" << k << " (resolution may be insufficient)\n";
            }
        }
    }

    if(options.runPhase) {
        const std::vector<std::string> targets = !options.surfaces.empty() ? options.surfaces : (options.quick
            ? std::vector<std::string>{"mixed_cubic","plane_x3_ky","crossing_planes_xy"}
            : std::vector<std::string>{"mixed_cubic","plane_x3_ky","crossing_planes_xy","whitney_umbrella","mixed_order_z","cubic_family"});
        const std::vector<std::array<double,3>> phases = {{{0.0,0.0,0.0}}, {{0.25,0.25,0.25}}, {{0.25,1.0/3.0,0.5}}};
        const std::vector<double> phaseKs = options.customKValues ? ks : (options.quick ? std::vector<double>{0.0,0.1} : std::vector<double>{-0.01,0.0,0.01});
        std::cout << "Grid-phase regression: translated the target surface by fractions of one cell.\n";
        for(const auto& id : targets) {
            const SurfaceCase* s=findSurface(id); if(!s) continue;
            for(double k : phaseKs) {
                std::vector<PhaseSnapshot> phaseSnapshots;
                for(const auto& phase : phases) {
                    const double dx=phase[0]*step, dy=phase[1]*step, dz=phase[2]*step;
                    const std::string expr=translatedExpression(s->expression,dx,dy,dz);
                    std::unique_ptr<Node> ast;
                    try { ast=parseExpression(expr); }
                    catch(const std::exception& e) { CHECK(false,"translated expression parse failed for "+id+": "+e.what()); continue; }
                    const auto mesh=Isosurface::build(ast.get(),k,resolution,kDomainHalfExtent);
                    const auto metrics=analyzeMesh(mesh,resolution,kDomainHalfExtent);
                    std::string meshRel;
                    // Always preserve phase variants of Mixed cubic; these are the most useful visual regression artifacts.
                    const bool exportThis=exports<kMaxMeshExports &&
                        (id=="mixed_cubic" || metrics.nonmanifoldEdges>0 || metrics.orientationConflicts>0 ||
                         metrics.interiorBoundaryEdges>0 || mesh.skippedDegenerateTriangles>0);
                    if(exportThis) {
                        std::ostringstream suffix; suffix<<"_k"<<safeId(number(k))<<"_phase"<<phase[0]<<"_"<<phase[1]<<"_"<<phase[2]<<".ply";
                        const fs::path p=options.out/"diagnostic_meshes"/(safeId(id)+suffix.str());
                        writeColoredPly(p,mesh,metrics); meshRel=fs::relative(p,options.out).generic_string(); ++exports;
                    }
                    writeMetricRow(phaseFile,*s,expr,k,resolution,phase[0],phase[1],phase[2],dx,dy,dz,metrics,mesh,meshRel);
                    writeEdgeRows(edgeFile,*s,k,resolution,phase[0],phase[1],phase[2],metrics);
                    phaseSnapshots.push_back({phase, metrics.triangles, metrics.boundaryEdges, metrics.interiorBoundaryEdges,
                        metrics.nonmanifoldEdges, metrics.orientationConflicts, metrics.nearDegenerateTriangles,
                        mesh.skippedDegenerateTriangles, metrics.minGradient, metrics.maxGradient});
                    if(exportThis) writeSkippedRows(skippedFile,*s,k,resolution,phase[0],phase[1],phase[2],mesh);
                    std::cout << "    " << id << " k=" << k << " phase=(" << phase[0] << ',' << phase[1] << ',' << phase[2]
                              << ") : boundary=" << metrics.boundaryEdges << " interior-boundary=" << metrics.interiorBoundaryEdges
                              << " nonmanifold=" << metrics.nonmanifoldEdges << " orientation-conflicts=" << metrics.orientationConflicts
                              << " skipped-degenerate=" << mesh.skippedDegenerateTriangles << "\n";
                }
                writePhaseVariationRow(phaseVariationFile,id,k,resolution,phaseSnapshots);
                if(!phaseSnapshots.empty()) {
                    const auto& first=phaseSnapshots.front();
                    bool phaseSensitive=false;
                    for(const auto& snap:phaseSnapshots) {
                        if(snap.interiorBoundary!=first.interiorBoundary || snap.nonmanifold!=first.nonmanifold ||
                           snap.orientation!=first.orientation || snap.skippedDegenerate!=first.skippedDegenerate) phaseSensitive=true;
                    }
                    std::cout << "      phase-sensitive topology=" << (phaseSensitive?"YES":"no (for tracked metrics)") << "\n";
                }
            }
        }
    }

    if(options.runInterpolation) {
        const std::vector<std::string> defaultTargets = {"plane_x3_ky","mixed_cubic","whitney_umbrella","monkey_saddle",
            "cubic_family","crossing_planes_xy","quadratic_crossing_unfolding","parabolic_cusp","higher_anisotropic_cusp"};
        const std::vector<std::string> targets = !options.interpolationSurfaces.empty() ? options.interpolationSurfaces :
            (!options.surfaces.empty() ? options.surfaces :
             (options.quick ? std::vector<std::string>{"plane_x3_ky","mixed_cubic","whitney_umbrella","cubic_family"} : defaultTargets));
        const std::vector<double> magnitudes = options.customInterpolationMagnitudes ? options.interpolationMagnitudes :
            (options.quick ? std::vector<double>{0.1} : std::vector<double>{0.1,0.01,0.001});
        std::cout << "Scalar-grid interpolation + reconstructed-mesh comparison:\\n";
        for(const auto& id : targets) {
            const SurfaceCase* s=findSurface(id);
            if(!s) { CHECK(false,"unknown interpolation surface ID: "+id); continue; }
            std::unique_ptr<Node> ast;
            try { ast=parseExpression(s->expression); }
            catch(const std::exception& e) { CHECK(false,"interpolation expression parse failed for "+id+": "+e.what()); continue; }
            for(double a : magnitudes) {
                const double kNeg=-a,kPos=a;
                const int N=resolution+1;
                const double h=kDomainHalfExtent, hStep=2.0*h/resolution;
                const auto idx=[N](int i,int j,int z){return (static_cast<std::size_t>(i)*N+j)*N+z;};
                std::vector<double> reconstructedSamples(static_cast<std::size_t>(N)*N*N);
                double maxFieldError=0.0, sumSquaredError=0.0;
                for(int i=0;i<N;++i) for(int j=0;j<N;++j) for(int z=0;z<N;++z) {
                    const double x=-h+i*hStep,y=-h+j*hStep,zz=-h+z*hStep;
                    const double fm=evalExpr<double>(ast.get(),x,y,zz,kNeg);
                    const double fp=evalExpr<double>(ast.get(),x,y,zz,kPos);
                    const double f0=evalExpr<double>(ast.get(),x,y,zz,0.0);
                    const double reconstructed=0.5*fm+0.5*fp;
                    const double error=std::fabs(reconstructed-f0);
                    maxFieldError=std::max(maxFieldError,error); sumSquaredError+=error*error;
                    reconstructedSamples[idx(i,j,z)]=reconstructed;
                }
                const std::size_t gridSamples=reconstructedSamples.size();
                const double rmsFieldError=std::sqrt(sumSquaredError/std::max<std::size_t>(1,gridSamples));
                const auto direct=Isosurface::build(ast.get(),0.0,resolution,h);
                const auto reconstructed=extractFromSampledField(ast.get(),0.0,resolution,h,reconstructedSamples);
                const auto directMetrics=analyzeMesh(direct,resolution,h);
                const auto reconstructedMetrics=analyzeMesh(reconstructed,resolution,h);
                const double directToReconstructed=sampledVertexToMeshDistance(direct,reconstructed);
                const double reconstructedToDirect=sampledVertexToMeshDistance(reconstructed,direct);
                const double symHausdorff=std::max(directToReconstructed,reconstructedToDirect);
                std::size_t normalsCompared=0,normalsExcluded=0;
                const double maxAngle=maxNormalAngleDegrees(direct,reconstructed,normalsCompared,normalsExcluded);
                const bool topologyMatches=directMetrics.triangles==reconstructedMetrics.triangles &&
                    directMetrics.weldedVertices==reconstructedMetrics.weldedVertices &&
                    directMetrics.uniqueEdges==reconstructedMetrics.uniqueEdges &&
                    directMetrics.boundaryEdges==reconstructedMetrics.boundaryEdges &&
                    directMetrics.interiorBoundaryEdges==reconstructedMetrics.interiorBoundaryEdges &&
                    directMetrics.nonmanifoldEdges==reconstructedMetrics.nonmanifoldEdges &&
                    directMetrics.orientationConflicts==reconstructedMetrics.orientationConflicts &&
                    directMetrics.nearDegenerateTriangles==reconstructedMetrics.nearDegenerateTriangles &&
                    direct.skippedDegenerateTriangles==reconstructed.skippedDegenerateTriangles;
                const double geometricTolerance=std::max(1e-7,hStep*1e-4);
                const bool geometryAgrees=!direct.vertices.empty() && !reconstructed.vertices.empty() && topologyMatches &&
                    std::isfinite(symHausdorff) && symHausdorff<=geometricTolerance;
                interpolationFile << csvQuote(s->id) << ',' << csvQuote(s->label) << ',' << csvQuote(s->expression) << ','
                    << kNeg << ',' << kPos << ',' << resolution << ',' << gridSamples << ',' << maxFieldError << ',' << rmsFieldError << ','
                    << directMetrics.triangles << ',' << reconstructedMetrics.triangles << ',' << direct.skippedDegenerateTriangles << ','
                    << reconstructed.skippedDegenerateTriangles << ',' << directMetrics.interiorBoundaryEdges << ','
                    << reconstructedMetrics.interiorBoundaryEdges << ',' << directMetrics.nonmanifoldEdges << ','
                    << reconstructedMetrics.nonmanifoldEdges << ',' << directMetrics.orientationConflicts << ','
                    << reconstructedMetrics.orientationConflicts << ',' << directToReconstructed << ',' << reconstructedToDirect << ','
                    << symHausdorff << ',' << normalsCompared << ',' << normalsExcluded << ',' << maxAngle << ','
                    << (topologyMatches?1:0) << ',' << (geometryAgrees?1:0) << '\n';
                const bool requiresExport=!geometryAgrees;
                if(requiresExport && exports+2<kMaxMeshExports) {
                    const fs::path dp=options.out/"diagnostic_meshes"/(safeId(id)+"_k0_direct_kpair"+safeId(number(a))+".ply");
                    const fs::path rp=options.out/"diagnostic_meshes"/(safeId(id)+"_k0_reconstructed_kpair"+safeId(number(a))+".ply");
                    writeColoredPly(dp,direct,directMetrics); writeColoredPly(rp,reconstructed,reconstructedMetrics); exports+=2;
                    writeEdgeRows(edgeFile,*s,0.0,resolution,0,0,0,directMetrics);
                    writeEdgeRows(edgeFile,*s,0.0,resolution,0,0,0,reconstructedMetrics);
                    writeSkippedRows(skippedFile,*s,0.0,resolution,0,0,0,direct);
                    writeSkippedRows(skippedFile,*s,0.0,resolution,0,0,0,reconstructed);
                }
                std::cout << "    " << id << " +/-" << a << ": max|field error|=" << maxFieldError
                          << " triangles direct/recon=" << directMetrics.triangles << '/' << reconstructedMetrics.triangles
                          << " approx symmetric Hausdorff=" << symHausdorff
                          << " max normal angle=" << maxAngle << " deg, geometry_agrees=" << (geometryAgrees?"yes":"NO") << "\\n";
            }
        }
    }

    std::cout << "\nWrote:\n  " << (options.out/"summary.csv")
              << "\n  " << (options.out/"phase_summary.csv")
              << "\n  " << (options.out/"phase_variation.csv")
              << "\n  " << (options.out/"interpolation_validation.csv")
              << "\n  " << (options.out/"singularity_assessment.csv")
              << "\n  " << (options.out/"edge_defects.csv")
              << "\n  " << (options.out/"skipped_degenerate_triangles.csv")
              << "\n  " << (options.out/"diagnostic_meshes") << "\n";
    std::cout << "Analytic assertion failures: " << g_failures << "\n";
    // Mesh defects are intentionally diagnostic and do not automatically fail
    // the test process. Analytic regression assertion failures do fail it.
    return g_failures==0 ? 0 : 1;
}
