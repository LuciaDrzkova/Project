#pragma once

// Research-grade local diagnostics for F_k(x,y,z)=F(x,y,z)+kG(x,y,z).
//
// This header is additive: it uses the project's Expr.h, Isosurface.h and
// SingularityConvergence.h APIs, but leaves the viewer and existing extractor
// unchanged. All numerical searches are explicitly labelled as non-certified.

#include "SingularityConvergence.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <numeric>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace iso_research {

using Point3 = iso_convergence::Point3;

struct Exponent3 { int x, y, z; };
inline const std::array<Exponent3, 20>& exponents() {
    static const std::array<Exponent3, 20> e{{
        {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1},
        {2,0,0}, {1,1,0}, {1,0,1}, {0,2,0}, {0,1,1}, {0,0,2},
        {3,0,0}, {2,1,0}, {2,0,1}, {1,2,0}, {1,1,1},
        {1,0,2}, {0,3,0}, {0,2,1}, {0,1,2}, {0,0,3}
    }};
    return e;
}
inline int exponentIndex(int x, int y, int z) {
    const auto& e = exponents();
    for (int i = 0; i < static_cast<int>(e.size()); ++i)
        if (e[static_cast<std::size_t>(i)].x == x &&
            e[static_cast<std::size_t>(i)].y == y &&
            e[static_cast<std::size_t>(i)].z == z) return i;
    return -1;
}
inline int factorial(int n) { return n < 2 ? 1 : (n == 2 ? 2 : 6); }

// A normalized multivariate Taylor jet through total degree 3.
// c[i,j,k] stores (partial^(i+j+k) f)/(i!j!k!).
struct Jet3 {
    std::array<double, 20> c{};
    bool smooth = true;
    Jet3() = default;
    explicit Jet3(double value) { c[0] = value; }
    static Jet3 variable(double value, int axis) {
        Jet3 q(value);
        if (axis >= 0 && axis < 3) q.c[static_cast<std::size_t>(1 + axis)] = 1.0;
        return q;
    }
    double value() const { return c[0]; }
};
inline Jet3 operator+(const Jet3& a, const Jet3& b) {
    Jet3 r; r.smooth = a.smooth && b.smooth;
    for (std::size_t i = 0; i < r.c.size(); ++i) r.c[i] = a.c[i] + b.c[i];
    return r;
}
inline Jet3 operator-(const Jet3& a, const Jet3& b) {
    Jet3 r; r.smooth = a.smooth && b.smooth;
    for (std::size_t i = 0; i < r.c.size(); ++i) r.c[i] = a.c[i] - b.c[i];
    return r;
}
inline Jet3 operator-(const Jet3& a) {
    Jet3 r; r.smooth = a.smooth;
    for (std::size_t i = 0; i < r.c.size(); ++i) r.c[i] = -a.c[i];
    return r;
}
inline Jet3 operator*(const Jet3& a, const Jet3& b) {
    Jet3 r; r.smooth = a.smooth && b.smooth;
    const auto& e = exponents();
    for (std::size_t i = 0; i < e.size(); ++i) {
        if (a.c[i] == 0.0) continue;
        for (std::size_t j = 0; j < e.size(); ++j) {
            if (b.c[j] == 0.0) continue;
            const int x = e[i].x + e[j].x;
            const int y = e[i].y + e[j].y;
            const int z = e[i].z + e[j].z;
            if (x + y + z > 3) continue;
            const int id = exponentIndex(x, y, z);
            if (id >= 0) r.c[static_cast<std::size_t>(id)] += a.c[i] * b.c[j];
        }
    }
    return r;
}
inline Jet3 composeUnary(const Jet3& a, double f0, double f1, double f2, double f3) {
    Jet3 u = a; u.c[0] = 0.0;
    Jet3 r(f0);
    r = r + u * Jet3(f1);
    r = r + (u * u) * Jet3(0.5 * f2);
    r = r + (u * u * u) * Jet3(f3 / 6.0);
    r.smooth = a.smooth;
    return r;
}
inline Jet3 reciprocal(const Jet3& a) {
    if (a.value() == 0.0 || !std::isfinite(a.value())) {
        Jet3 r(std::numeric_limits<double>::quiet_NaN()); r.smooth = false; return r;
    }
    const double x = a.value();
    return composeUnary(a, 1.0/x, -1.0/(x*x), 2.0/(x*x*x), -6.0/(x*x*x*x));
}
inline Jet3 operator/(const Jet3& a, const Jet3& b) { return a * reciprocal(b); }
inline Jet3 mc_sin(const Jet3& a) {
    const double x=a.value(); return composeUnary(a,std::sin(x),std::cos(x),-std::sin(x),-std::cos(x));
}
inline Jet3 mc_cos(const Jet3& a) {
    const double x=a.value(); return composeUnary(a,std::cos(x),-std::sin(x),-std::cos(x),std::sin(x));
}
inline Jet3 mc_tan(const Jet3& a) {
    const double x=a.value(), t=std::tan(x), sec2=1.0/(std::cos(x)*std::cos(x));
    return composeUnary(a,t,sec2,2.0*t*sec2,2.0*sec2*(sec2+2.0*t*t));
}
inline Jet3 mc_exp(const Jet3& a) {
    const double e=std::exp(a.value()); return composeUnary(a,e,e,e,e);
}
inline Jet3 mc_log(const Jet3& a) {
    const double x=a.value();
    if (!(x > 0.0)) { Jet3 r(std::numeric_limits<double>::quiet_NaN()); r.smooth=false; return r; }
    return composeUnary(a,std::log(x),1.0/x,-1.0/(x*x),2.0/(x*x*x));
}
inline Jet3 mc_sqrt(const Jet3& a) {
    const double x=a.value();
    if (!(x > 0.0)) {
        if (x == 0.0 && std::all_of(a.c.begin()+1,a.c.end(),[](double v){return v==0.0;})) return Jet3(0.0);
        Jet3 r(std::numeric_limits<double>::quiet_NaN()); r.smooth=false; return r;
    }
    const double r=std::sqrt(x);
    return composeUnary(a,r,0.5/r,-0.25/(x*r),0.375/(x*x*r));
}
inline Jet3 mc_abs(const Jet3& a) {
    if (a.value() == 0.0) { Jet3 r(std::fabs(a.value())); r.smooth=false; return r; }
    return a.value() > 0.0 ? a : -a;
}
inline Jet3 mc_pow(const Jet3& a, const Jet3& b) {
    bool exponentConstant = true;
    for (std::size_t i=1;i<b.c.size();++i) exponentConstant = exponentConstant && (b.c[i] == 0.0);
    if (!exponentConstant) {
        if (!(a.value() > 0.0)) { Jet3 r(std::numeric_limits<double>::quiet_NaN()); r.smooth=false; return r; }
        return mc_exp(b * mc_log(a));
    }
    const double p=b.value(), x=a.value();
    if (x == 0.0 && p < 0.0) { Jet3 r(std::numeric_limits<double>::quiet_NaN()); r.smooth=false; return r; }
    if (x < 0.0 && std::fabs(p-std::round(p)) > 1e-12) { Jet3 r(std::numeric_limits<double>::quiet_NaN()); r.smooth=false; return r; }
    const double f0=std::pow(x,p);
    const double f1=(p==0.0)?0.0:p*std::pow(x,p-1.0);
    const double f2=(p==0.0||p==1.0)?0.0:p*(p-1.0)*std::pow(x,p-2.0);
    const double f3=(p==0.0||p==1.0||p==2.0)?0.0:p*(p-1.0)*(p-2.0)*std::pow(x,p-3.0);
    return composeUnary(a,f0,f1,f2,f3);
}
inline Jet3 mc_min(const Jet3& a, const Jet3& b) {
    if (a.value()==b.value()) { Jet3 r=a; r.smooth=false; return r; }
    return a.value()<b.value()?a:b;
}
inline Jet3 mc_max(const Jet3& a, const Jet3& b) {
    if (a.value()==b.value()) { Jet3 r=a; r.smooth=false; return r; }
    return a.value()>b.value()?a:b;
}

struct DifferentialData {
    double value = std::numeric_limits<double>::quiet_NaN();
    Point3 gradient{};
    std::array<std::array<double,3>,3> hessian{};
    std::array<double,10> third{}; // xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz
    std::array<double,20> taylor{};
    bool smooth = false;
};
inline DifferentialData extractDifferential(const Jet3& j) {
    DifferentialData d; d.value=j.c[0]; d.smooth=j.smooth; d.taylor=j.c;
    d.gradient={j.c[1],j.c[2],j.c[3]};
    d.hessian[0][0]=2*j.c[4]; d.hessian[0][1]=d.hessian[1][0]=j.c[5];
    d.hessian[0][2]=d.hessian[2][0]=j.c[6]; d.hessian[1][1]=2*j.c[7];
    d.hessian[1][2]=d.hessian[2][1]=j.c[8]; d.hessian[2][2]=2*j.c[9];
    const auto& e=exponents();
    for (std::size_t i=10;i<20;++i) {
        const int mult=factorial(e[i].x)*factorial(e[i].y)*factorial(e[i].z);
        d.third[i-10]=j.c[i]*mult;
    }
    return d;
}
inline DifferentialData evaluateJet(const Node* ast, const Point3& p, double k) {
    const Jet3 x=Jet3::variable(p.x,0), y=Jet3::variable(p.y,1), z=Jet3::variable(p.z,2), kk(k);
    return extractDifferential(evalExpr<Jet3>(ast,x,y,z,kk));
}
inline double gradientNorm(const Point3& g) { return std::sqrt(g.x*g.x+g.y*g.y+g.z*g.z); }
inline double determinant3(const std::array<std::array<double,3>,3>& a) {
    return a[0][0]*(a[1][1]*a[2][2]-a[1][2]*a[2][1])
         - a[0][1]*(a[1][0]*a[2][2]-a[1][2]*a[2][0])
         + a[0][2]*(a[1][0]*a[2][1]-a[1][1]*a[2][0]);
}
inline std::array<double,3> symmetricEigenvalues(std::array<std::array<double,3>,3> a) {
    for (int sweep=0;sweep<24;++sweep) {
        int p=0,q=1;
        if (std::fabs(a[0][2])>std::fabs(a[p][q])) {p=0;q=2;}
        if (std::fabs(a[1][2])>std::fabs(a[p][q])) {p=1;q=2;}
        if (std::fabs(a[p][q])<1e-15) break;
        const double tau=(a[q][q]-a[p][p])/(2.0*a[p][q]);
        const double t=std::copysign(1.0,tau)/(std::fabs(tau)+std::sqrt(1.0+tau*tau));
        const double c=1.0/std::sqrt(1.0+t*t), s=t*c;
        const double app=a[p][p], aqq=a[q][q], apq=a[p][q];
        a[p][p]=app-t*apq; a[q][q]=aqq+t*apq; a[p][q]=a[q][p]=0.0;
        for (int r=0;r<3;++r) if(r!=p && r!=q) {
            const double arp=a[r][p], arq=a[r][q];
            a[r][p]=a[p][r]=c*arp-s*arq;
            a[r][q]=a[q][r]=s*arp+c*arq;
        }
    }
    std::array<double,3> e{{a[0][0],a[1][1],a[2][2]}};
    std::sort(e.begin(),e.end()); return e;
}
inline int firstNonzeroJetOrder(const std::array<double,20>& c) {
    double scale=0.0; for(double v:c) scale=std::max(scale,std::fabs(v));
    if (scale==0.0) return -1;
    if (scale==0.0) return -1;
    const double tol=512.0*std::numeric_limits<double>::epsilon()*scale;
    for (int i=0;i<20;++i) if(std::fabs(c[static_cast<std::size_t>(i)])>tol) {
        const auto e=exponents()[static_cast<std::size_t>(i)]; return e.x+e.y+e.z;
    }
    return -1;
}
inline double evalAbsScale(const Node* n, const Point3& p, double k) {
    if (!n) return 0.0;
    switch(n->type) {
        case Node::NUM: return std::fabs(n->num);
        case Node::VAR_X: return std::fabs(p.x);
        case Node::VAR_Y: return std::fabs(p.y);
        case Node::VAR_Z: return std::fabs(p.z);
        case Node::VAR_K: return std::fabs(k);
        case Node::NEG: return evalAbsScale(n->kids[0].get(),p,k);
        case Node::ADD: case Node::SUB:
            return evalAbsScale(n->kids[0].get(),p,k)+evalAbsScale(n->kids[1].get(),p,k);
        case Node::MUL:
            return evalAbsScale(n->kids[0].get(),p,k)*evalAbsScale(n->kids[1].get(),p,k);
        case Node::DIV: {
            const double d=evalAbsScale(n->kids[1].get(),p,k);
            return evalAbsScale(n->kids[0].get(),p,k)/std::max(d,1e-300);
        }
        case Node::POW: {
            const double a=evalAbsScale(n->kids[0].get(),p,k);
            const double b=evalExpr<double>(n->kids[1].get(),p.x,p.y,p.z,k);
            const double v=std::pow(std::max(a,1e-300),b);
            return std::isfinite(v)?v:std::numeric_limits<double>::max()/16.0;
        }
        case Node::CALL: {
            const double a=evalAbsScale(n->kids[0].get(),p,k);
            if(n->fname=="sin"||n->fname=="cos") return 1.0;
            if(n->fname=="exp") { const double v=std::exp(a); return std::isfinite(v)?v:std::numeric_limits<double>::max()/16.0; }
            if(n->fname=="sqrt") return std::sqrt(a);
            if(n->fname=="abs") return a;
            if(n->fname=="log"||n->fname=="ln") return std::max(1.0,std::fabs(std::log(std::max(a,1e-300))));
            if(n->kids.size()==2) return std::max(evalAbsScale(n->kids[0].get(),p,k),evalAbsScale(n->kids[1].get(),p,k));
            return std::max(1.0,a);
        }
    }
    return 1.0;
}

struct IsolationAssessment {
    std::string status = "UNKNOWN";
    double fieldValue = std::numeric_limits<double>::quiet_NaN();
    double gradientMagnitude = std::numeric_limits<double>::quiet_NaN();
    double hessianDeterminant = std::numeric_limits<double>::quiet_NaN();
    std::array<double,3> hessianEigenvalues{};
    std::string explanation;
};
inline IsolationAssessment assessReferenceSingularity(const Node* fAst, const Point3& p, double radius) {
    IsolationAssessment out;
    const DifferentialData d=evaluateJet(fAst,p,0.0);
    out.fieldValue=d.value; out.gradientMagnitude=gradientNorm(d.gradient);
    out.hessianDeterminant=determinant3(d.hessian);
    out.hessianEigenvalues=symmetricEigenvalues(d.hessian);
    const double scale=evalAbsScale(fAst,p,0.0);
    const double fTol=512.0*std::numeric_limits<double>::epsilon()*scale+1e-30;
    const double hScale=std::max({std::fabs(d.hessian[0][0]),std::fabs(d.hessian[0][1]),std::fabs(d.hessian[0][2]),
                                  std::fabs(d.hessian[1][1]),std::fabs(d.hessian[1][2]),std::fabs(d.hessian[2][2])});
    const double gTol=1e-10*std::max(hScale*std::max(radius,1e-6),1e-12);
    if (!d.smooth) {
        out.status="NOT_CERTIFIED_NONSMOOTH_OR_UNSUPPORTED";
        out.explanation="The expression is non-smooth or its third-order jet is unsupported at the probe point.";
        return out;
    }
    if (std::fabs(d.value)>fTol || out.gradientMagnitude>gTol) {
        out.status="PROBE_NOT_CONFIRMED_SINGULAR";
        out.explanation="The probe does not satisfy F(p)=0 and grad F(p)=0 within scale-aware floating-point tolerances.";
        return out;
    }
    const double detTol=512.0*std::numeric_limits<double>::epsilon()*std::max(hScale*hScale*hScale,1e-300);
    if (std::fabs(out.hessianDeterminant)>detTol) {
        out.status="NONDEGENERATE_SINGULAR_POINT_LOCAL_ISOLATION_FOLLOWS";
        out.explanation="The probe is numerically a singular point and det(H_F(p)) is nonzero. For exact input data, the inverse function theorem makes p an isolated zero of grad F, hence an isolated hypersurface singularity locally.";
        return out;
    }
    // A coarse witness scan may detect another singular point, but cannot prove
    // non-isolation. Degenerate Hessian cases remain explicitly inconclusive.
    int extraWitnesses=0;
    for(int i=-2;i<=2;++i) for(int j=-2;j<=2;++j) for(int l=-2;l<=2;++l) {
        const Point3 q{p.x+radius*0.5*i,p.y+radius*0.5*j,p.z+radius*0.5*l};
        if(iso_convergence::distance(q,p)<1e-7*std::max(1.0,radius)) continue;
        const auto v=evaluateJet(fAst,q,0.0);
        const double s=evalAbsScale(fAst,q,0.0);
        const double ft=512.0*std::numeric_limits<double>::epsilon()*s+1e-30;
        if(std::fabs(v.value)<=ft && gradientNorm(v.gradient)<=1e-10*std::max(1.0,hScale)) ++extraWitnesses;
    }
    if(extraWitnesses>0) {
        out.status="NONISOLATED_SINGULAR_WITNESSES_DETECTED_NOT_CERTIFIED";
        out.explanation="Additional grid witnesses approximately satisfy F=0 and grad F=0. This is evidence of a non-isolated singular locus, not a proof.";
    } else {
        out.status="DEGENERATE_HESSIAN_ISOLATION_INCONCLUSIVE";
        out.explanation="The Hessian is singular. A singular Hessian neither proves nor disproves isolatedness; higher-order/algebraic analysis is required.";
    }
    return out;
}

struct CriticalPoint {
    double k=0.0;
    Point3 point{};
    double fieldValue=std::numeric_limits<double>::quiet_NaN();
    double fieldScale=std::numeric_limits<double>::quiet_NaN();
    double fieldTolerance=std::numeric_limits<double>::quiet_NaN();
    double gradientMagnitude=std::numeric_limits<double>::quiet_NaN();
    double gradientTolerance=std::numeric_limits<double>::quiet_NaN();
    double newtonStepTolerance=std::numeric_limits<double>::quiet_NaN();
    double deduplicationTolerance=std::numeric_limits<double>::quiet_NaN();
    double hessianDeterminant=std::numeric_limits<double>::quiet_NaN();
    int iterations=0;
    bool probePointCandidate=false;
    std::string searchSource="MULTISTART_DAMPED_NEWTON";
    std::string status;
};

struct CriticalPointTolerances {
    double field = std::numeric_limits<double>::denorm_min() * 64.0;
    double gradient = std::numeric_limits<double>::denorm_min() * 64.0;
    double newtonStep = 0.0;
    double deduplication = 0.0;
};

inline double hessianMaxAbs(const std::array<std::array<double,3>,3>& h) {
    double m=0.0;
    for(const auto& row:h) for(double v:row) m=std::max(m,std::fabs(v));
    return m;
}

// These tolerances have different meanings and are intentionally not shared:
// field/gradient residuals test equations, step tolerance tests Newton progress,
// and deduplication tolerance decides whether two numerically found roots are
// the same candidate. No arbitrary 1e-30 field floor is added: that floor can
// incorrectly turn a small but nonzero critical value into a zero-level root.
inline CriticalPointTolerances criticalPointTolerances(
    const Node* ast,double k,const Point3& q,const Point3& center,double radius,
    const DifferentialData& d) {
    CriticalPointTolerances t;
    const double eps=std::numeric_limits<double>::epsilon();
    const double denormFloor=64.0*std::numeric_limits<double>::denorm_min();
    const double safeRadius=std::max(std::fabs(radius),1e-12);
    const double localLength=std::max(safeRadius,iso_convergence::distance(q,center));
    const double hScale=hessianMaxAbs(d.hessian);
    const double coordinateScale=std::max({std::fabs(q.x),std::fabs(q.y),std::fabs(q.z),safeRadius});
    const double fieldScale=evalAbsScale(ast,q,k);
    t.field=std::max(denormFloor,512.0*eps*fieldScale);
    t.gradient=std::max(denormFloor,512.0*eps*hScale*localLength);
    t.newtonStep=std::max(128.0*eps*std::max(coordinateScale,1e-12),
                          1e-12*safeRadius);
    t.deduplication=std::max(512.0*eps*std::max(coordinateScale,1e-12),
                              1e-10*safeRadius);
    return t;
}

inline double criticalPointQuality(const CriticalPoint& c) {
    const double gTol=std::max(c.gradientTolerance,64.0*std::numeric_limits<double>::denorm_min());
    const double fTol=std::max(c.fieldTolerance,64.0*std::numeric_limits<double>::denorm_min());
    const double gr=c.gradientMagnitude/gTol;
    const double fr=std::fabs(c.fieldValue)/fTol;
    return std::hypot(gr,fr);
}

inline void addOrImproveCriticalPoint(std::vector<CriticalPoint>& candidates,CriticalPoint candidate) {
    auto it=std::find_if(candidates.begin(),candidates.end(),[&](const CriticalPoint& existing) {
        const double tol=std::max(existing.deduplicationTolerance,candidate.deduplicationTolerance);
        return iso_convergence::distance(existing.point,candidate.point)<=tol;
    });
    if(it==candidates.end()) {
        candidates.push_back(std::move(candidate));
        return;
    }
    // Keep the candidate with the smallest normalized pair of residuals. If
    // quality ties, prefer direct evaluation of the selected probe point.
    const double oldQuality=criticalPointQuality(*it);
    const double newQuality=criticalPointQuality(candidate);
    if(newQuality<oldQuality ||
       (newQuality==oldQuality && candidate.probePointCandidate && !it->probePointCandidate)) {
        *it=std::move(candidate);
    }
}

inline CriticalPoint makeCriticalPoint(const Node* ast,double k,const Point3& q,
                                       const Point3& probe,double radius,int iterations,
                                       bool isProbePoint,const std::string& source) {
    const auto d=evaluateJet(ast,q,k);
    const auto tolerances=criticalPointTolerances(ast,k,q,probe,radius,d);
    CriticalPoint row;
    row.k=k; row.point=q; row.fieldValue=d.value;
    row.fieldScale=evalAbsScale(ast,q,k); row.fieldTolerance=tolerances.field;
    row.gradientMagnitude=gradientNorm(d.gradient); row.gradientTolerance=tolerances.gradient;
    row.newtonStepTolerance=tolerances.newtonStep;
    row.deduplicationTolerance=tolerances.deduplication;
    row.hessianDeterminant=determinant3(d.hessian); row.iterations=iterations;
    row.probePointCandidate=isProbePoint; row.searchSource=source;
    const bool stationary=row.gradientMagnitude<=row.gradientTolerance;
    const bool onLevel=std::fabs(row.fieldValue)<=row.fieldTolerance;
    if(stationary && onLevel) {
        row.status="NUMERICAL_SINGULAR_CANDIDATE_ON_ZERO_LEVEL_NOT_CERTIFIED";
    } else if(stationary) {
        row.status="STATIONARY_POINT_OFF_ZERO_LEVEL";
    } else {
        row.status="STATIONARY_RESIDUAL_ABOVE_SCALE_AWARE_TOLERANCE";
    }
    return row;
}

inline bool solve3x3(std::array<std::array<double,3>,3> a, Point3 b, Point3& x) {
    double m[3][4]={{a[0][0],a[0][1],a[0][2],b.x},
                    {a[1][0],a[1][1],a[1][2],b.y},
                    {a[2][0],a[2][1],a[2][2],b.z}};
    double matrixScale=0.0;
    for(int r=0;r<3;++r) for(int c=0;c<3;++c) matrixScale=std::max(matrixScale,std::fabs(m[r][c]));
    const double pivotTolerance=64.0*std::numeric_limits<double>::epsilon()*matrixScale;
    if(!(matrixScale>0.0) || !std::isfinite(matrixScale)) return false;
    for(int c=0;c<3;++c) {
        int pivot=c;
        for(int r=c+1;r<3;++r) if(std::fabs(m[r][c])>std::fabs(m[pivot][c])) pivot=r;
        if(std::fabs(m[pivot][c])<=pivotTolerance) return false;
        if(pivot!=c) for(int j=c;j<4;++j) std::swap(m[pivot][j],m[c][j]);
        const double div=m[c][c]; for(int j=c;j<4;++j) m[c][j]/=div;
        for(int r=0;r<3;++r) if(r!=c) { const double f=m[r][c]; for(int j=c;j<4;++j) m[r][j]-=f*m[c][j]; }
    }
    x={m[0][3],m[1][3],m[2][3]};
    return std::isfinite(x.x)&&std::isfinite(x.y)&&std::isfinite(x.z);
}

inline bool newtonStationary(const Node* ast,double k,const Point3& center,double radius,
                             Point3 seed,Point3& root,int& iterations,double& gradResidual) {
    Point3 x=seed;
    for(int it=0;it<60;++it) {
        const auto d=evaluateJet(ast,x,k);
        const double gn=gradientNorm(d.gradient);
        const auto tolerances=criticalPointTolerances(ast,k,x,center,radius,d);
        iterations=it;
        if(!std::isfinite(gn) || !d.smooth) return false;
        if(gn<=tolerances.gradient) { root=x; gradResidual=gn; return true; }
        Point3 delta;
        if(!solve3x3(d.hessian,d.gradient,delta)) return false;
        const double dn=iso_convergence::norm(delta);
        if(dn>radius*0.5 && dn>0.0) delta=delta*(radius*0.5/dn);
        bool accepted=false;
        double alpha=1.0;
        for(int ls=0;ls<16;++ls,alpha*=0.5) {
            const Point3 candidate=x-delta*alpha;
            if(iso_convergence::distance(candidate,center)>radius) continue;
            const auto trial=evaluateJet(ast,candidate,k);
            const double trialNorm=gradientNorm(trial.gradient);
            if(std::isfinite(trialNorm) && trial.smooth && trialNorm<gn) {
                x=candidate; accepted=true; break;
            }
        }
        // A small Newton step alone is not accepted as a stationary point: the
        // gradient equation must independently meet its own residual test.
        if(!accepted) return false;
    }
    const auto d=evaluateJet(ast,x,k);
    gradResidual=gradientNorm(d.gradient);
    root=x;
    const auto tolerances=criticalPointTolerances(ast,k,x,center,radius,d);
    return std::isfinite(gradResidual) && gradResidual<=tolerances.gradient;
}

inline std::vector<CriticalPoint> findCriticalPoints(const Node* ast,const std::vector<double>& kValues,
                                                      const Point3& probe,double radius) {
    std::vector<CriticalPoint> out;
    for(double k:kValues) {
        std::vector<CriticalPoint> candidates;

        // Always evaluate the selected probe exactly before Newton starts.
        // This prevents an approximate root from a neighbouring seed from
        // displacing an exact stationary probe during root deduplication.
        const auto probeJet=evaluateJet(ast,probe,k);
        const auto probeTol=criticalPointTolerances(ast,k,probe,probe,radius,probeJet);
        if(probeJet.smooth && gradientNorm(probeJet.gradient)<=probeTol.gradient) {
            addOrImproveCriticalPoint(candidates,
                makeCriticalPoint(ast,k,probe,probe,radius,0,true,"DIRECT_PROBE_EVALUATION"));
        }

        // Local 3x3x3 multi-start search. The exact center is intentionally
        // omitted because it was already evaluated above.
        for(int ix=-1;ix<=1;++ix) for(int iy=-1;iy<=1;++iy) for(int iz=-1;iz<=1;++iz) {
            if(ix==0 && iy==0 && iz==0) continue;
            const Point3 seed{probe.x+0.45*radius*ix,probe.y+0.45*radius*iy,probe.z+0.45*radius*iz};
            Point3 root; int iterations=0;
            double residual=std::numeric_limits<double>::infinity();
            if(!newtonStationary(ast,k,probe,radius,seed,root,iterations,residual)) continue;
            CriticalPoint candidate=makeCriticalPoint(ast,k,root,probe,radius,iterations,false,"MULTISTART_DAMPED_NEWTON");
            addOrImproveCriticalPoint(candidates,std::move(candidate));
        }

        if(candidates.empty()) {
            CriticalPoint row; row.k=k; row.status="NO_STATIONARY_POINT_FOUND_NOT_CERTIFIED"; out.push_back(row);
            continue;
        }
        std::sort(candidates.begin(),candidates.end(),[&](const CriticalPoint& a,const CriticalPoint& b) {
            return iso_convergence::distance(a.point,probe)<iso_convergence::distance(b.point,probe);
        });
        out.insert(out.end(),candidates.begin(),candidates.end());
    }
    return out;
}

struct Triangle { Point3 a,b,c,centroid,lo,hi; };
inline double component(const Point3& p,int axis) { return axis==0?p.x:(axis==1?p.y:p.z); }
inline Point3 minPoint(const Point3& a,const Point3& b) { return {std::min(a.x,b.x),std::min(a.y,b.y),std::min(a.z,b.z)}; }
inline Point3 maxPoint(const Point3& a,const Point3& b) { return {std::max(a.x,b.x),std::max(a.y,b.y),std::max(a.z,b.z)}; }
inline double pointBoxDistanceSquared(const Point3& p,const Point3& lo,const Point3& hi) {
    double s=0.0;
    for(int a=0;a<3;++a) { const double v=component(p,a), l=component(lo,a), h=component(hi,a);
        const double d=v<l?l-v:(v>h?v-h:0.0); s+=d*d; }
    return s;
}
inline double pointSegmentDistanceSquared(const Point3& p,const Point3& a,const Point3& b) {
    const Point3 ab=b-a; const double d=iso_convergence::normSquared(ab);
    const double t=d>0.0?std::clamp(iso_convergence::dot(p-a,ab)/d,0.0,1.0):0.0;
    return iso_convergence::normSquared(p-(a+ab*t));
}
inline double pointTriangleDistanceSquared(const Point3& p,const Triangle& t) {
    // Ericson's closest-point regions, with segment fallback for degenerate faces.
    const Point3 ab=t.b-t.a, ac=t.c-t.a, ap=p-t.a;
    const double d1=iso_convergence::dot(ab,ap), d2=iso_convergence::dot(ac,ap);
    if(d1<=0.0 && d2<=0.0) return iso_convergence::normSquared(ap);
    const Point3 bp=p-t.b; const double d3=iso_convergence::dot(ab,bp), d4=iso_convergence::dot(ac,bp);
    if(d3>=0.0 && d4<=d3) return iso_convergence::normSquared(bp);
    const double vc=d1*d4-d3*d2;
    if(vc<=0.0 && d1>=0.0 && d3<=0.0) { const double v=d1/(d1-d3); return iso_convergence::normSquared(p-(t.a+ab*v)); }
    const Point3 cp=p-t.c; const double d5=iso_convergence::dot(ab,cp), d6=iso_convergence::dot(ac,cp);
    if(d6>=0.0 && d5<=d6) return iso_convergence::normSquared(cp);
    const double vb=d5*d2-d1*d6;
    if(vb<=0.0 && d2>=0.0 && d6<=0.0) { const double w=d2/(d2-d6); return iso_convergence::normSquared(p-(t.a+ac*w)); }
    const double va=d3*d6-d5*d4;
    if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0) { const Point3 bc=t.c-t.b; const double w=(d4-d3)/((d4-d3)+(d5-d6)); return iso_convergence::normSquared(p-(t.b+bc*w)); }
    const double den=va+vb+vc;
    if(std::fabs(den)<1e-30) return std::min({pointSegmentDistanceSquared(p,t.a,t.b),pointSegmentDistanceSquared(p,t.b,t.c),pointSegmentDistanceSquared(p,t.c,t.a)});
    const double v=vb/den,w=vc/den;
    return iso_convergence::normSquared(p-(t.a+ab*v+ac*w));
}
inline std::vector<Triangle> meshTriangles(const IsosurfaceResult& mesh,const Point3& center,double radius) {
    std::vector<Triangle> out;
    out.reserve(mesh.vertices.size()/3);
    for(std::size_t i=0;i+2<mesh.vertices.size();i+=3) {
        auto p=[](const SurfaceVertex& v){return Point3{v.px,v.py,v.pz};};
        Triangle t; t.a=p(mesh.vertices[i]); t.b=p(mesh.vertices[i+1]); t.c=p(mesh.vertices[i+2]);
        t.centroid=(t.a+t.b+t.c)*(1.0/3.0);
        if(iso_convergence::distance(t.centroid,center)>radius) continue;
        t.lo=minPoint(t.a,minPoint(t.b,t.c)); t.hi=maxPoint(t.a,maxPoint(t.b,t.c)); out.push_back(t);
    }
    return out;
}
class TriangleBvh {
public:
    explicit TriangleBvh(std::vector<Triangle> triangles):triangles_(std::move(triangles)) {
        order_.resize(triangles_.size()); std::iota(order_.begin(),order_.end(),0);
        if(!order_.empty()) { nodes_.reserve(order_.size()*2); root_=build(0,order_.size()); }
    }
    bool empty() const { return triangles_.empty(); }
    double nearestDistance(const Point3& p) const {
        if(root_<0) return std::numeric_limits<double>::infinity();
        double best=std::numeric_limits<double>::infinity(); query(root_,p,best); return std::sqrt(best);
    }
private:
    struct BvhNode { Point3 lo,hi; int left=-1,right=-1; std::size_t begin=0,count=0; };
    std::vector<Triangle> triangles_; std::vector<std::size_t> order_; std::vector<BvhNode> nodes_; int root_=-1;
    int build(std::size_t begin,std::size_t end) {
        BvhNode n; n.begin=begin; n.count=end-begin;
        n.lo={std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
        n.hi={-n.lo.x,-n.lo.y,-n.lo.z};
        Point3 clo=n.lo, chi=n.hi;
        for(std::size_t i=begin;i<end;++i) { const auto& t=triangles_[order_[i]]; n.lo=minPoint(n.lo,t.lo); n.hi=maxPoint(n.hi,t.hi); clo=minPoint(clo,t.centroid); chi=maxPoint(chi,t.centroid); }
        const int id=static_cast<int>(nodes_.size()); nodes_.push_back(n);
        if(end-begin<=8) return id;
        const Point3 span=chi-clo; int axis=span.y>span.x?1:0; if(component(span,2)>component(span,axis)) axis=2;
        const std::size_t mid=begin+(end-begin)/2;
        std::nth_element(order_.begin()+static_cast<std::ptrdiff_t>(begin),order_.begin()+static_cast<std::ptrdiff_t>(mid),order_.begin()+static_cast<std::ptrdiff_t>(end),
            [&](std::size_t a,std::size_t b){return component(triangles_[a].centroid,axis)<component(triangles_[b].centroid,axis);});
        const int l=build(begin,mid), r=build(mid,end);
        nodes_[static_cast<std::size_t>(id)].left=l; nodes_[static_cast<std::size_t>(id)].right=r;
        nodes_[static_cast<std::size_t>(id)].count=0; return id;
    }
    void query(int id,const Point3& p,double& best) const {
        if(id<0) return;
        const auto& n=nodes_[static_cast<std::size_t>(id)];
        if(pointBoxDistanceSquared(p,n.lo,n.hi)>=best) return;
        if(n.left<0) { for(std::size_t i=n.begin;i<n.begin+n.count;++i) best=std::min(best,pointTriangleDistanceSquared(p,triangles_[order_[i]])); return; }
        const auto& a=nodes_[static_cast<std::size_t>(n.left)]; const auto& b=nodes_[static_cast<std::size_t>(n.right)];
        const double da=pointBoxDistanceSquared(p,a.lo,a.hi), db=pointBoxDistanceSquared(p,b.lo,b.hi);
        if(da<db) { query(n.left,p,best); query(n.right,p,best); } else { query(n.right,p,best); query(n.left,p,best); }
    }
};
inline std::vector<Point3> triangleInteriorSamples(const IsosurfaceResult& mesh,const Point3& center,double sampleRadius,int subdivisions=3) {
    std::vector<Point3> out;
    const int n=std::max(1,subdivisions);
    for(std::size_t i=0;i+2<mesh.vertices.size();i+=3) {
        const Point3 a{mesh.vertices[i].px,mesh.vertices[i].py,mesh.vertices[i].pz};
        const Point3 b{mesh.vertices[i+1].px,mesh.vertices[i+1].py,mesh.vertices[i+1].pz};
        const Point3 c{mesh.vertices[i+2].px,mesh.vertices[i+2].py,mesh.vertices[i+2].pz};
        const Point3 centroid=(a+b+c)*(1.0/3.0);
        if(iso_convergence::distance(centroid,center)>sampleRadius+1e-12) continue;
        for(int u=0;u<=n;++u) for(int v=0;v<=n-u;++v) {
            const double wa=static_cast<double>(u)/n, wb=static_cast<double>(v)/n, wc=1.0-wa-wb;
            const Point3 q=a*wa+b*wb+c*wc;
            if(iso_convergence::distance(q,center)<=sampleRadius) out.push_back(q);
        }
    }
    return out;
}
struct GeometrySample {
    double k=0.0;
    double meshToReference=std::numeric_limits<double>::quiet_NaN();
    double referenceToMesh=std::numeric_limits<double>::quiet_NaN();
    double symmetricEstimate=std::numeric_limits<double>::quiet_NaN();
    double meshToReferenceRms=std::numeric_limits<double>::quiet_NaN();
    double referenceToMeshRms=std::numeric_limits<double>::quiet_NaN();
    std::size_t sampleCountPerturbed=0, sampleCountReference=0;
    double interiorSamplingRadius=0.0;
    bool comparable=false;
};
inline std::pair<double,double> directedTriangleEstimate(const std::vector<Point3>& source,const TriangleBvh& target) {
    if(source.empty()||target.empty()) return {std::numeric_limits<double>::quiet_NaN(),std::numeric_limits<double>::quiet_NaN()};
    double maxd=0.0,sum=0.0;
    for(const auto& p:source) { const double d=target.nearestDistance(p); if(!std::isfinite(d)) return {NAN,NAN}; maxd=std::max(maxd,d); sum+=d*d; }
    return {maxd,std::sqrt(sum/static_cast<double>(source.size()))};
}
inline GeometrySample compareMeshes(double k,const IsosurfaceResult& mesh,const IsosurfaceResult& reference,
                                    const Point3& center,double radius,double gridSpacing,int subdivisions=3) {
    GeometrySample out; out.k=k;
    // Erode the requested ball by two grid widths. This avoids most of the
    // artificial mismatch caused by abruptly clipping triangle patches at the
    // ball boundary. The actual and interior radii are both recorded.
    out.interiorSamplingRadius=std::max(0.0,radius-2.0*gridSpacing);
    const auto meshTris=meshTriangles(mesh,center,radius+2.0*gridSpacing);
    const auto refTris=meshTriangles(reference,center,radius+2.0*gridSpacing);
    TriangleBvh meshBvh(meshTris), refBvh(refTris);
    const auto meshSamples=triangleInteriorSamples(mesh,center,out.interiorSamplingRadius,subdivisions);
    const auto refSamples=triangleInteriorSamples(reference,center,out.interiorSamplingRadius,subdivisions);
    out.sampleCountPerturbed=meshSamples.size(); out.sampleCountReference=refSamples.size();
    if(meshSamples.empty()||refSamples.empty()||meshBvh.empty()||refBvh.empty()) return out;
    const auto a=directedTriangleEstimate(meshSamples,refBvh);
    const auto b=directedTriangleEstimate(refSamples,meshBvh);
    out.meshToReference=a.first; out.meshToReferenceRms=a.second;
    out.referenceToMesh=b.first; out.referenceToMeshRms=b.second;
    if(std::isfinite(a.first)&&std::isfinite(b.first)) { out.symmetricEstimate=std::max(a.first,b.first); out.comparable=true; }
    return out;
}

struct TaylorRow { std::string name; DifferentialData data; int firstNonzeroOrder=-1; };
struct ResearchAnalysis {
    Point3 probe{};
    double requestedLocalRadius=0.0;
    std::string affineStatus="NOT_VERIFIED_AFFINE_IN_K";
    int gridResolution=0;
    double halfExtent=0.0, gridSpacing=0.0;
    int triangleSubdivisions=3;
    std::string fExpression="UNKNOWN",gExpression="UNKNOWN";
    TaylorRow fTaylor{"F",{}},gTaylor{"G",{}};
    IsolationAssessment isolation;
    std::vector<CriticalPoint> criticalPoints;
    std::vector<GeometrySample> geometry;
    std::string summary;
};
inline ResearchAnalysis analyze(const Node* familyAst,const iso_convergence::Options& options,
                                int triangleSubdivisions=3) {
    if(!familyAst) throw std::invalid_argument("Research analysis needs a parsed family expression.");
    if(options.kValues.empty()) throw std::invalid_argument("Research analysis needs explicit k values.");
    ResearchAnalysis out;
    out.probe=options.probe; out.requestedLocalRadius=options.localRadius;
    out.gridResolution=options.resolution; out.halfExtent=options.halfExtent;
    out.gridSpacing=2.0*options.halfExtent/options.resolution; out.triangleSubdivisions=triangleSubdivisions;
    const auto affine=iso_convergence::detail::decomposeAffineK(familyAst);
    if(!affine.valid) {
        out.summary="Taylor F/G analysis unavailable: the expression could not be structurally decomposed as F(x,y,z)+kG(x,y,z). The full k-dependent stationary search remains available.";
    } else {
        out.affineStatus="STRUCTURALLY_AFFINE_IN_K"; out.fExpression=affine.f; out.gExpression=affine.g;
        auto fAst=parseExpression(affine.f); auto gAst=parseExpression(affine.g);
        out.fTaylor.data=evaluateJet(fAst.get(),options.probe,0.0);
        out.fTaylor.firstNonzeroOrder=firstNonzeroJetOrder(out.fTaylor.data.taylor);
        out.gTaylor.data=evaluateJet(gAst.get(),options.probe,0.0);
        out.gTaylor.firstNonzeroOrder=firstNonzeroJetOrder(out.gTaylor.data.taylor);
        out.isolation=assessReferenceSingularity(fAst.get(),options.probe,options.localRadius);
    } 
    out.criticalPoints=findCriticalPoints(familyAst,options.kValues,options.probe,options.localRadius);
    const IsosurfaceResult reference=Isosurface::build(familyAst,0.0,options.resolution,options.halfExtent);
    out.geometry.reserve(options.kValues.size());
    const double h=2.0*options.halfExtent/options.resolution;
    for(double k:options.kValues) {
        IsosurfaceResult mesh=(k==0.0)?reference:Isosurface::build(familyAst,k,options.resolution,options.halfExtent);
        out.geometry.push_back(compareMeshes(k,mesh,reference,options.probe,options.localRadius,h,triangleSubdivisions));
    }
    std::ostringstream report;
    report << "Local isolation: " << out.isolation.status << ". ";
    if(out.affineStatus=="STRUCTURALLY_AFFINE_IN_K") {
        report << "Taylor jet orders at probe: first nonzero F order=" << out.fTaylor.firstNonzeroOrder
               << ", first nonzero G order=" << out.gTaylor.firstNonzeroOrder << ". ";
    }
    std::size_t singularCandidates=0, stationary=0;
    for(const auto& c:out.criticalPoints) { if(c.status.find("NUMERICAL_SINGULAR_CANDIDATE")!=std::string::npos) ++singularCandidates; if(c.status.find("NO_STATIONARY")!=0) ++stationary; }
    report << "Stationary-point search: " << stationary << " candidate rows, " << singularCandidates
           << " numerical on-zero-level singular candidates. Triangle-sampled geometry comparisons: "
           << out.geometry.size() << " k samples. All root-search labels and mesh distances are numerical diagnostics, not certified completeness or exact Hausdorff distances.";
    out.summary=report.str();
    return out;
}

inline std::string csvQuote(const std::string& s) {
    std::string r="\""; for(char c:s) { if(c=='\"') r+="\"\""; else r+=c; } return r+"\"";
}
inline std::string csvNumber(double v) { if(!std::isfinite(v)) return "NA"; std::ostringstream s; s<<std::setprecision(17)<<v; return s.str(); }
inline bool writeTaylorCsv(const std::string& path,const ResearchAnalysis& a,std::string* error=nullptr) {
    std::ofstream o(path); if(!o) { if(error)*error="Cannot open "+path; return false; }
    o<<"component,expression,first_nonzero_jet_order,smooth_jet,value,gx,gy,gz,hxx,hxy,hxz,hyy,hyz,hzz,"
       "hessian_determinant,eigenvalue_1,eigenvalue_2,eigenvalue_3,";
    static const char* names[20]={"c000","c100","c010","c001","c200","c110","c101","c020","c011","c002","c300","c210","c201","c120","c111","c102","c030","c021","c012","c003"};
    for(const char* n:names) o<<n<<',';
    o<<"affine_status,reference_isolation_status,probe_x,probe_y,probe_z\n";
    auto row=[&](const TaylorRow& t,const std::string& expr) {
        const auto& d=t.data; const auto eig=symmetricEigenvalues(d.hessian);
        o<<csvQuote(t.name)<<','<<csvQuote(expr)<<','<<t.firstNonzeroOrder<<','<<(d.smooth?1:0)<<','
         <<csvNumber(d.value)<<','<<csvNumber(d.gradient.x)<<','<<csvNumber(d.gradient.y)<<','<<csvNumber(d.gradient.z)<<','
         <<csvNumber(d.hessian[0][0])<<','<<csvNumber(d.hessian[0][1])<<','<<csvNumber(d.hessian[0][2])<<','
         <<csvNumber(d.hessian[1][1])<<','<<csvNumber(d.hessian[1][2])<<','<<csvNumber(d.hessian[2][2])<<','
         <<csvNumber(determinant3(d.hessian))<<','<<csvNumber(eig[0])<<','<<csvNumber(eig[1])<<','<<csvNumber(eig[2])<<',';
        for(double c:d.taylor)o<<csvNumber(c)<<',';
        o<<csvQuote(a.affineStatus)<<','<<csvQuote(a.isolation.status)<<','
         <<csvNumber(a.probe.x)<<','<<csvNumber(a.probe.y)<<','<<csvNumber(a.probe.z)<<"\n";
    };
    row(a.fTaylor,a.fExpression); row(a.gTaylor,a.gExpression); return static_cast<bool>(o);
}
inline bool writeCriticalPointsCsv(const std::string& path,const ResearchAnalysis& a,std::string* error=nullptr) {
    std::ofstream o(path); if(!o){if(error)*error="Cannot open "+path;return false;}
    o<<"k,critical_index,x,y,z,Fk_value,absolute_field_scale,field_zero_tolerance,gradient_norm,gradient_tolerance,"
       "newton_step_tolerance,deduplication_tolerance,hessian_determinant,newton_iterations,probe_point_candidate,"
       "search_source,status,search_method,reference_isolation_status,family_F,family_G\n";
    std::map<double,int> index;
    for(const auto& c:a.criticalPoints) {
        const int id=index[c.k]++;
        o<<csvNumber(c.k)<<','<<id<<','<<csvNumber(c.point.x)<<','<<csvNumber(c.point.y)<<','<<csvNumber(c.point.z)<<','
         <<csvNumber(c.fieldValue)<<','<<csvNumber(c.fieldScale)<<','<<csvNumber(c.fieldTolerance)<<','
         <<csvNumber(c.gradientMagnitude)<<','<<csvNumber(c.gradientTolerance)<<','
         <<csvNumber(c.newtonStepTolerance)<<','<<csvNumber(c.deduplicationTolerance)<<','
         <<csvNumber(c.hessianDeterminant)<<','<<c.iterations<<','<<(c.probePointCandidate?1:0)<<','
         <<csvQuote(c.searchSource)<<','<<csvQuote(c.status)<<','
         <<csvQuote("direct probe evaluation plus 26-start damped Newton on grad(F_k)=0; local and uncertified")<<','
         <<csvQuote(a.isolation.status)<<','<<csvQuote(a.fExpression)<<','<<csvQuote(a.gExpression)<<"\n";
    }
    return static_cast<bool>(o);
}
inline bool writeTriangleGeometryCsv(const std::string& path,const ResearchAnalysis& a,std::string* error=nullptr) {
    std::ofstream o(path); if(!o){if(error)*error="Cannot open "+path;return false;}
    o<<"k,mesh_to_reference_triangle_sampled_max,reference_to_mesh_triangle_sampled_max,symmetric_triangle_sampled_estimate,"
       "mesh_to_reference_rms,reference_to_mesh_rms,perturbed_source_samples,reference_source_samples,interior_sampling_radius,comparable,"
       "distance_method,requested_local_radius,grid_resolution_cells_per_axis,grid_spacing,domain_half_extent,triangle_subdivisions\n";
    for(const auto& g:a.geometry)
        o<<csvNumber(g.k)<<','<<csvNumber(g.meshToReference)<<','<<csvNumber(g.referenceToMesh)<<','<<csvNumber(g.symmetricEstimate)<<','
         <<csvNumber(g.meshToReferenceRms)<<','<<csvNumber(g.referenceToMeshRms)<<','<<g.sampleCountPerturbed<<','<<g.sampleCountReference<<','
         <<csvNumber(g.interiorSamplingRadius)<<','<<(g.comparable?1:0)<<','
         <<csvQuote("barycentric samples on source triangles; exact point-to-target-triangle distance accelerated by BVH; sampled estimate, not continuous Hausdorff")<<','
         <<csvNumber(a.requestedLocalRadius)<<','<<a.gridResolution<<','<<csvNumber(a.gridSpacing)<<','
         <<csvNumber(a.halfExtent)<<','<<a.triangleSubdivisions<<"\n";
    return static_cast<bool>(o);
}

} // namespace iso_research
