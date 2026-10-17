// Additive mesh-to-mesh interpolation experiment.
// This file does not modify or replace viewer functions. It consumes the same
// SurfaceCatalogue and Isosurface::build API as the application.
//
// For each parameterized surface and magnitude a > 0:
//   M- = mesh(F(x,-a)=0), M+ = mesh(F(x,+a)=0), M0 = mesh(F(x,0)=0).
// Since vertex indexing/connectivity can change across independently extracted
// meshes, it does NOT interpolate vertices by array index. Instead, each source
// mesh vertex is mapped to its closest point on the opposite endpoint mesh;
// halfway between them is the interpolated vertex. This is a deliberately
// explicit nearest-surface correspondence baseline. We do it in both directions
// and compare both interpolated meshes with the directly extracted M0 reference.
//
// The results are empirical diagnostics, not a universal proof of correspondence
// uniqueness. A correspondence can be ambiguous near a singularity or when the
// surfaces change topology; those ambiguities are part of what this test measures.

#include "Expr.h"
#include "Isosurface.h"
#include "SurfaceCatalogue.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cctype>
#include <cstdlib>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <sstream>
#include <string>
#include <stdexcept>
#include <utility>
#include <vector>

namespace fs = std::filesystem;

struct Vec3 {
    double x = 0.0, y = 0.0, z = 0.0;
};
static Vec3 operator+(const Vec3& a, const Vec3& b) { return {a.x+b.x,a.y+b.y,a.z+b.z}; }
static Vec3 operator-(const Vec3& a, const Vec3& b) { return {a.x-b.x,a.y-b.y,a.z-b.z}; }
static Vec3 operator*(const Vec3& a, double s) { return {a.x*s,a.y*s,a.z*s}; }
static Vec3 operator/(const Vec3& a, double s) { return {a.x/s,a.y/s,a.z/s}; }
static double dot(const Vec3& a, const Vec3& b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
static double norm2(const Vec3& a) { return dot(a,a); }
static double norm(const Vec3& a) { return std::sqrt(norm2(a)); }
static Vec3 minVec(const Vec3& a, const Vec3& b) { return {std::min(a.x,b.x),std::min(a.y,b.y),std::min(a.z,b.z)}; }
static Vec3 maxVec(const Vec3& a, const Vec3& b) { return {std::max(a.x,b.x),std::max(a.y,b.y),std::max(a.z,b.z)}; }
static double clamp01(double x) { return std::max(0.0, std::min(1.0, x)); }

struct TriangleMesh {
    std::vector<Vec3> vertices;
    std::vector<std::array<std::size_t,3>> faces;
};

struct VertexKey {
    std::int64_t x, y, z;
    bool operator<(const VertexKey& b) const {
        if (x != b.x) return x < b.x;
        if (y != b.y) return y < b.y;
        return z < b.z;
    }
};
static VertexKey makeKey(const Vec3& p, double tolerance) {
    return {static_cast<std::int64_t>(std::llround(p.x/tolerance)),
            static_cast<std::int64_t>(std::llround(p.y/tolerance)),
            static_cast<std::int64_t>(std::llround(p.z/tolerance))};
}
static Vec3 asVec(const SurfaceVertex& v) { return {v.px,v.py,v.pz}; }
static TriangleMesh weldTriangleSoup(const IsosurfaceResult& result, double tolerance) {
    TriangleMesh mesh;
    std::map<VertexKey,std::size_t> lookup;
    auto indexOf = [&](const SurfaceVertex& sv) {
        const Vec3 p=asVec(sv);
        const VertexKey key=makeKey(p,tolerance);
        const auto it=lookup.find(key);
        if(it!=lookup.end()) return it->second;
        const std::size_t id=mesh.vertices.size();
        mesh.vertices.push_back(p);
        lookup.emplace(key,id);
        return id;
    };
    for(std::size_t i=0;i+2<result.vertices.size();i+=3) {
        const std::size_t a=indexOf(result.vertices[i]);
        const std::size_t b=indexOf(result.vertices[i+1]);
        const std::size_t c=indexOf(result.vertices[i+2]);
        if(a==b || b==c || c==a) continue;
        mesh.faces.push_back({a,b,c});
    }
    return mesh;
}

static Vec3 closestPointOnTriangle(const Vec3& p, const Vec3& a, const Vec3& b, const Vec3& c) {
    // Real-Time Collision Detection, Christer Ericson, section 5.1.5.
    const Vec3 ab=b-a, ac=c-a, ap=p-a;
    const double d1=dot(ab,ap), d2=dot(ac,ap);
    if(d1<=0.0 && d2<=0.0) return a;
    const Vec3 bp=p-b;
    const double d3=dot(ab,bp), d4=dot(ac,bp);
    if(d3>=0.0 && d4<=d3) return b;
    const double vc=d1*d4-d3*d2;
    if(vc<=0.0 && d1>=0.0 && d3<=0.0) {
        const double v=d1/(d1-d3); return a+ab*v;
    }
    const Vec3 cp=p-c;
    const double d5=dot(ab,cp), d6=dot(ac,cp);
    if(d6>=0.0 && d5<=d6) return c;
    const double vb=d5*d2-d1*d6;
    if(vb<=0.0 && d2>=0.0 && d6<=0.0) {
        const double w=d2/(d2-d6); return a+ac*w;
    }
    const double va=d3*d6-d5*d4;
    if(va<=0.0 && (d4-d3)>=0.0 && (d5-d6)>=0.0) {
        const Vec3 bc=c-b;
        const double w=(d4-d3)/((d4-d3)+(d5-d6)); return b+bc*w;
    }
    const double denom=va+vb+vc;
    if(std::fabs(denom)<1e-30) {
        // Degenerate fallback: closest point on each segment.
        auto segment=[&](const Vec3& p0,const Vec3& p1) {
            const Vec3 d=p1-p0; const double l=norm2(d);
            const double t=l>0.0?clamp01(dot(p-p0,d)/l):0.0; return p0+d*t;
        };
        const Vec3 q0=segment(a,b),q1=segment(b,c),q2=segment(c,a);
        const double e0=norm2(p-q0),e1=norm2(p-q1),e2=norm2(p-q2);
        return e0<=e1 && e0<=e2?q0:(e1<=e2?q1:q2);
    }
    const double inv=1.0/denom;
    const double v=vb*inv, w=vc*inv;
    return a+ab*v+ac*w;
}

struct TriangleRef {
    Vec3 a,b,c,lo,hi,centroid;
    std::size_t faceIndex=0;
};
struct BvhNode {
    Vec3 lo,hi;
    std::size_t begin=0,count=0;
    int left=-1,right=-1;
};
static double pointAabbDistance2(const Vec3& p,const Vec3& lo,const Vec3& hi) {
    auto axis=[](double x,double a,double b) { return x<a?a-x:(x>b?x-b:0.0); };
    const double dx=axis(p.x,lo.x,hi.x),dy=axis(p.y,lo.y,hi.y),dz=axis(p.z,lo.z,hi.z);
    return dx*dx+dy*dy+dz*dz;
}
class TriangleBvh {
public:
    explicit TriangleBvh(const TriangleMesh& mesh) {
        tris_.reserve(mesh.faces.size());
        for(std::size_t i=0;i<mesh.faces.size();++i) {
            const auto f=mesh.faces[i];
            const Vec3 a=mesh.vertices[f[0]],b=mesh.vertices[f[1]],c=mesh.vertices[f[2]];
            const Vec3 cross={ (b.y-a.y)*(c.z-a.z)-(b.z-a.z)*(c.y-a.y),
                               (b.z-a.z)*(c.x-a.x)-(b.x-a.x)*(c.z-a.z),
                               (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x) };
            if(norm2(cross)<1e-28) continue;
            TriangleRef t; t.a=a;t.b=b;t.c=c;t.lo=minVec(a,minVec(b,c));t.hi=maxVec(a,maxVec(b,c));
            t.centroid=(a+b+c)/3.0;t.faceIndex=i; tris_.push_back(t);
        }
        order_.resize(tris_.size()); std::iota(order_.begin(),order_.end(),std::size_t{0});
        if(!order_.empty()) buildNode(0,order_.size());
    }
    bool empty() const { return order_.empty(); }
    Vec3 closestPoint(const Vec3& p) const {
        if(empty()) return p;
        double best=std::numeric_limits<double>::infinity(); Vec3 point=p;
        query(0,p,best,point); return point;
    }
    double distance(const Vec3& p) const { return norm(p-closestPoint(p)); }
private:
    std::vector<TriangleRef> tris_;
    std::vector<std::size_t> order_;
    std::vector<BvhNode> nodes_;
    static constexpr std::size_t leafSize_=8;
    int buildNode(std::size_t begin,std::size_t end) {
        BvhNode node;
        node.lo={std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity(),std::numeric_limits<double>::infinity()};
        node.hi={-node.lo.x,-node.lo.y,-node.lo.z};
        Vec3 clo=node.lo,chi=node.hi;
        for(std::size_t i=begin;i<end;++i) {
            const auto& t=tris_[order_[i]];
            node.lo=minVec(node.lo,t.lo);node.hi=maxVec(node.hi,t.hi);
            clo=minVec(clo,t.centroid);chi=maxVec(chi,t.centroid);
        }
        const int id=static_cast<int>(nodes_.size()); nodes_.push_back(node);
        const std::size_t count=end-begin;
        if(count<=leafSize_) { nodes_[id].begin=begin;nodes_[id].count=count; return id; }
        const Vec3 extent=chi-clo;
        int axis=(extent.y>extent.x)?1:0;
        if((axis==0?extent.x:extent.y)<extent.z) axis=2;
        const std::size_t mid=begin+count/2;
        auto component=[&](const TriangleRef& t)->double { return axis==0?t.centroid.x:(axis==1?t.centroid.y:t.centroid.z); };
        std::nth_element(order_.begin()+static_cast<std::ptrdiff_t>(begin),order_.begin()+static_cast<std::ptrdiff_t>(mid),order_.begin()+static_cast<std::ptrdiff_t>(end),
                         [&](std::size_t a,std::size_t b){return component(tris_[a])<component(tris_[b]);});
        const int left=buildNode(begin,mid),right=buildNode(mid,end);
        nodes_[id].left=left;nodes_[id].right=right;
        return id;
    }
    void query(int id,const Vec3& p,double& best,Vec3& point) const {
        const BvhNode& n=nodes_[static_cast<std::size_t>(id)];
        if(pointAabbDistance2(p,n.lo,n.hi)>best) return;
        if(n.left<0 && n.right<0) {
            for(std::size_t i=n.begin;i<n.begin+n.count;++i) {
                const auto& t=tris_[order_[i]];
                const Vec3 q=closestPointOnTriangle(p,t.a,t.b,t.c); const double d=norm2(p-q);
                if(d<best) {best=d;point=q;}
            }
            return;
        }
        const double dl=pointAabbDistance2(p,nodes_[static_cast<std::size_t>(n.left)].lo,nodes_[static_cast<std::size_t>(n.left)].hi);
        const double dr=pointAabbDistance2(p,nodes_[static_cast<std::size_t>(n.right)].lo,nodes_[static_cast<std::size_t>(n.right)].hi);
        if(dl<dr) {query(n.left,p,best,point);query(n.right,p,best,point);}
        else {query(n.right,p,best,point);query(n.left,p,best,point);}
    }
};

static TriangleMesh interpolateByClosestSurface(const TriangleMesh& source,const TriangleBvh& target) {
    TriangleMesh out; out.faces=source.faces; out.vertices.reserve(source.vertices.size());
    for(const Vec3& p:source.vertices) out.vertices.push_back((p+target.closestPoint(p))*0.5);
    return out;
}

struct TopologyMetrics {
    std::size_t vertices=0,edges=0,faces=0,boundaryEdges=0,nonmanifoldEdges=0,orientationConflicts=0,degenerateFaces=0;
    long long euler=0;
};
struct EdgeInfo { unsigned count=0; int orientationBalance=0; };
static TopologyMetrics topologyMetrics(const TriangleMesh& mesh,double areaEpsilon) {
    TopologyMetrics m; m.vertices=mesh.vertices.size();m.faces=mesh.faces.size();
    std::map<std::pair<std::size_t,std::size_t>,EdgeInfo> edges;
    for(const auto& f:mesh.faces) {
        const Vec3 a=mesh.vertices[f[0]],b=mesh.vertices[f[1]],c=mesh.vertices[f[2]];
        const Vec3 cr={ (b.y-a.y)*(c.z-a.z)-(b.z-a.z)*(c.y-a.y),
                        (b.z-a.z)*(c.x-a.x)-(b.x-a.x)*(c.z-a.z),
                        (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x) };
        if(0.5*norm(cr)<=areaEpsilon) ++m.degenerateFaces;
        for(int e=0;e<3;++e) {
            const std::size_t aId=f[static_cast<std::size_t>(e)],bId=f[static_cast<std::size_t>((e+1)%3)];
            const auto key=std::minmax(aId,bId); auto& info=edges[{key.first,key.second}];
            ++info.count; info.orientationBalance += (aId<bId)?1:-1;
        }
    }
    m.edges=edges.size();
    for(const auto& kv:edges) {
        const auto& e=kv.second;
        if(e.count==1) ++m.boundaryEdges;
        else if(e.count!=2) ++m.nonmanifoldEdges;
        else if(e.orientationBalance!=0) ++m.orientationConflicts;
    }
    m.euler=static_cast<long long>(m.vertices)-static_cast<long long>(m.edges)+static_cast<long long>(m.faces);
    return m;
}
struct DistanceStats { double mean=0,p95=0,max=0;std::size_t count=0; };
static DistanceStats pointToSurfaceStats(const std::vector<Vec3>& points,const TriangleBvh& target) {
    DistanceStats s; if(points.empty()||target.empty()) return s;
    std::vector<double> d;d.reserve(points.size()); double sum=0.0;
    for(const auto& p:points) { const double q=target.distance(p);d.push_back(q);sum+=q;s.max=std::max(s.max,q); }
    std::sort(d.begin(),d.end());s.count=d.size();s.mean=sum/static_cast<double>(d.size());
    const std::size_t i95=static_cast<std::size_t>(std::ceil(0.95*static_cast<double>(d.size())))-1;
    s.p95=d[std::min(i95,d.size()-1)];return s;
}
struct ScalarStats { double mean=0,p95=0,max=0;std::size_t count=0; };
static ScalarStats stats(std::vector<double> values) {
    ScalarStats s;if(values.empty())return s;
    double sum=0;for(double v:values){sum+=v;s.max=std::max(s.max,v);}std::sort(values.begin(),values.end());
    s.count=values.size();s.mean=sum/static_cast<double>(values.size());
    const std::size_t i95=static_cast<std::size_t>(std::ceil(0.95*static_cast<double>(values.size())))-1;
    s.p95=values[std::min(i95,values.size()-1)];return s;
}
struct NormalStats { double meanDeg=std::numeric_limits<double>::quiet_NaN(),maxDeg=std::numeric_limits<double>::quiet_NaN();std::size_t count=0; };
static Dual3 evalDualAt(const Node* ast,const Vec3& p) {
    return evalExpr<Dual3>(ast,Dual3(p.x,1,0,0),Dual3(p.y,0,1,0),Dual3(p.z,0,0,1),Dual3(0.0,0,0,0));
}
static NormalStats compareNormalsToReference(const TriangleMesh& candidate,const TriangleBvh& reference,const Node* ast) {
    NormalStats out;std::vector<double> angles;
    for(const Vec3& p:candidate.vertices) {
        const Vec3 q=reference.closestPoint(p);
        const Dual3 gp=evalDualAt(ast,p),gq=evalDualAt(ast,q);
        const double mp=std::sqrt(gp.dx*gp.dx+gp.dy*gp.dy+gp.dz*gp.dz);
        const double mq=std::sqrt(gq.dx*gq.dx+gq.dy*gq.dy+gq.dz*gq.dz);
        // Normals are undefined at singular/near-singular points. Omit those pairs.
        if(mp<1e-7||mq<1e-7||!std::isfinite(mp)||!std::isfinite(mq))continue;
        const double d=(gp.dx*gq.dx+gp.dy*gq.dy+gp.dz*gq.dz)/(mp*mq);
        const double angle=std::acos(std::max(-1.0,std::min(1.0,d)))*180.0/3.14159265358979323846;
        if(std::isfinite(angle))angles.push_back(angle);
    }
    const auto s=stats(angles);out.count=s.count;
    if(s.count>0) { out.meanDeg=s.mean;out.maxDeg=s.max; }
    return out;
}
static std::string csvQuote(const std::string& s) {
    std::string out="\"";for(char c:s){if(c=='\"')out+="\"\"";else out+=c;}return out+"\"";
}
static std::string safeId(std::string s) {
    for(char& c:s)if(!(std::isalnum(static_cast<unsigned char>(c))||c=='_'||c=='-'))c='_';return s;
}
static std::string kTag(double k) {
    std::ostringstream s;s<<std::fixed<<std::setprecision(4)<<k;std::string r=s.str();
    for(char& c:r){if(c=='-')c='m';else if(c=='.')c='p';}return r;
}
static bool writePly(const fs::path& path,const TriangleMesh& mesh) {
    std::ofstream f(path);if(!f)return false;
    f<<"ply\nformat ascii 1.0\nelement vertex "<<mesh.vertices.size()<<"\nproperty double x\nproperty double y\nproperty double z\nelement face "<<mesh.faces.size()<<"\nproperty list uchar int vertex_indices\nend_header\n";
    f<<std::setprecision(10);
    for(const auto& p:mesh.vertices)f<<p.x<<' '<<p.y<<' '<<p.z<<'\n';
    for(const auto& t:mesh.faces)f<<"3 "<<t[0]<<' '<<t[1]<<' '<<t[2]<<'\n';
    return static_cast<bool>(f);
}
static std::string behaviorName(iso_surfaces::Behavior b) {
    switch(b) {case iso_surfaces::Behavior::RegularForNonzeroK:return "regular_for_nonzero_k_expected";
        case iso_surfaces::Behavior::WhitneyUmbrella:return "whitney_umbrella_sign_dependent";
        case iso_surfaces::Behavior::PersistentOrigin:return "singularity_persists";}
    return "unknown";
}
static std::vector<std::string> splitStrings(const std::string& text) {
    std::vector<std::string> out;std::stringstream ss(text);std::string item;
    while(std::getline(ss,item,','))if(!item.empty())out.push_back(item);return out;
}
static std::vector<double> splitDoubles(const std::string& text) {
    std::vector<double> out;for(const auto& item:splitStrings(text)){const double v=std::stod(item);if(!(v>0.0)||!std::isfinite(v))throw std::runtime_error("k magnitudes must be positive finite numbers");out.push_back(v);}return out;
}
struct Options {
    int resolution=32;
    double halfExtent=1.0;
    fs::path out="mesh_to_mesh_interpolation_results";
    std::vector<double> magnitudes={0.1,0.01,0.001,0.0001};
    std::vector<std::string> surfaces={"plane_x3_ky","isolated_cusp_x","mixed_cubic","level_shift_cone","monkey_saddle","crossing_planes_xy"};
};
static Options parseOptions(int argc,char** argv) {
    Options o;
    for(int i=1;i<argc;++i){const std::string a=argv[i];
        auto next=[&]()->std::string{if(i+1>=argc)throw std::runtime_error("missing value after "+a);return argv[++i];};
        if(a=="--resolution")o.resolution=std::stoi(next());
        else if(a=="--domain-half-extent")o.halfExtent=std::stod(next());
        else if(a=="--out")o.out=next();
        else if(a=="--k-values")o.magnitudes=splitDoubles(next());
        else if(a=="--surfaces")o.surfaces=splitStrings(next());
        else if(a=="--help"||a=="-h"){
            std::cout<<"Usage: isosurface_mesh_interpolation_tests [--resolution N] [--domain-half-extent H] [--k-values 0.01,0.001] [--surfaces id,id] [--out DIR]\n";
            std::cout<<"Default magnitudes: 0.1,0.01,0.001,0.0001. Every magnitude tests both -k and +k against the direct k=0 mesh.\n";std::exit(0);
        }else throw std::runtime_error("unknown option: "+a);
    }
    if(o.resolution<8||o.resolution>256)throw std::runtime_error("resolution must be in [8,256]");
    if(!(o.halfExtent>0.0)||!std::isfinite(o.halfExtent))throw std::runtime_error("domain half extent must be positive finite");
    if(o.magnitudes.empty())throw std::runtime_error("provide at least one k magnitude");
    return o;
}

static void writeCsvHeader(std::ofstream& f) {
    f<<"surface_id,label,expression,analytic_expectation,k_magnitude,case_status,resolution,domain_half_extent,negative_triangles,positive_triangles,direct_zero_triangles,mid_from_negative_triangles,mid_from_positive_triangles,negative_skipped_degenerate,positive_skipped_degenerate,zero_skipped_degenerate,mid_from_negative_mean_abs_F0,mid_from_negative_p95_abs_F0,mid_from_negative_max_abs_F0,mid_from_positive_mean_abs_F0,mid_from_positive_p95_abs_F0,mid_from_positive_max_abs_F0,negative_to_zero_mean_dist,negative_to_zero_p95_dist,negative_to_zero_max_dist,zero_to_negative_mean_dist,zero_to_negative_p95_dist,zero_to_negative_max_dist,negative_symmetric_hausdorff_est,positive_to_zero_mean_dist,positive_to_zero_p95_dist,positive_to_zero_max_dist,zero_to_positive_mean_dist,zero_to_positive_p95_dist,zero_to_positive_max_dist,positive_symmetric_hausdorff_est,neg_mid_to_zero_mean_dist,neg_mid_to_zero_p95_dist,neg_mid_to_zero_max_dist,zero_to_neg_mid_mean_dist,zero_to_neg_mid_p95_dist,zero_to_neg_mid_max_dist,neg_mid_symmetric_hausdorff_est,pos_mid_to_zero_mean_dist,pos_mid_to_zero_p95_dist,pos_mid_to_zero_max_dist,zero_to_pos_mid_mean_dist,zero_to_pos_mid_p95_dist,zero_to_pos_mid_max_dist,pos_mid_symmetric_hausdorff_est,neg_mid_vs_negative_endpoint_error_ratio,pos_mid_vs_positive_endpoint_error_ratio,neg_mid_normal_pairs,neg_mid_mean_normal_angle_deg,neg_mid_max_normal_angle_deg,pos_mid_normal_pairs,pos_mid_mean_normal_angle_deg,pos_mid_max_normal_angle_deg,zero_boundary_edges,zero_nonmanifold_edges,zero_orientation_conflicts,negative_boundary_edges,negative_nonmanifold_edges,positive_boundary_edges,positive_nonmanifold_edges,neg_mid_boundary_edges,neg_mid_nonmanifold_edges,neg_mid_degenerate_faces,pos_mid_boundary_edges,pos_mid_nonmanifold_edges,pos_mid_degenerate_faces,neg_mid_ply,pos_mid_ply,direct_zero_ply\n";
}
static void writeDistanceFields(std::ofstream& f,const DistanceStats& a,const DistanceStats& b) {
    f<<a.mean<<','<<a.p95<<','<<a.max<<','<<b.mean<<','<<b.p95<<','<<b.max<<','<<std::max(a.max,b.max)<<',';
}
static bool findSurfaceById(const std::string& id,const iso_surfaces::SurfaceCase*& result) {
    for(const auto& s:iso_surfaces::catalogue())if(s.id==id){result=&s;return true;}return false;
}
static void writeUnavailableRow(std::ofstream& csv,const iso_surfaces::SurfaceCase& surface,
                                double magnitude,const std::string& status,int resolution,double halfExtent) {
    // CSV schema has 58 columns: fill stable key/configuration fields and leave
    // unavailable geometric measurements blank rather than silently dropping a case.
    csv<<csvQuote(surface.id)<<','<<csvQuote(surface.label)<<','<<csvQuote(surface.expression)<<','
       <<csvQuote(behaviorName(surface.behavior))<<','<<magnitude<<','<<status<<','<<resolution<<','<<halfExtent;
    for(int column=8;column<74;++column) csv<<',';
    csv<<'\n';
}

int main(int argc,char** argv) {
    try {
        const Options options=parseOptions(argc,argv);
        fs::create_directories(options.out);
        const fs::path csvPath=options.out/"mesh_to_mesh_interpolation.csv";
        std::ofstream csv(csvPath);
        if(!csv)throw std::runtime_error("cannot open output CSV: "+csvPath.string());
        csv<<std::setprecision(12);writeCsvHeader(csv);
        const double weldTolerance=std::max(1e-7,options.halfExtent*1e-7);
        const double areaEpsilon=std::pow((2.0*options.halfExtent/options.resolution),2)*1e-10;
        std::size_t successfulCases=0,unavailableCases=0,failedCases=0;
        std::cout<<"Mesh-to-mesh interpolation test: resolution="<<options.resolution
                 <<", domain=["<<-options.halfExtent<<','<<options.halfExtent<<"]^3, k magnitudes=";
        for(std::size_t i=0;i<options.magnitudes.size();++i){if(i)std::cout<<',';std::cout<<options.magnitudes[i];}
        std::cout<<"\nCorrespondence: each source vertex -> closest point on opposite endpoint mesh; interpolate halfway in both directions.\n";
        for(const std::string& id:options.surfaces) {
            const iso_surfaces::SurfaceCase* surface=nullptr;
            if(!findSurfaceById(id,surface)){std::cerr<<"ERROR: unknown catalogue surface ID: "<<id<<"\n";++failedCases;continue;}
            NodePtr ast;
            try{ast=parseExpression(surface->expression);}catch(const std::exception& e){std::cerr<<"ERROR parsing "<<id<<": "<<e.what()<<"\n";++failedCases;continue;}
            const IsosurfaceResult rawZero=Isosurface::build(ast.get(),0.0,options.resolution,options.halfExtent);
            const TriangleMesh zero=weldTriangleSoup(rawZero,weldTolerance);
            const TriangleBvh zeroBvh(zero);
            const TopologyMetrics zeroTop=topologyMetrics(zero,areaEpsilon);
            const fs::path zeroPly=options.out/(safeId(id)+"_direct_k0.ply");
            const bool zeroPlyOk=writePly(zeroPly,zero);
            std::cout<<"\n["<<id<<"] "<<surface->label<<"\n  direct k=0: triangles="<<zero.faces.size()<<", skipped="<<rawZero.skippedDegenerateTriangles
                     <<", boundary="<<zeroTop.boundaryEdges<<", nonmanifold="<<zeroTop.nonmanifoldEdges<<", Euler="<<zeroTop.euler<<"\n";
            for(const double magnitude:options.magnitudes) {
                const IsosurfaceResult rawNeg=Isosurface::build(ast.get(),-magnitude,options.resolution,options.halfExtent);
                const IsosurfaceResult rawPos=Isosurface::build(ast.get(),+magnitude,options.resolution,options.halfExtent);
                const TriangleMesh neg=weldTriangleSoup(rawNeg,weldTolerance),pos=weldTriangleSoup(rawPos,weldTolerance);
                std::string unavailableStatus;
                if(zero.faces.empty() || zeroBvh.empty()) unavailableStatus="EMPTY_DIRECT_ZERO_MESH";
                else if(neg.faces.empty()) unavailableStatus="EMPTY_NEGATIVE_ENDPOINT_MESH";
                else if(pos.faces.empty()) unavailableStatus="EMPTY_POSITIVE_ENDPOINT_MESH";
                if(!unavailableStatus.empty()) {
                    std::cout<<"  k=±"<<magnitude<<": "<<unavailableStatus<<" (M-="<<neg.faces.size()<<", M+="<<pos.faces.size()<<", M0="<<zero.faces.size()<<")\n";
                    writeUnavailableRow(csv,*surface,magnitude,unavailableStatus,options.resolution,options.halfExtent);
                    ++unavailableCases; continue;
                }
                const TriangleBvh negBvh(neg),posBvh(pos);
                if(negBvh.empty() || posBvh.empty()) {
                    writeUnavailableRow(csv,*surface,magnitude,"EMPTY_ENDPOINT_TRIANGLE_BVH",options.resolution,options.halfExtent);
                    ++unavailableCases; continue;
                }
                const TriangleMesh midNeg=interpolateByClosestSurface(neg,posBvh);
                const TriangleMesh midPos=interpolateByClosestSurface(pos,negBvh);
                const TriangleBvh midNegBvh(midNeg),midPosBvh(midPos);
                if(midNegBvh.empty() || midPosBvh.empty()) {
                    std::cout<<"  k=±"<<magnitude<<": midpoint collapsed to no non-degenerate triangles\n";
                    writeUnavailableRow(csv,*surface,magnitude,"COLLAPSED_MIDPOINT_MESH",options.resolution,options.halfExtent);
                    ++unavailableCases; continue;
                }
                const DistanceStats negEndpointToZero=pointToSurfaceStats(neg.vertices,zeroBvh);
                const DistanceStats zeroToNegEndpoint=pointToSurfaceStats(zero.vertices,negBvh);
                const DistanceStats posEndpointToZero=pointToSurfaceStats(pos.vertices,zeroBvh);
                const DistanceStats zeroToPosEndpoint=pointToSurfaceStats(zero.vertices,posBvh);
                const DistanceStats negToZero=pointToSurfaceStats(midNeg.vertices,zeroBvh);
                const DistanceStats zeroToNeg=pointToSurfaceStats(zero.vertices,midNegBvh);
                const DistanceStats posToZero=pointToSurfaceStats(midPos.vertices,zeroBvh);
                const DistanceStats zeroToPos=pointToSurfaceStats(zero.vertices,midPosBvh);
                const double negEndpointSym=std::max(negEndpointToZero.max,zeroToNegEndpoint.max);
                const double posEndpointSym=std::max(posEndpointToZero.max,zeroToPosEndpoint.max);
                const double negMidSym=std::max(negToZero.max,zeroToNeg.max);
                const double posMidSym=std::max(posToZero.max,zeroToPos.max);
                const double negImprovementRatio=negEndpointSym>1e-15?negMidSym/negEndpointSym:std::numeric_limits<double>::quiet_NaN();
                const double posImprovementRatio=posEndpointSym>1e-15?posMidSym/posEndpointSym:std::numeric_limits<double>::quiet_NaN();
                const NormalStats negNormals=compareNormalsToReference(midNeg,zeroBvh,ast.get());
                const NormalStats posNormals=compareNormalsToReference(midPos,zeroBvh,ast.get());
                auto fieldResidual=[&](const TriangleMesh& m) {
                    std::vector<double> values;values.reserve(m.vertices.size());
                    for(const Vec3& p:m.vertices) values.push_back(std::fabs(evalExpr<double>(ast.get(),p.x,p.y,p.z,0.0)));
                    return stats(std::move(values));
                };
                const ScalarStats negF=fieldResidual(midNeg),posF=fieldResidual(midPos);
                const auto negTop=topologyMetrics(neg,areaEpsilon),posTop=topologyMetrics(pos,areaEpsilon);
                const auto midNegTop=topologyMetrics(midNeg,areaEpsilon),midPosTop=topologyMetrics(midPos,areaEpsilon);
                const auto directTopology=zeroTop;
                const std::string tag=kTag(magnitude);
                const fs::path negPly=options.out/(safeId(id)+"_k"+tag+"_mid_from_negative.ply");
                const fs::path posPly=options.out/(safeId(id)+"_k"+tag+"_mid_from_positive.ply");
                const bool plysOk=writePly(negPly,midNeg)&&writePly(posPly,midPos);
                csv<<csvQuote(surface->id)<<','<<csvQuote(surface->label)<<','<<csvQuote(surface->expression)<<','<<csvQuote(behaviorName(surface->behavior))<<','
                   <<magnitude<<",OK,"<<options.resolution<<','<<options.halfExtent<<','
                   <<neg.faces.size()<<','<<pos.faces.size()<<','<<zero.faces.size()<<','<<midNeg.faces.size()<<','<<midPos.faces.size()<<','
                   <<rawNeg.skippedDegenerateTriangles<<','<<rawPos.skippedDegenerateTriangles<<','<<rawZero.skippedDegenerateTriangles<<','
                   <<negF.mean<<','<<negF.p95<<','<<negF.max<<','<<posF.mean<<','<<posF.p95<<','<<posF.max<<',';
                writeDistanceFields(csv,negEndpointToZero,zeroToNegEndpoint);
                writeDistanceFields(csv,posEndpointToZero,zeroToPosEndpoint);
                writeDistanceFields(csv,negToZero,zeroToNeg);
                writeDistanceFields(csv,posToZero,zeroToPos);
                csv<<negImprovementRatio<<','<<posImprovementRatio<<','
                   <<negNormals.count<<','<<negNormals.meanDeg<<','<<negNormals.maxDeg<<','<<posNormals.count<<','<<posNormals.meanDeg<<','<<posNormals.maxDeg<<','
                   <<directTopology.boundaryEdges<<','<<directTopology.nonmanifoldEdges<<','<<directTopology.orientationConflicts<<','
                   <<negTop.boundaryEdges<<','<<negTop.nonmanifoldEdges<<','<<posTop.boundaryEdges<<','<<posTop.nonmanifoldEdges<<','
                   <<midNegTop.boundaryEdges<<','<<midNegTop.nonmanifoldEdges<<','<<midNegTop.degenerateFaces<<','
                   <<midPosTop.boundaryEdges<<','<<midPosTop.nonmanifoldEdges<<','<<midPosTop.degenerateFaces<<','
                   <<csvQuote(negPly.string())<<','<<csvQuote(posPly.string())<<','<<csvQuote(zeroPly.string())<<'\n';
                std::cout<<"  k=±"<<magnitude<<": M- triangles="<<neg.faces.size()<<", M+ triangles="<<pos.faces.size()
                         <<", M0 triangles="<<zero.faces.size()<<"\n"
                         <<"    endpoint-to-M0 max distance: M-="<<negEndpointSym<<", M+="<<posEndpointSym<<"\n"
                         <<"    midpoint-from-M-: symmetric distance estimate="<<negMidSym<<", endpoint error ratio="<<negImprovementRatio
                         <<", mean|F0|="<<negF.mean<<", p95|F0|="<<negF.p95<<", max|F0|="<<negF.max
                         <<", normals compared="<<negNormals.count<<", mean angle="<<negNormals.meanDeg<<" deg\n"
                         <<"    midpoint-from-M+: symmetric distance estimate="<<posMidSym<<", endpoint error ratio="<<posImprovementRatio
                         <<", mean|F0|="<<posF.mean<<", p95|F0|="<<posF.p95<<", max|F0|="<<posF.max
                         <<", normals compared="<<posNormals.count<<", mean angle="<<posNormals.meanDeg<<" deg\n";
                if(!std::isfinite(negF.max)||!std::isfinite(posF.max)||!std::isfinite(negToZero.max)||!std::isfinite(posToZero.max)||!plysOk||!zeroPlyOk) {
                    std::cerr<<"  ERROR: non-finite measurements or PLY write failure for "<<id<<" k="<<magnitude<<"\n";++failedCases;
                } else ++successfulCases;
            }
        }
        csv.close();
        std::cout<<"\nWrote: "<<csvPath.string()<<"\n"
                 <<"Wrote direct k=0 reference meshes and midpoint PLY files to: "<<options.out.string()<<"\n"
                 <<"Completed cases: "<<successfulCases<<", unavailable geometry cases: "<<unavailableCases<<", fatal errors: "<<failedCases<<"\n"
                 <<"Interpretation: smaller midpoint-to-M0 distances and |F(x,0)| are better. Compare both directions and resolutions.\n"
                 <<"Normal angles omit points where either gradient magnitude is < 1e-7 because singular-point normals are undefined.\n";
        return failedCases==0?0:1;
    } catch(const std::exception& e) {
        std::cerr<<"ERROR: "<<e.what()<<"\n";return 2;
    }
}
