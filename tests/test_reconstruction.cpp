// Tests of src/analysis/Reconstruction.h: the weighted-cone patch on the ten
// surfaces of the thesis report, the field average, and the viewer entry point.
#include "analysis/Reconstruction.h"

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

int failures = 0;

void check(bool ok, const std::string& what) {
    if (!ok) {
        ++failures;
        std::printf("FAIL: %s\n", what.c_str());
    }
}

struct Case { const char* name; const char* f; int loops; };

// Largest distance from a mesh vertex in the box to {F = 0} (Newton projection).
double vertexError(const iso_reconstruct::Mesh& m, const iso_exact::Field& f, double box, double h) {
    double worst = 0.0;
    std::vector<char> used(m.v.size(), 0);
    for (const auto& t : m.t) used[t[0]] = used[t[1]] = used[t[2]] = 1;
    for (std::size_t i = 0; i < m.v.size(); ++i) {
        const auto& v = m.v[i];
        if (!used[i] || iso_reconstruct::infDistance(v, {0, 0, 0}) > box) continue;
        iso_geometry::Point3 p = v;
        if (iso_exact::detail::project(f, p, h)) worst = std::max(worst, iso_geometry::distance(p, v));
    }
    return worst;
}

} // namespace

int main() {
    // Expected number of boundary loops = number of components of the link
    // (x^2 = -y(y^2 - 3z^2) is positive on three sectors that meet only at P).
    const std::vector<Case> cases = {
        {"A1", "x^2+y^2-z^2+x*y*z", 2},       {"A2a", "x^2+y^2+z^3", 1},
        {"A2b", "x^2-y^2+z^3", 1},             {"A2c", "x^3+y^3+z^2+x^2", 1},
        {"D4a", "x^2+y^3+z^3", 1},             {"D4b", "x^2+y^3-3*y*z^2", 3},
        {"E6a", "x^3+y^3+z^3", 1},             {"E6b", "x^3+y^3-z^3", 1},
        {"E6c", "x^3+y^3+z^3+x*y*z", 1},       {"E6d", "x^3+2*y^3+3*z^3+2*x*y*z", 1},
    };
    const int res = 64;
    const double half = 0.4, h = 2.0 * half / res, box = 0.18;
    for (const Case& c : cases) {
        const NodePtr ast = parseExpression(std::string("(") + c.f + ")+k*(x)");
        iso_reconstruct::Settings s;
        s.eps = 1e-3;
        s.resolution = res;
        s.halfExtent = half;
        const auto direct = iso_reconstruct::reconstruct(ast.get(), iso_reconstruct::Method::FieldAverage, s);
        const auto cone = iso_reconstruct::reconstruct(ast.get(), iso_reconstruct::Method::SingularCone, s);
        const iso_exact::Field f(ast.get(), 0.0);
        const auto tDirect = iso_reconstruct::topology(direct.mesh, {0, 0, 0}, box, h);
        const auto tCone = iso_reconstruct::topology(cone.mesh, {0, 0, 0}, box, h);
        const double eDirect = vertexError(direct.mesh, f, 4.0 * h, h), eCone = vertexError(cone.mesh, f, 4.0 * h, h);
        std::printf("%-4s loops %d rings %d | direct: comp %zu chi %lld holes %zu | cone: comp %zu chi %lld holes %zu "
                    "nonmanifold %zu | vertex error near P: direct %.2e cone %.2e | %s\n",
                    c.name, cone.patch.loops, cone.patch.rings, tDirect.components, tDirect.euler, tDirect.holeEdges,
                    tCone.components, tCone.euler, tCone.holeEdges, tCone.nonManifoldEdges, eDirect, eCone,
                    cone.note.c_str());
        check(cone.local.weightsFound && cone.local.nondegenerate, std::string(c.name) + ": principal part");
        check(cone.patch.applied, std::string(c.name) + ": patch applied");
        check(cone.patch.loops == c.loops, std::string(c.name) + ": number of loops");
        check(cone.patch.openChains == 0 && cone.patch.projectionFailures == 0, std::string(c.name) + ": clean patch");
        check(tCone.contractible(), std::string(c.name) + ": cone patch is contractible in Q");
        check(eCone <= 1e-9, std::string(c.name) + ": patch vertices lie on F = 0");
    }

    // The field average is the grid of F: the same mesh as Marching Cubes of F.
    {
        const NodePtr ast = parseExpression("(x^3+y^3+z^3)+k*(x+y)");
        iso_reconstruct::Settings s;
        s.resolution = 32;
        s.halfExtent = 0.4;
        const auto avg = iso_reconstruct::reconstruct(ast.get(), iso_reconstruct::Method::FieldAverage, s);
        const auto direct = iso_reconstruct::fromResult(Isosurface::build(ast.get(), 0.0, 32, 0.4), 0.025);
        check(avg.mesh.t.size() == direct.t.size(), "field average equals the mesh of F");
    }

    // A point that is not singular: the cone patch is refused, the mesh of F is kept.
    {
        const NodePtr ast = parseExpression("(x^2+y^2+z^2-0.04)+k*(x)");
        iso_reconstruct::Settings s;
        s.resolution = 32;
        s.halfExtent = 0.4;
        const auto out = iso_reconstruct::reconstruct(ast.get(), iso_reconstruct::Method::SingularCone, s);
        check(!out.patch.applied && !out.note.empty() && !out.mesh.t.empty(), "non-singular P is refused");
    }

    // Every method runs through the viewer entry point.
    {
        const NodePtr ast = parseExpression("(x^2+y^2-z^2)+k*(1)");
        iso_reconstruct::Settings s;
        s.resolution = 32;
        s.halfExtent = 0.4;
        for (int m = 0; m < static_cast<int>(iso_reconstruct::Method::Count); ++m) {
            const auto out = iso_reconstruct::reconstruct(ast.get(), static_cast<iso_reconstruct::Method>(m), s);
            check(!out.mesh.t.empty(), std::string("method ") + iso_reconstruct::methodId(static_cast<iso_reconstruct::Method>(m)));
        }
    }

    if (failures) {
        std::printf("%d check(s) failed\n", failures);
        return 1;
    }
    std::printf("all reconstruction tests passed\n");
    return 0;
}
