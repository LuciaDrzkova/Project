#pragma once
// Shared source of truth for the viewer presets and headless singularity tests.
// All catalogue entries depend on k. K-independent geometry controls are intentionally excluded.
#include <string>
#include <vector>

namespace iso_surfaces {
enum class Behavior { RegularForNonzeroK, WhitneyUmbrella, PersistentOrigin };
struct SurfaceCase {
    std::string id;
    std::string label;
    std::string expression;
    Behavior behavior;
    std::string note;
    double defaultK = 0.0;
};
inline const std::vector<SurfaceCase>& catalogue() {
    static const std::vector<SurfaceCase> data = {
        {"plane_x3_ky", "Plane example: x^3-k*y (non-isolated)", "x*x*x-k*y", Behavior::RegularForNonzeroK, "At k=0 the whole plane x=0 is singular; for k!=0, dF/dy=-k.", 0},
        {"isolated_cusp_x", "Isolated cusp X: x^3+y^2+z^2+k*x", "x*x*x+y*y+z*z+k*x", Behavior::RegularForNonzeroK, "The cubic cusp has an isolated singular point at the origin when k=0.", 0},
        {"isolated_cusp_y", "Isolated cusp Y: x^2+y^3+z^2+k*y", "x*x+y*y*y+z*z+k*y", Behavior::RegularForNonzeroK, "Axis-permuted cubic cusp; the zero-level surface is regular for k!=0.", 0},
        {"isolated_cusp_z", "Isolated cusp Z: x^2+y^2+z^3+k*z", "x*x+y*y+z*z*z+k*z", Behavior::RegularForNonzeroK, "Axis-permuted cubic cusp; the zero-level surface is regular for k!=0.", 0},
        {"quartic_cusp_x", "Quartic cusp X: x^4+y^2+z^2+k*x", "x*x*x*x+y*y+z*z+k*x", Behavior::RegularForNonzeroK, "Quartic term in X; the k=0 fiber has an isolated singularity.", 0},
        {"quartic_cusp_y", "Quartic cusp Y: x^2+y^4+z^2+k*y", "x*x+y*y*y*y+z*z+k*y", Behavior::RegularForNonzeroK, "Axis-permuted quartic cusp.", 0},
        {"quartic_cusp_z", "Quartic cusp Z: x^2+y^2+z^4+k*z", "x*x+y*y+z*z*z*z+k*z", Behavior::RegularForNonzeroK, "Axis-permuted quartic cusp.", 0},
        {"higher_order_cusp", "Higher-order cusp: x^5+y^2+z^2+k*x", "x*x*x*x*x+y*y+z*z+k*x", Behavior::RegularForNonzeroK, "Fifth-order X term; the critical points move off the zero level for k!=0.", 0},
        {"mixed_cusp_y", "Mixed cusp Y: x^3+y^4+z^2+k*y", "x*x*x+y*y*y*y+z*z+k*y", Behavior::RegularForNonzeroK, "Mixed cubic/quartic orders with a Y-directed unfolding.", 0},
        {"mixed_cusp_z", "Mixed cusp Z: x^3+y^2+z^4+k*z", "x*x*x+y*y+z*z*z*z+k*z", Behavior::RegularForNonzeroK, "Mixed cubic/quartic orders with a Z-directed unfolding.", 0},
        {"mixed_cusp_x", "Mixed cusp X: x^2+y^3+z^4+k*x", "x*x+y*y*y+z*z*z*z+k*x", Behavior::RegularForNonzeroK, "Mixed orders with a linear X perturbation.", 0},
        {"mixed_order_z", "Mixed-order Z: x^4+y^3+z^2+k*z", "x*x*x*x+y*y*y+z*z+k*z", Behavior::RegularForNonzeroK, "The critical point is (0,0,-k/2), with F=-k^2/4 for k!=0.", 0},
        {"mixed_cubic", "Mixed cubic: x^3+y^3+z^2+k*z", "x*x*x+y*y*y+z*z+k*z", Behavior::RegularForNonzeroK, "The critical point is (0,0,-k/2), with F=-k^2/4 for k!=0.", 0},
        {"level_shift_cone", "Level-shift cone: x^2+y^2-z^2-k", "x*x+y*y-z*z-k", Behavior::RegularForNonzeroK, "A cone at k=0; nonzero k shifts the critical value away from zero.", 0},
        {"whitney_umbrella", "Whitney umbrella: x^2*z-y^2-k*z", "x*x*z-y*y-k*z", Behavior::WhitneyUmbrella, "For k>0, two additional singular points (±sqrt(k),0,0) lie on F=0.", 0},
        {"monkey_saddle", "Monkey saddle: x^3-3*x*y^2-k", "x*x*x-3*x*y*y-k", Behavior::RegularForNonzeroK, "The singular z-axis at k=0 is removed from the zero level by nonzero k.", 0},
        {"cubic_family", "Cubic family: x^3+y^3+z^3-3*k*x*y*z", "x*x*x+y*y*y+z*z*z-3*k*x*y*z", Behavior::PersistentOrigin, "The origin remains singular on F=0 for every k; this is an intentional counterexample.", 0},
        {"higher_anisotropic_cusp", "Higher anisotropic cusp: x^8+y^4+z^2-k*x", "pow(x,8)+pow(y,4)+z*z-k*x", Behavior::RegularForNonzeroK, "Anisotropic eighth-/fourth-/second-order unfolding.", 0.10000000000000001},
        {"mixed_even_cusp_x", "Mixed even cusp X: x^4+y^2+z^4-k*x", "x*x*x*x+y*y+z*z*z*z-k*x", Behavior::RegularForNonzeroK, "Anisotropic even-power singularity unfolded along X.", 0.10000000000000001},
        {"mixed_even_cusp_y", "Mixed even cusp Y: x^4+y^4+z^2-k*y", "x*x*x*x+y*y*y*y+z*z-k*y", Behavior::RegularForNonzeroK, "Anisotropic even-power singularity unfolded along Y.", 0.10000000000000001},
        {"mixed_even_cusp_z", "Mixed even cusp Z: x^2+y^4+z^4-k*z", "x*x+y*y*y*y+z*z*z*z-k*z", Behavior::RegularForNonzeroK, "Anisotropic even-power singularity unfolded along Z.", 0.10000000000000001},
        {"mixed_order_x", "Higher mixed-order cusp X: x^6+y^2+z^4-k*x", "pow(x,6)+y*y+z*z*z*z-k*x", Behavior::RegularForNonzeroK, "Sixth-/second-/fourth-order terms, with an X-directed unfolding.", 0.10000000000000001},
        {"mixed_order_y", "Higher mixed-order cusp Y: x^2+y^6+z^4-k*y", "x*x+pow(y,6)+z*z*z*z-k*y", Behavior::RegularForNonzeroK, "Second-/sixth-/fourth-order terms, with a Y-directed unfolding.", 0.10000000000000001},
        {"mixed_even_order_z", "Higher mixed-order cusp Z: x^4+y^4+z^2-k*z", "x*x*x*x+y*y*y*y+z*z-k*z", Behavior::RegularForNonzeroK, "A mixed-order cousin of the existing mixed-order Z preset.", 0.10000000000000001},
        {"monkey_saddle_z_unfolding", "Monkey saddle with z unfolding: x^3-3*x*y^2+z^2+k*z", "x*x*x-3*x*y*y+z*z+k*z", Behavior::RegularForNonzeroK, "A genuinely different unfolding from the original level-shifted monkey saddle.", 0},
        {"crossing_planes_xy", "Crossing planes XY: x*y-k*z", "x*y-k*z", Behavior::RegularForNonzeroK, "At nonzero k, dF/dz=-k everywhere, so the zero level is a graph.", 0},
        {"double_cone_x_shift", "Double cone X unfolding: x^2+y^2-z^2-k*x", "x*x+y*y-z*z-k*x", Behavior::RegularForNonzeroK, "Coordinate-directed unfolding of the double cone; initially shown at k=0.1.", 0.10000000000000001},
        {"double_cone_z_shift", "Double cone Z unfolding: x^2+y^2-z^2-k*z", "x*x+y*y-z*z-k*z", Behavior::RegularForNonzeroK, "Z-directed unfolding of the double cone; initially shown at k=0.1.", 0.10000000000000001},
        {"triple_plane_unfolding", "Three coordinate planes: x*y*z-k*(x+y+z)", "x*y*z-k*(x+y+z)", Behavior::RegularForNonzeroK, "Product singularity unfolding of the three coordinate planes.", 0},
        {"quadratic_crossing_unfolding", "Crossing quadratic sheets: x^2-y^2-k*z", "x*x-y*y-k*z", Behavior::RegularForNonzeroK, "For nonzero k, dF/dz=-k everywhere; the crossing opens into a graph.", 0.10000000000000001},
        {"parabolic_cusp", "Parabolic cusp: x^2+y^2-z^3-k*z", "x*x+y*y-z*z*z-k*z", Behavior::RegularForNonzeroK, "The k=0 fiber is a parabolic cusp; nonzero k changes its local geometry.", 0},
        {"radial_quartic_x", "Radial quartic unfolding: (x^2+y^2+z^2)^2-k*x", "pow(x*x+y*y+z*z,2)-k*x", Behavior::RegularForNonzeroK, "A radially symmetric fourth-degree field unfolded along X.", 0.10000000000000001},
        {"degree_six_isolated", "Degree-six isolated unfolding: x^6+y^6+z^6-k*x", "pow(x,6)+pow(y,6)+pow(z,6)-k*x", Behavior::RegularForNonzeroK, "A higher-degree positive-definite form; nonzero default k keeps the surface visible.", 0.10000000000000001},
        {"degree_six_mixed", "Degree-six mixed unfolding: x^8+y^6+z^2-k*z", "pow(x,8)+pow(y,6)+z*z-k*z", Behavior::RegularForNonzeroK, "An anisotropic high-order form unfolded along Z.", 0.10000000000000001},
    };
    return data;
}
} // namespace iso_surfaces
