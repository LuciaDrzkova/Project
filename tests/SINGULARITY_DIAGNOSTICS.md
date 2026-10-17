# Surface catalogue and test interpretation

## Catalogue (38 surfaces)

Every family is represented by its explicit expression in `test_singularity_diagnostics.cpp` so the test is reproducible without an implicit external preset table.

| ID | Expression | Expected analytic interpretation |
|---|---|---|
| `plane_x3_ky` | `x^3-k*y` | Regular zero level for nonzero k; the origin's gradient is `(0,-k,0)`. |
| `isolated_cusp_x` | `x^2+y^2+z^2-k*x` | Isolated point at k=0; regular small sphere for k!=0. |
| `isolated_cusp_y` | `x^2+y^2+z^2-k*y` | Axis-permuted isolated-point unfolding. |
| `isolated_cusp_z` | `x^2+y^2+z^2-k*z` | Axis-permuted isolated-point unfolding. |
| `quartic_cusp_x` | `x^4+y^2+z^2-k*x` | Higher-order isolated point; regular zero level for k!=0. |
| `quartic_cusp_y` | `x^2+y^4+z^2-k*y` | Axis-permuted quartic cusp. |
| `quartic_cusp_z` | `x^2+y^2+z^4-k*z` | Axis-permuted quartic cusp. |
| `higher_order_cusp` | `x^8+y^4+z^2-k*x` | High-order anisotropic isolated singularity. |
| `mixed_cusp_x` | `x^4+y^2+z^4-k*x` | Different even orders per coordinate. |
| `mixed_cusp_y` | `x^4+y^4+z^2-k*y` | Different even orders per coordinate. |
| `mixed_cusp_z` | `x^2+y^4+z^4-k*z` | Different even orders per coordinate. |
| `mixed_order_x` | `x^6+y^2+z^4-k*x` | Higher mixed-order unfolding. |
| `mixed_order_y` | `x^2+y^6+z^4-k*y` | Higher mixed-order unfolding. |
| `mixed_order_z` | `x^4+y^4+z^2-k*z` | Higher mixed-order unfolding. |
| `mixed_cubic` | `x^3+y^3+z^2+k*z` | For k!=0 the critical point `(0,0,-k/2)` has value `-k^2/4`, not zero. |
| `level_shift_cone` | `x^2+y^2-z^2-k` | Cone singularity at k=0; critical value `-k` for nonzero k. |
| `whitney_umbrella` | `x^2*z-y^2-k*z` | At k>0, two singular zero-level points `(±sqrt(k),0,0)` remain. |
| `monkey_saddle` | `x^3-3*x*y^2+z^2+k*z` | Origin singular at k=0; for nonzero k the critical value at `(0,0,-k/2)` is `-k^2/4`. |
| `cubic_family` | `x^3+y^3+z^3-3*k*x*y*z` | Origin singular for every k; it cannot satisfy a full singularity-removal requirement. |
| `crossing_planes_xy` | `x*y-k*z` | For k!=0, `dF/dz=-k` guarantees regularity. |
| `crossing_planes_yz` | `y*z-k*x` | Axis-permuted regular unfolding. |
| `crossing_planes_xz` | `x*z-k*y` | Axis-permuted regular unfolding. |
| `double_cone_x_shift` | `x^2+y^2-z^2-k*x` | The critical point's value is `-k^2/4`. |
| `double_cone_z_shift` | `x^2+y^2-z^2-k*z` | Critical value differs from zero for k!=0. |
| `triple_plane_unfolding` | `x*y*z-k*(x+y+z)` | Product singularity with linear unfolding; inspected numerically and by its critical equations. |
| `quadratic_crossing_unfolding` | `x^2-y^2-k*z` | For k!=0, `dF/dz=-k` guarantees regularity. |
| `parabolic_cusp` | `x^2+y^2-z^3-k*z` | Cusp-like noncompact surface; real critical points for k<0 have nonzero values. |
| `elliptic_cusp_x` | `x^3+y^2+z^2-k*x` | Cubic cusp unfolding. |
| `elliptic_cusp_y` | `x^2+y^3+z^2-k*y` | Axis-permuted cubic cusp. |
| `elliptic_cusp_z` | `x^2+y^2+z^3-k*z` | Axis-permuted cubic cusp. |
| `radial_quartic_x` | `(x^2+y^2+z^2)^2-k*x` | Radial quartic isolated-point unfolding. |
| `radial_quartic_z` | `(x^2+y^2+z^2)^2-k*z` | Axis-permuted radial quartic. |
| `degree_six_isolated` | `x^6+y^6+z^6-k*x` | High even-degree isolated singularity. |
| `degree_six_mixed` | `x^8+y^6+z^2-k*z` | High-order anisotropic case. |
| `sphere_control` | `x^2+y^2+z^2-0.25` | Smooth closed reference surface. |
| `torus_control` | `(x^2+y^2+z^2+0.2125)^2-1.21(x^2+y^2)` | Smooth polynomial torus, major radius 0.55 and tube radius 0.30. |
| `cylinder_control` | `x^2+y^2-0.25` | Smooth open surface after clipping to the box. |
| `hyperboloid_control` | `x^2+y^2-z^2-0.25` | Smooth clipped saddle-like reference surface. |

## Analytical critical-point witness assessment

`singularity_assessment.csv` records, for each evaluated surface/k pair, the catalog-specific closed-form critical-point candidates, their field values, gradient magnitudes, and whether any witness satisfies both F=0 and grad F=0 within scale-aware tolerances. The rules enumerate the real critical points for each listed family, with a representative point where the critical set is a line, axis, or circle. This is deliberately **not** a generic symbolic solver for arbitrary user expressions. A new surface added to the catalogue must also get a critical-point rule; otherwise the row is marked `UNASSESSED_NO_CRITICAL_POINT_RULE` and the analytic smoke check fails.

`EXPECTED_SINGULARITY_WITNESS_FOUND` means the expected singularity has an explicit critical-point witness. `NO_SINGULAR_WITNESS_ON_ZERO_LEVEL` means the catalogue's critical-point rule found no critical witness on F=0; inspect the rule text for the reason. `UNEXPECTED_SINGULARITY_WITNESS_FOUND` is a regression signal. The CSV also reports the origin separately so persistent-origin families are easy to identify.

## Why add grid-phase tests?

The original extractor samples a symmetric grid from `[-h,h]`. If exact zeros lie on grid vertices or grid-aligned curves, edge interpolation can collapse otherwise distinct edge crossings to the same location. To probe that, phase tests use the same expression with `x`, `y`, and `z` translated by fractions of a cell width. This shifts the surface relative to the same grid. The test records metrics; it does not assume that shifting must improve the result.

The Mixed cubic family is important because `F(t,-t,0)=0` for all t, including the grid-aligned line through the origin. For nonzero k, the gradient along that line is `(3t^2,3t^2,k)`, so it is regular. If anomalies change strongly as the surface is shifted by a fraction of a cell, exact-zero alignment is implicated.

## What each metric means

- `boundary_edges`: all edges incident to one triangle; many are expected on surfaces clipped by the domain box.
- `interior_boundary_edges`: edges incident to one triangle that do not appear to lie on the clipping box; these are more suspicious.
- `nonmanifold_edges`: edges incident to more than two triangles.
- `orientation_conflicts`: an edge incident to two triangles with the same directed orientation in both triangles.
- `near_degenerate_triangles`: output triangles whose area is very small relative to the square of a grid cell.
- `skipped_degenerate_triangles`: triangle candidates dropped by the existing production guard before entering the output mesh; `skipped_degenerate_triangles.csv` records all three candidate positions and area for focused/phase diagnostic runs.
- `phase_variation.csv`: range of each tracked metric over the phase variants; `topology_phase_sensitive=1` means interior-boundary, nonmanifold, orientation-conflict, or skipped-degenerate counts changed with phase.
- `min_gradient`: minimum analytic gradient magnitude at vertices of the emitted mesh. It is a mesh sample, not a certified minimum over the continuous surface.

No single mesh metric proves a mathematical singularity. Interpret these alongside the expression's analytic properties and the geometric resolution.
