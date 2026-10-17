# Shared parameterized surface catalogue and diagnostics

The viewer and `test_singularity_diagnostics.cpp` use `src/SurfaceCatalogue.h` as their single source of truth.
The catalogue contains 34 parameterized families. Every expression depends on `k`; k-independent control shapes are deliberately excluded.

| ID | Expression | Initial k | Expected behavior |
|---|---|---:|---|
| `plane_x3_ky` | `x*x*x-k*y` | `0` | At k=0 the whole plane x=0 is singular; for k!=0, dF/dy=-k. |
| `isolated_cusp_x` | `x*x*x+y*y+z*z+k*x` | `0` | The cubic cusp has an isolated singular point at the origin when k=0. |
| `isolated_cusp_y` | `x*x+y*y*y+z*z+k*y` | `0` | Axis-permuted cubic cusp; the zero-level surface is regular for k!=0. |
| `isolated_cusp_z` | `x*x+y*y+z*z*z+k*z` | `0` | Axis-permuted cubic cusp; the zero-level surface is regular for k!=0. |
| `quartic_cusp_x` | `x*x*x*x+y*y+z*z+k*x` | `0` | Quartic term in X; the k=0 fiber has an isolated singularity. |
| `quartic_cusp_y` | `x*x+y*y*y*y+z*z+k*y` | `0` | Axis-permuted quartic cusp. |
| `quartic_cusp_z` | `x*x+y*y+z*z*z*z+k*z` | `0` | Axis-permuted quartic cusp. |
| `higher_order_cusp` | `x*x*x*x*x+y*y+z*z+k*x` | `0` | Fifth-order X term; the critical points move off the zero level for k!=0. |
| `mixed_cusp_y` | `x*x*x+y*y*y*y+z*z+k*y` | `0` | Mixed cubic/quartic orders with a Y-directed unfolding. |
| `mixed_cusp_z` | `x*x*x+y*y+z*z*z*z+k*z` | `0` | Mixed cubic/quartic orders with a Z-directed unfolding. |
| `mixed_cusp_x` | `x*x+y*y*y+z*z*z*z+k*x` | `0` | Mixed orders with a linear X perturbation. |
| `mixed_order_z` | `x*x*x*x+y*y*y+z*z+k*z` | `0` | The critical point is (0,0,-k/2), with F=-k^2/4 for k!=0. |
| `mixed_cubic` | `x*x*x+y*y*y+z*z+k*z` | `0` | The critical point is (0,0,-k/2), with F=-k^2/4 for k!=0. |
| `level_shift_cone` | `x*x+y*y-z*z-k` | `0` | A cone at k=0; nonzero k shifts the critical value away from zero. |
| `whitney_umbrella` | `x*x*z-y*y-k*z` | `0` | For k>0, two additional singular points (±sqrt(k),0,0) lie on F=0. |
| `monkey_saddle` | `x*x*x-3*x*y*y-k` | `0` | The singular z-axis at k=0 is removed from the zero level by nonzero k. |
| `cubic_family` | `x*x*x+y*y*y+z*z*z-3*k*x*y*z` | `0` | The origin remains singular on F=0 for every k; this is an intentional counterexample. |
| `higher_anisotropic_cusp` | `pow(x,8)+pow(y,4)+z*z-k*x` | `0.1` | Anisotropic eighth-/fourth-/second-order unfolding. |
| `mixed_even_cusp_x` | `x*x*x*x+y*y+z*z*z*z-k*x` | `0.1` | Anisotropic even-power singularity unfolded along X. |
| `mixed_even_cusp_y` | `x*x*x*x+y*y*y*y+z*z-k*y` | `0.1` | Anisotropic even-power singularity unfolded along Y. |
| `mixed_even_cusp_z` | `x*x+y*y*y*y+z*z*z*z-k*z` | `0.1` | Anisotropic even-power singularity unfolded along Z. |
| `mixed_order_x` | `pow(x,6)+y*y+z*z*z*z-k*x` | `0.1` | Sixth-/second-/fourth-order terms, with an X-directed unfolding. |
| `mixed_order_y` | `x*x+pow(y,6)+z*z*z*z-k*y` | `0.1` | Second-/sixth-/fourth-order terms, with a Y-directed unfolding. |
| `mixed_even_order_z` | `x*x*x*x+y*y*y*y+z*z-k*z` | `0.1` | A mixed-order cousin of the existing mixed-order Z preset. |
| `monkey_saddle_z_unfolding` | `x*x*x-3*x*y*y+z*z+k*z` | `0` | A genuinely different unfolding from the original level-shifted monkey saddle. |
| `crossing_planes_xy` | `x*y-k*z` | `0` | At nonzero k, dF/dz=-k everywhere, so the zero level is a graph. |
| `double_cone_x_shift` | `x*x+y*y-z*z-k*x` | `0.1` | Coordinate-directed unfolding of the double cone; initially shown at k=0.1. |
| `double_cone_z_shift` | `x*x+y*y-z*z-k*z` | `0.1` | Z-directed unfolding of the double cone; initially shown at k=0.1. |
| `triple_plane_unfolding` | `x*y*z-k*(x+y+z)` | `0` | Product singularity unfolding of the three coordinate planes. |
| `quadratic_crossing_unfolding` | `x*x-y*y-k*z` | `0.1` | For nonzero k, dF/dz=-k everywhere; the crossing opens into a graph. |
| `parabolic_cusp` | `x*x+y*y-z*z*z-k*z` | `0` | The k=0 fiber is a parabolic cusp; nonzero k changes its local geometry. |
| `radial_quartic_x` | `pow(x*x+y*y+z*z,2)-k*x` | `0.1` | A radially symmetric fourth-degree field unfolded along X. |
| `degree_six_isolated` | `pow(x,6)+pow(y,6)+pow(z,6)-k*x` | `0.1` | A higher-degree positive-definite form; nonzero default k keeps the surface visible. |
| `degree_six_mixed` | `pow(x,8)+pow(y,6)+z*z-k*z` | `0.1` | An anisotropic high-order form unfolded along Z. |

## Running

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --target isosurface_lab isosurface_singularity_diagnostics
ctest --test-dir build --output-on-failure -R 'singularity_diagnostics_smoke'
./build/isosurface_singularity_diagnostics --resolution 48 --out singularity_diagnostics_artifacts
```

A small set of positive-even leading forms uses an initial `k=0.1` in the viewer, because the zero-level set at `k=0` would be only the origin and therefore look like an empty/point-like duplicate. Their tests still include `k=0` as the singular reference.

Analytical checks are catalog-specific witness checks, not a general symbolic proof engine. Mesh topology findings should be read separately from analytical singularity results.
