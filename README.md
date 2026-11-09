# Isosurface Singularity Lab

Interactive viewer and headless toolkit for studying singularities of
parameterised implicit surfaces

    S_k = { (x, y, z) : F(x, y, z, k) = 0 }

and how they unfold as the parameter `k` moves away from 0 (typically
`F_k = F + k·G`). The surface is extracted with Marching Cubes. Gradients are
exact (forward-mode automatic differentiation), and the analysis tools
measure how `S_k` approaches `S_0` near a chosen point.

## Build and run

Requirements: a C++17 compiler, CMake ≥ 3.18 and OpenGL 3.3. The viewer
downloads [pmp-library](https://github.com/pmp-library/pmp-library) 3.0 on
first configure.

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j
./build/isosurface_lab                      # the interactive viewer
ctest --test-dir build --output-on-failure  # all tests
```

Use `-DISOLAB_BUILD_APP=OFF` to build only the headless library, tests and
diagnostics tools (no pmp download, no OpenGL needed).

## What the viewer does

- **Function and k**: type any expression in `x, y, z, k` (operators
  `+ - * / ^`, functions `sin cos tan exp log/ln sqrt abs pow min max`) or
  pick one of 34 catalogued families, then drag `k`.
- **Surface**: coloured by `|∇F|` (purple = flat or singular, red = steep);
  optional gradient glyphs, sampling grid and cell-corner inspector.
- **2D slices**: grid-aligned cross-sections with marching squares, a slice
  gallery and full-precision CSV export.
- **Slice snapshots**: in the gallery, **Create snapshot...** picks three or
  more slices and a rectangle of grid nodes. It saves
  `snapshots/<F>_slices_<list>_u<a>-<b>_v<c>-<d>.svg` (node isovalues,
  contour, corner coordinates) plus a `.csv` with the same nodes.
- **Convergence analysis**: for a range of `k`, measures the point-to-mesh
  distance, local gradients, residuals and Hausdorff-type distances to the
  `k = 0` reference. It also runs Taylor-jet, isolation and stationary-point
  diagnostics. All local quantities are measured in the box
  `P + [−a, a]³`, `a = 0.45 × domain half-extent`. Results go to
  `experiments/run_<timestamp>/` (seven CSV files, see below).
- **Convergence classification**: **Classify convergence (k → 0)**, under the
  k-sweep, uses the same k selection (optionally log-spaced) and decides, for
  k < 0 and k > 0 separately, how `S_k → S_0` near P:
  - the rate `d_H ~ |k|^α` (linear, Hölder, or not converging);
  - what happens to P (it leaves `S_k`, or stays as a regular or singular
    point);
  - whether singular points of `S_k` are removed, persist at P, or move into P
    like `|k|^γ`, and their local type;
  - the exponent predicted from the Taylor jets at P;
  - the Milnor and Tjurina numbers of F at P and the essential coordinates
    of G (when F and G are polynomials of degree ≤ 3).

  Results go to `experiments/classify_<timestamp>/convergence_classification.csv`.
- **Batch**: runs the same analysis for `F + k·G_i` over 25 perturbations
  `G_i`. Each run also gets its classification, and the batch writes
  `classification_summary.csv` (one row per G and side of k = 0).

Long analyses run in the background, and the viewer only redraws when
something changes, so it uses no CPU while idle.

## Architecture

```
src/
├── core/        header-only math engine, no GUI dependencies
│   ├── Expr.h             parser → AST; templated evaluation (double, Dual3, jets)
│   ├── ExprProgram.h      AST compiled to a flat program; batched row evaluation
│   ├── Dual3.h            forward-mode AD number (value + exact gradient)
│   ├── Isosurface.h       marching cubes: sampling, zero tie-break, extraction
│   ├── IsosurfaceSlices.h grid-aligned 2D slices + CSV export
│   ├── SliceSnapshot.h    grid-window snapshot of selected slices (SVG + CSV)
│   ├── FileNames.h        file-system-safe names from expressions
│   ├── SurfaceCatalogue.h the 34 parameterised families (single source of truth)
│   ├── Grid.h             lattice coordinates shared by every component
│   ├── Geometry.h         points, closest points on triangles, KD-tree, BVH
│   ├── Parallel.h         fork-join parallelFor (nested calls run serially)
│   └── CsvFormat.h        CSV quoting and fast round-trip number formatting
├── analysis/    thesis diagnostics built on core
│   ├── SingularityConvergence.h  S_k → S_0 metrics, affine F + kG split, metadata
│   ├── SingularityResearch.h     Taylor jets, isolation, real Newton search, mesh-to-mesh comparison
│   ├── CubicPolynomial.h         exact polynomials of degree ≤ 3, real and complex evaluation
│   ├── LocalAlgebra.h            Milnor/Tjurina algebra, monomial basis, essential coordinates of G
│   ├── StationaryPoints.h        certified stationary points (complex Newton, Milnor-number check)
│   ├── ExactDistance.h           distances on the exact surfaces, extraction accuracy
│   ├── FamilyAnalysis.h          analyzeFamily(): all analyses on one shared mesh pass; CSV export
│   ├── DirectionSweep.h          map of linear perturbation directions, discriminant by bisection
│   ├── ConvergenceClassifier.h   k -> 0 classification: rates with CIs, singular points, prediction
│   └── SingularityBatch.h        25-family F + kG batch runner
└── app/         the PMP/ImGui viewer
    ├── main.cpp                entry point
    ├── IsosurfaceViewer.h/.cpp class declaration; frame loop, surface, picking, overlays
    ├── ViewerPanels.cpp        ImGui control panel and setup dialog
    ├── SliceGallery.cpp        2D slice state and gallery window
    ├── AnalysisJobs.cpp        background convergence analysis and batch
    ├── OverlayRenderer.h/.cpp  all raw OpenGL for helper geometry
    ├── BackgroundJob.h         worker-thread job polled from the UI thread
    └── ColorMap.h              colour ramp
tests/           unit/regression tests and headless diagnostics tools
```

Dependencies only point downwards: `app → analysis → core`. CMake exposes
`core` and `analysis` as the header-only target `isolab::core`, which every
test and tool links. The viewer adds pmp-library on top.

### Extraction pipeline (`core/Isosurface.h`)

1. **Sample**: `F` is compiled once (`ExprProgram`) and evaluated row by row
   at the lattice nodes `-h + 2h·i/N`, in parallel across x-slabs.
2. **Zero tie-break**: a node with `F = 0` exactly, or `|F|` below 1e-10 of
   its largest face neighbour (rounding noise), would make several edge
   crossings collapse onto it. If it has two or more face neighbours of
   opposite sign, it gets a positive bias of 5% of the smallest off-surface
   neighbour value. Grid-aligned zero sets with a single crossing per node,
   such as the plane `x = 0` of `x³` at `k = 0`, stay exact.
3. **March**: each cut edge is interpolated once, always from its lower to
   its upper node, and cached in a rolling two-plane table. Neighbouring
   cells therefore share bit-identical vertices. The result carries
   `vertexIds` so consumers can weld the mesh without hashing positions.
4. **Gradients**: each vertex normal and `|∇F|` comes from evaluating `F`
   with `Dual3` at the vertex.

### Convergence classification (`analysis/ConvergenceClassifier.h`)

**Prediction.** For `F + kG` at a singular point P, the classifier reads off:
- the weighted principal part of F (weights `w`, such that `P(t^w x) = t·P(x)`);
- the weighted order `d_w(G)` of G.

If `d_w(G) < 1`, the core of `S_k` lives on the inner scale
`λ = |k|^{1/(1−d_w(G))}`, which predicts
`α ∈ [w_min, w_max] / (1 − d_w(G))` (capped at 1). For a nondegenerate cubic
cone and G of order m this is `α = 1/(3−m)`.

**Measurement.** When F and G are polynomials of degree ≤ 3, `d_H` is the
local Hausdorff distance between the *exact* surfaces in the box
`P + [−a, a]³` (`analysis/ExactDistance.h`): mesh points are moved onto
`H_k = 0` by Newton's method, and point-to-surface distances are refined by
Newton's method on the Lagrange conditions. Extra seeds from nested finer
grids around P catch features smaller than a grid cell. Otherwise `d_H` is
measured between meshes and fitted only where `h/2 ≤ d_H ≤ 0.6·a`. Every
exponent comes with its standard error and a 95% confidence interval.

**Singular points.** With cubic F and G, the stationary points of `H_k` near P
come from a certified search (`analysis/StationaryPoints.h`): their
multiplicities must add up to the Milnor number of F at P, which guarantees
that none is missing.

### Analysis pipeline (`analysis/`)

`iso_research::analyzeFamily()` (`FamilyAnalysis.h`) extracts each `k` once,
in parallel, and runs the mesh diagnostics, the mesh-to-mesh comparison and
the exact-surface measurements on that same mesh. It writes:

| File | Content |
|---|---|
| `singularity_convergence.csv` | mesh diagnostics per k, with full metadata |
| `singularity_taylor.csv` | Taylor data of F and G at P |
| `singularity_critical_points.csv` | real Newton search (uncertified) |
| `singularity_triangle_geometry.csv` | mesh-to-mesh distances |
| `exact_geometry.csv` | exact local Hausdorff distance, dist(P, S_k), extraction error of the mesh |
| `stationary_points_certified.csv` | stationary points per k, with multiplicity, type and certificate |
| `local_algebra.csv` | Milnor/Tjurina numbers, basis, essential coordinates of G |

## Tests and diagnostics tools

| Executable | Purpose |
|---|---|
| `isosurface_lab_tests` | parser, AD gradients, marching-cubes topology (spheres, torus, cusp scaling) |
| `isosurface_slice_extraction_tests` | slice lattice and CSV columns |
| `isosurface_singularity_convergence_tests` | convergence trends and CSV metadata |
| `isosurface_singularity_research_tests` | jets, isolation, Newton search, triangle distances |
| `isosurface_singularity_batch_tests` | batch parsing and folder naming |
| `isosurface_convergence_classifier_tests` | predicted exponents (A1, A2, D4, Ẽ6) and end-to-end classification of known families |
| `isosurface_convergence_classification` | CLI: classification table for one F and a list of G, e.g. `--F "x^3+y^3+z^3" --G "1;x;x+y;x^2" --kmax 0.01 --out table.csv` |
| `isosurface_exact_analysis_tests` | Milnor/Tjurina numbers of the ten thesis surfaces, certified stationary points, exact distances, direction sweep |
| `isosurface_direction_sweep` | CLI: direction map of linear perturbations, e.g. `--F "x^3+y^3+z^3" --k -0.01 --out sweep` |
| `isosurface_thesis_experiments` | regenerates the complete experiment set of the thesis (about one hour): `isosurface_thesis_experiments out_dir` |
| `isosurface_deformation_criteria_tests`, `isosurface_unfolding_criterion_tests` | analytic singularity criteria |
| `isosurface_singularity_sweep_tests` | k sweep over the built-in families (writes `singularity_sweep_artifacts/`) |
| `isosurface_mesh_interpolation_tests` | mesh-to-mesh interpolation experiment |
| `isosurface_singularity_diagnostics` | full catalogue diagnostics CLI, see `tests/SINGULARITY_DIAGNOSTICS.md` and `tests/run_resolution_series.py` |

All measurements are numerical diagnostics, not proofs. Read the CSV
metadata columns before drawing conclusions from them.
