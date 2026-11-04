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
- **Convergence analysis**: for a range of `k`, measures the point-to-mesh
  distance, local gradients, residuals and Hausdorff-type distances to the
  `k = 0` reference. It also runs Taylor-jet, isolation and stationary-point
  diagnostics. Results go to `experiments/run_<timestamp>/` (four CSV files).
- **Batch**: runs the same analysis for `F + k·G_i` over 25 perturbations
  `G_i`.

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
│   ├── SurfaceCatalogue.h the 34 parameterised families (single source of truth)
│   ├── Grid.h             lattice coordinates shared by every component
│   ├── Geometry.h         points, point–triangle distance, KD-tree, BVH
│   ├── Parallel.h         fork-join parallelFor (nested calls run serially)
│   └── CsvFormat.h        CSV quoting and fast round-trip number formatting
├── analysis/    thesis diagnostics built on core
│   ├── SingularityConvergence.h  S_k → S_0 metrics, affine F + kG split, metadata
│   ├── SingularityResearch.h     Taylor jets, isolation, Newton search, mesh
│   │                             comparison; analyzeFamily() runs both together
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

### Analysis pipeline (`analysis/`)

`iso_research::analyzeFamily()` extracts each `k` once, in parallel, and runs
the convergence metrics and the research triangle-geometry comparison on that
same mesh. The stationary-point search also runs in parallel per `k`. All
results are written in the requested `k` order, independent of scheduling.

## Tests and diagnostics tools

| Executable | Purpose |
|---|---|
| `isosurface_lab_tests` | parser, AD gradients, marching-cubes topology (spheres, torus, cusp scaling) |
| `isosurface_slice_extraction_tests` | slice lattice and CSV columns |
| `isosurface_singularity_convergence_tests` | convergence trends and CSV metadata |
| `isosurface_singularity_research_tests` | jets, isolation, Newton search, triangle distances |
| `isosurface_singularity_batch_tests` | batch parsing and folder naming |
| `isosurface_deformation_criteria_tests`, `isosurface_unfolding_criterion_tests` | analytic singularity criteria |
| `isosurface_singularity_sweep_tests` | k sweep over the built-in families (writes `singularity_sweep_artifacts/`) |
| `isosurface_mesh_interpolation_tests` | mesh-to-mesh interpolation experiment |
| `isosurface_singularity_diagnostics` | full catalogue diagnostics CLI, see `tests/SINGULARITY_DIAGNOSTICS.md` and `tests/run_resolution_series.py` |

All measurements are numerical diagnostics, not proofs. Read the CSV
metadata columns before drawing conclusions from them.
