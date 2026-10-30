# Install the k -> 0 convergence button and range/count inputs

This installer is designed for the public `LuciaDrzkova/Project` C++ project. It does not connect to GitHub and does not modify the public repository. Run it against your own local checkout.

## Easiest way

1. Download `install_convergence.py`.
2. Put it in the root of your local `Project` folder (the folder containing `CMakeLists.txt`).
3. Run from that folder:

   ```bash
   python install_convergence.py
   ```

   Or pass the project path explicitly:

   ```bash
   python install_convergence.py /path/to/Project
   ```

4. Reconfigure and build:

   ```bash
   cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build -j
   ctest --test-dir build --output-on-failure -R singularity_convergence_regression
   ```

The script is self-contained: it installs `src/SingularityConvergence.h`, `tests/test_singularity_convergence.cpp`, inserts the new controls and button into `src/main.cpp`, and appends the test target to the existing `CMakeLists.txt`. It is intended to be re-runnable without duplicating UI code. Before changing existing files it creates one-time backups ending in `.before_convergence`.

## New UI

Inside **Gradient diagnostics**, the button **Analyze surface convergence** sits beside the k-sweep title. The controls below let you set:

- `k minimum`
- `k maximum`
- `Number of k samples` (3–101)
- mesh resolution

The program samples the selected interval uniformly, including both endpoints. If `k=0` lands on a sample, it is recorded in the CSV but excluded from positive/negative trend fits because `k=0` is the reference surface. If the selected interval does not cross zero, the report can only describe the sampled one-sided range.

The results are exported to `singularity_convergence.csv` in the program's working directory. The output gives point-to-mesh distance, sampled gradient magnitude, field residuals, and a local sampled-vertex Hausdorff proxy. These are numerical diagnostics, not exact Hausdorff distances or proofs. Repeat with higher grid resolution and inspect the CSV before drawing thesis-level conclusions. Large sample counts combined with high mesh resolution increase runtime.

## Notes

- The expression should represent a family affine in `k`, i.e. `F(x,y,z) + k*G(x,y,z)`, when interpreting results as a linear unfolding.
- The script checks that required source anchors match the current main file before editing. If it refuses to install due to a missing anchor, it leaves the project unchanged; share your local `src/main.cpp` and the installer can be adapted.
- No files in the public GitHub repository were changed by preparing this package.
