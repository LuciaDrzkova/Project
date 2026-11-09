# Headless tests and diagnostics tools. They depend only on isolab::core (no
# GUI), so they also build with -DISOLAB_BUILD_APP=OFF.

function(isolab_add_program target source)
  add_executable(${target} ${source})
  target_link_libraries(${target} PRIVATE isolab::core isolab::warnings)
endfunction()

# Parser, AD gradients and marching-cubes correctness/topology.
isolab_add_program(isosurface_lab_tests test_isosurface.cpp)
add_test(NAME isosurface_correctness COMMAND isosurface_lab_tests)

# Grid-aligned 2D cross-section extraction and CSV export.
isolab_add_program(isosurface_slice_extraction_tests test_isosurface_slices.cpp)
add_test(NAME isosurface_2d_slice_extraction
         COMMAND isosurface_slice_extraction_tests
                 --out ${CMAKE_BINARY_DIR}/isosurface_2d_slice_smoke_artifacts)

# Expanded diagnostics over the shared surface catalogue (also a CLI tool;
# see tests/run_resolution_series.py and tests/SINGULARITY_DIAGNOSTICS.md).
isolab_add_program(isosurface_singularity_diagnostics test_singularity_diagnostics.cpp)
add_test(NAME singularity_diagnostics_smoke
         COMMAND isosurface_singularity_diagnostics --quick --resolution 20
                 --out ${CMAKE_BINARY_DIR}/singularity_diagnostics_smoke_artifacts)

# Sweep the built-in surfaces through selected k values.
isolab_add_program(isosurface_singularity_sweep_tests test_singularity_sweep.cpp)
add_test(NAME singularity_k_sweep COMMAND isosurface_singularity_sweep_tests)
set_tests_properties(singularity_k_sweep PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

# Additive mesh-to-mesh interpolation experiment.
isolab_add_program(isosurface_mesh_interpolation_tests test_mesh_to_mesh_interpolation.cpp)
add_test(NAME mesh_to_mesh_interpolation_smoke
         COMMAND isosurface_mesh_interpolation_tests --resolution 20 --k-values 0.01,0.001
                 --surfaces plane_x3_ky,mixed_cubic,crossing_planes_xy
                 --out ${CMAKE_BINARY_DIR}/mesh_to_mesh_interpolation_smoke_artifacts)

# k -> 0 geometric/gradient convergence diagnostics.
isolab_add_program(isosurface_singularity_convergence_tests test_singularity_convergence.cpp)
add_test(NAME singularity_convergence_regression COMMAND isosurface_singularity_convergence_tests)

# Taylor jets, local isolation, stationary-point search, triangle-sampled geometry.
isolab_add_program(isosurface_singularity_research_tests test_singularity_research.cpp)
add_test(NAME singularity_research_analysis COMMAND isosurface_singularity_research_tests)

# Parser and naming helpers of the 25-case F + kG batch workflow.
isolab_add_program(isosurface_singularity_batch_tests test_singularity_batch.cpp)
add_test(NAME singularity_batch_helpers COMMAND isosurface_singularity_batch_tests)

# Analytic singularity criteria for deformation and unfolding families.
isolab_add_program(isosurface_deformation_criteria_tests test_deformation_criteria.cpp)
add_test(NAME deformation_criteria COMMAND isosurface_deformation_criteria_tests)

isolab_add_program(isosurface_unfolding_criterion_tests test_unfolding_criterion.cpp)
add_test(NAME unfolding_criterion COMMAND isosurface_unfolding_criterion_tests)

# Slice snapshot (rectangular grid window over selected slices -> SVG + CSV).
isolab_add_program(isosurface_slice_snapshot_tests test_slice_snapshot.cpp)
add_test(NAME slice_snapshot
         COMMAND isosurface_slice_snapshot_tests
                 --out ${CMAKE_BINARY_DIR}/slice_snapshot_test_artifacts)

# k -> 0 convergence classification (rates, singular-point behaviour, prediction).
isolab_add_program(isosurface_convergence_classifier_tests test_convergence_classifier.cpp)
add_test(NAME convergence_classifier COMMAND isosurface_convergence_classifier_tests)

# Headless classification tool (one F, a list of G): tables for the thesis.
isolab_add_program(isosurface_convergence_classification classify_families.cpp)

# Local algebra, certified stationary points, exact surface distances, direction sweep.
isolab_add_program(isosurface_exact_analysis_tests test_exact_analysis.cpp)
add_test(NAME exact_analysis COMMAND isosurface_exact_analysis_tests)

# Direction sweep of linear perturbations for one F.
isolab_add_program(isosurface_direction_sweep direction_sweep.cpp)

# Reproduces the complete experiment set of the thesis (chapter 5).
isolab_add_program(isosurface_thesis_experiments thesis_experiments.cpp)
