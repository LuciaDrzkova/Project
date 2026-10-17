cmake_minimum_required(VERSION 3.18)
project(IsosurfaceSingularityLab LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

include(FetchContent)
include(CTest)
enable_testing()

# Use the same PMP viewer/rendering setup as LuciaDrzkova/Isosurfaces (report).
set(PMP_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_TESTS OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_DOCS OFF CACHE BOOL "" FORCE)
set(PMP_BUILD_VIS ON CACHE BOOL "" FORCE)
set(PMP_INSTALL ON CACHE BOOL "" FORCE)
set(PMP_STRICT_COMPILATION OFF CACHE BOOL "" FORCE)
FetchContent_Declare(
  pmp-library
  GIT_REPOSITORY https://github.com/pmp-library/pmp-library.git
  GIT_TAG 3.0.0
)
FetchContent_MakeAvailable(pmp-library)

# Some pmp-library versions need the fetched GLFW include directory when
# compiling Dear ImGui's GLFW backend. Guard the optional target/path.
if(TARGET imgui)
  file(GLOB _pmp_glfw_dirs "${pmp-library_SOURCE_DIR}/external/glfw-*")
  list(LENGTH _pmp_glfw_dirs _pmp_glfw_count)
  if(_pmp_glfw_count GREATER 0)
    list(GET _pmp_glfw_dirs 0 _pmp_glfw_dir)
    target_include_directories(imgui PRIVATE "${_pmp_glfw_dir}/include")
  endif()
endif()

find_package(Threads REQUIRED)

# Keep the same executable name used by the original Project repository.
add_executable(isosurface_lab src/main.cpp)
target_include_directories(isosurface_lab PRIVATE src)
target_link_libraries(isosurface_lab PRIVATE pmp pmp_vis Threads::Threads)

# Math/parser/extraction tests are independent of the GUI and PMP viewer.
add_executable(isosurface_lab_tests tests/test_isosurface.cpp)
target_include_directories(isosurface_lab_tests PRIVATE src)
add_test(NAME isosurface_correctness COMMAND isosurface_lab_tests)

# Expanded, headless singularity/grid-phase diagnostics.
add_executable(isosurface_singularity_diagnostics tests/test_singularity_diagnostics.cpp)
target_include_directories(isosurface_singularity_diagnostics PRIVATE src)
add_test(NAME singularity_diagnostics_smoke
    COMMAND isosurface_singularity_diagnostics --quick --resolution 20 --out ${CMAKE_BINARY_DIR}/singularity_diagnostics_smoke_artifacts)

# Sweep every built-in surface through the selected k values and report
# nearby gradients, topology diagnostics, and winning +/- k candidates.
add_executable(isosurface_singularity_sweep_tests tests/test_singularity_sweep.cpp)
target_include_directories(isosurface_singularity_sweep_tests PRIVATE src)
add_test(NAME singularity_k_sweep COMMAND isosurface_singularity_sweep_tests)
set_tests_properties(singularity_k_sweep PROPERTIES WORKING_DIRECTORY "${CMAKE_SOURCE_DIR}")

