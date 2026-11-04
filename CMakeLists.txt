cmake_minimum_required(VERSION 3.18)
project(IsosurfaceSingularityLab VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "Build type" FORCE)
endif()

# Every executable (viewer, tests, diagnostics tools) lands directly in the
# build directory; tests/run_resolution_series.py and the docs rely on that.
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}")

option(ISOLAB_BUILD_APP "Build the interactive viewer (downloads pmp-library)" ON)
include(CTest) # provides BUILD_TESTING (ON by default)

# -----------------------------------------------------------------------------
# isolab::core - header-only math + analysis library (src/core, src/analysis).
# No GUI dependencies; shared by the viewer and every test/diagnostics tool.
# -----------------------------------------------------------------------------
find_package(Threads REQUIRED)

add_library(isolab_core INTERFACE)
add_library(isolab::core ALIAS isolab_core)
target_include_directories(isolab_core INTERFACE "${CMAKE_CURRENT_SOURCE_DIR}/src")
target_compile_features(isolab_core INTERFACE cxx_std_17)
target_link_libraries(isolab_core INTERFACE Threads::Threads)

# Floating-point std::to_chars makes the large CSV exports ~4x faster. libc++
# provides it without advertising __cpp_lib_to_chars, so probe for it.
include(CheckCXXSourceCompiles)
check_cxx_source_compiles([[
  #include <charconv>
  int main() {
    char buffer[40];
    return std::to_chars(buffer, buffer + 40, 0.1, std::chars_format::general, 17).ptr == buffer;
  }
]] ISOLAB_HAVE_FLOAT_TO_CHARS)
if(ISOLAB_HAVE_FLOAT_TO_CHARS)
  target_compile_definitions(isolab_core INTERFACE ISO_HAVE_FLOAT_TO_CHARS=1)
endif()

# Strict warnings for this project's own targets only (not for dependencies).
add_library(isolab_warnings INTERFACE)
add_library(isolab::warnings ALIAS isolab_warnings)
target_compile_options(isolab_warnings INTERFACE
  $<$<CXX_COMPILER_ID:GNU,Clang,AppleClang>:-Wall -Wextra -Wpedantic>
  $<$<CXX_COMPILER_ID:MSVC>:/W4>)

# -----------------------------------------------------------------------------
# isosurface_lab - the interactive PMP/ImGui viewer (src/app).
# -----------------------------------------------------------------------------
if(ISOLAB_BUILD_APP)
  include(FetchContent)
  # Same PMP viewer/rendering setup as LuciaDrzkova/Isosurfaces (report).
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

  # Treat dependency headers as system headers so their warnings (e.g. in
  # imgui.h) do not drown out warnings in our own code.
  foreach(_dependency IN ITEMS pmp pmp_vis imgui)
    if(TARGET ${_dependency})
      get_target_property(_includes ${_dependency} INTERFACE_INCLUDE_DIRECTORIES)
      if(_includes)
        set_target_properties(${_dependency} PROPERTIES INTERFACE_SYSTEM_INCLUDE_DIRECTORIES "${_includes}")
      endif()
    endif()
  endforeach()

  add_executable(isosurface_lab
    src/app/main.cpp
    src/app/IsosurfaceViewer.cpp
    src/app/ViewerPanels.cpp
    src/app/SliceGallery.cpp
    src/app/AnalysisJobs.cpp
    src/app/OverlayRenderer.cpp
  )
  target_include_directories(isosurface_lab PRIVATE src/app)

  # Embed the source directory through a generated header (experiments/ lives
  # next to the sources). A header rather than a compile definition keeps
  # paths with spaces intact.
  set(ISOSURFACE_PROJECT_ROOT_VALUE "${CMAKE_SOURCE_DIR}")
  string(REPLACE "\\" "\\\\" ISOSURFACE_PROJECT_ROOT_VALUE "${ISOSURFACE_PROJECT_ROOT_VALUE}")
  string(REPLACE "\"" "\\\"" ISOSURFACE_PROJECT_ROOT_VALUE "${ISOSURFACE_PROJECT_ROOT_VALUE}")
  set(ISOSURFACE_GENERATED_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/generated")
  configure_file(src/app/ProjectPaths.h.in "${ISOSURFACE_GENERATED_INCLUDE_DIR}/ProjectPaths.h" @ONLY)
  target_include_directories(isosurface_lab PRIVATE "${ISOSURFACE_GENERATED_INCLUDE_DIR}")

  target_link_libraries(isosurface_lab PRIVATE isolab::core isolab::warnings pmp pmp_vis)
endif()

# -----------------------------------------------------------------------------
# Tests and headless diagnostics tools.
# -----------------------------------------------------------------------------
if(BUILD_TESTING)
  add_subdirectory(tests)
endif()
