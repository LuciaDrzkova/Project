cmake_minimum_required(VERSION 3.15)
project(IsosurfaceSingularityLab CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

include(FetchContent)

# ---- GLFW (windowing + GL context) ----
# If your system already has a dev package (Ubuntu/Debian:
#   sudo apt install libglfw3-dev libglu1-mesa-dev
# macOS:  brew install glfw
# ), CMake will use that automatically; FetchContent below is only a
# fallback so the project still builds with nothing preinstalled.
find_package(glfw3 QUIET)
if(NOT glfw3_FOUND)
  set(GLFW_BUILD_EXAMPLES OFF CACHE BOOL "" FORCE)
  set(GLFW_BUILD_TESTS OFF CACHE BOOL "" FORCE)
  set(GLFW_BUILD_DOCS OFF CACHE BOOL "" FORCE)
  set(GLFW_INSTALL OFF CACHE BOOL "" FORCE)
  FetchContent_Declare(
    glfw
    GIT_REPOSITORY https://github.com/glfw/glfw.git
    GIT_TAG 3.4
    GIT_SHALLOW TRUE
  )
  FetchContent_MakeAvailable(glfw)
endif()

# ---- Dear ImGui (core + GLFW/OpenGL2 backends) ----
FetchContent_Declare(
  imgui
  GIT_REPOSITORY https://github.com/ocornut/imgui.git
  GIT_TAG v1.91.0
  GIT_SHALLOW TRUE
)
FetchContent_MakeAvailable(imgui)

add_library(imgui STATIC
  ${imgui_SOURCE_DIR}/imgui.cpp
  ${imgui_SOURCE_DIR}/imgui_draw.cpp
  ${imgui_SOURCE_DIR}/imgui_tables.cpp
  ${imgui_SOURCE_DIR}/imgui_widgets.cpp
  ${imgui_SOURCE_DIR}/imgui_demo.cpp
  ${imgui_SOURCE_DIR}/backends/imgui_impl_glfw.cpp
  ${imgui_SOURCE_DIR}/backends/imgui_impl_opengl2.cpp
)
target_include_directories(imgui PUBLIC
  ${imgui_SOURCE_DIR}
  ${imgui_SOURCE_DIR}/backends
)
target_link_libraries(imgui PUBLIC glfw)

find_package(OpenGL REQUIRED)

add_executable(isosurface_lab
  src/main.cpp
)
target_link_libraries(isosurface_lab PRIVATE imgui glfw OpenGL::GL)

if(APPLE)
  target_link_libraries(isosurface_lab PRIVATE "-framework Cocoa" "-framework IOKit" "-framework CoreVideo")
  # We intentionally use the legacy/compatibility OpenGL API (glBegin/glEnd,
  # the fixed-function matrix stack) so the app needs no GL function loader.
  # Apple deprecated that API in 10.14 but it still works; silence the
  # resulting warning spam rather than pretend otherwise.
  target_compile_definitions(isosurface_lab PRIVATE GL_SILENCE_DEPRECATION)
endif()
if(UNIX AND NOT APPLE)
  find_package(Threads REQUIRED)
  target_link_libraries(isosurface_lab PRIVATE Threads::Threads)
endif()

# ---- correctness tests (no GL/windowing deps: parser + isosurface math only) ----
# Checks: a sphere's extracted vertices lie on the sphere with the correct
# analytic (AD) radial gradient, the cuspidal family x^3-ky has zero
# gradient everywhere at k=0 and nonzero off it, the mesh is watertight
# (every edge shared by exactly 2 triangles), and parser error handling.
enable_testing()
add_executable(isosurface_lab_tests EXCLUDE_FROM_ALL tests/test_isosurface.cpp)
target_include_directories(isosurface_lab_tests PRIVATE src)
add_test(NAME isosurface_correctness COMMAND isosurface_lab_tests)
