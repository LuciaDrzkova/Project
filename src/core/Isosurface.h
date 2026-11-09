// Isosurface.h
//
// Extracts the zero level set {F(x,y,z;k) = 0} of a user-supplied scalar
// field on a regular grid using the classic 256-case Marching Cubes
// algorithm (Lorensen & Cline 1987). Each grid cell (cube) has 8 corners;
// the sign of F at each corner (inside/outside) selects one of 256 cases,
// and a precomputed triangle table (kTriTable) says which of the cube's 12
// edges to connect into triangles for that case.
//
// Pipeline (Isosurface::build):
//   1. sampleField  - F is evaluated (value only, by the compiled ExprProgram,
//                     row by row and in parallel for larger grids) at every
//                     node -h + 2h*i/N of the lattice.
//   2. zero tie-break - nodes on F=0, exactly or within rounding noise, get a
//                     deterministic positive bias (see applyZeroTieBreak).
//   3. extract      - marching cubes. Every cut grid edge is interpolated
//                     exactly once, always from its lower to its upper node,
//                     and cached in a rolling two-plane edge table, so cells
//                     sharing an edge share the identical vertex. The vertex
//                     gradient/normal comes from re-evaluating F with Dual3
//                     (exact forward-mode AD) at the interpolated point.
//
// Known trade-off: a handful of the 256 cases are "ambiguous" on a cube face
// (two diagonally-opposite corners inside, the other two outside) where the
// standard table's choice of diagonal can, in principle, leave a crack
// against a neighbouring cube. This is rare for smooth, well-resolved fields
// and is not resolved here (no asymptotic decider); the watertightness tests
// in tests/test_isosurface.cpp are the check that would catch it.
#pragma once
#include "Expr.h"
#include "ExprProgram.h"
#include "Grid.h"
#include "Parallel.h"

#include <algorithm>
#include <array>
#include <chrono>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

struct SurfaceVertex {
    float px, py, pz;      // position
    float nx, ny, nz;      // unit normal = normalized gradient
    float gradMag;         // |grad F| at this point (un-normalized)
};

struct SkippedTriangleGeometry {
    float ax, ay, az, bx, by, bz, cx, cy, cz;
};

struct IsosurfaceResult {
    // Non-indexed triangle soup: three corners per triangle.
    std::vector<SurfaceVertex> vertices;
    // Per corner (parallel to `vertices`): the id of the grid-edge crossing the
    // corner came from. Equal ids mean the same welded vertex; ids are dense in
    // [0, uniqueVertexCount). Empty for hand-assembled results.
    std::vector<std::uint32_t> vertexIds;
    std::size_t uniqueVertexCount = 0;
    double buildMillis = 0.0;
    double maxGradMag = 1e-9; // for color normalization
    std::size_t skippedDegenerateTriangles = 0; // triangles dropped by the degenerate-area guard
    std::size_t exactZeroGridSamples = 0;       // nodes on F=0 (exactly or numerically) before tie-breaking
    std::vector<SkippedTriangleGeometry> skippedDegenerateGeometry; // positions of dropped triangle candidates
};

// Cube corner numbering/positions (standard Marching Cubes convention,
// matching the kTriTable/kEdgeCorners data below):
//   0:(0,0,0) 1:(1,0,0) 2:(1,1,0) 3:(0,1,0)
//   4:(0,0,1) 5:(1,0,1) 6:(1,1,1) 7:(0,1,1)
inline constexpr int kCubeCornerOffset[8][3] = {
    {0,0,0}, {1,0,0}, {1,1,0}, {0,1,0},
    {0,0,1}, {1,0,1}, {1,1,1}, {0,1,1},
};

// The 12 cube edges, each given as the pair of corner indices it connects
// (using the corner numbering above). Edge i is where kTriTable's entries
// of value i refer to.
inline constexpr int kEdgeCorners[12][2] = {
    {0,1}, {1,2}, {2,3}, {3,0}, {4,5}, {5,6}, {6,7}, {7,4}, {0,4}, {1,5}, {2,6}, {3,7}
};

// Classic Marching Cubes triangle table: for each of the 256 possible
// inside/outside corner configurations (bit i set <=> corner i is
// "inside", i.e. F<0), up to 5 triangles (15 edge indices, -1 terminated)
// describing how to connect edge-crossing points into triangles.
//
// This is the well-known public-domain table originally attributed to
// Cory Bloyd / Paul Bourke's "Polygonising a scalar field" article. This
// copy was fetched from libigl (github.com/libigl/libigl,
// include/igl/marching_cubes_tables.h, MIT licensed) via raw.githubusercontent.com
// and mechanically reformatted (not hand-transcribed) to avoid transcription
// errors, and verified by the topology tests in tests/test_isosurface.cpp.
inline constexpr std::int8_t kTriTable[256][16] = {
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 1, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 8, 3, 9, 8, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 3, 1, 2,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 9, 2,10, 0, 2, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 2, 8, 3, 2,10, 8,10, 9, 8,-1,-1,-1,-1,-1,-1,-1},
    { 3,11, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0,11, 2, 8,11, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 9, 0, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1,11, 2, 1, 9,11, 9, 8,11,-1,-1,-1,-1,-1,-1,-1},
    { 3,10, 1,11,10, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0,10, 1, 0, 8,10, 8,11,10,-1,-1,-1,-1,-1,-1,-1},
    { 3, 9, 0, 3,11, 9,11,10, 9,-1,-1,-1,-1,-1,-1,-1},
    { 9, 8,10,10, 8,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 3, 0, 7, 3, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 1, 9, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 1, 9, 4, 7, 1, 7, 3, 1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,10, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 3, 4, 7, 3, 0, 4, 1, 2,10,-1,-1,-1,-1,-1,-1,-1},
    { 9, 2,10, 9, 0, 2, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1},
    { 2,10, 9, 2, 9, 7, 2, 7, 3, 7, 9, 4,-1,-1,-1,-1},
    { 8, 4, 7, 3,11, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {11, 4, 7,11, 2, 4, 2, 0, 4,-1,-1,-1,-1,-1,-1,-1},
    { 9, 0, 1, 8, 4, 7, 2, 3,11,-1,-1,-1,-1,-1,-1,-1},
    { 4, 7,11, 9, 4,11, 9,11, 2, 9, 2, 1,-1,-1,-1,-1},
    { 3,10, 1, 3,11,10, 7, 8, 4,-1,-1,-1,-1,-1,-1,-1},
    { 1,11,10, 1, 4,11, 1, 0, 4, 7,11, 4,-1,-1,-1,-1},
    { 4, 7, 8, 9, 0,11, 9,11,10,11, 0, 3,-1,-1,-1,-1},
    { 4, 7,11, 4,11, 9, 9,11,10,-1,-1,-1,-1,-1,-1,-1},
    { 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 9, 5, 4, 0, 8, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 5, 4, 1, 5, 0,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 8, 5, 4, 8, 3, 5, 3, 1, 5,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,10, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 3, 0, 8, 1, 2,10, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1},
    { 5, 2,10, 5, 4, 2, 4, 0, 2,-1,-1,-1,-1,-1,-1,-1},
    { 2,10, 5, 3, 2, 5, 3, 5, 4, 3, 4, 8,-1,-1,-1,-1},
    { 9, 5, 4, 2, 3,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0,11, 2, 0, 8,11, 4, 9, 5,-1,-1,-1,-1,-1,-1,-1},
    { 0, 5, 4, 0, 1, 5, 2, 3,11,-1,-1,-1,-1,-1,-1,-1},
    { 2, 1, 5, 2, 5, 8, 2, 8,11, 4, 8, 5,-1,-1,-1,-1},
    {10, 3,11,10, 1, 3, 9, 5, 4,-1,-1,-1,-1,-1,-1,-1},
    { 4, 9, 5, 0, 8, 1, 8,10, 1, 8,11,10,-1,-1,-1,-1},
    { 5, 4, 0, 5, 0,11, 5,11,10,11, 0, 3,-1,-1,-1,-1},
    { 5, 4, 8, 5, 8,10,10, 8,11,-1,-1,-1,-1,-1,-1,-1},
    { 9, 7, 8, 5, 7, 9,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 9, 3, 0, 9, 5, 3, 5, 7, 3,-1,-1,-1,-1,-1,-1,-1},
    { 0, 7, 8, 0, 1, 7, 1, 5, 7,-1,-1,-1,-1,-1,-1,-1},
    { 1, 5, 3, 3, 5, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 9, 7, 8, 9, 5, 7,10, 1, 2,-1,-1,-1,-1,-1,-1,-1},
    {10, 1, 2, 9, 5, 0, 5, 3, 0, 5, 7, 3,-1,-1,-1,-1},
    { 8, 0, 2, 8, 2, 5, 8, 5, 7,10, 5, 2,-1,-1,-1,-1},
    { 2,10, 5, 2, 5, 3, 3, 5, 7,-1,-1,-1,-1,-1,-1,-1},
    { 7, 9, 5, 7, 8, 9, 3,11, 2,-1,-1,-1,-1,-1,-1,-1},
    { 9, 5, 7, 9, 7, 2, 9, 2, 0, 2, 7,11,-1,-1,-1,-1},
    { 2, 3,11, 0, 1, 8, 1, 7, 8, 1, 5, 7,-1,-1,-1,-1},
    {11, 2, 1,11, 1, 7, 7, 1, 5,-1,-1,-1,-1,-1,-1,-1},
    { 9, 5, 8, 8, 5, 7,10, 1, 3,10, 3,11,-1,-1,-1,-1},
    { 5, 7, 0, 5, 0, 9, 7,11, 0, 1, 0,10,11,10, 0,-1},
    {11,10, 0,11, 0, 3,10, 5, 0, 8, 0, 7, 5, 7, 0,-1},
    {11,10, 5, 7,11, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 3, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 9, 0, 1, 5,10, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 8, 3, 1, 9, 8, 5,10, 6,-1,-1,-1,-1,-1,-1,-1},
    { 1, 6, 5, 2, 6, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 6, 5, 1, 2, 6, 3, 0, 8,-1,-1,-1,-1,-1,-1,-1},
    { 9, 6, 5, 9, 0, 6, 0, 2, 6,-1,-1,-1,-1,-1,-1,-1},
    { 5, 9, 8, 5, 8, 2, 5, 2, 6, 3, 2, 8,-1,-1,-1,-1},
    { 2, 3,11,10, 6, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {11, 0, 8,11, 2, 0,10, 6, 5,-1,-1,-1,-1,-1,-1,-1},
    { 0, 1, 9, 2, 3,11, 5,10, 6,-1,-1,-1,-1,-1,-1,-1},
    { 5,10, 6, 1, 9, 2, 9,11, 2, 9, 8,11,-1,-1,-1,-1},
    { 6, 3,11, 6, 5, 3, 5, 1, 3,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8,11, 0,11, 5, 0, 5, 1, 5,11, 6,-1,-1,-1,-1},
    { 3,11, 6, 0, 3, 6, 0, 6, 5, 0, 5, 9,-1,-1,-1,-1},
    { 6, 5, 9, 6, 9,11,11, 9, 8,-1,-1,-1,-1,-1,-1,-1},
    { 5,10, 6, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 3, 0, 4, 7, 3, 6, 5,10,-1,-1,-1,-1,-1,-1,-1},
    { 1, 9, 0, 5,10, 6, 8, 4, 7,-1,-1,-1,-1,-1,-1,-1},
    {10, 6, 5, 1, 9, 7, 1, 7, 3, 7, 9, 4,-1,-1,-1,-1},
    { 6, 1, 2, 6, 5, 1, 4, 7, 8,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2, 5, 5, 2, 6, 3, 0, 4, 3, 4, 7,-1,-1,-1,-1},
    { 8, 4, 7, 9, 0, 5, 0, 6, 5, 0, 2, 6,-1,-1,-1,-1},
    { 7, 3, 9, 7, 9, 4, 3, 2, 9, 5, 9, 6, 2, 6, 9,-1},
    { 3,11, 2, 7, 8, 4,10, 6, 5,-1,-1,-1,-1,-1,-1,-1},
    { 5,10, 6, 4, 7, 2, 4, 2, 0, 2, 7,11,-1,-1,-1,-1},
    { 0, 1, 9, 4, 7, 8, 2, 3,11, 5,10, 6,-1,-1,-1,-1},
    { 9, 2, 1, 9,11, 2, 9, 4,11, 7,11, 4, 5,10, 6,-1},
    { 8, 4, 7, 3,11, 5, 3, 5, 1, 5,11, 6,-1,-1,-1,-1},
    { 5, 1,11, 5,11, 6, 1, 0,11, 7,11, 4, 0, 4,11,-1},
    { 0, 5, 9, 0, 6, 5, 0, 3, 6,11, 6, 3, 8, 4, 7,-1},
    { 6, 5, 9, 6, 9,11, 4, 7, 9, 7,11, 9,-1,-1,-1,-1},
    {10, 4, 9, 6, 4,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4,10, 6, 4, 9,10, 0, 8, 3,-1,-1,-1,-1,-1,-1,-1},
    {10, 0, 1,10, 6, 0, 6, 4, 0,-1,-1,-1,-1,-1,-1,-1},
    { 8, 3, 1, 8, 1, 6, 8, 6, 4, 6, 1,10,-1,-1,-1,-1},
    { 1, 4, 9, 1, 2, 4, 2, 6, 4,-1,-1,-1,-1,-1,-1,-1},
    { 3, 0, 8, 1, 2, 9, 2, 4, 9, 2, 6, 4,-1,-1,-1,-1},
    { 0, 2, 4, 4, 2, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 8, 3, 2, 8, 2, 4, 4, 2, 6,-1,-1,-1,-1,-1,-1,-1},
    {10, 4, 9,10, 6, 4,11, 2, 3,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 2, 2, 8,11, 4, 9,10, 4,10, 6,-1,-1,-1,-1},
    { 3,11, 2, 0, 1, 6, 0, 6, 4, 6, 1,10,-1,-1,-1,-1},
    { 6, 4, 1, 6, 1,10, 4, 8, 1, 2, 1,11, 8,11, 1,-1},
    { 9, 6, 4, 9, 3, 6, 9, 1, 3,11, 6, 3,-1,-1,-1,-1},
    { 8,11, 1, 8, 1, 0,11, 6, 1, 9, 1, 4, 6, 4, 1,-1},
    { 3,11, 6, 3, 6, 0, 0, 6, 4,-1,-1,-1,-1,-1,-1,-1},
    { 6, 4, 8,11, 6, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 7,10, 6, 7, 8,10, 8, 9,10,-1,-1,-1,-1,-1,-1,-1},
    { 0, 7, 3, 0,10, 7, 0, 9,10, 6, 7,10,-1,-1,-1,-1},
    {10, 6, 7, 1,10, 7, 1, 7, 8, 1, 8, 0,-1,-1,-1,-1},
    {10, 6, 7,10, 7, 1, 1, 7, 3,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2, 6, 1, 6, 8, 1, 8, 9, 8, 6, 7,-1,-1,-1,-1},
    { 2, 6, 9, 2, 9, 1, 6, 7, 9, 0, 9, 3, 7, 3, 9,-1},
    { 7, 8, 0, 7, 0, 6, 6, 0, 2,-1,-1,-1,-1,-1,-1,-1},
    { 7, 3, 2, 6, 7, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 2, 3,11,10, 6, 8,10, 8, 9, 8, 6, 7,-1,-1,-1,-1},
    { 2, 0, 7, 2, 7,11, 0, 9, 7, 6, 7,10, 9,10, 7,-1},
    { 1, 8, 0, 1, 7, 8, 1,10, 7, 6, 7,10, 2, 3,11,-1},
    {11, 2, 1,11, 1, 7,10, 6, 1, 6, 7, 1,-1,-1,-1,-1},
    { 8, 9, 6, 8, 6, 7, 9, 1, 6,11, 6, 3, 1, 3, 6,-1},
    { 0, 9, 1,11, 6, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 7, 8, 0, 7, 0, 6, 3,11, 0,11, 6, 0,-1,-1,-1,-1},
    { 7,11, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 7, 6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 3, 0, 8,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 1, 9,11, 7, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 8, 1, 9, 8, 3, 1,11, 7, 6,-1,-1,-1,-1,-1,-1,-1},
    {10, 1, 2, 6,11, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,10, 3, 0, 8, 6,11, 7,-1,-1,-1,-1,-1,-1,-1},
    { 2, 9, 0, 2,10, 9, 6,11, 7,-1,-1,-1,-1,-1,-1,-1},
    { 6,11, 7, 2,10, 3,10, 8, 3,10, 9, 8,-1,-1,-1,-1},
    { 7, 2, 3, 6, 2, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 7, 0, 8, 7, 6, 0, 6, 2, 0,-1,-1,-1,-1,-1,-1,-1},
    { 2, 7, 6, 2, 3, 7, 0, 1, 9,-1,-1,-1,-1,-1,-1,-1},
    { 1, 6, 2, 1, 8, 6, 1, 9, 8, 8, 7, 6,-1,-1,-1,-1},
    {10, 7, 6,10, 1, 7, 1, 3, 7,-1,-1,-1,-1,-1,-1,-1},
    {10, 7, 6, 1, 7,10, 1, 8, 7, 1, 0, 8,-1,-1,-1,-1},
    { 0, 3, 7, 0, 7,10, 0,10, 9, 6,10, 7,-1,-1,-1,-1},
    { 7, 6,10, 7,10, 8, 8,10, 9,-1,-1,-1,-1,-1,-1,-1},
    { 6, 8, 4,11, 8, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 3, 6,11, 3, 0, 6, 0, 4, 6,-1,-1,-1,-1,-1,-1,-1},
    { 8, 6,11, 8, 4, 6, 9, 0, 1,-1,-1,-1,-1,-1,-1,-1},
    { 9, 4, 6, 9, 6, 3, 9, 3, 1,11, 3, 6,-1,-1,-1,-1},
    { 6, 8, 4, 6,11, 8, 2,10, 1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,10, 3, 0,11, 0, 6,11, 0, 4, 6,-1,-1,-1,-1},
    { 4,11, 8, 4, 6,11, 0, 2, 9, 2,10, 9,-1,-1,-1,-1},
    {10, 9, 3,10, 3, 2, 9, 4, 3,11, 3, 6, 4, 6, 3,-1},
    { 8, 2, 3, 8, 4, 2, 4, 6, 2,-1,-1,-1,-1,-1,-1,-1},
    { 0, 4, 2, 4, 6, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 9, 0, 2, 3, 4, 2, 4, 6, 4, 3, 8,-1,-1,-1,-1},
    { 1, 9, 4, 1, 4, 2, 2, 4, 6,-1,-1,-1,-1,-1,-1,-1},
    { 8, 1, 3, 8, 6, 1, 8, 4, 6, 6,10, 1,-1,-1,-1,-1},
    {10, 1, 0,10, 0, 6, 6, 0, 4,-1,-1,-1,-1,-1,-1,-1},
    { 4, 6, 3, 4, 3, 8, 6,10, 3, 0, 3, 9,10, 9, 3,-1},
    {10, 9, 4, 6,10, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 9, 5, 7, 6,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 3, 4, 9, 5,11, 7, 6,-1,-1,-1,-1,-1,-1,-1},
    { 5, 0, 1, 5, 4, 0, 7, 6,11,-1,-1,-1,-1,-1,-1,-1},
    {11, 7, 6, 8, 3, 4, 3, 5, 4, 3, 1, 5,-1,-1,-1,-1},
    { 9, 5, 4,10, 1, 2, 7, 6,11,-1,-1,-1,-1,-1,-1,-1},
    { 6,11, 7, 1, 2,10, 0, 8, 3, 4, 9, 5,-1,-1,-1,-1},
    { 7, 6,11, 5, 4,10, 4, 2,10, 4, 0, 2,-1,-1,-1,-1},
    { 3, 4, 8, 3, 5, 4, 3, 2, 5,10, 5, 2,11, 7, 6,-1},
    { 7, 2, 3, 7, 6, 2, 5, 4, 9,-1,-1,-1,-1,-1,-1,-1},
    { 9, 5, 4, 0, 8, 6, 0, 6, 2, 6, 8, 7,-1,-1,-1,-1},
    { 3, 6, 2, 3, 7, 6, 1, 5, 0, 5, 4, 0,-1,-1,-1,-1},
    { 6, 2, 8, 6, 8, 7, 2, 1, 8, 4, 8, 5, 1, 5, 8,-1},
    { 9, 5, 4,10, 1, 6, 1, 7, 6, 1, 3, 7,-1,-1,-1,-1},
    { 1, 6,10, 1, 7, 6, 1, 0, 7, 8, 7, 0, 9, 5, 4,-1},
    { 4, 0,10, 4,10, 5, 0, 3,10, 6,10, 7, 3, 7,10,-1},
    { 7, 6,10, 7,10, 8, 5, 4,10, 4, 8,10,-1,-1,-1,-1},
    { 6, 9, 5, 6,11, 9,11, 8, 9,-1,-1,-1,-1,-1,-1,-1},
    { 3, 6,11, 0, 6, 3, 0, 5, 6, 0, 9, 5,-1,-1,-1,-1},
    { 0,11, 8, 0, 5,11, 0, 1, 5, 5, 6,11,-1,-1,-1,-1},
    { 6,11, 3, 6, 3, 5, 5, 3, 1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,10, 9, 5,11, 9,11, 8,11, 5, 6,-1,-1,-1,-1},
    { 0,11, 3, 0, 6,11, 0, 9, 6, 5, 6, 9, 1, 2,10,-1},
    {11, 8, 5,11, 5, 6, 8, 0, 5,10, 5, 2, 0, 2, 5,-1},
    { 6,11, 3, 6, 3, 5, 2,10, 3,10, 5, 3,-1,-1,-1,-1},
    { 5, 8, 9, 5, 2, 8, 5, 6, 2, 3, 8, 2,-1,-1,-1,-1},
    { 9, 5, 6, 9, 6, 0, 0, 6, 2,-1,-1,-1,-1,-1,-1,-1},
    { 1, 5, 8, 1, 8, 0, 5, 6, 8, 3, 8, 2, 6, 2, 8,-1},
    { 1, 5, 6, 2, 1, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 3, 6, 1, 6,10, 3, 8, 6, 5, 6, 9, 8, 9, 6,-1},
    {10, 1, 0,10, 0, 6, 9, 5, 0, 5, 6, 0,-1,-1,-1,-1},
    { 0, 3, 8, 5, 6,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {10, 5, 6,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {11, 5,10, 7, 5,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {11, 5,10,11, 7, 5, 8, 3, 0,-1,-1,-1,-1,-1,-1,-1},
    { 5,11, 7, 5,10,11, 1, 9, 0,-1,-1,-1,-1,-1,-1,-1},
    {10, 7, 5,10,11, 7, 9, 8, 1, 8, 3, 1,-1,-1,-1,-1},
    {11, 1, 2,11, 7, 1, 7, 5, 1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 3, 1, 2, 7, 1, 7, 5, 7, 2,11,-1,-1,-1,-1},
    { 9, 7, 5, 9, 2, 7, 9, 0, 2, 2,11, 7,-1,-1,-1,-1},
    { 7, 5, 2, 7, 2,11, 5, 9, 2, 3, 2, 8, 9, 8, 2,-1},
    { 2, 5,10, 2, 3, 5, 3, 7, 5,-1,-1,-1,-1,-1,-1,-1},
    { 8, 2, 0, 8, 5, 2, 8, 7, 5,10, 2, 5,-1,-1,-1,-1},
    { 9, 0, 1, 5,10, 3, 5, 3, 7, 3,10, 2,-1,-1,-1,-1},
    { 9, 8, 2, 9, 2, 1, 8, 7, 2,10, 2, 5, 7, 5, 2,-1},
    { 1, 3, 5, 3, 7, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 7, 0, 7, 1, 1, 7, 5,-1,-1,-1,-1,-1,-1,-1},
    { 9, 0, 3, 9, 3, 5, 5, 3, 7,-1,-1,-1,-1,-1,-1,-1},
    { 9, 8, 7, 5, 9, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 5, 8, 4, 5,10, 8,10,11, 8,-1,-1,-1,-1,-1,-1,-1},
    { 5, 0, 4, 5,11, 0, 5,10,11,11, 3, 0,-1,-1,-1,-1},
    { 0, 1, 9, 8, 4,10, 8,10,11,10, 4, 5,-1,-1,-1,-1},
    {10,11, 4,10, 4, 5,11, 3, 4, 9, 4, 1, 3, 1, 4,-1},
    { 2, 5, 1, 2, 8, 5, 2,11, 8, 4, 5, 8,-1,-1,-1,-1},
    { 0, 4,11, 0,11, 3, 4, 5,11, 2,11, 1, 5, 1,11,-1},
    { 0, 2, 5, 0, 5, 9, 2,11, 5, 4, 5, 8,11, 8, 5,-1},
    { 9, 4, 5, 2,11, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 2, 5,10, 3, 5, 2, 3, 4, 5, 3, 8, 4,-1,-1,-1,-1},
    { 5,10, 2, 5, 2, 4, 4, 2, 0,-1,-1,-1,-1,-1,-1,-1},
    { 3,10, 2, 3, 5,10, 3, 8, 5, 4, 5, 8, 0, 1, 9,-1},
    { 5,10, 2, 5, 2, 4, 1, 9, 2, 9, 4, 2,-1,-1,-1,-1},
    { 8, 4, 5, 8, 5, 3, 3, 5, 1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 4, 5, 1, 0, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 8, 4, 5, 8, 5, 3, 9, 0, 5, 0, 3, 5,-1,-1,-1,-1},
    { 9, 4, 5,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4,11, 7, 4, 9,11, 9,10,11,-1,-1,-1,-1,-1,-1,-1},
    { 0, 8, 3, 4, 9, 7, 9,11, 7, 9,10,11,-1,-1,-1,-1},
    { 1,10,11, 1,11, 4, 1, 4, 0, 7, 4,11,-1,-1,-1,-1},
    { 3, 1, 4, 3, 4, 8, 1,10, 4, 7, 4,11,10,11, 4,-1},
    { 4,11, 7, 9,11, 4, 9, 2,11, 9, 1, 2,-1,-1,-1,-1},
    { 9, 7, 4, 9,11, 7, 9, 1,11, 2,11, 1, 0, 8, 3,-1},
    {11, 7, 4,11, 4, 2, 2, 4, 0,-1,-1,-1,-1,-1,-1,-1},
    {11, 7, 4,11, 4, 2, 8, 3, 4, 3, 2, 4,-1,-1,-1,-1},
    { 2, 9,10, 2, 7, 9, 2, 3, 7, 7, 4, 9,-1,-1,-1,-1},
    { 9,10, 7, 9, 7, 4,10, 2, 7, 8, 7, 0, 2, 0, 7,-1},
    { 3, 7,10, 3,10, 2, 7, 4,10, 1,10, 0, 4, 0,10,-1},
    { 1,10, 2, 8, 7, 4,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 9, 1, 4, 1, 7, 7, 1, 3,-1,-1,-1,-1,-1,-1,-1},
    { 4, 9, 1, 4, 1, 7, 0, 8, 1, 8, 7, 1,-1,-1,-1,-1},
    { 4, 0, 3, 7, 4, 3,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 4, 8, 7,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 9,10, 8,10,11, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 3, 0, 9, 3, 9,11,11, 9,10,-1,-1,-1,-1,-1,-1,-1},
    { 0, 1,10, 0,10, 8, 8,10,11,-1,-1,-1,-1,-1,-1,-1},
    { 3, 1,10,11, 3,10,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 2,11, 1,11, 9, 9,11, 8,-1,-1,-1,-1,-1,-1,-1},
    { 3, 0, 9, 3, 9,11, 1, 2, 9, 2,11, 9,-1,-1,-1,-1},
    { 0, 2,11, 8, 0,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 3, 2,11,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 2, 3, 8, 2, 8,10,10, 8, 9,-1,-1,-1,-1,-1,-1,-1},
    { 9,10, 2, 0, 9, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 2, 3, 8, 2, 8,10, 0, 1, 8, 1,10, 8,-1,-1,-1,-1},
    { 1,10, 2,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 1, 3, 8, 9, 1, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 9, 1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    { 0, 3, 8,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
    {-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1,-1},
};

class Isosurface {
public:
    // Positive bias given to an on-surface node, as a fraction of the smallest
    // finite, nonzero, off-surface face-neighbour |F| (so it stays below every
    // such neighbour magnitude: at most a small fraction of a cell).
    static constexpr double kZeroBiasFraction = 5e-2;
    // A node is "on the surface" when F is exactly 0, or when |F| is at most
    // this fraction of the largest finite face-neighbour |F|. Such a value is
    // rounding noise (e.g. x^2+y^2+z^2-0.25 at the lattice point (1/3,1/3,1/6)
    // evaluates to ~1e-17): its edge crossings would all collapse onto the
    // node and emit zero-area slivers, exactly like a true zero.
    static constexpr double kNumericalZeroTolerance = 1e-10;
    // Largest grid resolution accepted by the analyses (cells per axis). The
    // extractor itself has no limit; 512^3 cells need about 1 GB of samples.
    static constexpr int kMaxAnalysisResolution = 512;
    // Below this many samples, threading costs more than it saves.
    static constexpr std::size_t kParallelSampleThreshold = 32768;

    // domainHalfExtent: grid spans [-h,h]^3. res: number of cells per axis.
    static IsosurfaceResult build(const Node* ast, double k, int res, double domainHalfExtent) {
        const auto t0 = std::chrono::steady_clock::now();
        res = std::max(res, 2);
        const ExprProgram program(ast);
        IsosurfaceResult result = extract(program, k, res, domainHalfExtent,
                                          sampleField(program, k, res, domainHalfExtent));
        result.buildMillis = elapsedMillis(t0);
        return result;
    }

    // F at every lattice node, x-major with z fastest: index (i*N + j)*N + l.
    // Each z-row is evaluated as one batch by the compiled program.
    static std::vector<double> sampleField(const ExprProgram& program, double k, int res, double h) {
        res = std::max(res, 2);
        const int N = res + 1;
        const std::size_t n = static_cast<std::size_t>(N);
        const std::vector<double> coords = axisCoordinates(res, h);
        std::vector<double> values(n * n * n);
        const unsigned threads = values.size() >= kParallelSampleThreshold ? 0u : 1u;
        iso_parallel::parallelFor(n, [&](std::size_t i) {
            std::vector<double> scratch;
            for (std::size_t j = 0; j < n; ++j)
                program.evaluate(ExprProgram::Lane::shared(coords[i]), ExprProgram::Lane::shared(coords[j]),
                                 ExprProgram::Lane::array(coords.data()), k, n,
                                 values.data() + (i * n + j) * n, scratch);
        }, threads);
        return values;
    }
    static std::vector<double> sampleField(const Node* ast, double k, int res, double h) {
        return sampleField(ExprProgram(ast), k, res, h);
    }

    // Marching cubes over caller-supplied lattice samples (layout as returned
    // by sampleField). The zero tie-break is applied first, exactly as in
    // build(); the expression and `k` are only used for the analytic vertex
    // gradients.
    static IsosurfaceResult extract(const Node* ast, double k, int res, double h,
                                    std::vector<double> values) {
        return extract(ExprProgram(ast), k, res, h, std::move(values));
    }
    static IsosurfaceResult extract(const ExprProgram& program, double k, int res, double h,
                                    std::vector<double> values) {
        const auto t0 = std::chrono::steady_clock::now();
        res = std::max(res, 2);
        const int N = res + 1;
        const std::size_t plane = static_cast<std::size_t>(N) * N;
        if (values.size() != plane * N)
            throw std::invalid_argument("Isosurface::extract: sample count does not match (res+1)^3.");

        IsosurfaceResult result;
        result.exactZeroGridSamples = applyZeroTieBreak(values, N);

        const std::vector<double> coords = axisCoordinates(res, h);
        const Dual3 kd(k, 0.0, 0.0, 0.0);
        const std::size_t stride[3] = {plane, static_cast<std::size_t>(N), 1};
        std::vector<SurfaceVertex> unique;

        // Interpolates the crossing on the edge from lattice node `lower` (index
        // i,j,l) to its +1 neighbour along `axis`, then evaluates the exact
        // gradient there. Each edge is computed once, so the result does not
        // depend on which neighbouring cell asked first.
        auto cutEdge = [&](std::size_t lower, int axis, const int (&node)[3]) {
            const double fa = values[lower];
            const double fb = values[lower + stride[axis]];
            const double t = fa / (fa - fb); // fa and fb lie on opposite sides of 0
            double p[3] = {coords[static_cast<std::size_t>(node[0])],
                           coords[static_cast<std::size_t>(node[1])],
                           coords[static_cast<std::size_t>(node[2])]};
            const double upper = coords[static_cast<std::size_t>(node[axis] + 1)];
            p[axis] += t * (upper - p[axis]);
            const Dual3 r = program.evaluate<Dual3>(Dual3(p[0], 1, 0, 0), Dual3(p[1], 0, 1, 0),
                                                    Dual3(p[2], 0, 0, 1), kd);
            SurfaceVertex v{};
            v.px = static_cast<float>(p[0]);
            v.py = static_cast<float>(p[1]);
            v.pz = static_cast<float>(p[2]);
            const double gm = std::sqrt(r.dx * r.dx + r.dy * r.dy + r.dz * r.dz);
            v.gradMag = static_cast<float>(gm);
            if (gm > 1e-12) {
                v.nx = static_cast<float>(r.dx / gm);
                v.ny = static_cast<float>(r.dy / gm);
                v.nz = static_cast<float>(r.dz / gm);
            } else {
                v.nx = 0.0f; v.ny = 0.0f; v.nz = 1.0f;
            }
            return v;
        };

        // Degenerate-triangle guard: reject any triangle whose float area is
        // ~0 right before emitting it (a value epsilon would underflow back to
        // the same float32 position). With the tie-break above this should be
        // rare; rejected candidates are counted and kept for diagnostics.
        const float step = static_cast<float>(2.0 * h / res);
        const float degenerateAreaEps = step * step * 1e-8f;
        auto emitTri = [&](std::uint32_t ia, std::uint32_t ib, std::uint32_t ic) {
            const SurfaceVertex& a = unique[ia];
            const SurfaceVertex& b = unique[ib];
            const SurfaceVertex& c = unique[ic];
            const float ex1 = b.px - a.px, ey1 = b.py - a.py, ez1 = b.pz - a.pz;
            const float ex2 = c.px - a.px, ey2 = c.py - a.py, ez2 = c.pz - a.pz;
            const float fnx = ey1 * ez2 - ez1 * ey2;
            const float fny = ez1 * ex2 - ex1 * ez2;
            const float fnz = ex1 * ey2 - ey1 * ex2;
            if (fnx * fnx + fny * fny + fnz * fnz < degenerateAreaEps * degenerateAreaEps) {
                ++result.skippedDegenerateTriangles;
                result.skippedDegenerateGeometry.push_back({a.px, a.py, a.pz, b.px, b.py, b.pz, c.px, c.py, c.pz});
                return;
            }
            // Orient each face to agree with the analytic gradient.
            const float dot = fnx * (a.nx + b.nx + c.nx) + fny * (a.ny + b.ny + c.ny) +
                              fnz * (a.nz + b.nz + c.nz);
            if (dot < 0.0f) std::swap(ib, ic);
            for (std::uint32_t id : {ia, ib, ic}) {
                result.vertices.push_back(unique[id]);
                result.vertexIds.push_back(id);
            }
            result.maxGradMag = std::max(result.maxGradMag,
                static_cast<double>(std::max({a.gradMag, b.gradMag, c.gradMag})));
        };

        // Rolling edge cache for the cells of one x-slab [i, i+1]:
        //   xEdges[j*N+l]              x-directed edge (i,j,l)-(i+1,j,l)
        //   lowPlane/highPlane[a*N*N+j*N+l]  y- (a=0) or z- (a=1) directed edge
        //                              starting at (i,j,l) / (i+1,j,l)
        // Moving to the next slab turns the high plane into the low plane.
        constexpr std::int32_t kNone = -1;
        std::vector<std::int32_t> xEdges(plane, kNone);
        std::vector<std::int32_t> lowPlane(2 * plane, kNone);
        std::vector<std::int32_t> highPlane(2 * plane, kNone);
        const auto& edges = edgeTable();

        for (int i = 0; i < res; ++i) {
            if (i > 0) {
                std::swap(lowPlane, highPlane);
                std::fill(highPlane.begin(), highPlane.end(), kNone);
                std::fill(xEdges.begin(), xEdges.end(), kNone);
            }
            for (int j = 0; j < res; ++j) {
                for (int l = 0; l < res; ++l) {
                    const std::size_t base = (static_cast<std::size_t>(i) * N + j) * N + l;
                    // bit c set <=> corner c is "inside" (F < 0).
                    int cubeIndex = 0;
                    for (int c = 0; c < 8; ++c) {
                        const std::size_t id = base + kCubeCornerOffset[c][0] * plane +
                                               kCubeCornerOffset[c][1] * stride[1] +
                                               static_cast<std::size_t>(kCubeCornerOffset[c][2]);
                        if (values[id] < 0.0) cubeIndex |= (1 << c);
                    }
                    if (cubeIndex == 0 || cubeIndex == 255) continue; // entirely outside / inside

                    auto edgeVertex = [&](int e) -> std::uint32_t {
                        const EdgeInfo& edge = edges[static_cast<std::size_t>(e)];
                        const int node[3] = {i + edge.lower[0], j + edge.lower[1], l + edge.lower[2]};
                        const std::size_t inPlane = static_cast<std::size_t>(node[1]) * N +
                                                    static_cast<std::size_t>(node[2]);
                        std::int32_t& slot = edge.axis == 0
                            ? xEdges[inPlane]
                            : (edge.lower[0] == 0 ? lowPlane : highPlane)[(edge.axis - 1) * plane + inPlane];
                        if (slot == kNone) {
                            const std::size_t lower = (static_cast<std::size_t>(node[0]) * N + node[1]) * N + node[2];
                            slot = static_cast<std::int32_t>(unique.size());
                            unique.push_back(cutEdge(lower, edge.axis, node));
                        }
                        return static_cast<std::uint32_t>(slot);
                    };

                    const std::int8_t* tri = kTriTable[cubeIndex];
                    for (int t = 0; tri[t] != -1; t += 3)
                        emitTri(edgeVertex(tri[t]), edgeVertex(tri[t + 1]), edgeVertex(tri[t + 2]));
                }
            }
        }

        result.uniqueVertexCount = unique.size();
        result.buildMillis = elapsedMillis(t0);
        return result;
    }

    // Finds lattice nodes that lie on F=0 (exactly, or within rounding noise;
    // see kNumericalZeroTolerance) and moves those whose edge crossings would
    // collapse deterministically to the positive side. Without this, several
    // cube edges meeting at such a node report (almost) the same crossing
    // point and the table emits zero-area triangles. The bias is computed
    // from the original values of off-surface neighbours only. This is a
    // local grid tie-break (a fraction of a cell), not a change of the
    // analytic field: vertex gradients and all residuals still use the
    // original expression. Returns the number of on-surface nodes found.
    static std::size_t applyZeroTieBreak(std::vector<double>& f, int N) {
        const std::size_t plane = static_cast<std::size_t>(N) * N;
        const std::size_t n = static_cast<std::size_t>(N);
        double maxAbs = 0.0;
        for (double v : f)
            if (std::isfinite(v)) maxAbs = std::max(maxAbs, std::fabs(v));

        auto forEachNeighbour = [&](std::size_t id, auto&& visit) {
            const std::size_t i = id / plane, j = (id / n) % n, l = id % n;
            if (i > 0) visit(id - plane);
            if (i + 1 < n) visit(id + plane);
            if (j > 0) visit(id - n);
            if (j + 1 < n) visit(id + n);
            if (l > 0) visit(id - 1);
            if (l + 1 < n) visit(id + 1);
        };

        // Pass 1: find on-surface nodes. The global prefilter is a necessary
        // condition (every local scale is <= maxAbs) and skips almost all nodes.
        std::vector<std::size_t> zeros;
        const double prefilter = kNumericalZeroTolerance * maxAbs;
        for (std::size_t id = 0; id < f.size(); ++id) {
            const double magnitude = std::fabs(f[id]);
            if (f[id] != 0.0) {
                if (!(magnitude <= prefilter)) continue; // also skips NaN
                double scale = 0.0;
                forEachNeighbour(id, [&](std::size_t nb) {
                    if (std::isfinite(f[nb])) scale = std::max(scale, std::fabs(f[nb]));
                });
                if (!(magnitude <= kNumericalZeroTolerance * scale)) continue;
            }
            zeros.push_back(id); // ascending, so binary_search works below
        }
        if (zeros.empty()) return 0;

        // Pass 2: compute every bias from the untouched values, then apply.
        // Only a node whose crossings would collapse is moved: that needs two
        // or more face neighbours on the opposite side (F<0 vs F>=0), since
        // each such neighbour cuts one edge exactly at the node. With at most
        // one, the node's single crossing is exact and harmless; leaving it
        // keeps grid-aligned zero sets (e.g. the plane x=0 of x^3 at k=0, or
        // the coordinate planes of x*y*z) exactly on their analytic locus.
        std::vector<double> biased(zeros.size());
        for (std::size_t z = 0; z < zeros.size(); ++z) {
            const std::size_t id = zeros[z];
            const bool inside = f[id] < 0.0;
            int opposite = 0;
            double localScale = std::numeric_limits<double>::infinity();
            forEachNeighbour(id, [&](std::size_t nb) {
                const double v = f[nb];
                if (!std::isfinite(v)) return;
                if ((v < 0.0) != inside) ++opposite;
                if (v != 0.0 && !std::binary_search(zeros.begin(), zeros.end(), nb))
                    localScale = std::min(localScale, std::fabs(v));
            });
            biased[z] = f[id];
            if (opposite < 2) continue;
            if (!std::isfinite(localScale)) localScale = maxAbs;
            if (localScale > 0.0 && std::isfinite(localScale)) {
                const double bias = localScale * kZeroBiasFraction;
                biased[z] = bias > 0.0 ? bias : std::nextafter(0.0, 1.0);
            }
        }
        for (std::size_t z = 0; z < zeros.size(); ++z) f[zeros[z]] = biased[z];
        return zeros.size();
    }

private:
    struct EdgeInfo {
        int axis;     // 0 = x, 1 = y, 2 = z
        int lower[3]; // cube-corner offset of the edge's lower endpoint
    };

    static const std::array<EdgeInfo, 12>& edgeTable() {
        static const std::array<EdgeInfo, 12> table = [] {
            std::array<EdgeInfo, 12> out{};
            for (int e = 0; e < 12; ++e) {
                const int* a = kCubeCornerOffset[kEdgeCorners[e][0]];
                const int* b = kCubeCornerOffset[kEdgeCorners[e][1]];
                EdgeInfo& info = out[static_cast<std::size_t>(e)];
                info.axis = a[0] != b[0] ? 0 : (a[1] != b[1] ? 1 : 2);
                for (int c = 0; c < 3; ++c) info.lower[c] = std::min(a[c], b[c]);
            }
            return out;
        }();
        return table;
    }

    static std::vector<double> axisCoordinates(int res, double h) {
        std::vector<double> coords(static_cast<std::size_t>(res) + 1);
        for (int i = 0; i <= res; ++i) coords[static_cast<std::size_t>(i)] = gridCoordinate(i, res, h);
        return coords;
    }

    static double elapsedMillis(std::chrono::steady_clock::time_point start) {
        return std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
    }
};
