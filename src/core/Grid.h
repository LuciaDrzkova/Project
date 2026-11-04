#pragma once

// Single definition of the sampling lattice shared by the 3D extractor, the 2D
// slice exporter and the viewer overlays. Grid node i of an N-cell axis over
// [-h, h] sits at -h + 2h*i/N. Computing it this way (rather than accumulating
// i*step) keeps exactly representable nodes such as 0, +/-h/2 and +/-h exact,
// so every component evaluates F at bit-identical coordinates.

#include <algorithm>

inline double gridCoordinate(int index, int resolution, double halfExtent) {
    return -halfExtent + (2.0 * halfExtent * static_cast<double>(index) /
                          static_cast<double>(std::max(1, resolution)));
}
