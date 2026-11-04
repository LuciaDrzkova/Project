#pragma once

// The viewer's single color ramp (dark purple -> blue -> green -> yellow -> red),
// used for |grad F| on the surface, gradient glyphs and slice-index contours.

#include <pmp/visualization/gl.h>
#include <pmp/types.h>

#include <algorithm>

namespace iso_app {

inline pmp::Color gradient_color(float t)
{
    t = std::clamp(t, 0.0f, 1.0f);
    static constexpr float stops[5][3] = {
        {0.15f, 0.05f, 0.45f},
        {0.10f, 0.45f, 0.75f},
        {0.25f, 0.75f, 0.35f},
        {0.95f, 0.80f, 0.15f},
        {0.90f, 0.15f, 0.15f}};
    const float scaled = t * 4.0f;
    const int i0 = std::min(3, static_cast<int>(scaled));
    const float f = scaled - static_cast<float>(i0);
    return pmp::Color(stops[i0][0] + f * (stops[i0 + 1][0] - stops[i0][0]),
                      stops[i0][1] + f * (stops[i0 + 1][1] - stops[i0][1]),
                      stops[i0][2] + f * (stops[i0 + 1][2] - stops[i0][2]));
}

} // namespace iso_app
