#pragma once

// Draws the viewer's helper geometry (sampling grid, selected cell, gradient
// glyphs, markers, slice contours) on top of the PMP surface. Owns every raw
// OpenGL object it uses: one VAO/VBO per layer plus a flat-color shader.
//
// Usage per frame: set() any layer whose geometry changed (uploaded lazily),
// then begin(), draw() the visible layers, end().

#include <pmp/mat_vec.h>
#include <pmp/types.h>
#include <pmp/visualization/gl.h>
#include <pmp/visualization/shader.h>

#include <array>
#include <cstddef>
#include <vector>

namespace iso_app {

struct OverlayVertex
{
    float x, y, z;
    float r, g, b;
};

inline OverlayVertex overlay_vertex(double x, double y, double z, float r, float g, float b)
{
    return {static_cast<float>(x), static_cast<float>(y), static_cast<float>(z), r, g, b};
}

inline OverlayVertex overlay_vertex(double x, double y, double z, const pmp::Color& c)
{
    return overlay_vertex(x, y, z, c[0], c[1], c[2]);
}

class OverlayRenderer
{
public:
    enum Layer : std::size_t
    {
        Grid,
        CellEdges,
        CellCorners,
        GradientGlyphs,
        PointOfInterest,
        SelectedGridVertex,
        SliceContours,
        SlicePlane,
        SliceSelectedPoint,
        LayerCount
    };

    enum class Depth
    {
        Test,      // regular depth test
        LessEqual, // for geometry lying exactly on the surface
        Off        // always on top
    };

    struct Style
    {
        GLenum primitive = GL_LINES;
        float alpha = 1.0f;
        float size = 1.0f; // line width or point size, clamped to driver support
        Depth depth = Depth::Test;
    };

    OverlayRenderer() = default;
    OverlayRenderer(const OverlayRenderer&) = delete;
    OverlayRenderer& operator=(const OverlayRenderer&) = delete;

    bool initialized() const { return initialized_; }
    void initialize(); // requires a current GL context
    void release();    // requires a current GL context

    void set(Layer layer, std::vector<OverlayVertex> vertices);
    bool empty(Layer layer) const { return layers_[layer].vertices.empty(); }

    void begin(const pmp::mat4& projection, const pmp::mat4& modelview);
    void draw(Layer layer, const Style& style);
    void end();

private:
    struct Buffer
    {
        GLuint vao = 0;
        GLuint vbo = 0;
        std::vector<OverlayVertex> vertices;
        bool dirty = false;
    };

    float line_width(float requested) const;
    float point_size(float requested) const;

    std::array<Buffer, LayerCount> layers_{};
    pmp::Shader shader_;
    bool initialized_ = false;
    // Many Core Profile drivers (macOS, and Mesa on Linux) report a generous
    // GL_ALIASED_LINE_WIDTH_RANGE but then raise GL_INVALID_VALUE for any
    // glLineWidth other than 1.0. initialize() therefore actually tries a
    // wide line once and falls back to 1px if the driver rejects it.
    bool wide_lines_supported_ = false;
    std::array<GLfloat, 2> line_width_range_{1.0f, 1.0f};
    std::array<GLfloat, 2> point_size_range_{1.0f, 1.0f};
};

} // namespace iso_app
