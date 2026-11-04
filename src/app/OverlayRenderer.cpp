#include "OverlayRenderer.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <utility>

namespace iso_app {

namespace {

constexpr const char* kVertexShader = R"GLSL(
    #version 330 core
    layout(location = 0) in vec3 vertex;
    layout(location = 1) in vec3 vertex_color;
    uniform mat4 projection_matrix;
    uniform mat4 modelview_matrix;
    uniform float point_size;
    out vec3 color;
    void main()
    {
        gl_Position = projection_matrix * modelview_matrix * vec4(vertex, 1.0);
        gl_PointSize = point_size;
        color = vertex_color;
    }
)GLSL";

constexpr const char* kFragmentShader = R"GLSL(
    #version 330 core
    in vec3 color;
    uniform float alpha;
    out vec4 frag_color;
    void main() { frag_color = vec4(color, alpha); }
)GLSL";

void sanitize_range(std::array<GLfloat, 2>& range)
{
    if (!std::isfinite(range[0]) || !std::isfinite(range[1]) || range[0] < 1.0f || range[1] < range[0])
        range = {1.0f, 1.0f};
}

} // namespace

void OverlayRenderer::initialize()
{
    if (initialized_)
        return;

    glGetFloatv(GL_ALIASED_LINE_WIDTH_RANGE, line_width_range_.data());
    glGetFloatv(GL_ALIASED_POINT_SIZE_RANGE, point_size_range_.data());
    // The range query alone is not trustworthy (see the header): clear any
    // pending error, try the widest width this viewer ever requests, and
    // check whether the driver really accepted it.
    while (glGetError() != GL_NO_ERROR) {}
    glLineWidth(std::clamp(3.5f, line_width_range_[0], line_width_range_[1]));
    wide_lines_supported_ = (glGetError() == GL_NO_ERROR);
    glLineWidth(1.0f);
    glGetError(); // a reset to 1.0 cannot meaningfully fail; discard
    sanitize_range(line_width_range_);
    sanitize_range(point_size_range_);

    shader_.source(kVertexShader, kFragmentShader);
    for (Buffer& layer : layers_)
    {
        glGenVertexArrays(1, &layer.vao);
        glGenBuffers(1, &layer.vbo);
        glBindVertexArray(layer.vao);
        glBindBuffer(GL_ARRAY_BUFFER, layer.vbo);
        glBufferData(GL_ARRAY_BUFFER, 0, nullptr, GL_DYNAMIC_DRAW);
        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(OverlayVertex),
                              reinterpret_cast<void*>(offsetof(OverlayVertex, x)));
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(OverlayVertex),
                              reinterpret_cast<void*>(offsetof(OverlayVertex, r)));
        glEnableVertexAttribArray(1);
        layer.dirty = true;
    }
    glBindVertexArray(0);
    initialized_ = true;
}

void OverlayRenderer::release()
{
    for (Buffer& layer : layers_)
    {
        if (layer.vbo) glDeleteBuffers(1, &layer.vbo);
        if (layer.vao) glDeleteVertexArrays(1, &layer.vao);
        layer.vbo = 0;
        layer.vao = 0;
    }
    initialized_ = false;
}

void OverlayRenderer::set(Layer layer, std::vector<OverlayVertex> vertices)
{
    layers_[layer].vertices = std::move(vertices);
    layers_[layer].dirty = true;
}

void OverlayRenderer::begin(const pmp::mat4& projection, const pmp::mat4& modelview)
{
    shader_.use();
    shader_.set_uniform("projection_matrix", projection);
    shader_.set_uniform("modelview_matrix", modelview);
    shader_.set_uniform("point_size", point_size(1.0f));
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glDepthMask(GL_FALSE);
    glEnable(GL_PROGRAM_POINT_SIZE);
}

void OverlayRenderer::draw(Layer layer, const Style& style)
{
    Buffer& buffer = layers_[layer];
    if (buffer.vertices.empty())
        return;
    if (buffer.dirty)
    {
        glBindBuffer(GL_ARRAY_BUFFER, buffer.vbo);
        glBufferData(GL_ARRAY_BUFFER, static_cast<GLsizeiptr>(buffer.vertices.size() * sizeof(OverlayVertex)),
                     buffer.vertices.data(), GL_DYNAMIC_DRAW);
        glBindBuffer(GL_ARRAY_BUFFER, 0);
        buffer.dirty = false;
    }

    shader_.set_uniform("alpha", style.alpha);
    if (style.primitive == GL_LINES)
        glLineWidth(line_width(style.size));
    else if (style.primitive == GL_POINTS)
        shader_.set_uniform("point_size", point_size(style.size));

    GLint previous_depth_func = GL_LESS;
    const GLboolean depth_was_enabled = glIsEnabled(GL_DEPTH_TEST);
    if (style.depth == Depth::LessEqual)
    {
        glGetIntegerv(GL_DEPTH_FUNC, &previous_depth_func);
        glDepthFunc(GL_LEQUAL);
    }
    else if (style.depth == Depth::Off)
    {
        glDisable(GL_DEPTH_TEST);
    }

    glBindVertexArray(buffer.vao);
    glDrawArrays(style.primitive, 0, static_cast<GLsizei>(buffer.vertices.size()));

    if (style.depth == Depth::LessEqual)
        glDepthFunc(static_cast<GLenum>(previous_depth_func));
    else if (style.depth == Depth::Off && depth_was_enabled)
        glEnable(GL_DEPTH_TEST);
}

void OverlayRenderer::end()
{
    glBindVertexArray(0);
    glDepthMask(GL_TRUE);
    glDisable(GL_BLEND);
    glDisable(GL_PROGRAM_POINT_SIZE);
    shader_.disable();
}

float OverlayRenderer::line_width(float requested) const
{
    if (!wide_lines_supported_)
        return 1.0f;
    return std::clamp(requested, line_width_range_[0], line_width_range_[1]);
}

float OverlayRenderer::point_size(float requested) const
{
    return std::clamp(requested, point_size_range_[0], point_size_range_[1]);
}

} // namespace iso_app
