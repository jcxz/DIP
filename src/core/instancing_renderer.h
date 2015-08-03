#ifndef INSTANCING_RENDERER_H
#define INSTANCING_RENDERER_H

#include "core/base_renderer.h"
#include "utils/geom.h"


class InstancingRenderer : public BaseRenderer
{
  public:
    InstancingRenderer(void)
      : BaseRenderer()
      , m_prog_particle_colors()
      , m_prog_uniform_color()
      , m_particle_geom()
    { }

  protected:
    virtual bool reset_impl(int /* w */, int /* h */) override;
    virtual void render_impl(const QQuaternion & rotation,
                             const QVector3D & scale,
                             const QVector3D & translation,
                             const QQuaternion & /* camera_rotation */,
                             GLuint vbo_positions,
                             GLuint vbo_colors,
                             size_t part_cnt,
                             const QVector3D & /* volume_size */) override;

  private:
    // OpenGL shaders
    QOpenGLShaderProgram m_prog_particle_colors;
    QOpenGLShaderProgram m_prog_uniform_color;

    // vertex buffers with sphere geometry (or the geometry of objects that will represent particles)
    utils::geom::Model m_particle_geom;
};

#endif // INSTANCING_RENDERER_H
