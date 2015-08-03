#ifndef POINT_SPRITE_RENDERER_H
#define POINT_SPRITE_RENDERER_H

#include "core/base_renderer.h"


class PointSpriteRenderer : public BaseRenderer
{
  public:
    PointSpriteRenderer(void)
      : BaseRenderer()
      , m_prog()
      , m_vao()
      , m_w(0.0f)
      , m_h(0.0f)
    {
      if (!init())
        throw std::runtime_error("Failed to initialize PointSpriteRenderer");
    }

  protected:
    virtual bool resize_impl(int w, int h) override;
    virtual bool reset_impl(int w, int h) override;
    virtual void render_impl(const QQuaternion & rotation,
                             const QVector3D & scale,
                             const QVector3D & translation,
                             const QQuaternion & /* camera_rotation */,
                             GLuint vbo_positions,
                             GLuint vbo_colors,
                             size_t part_cnt,
                             const QVector3D & /* volume_size */) override;

  private:
    bool init(void);

  private:
    // OpenGL shaders
    QOpenGLShaderProgram m_prog;
    QOpenGLVertexArrayObject m_vao;
    GLfloat m_w;
    GLfloat m_h;
};

#endif // POINT_SPRITE_RENDERER_H
