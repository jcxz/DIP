#include "core/point_sprite_renderer.h"
#include "utils/ogl.h"
#include "utils/debug.h"

// include OpenCL definitions in crossplatform manner
//#include <boost/compute/cl.hpp>


//static constexpr float particle_radius = 10.0f; //5.0f;



bool PointSpriteRenderer::resize_impl(int w, int h)
{
  m_w = w;
  m_h = h;
  return true;
}


bool PointSpriteRenderer::reset_impl(int w, int h)
{
  m_w = w;
  m_h = h;
  return true;
}


void PointSpriteRenderer::render_impl(const QQuaternion & rotation,
                                      const QVector3D & scale,
                                      const QVector3D & translation,
                                      const QQuaternion & /* camera_rotation */,
                                      GLuint vbo_positions,
                                      GLuint vbo_colors,
                                      size_t part_cnt,
                                      const QVector3D & /* volume_size */)
{
  // render the particle system
  OGLF->glEnable(GL_DEPTH_TEST);

  m_vao.bind();

  GLuint shader_id = m_prog.programId();

  OGLF->glUseProgram(shader_id);

  // calculate scene transformation
  QMatrix4x4 mv;

  mv.translate(translation);
  mv.rotate(rotation);
  mv.scale(scale);

  //QMatrix4x4 mvp(m_proj * mv);

  // set transformation matrices
  OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "proj"), 1, GL_FALSE, m_proj.constData());
  OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "mv"), 1, GL_FALSE, mv.constData());
  //OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "mvp"), 1, GL_FALSE, mvp.constData());
  OGLF->glUniform2f(OGLF->glGetUniformLocation(shader_id, "screen_size"), m_w, m_h);

  // set light parameters
  OGLF->glUniform3f(OGLF->glGetUniformLocation(shader_id, "light_pos"),
                    m_light_pos.x(), m_light_pos.y(), m_light_pos.z());
  OGLF->glUniform3f(OGLF->glGetUniformLocation(shader_id, "light_col_a"),
                    m_light_ambient_col.x(), m_light_ambient_col.y(), m_light_ambient_col.z());
  OGLF->glUniform3f(OGLF->glGetUniformLocation(shader_id, "light_col_d"),
                    m_light_diffuse_col.x(), m_light_diffuse_col.y(), m_light_diffuse_col.z());

  // set attributes
  if (m_use_uniform_color)
  {
    OGLF->glUniform1i(OGLF->glGetUniformLocation(shader_id, "use_uniform_color"), 1);
    OGLF->glUniform3f(OGLF->glGetUniformLocation(shader_id, "particle_uniform_color"),
                      m_particle_col.x(), m_particle_col.y(), m_particle_col.z());
    OGLF->glDisableVertexAttribArray(5);
  }
  else
  {
    OGLF->glUniform1i(OGLF->glGetUniformLocation(shader_id, "use_uniform_color"), 0);
    OGLF->glBindBuffer(GL_ARRAY_BUFFER, vbo_colors);
    OGLF->glEnableVertexAttribArray(5);
    OGLF->glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, 0, (void *) (0));
  }

  OGLF->glBindBuffer(GL_ARRAY_BUFFER, vbo_positions);
  OGLF->glEnableVertexAttribArray(4);
  OGLF->glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 0, (void *) (0));

  // render points
  OGLF->glEnable(GL_POINT_SPRITE);
  //OGLF->glPointSize(particle_radius * 2.0f);
  OGLF->glEnable(GL_PROGRAM_POINT_SIZE);

  OGLF->glDrawArrays(GL_POINTS, 0, part_cnt);

  // clean up
  OGLF->glDisable(GL_PROGRAM_POINT_SIZE);
  OGLF->glDisable(GL_POINT_SPRITE);

  OGLF->glBindBuffer(GL_ARRAY_BUFFER, 0);

  OGLF->glUseProgram(0);

  m_vao.release();
}


bool PointSpriteRenderer::init(void)
{
  // compile shaders
  if (!utils::ogl::buildShaderProgram(m_prog,
                                      ":/src/opengl/point_sprite_renderer.vert",
                                      ":/src/opengl/point_sprite_renderer.frag"))
  {
    ERRORM("Failed to compile shaders for point sprite program");
    return false;
  }

  // set default light parameters
  m_light_pos.setX(-10.0f);
  m_light_pos.setY( 10.0f);
  m_light_pos.setZ( 15.0f);

  m_light_ambient_col.setX(0.8f);
  m_light_ambient_col.setY(0.8f);
  m_light_ambient_col.setZ(0.8f);

  m_light_diffuse_col.setX(1.0f);
  m_light_diffuse_col.setY(1.0f);
  m_light_diffuse_col.setZ(1.0f);

  // uniform color program
  GLuint uc_prog = m_prog.programId();

  OGLF->glUseProgram(uc_prog);
  //OGLF->glUniform3f(OGLF->glGetUniformLocation(uc_prog, "light_pos"), -10.0f, 10.0f, 15.0f);
  //OGLF->glUniform3f(OGLF->glGetUniformLocation(uc_prog, "light_col_a"), 0.8f, 0.8f, 0.8f);
  //OGLF->glUniform3f(OGLF->glGetUniformLocation(uc_prog, "light_col_d"), 1.0f, 1.0f, 1.0f);
  OGLF->glUniform3f(OGLF->glGetUniformLocation(uc_prog, "light_col_s"), 1.0f, 1.0f, 1.0f);

  OGLF->glUseProgram(0);

  // Create vertex array object
  return m_vao.create();
}
