#include "core/instancing_renderer.h"
#include "utils/ogl.h"
#include "utils/debug.h"

// include OpenCL definitions in crossplatform manner
//#include <boost/compute/cl.hpp>



bool InstancingRenderer::reset_impl(int /* w */, int /* h */)
{
  /* compile shaders */
  if (!utils::ogl::buildShaderProgram(m_prog_particle_colors,
                                      ":/src/opengl/ps_particle_colors.vert",
                                      ":/src/opengl/ps_particle_colors.frag"))
  {
    ERRORM("Failed to compile shaders for particle colors program");
    return false;
  }

  if (!utils::ogl::buildShaderProgram(m_prog_uniform_color,
                                      ":/src/opengl/ps_uniform_color.vert",
                                      ":/src/opengl/ps_uniform_color.frag"))
  {
    ERRORM("Failed to compile shaders for uniform color program");
    return false;
  }

  /* per particle colors program */
  GLuint pc_prog = m_prog_particle_colors.programId();

  OGLF->glUseProgram(pc_prog);
  //OGLF->glUniform3f(glGetUniformLocation(prog, "light_pos"), 0.0f, 0.0f, -1.0f);     // shader assumes the position is normalized
  //OGLF->glUniform3f(glGetUniformLocation(prog, "light_pos"), -1.0f, 1.0f, -1.0f);
  OGLF->glUniform3f(OGLF->glGetUniformLocation(pc_prog, "light_pos"), -10.0f, 10.0f, 15.0f);
  OGLF->glUniform3f(OGLF->glGetUniformLocation(pc_prog, "light_col_a"), 0.8f, 0.8f, 0.8f);
  OGLF->glUniform3f(OGLF->glGetUniformLocation(pc_prog, "light_col_d"), 1.0f, 1.0f, 1.0f);
  OGLF->glUniform3f(OGLF->glGetUniformLocation(pc_prog, "light_col_s"), 1.0f, 1.0f, 1.0f);

  /* uniform color program */
  GLuint uc_prog = m_prog_uniform_color.programId();

  OGLF->glUseProgram(uc_prog);
  OGLF->glUniform3f(OGLF->glGetUniformLocation(uc_prog, "light_pos"), -10.0f, 10.0f, 15.0f);
  OGLF->glUniform3f(OGLF->glGetUniformLocation(uc_prog, "light_col_a"), 0.8f, 0.8f, 0.8f);
  OGLF->glUniform3f(OGLF->glGetUniformLocation(uc_prog, "light_col_d"), 1.0f, 1.0f, 1.0f);
  OGLF->glUniform3f(OGLF->glGetUniformLocation(uc_prog, "light_col_s"), 1.0f, 1.0f, 1.0f);
  //OGLF->glUniform3f(OGLF->glGetUniformLocation(uc_prog, "particle_col"), 0.5f, 0.5f, 1.0f);

  OGLF->glUseProgram(0);

  /* load models */
  if (!utils::geom::genSphere(m_particle_geom))
  //if (!geom::genPrism(m_particle_geom))
  {
    ERRORM("Failed to generate sphere model");
    return false;
  }

  return true;
}


void InstancingRenderer::render_impl(const QQuaternion & rotation,
                                     const QVector3D & scale,
                                     const QVector3D & translation,
                                     const QQuaternion & /* camera_rotation */,
                                     GLuint vbo_positions,
                                     GLuint vbo_colors,
                                     size_t part_cnt,
                                     const QVector3D & /* volume_size */)
{
  // calculate scene transformation
  QMatrix4x4 mv;

  mv.translate(translation);
  mv.rotate(rotation);
  mv.scale(scale);

  // render the particle system
  OGLF->glEnable(GL_DEPTH_TEST);

  OGLF->glBindVertexArray(m_particle_geom.vao);

  GLuint shader_id = 0;

  if (m_use_uniform_color)
  {
    shader_id = m_prog_uniform_color.programId();
    OGLF->glUseProgram(shader_id);
    OGLF->glUniform3f(OGLF->glGetUniformLocation(shader_id, "particle_col"),
                      m_particle_col.x(), m_particle_col.y(), m_particle_col.z());
  }
  else
  {
    shader_id = m_prog_particle_colors.programId();
    OGLF->glUseProgram(shader_id);

    OGLF->glBindBuffer(GL_ARRAY_BUFFER, vbo_colors);
    OGLF->glEnableVertexAttribArray(5);
    // for performance reasons the position is 4 component array
    OGLF->glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, 0, (void *) (0));
    //OGLF->glVertexAttribPointer(5, 4, GL_FLOAT, GL_FALSE, sizeof(cl_float4), (void *) (0));
    OGLF->glVertexAttribDivisor(5, 1);
  }

  OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "proj"), 1, GL_FALSE, m_proj.constData());
  OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "mv"), 1, GL_FALSE, mv.constData());

  /* calculate the normal matrix (assume only rotations, i.e only orthogonal matrices) */
  QMatrix3x3 mv_normal = mv.toGenericMatrix<3, 3>();   // otherwise add glm::transpose(glm::inverse(mv));
  OGLF->glUniformMatrix3fv(OGLF->glGetUniformLocation(shader_id, "mv_normal"), 1, GL_FALSE, mv_normal.constData());

  OGLF->glBindBuffer(GL_ARRAY_BUFFER, vbo_positions);
  OGLF->glEnableVertexAttribArray(4);
  // for performance reasons the position is 4 component array
  OGLF->glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 0, (void *) (0));
  //OGLF->glVertexAttribPointer(4, 3, GL_FLOAT, GL_FALSE, sizeof(cl_float4), (void *) (0));
  OGLF->glVertexAttribDivisor(4, 1);

  OGLF->glDrawArraysInstanced(m_particle_geom.mode, 0, m_particle_geom.count, part_cnt);

  OGLF->glBindBuffer(GL_ARRAY_BUFFER, 0);

  OGLF->glBindVertexArray(0);

  OGLF->glUseProgram(0);
}
