#include "core/base_renderer.h"
#include "utils/ogl.h"
#include "utils/debug.h"


// forward declaration
namespace {

extern const float  g_skybox_cube_coords[];
extern const size_t g_skybox_cube_coords_size;
extern const size_t g_skybox_cube_coords_vertex_count;

} // End of private namespace



bool BaseRenderer::setSkyBoxTextures(const char *posx, const char *negx,
                                     const char *posy, const char *negy,
                                     const char *posz, const char *negz)
{
  if (!utils::ogl::loadSkyBoxTexture(m_tex_sky_box, posx, negx, posy, negy, posz, negz))
  {
    ERRORM("Failed to load Sky Box textures from: "
           << posx << "\n" << negx << "\n"
           << posy << "\n" << negy << "\n"
           << posz << "\n" << negz);
    return false;
  }
  return true;
}


bool BaseRenderer::resize(int w, int h)
{
  glViewport(0, 0, w, h);
  setPerspectiveProjection(w, h);
  return resize_impl(w, h);
}


void BaseRenderer::render(const QQuaternion & rotation,
                          const QVector3D & scale,
                          const QVector3D & translation,
                          const QQuaternion & camera_rotation,
                          GLuint vbo_positions, GLuint vbo_colors,
                          size_t part_cnt, const QVector3D & volume_size)
{
  m_stats.begin("base_renderer_total_frame_time");

  if (m_draw_sky_box)
  {
    // the skybox covers up pretty much all of the screen,
    // so we can save some instructions by clearing only the depth buffer
    OGLF->glClear(GL_DEPTH_BUFFER_BIT);
    //renderSkyBox(rotation, scale, translation);
    renderSkyBox(camera_rotation, scale, translation);
  }
  else
  {
    OGLF->glClearColor(m_clear_col.x(), m_clear_col.y(),
                       m_clear_col.z(), m_clear_col.w());
    OGLF->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  }

  if (m_draw_bounding_volume)
  {
    renderBBox(rotation, scale, translation, volume_size);
  }

  //renderUniformGrid(m_mixing ? QQuaternion() : rotation,
  //                  scale, translation,
  //                  32, 32, 32);

  render_impl(m_mixing ? QQuaternion() : rotation,
              scale, translation, camera_rotation,
              vbo_positions, vbo_colors, part_cnt,
              volume_size);

  m_stats.end("base_renderer_total_frame_time");
}


void BaseRenderer::renderUniformGrid(const QQuaternion & rotation,
                                     const QVector3D & scale,
                                     const QVector3D & translation,
                                     const size_t grid_w,
                                     const size_t grid_h,
                                     const size_t grid_d)
{
  qDebug() << __PRETTY_FUNCTION__;

  QMatrix4x4 mv;

  mv.translate(translation);
  mv.rotate(rotation);
  mv.scale(scale);

  m_prog_uniform_grid.bind();

  m_prog_uniform_grid.setUniformValue("proj", m_proj);
  m_prog_uniform_grid.setUniformValue("mv", mv);
  m_prog_uniform_grid.setUniformValue("cell_size", 8.0f);

  m_vao_uniform_grid.bind();

  OGLF->glEnable(GL_DEPTH_TEST);
  //OGLF->glDisable(GL_DEPTH_TEST);
  //OGLF->glDisable(GL_CULL_FACE);

  // x-direction
  // grid_d * grid_w

  // y-direction
  // grid_w * grid_d
  m_prog_uniform_grid.setUniformValue("len", (float) grid_h);
  m_prog_uniform_grid.setUniformValue("w", (int) grid_w);
  OGLF->glDrawArrays(GL_LINES, 0, grid_w * grid_d * 2);
  //OGLF->glDrawArrays(GL_TRIANGLES, 0, 5);

  // z-direction
  // grid_w * grid_h

  OGLF->glDisable(GL_DEPTH_TEST);

  m_vao_uniform_grid.release();

  m_prog_uniform_grid.release();
}


void BaseRenderer::renderBBox(const QQuaternion & rotation,
                              const QVector3D & scale,
                              const QVector3D & translation,
                              const QVector3D & size)
{
  QMatrix4x4 mv;

  mv.translate(translation);
  //mv.translate(0.0f, 0.0f, -1.0f);
  mv.rotate(rotation);
  mv.scale(scale);
  //mv.scale(1.1f, 1.1f, 1.1f);   // to make the bounding box slightly larger than its contents

  m_prog_bounding_volume.bind();

  m_prog_bounding_volume.setUniformValue("proj", m_proj);
  m_prog_bounding_volume.setUniformValue("mv", mv);
  //m_prog_bounding_volume.setUniformValue("dimensions", m_volume_size);
  m_prog_bounding_volume.setUniformValue("dimensions", size);

  m_vao_bounding_volume.bind();
  OGLF->glEnable(GL_DEPTH_TEST);
  OGLF->glDrawArrays(GL_LINES, 0, 24);
  OGLF->glDisable(GL_DEPTH_TEST);
  m_vao_bounding_volume.release();

  m_prog_bounding_volume.release();
}


void BaseRenderer::renderSkyBox(const QQuaternion & rotation,
                                const QVector3D & /* scale */,
                                const QVector3D & translation)
{
#if 0
  QMatrix4x4 mv;
  //mv.translate(0.0f, 0.0f, -32.0f);
  //mv.translate(translation);
  //mv.translate(0.0f, 0.0f, -40.0f);
  //mv.translate(0.0f, 0.0f, -80.0f);
  mv.rotate(rotation);
  mv.scale(scale);
#elif 0
  QMatrix4x4 mv;
  mv.translate(translation);
  mv.scale(scale);
  mv.rotate(rotation);
#else
  QMatrix4x4 mv;
  //mv.translate(translation);
  //mv.translate(0.0f, 0.0f, translation.z());
  //float sz = 100.0f;
  //mv.translate(0.0f, 0.0f, -sz);
  mv.rotate(rotation);
#endif

  QMatrix4x4 mvp(m_proj * mv);

  OGLF->glActiveTexture(GL_TEXTURE0);
  OGLF->glBindTexture(GL_TEXTURE_CUBE_MAP, m_tex_sky_box.textureId());

  m_prog_sky_box.bind();
  m_prog_sky_box.setUniformValue("mvp", mvp);
  m_prog_sky_box.setUniformValue("sky_box_cube_size", SKY_BOX_SIZE_HALF); //100.0f); //40.0f); // 10.0f
  //m_prog_sky_box.setUniformValue("sky_box_cube_size", sz);
  m_prog_sky_box.setUniformValue("tex_sky_box_cubemap", 0);

  m_vao_sky_box_cube.bind();
  OGLF->glDepthMask(GL_FALSE);
  OGLF->glEnable(GL_CULL_FACE);
  OGLF->glDrawArrays(GL_TRIANGLES, 0, g_skybox_cube_coords_vertex_count);
  OGLF->glDisable(GL_CULL_FACE);
  OGLF->glDepthMask(GL_TRUE);
  m_vao_sky_box_cube.release();

  m_prog_sky_box.release();
}


bool BaseRenderer::init(void)
{
  OGLF->glGetIntegerv(GL_FRAMEBUFFER_BINDING, (GLint *) &m_default_fbo);
  DBGM("default frame buffer object id: " << m_default_fbo);
  return initBoundingBox() && initUniformGrid() && initSkyBox();
}


bool BaseRenderer::initBoundingBox(void)
{
  if (!utils::ogl::buildShaderProgram(m_prog_bounding_volume,
                                      ":/src/opengl/ps_bounding_volume.vert",
                                      ":/src/opengl/ps_bounding_volume.frag"))
  {
    ERRORM("Failed to compile shaders for bounding volume program");
    return false;
  }

  m_prog_bounding_volume.bind();
  m_prog_bounding_volume.setUniformValue("col", QVector3D(1.0f, 0.0f, 0.0f));
  m_prog_bounding_volume.release();

  if (!m_vao_bounding_volume.create())
  {
    ERRORM("Failed to create Vertex Array Object for bounding volume program in Base Renderer");
    return false;
  }

  return true;
}


bool BaseRenderer::initUniformGrid(void)
{
  if (!utils::ogl::buildShaderProgram(m_prog_uniform_grid,
                                      ":/src/opengl/uniform_grid.vert",
                                      ":/src/opengl/uniform_grid.frag"))
  {
    ERRORM("Failed to compile shaders for uniform grid program");
    return false;
  }

  m_prog_uniform_grid.bind();
  //m_prog_uniform_grid.setUniformValue("col", QVector3D(1.0f, 0.0f, 0.0f));
  m_prog_uniform_grid.setUniformValue("col", QVector3D(0.0f, 1.0f, 0.0f));
  m_prog_uniform_grid.release();

  if (!m_vao_uniform_grid.create())
  {
    ERRORM("Failed to create Vertex Array Object uniform grid program in Base Renderer");
    return false;
  }

  return true;
}


bool BaseRenderer::initSkyBox(void)
{
  //setSkyBoxTextures(":/data/skybox/debug_posx.png",
  //                  ":/data/skybox/debug_negx.png",
  //                  ":/data/skybox/debug_posy.png",
  //                  ":/data/skybox/debug_negy.png",
  //                  ":/data/skybox/debug_posz.png",
  //                  ":/data/skybox/debug_negz.png");

  setSkyBoxTextures(":/data/skybox/debug_large_posx.png",
                    ":/data/skybox/debug_large_negx.png",
                    ":/data/skybox/debug_large_posy.png",
                    ":/data/skybox/debug_large_negy.png",
                    ":/data/skybox/debug_large_posz.png",
                    ":/data/skybox/debug_large_negz.png");

  // create buffer object for the cube geometry
  m_vbo_sky_box_cube.setUsagePattern(QOpenGLBuffer::StaticDraw);
  if ((!m_vao_sky_box_cube.create()) || (!m_vbo_sky_box_cube.create()))
  {
    ERRORM("Failed to create vertex buffer object or "
           "vertex array object for sky box cube");
    return false;
  }

  m_vao_sky_box_cube.bind();

  m_vbo_sky_box_cube.bind();
  m_vbo_sky_box_cube.allocate(g_skybox_cube_coords, g_skybox_cube_coords_size);

  // initialize the skybox shader program
  if (!utils::ogl::buildShaderProgram(m_prog_sky_box,
                                      ":/src/opengl/skybox.vert",
                                      ":/src/opengl/skybox.frag"))
  {
    ERRORM("Failed to build SkyBox shader program");
    return false;
  }

  // initialize attributes and parameters that do not change much
  // while the application is running
  m_prog_sky_box.bind();

  int attr_pos = m_prog_sky_box.attributeLocation("pos");
  OGLF->glEnableVertexAttribArray(attr_pos);
  OGLF->glVertexAttribPointer(attr_pos, 3, GL_FLOAT, GL_FALSE, 0, (void *) 0);

  m_prog_sky_box.release();

  m_vao_sky_box_cube.release();

  return true;
}


namespace {

const float g_skybox_cube_coords[] = {
  -1.0f,  1.0f, -1.0f,
  -1.0f, -1.0f, -1.0f,
   1.0f, -1.0f, -1.0f,
   1.0f, -1.0f, -1.0f,
   1.0f,  1.0f, -1.0f,
  -1.0f,  1.0f, -1.0f,

  -1.0f, -1.0f,  1.0f,
  -1.0f, -1.0f, -1.0f,
  -1.0f,  1.0f, -1.0f,
  -1.0f,  1.0f, -1.0f,
  -1.0f,  1.0f,  1.0f,
  -1.0f, -1.0f,  1.0f,

   1.0f, -1.0f, -1.0f,
   1.0f, -1.0f,  1.0f,
   1.0f,  1.0f,  1.0f,
   1.0f,  1.0f,  1.0f,
   1.0f,  1.0f, -1.0f,
   1.0f, -1.0f, -1.0f,

  -1.0f, -1.0f,  1.0f,
  -1.0f,  1.0f,  1.0f,
   1.0f,  1.0f,  1.0f,
   1.0f,  1.0f,  1.0f,
   1.0f, -1.0f,  1.0f,
  -1.0f, -1.0f,  1.0f,

  -1.0f,  1.0f, -1.0f,
   1.0f,  1.0f, -1.0f,
   1.0f,  1.0f,  1.0f,
   1.0f,  1.0f,  1.0f,
  -1.0f,  1.0f,  1.0f,
  -1.0f,  1.0f, -1.0f,

  -1.0f, -1.0f, -1.0f,
  -1.0f, -1.0f,  1.0f,
   1.0f, -1.0f, -1.0f,
   1.0f, -1.0f, -1.0f,
  -1.0f, -1.0f,  1.0f,
   1.0f, -1.0f,  1.0f
};

const size_t g_skybox_cube_coords_size = sizeof(g_skybox_cube_coords);
const size_t g_skybox_cube_coords_vertex_count = sizeof(g_skybox_cube_coords) /
                                                 (sizeof(g_skybox_cube_coords[0]) * 3);

} // End of private namespace
