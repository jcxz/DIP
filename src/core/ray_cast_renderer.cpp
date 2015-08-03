#include "core/ray_cast_renderer.h"
#include "utils/ogl.h"
#include "utils/ocl.h"
#include "utils/debug.h"

#include <QMatrix4x4>
#include <QFile>


#define F() qDebug() << __PRETTY_FUNCTION__;


namespace {

const GLfloat g_cube_vertices[][3] = {
  { -1.0f, -1.0f, -1.0f },  // 0
  { -1.0f, -1.0f,  1.0f },  // 1
  { -1.0f,  1.0f, -1.0f },  // 2
  { -1.0f,  1.0f,  1.0f },  // 3
  {  1.0f, -1.0f, -1.0f },  // 4
  {  1.0f, -1.0f,  1.0f },  // 5
  {  1.0f,  1.0f, -1.0f },  // 6
  {  1.0f,  1.0f,  1.0f }   // 7
};

const int g_cube_vertices_cnt = sizeof(g_cube_vertices) / sizeof(g_cube_vertices[0]);

const GLuint g_cube_indices[] = {
  // predna stena
  1, 5, 7,
  1, 7, 3,

  // zadna stena
  2, 4, 0,
  2, 6, 4,

  // prava stena
  5, 4, 6,
  5, 6, 7,

  // lava stena
  0, 1, 3,
  0, 3, 2,

  // horna stena
  3, 7, 6,
  3, 6, 2,

  // dolna stena
  0, 4, 5,
  0, 5, 1
};

const int g_cube_indices_cnt = sizeof(g_cube_indices) / sizeof(g_cube_indices[0]);

} // End of private namespace


RayCastRenderer::~RayCastRenderer(void)
{
  F();

  OGLF->glDeleteTextures(1, &m_color_attach);
  // pretoze Qt pouziva fbo cislo 1
  if ((m_fbo != 0) && (m_fbo != m_default_fbo))
  {
    OGLF->glDeleteFramebuffers(1, &m_fbo);
  }
}


bool RayCastRenderer::resize_impl(int w, int h)
{
  F();
  return resetFramebuffer(w, h);
}


bool RayCastRenderer::reset_impl(int w, int h)
{
  F();

  if ((m_data_width <= 0) || (m_data_height <= 0) || (m_data_depth <= 0))
  {
    ERRORM("RayCastRenderer: grid size is not set");
    return false;
  }

  if ((m_cell_starts_buf.get() == nullptr) || (m_cell_ends_buf.get() == nullptr))
  {
    ERRORM("RayCastRenderer: cell_starts or cell_ends buffer is not set");
    return false;
  }

  return resetDataTexture() &&
      resetBuffers() &&
      resetFramebuffer(w, h) &&
      resetTransferFunctionTexture();
}


bool RayCastRenderer::initGL(void)
{
  F();

  // kompilacia shaderov
  m_prog_gen_back_faces.addShaderFromSourceFile(QOpenGLShader::Vertex, ":/src/opengl/ray_cast_renderer.vert");
  m_prog_gen_back_faces.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/src/opengl/ray_cast_renderer_back_faces.frag");
  m_prog_gen_back_faces.link();

  //*** Kompilacia shader programu, ktory zobrazuje frame bufffer na obrazovku
  m_prog_ray_cast.addShaderFromSourceFile(QOpenGLShader::Vertex,   ":/src/opengl/ray_cast_renderer.vert");
  m_prog_ray_cast.addShaderFromSourceFile(QOpenGLShader::Fragment, ":/src/opengl/ray_cast_renderer.frag");
  m_prog_ray_cast.link();

  //*** Vytvorenie OpenGL bufferov
  m_vao.create();
  m_vao.bind();

  // tieto 2 prikazy v skutocnosti ani nie su potrebne,
  // pretoze StaticDraw je pre QOpenGLBuffer defaultny
  m_vbo_cube.setUsagePattern(QOpenGLBuffer::StaticDraw);
  m_ibo_cube.setUsagePattern(QOpenGLBuffer::StaticDraw);

  m_vbo_cube.create();
  m_ibo_cube.create();

  m_ibo_cube.bind();
  m_ibo_cube.allocate(g_cube_indices, sizeof(g_cube_indices));

  m_vbo_cube.bind();
  m_vbo_cube.allocate(g_cube_vertices, sizeof(g_cube_vertices));

  // nastavenie arributov
  m_prog_gen_back_faces.bind();
  int attr_pos = m_prog_gen_back_faces.attributeLocation("pos");

  OGLF->glEnableVertexAttribArray(attr_pos);
  OGLF->glVertexAttribPointer(attr_pos, 3, GL_FLOAT, GL_FALSE, sizeof(g_cube_vertices[0]), (void *) 0);

  // odbindovanie VAO a buffer objectu
  OGLF->glBindVertexArray(0);
  OGLF->glBindBuffer(GL_ARRAY_BUFFER, 0);
  OGLF->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);

  OGLF->glUseProgram(0);

  return true;
}


bool RayCastRenderer::initCL(const boost::compute::context & ctx,
                             const boost::compute::command_queue & queue)
{
  F();

  try
  {
    // build opencl program
    QFile f(":/src/opencl/ray_cast_renderer.cl");
    if (!f.open(QFile::ReadOnly)) return false;

    m_prog = boost::compute::program::create_with_source(f.readAll().toStdString(), ctx);
    m_prog.build();

    // print build log in case there are any warnings
    std::string log(m_prog.build_log());
    if (!log.empty())
    {
      WARNM("---------------------- Build log ---------------------------------");
      WARNM(log);
      WARNM("------------------- End of Build log -----------------------------");
    }

    // create kernels
    m_calculate_particle_cnts_kernel = m_prog.create_kernel("ray_cast_renderer_calculate_particle_cnts");
    utils::ocl::printKernelInfo(m_calculate_particle_cnts_kernel);

    m_calculate_gradients_kernel = m_prog.create_kernel("ray_cast_renderer_calculate_gradients");
    utils::ocl::printKernelInfo(m_calculate_gradients_kernel);

    m_calculate_normals_kernel = m_prog.create_kernel("ray_cast_renderer_calculate_normals");
    utils::ocl::printKernelInfo(m_calculate_normals_kernel);

    m_queue = queue;
    m_cl_ctx = ctx;
  }
  catch (const boost::compute::opencl_error & e)
  {
    ERRORM(e.what());
    if (e.error_code() == CL_BUILD_PROGRAM_FAILURE) ERRORM(m_prog.build_log());
    return false;
  }
  catch (const std::exception & e)
  {
    ERRORM("An unexpected error occured during SPHUniformGrid initialization: " << e.what());
    return false;
  }

  return true;
}


bool RayCastRenderer::resetFramebuffer(int w, int h)
{
  F();

  // alokacia textury
  if (m_color_attach == 0) OGLF->glGenTextures(1, &m_color_attach);

  OGLF->glBindTexture(GL_TEXTURE_2D, m_color_attach);

  OGLF->glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_UNSIGNED_BYTE, nullptr);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

  OGLF->glBindTexture(GL_TEXTURE_2D, 0);

  // Nastavenie framebufferu
  if (m_fbo == 0) OGLF->glGenFramebuffers(1, &m_fbo);

  OGLF->glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);

  OGLF->glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_color_attach, 0);

  // Nastavenie zoznamu attachementov, do ktorych sa bude kreslit
  GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT0 };
  OGLF->glDrawBuffers(1, draw_buffers);

  if (OGLF->glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
  {
    WARNM("**initFramebuffer: Failed to create frame buffer: frame buffer is not complete");
    return false;
  }

  OGLF->glBindFramebuffer(GL_FRAMEBUFFER, m_default_fbo);

  return true;
}


bool RayCastRenderer::resetDataTexture(void)
{
  F();

  // nahratie dat do textury
  if (!m_data.create())
  {
    ERRORM("Failed to create OpenGL 3D texture for volume data");
    return false;
  }

  OGLF->glBindTexture(GL_TEXTURE_3D, m_data.textureId());

  //OGLF->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
  //OGLF->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);
  //OGLF->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_BORDER);
  OGLF->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  OGLF->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  OGLF->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  OGLF->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  OGLF->glTexParameteri(GL_TEXTURE_3D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

  OGLF->glTexImage3D(GL_TEXTURE_3D, 0, GL_RGBA,
                     m_data_width, m_data_height, m_data_depth,
                     0, GL_RGBA, GL_FLOAT, nullptr);

  OGLF->glBindTexture(GL_TEXTURE_3D, 0);

  return true;
}


bool RayCastRenderer::resetTransferFunctionTexture(void)
{
  F();

  // nahratie dat do textury
  if (!m_transfer_func.create())
  {
    ERRORM("Failed to create OpenGL 1D texture for transfer function");
    return false;
  }

  static constexpr int n = 1024;
  static constexpr float factor = 255.0f / float(n);

  unsigned char *pixels = new unsigned char[n * 4];

#define BLUE_FOAM_CLEAR_WATER

  for (int i = 0; i < n; ++i)
  {
#ifdef BLUE_FOAM_CLEAR_WATER
    //if (i < 5)
    //if (i < -1)
    if (i > 256)
    {
      pixels[i * 4 + 0] = 0;  // red
      pixels[i * 4 + 1] = 0;  // green
      pixels[i * 4 + 2] = i * factor;  // blue
      pixels[i * 4 + 3] = 0;
    }
    else
    {
      pixels[i * 4 + 0] = 0;  // red
      pixels[i * 4 + 1] = 0;  // green
      pixels[i * 4 + 2] = i * factor;  // blue
      //pixels[i * 4 + 3] = i * factor; //255;
      pixels[i * 4 + 3] = i * factor; //255;
    }
#elif YELLOW_GREEN_WATER
    if (i > 256)
    {
      pixels[i * 4 + 0] = 0;  // red
      pixels[i * 4 + 1] = 0;  // green
      pixels[i * 4 + 2] = i * factor;  // blue
      pixels[i * 4 + 3] = 0;
    }
    else
    {
      pixels[i * 4 + 0] = i * i * factor;  // red
      pixels[i * 4 + 1] = i * i * factor;  // green
      pixels[i * 4 + 2] = i * factor;  // blue
      pixels[i * 4 + 3] = i * factor; //255;
    }
#elif GRAY_WATER
    if (i > 256)
    {
      pixels[i * 4 + 0] = 0;  // red
      pixels[i * 4 + 1] = 0;  // green
      pixels[i * 4 + 2] = i * factor;  // blue
      pixels[i * 4 + 3] = 0;
    }
    else
    {
      pixels[i * 4 + 0] = i * i * factor;  // red
      pixels[i * 4 + 1] = i * i * factor;  // green
      pixels[i * 4 + 2] = i * i * factor;  // blue
      pixels[i * 4 + 3] = i * factor; //255;
    }
#else
    if (i > 256)
    {
      pixels[i * 4 + 0] = 0;  // red
      pixels[i * 4 + 1] = 0;  // green
      pixels[i * 4 + 2] = i * factor;  // blue
      pixels[i * 4 + 3] = 10 / i;
    }
    else
    {
      pixels[i * 4 + 0] = i * factor;  // red
      pixels[i * 4 + 1] = i * factor;  // green
      pixels[i * 4 + 2] = (256 / i) * factor;  // blue
      pixels[i * 4 + 3] = i * factor; //255;
    }
#endif
  }

  m_transfer_func.bind();

  OGLF->glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);
  OGLF->glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  OGLF->glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);

  OGLF->glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, n, 0, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

  OGLF->glBindTexture(GL_TEXTURE_1D, 0);

  delete [] pixels;

  return true;
}


// This function has to be called after all the OpenGL textures have been created
bool RayCastRenderer::resetBuffers(void)
{
  F();

  size_t cell_cnt = m_data_width * m_data_height * m_data_depth;

  try
  {
    m_cell_cnts_buf = boost::compute::buffer(m_cl_ctx, cell_cnt * sizeof(cl_uint));
    m_grads_x_buf = boost::compute::buffer(m_cl_ctx, cell_cnt * sizeof(cl_float));
    m_grads_y_buf = boost::compute::buffer(m_cl_ctx, cell_cnt * sizeof(cl_float));
    m_grads_z_buf = boost::compute::buffer(m_cl_ctx, cell_cnt * sizeof(cl_float));
    m_data_cl_img = boost::compute::opengl_texture(m_cl_ctx, GL_TEXTURE_3D, 0, m_data.textureId());
  }
  catch (const boost::compute::opencl_error & e)
  {
    ERRORM(e.what());
    return false;
  }
  catch (const std::exception & e)
  {
    ERRORM("An unexpected error occured during SPHUniformGrid initialization: " << e.what());
    return false;
  }

  return true;
}


void RayCastRenderer::calcNormals(void)
{
  //F();

  // synchronize with OpenGL
  utils::ocl::GLSyncHandler sync(m_queue, 1, &m_data_cl_img.get());
  if (!sync) return;

  size_t global_work_size[3] = { m_data_width, m_data_height, m_data_depth };
  size_t cell_cnt = m_data_width * m_data_height * m_data_depth;
  cl_int err = CL_SUCCESS;

  //*** calculate particle counts
  utils::ocl::KernelArgs(m_calculate_particle_cnts_kernel,
                         "ray_cast_renderer_calculate_particle_cnts")
      .arg(m_cell_starts_buf)
      .arg(m_cell_ends_buf)
      .arg(m_cell_cnts_buf);
  err = clEnqueueNDRangeKernel(m_queue, m_calculate_particle_cnts_kernel, 1,
                               nullptr, &cell_cnt, nullptr,
                               0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    WARNM("RayCastRenderer: Failed to enqueue ray_cast_renderer_calculate_particle_cnts kernel"
          << boost::compute::opencl_error::to_string(err)
          << "(" << err << ")");
  }

  //*** calculate gradients
  utils::ocl::KernelArgs(m_calculate_gradients_kernel,
                         "ray_cast_renderer_calculate_gradients")
      .arg(m_cell_cnts_buf)
      .arg(m_grads_x_buf)
      .arg(m_grads_y_buf)
      .arg(m_grads_z_buf)
      .arg<cl_uint>(m_data_width)
      .arg<cl_uint>(m_data_height)
      .arg<cl_uint>(m_data_depth);
  err = clEnqueueNDRangeKernel(m_queue, m_calculate_gradients_kernel, 3,
                               nullptr, global_work_size, nullptr,
                               0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    WARNM("RayCastRenderer: Failed to enqueue ray_cast_renderer_calculate_gradients kernel"
          << boost::compute::opencl_error::to_string(err)
          << "(" << err << ")");
  }

  //*** calculate normals
  utils::ocl::KernelArgs(m_calculate_normals_kernel,
                         "ray_cast_renderer_calculate_normals")
      .arg(m_cell_cnts_buf)
      .arg(m_grads_x_buf)
      .arg(m_grads_y_buf)
      .arg(m_grads_z_buf)
      .arg(m_data_cl_img)
      .arg<cl_uint>(m_data_width)
      .arg<cl_uint>(m_data_height)
      .arg<cl_uint>(m_data_depth);
  err = clEnqueueNDRangeKernel(m_queue, m_calculate_normals_kernel, 3,
                               nullptr, global_work_size, nullptr,
                               0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    WARNM("RayCastRenderer: Failed to enqueue ray_cast_renderer_calculate_normals kernel"
          << boost::compute::opencl_error::to_string(err)
          << "(" << err << ")");
  }

#if 0
  cl_float4 fill_color = { 1.0f, 1.0f, 1.0f, 1.0f };
  size_t origin[3] = { 0, 0, 0 };
  size_t region[3] = { m_data_width, m_data_height, m_data_depth };
  err = clEnqueueFillImage(m_queue, m_data_cl_img.get(),
                           &fill_color, origin, region,
                           0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    WARNM("RayCastRenderer: Failed to fill image"
          << boost::compute::opencl_error::to_string(err)
          << "(" << err << ")");
  }
#endif
}


void RayCastRenderer::render_impl(const QQuaternion & rotation,
                                  const QVector3D & scale,
                                  const QVector3D & translation,
                                  const QQuaternion & /* camera_rotation */,
                                  GLuint /* vbo_positions */,
                                  GLuint /* vbo_colors */,
                                  size_t /* part_cnt */,
                                  const QVector3D & volume_size)
                                  //float peel_depth,
                                  //int detail
{
  if (!isDataValid())
  {
    WARNM("Not rendering because: Failed to create 3D texture for data");
    return;
  }

  // Transformacna model-view matica
  QMatrix4x4 mv;

  mv.translate(0.0f, 0.0f, -1.0f);
  mv.translate(translation);
  mv.rotate(rotation);
  mv.scale(scale);

  // Uprava rozmerov volumetrickych dat, tak aby presne sedeli na jednotkovu kocku
  float max_size = std::max(volume_size.x(), std::max(volume_size.y(), volume_size.z()));
  mv.scale((volume_size.x() / max_size),
           (volume_size.y() / max_size),
           (volume_size.z() / max_size));

  mv.scale(volume_size.x(), volume_size.y(), volume_size.z());
  //mv.scale(volume_size.x() / 2.0f, volume_size.y() / 2.0f, volume_size.z() / 2.0f);
  //mv.scale(volume_size.x() / 4.0f, volume_size.y() / 4.0f, volume_size.z() / 4.0f);
  //mv.scale(volume_size.x() / 8.0f, volume_size.y() / 8.0f, volume_size.z() / 8.0f);

  // nabindovanie textur, geometrie, povolenie cullingu, depth testovania a HW blendingu (robi sa v shaderi)
  OGLF->glActiveTexture(GL_TEXTURE0);
  OGLF->glBindTexture(GL_TEXTURE_1D, m_transfer_func.textureId());
  OGLF->glActiveTexture(GL_TEXTURE1);
  OGLF->glBindTexture(GL_TEXTURE_2D, m_color_attach);
  OGLF->glActiveTexture(GL_TEXTURE2);
  OGLF->glBindTexture(GL_TEXTURE_3D, m_data.textureId());

  m_vao.bind();

  OGLF->glEnable(GL_DEPTH_TEST);
  OGLF->glDepthMask(GL_FALSE);
  OGLF->glEnable(GL_CULL_FACE);
  // blending je potreba povolit, aby neboli casti ciar z bounding box-u prekryte
  // farbou vzduchu vo volumetrickych datach
  OGLF->glEnable(GL_BLEND);
  //OGLF->glDisable(GL_BLEND);

  // Kreslenie sceny do framebufferu (vygenerovanie exit pointov)s
  OGLF->glBindFramebuffer(GL_FRAMEBUFFER, m_fbo);
  OGLF->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  OGLF->glCullFace(GL_FRONT);

  m_prog_gen_back_faces.bind();
  m_prog_gen_back_faces.setUniformValue("mvp", m_proj * mv);

  OGLF->glDrawElements(GL_TRIANGLES, g_cube_indices_cnt, GL_UNSIGNED_INT, nullptr);

  // Ray-casting
  OGLF->glBindFramebuffer(GL_FRAMEBUFFER, m_default_fbo);
  // Toto uz netreba, pretoze defaultny frame buffer clearuje base class
  // a pripadne este aj kresli bounding box ak treba
  //OGLF->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  OGLF->glCullFace(GL_BACK);

  int detail = 1000; //5000;
  float default_step = 1.0f / max_size;
  float step;

  if (detail <= 0)
    step = default_step;
  else
    step = 1.0f / float(detail);

  //qDebug() << "step=" << step << ", default_step=" << default_step << ", correction factor=" << (step / default_step);

  m_prog_ray_cast.bind();
  m_prog_ray_cast.setUniformValue("mvp", m_proj * mv);
  m_prog_ray_cast.setUniformValue("step", step);
  //m_prog_ray_cast.setUniformValue("offset", peel_depth);
  m_prog_ray_cast.setUniformValue("offset", 0.0f);
  m_prog_ray_cast.setUniformValue("alpha_correction_factor", step / default_step);
  m_prog_ray_cast.setUniformValue("tex_transfer_func", 0);
  m_prog_ray_cast.setUniformValue("tex_back_faces", 1);
  m_prog_ray_cast.setUniformValue("tex_volume_data", 2);
  //m_prog_ray_cast.setUniformValue("use_tf", m_use_transfer_func);

  if (m_use_lighting && m_use_transfer_func)
  {
    m_prog_ray_cast.setUniformValue("method", 1);
    m_prog_ray_cast.setUniformValue("La", m_light_ambient_col);
    m_prog_ray_cast.setUniformValue("Ld", m_light_diffuse_col);
    m_prog_ray_cast.setUniformValue("light_pos", m_light_pos);
  }
  else if (m_use_transfer_func)
  {
    m_prog_ray_cast.setUniformValue("method", 2);
  }
  else
  {
    m_prog_ray_cast.setUniformValue("method", 3);
  }

  OGLF->glEnable(GL_DEPTH_TEST);
  OGLF->glDrawElements(GL_TRIANGLES, g_cube_indices_cnt, GL_UNSIGNED_INT, nullptr);

  // Odbindovanie programu, geometrie, textur a zakazanie cullingu
  OGLF->glUseProgram(0);

  OGLF->glDisable(GL_BLEND);
  OGLF->glDisable(GL_CULL_FACE);
  OGLF->glDepthMask(GL_TRUE);
  OGLF->glDisable(GL_DEPTH_TEST);

  OGLF->glBindVertexArray(0);

  OGLF->glActiveTexture(GL_TEXTURE2);
  OGLF->glBindTexture(GL_TEXTURE_3D, 0);
  OGLF->glActiveTexture(GL_TEXTURE1);
  OGLF->glBindTexture(GL_TEXTURE_2D, 0);
  OGLF->glActiveTexture(GL_TEXTURE0);
  OGLF->glBindTexture(GL_TEXTURE_1D, 0);
}
