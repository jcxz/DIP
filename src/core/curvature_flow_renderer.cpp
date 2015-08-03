#include "core/curvature_flow_renderer.h"
#include "utils/ogl.h"
#include "utils/debug.h"



CurvatureFlowRenderer::~CurvatureFlowRenderer(void)
{
  // pretoze Qt pouziva fbo cislo 1
  if ((m_fbo_depth != 0) && (m_fbo_depth != m_default_fbo))
  {
    OGLF->glDeleteFramebuffers(1, &m_fbo_depth);
  }

  // delete textures and renderbuffers attached to depth FBO
  if (m_tex_depth[0] == 0) OGLF->glGenTextures(2, m_tex_depth);
  if (m_rbo_depth == 0) OGLF->glGenRenderbuffers(1, &m_rbo_depth);

  // pretoze Qt pouziva fbo cislo 1
  if ((m_fbo_thickness != 0) && (m_fbo_thickness != m_default_fbo))
  {
    OGLF->glDeleteFramebuffers(1, &m_fbo_thickness);
  }
}


bool CurvatureFlowRenderer::resize_impl(int w, int h)
{
  m_w = w;
  m_h = h;
  return resetDepthFBO(w, h) && resetThicknessFBO(w, h);
}


bool CurvatureFlowRenderer::reset_impl(int w, int h)
{
  m_w = w;
  m_h = h;
  return resetDepthFBO(w, h) && resetThicknessFBO(w, h);
}


void CurvatureFlowRenderer::render_impl(const QQuaternion & rotation,
                                        const QVector3D & scale,
                                        const QVector3D & translation,
                                        const QQuaternion & camera_rotation,
                                        GLuint vbo_positions,
                                        GLuint /* vbo_colors */,
                                        size_t part_cnt,
                                        const QVector3D & /* volume_size */)
{
  m_stats.begin("curvature_flow_renderer_total_time");

  // calculate scene transformation
  QMatrix4x4 mv;

  mv.translate(translation);
  mv.rotate(rotation);
  mv.scale(scale);

  m_vao.bind();

  // set up vertex buffer with particle positions
  OGLF->glBindBuffer(GL_ARRAY_BUFFER, vbo_positions);
  OGLF->glEnableVertexAttribArray(4);
  OGLF->glVertexAttribPointer(4, 4, GL_FLOAT, GL_FALSE, 0, (void *) (0));

  // thickness pass
  m_stats.begin("Thickness pass");
  OGLF->glBindFramebuffer(GL_FRAMEBUFFER, m_fbo_thickness);
  OGLF->glViewport(0, 0, m_w / m_thickness_scale, m_h / m_thickness_scale);
  OGLF->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  //setPerspectiveProjection(m_w / m_thickness_scale, m_h / m_thickness_scale);
  renderThickness(mv, part_cnt);
  //setPerspectiveProjection(m_w, m_h);
  OGLF->glViewport(0, 0, m_w, m_h);
  m_stats.end("Thickness pass");

  OGLF->glBindFramebuffer(GL_FRAMEBUFFER, m_fbo_depth);
  OGLF->glClearColor(0.0f, 0.0f, 0.0f, 0.0f);

  // depth pass
  m_stats.begin("Depth pass");
  OGLF->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  renderDepth(mv, part_cnt);
  m_stats.end("Depth pass");

  // depth smoothing pass (uses the fbo_depth)
  if (m_display_mode != MODE_Depth)
  {
    m_stats.begin("Smoothing pass");
    smoothDepth();
    m_stats.end("Smoothing pass");
  }

  // final composition pass
  m_stats.begin("Composition pass");
  OGLF->glBindFramebuffer(GL_FRAMEBUFFER, m_default_fbo);
  QMatrix4x4 camera_mv;
  camera_mv.rotate(camera_rotation);
  renderComposed(mv, camera_mv);
  m_stats.end("Composition pass");

  OGLF->glUseProgram(0);

  OGLF->glBindBuffer(GL_ARRAY_BUFFER, 0);

  m_vao.release();

  m_stats.end("curvature_flow_renderer_total_time");
}


void CurvatureFlowRenderer::renderThickness(const QMatrix4x4 & mv, size_t part_cnt)
{
  // pick shader program
  GLuint shader_id = m_prog_thickness.programId();
  OGLF->glUseProgram(shader_id);

  QMatrix4x4 proj;
  proj.perspective(30.0f, ((float) (m_w / m_thickness_scale)) / ((float) (m_h / m_thickness_scale)), 5.0f, 1000.0f);

  // set uniforms
  OGLF->glUniform1i(OGLF->glGetUniformLocation(shader_id, "use_uniform_color"), 1);
  OGLF->glUniform3f(OGLF->glGetUniformLocation(shader_id, "particle_uniform_color"),
                    m_particle_col.x(), m_particle_col.y(), m_particle_col.z());
  //OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "proj"), 1, GL_FALSE, m_proj.constData());
  OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "proj"), 1, GL_FALSE, proj.constData());
  OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "mv"), 1, GL_FALSE, mv.constData());
  OGLF->glUniform2f(OGLF->glGetUniformLocation(shader_id, "screen_size"),
                    m_w / m_thickness_scale, m_h / m_thickness_scale);

  // set rendering options
  OGLF->glEnable(GL_POINT_SPRITE);
  //OGLF->glPointSize(particle_radius * 2.0f);
  OGLF->glEnable(GL_PROGRAM_POINT_SIZE);
  OGLF->glEnable(GL_BLEND);
  OGLF->glBlendFunc(GL_ONE, GL_ONE);
  OGLF->glDisable(GL_DEPTH_TEST);

  // render point sprites
  OGLF->glDrawArrays(GL_POINTS, 0, part_cnt);

  // revert rendering options
  OGLF->glEnable(GL_DEPTH_TEST);
  //OGLF->glBlendFunc(GL_ONE, GL_ZERO);   // the default blending function according to manual
  OGLF->glDisable(GL_BLEND);
  OGLF->glDisable(GL_PROGRAM_POINT_SIZE);
  OGLF->glDisable(GL_POINT_SPRITE);
}


void CurvatureFlowRenderer::renderDepth(const QMatrix4x4 & mv, size_t part_cnt)
{
  // pick shader program
  GLuint shader_id = m_prog_depth.programId();
  OGLF->glUseProgram(shader_id);

  // set uniforms
  OGLF->glUniform1i(OGLF->glGetUniformLocation(shader_id, "use_uniform_color"), 1);
  OGLF->glUniform3f(OGLF->glGetUniformLocation(shader_id, "particle_uniform_color"),
                    m_particle_col.x(), m_particle_col.y(), m_particle_col.z());
  OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "proj"), 1, GL_FALSE, m_proj.constData());
  OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "mv"), 1, GL_FALSE, mv.constData());
  OGLF->glUniform2f(OGLF->glGetUniformLocation(shader_id, "screen_size"), m_w, m_h);

  // set rendering options
  OGLF->glEnable(GL_POINT_SPRITE);
  //OGLF->glPointSize(particle_radius * 2.0f);
  OGLF->glEnable(GL_PROGRAM_POINT_SIZE);
  OGLF->glEnable(GL_DEPTH_TEST);
  OGLF->glDisable(GL_BLEND);

  // render point sprites
  OGLF->glDrawArrays(GL_POINTS, 0, part_cnt);

  // revert rendering options
  OGLF->glEnable(GL_BLEND);
  OGLF->glDisable(GL_DEPTH_TEST);
  OGLF->glDisable(GL_PROGRAM_POINT_SIZE);
  OGLF->glDisable(GL_POINT_SPRITE);
}


void CurvatureFlowRenderer::smoothDepth(void)
{
  // pick shader program
  GLuint shader_id = m_prog_smooth.programId();
  OGLF->glUseProgram(shader_id);

  // set uniforms
  OGLF->glUniform1i(OGLF->glGetUniformLocation(shader_id, "tex_particle_depth"), 0);
  // OGLF->glUniform2f(OGLF->glGetUniformLocation(shader_id, "screen_size"), m_w, m_h);
  // OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "proj"), 1, GL_FALSE, m_proj.constData());
  // //OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "proj"), 1, GL_FALSE, m_proj.inverted().constData());

  OGLF->glUniform2f(OGLF->glGetUniformLocation(shader_id, "dx"), (GLfloat) (1.0f / m_w), (GLfloat) (0.0f));
  OGLF->glUniform2f(OGLF->glGetUniformLocation(shader_id, "dy"), (GLfloat) (0.0f), (GLfloat) (1.0f / m_h));

  GLfloat cx = -2.0f / (m_w * m_proj(0, 0));
  GLfloat cy = -2.0f / (m_h * m_proj(1, 1));
  OGLF->glUniform1f(OGLF->glGetUniformLocation(shader_id, "cx"), cx);
  OGLF->glUniform1f(OGLF->glGetUniformLocation(shader_id, "cy"), cy);
  OGLF->glUniform1f(OGLF->glGetUniformLocation(shader_id, "cx2"), cx * cx);
  OGLF->glUniform1f(OGLF->glGetUniformLocation(shader_id, "cy2"), cy * cy);
  OGLF->glUniform1f(OGLF->glGetUniformLocation(shader_id, "cx2cy2"), cx * cx * cy * cy);

  // smooth the generated depth
  OGLF->glDisable(GL_DEPTH_TEST);
  OGLF->glDisable(GL_BLEND);

  // smoothing loop
  for (int i = 0; i < m_smooth_iter_cnt; ++i)
  {
    // activate and bind the previous depth texture
    OGLF->glActiveTexture(GL_TEXTURE0);
    OGLF->glBindTexture(GL_TEXTURE_2D, activeDepthTexture());

    // setup which color attachement will be written into
    OGLF->glDrawBuffer(drawBufferAttachement());  // set which attachements will be used for drawing

    // draw a screen space quad to run the smoothing shader on every pixel of the screen
    OGLF->glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    // switch buffers
    swapActiveDepthTextures();
  }

  OGLF->glEnable(GL_BLEND);
  OGLF->glEnable(GL_DEPTH_TEST);
}


void CurvatureFlowRenderer::renderComposed(const QMatrix4x4 & mv, const QMatrix4x4 & camera_rot)
{
  // this function renders a screen space quad and composes the results
  // of previous passes into final result

  // bind necessary textures
  OGLF->glActiveTexture(GL_TEXTURE0);
  OGLF->glBindTexture(GL_TEXTURE_2D, activeDepthTexture());
  //OGLF->glBindTexture(GL_TEXTURE_2D, inactiveDepthTexture());

  OGLF->glActiveTexture(GL_TEXTURE1);
  OGLF->glBindTexture(GL_TEXTURE_2D, m_tex_thickness.textureId());

  OGLF->glActiveTexture(GL_TEXTURE2);
  bindSkyBoxTexture();

  // pick shader program
  GLuint shader_id = m_prog_compose.programId();
  OGLF->glUseProgram(shader_id);

  // set uniforms
  OGLF->glUniform1i(OGLF->glGetUniformLocation(shader_id, "tex_particle_depth"), 0);
  OGLF->glUniform1i(OGLF->glGetUniformLocation(shader_id, "tex_particle_thickness"), 1);
  OGLF->glUniform1i(OGLF->glGetUniformLocation(shader_id, "tex_sky_box_cubemap"), 2);
  OGLF->glUniform2f(OGLF->glGetUniformLocation(shader_id, "screen_size"), m_w, m_h);
  OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "proj"), 1, GL_FALSE, m_proj.constData());
  //OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "proj"), 1, GL_FALSE, m_proj.inverted().constData());
  OGLF->glUniformMatrix4fv(OGLF->glGetUniformLocation(shader_id, "mv"), 1, GL_FALSE, mv.constData());
  OGLF->glUniformMatrix3fv(OGLF->glGetUniformLocation(shader_id, "camera_mv"), 1, GL_FALSE,
                           camera_rot.toGenericMatrix<3,3>().constData());

  // ak bude posledny stlpec a posledny riadok mv matice [0,0,0,1], tak potom by malo platit,
  // ze mat3(inv(mv)) == inv(mat3(mv)), pre potvrdenie pozri vzorce tu:
  // http://www.cg.info.hiroshima-cu.ac.jp/~miyazaki/knowledge/teche23.html
  OGLF->glUniformMatrix3fv(OGLF->glGetUniformLocation(shader_id, "mv_inv"), 1, GL_FALSE,
                           mv.inverted().toGenericMatrix<3,3>().constData());

  //qDebug() << "proj=" << m_proj;

  OGLF->glUniform1f(OGLF->glGetUniformLocation(shader_id, "cx"), (GLfloat) (-2.0f / (m_w * m_proj(0, 0))));
  OGLF->glUniform1f(OGLF->glGetUniformLocation(shader_id, "cy"), (GLfloat) (-2.0f / (m_h * m_proj(1, 1))));
  //OGLF->glUniform1f(OGLF->glGetUniformLocation(shader_id, "cx"), (GLfloat) (2.0f / (m_w * m_proj(0, 0))));
  //OGLF->glUniform1f(OGLF->glGetUniformLocation(shader_id, "cy"), (GLfloat) (2.0f / (m_h * m_proj(1, 1))));
  OGLF->glUniform2f(OGLF->glGetUniformLocation(shader_id, "dx"), (GLfloat) (1.0f / m_w), (GLfloat) (0.0f));
  OGLF->glUniform2f(OGLF->glGetUniformLocation(shader_id, "dy"), (GLfloat) (0.0f), (GLfloat) (1.0f / m_h));
  OGLF->glUniform1f(OGLF->glGetUniformLocation(shader_id, "proj_00_inv"), (GLfloat) (1.0f / m_proj(0, 0)));
  OGLF->glUniform1f(OGLF->glGetUniformLocation(shader_id, "proj_11_inv"), (GLfloat) (1.0f / m_proj(1, 1)));

  OGLF->glUniform1i(OGLF->glGetUniformLocation(shader_id, "display_mode"), m_display_mode);

  // render the quad
  OGLF->glEnable(GL_BLEND);
  OGLF->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
  //OGLF->glEnable(GL_DEPTH_TEST);
  //OGLF->glDepthMask(GL_FALSE);
  //OGLF->glDisable(GL_BLEND);
  OGLF->glDisable(GL_DEPTH_TEST);
  OGLF->glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
  //OGLF->glDepthMask(GL_TRUE);
  //OGLF->glDisable(GL_DEPTH_TEST);
  //OGLF->glDisable(GL_BLEND);

  // unbind the textures
  OGLF->glActiveTexture(GL_TEXTURE1);
  OGLF->glBindTexture(GL_TEXTURE_2D, 0);

  OGLF->glActiveTexture(GL_TEXTURE0);
  OGLF->glBindTexture(GL_TEXTURE_2D, 0);
}


bool CurvatureFlowRenderer::resetDepthFBO(int w, int h)
{
  // TODO: 2 rozne textury a 1 render buffer pre depth fbo
  // je mozno az prilis vela.
  // Ako optimalizaciu mozno skusit len dve floatove textury
  // a v depth shaderi zapisovat len do gl_FragDepth s tym, ze
  // pri kazdom prehodeni textur pri smoothovani bude treba
  // zavolat glFramebufferTexture2D nastavit spravny attachement.

  // allocate render buffer
  if (m_rbo_depth == 0) OGLF->glGenRenderbuffers(1, &m_rbo_depth);
  OGLF->glBindRenderbuffer(GL_RENDERBUFFER, m_rbo_depth);
  OGLF->glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH_COMPONENT24, w, h);

  // allocate textures
  if (m_tex_depth[0] == 0) OGLF->glGenTextures(2, m_tex_depth);

  OGLF->glBindTexture(GL_TEXTURE_2D, m_tex_depth[0]);
  // GL_R32F by mi malo dat floatove textury, ktore nie su clampovane do [0,1]
  OGLF->glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, nullptr);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  //OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);

  OGLF->glBindTexture(GL_TEXTURE_2D, m_tex_depth[1]);
  OGLF->glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, nullptr);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
  //OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAX_LEVEL, 0);

  OGLF->glBindTexture(GL_TEXTURE_2D, 0);

  // Setup framebuffer
  if (m_fbo_depth == 0) OGLF->glGenFramebuffers(1, &m_fbo_depth);

  OGLF->glBindFramebuffer(GL_FRAMEBUFFER, m_fbo_depth);
  OGLF->glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_tex_depth[0], 0);
  OGLF->glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT1, GL_TEXTURE_2D, m_tex_depth[1], 0);
  OGLF->glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_ATTACHMENT, GL_RENDERBUFFER, m_rbo_depth);

  //OGLF->glDrawBuffer(GL_COLOR_ATTACHMENT0);  // set which attachements will be drawn to
  OGLF->glDrawBuffer(drawBufferAttachement());  // set which attachements will be drawn to

  if (OGLF->glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
  {
    ERRORM("Failed to create frame buffer: frame buffer is not complete");
    return false;
  }

  OGLF->glBindFramebuffer(GL_FRAMEBUFFER, m_default_fbo);

  return true;
}


bool CurvatureFlowRenderer::resetThicknessFBO(int w, int h)
{
  // alokacia textury
  if (!m_tex_thickness.create())
  {
    ERRORM("Failed to create texture for thickness FBO");
    return false;
  }

  m_tex_thickness.bind();

  //OGLF->glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, w, h, 0, GL_RGBA, GL_FLOAT, nullptr);
  // GL_R32F by mi malo dat floatove textury, ktore nie su clampovane do [0,1]
  OGLF->glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F,
                     w / m_thickness_scale, h / m_thickness_scale,
                     0, GL_RED, GL_FLOAT, nullptr);
  //OGLF->glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, w, h, 0, GL_RED, GL_FLOAT, nullptr);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  //OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_BORDER);  // ???
  //OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_BORDER);  // ???

  m_tex_thickness.release();

  // Nastavenie framebufferu
  if (m_fbo_thickness == 0) OGLF->glGenFramebuffers(1, &m_fbo_thickness);

  OGLF->glBindFramebuffer(GL_FRAMEBUFFER, m_fbo_thickness);

  OGLF->glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, m_tex_thickness.textureId(), 0);

  // Nastavenie zoznamu attachementov, do ktorych sa bude kreslit
  GLenum draw_buffers[] = { GL_COLOR_ATTACHMENT0 };
  OGLF->glDrawBuffers(1, draw_buffers);

  //glDrawBuffer(GL_NONE);

  if (OGLF->glCheckFramebufferStatus(GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE)
  {
    ERRORM("Failed to create frame buffer: frame buffer is not complete");
    return false;
  }

  OGLF->glBindFramebuffer(GL_FRAMEBUFFER, m_default_fbo);

  return true;
}


bool CurvatureFlowRenderer::init(void)
{
  // depth program
  if (!utils::ogl::buildShaderProgram(m_prog_depth,
                                      ":/src/opengl/curvature_flow_renderer.vert",
                                      ":/src/opengl/curvature_flow_renderer_depth.frag"))
  {
    ERRORM("Failed to compile depth pass program");
    return false;
  }

  // thickness program
  if (!utils::ogl::buildShaderProgram(m_prog_thickness,
                                      ":/src/opengl/curvature_flow_renderer.vert",
                                      ":/src/opengl/curvature_flow_renderer_thickness.frag"))
  {
    ERRORM("Failed to compile thickness pass program");
    return false;
  }

  // program to smooth out depth buffer
  if (!utils::ogl::buildShaderProgram(m_prog_smooth,
                                      ":/src/opengl/curvature_flow_renderer_quad.vert",
                                      ":/src/opengl/curvature_flow_renderer_curvatureflow.frag"))
  {
    ERRORM("Failed to compile curvature flow smoothing pass program");
    return false;
  }

  // program to combine previous stages
  if (!utils::ogl::buildShaderProgram(m_prog_compose,
                                      ":/src/opengl/curvature_flow_renderer_quad.vert",
                                      ":/src/opengl/curvature_flow_renderer_compose.frag"))
  {
    ERRORM("Failed to compile final composition program");
    return false;
  }

  // Create vertex array object
  return m_vao.create();
}
