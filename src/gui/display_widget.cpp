#include "gui/display_widget.h"

#include <QWheelEvent>
#include <QKeyEvent>
#include <QMatrix4x4>
#include <QPainter>
#include <QApplication>



DisplayWidget::~DisplayWidget(void)
{
  makeCurrent(); // to have an active OpenGL context in renderer destructor

  OGLF->glDeleteQueries(1, &m_gl_timer_id);
}


QString DisplayWidget::particleSystemTypeToString(ParticleSystemType t)
{
  switch (t)
  {
    case PS_Debug:          return tr("Debugging particle system");
    case PS_SPHNaive:       return tr("Naive SPH simulator");
    case PS_SPHOptimized:   return tr("Optimized SPH simulator");
    case PS_SPHUniformGrid: return tr("Uniform Grid SPH simulator");
    default:                return tr("Unknown particle system");
  }

  return tr("Unknown particle system");
}


QString DisplayWidget::rendererTypeToString(RendererType t)
{
  switch (t)
  {
    case REN_Debug:         return tr("Debugging renderer");
    case REN_Instancing:    return tr("Instancing renderer");
    case REN_PointSprite:   return tr("Point Sprite renderer");
    case REN_CurvatureFlow: return tr("Screen Space Fluid renderer with Curvature Flow");
    case REN_Texture:       return tr("Texture based renderer");
    case REN_RayCast:       return tr("Ray casting based renderer");
    default:                return tr("Unknown renderer");
  }

  return tr("Unknown renderer");
}


void DisplayWidget::setSkyBox(SkyBoxType type)
{
  m_cur_skybox_type = type;

  switch (type)
  {
    case SKY_Forest:
      m_cur_ren->setSkyBoxTextures(":/data/skybox/forest_posx.jpg",
                                   ":/data/skybox/forest_negx.jpg",
                                   ":/data/skybox/forest_posy.jpg",
                                   ":/data/skybox/forest_negy.jpg",
                                   ":/data/skybox/forest_posz.jpg",
                                   ":/data/skybox/forest_negz.jpg");
      break;

    case SKY_PalmTrees:
      m_cur_ren->setSkyBoxTextures(":/data/skybox/palmtrees_posx.jpg",
                                   ":/data/skybox/palmtrees_negx.jpg",
                                   ":/data/skybox/palmtrees_posy.jpg",
                                   ":/data/skybox/palmtrees_negy.jpg",
                                   ":/data/skybox/palmtrees_posz.jpg",
                                   ":/data/skybox/palmtrees_negz.jpg");
      break;

    case SKY_Clouds:
      m_cur_ren->setSkyBoxTextures(":/data/skybox/clouds_posx.bmp",
                                   ":/data/skybox/clouds_negx.bmp",
                                   ":/data/skybox/clouds_posy.bmp",
                                   ":/data/skybox/clouds_negy.bmp",
                                   ":/data/skybox/clouds_posz.bmp",
                                   ":/data/skybox/clouds_negz.bmp");
      break;

    case SKY_Debug:
      m_cur_ren->setSkyBoxTextures(":/data/skybox/debug_posx.png",
                                   ":/data/skybox/debug_negx.png",
                                   ":/data/skybox/debug_posy.png",
                                   ":/data/skybox/debug_negy.png",
                                   ":/data/skybox/debug_posz.png",
                                   ":/data/skybox/debug_negz.png");
      break;

    case SKY_DebugLarge:
      m_cur_ren->setSkyBoxTextures(":/data/skybox/debug_large_posx.png",
                                   ":/data/skybox/debug_large_negx.png",
                                   ":/data/skybox/debug_large_posy.png",
                                   ":/data/skybox/debug_large_negy.png",
                                   ":/data/skybox/debug_large_posz.png",
                                   ":/data/skybox/debug_large_negz.png");
      break;

    default:
      WARNM("Unknown SkyBox type: " << int(type));
      return;
  }
}


void DisplayWidget::toggleParticleSystem(void)
{
  setParticleSystem_impl(toggleParticleSystemType());
  update();
}


void DisplayWidget::setParticleSystem(ParticleSystemType type)
{
  setParticleSystem_impl(type);
  update();
}


void DisplayWidget::setParticleSystem_impl(ParticleSystemType type)
{
  makeCurrent();

  switch (type)
  {
    case PS_Debug:
      qDebug() << "Using debug particle system";
      m_cur_ps = m_debug_system.get();
      break;

    case PS_SPHNaive:
      qDebug() << "Using naive sph fluid particle system";
      m_cur_ps = m_sph_naive.get();
      break;

    case PS_SPHOptimized:
      qDebug() << "Using optimized sph fluid particle system";
      m_cur_ps = m_sph_optimized.get();
      break;

    case PS_SPHUniformGrid:
      qDebug() << "Using sph fluid particle system with uniform grid";
      m_cur_ps = m_sph_uniform_grid.get();
      break;

    default:
      qDebug() << "Unknown particle system: " << type;
      // the return here is very important, otherwise
      // the code under this switch will break
      return;
  }

  m_cur_ps_type = type;

  particleSystemReset();
}


void DisplayWidget::toggleRenderer(void)
{
  setRenderer_impl(toggleRendererType());
}


void DisplayWidget::setRenderer(RendererType type)
{
  setRenderer_impl(type);
  update();
}


void DisplayWidget::setRenderer_impl(RendererType type)
{
  makeCurrent();

  switch (type)
  {
    case REN_Debug:
      qDebug() << "Using debug renderer";
      break;

    case REN_Instancing:
      qDebug() << "Using instancing renderer";
      m_cur_ren = m_ren_instancing.get();
      break;

    case REN_PointSprite:
      qDebug() << "Using point sprite renderer";
      m_cur_ren = m_ren_point_sprite.get();
      break;

    case REN_CurvatureFlow:
      qDebug() << "Using Screen Space renderer with Curvature Flow";
      m_cur_ren = m_ren_curvature_flow.get();
      break;

    case REN_Texture:
      qDebug() << "Using texture renderer";
      break;

    case REN_RayCast:
      qDebug() << "Using raycasting renderer";
      m_cur_ren = m_ren_ray_cast.get();
      break;

    default:
      qDebug() << "Unknown renderer: " << type;
      // the return here is very important, otherwise
      // the code under this switch will break
      return;
  }

  m_cur_ren_type = type;

  rendererReset();
}


io::VoxelMesh *DisplayWidget::loadVoxelMesh(io::VoxelMesh & mesh)
{
  const size_t w = m_user_grid_w;
  const size_t h = m_user_grid_h;
  const size_t d = m_user_grid_d;
  QString voxel_mesh_file(":/data/voxel_meshes/%1_%2.raw");

  switch (m_reset_mode)
  {
    case RM_ArmadilloFilled:
      voxel_mesh_file = voxel_mesh_file.arg("armadillo_full");
      break;

    case RM_ArmadilloHollow:
      voxel_mesh_file = voxel_mesh_file.arg("armadillo_hollow");
      break;

    case RM_BunnyFilled:
      voxel_mesh_file = voxel_mesh_file.arg("bunny_full");
      break;

    case RM_BunnyHollow:
      voxel_mesh_file = voxel_mesh_file.arg("bunny_hollow");
      break;

    case RM_Prism:
      m_sph_params->setParticleCount(m_user_particle_count);
      return nullptr;

    default:
      WARNM("Unknown particle position reset mode");
      return nullptr;
  }

  bool ret = true;

  if ((w >= 128) && (h >= 128) && (d >= 128))
    ret = mesh.loadFromRaw(voxel_mesh_file.arg("256x256x256").toStdString().c_str(), 256, 256, 256);
  else if ((w >= 64) && (h >= 64) && (d >= 64))
    ret = mesh.loadFromRaw(voxel_mesh_file.arg("128x128x128").toStdString().c_str(), 128, 128, 128);
  else if ((w >= 32) && (h >= 32) && (d >= 32))
    ret = mesh.loadFromRaw(voxel_mesh_file.arg("64x64x64").toStdString().c_str(), 64, 64, 64);
  else
  {
    WARNM("No suitable voxel mesh on scale " << w << "x" << h << "x" << d);
    return nullptr;
  }

  if (!ret)
  {
    WARNM("Failed to load the voxel mesh from " << voxel_mesh_file.toStdString());
    return nullptr;
  }

  //mesh.dump(std::cout);

  return &mesh;
}


bool DisplayWidget::particleSystemReset(void)
{
  io::VoxelMesh mesh;
  if (!m_cur_ps->reset(loadVoxelMesh(mesh)))
  {
    emit error(tr("Failed to initialize particle system"));
    return false;
  }

  // set renderer coloring mode
  if (m_cur_ren) m_cur_ren->setUseUniformColor(m_cur_ps_type != PS_Debug);

  // set grid and set volume size
  m_ren_ray_cast->setGrid(m_sph_uniform_grid->cellStartsBuffer(),
                          m_sph_uniform_grid->cellEndsBuffer());

  m_ren_ray_cast->setGridSize(m_sph_params->gridWidth(),
                              m_sph_params->gridHeight(),
                              m_sph_params->gridDepth());

  return true;
}


bool DisplayWidget::rendererReset(void)
{
  // set grid and set volume size
  m_ren_ray_cast->setGrid(m_sph_uniform_grid->cellStartsBuffer(),
                          m_sph_uniform_grid->cellEndsBuffer());

  m_ren_ray_cast->setGridSize(m_sph_params->gridWidth(),
                              m_sph_params->gridHeight(),
                              m_sph_params->gridDepth());

  // reset the renderer
  if (!m_cur_ren->reset(width(), height()))
  {
    emit error(tr("Failed to initialize renderer"));
    return false;
  }

  //m_cur_ren->setUseUniformColor(false);
  m_cur_ren->setUseUniformColor(m_cur_ps_type != PS_Debug);
  m_cur_ren->setUseLighting(true);
  //m_cur_ren->setLightPosition(QVector3D(0.0f, 300.0f, 300.0f));
  //m_cur_ren->setLightPosition(QVector3D(-300.0f, -300.0f, -300.0f));
  m_cur_ren->setLightPosition(QVector3D(300.0f, 300.0f, 300.0f));
  //m_cur_ren->setLightPosition(QVector3D(300.0f, 300.0f, -300.0f));
  m_cur_ren->setLightAmbientColor(QVector3D(0.2f, 0.2f, 0.2f));
  m_cur_ren->setLightDiffuseColor(QVector3D(0.8f, 0.8f, 0.8f));

  //OGLF->glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  //OGLF->glClearColor(0.8f, 0.8f, 0.8f, 1.0f);
  //OGLF->glClearColor(0.5f, 0.0f, 0.0f, 1.0f);
  m_cur_ren->setClearColor(QVector4D(0.8f, 0.8f, 0.8f, 1.0f));
  m_cur_ren->setDrawSkyBox(true);

  setSkyBox(m_cur_skybox_type);

  return true;
}


bool DisplayWidget::simulationReset(void)
{
#if 0
  if (m_mesh_armadillo.isValid())
  {
    if (!m_cur_ps->reset(m_mesh_armadillo))
    {
      emit error(tr("Failed to reinitialize particle system"));
      return false;
    }
  }
  else
  {
    //if (!m_cur_ps->reset(DEFAULT_PARTICLE_COUNT))
    //if (!m_cur_ps->reset(m_particle_count))
    if (!m_cur_ps->reset())
    {
      emit error(tr("Failed to reinitialize particle system"));
      return false;
    }
  }

  // set grid and set volume size
  m_ren_ray_cast->setGrid(m_sph_uniform_grid->cellStartsBuffer(),
                          m_sph_uniform_grid->cellEndsBuffer());

  m_ren_ray_cast->setGridSize(m_sph_params->gridWidth(),
                              m_sph_params->gridHeight(),
                              m_sph_params->gridDepth());

  if (!m_cur_ren->reset(width(), height()))
  {
    emit error(tr("Failed to initialize renderer"));
    return false;
  }

  // set renderer coloring mode
  if (m_cur_ren) m_cur_ren->setUseUniformColor(m_cur_ps_type != PS_Debug);
  return true;
#else
  return particleSystemReset() && rendererReset();
#endif
}


void DisplayWidget::resetXYPanning(void)
{
  m_pan_pos = QPoint();
  m_pan_dx_prev = 0.0f;
  m_pan_dy_prev = 0.0f;
  m_pan_dx = 0.0f;
  m_pan_dy = 0.0f;
}


void DisplayWidget::resetZTranslation(void)
{
  m_z_dist_pos = 0;
  m_z_dist = 0.0f;
  m_z_dist_prev = 0.0f;
}


void DisplayWidget::setDetail(int level)
{
  m_detail = level;
  update();
}


void DisplayWidget::setAutoSubsampling(bool enabled)
{
  m_auto_subsampling = enabled;
  update();
}


void DisplayWidget::setHighQuality(void)
{
  m_high_quality = true;
  update();
}


void DisplayWidget::setLowQuality(void)
{
  m_high_quality = false;
  update();
}


void DisplayWidget::setDisplayBBox(bool enabled)
{
  assert(m_cur_ren != nullptr);
  m_cur_ren->setDrawBoundingVolume(enabled);
  update();
}


void DisplayWidget::setUseShading(bool enabled)
{
  assert(m_cur_ren != nullptr);
  m_cur_ren->setUseLighting(enabled);
  update();
}


void DisplayWidget::setLightPosition(const QVector3D & light_pos)
{
  assert(m_cur_ren != nullptr);
  m_cur_ren->setLightPosition(light_pos);
  update();
}


void DisplayWidget::setLightAmbientColor(const QVector3D & ambient_col)
{
  assert(m_cur_ren != nullptr);
  m_cur_ren->setLightAmbientColor(ambient_col);
  update();
}


void DisplayWidget::setLightDiffuseColor(const QVector3D & diffuse_col)
{
  assert(m_cur_ren != nullptr);
  m_cur_ren->setLightDiffuseColor(diffuse_col);
  update();
}


void DisplayWidget::initializeGL(void)
{
  qDebug() << __PRETTY_FUNCTION__;

  // inicializacia OpenGL Funkcii
  //OGLF->initializeOpenGLFunctions();

  // inicializacia logovania sprav z debug kontextu
  m_logger.init();

  // initialize OpenCL context
  if (!utils::ocl::initCLGLContext(m_cl_ctx))
  {
    emit error(tr("Failed to initialize OpenCL context"));
    return;
  }

  // initialize SPH parameters and set some reasonable defaults
  m_sph_params.reset(new SPHParams(m_cl_ctx));
  //m_sph_params->setParticleCount(m_particle_count);
  m_sph_params->setParticleCount(DEFAULT_PARTICLE_COUNT);
  //m_sph_params->setGridSize({ 512, 512, 512, 1});
  //m_sph_params->setGridSize({ 512, 128, 512, 1});
  //m_sph_params->setGridSize({ 512, 64, 512, 1});
  //m_sph_params->setGridSize({ 256, 256, 256, 1});
  //m_sph_params->setGridSize({ 256, 32, 256, 1});
  //m_sph_params->setGridSize({ 128, 128, 128, 1});
  //m_sph_params->setGridSize({ 128, 64, 128, 1});
  //m_sph_params->setGridSize({ 128, 32, 128, 1});
  //m_sph_params->setGridSize({ 64, 128, 64, 1});
  //m_sph_params->setGridSize({ 64, 32, 64, 1});
  //m_sph_params->setGridSize({ 32, 32, 32, 1});
  //m_sph_params->setGridSize({ 32, 128, 32, 1});
  //m_sph_params->setGridSize({ 32, 64, 32, 1});
  m_sph_params->setGridSize({ DEFAULT_GRID_W, DEFAULT_GRID_H, DEFAULT_GRID_D, 1});

  // override default SPH parameter settings by arguments from command line
  m_sph_params->parseCmdArgs(QApplication::arguments());

  m_user_grid_w = m_sph_params->gridWidth();
  m_user_grid_h = m_sph_params->gridHeight();
  m_user_grid_d = m_sph_params->gridDepth();
  m_user_particle_count = m_sph_params->particleCount();

  // create particle systems
  m_sph_naive.reset(new SPHNaive(m_cl_ctx, m_sph_params.get()));
  m_sph_optimized.reset(new SPHOptimized(m_cl_ctx, m_sph_params.get()));
  m_sph_uniform_grid.reset(new SPHUniformGrid(m_cl_ctx, m_sph_params.get()));
  m_debug_system.reset(new DebugParticleSystem(m_cl_ctx, m_sph_params.get()));

  // create renderers
  m_ren_instancing.reset(new InstancingRenderer);
  m_ren_point_sprite.reset(new PointSpriteRenderer);
  m_ren_curvature_flow.reset(new CurvatureFlowRenderer);
  m_ren_ray_cast.reset(new RayCastRenderer(m_cl_ctx, m_sph_uniform_grid->clCommandQueue()));

  // set default settings
  setParticleSystem_impl(PS_SPHUniformGrid);
  //setRenderer_impl(REN_RayCast);
  setRenderer_impl(REN_PointSprite);
  //setRenderer_impl(REN_CurvatureFlow);

#if 0
  //m_cur_ren->setUseUniformColor(false);
  m_cur_ren->setUseUniformColor(m_cur_ps_type != PS_Debug);
  m_cur_ren->setUseLighting(true);
  //m_cur_ren->setLightPosition(QVector3D(0.0f, 300.0f, 300.0f));
  //m_cur_ren->setLightPosition(QVector3D(-300.0f, -300.0f, -300.0f));
  m_cur_ren->setLightPosition(QVector3D(300.0f, 300.0f, 300.0f));
  //m_cur_ren->setLightPosition(QVector3D(300.0f, 300.0f, -300.0f));
  m_cur_ren->setLightAmbientColor(QVector3D(0.2f, 0.2f, 0.2f));
  m_cur_ren->setLightDiffuseColor(QVector3D(0.8f, 0.8f, 0.8f));

  //OGLF->glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
  //OGLF->glClearColor(0.8f, 0.8f, 0.8f, 1.0f);
  //OGLF->glClearColor(0.5f, 0.0f, 0.0f, 1.0f);
  m_cur_ren->setClearColor(QVector4D(0.8f, 0.8f, 0.8f, 1.0f));
  m_cur_ren->setDrawSkyBox(true);
#endif
  m_sim_frame_num = 0;

  OGLF->glGenQueries(1, &m_gl_timer_id);

  m_initialized = true;

  startTimer(33);
}


void DisplayWidget::paintGL(void)
{
  if (!m_initialized)
  {
    WARNM("Not rendering because DisplayWidget was not successfully initialized");
    return;
  }

#if 1
  auto start = std::chrono::high_resolution_clock::now();

  // aktualizacia simulacie
  if ((!m_pause) || (m_step_frame))
  {
    if (m_mixing)
      m_sph_uniform_grid->rotate(m_track_ball.getRotation());
    else
      m_sph_uniform_grid->rotate(QQuaternion());

    m_cur_ps->update();
    if (!m_freeze_timers) m_sim_time_ms = m_cur_ps->frameTime();
    ++m_sim_frame_num;
    if ((m_cur_ren == m_ren_ray_cast.get()) && (m_calc_normals))
    {
      m_ren_ray_cast->calcNormals();
    }
  }

  // vykreslenie sceny
  GLuint timer_ren_available = 1;

  if (m_frame_num > 0)
  {
    OGLF->glGetQueryObjectuiv(m_gl_timer_id, GL_QUERY_RESULT_AVAILABLE, &timer_ren_available);
  }

  if (timer_ren_available)
  {
    if ((m_frame_num > 0) && (m_freeze_timers == false))
    {
      GLuint elapsed_ns = 0;
      OGLF->glGetQueryObjectuiv(m_gl_timer_id, GL_QUERY_RESULT, &elapsed_ns);
      m_ren_time_ms = elapsed_ns / 1000000.0f;
    }

    OGLF->glBeginQuery(GL_TIME_ELAPSED, m_gl_timer_id);
  }

  //OGLF->glClearColor(0.8f, 0.8f, 0.8f, 1.0f);
  //OGLF->glClearColor(0.0f, 0.0f, 0.0f, 0.0f);
  //OGLF->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  m_cur_ren->setDefaultFBO(defaultFramebufferObject());
  m_cur_ren->setMixing(m_mixing);

  cl_float3 bvs = m_sph_params->boundingVolumeSize();
  m_cur_ren->render(m_track_ball.getRotation(),
                    //QVector3D(1.0f, 1.0f, 1.0f),
                    QVector3D(m_scale, m_scale, m_scale),
                    //QVector3D(0.0f, 0.0f, -32.0f),
                    //QVector3D(0.0f, 0.0f, -64),
                    QVector3D(m_pan_dx_prev + m_pan_dx,
                              m_pan_dy_prev + m_pan_dy,
                              m_z_dist_prev + m_z_dist - 64.0f),
                    m_track_ball_camera.getRotation(),
                    m_cur_ps->positionsVBO(),
                    m_cur_ps->colorsVBO(),
                    m_sph_params->particleCount(),
                    QVector3D(bvs.s[0], bvs.s[1], bvs.s[2]));

  /* display messages */
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_CULL_FACE);

  QPainter painter(this);

  //painter.beginNativePainting();
  drawHUD(painter);
  //painter.endNativePainting();

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

  if (timer_ren_available) OGLF->glEndQuery(GL_TIME_ELAPSED);

  if ((m_display_info) && (m_freeze_timers == false))
  {
    glFinish();
    auto end = std::chrono::high_resolution_clock::now();
    m_total_time_ms = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
  }

  m_step_frame = false;
  ++m_frame_num;
#elif 0
  // inicializacia rendereru (renderer treba inicializovat az vtedy ked uz
  // existuje platny OpenGL kontext.
  // A navyse medzi renderermi sa da aj prepinat ...)
  if (m_renderer_changed)
  {
    qDebug() << "width=" << width() << ", height=" << height();
    if (!m_renderer->reset(width(), height()))
    {
      emit error(tr("Failed to initialize renderer"));
      return;
    }

    m_renderer->setPerspectiveProjection(width(), height());

    if (m_renderer_type == RayCastRenderer)
    {
      static_cast<RayCastVolumeRenderer *>(m_renderer.get())->setDefaultFBO(defaultFramebufferObject());
    }

    m_renderer_changed = false;
  }

  // nastavenie properties rendereru, ktore sa mozu menit kazdy frame
  m_renderer->setRenderBBox(m_display_bbox);
  m_renderer->setUseLighting(m_use_shading);
  m_renderer->setLightPosition(m_light_pos);
  m_renderer->setLightAmbientColor(m_light_ambient_col);
  m_renderer->setLightDiffuseColor(m_light_diffuse_col);

  // renderovanie
  if ((m_high_quality) || (m_auto_subsampling == false))
  {
    m_renderer->render(m_track_ball.getRotation(),
                       QVector3D(m_scale, m_scale, m_scale),
                       m_transl,
                       m_detail);
  }
  else
  {
    m_renderer->renderPreview(m_track_ball.getRotation(),
                              QVector3D(m_scale, m_scale, m_scale),
                              m_transl);
  }
#elif 0
  auto start = std::chrono::high_resolution_clock::now();

  // aktualizacia simulacie
  if ((!m_pause) || (m_step_frame))
  {
    ++m_frame_num;
    m_cur_ps->update();
  }

  // vypocet transformacie sceny
  QMatrix4x4 mv;
  mv.translate(0.0f, 0.0f, -32.0f);
  mv.rotate(m_track_ball.getRotation());
  mv.scale(m_scale);

  QMatrix4x4 proj;
  proj.perspective(30.0f, ((float) width()) / ((float) height()), 0.01f, 1000.0f);

  // vykreslenie sceny
  OGLF->glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  m_cur_ps->render(mv, proj);

  /* display messages */
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_CULL_FACE);

  QPainter painter(this);

  //painter.beginNativePainting();
  drawHUD(painter);
  //painter.endNativePainting();

  glEnable(GL_DEPTH_TEST);
  glEnable(GL_CULL_FACE);

  if (m_display_info)
  {
    glFinish();

    auto end = std::chrono::high_resolution_clock::now();
    double t = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    //std::cout << "t=" << t << std::endl;

    //std::time_t ttp_start = std::chrono::high_resolution_clock::to_time_t(start);
    //std::cout << "start: " << std::ctime(&ttp_start) << std::endl;
    //std::cout << "start: " << ttp_start << std::endl;

    //std::time_t ttp_end = std::chrono::high_resolution_clock::to_time_t(end);
    //std::cout << "end: " << std::ctime(&ttp_end) << std::endl;
    //std::cout << "end: " << ttp_end << std::endl;

    m_fps = 1000.0 / t;
  }

  m_step_frame = false;
#endif
}


void DisplayWidget::resizeGL(int w, int h)
{
  qDebug() << __PRETTY_FUNCTION__;
  if (m_cur_ren != nullptr)
  {
    if (!m_cur_ren->resize(w, h))
    {
      WARNM("Failed to resize render");
    }
  }
}


void DisplayWidget::keyPressEvent(QKeyEvent *event)
{
  if (!m_initialized)
  {
    WARNM("Not handling keyboard because DisplayWidget was not successfully initialized");
    return QOpenGLWidget::keyPressEvent(event);
  }

  switch (event->key())
  {
    //case Qt::Key_Shift: m_shift_pressed = true; break;

    case Qt::Key_4:     setSkyBox(SKY_Debug);                             break;
    case Qt::Key_5:     setSkyBox(SKY_DebugLarge);                        break;
    case Qt::Key_7:     setSkyBox(SKY_Forest);                            break;
    case Qt::Key_8:     setSkyBox(SKY_PalmTrees);                         break;
    case Qt::Key_9:     setSkyBox(SKY_Clouds);                            break;
    //case Qt::Key_D:     m_sph_params->toggleDrain();               break;
    case Qt::Key_F:     m_sph_params->toggleFountain();                   break;
    case Qt::Key_W:     m_sph_params->emitWave();                         break;
    case Qt::Key_H:     m_display_help = !m_display_help;                 break;
    case Qt::Key_I:     m_display_info = !m_display_info;                 break;
    case Qt::Key_Space: m_pause = !m_pause;                               break;
    case Qt::Key_N:     m_step_frame = true;                              break;
    case Qt::Key_Z:     m_freeze_timers = !m_freeze_timers;               break;
    case Qt::Key_G:     m_calc_normals = !m_calc_normals;                 break;
    case Qt::Key_M:     m_mixing = !m_mixing;                             break;
    //case Qt::Key_C:     m_ren_curvature_flow->toggleDisplayMode(); break;
    case Qt::Key_C:     m_ren_curvature_flow->toggleCurvatureFlowSteps(); break;

    case Qt::Key_P:
        if (m_ren_curvature_flow->displayMode() ==
            CurvatureFlowRenderer::MODE_PureRefraction)
        {
          std::cerr << "Other display modes" << std::endl;
          m_ren_curvature_flow->setDisplayMode(
                CurvatureFlowRenderer::MODE_Lighting);
        }
        else
        {
          std::cerr << "Pure reflection mode" << std::endl;
          m_ren_curvature_flow->setDisplayMode(
                CurvatureFlowRenderer::MODE_PureRefraction);
        }
      break;

    case Qt::Key_O:
      m_track_ball.reset();
      m_scale = DEFAULT_SCALE;
      m_track_ball_camera.reset();
      resetXYPanning();
      resetZTranslation();
      break;

    case Qt::Key_Left:
      m_track_ball.push(QPointF( 0.00f, 0.0f));
      m_track_ball.move(QPointF(-0.05f, 0.0f));
      m_track_ball_camera.push(QPointF( 0.00f, 0.0f));
      m_track_ball_camera.move(QPointF(-0.05f, 0.0f));
      break;

    case Qt::Key_Right:
      m_track_ball.push(QPointF(0.00f, 0.0f));
      m_track_ball.move(QPointF(0.05f, 0.0f));
      m_track_ball_camera.push(QPointF(0.00f, 0.0f));
      m_track_ball_camera.move(QPointF(0.05f, 0.0f));
      break;

    case Qt::Key_Down:
      m_track_ball.push(QPointF(0.0f,  0.00f));
      m_track_ball.move(QPointF(0.0f, -0.05f));
      m_track_ball_camera.push(QPointF(0.0f,  0.00f));
      m_track_ball_camera.move(QPointF(0.0f, -0.05f));
      break;

    case Qt::Key_Up:
      m_track_ball.push(QPointF(0.0f, 0.00f));
      m_track_ball.move(QPointF(0.0f, 0.05f));
      m_track_ball_camera.push(QPointF(0.0f, 0.00f));
      m_track_ball_camera.move(QPointF(0.0f, 0.05f));
      break;

    case Qt::Key_T:
      setParticleSystem_impl(toggleParticleSystemType());
      break;

    case Qt::Key_V:
      setRenderer_impl(toggleRendererType());
      m_cur_ren->resize(width(), height());
      break;

    case Qt::Key_S:
      if (event->modifiers() & Qt::ControlModifier) m_debug_system->toggleSpiral();
      break;

    case Qt::Key_B:
      INFOM("Bounding volume: " << (m_cur_ren->toggleDrawBoundingVolume() ? "on" : "off"));
      break;

    case Qt::Key_R:
      //simulationReset();
      particleSystemReset();
      break;

    case Qt::Key_1:
      m_reset_mode = RM_Prism;
      particleSystemReset();
      break;

    case Qt::Key_2:
      m_reset_mode = (event->modifiers() & Qt::ControlModifier) ? RM_ArmadilloFilled : RM_ArmadilloHollow;
      particleSystemReset();
      break;

    case Qt::Key_3:
      m_reset_mode = (event->modifiers() & Qt::ControlModifier) ? RM_BunnyFilled : RM_BunnyHollow;
      particleSystemReset();
      break;
  }

  update();

  return QOpenGLWidget::keyPressEvent(event);
}


void DisplayWidget::keyReleaseEvent(QKeyEvent *event)
{
  return QOpenGLWidget::keyReleaseEvent(event);
}


void DisplayWidget::mousePressEvent(QMouseEvent *event)
{
  if (event->button() == Qt::LeftButton)
  {
    if (event->modifiers() == Qt::ShiftModifier)
    {
      // rotate with camera when shift is pressed
      m_track_ball_camera.push(event->pos(), width(), height());
    }
    else if (event->modifiers() == Qt::ControlModifier)
    {
      m_track_ball.push(event->pos(), width(), height());
      m_track_ball_camera.push(event->pos(), width(), height());
    }
    else
      m_track_ball.push(event->pos(), width(), height());

    m_high_quality = false;
    event->accept();
    update();
  }
  else if (event->buttons() & Qt::MiddleButton)
  {
    m_pan_pos = event->pos();
  }
  else if (event->buttons() & Qt::RightButton)
  {
    m_z_dist_pos = event->y();
  }

  return QOpenGLWidget::mousePressEvent(event);
}


void DisplayWidget::mouseReleaseEvent(QMouseEvent *event)
{
  m_high_quality = true;

  // remember the old direction
  m_pan_dx_prev = m_pan_dx_prev + m_pan_dx;
  m_pan_dy_prev = m_pan_dy_prev + m_pan_dy;
  // and clear the new direction (necessary in case that paintEvent
  // would occur before mousePressEvent)
  m_pan_dx = 0.0f;
  m_pan_dy = 0.0f;

  m_z_dist_prev = m_z_dist_prev + m_z_dist;
  m_z_dist = 0.0f;

  update();

  return QOpenGLWidget::mouseReleaseEvent(event);
}


void DisplayWidget::mouseMoveEvent(QMouseEvent *event)
{
  if (event->buttons() & Qt::LeftButton)
  {
    if (event->modifiers() == Qt::ShiftModifier)
    {
      // rotate with camera when shift is pressed
      m_track_ball_camera.move(event->pos(), width(), height());
    }
    else if (event->modifiers() == Qt::ControlModifier)
    {
      m_track_ball.move(event->pos(), width(), height());
      m_track_ball_camera.move(event->pos(), width(), height());
    }
    else
      m_track_ball.move(event->pos(), width(), height());

    m_high_quality = false;
    event->accept();
    update();
  }
  else if (event->buttons() & Qt::MiddleButton)
  {
    float aspect = float(width()) / float(height());
    m_pan_dx = 20.0f * aspect * (float(event->x() - m_pan_pos.x()) / float(width()));
    // -event->y() because of different coordinate systems used in OpenGL and in Qt
    m_pan_dy = 20.0f * (float(m_pan_pos.y() - event->y()) / float(height()));

    qDebug() << "m_pan_dx=" << m_pan_dx
             << "m_pan_dy=" << m_pan_dy
             << "m_pan_prev_dx=" << m_pan_dx_prev
             << "m_pan_prev_dy=" << m_pan_dy_prev;

    event->accept();
    update();
  }
  else if (event->buttons() & Qt::RightButton)
  {
    m_z_dist = (event->y() - m_z_dist_pos) * 0.2f;
    qDebug() << "m_z_dist=" << m_z_dist << "m_z_dist_prev=" << m_z_dist_prev;
    event->accept();
    update();
  }

  return QOpenGLWidget::mouseMoveEvent(event);
}


void DisplayWidget::wheelEvent(QWheelEvent *event)
{
  float delta = float(event->delta()) * 0.0001f;

  m_scale += delta;

  if (m_scale >= 2.0f) m_scale = 2.0f;
  else if (m_scale <= 0.000001f) m_scale = 0.000001f;

  qDebug() << "scale=" << m_scale;

  update();

  return QOpenGLWidget::wheelEvent(event);
}


void DisplayWidget::timerEvent(QTimerEvent * /* event */)
{
  update();
}


void DisplayWidget::drawHUD(QPainter & painter)
{
  int text_y_pos = 20;

  painter.setRenderHint(QPainter::TextAntialiasing);
  //painter.setPen(Qt::white);
  painter.setPen(QColor(0x80, 0x80, 0x80));
  painter.setFont(QFont("Arial", SMALL_FONT_SIZE));

  // draw basic commands
  painter.drawText(10, text_y_pos, tr("Press 'H' to display help"));
  text_y_pos += (SMALL_FONT_SIZE + LINE_SPACING_SIZE);

  painter.drawText(10, text_y_pos, tr("Press 'I' to display status information"));
  text_y_pos += (SMALL_FONT_SIZE + LINE_SPACING_SIZE);

  painter.setFont(QFont("Arial", MEDIUM_FONT_SIZE));

  // draw info text
  text_y_pos = displayInfo(painter, text_y_pos, MEDIUM_FONT_SIZE);

  // draw help text
  displayHelp(painter, text_y_pos, MEDIUM_FONT_SIZE);
}


int DisplayWidget::displayInfo(QPainter & painter, int height, int font_size)
{
  static constexpr int offset_x = 10;

  if (!m_display_info) return height;

  height += (font_size + LINE_SPACING_SIZE);  // add an empty line

  painter.drawText(offset_x, height, currentParticleSystemName());
  height += (font_size + LINE_SPACING_SIZE);

  painter.drawText(offset_x, height, tr("Number of particles: %1").arg(m_sph_params->particleCount()));
  height += (font_size + LINE_SPACING_SIZE);

  painter.drawText(offset_x, height, tr("Grid Size: %1x%2x%3").arg(m_sph_params->gridWidth())
                                                              .arg(m_sph_params->gridHeight())
                                                              .arg(m_sph_params->gridDepth()));
  height += (font_size + LINE_SPACING_SIZE);

  painter.drawText(offset_x, height, currentRendererName());
  height += (font_size + LINE_SPACING_SIZE);

  painter.drawText(offset_x, height, tr("Frame: %1").arg(m_sim_frame_num));
  height += (font_size + LINE_SPACING_SIZE);

  painter.drawText(offset_x, height, tr("FPS: %1").arg(1000.0f / m_total_time_ms));
  height += (font_size + LINE_SPACING_SIZE);

  painter.drawText(offset_x, height, tr("Frame time (ms): %1").arg(m_total_time_ms));
  height += (font_size + LINE_SPACING_SIZE);

  painter.drawText(offset_x, height, tr("Simulation FPS: %1").arg(1000.0f / m_sim_time_ms));
  height += (font_size + LINE_SPACING_SIZE);

  painter.drawText(offset_x, height, tr("Simulation time (ms): %1").arg(m_sim_time_ms));
  height += (font_size + LINE_SPACING_SIZE);

  painter.drawText(offset_x, height, tr("Rendering FPS: %1").arg(1000.0f / m_ren_time_ms));
  height += (font_size + LINE_SPACING_SIZE);

  painter.drawText(offset_x, height, tr("Rendering time (ms): %1").arg(m_ren_time_ms));
  height += (font_size + LINE_SPACING_SIZE);

  return height;
}


int DisplayWidget::displayHelp(QPainter & painter, int height, int font_size)
{
  if (!m_display_help) return height;

  static QString help_strings[] = {
    tr("Press C to display the steps of Curvature Flow Renderer"),
    tr("Press P to display Curvature Flow renderering with refractions"),
    //tr("Press D to toggle drain effect On/Off"),
    tr("Press F to toggle fountain effect On/Off"),
    tr("Press W to emit wave"),
    tr("Press H to toggle On/Off this help message"),
    tr("Press I to toggle On/Off status information display"),
    tr("Press B to show/hide bounding volume box"),
    tr("Press R to restart simulation"),
    tr("Press O to reset all rotations, scale and translations"),
    tr("Press T to toggle simulation engine"),
    tr("Press V to toggle visualization engine"),
    tr("Press N to advance simulation by one frame (has to be paused first)"),
    tr("Press Z to freeze performance counters"),
    tr("Press SPACE BAR to pause/restart simulation"),
    tr("Press 1 to restart without using any voxel model"),
    tr("Press 2 to restart using armadillo voxel model"),
    tr("Press 3 to restart using bunny voxel model"),
    tr("Press 7 to use forest sky box"),
    tr("Press 8 to use beach with palm trees sky box"),
    tr("Press 9 to use clouds sky box")
  };

  static const int n = sizeof(help_strings) / sizeof(help_strings[0]);

  static constexpr int offset_x = 10;

  height += (font_size + LINE_SPACING_SIZE);  // add an empty line

  for (unsigned int i = 0; i < n; ++i)
  {
    painter.drawText(offset_x, height, help_strings[i]);
    height += (font_size + LINE_SPACING_SIZE);
  }

  return height;
}
