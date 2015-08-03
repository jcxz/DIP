/*
   TODO:
     * turn on and off some effects every other frame in the runParticleSystem
       method to vary the load a bit.
 */

#include "test/benchmark.h"
#include "utils/debug.h"
#include "utils/stats.h"
#include "utils/timer.h"

#include <chrono>
#include <unordered_map>

#include <QDir>
#include <QFile>
#include <QDateTime>
#include <QOpenGLTimerQuery>



namespace {

double calcMedian(const std::vector<double> & frame_times)
{
  // copy the data for sorting, to do not destroy the original order of timings
  std::vector<double> tmp(frame_times);
  std::sort(tmp.begin(), tmp.end());
  return tmp[tmp.size() / 2];
}


void printStatsHelper(std::ostream & os, const double median,
                      const test::BenchMark::tStats & stats)
{
  os << "  min frame time (ms)    : "  << (stats.min())          << std::endl;
  os << "  max frame time (ms)    : "  << (stats.max())          << std::endl;
  os << "  avg frame time (ms)    : "  << (stats.avg())          << std::endl;
  os << "  median frame time (ms) : "  << (median)               << std::endl;
  os << "  max fps                : "  << (1000.0 / stats.min()) << std::endl;
  os << "  min fps                : "  << (1000.0 / stats.max()) << std::endl;
  os << "  avg fps                : "  << (1000.0 / stats.avg()) << std::endl;
  os << "  median fps             : "  << (1000.0 / median)      << std::endl;
}

void printStatsAsCsvHelper(QTextStream & s, const double median,
                           const test::BenchMark::tStats & stats)
{
  s << (stats.min())          << ";"
    << (stats.max())          << ";"
    << (stats.avg())          << ";"
    << (median)               << ";"
    << (1000.0 / stats.max()) << ";"
    << (1000.0 / stats.min()) << ";"
    << (1000.0 / stats.avg()) << ";"
    << (1000.0 / median)      << ";";
}

} // End of private namespace



namespace test {

bool BenchMark::run(const QString & dir_name)
{
  ParticleSystem *ps = currentParticleSystem();
  BaseRenderer *ren = currentRenderer();

  // initialize particle simulator
  if (!ps->reset())
  {
    ERRORM("runParticleSystem: Failed to reset fluid particle system");
    return false;
  }

  ps->clearPerformanceCounters();

  // initialize renderer
  m_ren_ray_cast->setGrid(m_sph_uniform_grid->cellStartsBuffer(),
                          m_sph_uniform_grid->cellEndsBuffer());

  m_ren_ray_cast->setGridSize(m_sph_params->gridWidth(),
                              m_sph_params->gridHeight(),
                              m_sph_params->gridDepth());

  if (!ren->reset(m_viewport_w, m_viewport_h))
  {
    ERRORM("Failed to initialize renderer");
    return false;
  }

  if (!ren->resize(m_viewport_w, m_viewport_h))
  {
    ERRORM("Failed to resize renderer");
    return false;
  }

  ren->setUseUniformColor(true);
  ren->setUseLighting(true);
  ren->setLightPosition(QVector3D(300.0f, 300.0f, 300.0f));
  ren->setLightAmbientColor(QVector3D(0.2f, 0.2f, 0.2f));
  ren->setLightDiffuseColor(QVector3D(0.8f, 0.8f, 0.8f));
  ren->setClearColor(QVector4D(0.8f, 0.8f, 0.8f, 1.0f));
  ren->setDrawSkyBox(true);
  ren->setMixing(false);

  ren->setSkyBoxTextures(":/data/skybox/forest_posx.jpg",
                         ":/data/skybox/forest_negx.jpg",
                         ":/data/skybox/forest_posy.jpg",
                         ":/data/skybox/forest_negy.jpg",
                         ":/data/skybox/forest_posz.jpg",
                         ":/data/skybox/forest_negz.jpg");

  ren->clearPerformanceCounters();

  // create OpenGL timer query
  QOpenGLTimeMonitor time_monitor;
  time_monitor.setSampleCount(2);
  if (!time_monitor.create())
  {
    ERRORM("Failed to create OpenGL time monitor");
    return false;
  }

  // create frame buffer for offscreen rendering
  QOpenGLFramebufferObject fbo(m_viewport_w, m_viewport_h);
  if (!fbo.isValid())
  {
    ERRORM("Failed to create frame buffer, frame buffer is ot valid");
    return false;
  }

  fbo.bind();

  ren->setDefaultFBO(fbo.handle());

  // clear the statistics counters
  m_sim_frame_times.clear();
  m_ren_frame_times.clear();
  m_total_frame_times.clear();

  m_sim_stats.clear();
  m_ren_stats.clear();
  m_total_stats.clear();

  // run the simulation
  constexpr float scale = 0.03f; //0.1f;

  for (int i = 0; i < m_frame_cnt; ++i)
  {
#if 0
    auto start = std::chrono::high_resolution_clock::now();
    sph.updateBlocking();
    auto end = std::chrono::high_resolution_clock::now();
    double t = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
    stats.add(t);
#elif 0
    t.start();
    sph.updateBlocking();
    t.stop();
    stats.add(utils::time::elapsedMilliSeconds(t));
#else
    ps->update();
#endif

    time_monitor.recordSample();

    if (ren == m_ren_ray_cast.get()) m_ren_ray_cast->calcNormals();

    cl_float3 bvs = m_sph_params->boundingVolumeSize();
    ren->render(QQuaternion(),                            // rotation
                QVector3D(scale, scale, scale),           // scale
                //QVector3D(0.0f, 0.0f, -0.5f * bvs.s[2]),  // translation
                QVector3D(0.0f, 0.0f, -6.0f * scale * bvs.s[2]),  // translation
                QQuaternion(),                            // camera rotation
                ps->positionsVBO(),
                ps->colorsVBO(),
                m_sph_params->particleCount(),
                QVector3D(bvs.s[0], bvs.s[1], bvs.s[2]));

    time_monitor.recordSample();

    // record the simulation time (in milliseconds)
    double t_sim = ps->waitForFrameTime();

    // record the rendering time (in milliseconds)
    double t_ren = time_monitor.waitForIntervals()[0] / 1000000.0;

    time_monitor.reset();

    // record statistics
    m_sim_frame_times.push_back(t_sim);
    m_ren_frame_times.push_back(t_ren);
    m_total_frame_times.push_back(t_sim + t_ren);

    m_sim_stats.add(t_sim);
    m_ren_stats.add(t_ren);
    m_total_stats.add(t_sim + t_ren);

    if ((!dir_name.isNull()) && (m_save_intensity > 0) && ((i % m_save_intensity) == 0))
    {
      QImage frame(fbo.toImage());
      if (frame.isNull())
      {
        WARNM("Failed to grab frame #" << i);
      }

      QString fname(QString("%1/frame_%2.jpg").arg(dir_name).arg(i));
      if (!frame.save(fname))
      {
        WARNM("Failed to save frame #" << i << " to file " << fname.toStdString());
      }
    }
  }

  fbo.release();

  // recalculate median
  m_sim_median   = calcMedian(m_sim_frame_times);
  m_ren_median   = calcMedian(m_ren_frame_times);
  m_total_median = calcMedian(m_total_frame_times);

  return true;
}


void BenchMark::printStats(std::ostream & os)
{
  os << "Statistics:" << std::endl;
  os << "  Frames executed        : "  << (m_frame_cnt)                   << std::endl;
  os << "  Number of particles    : "  << (m_sph_params->particleCount()) << std::endl;
  os << "  Grid size              : [" << (m_sph_params->gridWidth())     << "x"
                                       << (m_sph_params->gridHeight())    << "x"
                                       << (m_sph_params->gridDepth())     << "]" << std::endl;

  os << "Simulation:" << std::endl;
  printStatsHelper(os, m_sim_median, m_sim_stats);

  os << "Rendering:" << std::endl;
  printStatsHelper(os, m_ren_median, m_ren_stats);
}


bool BenchMark::saveStats(const QString & filename)
{
  std::ofstream file(filename.toStdString());
  if (!file.is_open())
  {
    ERRORM("Failed to create file: " << filename.toStdString());
    return false;
  }

  printStats(file);

  return true;
}


bool BenchMark::saveHistogramAsCsv(const QString & filename)
{
  // create the output file
  QFile f(filename);
  if (!f.open(QFile::WriteOnly))
  {
    ERRORM("Failed to create file named: " << filename.toStdString());
    return false;
  }

  // build the histogram
  std::unordered_map<double, int> hist;
  for (const double time : m_sim_frame_times)
  {
    double tv = std::trunc(time * 1000.0f) / 1000.0f;  // round the time value to 3 decimal places
    ++hist[tv];
  }

  // write the histogram to file
  QTextStream s(&f);
  s << "frameTime;count\n";
  for (const auto & it : hist)
  {
    s << it.first << ";" << it.second << "\n";
  }

  return true;
}


bool BenchMark::saveTimelineAsCsv(const QString & filename)
{
  // create the output file
  QFile f(filename);
  if (!f.open(QFile::WriteOnly))
  {
    ERRORM("Failed to create file named: " << filename.toStdString());
    return false;
  }

  // write the timeline to file
  QTextStream s(&f);
  s << "frameNumber;simulationTime;renderingTime;totalFrameTime\n";
  for (int i = 0; i < m_frame_cnt; ++i)
  {
    double t_sim = m_sim_frame_times[i];
    double t_ren = m_ren_frame_times[i];
    double t_total = t_sim + t_ren;
    s << i << ";" << t_sim << ";" << t_ren << ";" << t_total << "\n";
  }

  return true;
}


bool BenchMark::saveLastSimulationPerformanceCounters(const QString & filename)
{
  std::ofstream file(filename.toStdString());
  if (!file.is_open())
  {
    ERRORM("Failed to create file: " << filename.toStdString());
    return false;
  }

  currentParticleSystem()->performanceCounters().printAsTable(file);

  return true;
}


bool BenchMark::saveLastRenderingPerformanceCounters(const QString & filename)
{
  std::ofstream file(filename.toStdString());
  if (!file.is_open())
  {
    ERRORM("Failed to create file: " << filename.toStdString());
    return false;
  }

  currentRenderer()->performanceCounters().printAsTable(file);

  return true;
}


bool BenchMark::init(void)
{
  // create OpenCL context
  if (!utils::ocl::initCLGLContext(m_cl_ctx))
  {
    ERRORM("runParticleSystem: Failed to create OpenGL/OpenCL shared context");
    return false;
  }

  // initialize simulation parameters
  m_sph_params.reset(new SPHParams(m_cl_ctx));
  m_sph_params->setParticleCount(DEFAULT_PARTICLE_COUNT);
  m_sph_params->setGridSize(DEFAULT_GRID_WIDTH,
                            DEFAULT_GRID_HEIGHT,
                            DEFAULT_GRID_DEPTH);

  // create simulators
  m_sph_naive.reset(new SPHNaive(m_cl_ctx, m_sph_params.get()));
  m_sph_naive->setPrintStatsOnExit(false);

  m_sph_optimized.reset(new SPHOptimized(m_cl_ctx, m_sph_params.get()));
  m_sph_optimized->setPrintStatsOnExit(false);

  m_sph_uniform_grid.reset(new SPHUniformGrid(m_cl_ctx, m_sph_params.get()));
  m_sph_uniform_grid->setPrintStatsOnExit(false);

  // create renderers
  m_ren_instancing.reset(new InstancingRenderer);
  m_ren_instancing->setPrintStatsOnExit(false);

  m_ren_point_sprite.reset(new PointSpriteRenderer);
  m_ren_point_sprite->setPrintStatsOnExit(false);

  m_ren_curvature_flow.reset(new CurvatureFlowRenderer);
  m_ren_curvature_flow->setPrintStatsOnExit(false);

  m_ren_ray_cast.reset(new RayCastRenderer(m_cl_ctx, m_sph_uniform_grid->clCommandQueue()));
  m_ren_ray_cast->setPrintStatsOnExit(false);

  return true;
}


ParticleSystem *BenchMark::currentParticleSystem(void) const
{
  switch (m_cur_ps_type)
  {
    case PS_SPHUniformGrid: return m_sph_uniform_grid.get();
    case PS_SPHOptimized:   return m_sph_optimized.get();
    case PS_SPHNaive:       return m_sph_naive.get();
    case PS_None:           return nullptr;
  }

  return nullptr;
}


BaseRenderer *BenchMark::currentRenderer(void) const
{
  switch (m_cur_ren_type)
  {
    case REN_Instancing:    return m_ren_instancing.get();
    case REN_PointSprite:   return m_ren_point_sprite.get();
    case REN_CurvatureFlow: return m_ren_curvature_flow.get();
    case REN_RayCast:       return m_ren_ray_cast.get();
    case REN_None:          return nullptr;
  }

  return nullptr;
}


bool runBenchmarks(const size_t particle_counts[], const int particle_counts_n,
                   const size_t grid_sizes[][3],   const int grid_sizes_n,
                   const QString & dir_name)
{
  QString path(utils::misc::makePathFromCurrentDate(dir_name));
  QFile f(path + "/all_stats.csv");
  if (!f.open(QFile::WriteOnly))
  {
    ERRORM("Failed to create file named: " << QString(path + "/all_stats.csv").toStdString());
    return false;
  }

  QTextStream s(&f);
  s << "particleCount;gridWidth;gridHeight;gridDepth;"
       "simMinTime;simMaxTime;simAvgTime;simMedianTime;simMinFps;simMaxFps;simAvgFps;simMedianFps;"
       "renMinTime;renMaxTime;renAvgTime;renMedianTime;renMinFps;renMaxFps;renAvgFps;renMedianFps;"
       "totalMinTime;totalMaxTime;totalAvgTime;totalMedianTime;totalMinFps;totalMaxFps;totalAvgFps;totalMedianFps;\n";

  // initialize benchmark
  test::BenchMark b;
  //b.setSaveIntensity(100);

  // run benchmark
  for (int j = 0; j < particle_counts_n; ++j)
  {
    size_t part_cnt = particle_counts[j];
    QString filename("%1/stats_%2.csv");

    QFile f2(filename.arg(path).arg(part_cnt));
    if (!f2.open(QFile::WriteOnly))
    {
      ERRORM("Failed to create file named: " << filename.arg(path).arg(part_cnt).toStdString());
      return false;
    }

    QTextStream s2(&f2);
    s2 << "Počet častíc;Rozmery mriežky;Simulácia;Renderovanie;\n";

    for (int i = 0; i < grid_sizes_n; ++i)
    {
      size_t grid_w   = grid_sizes[i][0];
      size_t grid_h   = grid_sizes[i][1];
      size_t grid_d   = grid_sizes[i][2];

      std::cerr << "=====================================================================";
      std::cerr << "Running " << part_cnt << " " << grid_w << "x" << grid_h << "x" << grid_d;

      // set parameters
      SPHParams & params = b.sphParamsRef();
      params.setParticleCount(part_cnt);
      params.setGridSize(grid_w, grid_h, grid_d);

      // run simulation
      b.run(path);

      // write the results
      const test::BenchMark::tStats & sim_stats = b.simStats();
      const test::BenchMark::tStats & ren_stats = b.renStats();

      s << part_cnt << ";" << grid_w << ";" << grid_h << ";" << grid_d << ";";

      printStatsAsCsvHelper(s, b.simMedian(),   sim_stats);
      printStatsAsCsvHelper(s, b.renMedian(),   ren_stats);
      printStatsAsCsvHelper(s, b.totalMedian(), b.totalStats());

      s << "\n";

      s2 << part_cnt << ";" << grid_w << "x" << grid_h << "x" << grid_d << ";"
         << sim_stats.avg() << ";" << ren_stats.avg() << ";\n";

      // save results for inidividual runs
      b.saveStats(QString("%1/stats_%2_%3x%4x%5.txt")
                  .arg(path).arg(part_cnt)
                  .arg(grid_w).arg(grid_h).arg(grid_d));

      b.saveHistogramAsCsv(QString("%1/histogram_%2_%3x%4x%5.csv")
                           .arg(path).arg(part_cnt)
                           .arg(grid_w).arg(grid_h).arg(grid_d));

      b.saveTimelineAsCsv(QString("%1/timeline_%2_%3x%4x%5.csv")
                          .arg(path).arg(part_cnt)
                          .arg(grid_w).arg(grid_h).arg(grid_d));

      b.saveLastSimulationPerformanceCounters(
            QString("%1/perf_counters_sim_%2_%3x%4x%5.txt")
            .arg(path).arg(part_cnt)
            .arg(grid_w).arg(grid_h).arg(grid_d));

      b.saveLastRenderingPerformanceCounters(
            QString("%1/perf_counters_ren_%2_%3x%4x%5.txt")
            .arg(path).arg(part_cnt)
            .arg(grid_w).arg(grid_h).arg(grid_d));
    }
  }

  return true;
}

} // End of namespace test
