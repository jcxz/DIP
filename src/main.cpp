#include "gui/display_widget.h"
#include "test/benchmark.h"
#include "test/sph_test.h"
#include "utils/debug.h"
#include "utils/misc.h"

#include <QApplication>
#include <QSurfaceFormat>
#include <QOffscreenSurface>
#include <QMessageBox>



static bool initOffscreen(QOpenGLContext & ctx, QOffscreenSurface & surf)
{
  if (!ctx.create())
  {
    ERRORM("Failed to create OpenGL context");
    return false;
  }

  surf.setFormat(ctx.format());
  surf.create();

  ctx.makeCurrent(&surf);

  return true;
}


static int runRenderingBenchmark(void)
{
  QOpenGLContext ctx;
  QOffscreenSurface surf;
  if (!initOffscreen(ctx, surf)) return 1;

  // initialize benchmark
  test::BenchMark b;

  if (!b.sphParamsRef().parseCmdArgs(QApplication::arguments()))
  {
    ERRORM("Failed to parse command line arguments");
    return 1;
  }

  //b.setSaveIntensity(10);
  b.setSaveIntensity(100);

  // initialize the output path
  QString path(utils::misc::makePathFromCurrentDate("../benchmark_rendering"));

  // run the benchmark with varying number of curvatureflow iterations
  //for (int i = 10; i < 200; i += 10)
  int i = 50;
  {
    QString dir(path + "/" + QString::number(i));
    if (!QDir().mkpath(dir))
    {
      WARNM("Failed to create output folder " << dir.toStdString());
      //continue;
    }

    // run benchmark
    b.setCurvatureFlowIterations(i);
    b.run(dir);

    // save results
    b.saveStats(dir + "/stats.txt");
    b.saveHistogramAsCsv(dir + "/histogram.csv");
    b.saveTimelineAsCsv(dir + "/timeline.csv");
    b.saveLastSimulationPerformanceCounters(dir + "/perf_counters_sim.txt");
    b.saveLastRenderingPerformanceCounters(dir + "/perf_counters_ren.txt");
  }

  return 0;
}


static int runComplexBenchmark(void)
{
  QOpenGLContext ctx;
  QOffscreenSurface surf;
  if (!initOffscreen(ctx, surf)) return 1;

  // initialize benchmark
  static constexpr size_t particle_counts[] = {
#if 0
  #ifdef SMALL_BENCHMARK  // intended for notebook, or other low performance device
    16384, 32768, 65536, 131072, 262144
  #else
    131072, 262144, 524288, 1048576, 2097152, 4194304
  #endif
#else
    131072, 262144, 524288
#endif
  };

  static constexpr int particle_counts_n = sizeof(particle_counts) /
                                           sizeof(particle_counts[0]);

  static constexpr size_t grid_sizes[][3] = {
    {  32,  32,  32 },
    {  64,  64,  64 },
    { 128, 128, 128 },
    { 256, 256, 256 },
#if 0
    { 512, 512, 512 },

    {  64,  32,  64 },
    { 128,  32, 128 },
    { 128,  64, 128 },
    { 256,  32, 256 },
    { 256,  64, 256 },
    { 256, 128, 256 },
    { 512,  32, 512 },
    { 512,  64, 512 },
    { 512, 128, 512 },
    { 512, 256, 512 },
#endif
  };

  static constexpr int grid_sizes_n = sizeof(grid_sizes) /
                                      sizeof(grid_sizes[0]);

  test::runBenchmarks(particle_counts, particle_counts_n,
                      grid_sizes, grid_sizes_n,
                      "../benchmark_complete");

  return 0;
}


static int runVerification(void)
{
  using namespace test;

  QOpenGLContext ctx;
  QOffscreenSurface surface;
  if (!initOffscreen(ctx, surface)) return 1;

  //SPHOptimized_vs_SPHOptimized t(2000);
  //SPHUniformGrid_vs_SPHUniformGrid t(2000);
  SPHUniformGrid_vs_SPHOptimized t(2000);
  //SPHUniformGrid_vs_SPHNaive t(2000);
  //SPHOptimized_vs_SPHNaive t(2000);
  if (!t.run(std::cerr, 10, 1)) return 1;

  return 0;
}


static int runDebug(void)
{
  QOpenGLContext gl_ctx;
  QOffscreenSurface surface;
  if (!initOffscreen(gl_ctx, surface)) return 1;

  //boost::compute::device dev;
  //if (!utils::ocl::findGPUDevice(dev)) return 1;
  //boost::compute::context cl_ctx(dev);

  boost::compute::context cl_ctx;
  if (!utils::ocl::initCLGLContext(cl_ctx)) return 1;

  SPHParams params(cl_ctx);
  params.setSeed(10000);
  params.setParticleCount(100);
  params.setGridSize(2, 2, 2);

  SPHUniformGrid sph(cl_ctx, &params);
  sph.setPrintStatsOnExit(false);
  if (!sph.reset()) return 1;

  for (int i = 0; i < 10; ++i) sph.update();

  return 0;
}


static void displayErrorMessage(const QString & msg)
{
  QMessageBox::critical(nullptr, QObject::tr("Error"), msg);
}



int main(int argc, char *argv[])
{
  QApplication app(argc, argv);

  // run benchmark or verification instead of gui if requested
  QStringList args = QApplication::arguments();
  if (args.contains("-benchmark-rendering")) return runRenderingBenchmark();
  else if (args.contains("-benchmark")) return runComplexBenchmark();
  else if (args.contains("-verify")) return runVerification();
  else if (args.contains("-debug")) return runDebug();

  // Setup default OpenGL format (core profil version 3.3)
  QSurfaceFormat format;
  //format.setVersion(3, 3);
  format.setProfile(QSurfaceFormat::CoreProfile);
  //format.setProfile(QSurfaceFormat::CompatibilityProfile);
  format.setOption(QSurfaceFormat::DebugContext);
  QSurfaceFormat::setDefaultFormat(format);

  DisplayWidget w;

  // Note: this syntax is only supported by Qt 5, Qt 4 will not munch this
  QObject::connect(&w, &DisplayWidget::error, displayErrorMessage);

  w.show();

  return app.exec();
}
