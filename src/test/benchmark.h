#ifndef BENCHMARK_H
#define BENCHMARK_H

#include "core/sph_params.h"
#include "core/sph_naive.h"
#include "core/sph_optimized.h"
#include "core/sph_uniform_grid.h"
#include "core/instancing_renderer.h"
#include "core/point_sprite_renderer.h"
#include "core/curvature_flow_renderer.h"
#include "core/ray_cast_renderer.h"

#include <vector>
#include <stdexcept>
#include <memory>
#include <QString>


namespace test {

class BenchMark
{
  public:
    enum RendererType {
      REN_Instancing,
      REN_PointSprite,
      REN_CurvatureFlow,
      REN_RayCast,
      REN_None
    };

    enum ParticleSystemType {
      PS_SPHUniformGrid,
      PS_SPHOptimized,
      PS_SPHNaive,
      PS_None
    };

  public:
    typedef std::vector<double> tFrameTimesVec;
    typedef utils::stats::Statistics<double> tStats;

  private:
    // benchmark defaults
    static constexpr int DEFAULT_FRAME_COUNT = 1000;
    static constexpr int DEFAULT_SAVE_INTENSITY = 0;

    // simulation defaults
    static constexpr int DEFAULT_PARTICLE_COUNT = 60000;
    static constexpr int DEFAULT_GRID_WIDTH  = 64;
    static constexpr int DEFAULT_GRID_HEIGHT = 64;
    static constexpr int DEFAULT_GRID_DEPTH  = 64;

    // rendering defaults
    static constexpr int DEFAULT_VIEWPORT_W = 1920;
    static constexpr int DEFAULT_VIEWPORT_H = 1080;

  public:
    BenchMark(void)
      : m_cl_ctx()
      , m_sph_params(nullptr)
      , m_sph_naive(nullptr)
      , m_sph_optimized(nullptr)
      , m_sph_uniform_grid(nullptr)
      , m_cur_ps_type(PS_SPHUniformGrid)
      , m_ren_instancing(nullptr)
      , m_ren_point_sprite(nullptr)
      , m_ren_curvature_flow(nullptr)
      , m_ren_ray_cast(nullptr)
      , m_cur_ren_type(REN_CurvatureFlow)
      , m_sim_frame_times()
      , m_ren_frame_times()
      , m_total_frame_times()
      , m_frame_cnt(DEFAULT_FRAME_COUNT)
      , m_save_intensity(DEFAULT_SAVE_INTENSITY)
      , m_viewport_w(DEFAULT_VIEWPORT_W)
      , m_viewport_h(DEFAULT_VIEWPORT_H)
      , m_sim_stats()
      , m_ren_stats()
      , m_total_stats()
      , m_sim_median(0.0)
      , m_ren_median(0.0)
      , m_total_median(0.0)
    {
      if (!init())
        throw std::runtime_error("Failed to initialize OpenCL/OpenGL context in Benchmark");
    }

    const tFrameTimesVec & simulationTimesRef(void) const { return m_sim_frame_times; }
    const tFrameTimesVec & renderingTimesRef(void) const { return m_ren_frame_times; }
    const tFrameTimesVec & totalTimesRef(void) const { return m_total_frame_times; }

    const tStats & simStats(void) const { return m_sim_stats; }
    double simMedian(void) const { return m_sim_median; }

    const tStats & renStats(void) const { return m_ren_stats; }
    double renMedian(void) const { return m_ren_median; }

    const tStats & totalStats(void) const { return m_total_stats; }
    double totalMedian(void) const { return m_total_median; }

    // This property determines how often will the output of rendering be saved.
    // A number that is less than or equal to 0 means the output will never be
    // saved, number 1 means that every frame will be saved, number 50 means that
    // only every 50th frame will be saved to a file on disk, ...
    int saveIntensity(void) const { return m_save_intensity; }
    void setSaveIntensity(int intensity) { m_save_intensity = intensity; }

    int frameCount(void) const { return m_frame_cnt; }
    void setFrameCount(const int frame_cnt) { m_frame_cnt = frame_cnt; }

    void setRenderer(RendererType type) { m_cur_ren_type = type; }
    void setParticleSystem(ParticleSystemType type) { m_cur_ps_type = type; }

    const SPHParams & sphParamsRef(void) const
    { assert(m_sph_params.get() != nullptr); return *m_sph_params; }

    SPHParams & sphParamsRef(void)
    { assert(m_sph_params.get() != nullptr); return *m_sph_params; }

    void setCurvatureFlowIterations(const int n)
    {
      assert(m_ren_curvature_flow.get() != nullptr);
      return m_ren_curvature_flow->setSmoothingIterationCount(n);
    }

    // Functions to run start benchmarking a given simulation type.
    // When dir_name is not null, then this string will be used as a
    // directory where inidividual output frames will be saved as images.
    bool run(const QString & dir_name = QString());

    void printStats(std::ostream & os);
    bool saveStats(const QString & filename);
    bool saveHistogramAsCsv(const QString & filename);
    bool saveTimelineAsCsv(const QString & filename);

    bool saveLastSimulationPerformanceCounters(const QString & filename);
    bool saveLastRenderingPerformanceCounters(const QString & filename);

  private:
    bool init(void);
    ParticleSystem *currentParticleSystem(void) const;
    BaseRenderer *currentRenderer(void) const;

  private:
    boost::compute::context m_cl_ctx;

    // fluid simulators
    std::unique_ptr<SPHParams> m_sph_params;
    std::unique_ptr<SPHNaive> m_sph_naive;
    std::unique_ptr<SPHOptimized> m_sph_optimized;
    std::unique_ptr<SPHUniformGrid> m_sph_uniform_grid;
    ParticleSystemType m_cur_ps_type;

    // fluid renderers
    std::unique_ptr<InstancingRenderer> m_ren_instancing;
    std::unique_ptr<PointSpriteRenderer> m_ren_point_sprite;
    std::unique_ptr<CurvatureFlowRenderer> m_ren_curvature_flow;
    std::unique_ptr<RayCastRenderer> m_ren_ray_cast;
    RendererType m_cur_ren_type;

    // benchmark settings
    tFrameTimesVec m_sim_frame_times;    // a vector with meassured frame times for the simulation part
    tFrameTimesVec m_ren_frame_times;    // a vector with meassured frame times for the rendering part
    tFrameTimesVec m_total_frame_times;  // a vector with meassured frame times for the rendering part
    int m_frame_cnt;
    int m_save_intensity;
    int m_viewport_w;
    int m_viewport_h;

    // a cache for summary statistics
    tStats m_sim_stats;
    tStats m_ren_stats;
    tStats m_total_stats;
    double m_sim_median;
    double m_ren_median;
    double m_total_median;
};


bool runBenchmarks(const size_t particle_counts[], const int particle_counts_n,
                   const size_t grid_sizes[][3],   const int grid_sizes_n,
                   const QString & dir_name =  QString("."));

} // End of namespace test

#endif // BENCHMARK_H
