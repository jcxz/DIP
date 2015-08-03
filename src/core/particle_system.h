#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include "utils/ocl.h"
#include "utils/misc.h"
#include "core/sph_params.h"

#include <boost/compute/context.hpp>
#include <boost/compute/command_queue.hpp>
#include <stdexcept>
#include <array>


namespace io { class VoxelMesh; }


/**
 * Common base class for all particle systems
 */
class ParticleSystem
{
  public:
    ParticleSystem(const boost::compute::context & ctx, SPHParams *params)
      : m_cl_ctx(ctx)
      , m_queue()
      , m_sim_params(params)
      , m_particle_pos_buf(ctx)
      , m_particle_col_buf(ctx)
      , m_stats()
      , m_cl_timer()
      , m_print_stats(true)
    {
      assert(m_sim_params != nullptr);
      m_queue = boost::compute::command_queue(m_cl_ctx, m_cl_ctx.get_device(),
                                              boost::compute::command_queue::enable_profiling);
      m_cl_timer.setCommandQueue(m_queue);
    }

    virtual ~ParticleSystem(void)
    {
      if (m_print_stats)
      {
        std::cerr << "Simulation performance statistics:\n" <<  m_stats << std::endl;
      }
    }

    const boost::compute::context & clContext(void) const { return m_cl_ctx; }
    const boost::compute::command_queue & clCommandQueue(void) const { return m_queue; }

    GLuint positionsVBO(void) const { return m_particle_pos_buf.getGLID(); }
    GLuint colorsVBO(void) const { return m_particle_col_buf.getGLID(); }

    bool isFrameTimeValid(void) const { return m_cl_timer.isTimeValid(); }
    double frameTime(void) const { return m_cl_timer.elapsedMiliseconds(); }
    double waitForFrameTime(void) { m_cl_timer.waitForValid(); return frameTime(); }

    void setPrintStatsOnExit(bool enabled = true) { m_print_stats = enabled; }
    void clearPerformanceCounters(void) { m_stats.clear(); }
    const utils::ocl::PerfStats & performanceCounters(void) const { return m_stats; }

    void setSimulationParams(SPHParams *params) { m_sim_params = params; }

    // reset the particle system
    // initializes buffers and shared data
    bool reset(void) { return reset_impl(nullptr); }
    bool reset(const io::VoxelMesh & mesh) { return reset_impl(&mesh); }
    bool reset(const io::VoxelMesh *mesh) { return reset_impl(mesh); }
    bool reset(SPHParams *params) { m_sim_params = params; return reset_impl(nullptr); }

    // recalculate the particle system
    //void update(float time_step = 1.0f) { return update_impl(time_step); }
    void update(float time_step = 0.003f) { return update_impl(time_step); }

    //void updateBlocking(float time_step = 1.0f)
    void updateBlocking(float time_step = 0.003f)
    {
      update(time_step);
      m_queue.finish();
    }

  protected:
    virtual bool reset_impl(const io::VoxelMesh *mesh) = 0;
    virtual void update_impl(float time_step) = 0;

  protected:
    // kernel execution helpers
    bool runKernelHelper(cl_kernel kern, const char *name,
                         const size_t *global_work_size,
                         const size_t *local_work_size);

    bool runKernelHelper(const size_t *global_work_size, const size_t *local_work_size,
                         const utils::ocl::KernelArgs & kern_args)
    {
      return kern_args &&
          runKernelHelper(kern_args.kernel(), kern_args.kernelName(),
                          global_work_size, local_work_size);
    }

    bool runKernel(size_t global_work_size,
                   const utils::ocl::KernelArgs & kern_args)
    { return runKernelHelper(&global_work_size, nullptr, kern_args); }

    bool runKernel(size_t global_work_size, size_t local_work_size,
                   const utils::ocl::KernelArgs & kern_args)
    { return runKernelHelper(&global_work_size, &local_work_size, kern_args); }

    bool runKernel(boost::compute::kernel & kern, const char *name,
                   size_t global_work_size)
    { return runKernelHelper(kern, name, &global_work_size, nullptr); }

    bool runKernel(boost::compute::kernel & kern, const char *name,
                   size_t global_work_size, size_t local_work_size)
    { return runKernelHelper(kern, name, &global_work_size, &local_work_size); }

    // OpenGL synchronization
    template <const int N>
    using tGLBufferArray = typename std::array<utils::misc::ConstReference<utils::ocl::GLBuffer>, N>;

    template <const int N>
    //bool acquireGL(const std::array<const utils::ocl::GLBuffer &, N> & bufs)
    //bool acquireGL(const std::array<const utils::ocl::GLBuffer *, N> & bufs)
    bool acquireGL(const tGLBufferArray<N> & bufs)
    {
      // create a list of memory objects to be waited on
      cl_mem mems[N];
      //for (int i = 0; i < N; ++i) mems[i] = bufs[i].getCLID();
      for (int i = 0; i < N; ++i) mems[i] = bufs[i].get().getCLID();
      //for (int i = 0; i < N; ++i) mems[i] = bufs[i]->getCLID();

      // wait for OpenGL to finish rendering
      glFinish();

      // acquire access to the shared vertex buffer object
      cl_int err = clEnqueueAcquireGLObjects(m_queue, N, mems, 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        WARNM("Failed to acquire an exclusive access to one of OpenGL's vertex buffer objects: "
              << boost::compute::opencl_error::to_string(err)
              << "(" << err << ")");
        return false;
      }

      return true;
    }

    template <const int N>
    //bool releaseGL(const std::array<const utils::ocl::GLBuffer &, N> & bufs)
    bool releaseGL(const tGLBufferArray<N> & bufs)
    {
      // create a list of memory objects to be waited on
      cl_mem mems[N];
      //for (int i = 0; i < N; ++i) mems[i] = bufs[i].getCLID();
      for (int i = 0; i < N; ++i) mems[i] = bufs[i].get().getCLID();

      // unlock the vertex buffer object, so that OpenGL can continue using it
      cl_int err = clEnqueueReleaseGLObjects(m_queue, N, mems, 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        WARNM("Failed to release an exclusive access to one of OpenGL's vertex buffer objects: "
              << boost::compute::opencl_error::to_string(err)
              << "(" << err << ")");
        return false;
      }

      // wait for OpenCL to finish processing
      clFinish(m_queue);

      return true;
    }

  protected:
    // OpenCL context data
    boost::compute::context m_cl_ctx;       // OpenCL context
    boost::compute::command_queue m_queue;  // OpenCL command queue

    // simulation parameters
    SPHParams *m_sim_params;

    // memory objects with particle data
    utils::ocl::GLBuffer m_particle_pos_buf;  // a buffer with particle positions (shared with OpenGL)
    utils::ocl::GLBuffer m_particle_col_buf;  // a buffer with particle colors (shared with OpenGL)

    // statistics and performance counters
    utils::ocl::PerfStats m_stats;
    utils::ocl::Timer m_cl_timer;

  private:
    bool m_print_stats;           // whether printing statistics in destructor is enabled
};

#endif // PARTICLE_SYSTEM_H
