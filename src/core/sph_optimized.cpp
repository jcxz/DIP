#include "core/particle_system_internal.h"
#include "core/sph_optimized.h"
#include "utils/macros.h"
#include "utils/debug.h"
#include "utils/ocl.h"

#include <QFile>
#include <ctime>
#include <boost/compute/algorithm/sort.hpp>
#include <boost/compute/algorithm/reduce.hpp>

//#define DEBUG_PARAMS
#define BLOCK_SIZE 256 //192 //128 //96 //32 //16 //512



bool SPHOptimized::init(void)
{
  try
  {
    // build opencl program
#ifdef DEBUG_SOURCE
    QFile f(":/src/opencl/sph_optimized.cl");
    if (!f.open(QFile::ReadOnly)) return false;

    m_sph_prog = boost::compute::program::create_with_source(f.readAll().toStdString(), m_cl_ctx);
    m_sph_prog.build("-I../DIP/src/core -DDEBUG_SOURCE");
#else
    QFile f_cl(":/src/opencl/sph_optimized.cl");
    QFile f_h(":/src/core/sph_ocl_common.h");
    if ((!f_cl.open(QFile::ReadOnly)) || (!f_h.open(QFile::ReadOnly)))
    {
      ERRORM("Failed to open OpenCL program source files");
      return false;
    }

    QByteArray src(f_h.readAll());
    src.append(f_cl.readAll());

    m_sph_prog = boost::compute::program::create_with_source(src.toStdString(), m_cl_ctx);
    m_sph_prog.build();
#endif

    // print build log in case there are any warnings
    std::string log(m_sph_prog.build_log());
    if (!log.empty())
    {
      WARNM("---------------------- Build log ---------------------------------");
      WARNM(log);
      WARNM("------------------- End of Build log -----------------------------");
    }

    // create kernels
    m_sph_reset_kernel = m_sph_prog.create_kernel("sph_optimized_reset");
    utils::ocl::printKernelInfo(m_sph_reset_kernel);

    m_sph_compute_step_kernel = m_sph_prog.create_kernel("sph_optimized_compute_step");
    utils::ocl::printKernelInfo(m_sph_compute_step_kernel);

    m_sph_compute_force_kernel = m_sph_prog.create_kernel("sph_optimized_compute_force");
    utils::ocl::printKernelInfo(m_sph_compute_force_kernel);

    m_sph_compute_pressure_kernel = m_sph_prog.create_kernel("sph_optimized_compute_pressure");
    utils::ocl::printKernelInfo(m_sph_compute_pressure_kernel);
  }
  catch (const boost::compute::opencl_error & e)
  {
    ERRORM("SPHOptimized: An OpenCL error occured: " << e.what());
    if (e.error_code() == CL_BUILD_PROGRAM_FAILURE) ERRORM(m_sph_prog.build_log());
    return false;
  }
  catch (const std::exception & e)
  {
    ERRORM("SPHOptimized: An unexpected error occured during SPHOptimized initialization: "
           << e.what());
    return false;
  }

  return true;
}


bool SPHOptimized::reset_impl(const io::VoxelMesh * /* mesh */)
{
  size_t part_cnt = m_sim_params->particleCount();

  // calculate optimal global and local work sizes and buffer padding
  size_t block_size = BLOCK_SIZE;
  size_t grid_size = (part_cnt + block_size - 1) / block_size;
  size_t local_work_size = block_size;
  size_t global_work_size = grid_size * block_size;
  size_t padding = global_work_size - part_cnt;

  INFOM("Global Work Size      : " << global_work_size);
  INFOM("Local Work Size       : " << local_work_size);
  INFOM("Memory Object padding : " << padding);

  /* initialize OpenGL shared buffer with particle positions */
  if (!m_particle_pos_buf.bufferData(nullptr, (part_cnt + padding) * sizeof(cl_float4),
                                     utils::ocl::GLBuffer::WRITE_ONLY))
  {
    ERRORM("SPHOptimized: Failed to initialize position GLBuffer");
    return false;
  }

  /* allocate buffers on GPU */
  try
  {
    m_velocity_buf = boost::compute::buffer(m_cl_ctx, (part_cnt + padding)* sizeof(cl_float4));
    m_prev_velocity_buf = boost::compute::buffer(m_cl_ctx, (part_cnt + padding)* sizeof(cl_float4));
    m_pressure_buf = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float));
    m_density_buf = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float));
    m_force_buf = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float4));
  }
  catch (const boost::compute::opencl_error & e)
  {
    ERRORM("SPHOptimized: An OpenCL error occured: " << e.what());
    return false;
  }
  catch (const std::exception & e)
  {
    ERRORM("SPHOptimized: An unexpected error occured during SPHOptimized initialization: "
           << e.what());
    return false;
  }

  /* set kernel parameters that do not change once the simulation is prepared */
  m_sim_params->setTime(0.0f);
  m_sim_params->clearEffectFlags();

  INFOM(*m_sim_params);

  /* compute pressure kernel's arguments */
  if (!utils::ocl::KernelArgs(m_sph_compute_pressure_kernel, "sph_optimized_compute_pressure")
            .arg(m_particle_pos_buf.getCLID())
            .arg(m_density_buf)
            .arg(m_pressure_buf)
            .arg(m_sim_params->buffer()))
  {
    return false;
  }

  /* compute force kernel's arguments */
  if (!utils::ocl::KernelArgs(m_sph_compute_force_kernel, "sph_optimized_compute_force")
            .arg(m_particle_pos_buf.getCLID())
            .arg(m_density_buf)
            .arg(m_pressure_buf)
            .arg(m_force_buf)
            .arg(m_velocity_buf)
            .arg(m_sim_params->buffer()))
  {
    return false;
  }

  /* compute step kernel's arguments */
  if (!utils::ocl::KernelArgs(m_sph_compute_step_kernel, "sph_optimized_compute_step")
            .arg(m_particle_pos_buf.getCLID())
            .arg(m_force_buf)
            .arg(m_velocity_buf)
            .arg(m_prev_velocity_buf)
            .arg(m_sim_params->buffer()))
  {
    return false;
  }

  m_stats.addEvent(m_sim_params->upload(m_queue), "Update SPHOptimized params");

  // run the reset kernel to initialize particle data
#if 1
  cl_mem buffers[] = { m_particle_pos_buf.getCLID() };
  utils::ocl::GLSyncHandler sync(m_queue, sizeof(buffers) / sizeof(buffers[0]), buffers);
  if (!sync) return false;
#else
  utils::ocl::GLSyncHandler sync(m_queue, m_particle_pos_buf.getCLID());
  if (!sync) return false;
#endif

  // run the initialization kernel
  return deviceCall(m_sph_reset_kernel, "sph_optimized_reset",
                    global_work_size, local_work_size,
                    m_particle_pos_buf.getCLID(),
                    m_velocity_buf,
                    m_prev_velocity_buf,
                    m_pressure_buf,
                    m_density_buf,
                    m_force_buf,
                    m_sim_params->buffer());
}


void SPHOptimized::update_impl(float time_step)
{
  if (m_cl_timer.isReady()) m_cl_timer.start();

  m_stats.addEvent(m_sim_params->upload(m_queue), "Update SPHOptimized params");

#ifdef DEBUG_PARAMS
  m_sim_params->enqueuePrintGPUBuffer(m_queue);
#endif

  // synchronise with OpenGL
  cl_mem buffers[] = { m_particle_pos_buf.getCLID() };

  utils::ocl::GLSyncHandler sync(m_queue, sizeof(buffers) / sizeof(buffers[0]), buffers);
  if (!sync) return;

  // calculate global work size based on local work size
  size_t part_cnt = m_sim_params->particleCount();
  size_t block_size = BLOCK_SIZE;
  size_t grid_size = (part_cnt + block_size - 1) / block_size;
  size_t local_work_size = block_size;
  size_t global_work_size = grid_size * block_size;

  // compute pressure
  runKernel(m_sph_compute_pressure_kernel, "sph_optimized_compute_pressure",
            global_work_size, local_work_size);

  // compute force
  runKernel(m_sph_compute_force_kernel, "sph_optimized_compute_force",
            global_work_size, local_work_size);

  // integrate
  runKernel(m_sph_compute_step_kernel, "sph_optimized_compute_step",
            global_work_size, local_work_size);

  // advance simulation time
  m_sim_params->advanceSimulation(time_step);

  if (m_cl_timer.isReady()) m_cl_timer.stop();
}
