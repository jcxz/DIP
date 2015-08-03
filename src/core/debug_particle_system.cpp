#include "core/debug_particle_system.h"
#include "utils/macros.h"
#include "utils/debug.h"

#include <QFile>
#include <ctime>



bool DebugParticleSystem::init(void)
{
  try
  {
    /* create gen_rand_particles program and kernel */
    {
      QFile f(":/src/opencl/gen_rand_particles.cl");
      if (!f.open(QFile::ReadOnly)) return false;
      m_gen_rand_particles_prog = boost::compute::program::build_with_source(f.readAll().toStdString(), m_cl_ctx);
      m_gen_rand_particles_kernel = m_gen_rand_particles_prog.create_kernel("gen_part_positions");
    }

    /* create polar_spiral program and kernel */
    {
      QFile f(":/src/opencl/polar_spiral.cl");
      if (!f.open(QFile::ReadOnly)) return false;
      m_polar_spiral_prog = boost::compute::program::build_with_source(f.readAll().toStdString(), m_cl_ctx);
      m_polar_spiral_kernel = m_polar_spiral_prog.create_kernel("polar_spiral");
    }
  }
  catch (const boost::compute::opencl_error & e)
  {
    ERRORM("DebugParticleSystem: An OpenCL error occured: " << e.what());
    return false;
  }
  catch (const std::exception & e)
  {
    ERRORM("DebugParticleSystem: An unexpected error occured during ParticleSystem initialization: "
           << e.what());
    return false;
  }

  return true;
}


bool DebugParticleSystem::reset_impl(const io::VoxelMesh * /* mesh */)
{
  INFOM("DebugParticleSystem: Initializing Test Simulation data");

  unsigned int part_num = m_sim_params->particleCount();

  /* allocate the buffers to hold particle positions and colors */
  if (!m_particle_pos_buf.bufferData(nullptr, part_num * sizeof(cl_float4), utils::ocl::GLBuffer::WRITE_ONLY))
  {
    ERRORM("DebugParticleSystem: Failed to initialize position GLBuffer");
    return false;
  }

  if (!m_particle_col_buf.bufferData(nullptr, part_num * sizeof(cl_float4), utils::ocl::GLBuffer::WRITE_ONLY))
  {
    ERRORM("DebugParticleSystem: Failed to initialize color GLBuffer");
    return false;
  }

  /* initialize kernel arguments */
  if (!utils::ocl::KernelArgs(m_gen_rand_particles_kernel, "m_gen_rand_particles_kernel")
            .arg(m_particle_pos_buf.getCLID())
            .arg(m_particle_col_buf.getCLID()))
            //.arg((cl_ulong) (time(nullptr))))
  {
    return false;
  }

  cl_float3 position = { 0.0f, 0.0f, 0.0f };
  if (!utils::ocl::KernelArgs(m_polar_spiral_kernel, "m_polar_spiral_kernel")
            .arg(m_particle_pos_buf.getCLID())
            .arg(m_particle_col_buf.getCLID())
            .arg(position))
  {
    return false;
  }

  m_sim_params->setTime(0.0f);

  return true;
}


void DebugParticleSystem::update_impl(float time_step)
{
  if (m_cl_timer.isReady()) m_cl_timer.start();

  size_t part_cnt = m_sim_params->particleCount();
  float sim_time = m_sim_params->time();

  if (m_spiral)
    m_polar_spiral_kernel.set_arg(3, cl_ulong(time(nullptr) + sim_time));
  else
    m_gen_rand_particles_kernel.set_arg(2, cl_ulong(time(nullptr) + sim_time));

  cl_mem buffers[] = { m_particle_pos_buf.getCLID(), m_particle_col_buf.getCLID() };

  utils::ocl::GLSyncHandler sync(m_queue, sizeof(buffers) / sizeof(buffers[0]), buffers);
  if (!sync) return;

  cl_int err = CL_SUCCESS;

  if (m_spiral)
  {
    err = clEnqueueNDRangeKernel(m_queue, m_polar_spiral_kernel, 1,
                                 nullptr, &part_cnt, nullptr,
                                 0, nullptr, m_stats.event("polar_spiral"));
  }
  else
  {
    err = clEnqueueNDRangeKernel(m_queue, m_gen_rand_particles_kernel, 1,
                                 nullptr, &part_cnt, nullptr,
                                 0, nullptr, m_stats.event("gen_part_positions"));
  }

  if (err != CL_SUCCESS)
  {
    WARNM("DebugParticleSystem: Failed to enqueue test simulation kernel: "
          << utils::ocl::errorToStr(err));
  }

  m_sim_params->advanceTime(time_step);  // 3.0f;

  if (m_cl_timer.isReady()) m_cl_timer.stop();
}
