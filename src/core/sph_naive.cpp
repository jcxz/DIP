#include "core/particle_system_internal.h"
#include "core/sph_naive.h"
#include "utils/macros.h"
#include "utils/debug.h"

#include <QFile>
#include <ctime>



bool SPHNaive::init(void)
{
  try
  {
    QFile f(":/src/opencl/sph_naive.cl");
    if (!f.open(QFile::ReadOnly)) return false;
    //m_sph_prog = boost::compute::program::build_with_source(f.readAll().toStdString(), m_cl_ctx);
    m_sph_prog = boost::compute::program::create_with_source(f.readAll().toStdString(), m_cl_ctx);
    m_sph_prog.build();

    m_sph_reset_kernel = m_sph_prog.create_kernel("sph_naive_reset");
    m_sph_compute_step_kernel = m_sph_prog.create_kernel("sph_naive_compute_step");
    m_sph_compute_force_kernel = m_sph_prog.create_kernel("sph_naive_compute_force");
    m_sph_compute_pressure_kernel = m_sph_prog.create_kernel("sph_naive_compute_pressure");
  }
  catch (const boost::compute::opencl_error & e)
  {
    ERRORM("SPHNaive: An OpenCL error occured: " << e.what());
    if (e.error_code() == CL_BUILD_PROGRAM_FAILURE) ERRORM(m_sph_prog.build_log());
    return false;
  }
  catch (const std::exception & e)
  {
    ERRORM("SPHNaive: An unexpected error occured during SPHNaive initialization: "
           << e.what());
    return false;
  }

  return true;
}


bool SPHNaive::reset_impl(const io::VoxelMesh * /* mesh */)
{
  size_t part_cnt = m_sim_params->particleCount();
  m_sim_params->setTime(0.0f);
  m_sim_params->clearEffectFlags();

  // initialize OpenGL shared buffer with particle positions
  if (!m_particle_pos_buf.bufferData(nullptr, part_cnt * sizeof(cl_float4),
                                     utils::ocl::GLBuffer::WRITE_ONLY))
  {
    ERRORM("SPHNaive: Failed to initialize position GLBuffer");
    return false;
  }

  // allocate buffers on GPU
  try
  {
    m_velocity_buf = boost::compute::buffer(m_cl_ctx, part_cnt * sizeof(cl_float4));
    m_prev_velocity_buf = boost::compute::buffer(m_cl_ctx, part_cnt * sizeof(cl_float4));
    m_pressure_buf = boost::compute::buffer(m_cl_ctx, part_cnt * sizeof(cl_float));
    m_density_buf = boost::compute::buffer(m_cl_ctx, part_cnt * sizeof(cl_float));
    m_force_buf = boost::compute::buffer(m_cl_ctx, part_cnt * sizeof(cl_float4));
  }
  catch (const boost::compute::opencl_error & e)
  {
    ERRORM("SPHNaive: An OpenCL error occured: " << e.what());
    return false;
  }
  catch (const std::exception & e)
  {
    ERRORM("SPHNaive: An unexpected error occured during SPHNaive initialization: "
           << e.what());
    return false;
  }

  // set kernel parameters that do not change once the simulation is prepared

#if 0
#define SIM_SCALE ((cl_float) (0.004f))
#define SMOOTH_RADIUS ((cl_float) (0.01f))
#define RADIUS2 ((cl_float) ((SMOOTH_RADIUS) * (SMOOTH_RADIUS)))

// m = Ro * (V / n)
// Ro   ... hustota tekutiny
// V    ... objem
// n    ... pocet castic
#define MASS ((cl_float) (0.00020543f))
#define POLYKERN ((cl_float) (315.0f / (64.0f * 3.141592 * pow(SMOOTH_RADIUS, 9))))
#define RESTDENSITY ((cl_float) (600.0f))
#define INTSTIFFNESS ((cl_float) (1.0f))
#define MASS_POLYKERN ((cl_float) ((MASS) * (POLYKERN)))

#define VISCOSITY ((cl_float) (0.2f))
#define LAPKERN ((cl_float) (45.0f / (3.141592 * pow(SMOOTH_RADIUS, 6))))
#define VTERM ((cl_float) ((LAPKERN) * (VISCOSITY)))
#define SPIKEYKERN ((cl_float) (-45.0f / (3.141592 * pow(SMOOTH_RADIUS, 6))))
#define SPIKEYKERN_HALF ((cl_float) ((SPIKEYKERN) * (-0.5f)))

#define SLOPE ((cl_float) (0.0f))        // ???
#define LEFTWAVE ((cl_float) (0.0f))     // ???
#define RIGHTWAVE ((cl_float) (0.0f))   // ???
#define DELTATIME ((cl_float) (.003f))   // ???
#define LIMIT ((cl_float) (200.0f))
#define EXTSTIFFNESS ((cl_float) (10000.0f))
#define EXTDAMPING ((cl_float) (256.0f))
//#define RADIUS ((cl_float) (0.004f))
#define RADIUS ((cl_float) (1.0f / 96.0f))
#else



#define SIM_SCALE ((tFloat) (0.005f))
#define SMOOTH_RADIUS ((tFloat) (0.01f))
#define RADIUS2 ((tFloat) ((SMOOTH_RADIUS) * (SMOOTH_RADIUS)))

// m = Ro * (V / n)
// Ro   ... hustota tekutiny
// V    ... objem
// n    ... pocet castic
#define MASS ((tFloat) (0.00020543f))
#define POLYKERN ((tFloat) (315.0f / (64.0f * 3.141592 * pow(SMOOTH_RADIUS, 9))))
#define RESTDENSITY ((tFloat) (600.0f))
#define INTSTIFFNESS ((tFloat) (1.5f))
#define MASS_POLYKERN ((tFloat) ((MASS) * (POLYKERN)))

#define VISCOSITY ((tFloat) (0.35f))
#define LAPKERN ((tFloat) (45.0f / (3.141592 * pow(SMOOTH_RADIUS, 6))))
#define VTERM ((tFloat) ((LAPKERN) * (VISCOSITY)))
#define SPIKEYKERN ((tFloat) (-45.0f / (3.141592 * pow(SMOOTH_RADIUS, 6))))
#define SPIKEYKERN_HALF ((tFloat) ((SPIKEYKERN) * (-0.5f)))

#define SLOPE ((tFloat) (0.0f))        // ???
#define LEFTWAVE ((tFloat) (0.0f))     // ???
#define RIGHTWAVE ((tFloat) (0.0f))   // ???
#define DELTATIME ((tFloat) (0.003f))   // ???
#define LIMIT ((tFloat) (150.0f))
#define EXTSTIFFNESS ((tFloat) (50000.0f))
#define EXTDAMPING ((tFloat) (100.0f))

#define RADIUS ((tFloat) (1.0f / 96.0f))
//#define RADIUS ((tFloat) (0.02f))
#endif

  // compute pressure kernel's arguments
  if (!utils::ocl::KernelArgs(m_sph_compute_pressure_kernel, "sph_naive_compute_pressure")
            .arg(m_particle_pos_buf.getCLID())
            .arg(m_density_buf)
            .arg(m_pressure_buf)
            .arg(SIM_SCALE)
            //.arg(SMOOTH_RADIUS)
            .arg(RADIUS2)
            //.arg(MASS)
            //.arg(POLYKERN)
            .arg(MASS_POLYKERN)
            .arg(RESTDENSITY)
            .arg(INTSTIFFNESS)
            .arg((unsigned int) (part_cnt)))
  {
    return false;
  }

  // compute force kernel's arguments
  if (!utils::ocl::KernelArgs(m_sph_compute_force_kernel, "sph_naive_compute_force")
            .arg(m_particle_pos_buf.getCLID())
            .arg(m_density_buf)
            .arg(m_pressure_buf)
            .arg(m_force_buf)
            .arg(m_velocity_buf)
            .arg(SIM_SCALE)
            .arg(SMOOTH_RADIUS)
            .arg(RADIUS2)
            //.arg(VISCOSITY)
            //.arg(LAPKERN)
            .arg(VTERM)
            .arg(SPIKEYKERN_HALF)
            .arg((unsigned int) (part_cnt)))
  {
    return false;
  }

  // compute step kernel's arguments
  if (!utils::ocl::KernelArgs(m_sph_compute_step_kernel, "sph_naive_compute_step")
            .arg(m_particle_pos_buf.getCLID())
            .arg(m_force_buf)
            .arg(m_velocity_buf)
            .arg(m_prev_velocity_buf)
            .arg(SLOPE)
            .arg(LEFTWAVE)
            .arg(RIGHTWAVE)
            .arg(DELTATIME)
            .arg(LIMIT)
            .arg(EXTSTIFFNESS)
            .arg(EXTDAMPING)
            .arg(RADIUS)
            .arg(m_sim_params->volumeMin())
            .arg(m_sim_params->volumeMax())
            .arg(SIM_SCALE)
            .arg(MASS))
  {
    return false;
  }

  // run the reset kernel to initialize particle data
#if 1
  cl_mem buffers[] = { m_particle_pos_buf.getCLID() };
  utils::ocl::GLSyncHandler sync(m_queue, sizeof(buffers) / sizeof(buffers[0]), buffers);
  if (!sync) return false;
#else
  utils::ocl::GLSyncHandler sync(m_queue, m_particle_pos_buf.getCLID());
  if (!sync) return false;
#endif

  return deviceCall2(m_sph_reset_kernel, "sph_naive_reset",
                     part_cnt,
                     m_particle_pos_buf.getCLID(),
                     m_velocity_buf,
                     m_prev_velocity_buf,
                     m_pressure_buf,
                     m_density_buf,
                     m_force_buf,
                     m_sim_params->volumeMin(),
                     m_sim_params->volumeMax(),
                     m_sim_params->seed());
}


void SPHNaive::update_impl(float time_step)
{
  if (m_cl_timer.isReady()) m_cl_timer.start();

  // set kernel arguments that change every frame
  m_sph_compute_step_kernel.set_arg(16, (cl_float) (m_sim_params->time()));
  m_sph_compute_step_kernel.set_arg(17, (cl_uint) (m_sim_params->effectFlags()));
  //m_sph_compute_step_kernel.setArg(18, calcGravitationVector(m_rx, m_ry));

  // synchronise with OpenGL
  cl_mem buffers[] = { m_particle_pos_buf.getCLID() };

  utils::ocl::GLSyncHandler sync(m_queue, sizeof(buffers) / sizeof(buffers[0]), buffers);
  if (!sync) return;

  size_t part_cnt = m_sim_params->particleCount();

  // compute pressure
  runKernel(m_sph_compute_pressure_kernel, "sph_naive_compute_pressure", part_cnt);

  // compute force
  runKernel(m_sph_compute_force_kernel, "sph_naive_compute_force", part_cnt);

  // integrate
  runKernel(m_sph_compute_step_kernel, "sph_naive_compute_step", part_cnt);

  // advance simulation time
  m_sim_params->advanceSimulation(time_step);

  if (m_cl_timer.isReady()) m_cl_timer.stop();
}
