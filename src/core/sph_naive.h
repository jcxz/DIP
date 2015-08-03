#ifndef SPH_NAIVE_H
#define SPH_NAIVE_H

#include "core/fluid_particle_system.h"

#include <boost/compute/program.hpp>
#include <boost/compute/kernel.hpp>
#include <boost/compute/buffer.hpp>

namespace test { class SPHNaiveTestAdapter; }


class SPHNaive : public FluidParticleSystem
{
  public:
    SPHNaive(const boost::compute::context & ctx, SPHParams *params)
      : FluidParticleSystem(ctx, params)
      , m_sph_prog()
      , m_sph_reset_kernel()
      , m_sph_compute_step_kernel()
      , m_sph_compute_force_kernel()
      , m_sph_compute_pressure_kernel()
      , m_velocity_buf()
      , m_pressure_buf()
      , m_density_buf()
      , m_force_buf()
      , m_prev_velocity_buf()
    {
      // initialize OpenCL context, compile kernels
      if (!init())
      {
        throw std::runtime_error("Failed to construct SPHNaive: OpenCL initialization failed");
      }
    }

  protected:
    virtual bool reset_impl(const io::VoxelMesh * /* mesh */) override;
    virtual void update_impl(float time_step) override;

  private:
    // initializes the OpenCL program and kernel for SPH simulation
    bool init(void);

  private:
    // OpenCL programs
    boost::compute::program m_sph_prog;

    // OpenCL kernels
    boost::compute::kernel m_sph_reset_kernel;             // a reset kernel for the SPH simulation
    boost::compute::kernel m_sph_compute_step_kernel;      // a kernel to compute a single SPH step
    boost::compute::kernel m_sph_compute_force_kernel;     // kernel for computing forces
    boost::compute::kernel m_sph_compute_pressure_kernel;  // kernel for computing the pressure inside of the fluid

    // buffers for SPH simulation
    boost::compute::buffer m_velocity_buf;
    boost::compute::buffer m_pressure_buf;
    boost::compute::buffer m_density_buf;
    boost::compute::buffer m_force_buf;
    boost::compute::buffer m_prev_velocity_buf;

    // give the test adapter class access to the internals
    friend class test::SPHNaiveTestAdapter;
};

#endif // SPH_NAIVE_H
