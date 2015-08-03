#ifndef DEBUG_PARTICLE_SYSTEM_H
#define DEBUG_PARTICLE_SYSTEM_H

#include "core/particle_system.h"

#include <boost/compute/program.hpp>
#include <boost/compute/kernel.hpp>


/**
 * A dummy debugging particle system
 */
class DebugParticleSystem : public ParticleSystem
{
  public:
    DebugParticleSystem(const boost::compute::context & ctx, SPHParams *params)
      : ParticleSystem(ctx, params)
      , m_polar_spiral_prog()
      , m_gen_rand_particles_prog()
      , m_polar_spiral_kernel()
      , m_gen_rand_particles_kernel()
      , m_spiral(true)
    {
      // initialize OpenCL context, compile kernels
      if (!init())
      {
        throw std::runtime_error("Failed to construct TestSystem: OpenCL initialization failed");
      }
    }

    bool spiral(void) const { return m_spiral; }
    void setSpiral(bool spiral = true) { m_spiral = spiral; }
    bool toggleSpiral(void) { return m_spiral = !m_spiral; }


  protected:
    virtual bool reset_impl(const io::VoxelMesh * /* mesh */) override;
    virtual void update_impl(float time_step) override;

  private:
    // initializes OpenCL context
    bool init(void);

  private:
    // OpenCL programs
    boost::compute::program m_polar_spiral_prog;
    boost::compute::program m_gen_rand_particles_prog;

    // OpenCL kernels
    boost::compute::kernel m_polar_spiral_kernel;        // a kernel to generate polar spiral
    boost::compute::kernel m_gen_rand_particles_kernel;  // testing kernel

    bool m_spiral;  // whether to generate archimedean spiral or just random positions
};

#endif // DEBUG_PARTICLE_SYSTEM_H
