#ifndef FLUID_PARTICLE_SYSTEM_H
#define FLUID_PARTICLE_SYSTEM_H

#include "core/particle_system.h"

#include <boost/compute/program.hpp>
#include <boost/compute/kernel.hpp>
#include <boost/compute/buffer.hpp>


class FluidParticleSystem : public ParticleSystem
{
  public:
    FluidParticleSystem(const boost::compute::context & ctx, SPHParams *params)
      : ParticleSystem(ctx, params)
      , m_rx(0)
      , m_ry(0)
    { }

    void setRotation(float rx, float ry) { m_rx = rx; m_ry = ry; }

  protected:
    cl_float4 calcGravitationVector(int rx, int ry);

  protected:
    int m_rx;
    int m_ry;
};

#endif // FLUID_PARTICLE_SYSTEM_H
