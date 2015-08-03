#ifndef SPH_OPTIMIZED_TEST_ADAPTER_H
#define SPH_OPTIMIZED_TEST_ADAPTER_H

#include "test/sph_test_adapter.h"
#include "core/sph_optimized.h"

#include <memory>


namespace test {

class SPHOptimizedTestAdapter : public SPHTestAdapter<SPHOptimizedTestAdapter>
{
    friend class SPHTestAdapter;

  public:
    SPHOptimizedTestAdapter(void)
      : m_ps(nullptr)
      , m_pos_ptr(nullptr)
      , m_density_ptr(nullptr)
      , m_pressure_ptr(nullptr)
      , m_forces_ptr(nullptr)
      , m_vel_ptr(nullptr)
      , m_prevvel_ptr(nullptr)
    { }

    ~SPHOptimizedTestAdapter(void)
    { endTraverse_impl(); }

  private:
    bool init_impl(const boost::compute::context & ctx, SPHParams *params)
    {
      m_ps.reset(new SPHOptimized(ctx, params));
      m_ps->setPrintStatsOnExit(false); // disable printing performance statistics
      return m_ps->reset();
    }

    bool update_impl(void) { assert(m_ps.get() != nullptr); m_ps->update(); return true; }

    bool beginTraverse_impl(void);

    cl_float4 pos_impl(int idx) const { assert(m_pos_ptr != nullptr); return m_pos_ptr[idx]; }
    cl_float density_impl(int idx) const { assert(m_density_ptr != nullptr); return m_density_ptr[idx]; }
    cl_float pressure_impl(int idx) const { assert(m_pressure_ptr != nullptr); return m_pressure_ptr[idx]; }
    cl_float4 forces_impl(int idx) const { assert(m_forces_ptr != nullptr); return m_forces_ptr[idx]; }
    cl_float4 vel_impl(int idx) const { assert(m_vel_ptr != nullptr); return m_vel_ptr[idx]; }
    cl_float4 prevvel_impl(int idx) const { assert(m_prevvel_ptr != nullptr); return m_prevvel_ptr[idx]; }

    bool endTraverse_impl(void);

    bool dumpParameters_impl(std::ostream & os) const
    {
      assert(m_ps.get() != nullptr);
      os << "SPHOptimized simulation parameters" << std::endl;
      os << *m_ps->m_sim_params << std::endl;
      return true;
    }

    size_t particleCount_impl(void) const { return m_ps->m_sim_params->particleCount(); }

    bool dumpParticle_impl(std::ostream & /* os */, int /* i */) const;

  private:
    //cl_float4 *mapVelocityBuf(void)
    //{ return (cl_float4 *) m_queue.enqueue_map_buffer(m_velocity_buf, CL_MAP_READ, 0, m_num_particles * sizeof(T)); }

    cl_float4 *mapVelocityBuf(void)     { return mapBuffer<cl_float4>(m_ps->m_queue, m_ps->m_velocity_buf,      particleCount_impl()); }
    cl_float  *mapPressureBuf(void)     { return mapBuffer<cl_float> (m_ps->m_queue, m_ps->m_pressure_buf,      particleCount_impl()); }
    cl_float  *mapDensityBuf(void)      { return mapBuffer<cl_float> (m_ps->m_queue, m_ps->m_density_buf,       particleCount_impl()); }
    cl_float4 *mapForceBuf(void)        { return mapBuffer<cl_float4>(m_ps->m_queue, m_ps->m_force_buf,         particleCount_impl()); }
    cl_float4 *mapPrevVelocityBuf(void) { return mapBuffer<cl_float4>(m_ps->m_queue, m_ps->m_prev_velocity_buf, particleCount_impl()); }
    cl_float4 *mapPositionBuf(void)     { return mapBuffer<cl_float4>(m_ps->m_queue, m_ps->m_particle_pos_buf,  particleCount_impl()); }

    void unmapVelocityBuf(const cl_float4 *ptr)     { return unmapBuffer(m_ps->m_queue, m_ps->m_velocity_buf,      ptr); }
    void unmapPressureBuf(const cl_float *ptr)      { return unmapBuffer(m_ps->m_queue, m_ps->m_pressure_buf,      ptr); }
    void unmapDensityBuf(const cl_float *ptr)       { return unmapBuffer(m_ps->m_queue, m_ps->m_density_buf,       ptr); }
    void unmapForceBuf(const cl_float4 *ptr)        { return unmapBuffer(m_ps->m_queue, m_ps->m_force_buf,         ptr); }
    void unmapPrevVelocityBuf(const cl_float4 *ptr) { return unmapBuffer(m_ps->m_queue, m_ps->m_prev_velocity_buf, ptr); }
    void unmapPositionBuf(const cl_float4 *ptr)     { return unmapBuffer(m_ps->m_queue, m_ps->m_particle_pos_buf,  ptr); }

  private:
    SPHOptimizedTestAdapter(const SPHOptimizedTestAdapter & );
    SPHOptimizedTestAdapter & operator=(const SPHOptimizedTestAdapter & );

  private:
    std::unique_ptr<SPHOptimized> m_ps;
    cl_float4 *m_pos_ptr;
    cl_float  *m_density_ptr;
    cl_float  *m_pressure_ptr;
    cl_float4 *m_forces_ptr;
    cl_float4 *m_vel_ptr;
    cl_float4 *m_prevvel_ptr;
};

} // End of namespace test

#endif // SPH_OPTIMIZED_TEST_ADAPTER_H
