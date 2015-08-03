#include "test/sph_naive_test_adapter.h"



namespace test {

bool SPHNaiveTestAdapter::beginTraverse_impl(void)
{
  SPHNaive *ps = m_ps.get();
  assert(ps != nullptr);

  m_pos_ptr = mapPositionBuf();
  m_density_ptr = mapDensityBuf();
  m_pressure_ptr = mapPressureBuf();
  m_forces_ptr = mapForceBuf();
  m_vel_ptr = mapVelocityBuf();
  m_prevvel_ptr = mapPrevVelocityBuf();

  return true;
}


bool SPHNaiveTestAdapter::endTraverse_impl(void)
{
  SPHNaive *ps = m_ps.get();
  assert(ps != nullptr);

  unmapPositionBuf(m_pos_ptr);
  unmapDensityBuf(m_density_ptr);
  unmapPressureBuf(m_pressure_ptr);
  unmapForceBuf(m_forces_ptr);
  unmapVelocityBuf(m_vel_ptr);
  unmapPrevVelocityBuf(m_prevvel_ptr);

  return true;
}


bool SPHNaiveTestAdapter::dumpParticle_impl(std::ostream & /* os */, int /* i */) const
{
  return true;
}

} // End of namespace test
