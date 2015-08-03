#include "test/sph_test.h"
#include "test/sph_test_adapter.h"
#include "utils/ocl.h"

#include <limits>
#include <iomanip>


#define DECL_DIFF_FUNCTOR(name, func, type) \
  template <typename TestAdapter1, typename TestAdapter2> \
  struct name ## Diff \
  { \
    name ## Diff(TestAdapter1 & sim1, TestAdapter2 & sim2) \
      : m_sim1(sim1), m_sim2(sim2) \
    { } \
    \
    type diff(int i) const \
    { return abs(m_sim1.func(i) - m_sim2.func(i)); } \
    \
    void print(std::ostream & os, int i) const \
    { os << m_sim1.func(i) << " <-> " << m_sim2.func(i) << std::endl; } \
    \
    const TestAdapter1 & m_sim1; \
    const TestAdapter2 & m_sim2; \
  };




namespace {

inline cl_float4 operator-(const cl_float4 & v1, const cl_float4 & v2)
{
  return { (v1.s[0] - v2.s[0]), (v1.s[1] - v2.s[1]),
           (v1.s[2] - v2.s[2]), (v1.s[3] - v2.s[3]) };
}

inline cl_float4 abs(const cl_float4 & v)
{
  return { fabsf(v.s[0]), fabsf(v.s[1]), fabsf(v.s[2]), fabsf(v.s[3]) };
}

inline cl_float abs(const cl_float v)
{
  return fabsf(v);
}


// Position attribute of particle
DECL_DIFF_FUNCTOR(Pos, pos, cl_float4)

// Density attribute of particle
DECL_DIFF_FUNCTOR(Density, density, cl_float)

// Pressure attribute of particle
DECL_DIFF_FUNCTOR(Pressure, pressure, cl_float)

// Forces attribute of particle
DECL_DIFF_FUNCTOR(Forces, forces, cl_float4)

// Velocity attribute of particle
DECL_DIFF_FUNCTOR(Vel, vel, cl_float4)

// Previous velocity attribute of particle
DECL_DIFF_FUNCTOR(PrevVel, prevvel, cl_float4)

}


namespace test {

template <typename TestAdapter1, typename TestAdapter2>
template <typename DiffFunc>
bool SPHTest<TestAdapter1, TestAdapter2>::cmpAndPrintFloat
     (std::ostream & os, const unsigned int max_diffs, const DiffFunc & f)
{
  cl_float sum = 0.0f;
  cl_float min = std::numeric_limits<cl_float>::max();
  cl_float max = std::numeric_limits<cl_float>::min();
  cl_uint count = 0;
  int part_cnt = m_params->particleCount();

  for (int i = 0; i < part_cnt; ++i)
  {
    cl_float diff = f.diff(i);

    if (diff > m_max_tolerance)
    {
      sum += diff;
      if (diff > max) max = diff;
      if (diff < min) min = diff;
      ++count;
      if (count <= max_diffs) f.print(os, i);
    }
  }

  if (count == 0)
  {
    os << "Results are exactly equal" << std::endl;
  }
  else
  {
    os << std::fixed << std::setprecision(10)
       << "Results VARY: Avg=" << (sum / ((cl_float) count))
       << ", Sum=" << sum << ", Min=" << min << ", Max=" << max
       << ", Count=" << count << std::endl;
  }

  return true;
}


template <typename TestAdapter1, typename TestAdapter2>
template <typename DiffFunc>
bool SPHTest<TestAdapter1, TestAdapter2>::cmpAndPrintFloat4
     (std::ostream & os, const unsigned int max_diffs, const DiffFunc & f)
{
  static constexpr cl_float MAX = std::numeric_limits<cl_float>::max();
  static constexpr cl_float MIN = std::numeric_limits<cl_float>::min();

  cl_float4 sum = { 0.0f, 0.0f, 0.0f, 0.0f };
  cl_float4 min = { MAX, MAX, MAX, MAX };
  cl_float4 max = { MIN, MIN, MIN, MIN };
  cl_uint4 count = { 0, 0, 0, 0 };
  unsigned int prints = 0;
  int part_cnt = m_params->particleCount();

  for (int j = 0; j < part_cnt; ++j)
  {
    cl_float4 diff = f.diff(j);

    int b = 0;

    for (int i = 0; i < 4; ++i)
    {
      if (diff.s[i] > m_max_tolerance)
      {
        sum.s[i] += diff.s[i];
        if (diff.s[i] > max.s[i]) max.s[i] = diff.s[i];
        if (diff.s[i] < min.s[i]) min.s[i] = diff.s[i];
        ++count.s[i];
        b = 1;
      }
    }

    prints += b;
    if ((b) && (prints <= max_diffs)) f.print(os, j);
  }

  for (int i = 0; i < 4; ++i)
  {
    if (count.s[i] == 0)
    {
      os << "Results for the " << i << ". component are exactly equal" << std::endl;
    }
    else
    {
      os << std::fixed << std::setprecision(10)
         << "Results for the " << i << ". component VARY: Avg=" << (sum.s[i] / ((cl_float) count.s[i]))
         << ", Sum=" << sum.s[i] << ", Min=" << min.s[i] << ", Max=" << max.s[i]
         << ", Count=" << count.s[i] << std::endl;
    }
  }

  return true;
}


template <typename TestAdapter1, typename TestAdapter2>
bool SPHTest<TestAdapter1, TestAdapter2>::cmpAndPrintParticles
     (std::ostream & os, const unsigned int max_print)
{
  if ((!m_sim1.beginTraverse()) || (!m_sim2.beginTraverse()))
  {
    WARNM("Failed to start traversing results");
    return false;
  }

  m_sim1.sort();
  m_sim2.sort();

  bool ret = true;

  os << "Position component:" << std::endl;
  ret &= cmpAndPrintFloat4(os, max_print, PosDiff<TestAdapter1, TestAdapter2>(m_sim1, m_sim2));

  os << "Density component:" << std::endl;
  ret &= cmpAndPrintFloat (os, max_print, DensityDiff<TestAdapter1, TestAdapter2>(m_sim1, m_sim2));

  os << "Pressure component:" << std::endl;
  ret &= cmpAndPrintFloat (os, max_print, PressureDiff<TestAdapter1, TestAdapter2>(m_sim1, m_sim2));

  os << "Force component:" << std::endl;
  ret &= cmpAndPrintFloat4(os, max_print, ForcesDiff<TestAdapter1, TestAdapter2>(m_sim1, m_sim2));

  os << "Velocity component:" << std::endl;
  ret &= cmpAndPrintFloat4(os, max_print, VelDiff<TestAdapter1, TestAdapter2>(m_sim1, m_sim2));

  os << "Previous velocity component:" << std::endl;
  ret &= cmpAndPrintFloat4(os, max_print, PrevVelDiff<TestAdapter1, TestAdapter2>(m_sim1, m_sim2));

  if ((!m_sim1.endTraverse()) || (!m_sim2.endTraverse()))
  {
    WARNM("Failed to end traversing results");
    return false;
  }


  return ret;
}


template <typename TestAdapter1, typename TestAdapter2>
bool SPHTest<TestAdapter1, TestAdapter2>::reset(int nparts, int grid_size)
{
  if (!utils::ocl::initCLGLContext(m_cl_ctx)) return false;

  m_params.reset(new SPHParams(m_cl_ctx));
  //m_params->setSeed(10000);
  m_params->setParticleCount(nparts);
  m_params->setGridSize(grid_size, grid_size, grid_size);

  return m_sim1.init(m_cl_ctx, m_params.get()) &&
         m_sim2.init(m_cl_ctx, m_params.get());
}


template <typename TestAdapter1, typename TestAdapter2>
bool SPHTest<TestAdapter1, TestAdapter2>::run
     (std::ostream & os, const int max_print, const int niters)
{
  // print simulation parameters
  //m_sim1.dumpParameters(os);
  //m_sim2.dumpParameters(os);
  os << "Simulation parameters:\n" << *m_params << std::endl;

  // compare data before running the any simulation
  os << "================================================================" << std::endl;
  os << "Comparison of initial states:" << std::endl;
  os << "================================================================" << std::endl;

  cmpAndPrintParticles(os, max_print);

  // run both simulations
  for (int i = 0; i < niters; ++i)
  {
    m_sim1.update();
    m_sim2.update();
  }

  // compare particle positions on their outputs
  os << "================================================================" << std::endl;
  os << "Comparison of computed states:" << std::endl;
  os << "================================================================" << std::endl;

  cmpAndPrintParticles(os, max_print);

  return true;
}

// explicit instantiation of the usefull tests
template class SPHTest<SPHOptimizedTestAdapter, SPHOptimizedTestAdapter>;
template class SPHTest<SPHUniformGridTestAdapter, SPHUniformGridTestAdapter>;
template class SPHTest<SPHNaiveTestAdapter, SPHNaiveTestAdapter>;

template class SPHTest<SPHUniformGridTestAdapter, SPHOptimizedTestAdapter>;
template class SPHTest<SPHUniformGridTestAdapter, SPHNaiveTestAdapter>;
template class SPHTest<SPHOptimizedTestAdapter, SPHNaiveTestAdapter>;

} // End of namespace test
