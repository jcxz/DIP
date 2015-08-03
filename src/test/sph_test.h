#ifndef SPH_TEST_H
#define SPH_TEST_H

#include "test/sph_optimized_test_adapter.h"
#include "test/sph_uniform_grid_test_adapter.h"
#include "test/sph_naive_test_adapter.h"

#include <stdexcept>


namespace test {

template <typename TestAdapter1, typename TestAdapter2>
class SPHTest
{
  private:
    static constexpr float DEF_TOLERANCE = 0.0f;

  public:
    explicit SPHTest(int nparts, int grid_size = 32)
      : m_cl_ctx()
      , m_params(nullptr)
      , m_max_tolerance(DEF_TOLERANCE)
      , m_sim1()
      , m_sim2()
    {
      if (!reset(nparts, grid_size))
        throw std::runtime_error("Failed to initialize simulation engines");
    }

    void setMaxTolerance(float tol) { m_max_tolerance = tol; }
    bool reset(int nparts, int grid_size);
    bool run(std::ostream & os, const int max_print = 1, const int niters = 1);

  private:
    template <typename DiffFunc>
    bool cmpAndPrintFloat(std::ostream & os, const unsigned int n, const DiffFunc & f);

    template <typename DiffFunc>
    bool cmpAndPrintFloat4(std::ostream & os, const unsigned int n, const DiffFunc & f);

    bool cmpAndPrintParticles(std::ostream & os, const unsigned int max_print);

  private:
    boost::compute::context m_cl_ctx;
    std::unique_ptr<SPHParams> m_params;
    float m_max_tolerance;
    TestAdapter1 m_sim1;
    TestAdapter2 m_sim2;
};

// Definition of available tests
typedef SPHTest<SPHOptimizedTestAdapter, SPHOptimizedTestAdapter>     SPHOptimized_vs_SPHOptimized;
typedef SPHTest<SPHUniformGridTestAdapter, SPHUniformGridTestAdapter> SPHUniformGrid_vs_SPHUniformGrid;
typedef SPHTest<SPHNaiveTestAdapter, SPHNaiveTestAdapter>             SPHnaive_vs_SPHNaive;

typedef SPHTest<SPHUniformGridTestAdapter, SPHOptimizedTestAdapter>   SPHUniformGrid_vs_SPHOptimized;
typedef SPHTest<SPHUniformGridTestAdapter, SPHNaiveTestAdapter>       SPHUniformGrid_vs_SPHNaive;
typedef SPHTest<SPHOptimizedTestAdapter, SPHNaiveTestAdapter>         SPHOptimized_vs_SPHNaive;

} // End of namespace test

#endif // SPH_TEST_H
