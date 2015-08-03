#include "core/particle_system.h"
#include "utils/macros.h"
#include "utils/debug.h"

#include <boost/compute/interop/opengl/context.hpp>

#include <ctime>



bool ParticleSystem::runKernelHelper(cl_kernel kern, const char *name,
                                     const size_t *global_work_size,
                                     const size_t *local_work_size)
{
  cl_int err = clEnqueueNDRangeKernel(m_queue, kern, 1,
                                      nullptr, global_work_size, local_work_size,
                                      0, nullptr, m_stats.event(name));
  if (err != CL_SUCCESS)
  {
    WARNM("Failed to enqueue " << name << " kernel: " <<
          boost::compute::opencl_error::to_string(err) <<
          "(" << err << ")");
    return false;
  }

  return true;
}
