// TODO: blah blah blah
//   All this test adapter and all this test code should be rather
//   wasted and banned for eternity from the whole code base, because
//   it is just a pile of stinky crap created over a late night
//   delusion of pure genius.

#ifndef SPH_TEST_ADAPTER_H
#define SPH_TEST_ADAPTER_H

#include "utils/ocl.h"

#include <boost/compute/command_queue.hpp>


class SPHParams;

namespace test {

// This class defines an interface that all other test adapeters should adhere to
template <typename DerivedClass>
class SPHTestAdapter
{
  public:
    bool init(const boost::compute::context & ctx, SPHParams *params)
    { return static_cast<DerivedClass *>(this)->init_impl(ctx, params); }

    bool update(void)
    { return static_cast<DerivedClass *>(this)->update_impl(); }

    bool beginTraverse(void)
    { return static_cast<DerivedClass *>(this)->beginTraverse_impl(); }

    cl_float4 pos(int idx) const      { return static_cast<const DerivedClass *>(this)->pos_impl(idx); }
    cl_float  density(int idx) const  { return static_cast<const DerivedClass *>(this)->density_impl(idx); }
    cl_float  pressure(int idx) const { return static_cast<const DerivedClass *>(this)->pressure_impl(idx); }
    cl_float4 forces(int idx) const   { return static_cast<const DerivedClass *>(this)->forces_impl(idx); }
    cl_float4 vel(int idx) const      { return static_cast<const DerivedClass *>(this)->vel_impl(idx); }
    cl_float4 prevvel(int idx) const  { return static_cast<const DerivedClass *>(this)->prevvel_impl(idx); }

    bool endTraverse(void)
    { return static_cast<DerivedClass *>(this)->endTraverse_impl(); }

    bool dumpParameters(std::ostream & os) const
    { return static_cast<const DerivedClass *>(this)->dumpParameters_impl(os); }

    bool dumpParticle(std::ostream & os, int i) const
    { return static_cast<const DerivedClass *>(this)->dumpParticle_impl(os, i); }

    void sort(void)
    {
      std::cerr << __PRETTY_FUNCTION__ << std::endl;

      DerivedClass *dc = static_cast<DerivedClass *>(this);

      assert(dc->m_ps.get() != nullptr);

      assert(dc->m_pos_ptr      != nullptr);
      assert(dc->m_density_ptr  != nullptr);
      assert(dc->m_pressure_ptr != nullptr);
      assert(dc->m_forces_ptr   != nullptr);
      assert(dc->m_vel_ptr      != nullptr);
      assert(dc->m_prevvel_ptr  != nullptr);

      //auto f = [](const cl_float4 & lhs, const cl_float4 & rhs) {
      //  return (lhs.s[0] < rhs.s[0]) &&
      //         (lhs.s[1] < rhs.s[1]) &&
      //         (lhs.s[2] < rhs.s[2]) &&
      //         (lhs.s[3] < rhs.s[3]);
      //  //return memcmp(lhs.s, rhs.s, sizeof(cl_float4)) < 0;
      //  //return memcmp(lhs.s, rhs.s, sizeof(cl_float4)) >= 0;
      //};

      //std::sort(dc->m_pos_ptr,      dc->m_pos_ptr      +  dc->particleCount_impl(), f);
      //std::sort(dc->m_density_ptr,  dc->m_density_ptr  +  dc->particleCount_impl());
      //std::sort(dc->m_pressure_ptr, dc->m_pressure_ptr +  dc->particleCount_impl());
      //std::sort(dc->m_forces_ptr,   dc->m_forces_ptr   +  dc->particleCount_impl(), f);
      //std::sort(dc->m_vel_ptr,      dc->m_vel_ptr      +  dc->particleCount_impl(), f);
      //std::sort(dc->m_prevvel_ptr,  dc->m_prevvel_ptr  +  dc->particleCount_impl(), f);

      std::sort((float *) dc->m_pos_ptr,      ((float *) dc->m_pos_ptr      +  dc->particleCount_impl() * 4));
      std::sort((float *) dc->m_density_ptr,  ((float *) dc->m_density_ptr  +  dc->particleCount_impl()));
      std::sort((float *) dc->m_pressure_ptr, ((float *) dc->m_pressure_ptr +  dc->particleCount_impl()));
      std::sort((float *) dc->m_forces_ptr,   ((float *) dc->m_forces_ptr   +  dc->particleCount_impl() * 4));
      std::sort((float *) dc->m_vel_ptr,      ((float *) dc->m_vel_ptr      +  dc->particleCount_impl() * 4));
      std::sort((float *) dc->m_prevvel_ptr,  ((float *) dc->m_prevvel_ptr  +  dc->particleCount_impl() * 4));
    }

  protected:
    template <typename T>
    T *mapBuffer(boost::compute::command_queue & queue, const boost::compute::buffer & buf,
                 size_t n, cl_map_flags flags = CL_MAP_READ)
    {
      return (T *) queue.enqueue_map_buffer(buf, flags, 0, n * sizeof(T));
    }

    template <typename T>
    void unmapBuffer(boost::compute::command_queue & queue, const boost::compute::buffer & buf, T *ptr)
    {
      queue.enqueue_unmap_buffer(buf, (void *) ptr);
    }

    template <typename T>
    T *mapBuffer(boost::compute::command_queue & queue, const utils::ocl::GLBuffer & buf,
                 size_t n, cl_map_flags flags = CL_MAP_READ)
    {
      cl_int err = CL_SUCCESS;
      T *p = (T *) clEnqueueMapBuffer(queue, buf.getCLID(), CL_TRUE,
                                      flags, 0, n * sizeof(T),
                                      0, nullptr, nullptr, &err);
      if (p == nullptr)
      {
        WARNM("Failed to map GL Buffer: " <<
              boost::compute::opencl_error::to_string(err) <<
              "(" << err << ")");
      }
      return p;
    }

    template <typename T>
    void unmapBuffer(boost::compute::command_queue & queue, const utils::ocl::GLBuffer & buf, T *ptr)
    {
      cl_int err = clEnqueueUnmapMemObject(queue, buf.getCLID(), (void *) ptr, 0, nullptr, nullptr);
      if (err != CL_SUCCESS)
      {
        WARNM("Failed to unmap GL Buffer: " <<
              boost::compute::opencl_error::to_string(err) <<
              "(" << err << ")");
      }
    }
};

}

#endif // SPH_TEST_ADAPTER_H
