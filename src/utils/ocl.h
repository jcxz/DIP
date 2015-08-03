#ifndef OCL_LIB_H
#define OCL_LIB_H

//#define WIN32_LEAN_AND_MEAN
//#define NOMINMAX               // disable min, max macros from Windows headers

#include "utils/ogl.h"
#include "utils/stats.h"
#include "utils/debug.h"

#include <boost/compute/kernel.hpp>
#include <boost/compute/command_queue.hpp>
#include <boost/compute/event.hpp>
#include <boost/compute/interop/opengl/cl_gl.hpp>

#include <stdexcept>
#include <unordered_map>
#include <memory>


namespace utils {

namespace ocl {

///////////////////////////////////////////////////////////////////////////////
// Error handling

const char *errorToStr(cl_int err);

/** OpenCL exception class */
struct Exception : std::runtime_error
{
  explicit Exception(cl_int err) throw()
    : std::runtime_error(errorToStr(err))
  {
  }

  explicit Exception(const std::string & msg, cl_int err) throw()
    : std::runtime_error(msg + ": " + errorToStr(err))
  {
  }
};

///////////////////////////////////////////////////////////////////////////////
// Device and platform management

/**
 * This function will attempt to determine the device and platform id of the device
 * that is used by OpenGL. This function should be called only after a valid OpenGL
 * context is created.
 *
 * @param device a pointer to cl_device_id variable that will contain the selected device
 * @param platform a pointer to cl_platform_id variable that will contain the selected platform
 *
 * @return true on success, false in case no interoperable device is found.
 * Note however that false will be returned even in cases, when necessary extensions are
 * not supported by OpenCL implementation, so this does not necesarilly mean that no
 * device is capable of sharing between OpenGL and OpenCL is present.
 */
bool selectGLDeviceAndPlatform(cl_device_id *device, cl_platform_id *platform);

/**
 * A function to select appropriate platform and device id-s
 *
 * @param device a pointer to cl_device_id variable that will contain the selected device
 * @param platform a pointer to cl_platform_id variable that will contain the selected platform
 * @param dev_type the preferred type of device
 *
 * @return true on success, false otherwise
 */
bool selectPlatformAndDevice(cl_device_id *device, cl_platform_id *platform,
                             cl_device_type dev_type = CL_DEVICE_TYPE_GPU);

/**
 * @brief findGPUDevice will fill the dev parameter with the first OpenCL
 * capable gpu device and return true. In case no such device is present
 * in the system, then this function will return false and fill the dev
 * parameter the default OpenCL device.
 * @param dev an output parameter for the located gpu device.
 * @return true on success, false if no gpu device was found.
 */
bool findGPUDevice(boost::compute::device & dev);

///////////////////////////////////////////////////////////////////////////////
// Kernel management

/**
 * A class to initialize kernel arguments
 */
class KernelArgs
{
  public:
    template <typename T>
    struct LocalArg
    {
      size_t size;
      explicit LocalArg(size_t n) : size(n * sizeof(T)) { }
    };

  public:
    explicit KernelArgs(boost::compute::kernel & kernel, const char *kernel_name = "")
      : m_kernel(kernel.get()), m_next_arg(0),
        m_err(CL_SUCCESS), m_kernel_name(kernel_name)
    {
      assert(m_kernel != nullptr);
    }

    explicit KernelArgs(cl_kernel kernel, const char *kernel_name = "")
      : m_kernel(kernel), m_next_arg(0),
        m_err(CL_SUCCESS), m_kernel_name(kernel_name)
    {
      assert(m_kernel != nullptr);
    }

    template <typename T>
    KernelArgs & arg(T param)
    { return setArgNext(sizeof(T), &param); }

    template <typename T>
    KernelArgs & arg(T param, cl_uint index)
    { return setArgIndex(sizeof(T), &param, index); }

    template <typename T>
    KernelArgs & arg(LocalArg<T> param)
    { return setArgNext(param.size, nullptr); }

    template <typename T>
    KernelArgs & arg(LocalArg<T> param, cl_uint index)
    { return setArgIndex(param.size, nullptr, index); }

    template <typename T>
    KernelArgs & operator,(T param)
    { return setArgNext(sizeof(T), &param); }

    template <typename T>
    KernelArgs & operator,(LocalArg<T> param)
    { return setArgNext(param.size, nullptr); }

    bool hasError(void) const { return m_err != CL_SUCCESS; }
    cl_int error(void) const { return m_err; }
    const char *errorString(void) const { return errorToStr(m_err); }
    operator bool(void) const { return m_err == CL_SUCCESS; }

    cl_int argIndex(void) const { return m_next_arg - 1; }

    const char *kernelName(void) const { return m_kernel_name; }
    cl_kernel kernel(void) const { return m_kernel; }

  private:
    KernelArgs & setArgNext(size_t size, const void *value)
    {
      if (m_err != CL_SUCCESS) return *this;
      m_err = clSetKernelArg(m_kernel, m_next_arg++, size, value);
      if (m_err != CL_SUCCESS)
      {
        WARNM("Failed to set argument number "
              << (m_next_arg - 1) << "of " << m_kernel_name
              << "kernel : " << errorToStr(m_err));
      }
      return *this;
    }

    KernelArgs & setArgIndex(size_t size, const void *value, cl_uint index)
    {
      if (m_err != CL_SUCCESS) return *this;
      m_next_arg = index;
      m_err = clSetKernelArg(m_kernel, m_next_arg++, size, value);
      if (m_err != CL_SUCCESS)
      {
        WARNM("Failed to set argument number "
              << (m_next_arg - 1) << "of " << m_kernel_name
              << "kernel : " << errorToStr(m_err));
      }
      return *this;
    }

  private:
    cl_kernel m_kernel;          // the OpenCL kernel to be setup
    cl_uint m_next_arg;          // the number of the next argument that will be initialized
    cl_int m_err;                // whether any error occured while processing arguments
    const char *m_kernel_name;   // kernel's string name (for better diagnostic messages)
};

///////////////////////////////////////////////////////////////////////////////
// OpenGL interoperability

bool initCLGLContext(boost::compute::context & ctx);

/** Class representing a buffer that is shared between OpenCL and OpenGL */
class GLBuffer
{
  public:
    enum AccessType {
      READ_ONLY  = CL_MEM_READ_ONLY,
      WRITE_ONLY = CL_MEM_WRITE_ONLY,
      READ_WRITE = CL_MEM_READ_WRITE
    };

  public:
    explicit GLBuffer(cl_context ctx = nullptr)
      : m_vbo(0),
        m_mem(nullptr),
        m_ctx(ctx)
    {
      OGLF->glGenBuffers(1, &m_vbo);
      //GLenum err = OGLF->glGetError();
      //if (err != GL_NO_ERROR) throw Exception("Failed to construct GLBuffer", err);
    }

    GLBuffer(GLBuffer && other)
    {
      m_vbo = other.m_vbo;
      other.m_vbo = 0;
      m_mem = other.m_mem;
      other.m_mem = nullptr;
    }

    ~GLBuffer(void)
    {
      clReleaseMemObject(m_mem);
      OGLF->glDeleteBuffers(1, &m_vbo);
    }

    GLBuffer & operator=(GLBuffer && other)
    {
      swap(other);
      return *this;
    }

    GLuint getGLID(void) const { return m_vbo; }
    cl_mem getCLID(void) const { return m_mem; }
    void setCLContext(cl_context ctx) { m_ctx = ctx; }

    void swap(GLBuffer & other)
    {
      GLuint tmp_vbo = other.m_vbo;
      other.m_vbo = m_vbo;
      m_vbo = tmp_vbo;

      cl_mem tmp_mem = other.m_mem;
      other.m_mem = m_mem;
      m_mem = tmp_mem;
    }

    bool bufferData(const GLvoid *data,               // a pointer to data that will be stored in the buffer
                    GLsizeiptr size,                  // the size of the data that will be stored (in bytes)
                    AccessType at = READ_WRITE,       // the access type of OpenCL's memory object
                    GLenum usage = GL_DYNAMIC_DRAW);  // the way the buffer will be utilized (assume, that the contents
                                                      // of the buffer will be used in rendering, but that they will be
                                                      // changed often by an OpenCL kernel)

  private:
    GLBuffer(const GLBuffer & );
    GLBuffer & operator=(const GLBuffer & );

  private:
    GLuint m_vbo;      /// vertex buffer object of OpenGL
    cl_mem m_mem;      /// OpenCL memory object
    cl_context m_ctx;  /// OpenCL context to which the buffer is bound (not owned by this class)
};

/**
 * A class to synchronize OpenGL and OpenCL memory objects
 */
#if 1
class GLSyncHandler
{
  public:
#if 0
    GLSyncHandler(cl_command_queue queue, cl_mem buffer)
      : m_mem(buffer), m_queue(queue), m_num_buffers(1),
        m_buffers(&m_mem), m_err(CL_SUCCESS)
    {
      acquireGLObjects();
    }
#endif
    GLSyncHandler(cl_command_queue queue, cl_uint num_buffers, const cl_mem *buffers)
      : m_queue(queue), m_num_buffers(num_buffers),
        m_buffers(buffers), m_err(CL_SUCCESS)
    {
      acquireGLObjects();
    }

    ~GLSyncHandler(void)
    {
      /* unlock the vertex buffer object, so that OpenGL can continue using it */
      m_err = clEnqueueReleaseGLObjects(m_queue, m_num_buffers, m_buffers, 0, nullptr, nullptr);
      if (m_err != CL_SUCCESS)
      {
        WARNM("Failed to release an exclusive access to one of OpenGL's vertex buffer objects: "
              << errorToStr(m_err));
      }

      /* wait for OpenCL to finish processing */
      clFinish(m_queue);
    }

    bool hasError(void) const { return m_err != CL_SUCCESS; }
    operator bool(void) const { return m_err == CL_SUCCESS; }

  private:
    void acquireGLObjects(void)
    {
      assert(m_queue != nullptr);
      assert(m_buffers != nullptr);

      /* wait for OpenGL to finish rendering */
      glFinish();

      /* acquire access to the shared vertex buffer object */
      m_err = clEnqueueAcquireGLObjects(m_queue, m_num_buffers, m_buffers, 0, nullptr, nullptr);
      if (m_err != CL_SUCCESS)
      {
        WARNM("Failed to acquire an exclusive access to one of OpenGL's vertex buffer objects: "
              << errorToStr(m_err));
      }
    }

  private:
    // just disable copying for now
    GLSyncHandler(const GLSyncHandler & );
    GLSyncHandler & operator=(const GLSyncHandler & );

  private:
    //cl_mem m_mem;
    cl_command_queue m_queue;
    cl_uint m_num_buffers;
    const cl_mem *m_buffers;
    cl_int m_err;
};
#else
class GLSyncGuard
{
  public:
    explicit GLSyncGuard(QCLContextGL & ctx, const QCLMemoryObject & mem)
      : m_p_cl_ctx(&ctx), m_p_mem(&mem)
    { ctx.acquire(mem).waitForFinished(); }

    ~GLSyncGuard(void)
    { m_p_cl_ctx->release(*m_p_mem).waitForFinished(); }

  private:
    QCLContextGL *m_p_cl_ctx;
    const QCLMemoryObject *m_p_mem;
};
#endif

///////////////////////////////////////////////////////////////////////////////
// Memory management

// buffer writing
template <typename T>
inline bool writeBuffer(boost::compute::command_queue & queue,
                        const boost::compute::buffer & buf,
                        const T *host_data, const size_t n)
{
  try
  {
    queue.enqueue_write_buffer(buf, 0, n * sizeof(T), (const void *) host_data);
  }
  catch (boost::compute::opencl_error & err)
  {
    ERRORM("Failed to write to boost::compute::buffer: " <<
           err.error_string() << "(" << err.error_code() << ")");
    return false;
  }
  return true;
}

template <typename T>
inline bool writeBuffer(boost::compute::command_queue & queue,
                        const utils::ocl::GLBuffer & buf,
                        const T *host_data, const size_t n)
{
  cl_int err = clEnqueueWriteBuffer(queue, buf.getCLID(), CL_TRUE,
                                    0, n * sizeof(T), (const void *) host_data,
                                    0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    ERRORM("Failed to write to GL Buffer: " <<
          boost::compute::opencl_error::to_string(err) <<
          "(" << err << ")");
    return false;
  }
  return true;
}

// buffer clearing
template<typename T>
inline bool zeroBuffer(boost::compute::command_queue & queue,
                       const boost::compute::buffer & buf,
                       const size_t n)
{
  try
  {
    unsigned char pattern = 0;
    queue.enqueue_fill_buffer(buf, &pattern, sizeof(pattern), 0, n * sizeof(T));
  }
  catch (boost::compute::opencl_error & err)
  {
    ERRORM("Failed to zero out boost::compute::buffer: " <<
          err.error_string() << "(" << err.error_code() << ")");
    return false;
  }
  return true;
}

template <typename T>
inline bool zeroBuffer(boost::compute::command_queue & queue,
                       const utils::ocl::GLBuffer & buf,
                       const size_t n)
{
  unsigned char pattern = 0;
  cl_int err = clEnqueueFillBuffer(queue, buf.getCLID(),
                                   &pattern, sizeof(pattern),
                                   0, n * sizeof(T),
                                   0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    ERRORM("Failed to zero out GL Buffer: " <<
          boost::compute::opencl_error::to_string(err) <<
          "(" << err << ")");
    return false;
  }
  return true;
}

// buffer mapping
template <typename T>
inline T *mapBuffer(boost::compute::command_queue & queue,
                    const boost::compute::buffer & buf,
                    const size_t n,const cl_map_flags flags = CL_MAP_READ)
{
  return (T *) queue.enqueue_map_buffer(buf, flags, 0, n * sizeof(T));
}

template <typename T>
inline void unmapBuffer(boost::compute::command_queue & queue,
                        const boost::compute::buffer & buf, T *ptr)
{
  queue.enqueue_unmap_buffer(buf, (void *) ptr);
}

template <typename T>
inline T *mapBuffer(boost::compute::command_queue & queue,
                    const utils::ocl::GLBuffer & buf,
                    const size_t n, const cl_map_flags flags = CL_MAP_READ)
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
inline void unmapBuffer(boost::compute::command_queue & queue,
                        const utils::ocl::GLBuffer & buf, T *ptr)
{
  cl_int err = clEnqueueUnmapMemObject(queue, buf.getCLID(), (void *) ptr, 0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    WARNM("Failed to unmap GL Buffer: " <<
          boost::compute::opencl_error::to_string(err) <<
          "(" << err << ")");
  }
}

///////////////////////////////////////////////////////////////////////////////
// Information display

void printDeviceInfo(const boost::compute::device & dev);
void printDeviceInfo(const std::vector<boost::compute::device> & devs);
void printKernelInfo(boost::compute::kernel & kern);

template <typename T, typename BufType>
void dumpBuffer(boost::compute::command_queue & queue,
                const BufType & buf, std::ostream & os,
                size_t n, int max_print = 10)
{
  const T *p_data = mapBuffer<T>(queue, buf, n);
  if (p_data == nullptr)
  {
    WARNM("dumpBuffer: Failed to map buffer to host for dumping");
    return;
  }
  utils::debug::printArray1D(p_data, std::min(max_print, (int) n), os);
  unmapBuffer(queue, buf, p_data);
}

///////////////////////////////////////////////////////////////////////////////
// Performance counters and profiling

// A specialized OpenCL timer
// The class can be used to measure the time it took to execute an arbitrary number
// of OpenCL commands.
// The code works by injecting special marker commands to the provided command queue,
// which are then queried for the amount of time they spent in the queue.
class Timer
{
  public:
    explicit Timer(const boost::compute::command_queue & queue = boost::compute::command_queue())
      : m_queue(queue)
      , m_start()
      , m_end()
      , m_start_time(0.0f)
      , m_end_time(0.0f)
    { }

    void setCommandQueue(const boost::compute::command_queue & queue) { m_queue = queue; }

    void start(void)
    {
      m_start = m_queue.enqueue_marker();
      m_start.set_callback(cl_callback, CL_COMPLETE, (void *) this);
    }

    void stop(void)
    {
      m_end = m_queue.enqueue_marker();
      m_end.set_callback(cl_callback, CL_COMPLETE, (void *) this);
    }

    // Whether the timer is ready to be used, i.e whether the start and stop methods can be
    // called safely, i.e whether the previous measuring has already completed or the timer was not used yet
    bool isReady(void) const { return (m_end == nullptr) || (m_end.status() == CL_COMPLETE); }

    // This method is used when the timer is already finished, or whether it is still in progress
    // (That is both start() and stop() have already been called, but the measure block has not completed yet)
    bool isTimeValid(void) const { return (m_end != nullptr) && (m_end.status() == CL_COMPLETE); }

    void waitForValid(void) { if (m_end) m_end.wait(); }
    bool isStarted(void) const { return (m_start != nullptr) && (m_start.status() == CL_COMPLETE); }
    double elapsedMiliseconds(void) const { return (m_end_time - m_start_time) / 1000000.0f; }

  private:
    static void CL_CALLBACK cl_callback(cl_event ev, cl_int /* status */, void *timer);

  private:
    boost::compute::command_queue m_queue;
    boost::compute::event m_start;
    boost::compute::event m_end;
    double m_start_time;
    double m_end_time;
};


// A class to represent a single performance record
class PerfStatsRecord
{
  public:
    PerfStatsRecord(void)
      : m_count(0)
      , m_queue_time()
      , m_submit_time()
      , m_exec_time()
    {
    }

    void add(cl_ulong queue_time, cl_ulong submit_time, cl_ulong exec_time)
    {
      m_queue_time.add(queue_time, m_count);
      m_submit_time.add(submit_time, m_count);
      m_exec_time.add(exec_time, m_count);
      m_count++;
    }

    void print(const std::string & name, std::ostream & os);

  private:
    unsigned int m_count;                              /// the number of times this record has been recorded
    utils::stats::StatisticsBase<cl_ulong> m_queue_time;   /// aggregated statistical records for in queue time
    utils::stats::StatisticsBase<cl_ulong> m_submit_time;  /// aggregated statistical records for submit time
    utils::stats::StatisticsBase<cl_ulong> m_exec_time;    /// aggregated statistical records for execution time
};


class PerfStats
{
  private:
    struct EventProxy
    {
      PerfStatsRecord *m_stats_rec;
      cl_event m_event;

      EventProxy(PerfStatsRecord *stats_rec)
        : m_stats_rec(stats_rec), m_event(nullptr)
      { }

      ~EventProxy(void)
      {
        if (m_event != nullptr)
        {
          clSetEventCallback(m_event, CL_COMPLETE, cl_callback, (void *) m_stats_rec);
          clReleaseEvent(m_event);
        }
      }

      operator cl_event * (void) { return &m_event; }
    };

    struct EventProxyPtr
    {
      PerfStatsRecord *m_stats_rec;
      cl_event *m_event;

      EventProxyPtr(cl_event *event, PerfStatsRecord *stats_rec)
        : m_stats_rec(stats_rec), m_event(event)
      { }

      ~EventProxyPtr(void)
      {
        if ((m_event != nullptr) && (*m_event != nullptr))
        {
          clSetEventCallback(*m_event, CL_COMPLETE, cl_callback, (void *) m_stats_rec);
        }
      }

      operator cl_event * (void) { return m_event; }
    };

    struct EventProxyRef
    {
      PerfStatsRecord *m_stats_rec;
      boost::compute::event & m_event;

      EventProxyRef(boost::compute::event & event, PerfStatsRecord *stats_rec)
        : m_stats_rec(stats_rec), m_event(event)
      { }

      ~EventProxyRef(void)
      {
        m_event.set_callback(cl_callback, CL_COMPLETE, (void *) m_stats_rec);
      }

      operator cl_event * (void) { return &m_event.get(); }
    };

    static void CL_CALLBACK cl_callback(cl_event ev, cl_int /* status */, void *stats_rec);

  private:
    typedef std::unordered_map<std::string, std::unique_ptr<PerfStatsRecord>> tContainer;

  public:
    enum SortOption {
      SortByName,
      SortByType,
      SortByCount,
      SortByQueueTime,
      SortBySubmitTime,
      SortByExecTime,
      SortByTotalTime,
    };

    enum FilterFlags {
      FilterNone  = 0,
      FilterExec  = (1 << 0),
      FilterMem   = (1 << 1),
      FilterSync  = (1 << 2),
      FilterOther = (1 << 3),
      FilterAll   = FilterExec | FilterMem | FilterSync | FilterOther
    };

  public:
    PerfStats(void) : m_stats() { }

    void clear(void) { return m_stats.clear(); }

    //bool beginFrame(const std::string & name, boost::compute::command_queue & queue)
    //{
    //}

    //bool endFrame(const std::string & name, boost::compute::command_queue & queue)
    //{
    //}

    //void addTimer(const std::string & name, const Timer & t)
    //{
    //  PerfStatsRecord *rec = insertStat(name);
    //  rec->add(t.elapsedMiliseconds());
    //}

    void addEvent(boost::compute::event & ev, const std::string & name)
    {
      PerfStatsRecord *rec = insertStat(name);
      ev.set_callback(cl_callback, CL_COMPLETE, (void *) rec);
    }

    void addEvent(cl_event ev, const std::string & name)
    {
      PerfStatsRecord *rec = insertStat(name);
      if (ev != nullptr) clSetEventCallback(ev, CL_COMPLETE, cl_callback, (void *) rec);
    }

    EventProxyRef event(boost::compute::event & event, const std::string & name)
    { return EventProxyRef(event, insertStat(name)); }

    EventProxyPtr event(cl_event *event, const std::string & name)
    { return EventProxyPtr(event, insertStat(name)); }

    EventProxy event(const std::string & name)
    { return EventProxy(insertStat(name)); }

    void print(std::ostream & os, SortOption so = SortByExecTime, FilterFlags fflags = FilterAll) const;

    void printAsTable(std::ostream & os) const
    {
      for (auto & it : m_stats)
      {
        it.second->print(it.first, os);
        os << std::endl;
      }
    }

    friend std::ostream & operator<<(std::ostream & os, const PerfStats & stats)
    {
      //stats.print(os);
      stats.printAsTable(os);
      return os;
    }

  private:
    PerfStatsRecord *insertStat(const std::string & name)
    {
      tContainer::iterator it = m_stats.find(name);
      if (it == m_stats.end())
      {
        PerfStatsRecord *ptr = new PerfStatsRecord;
        m_stats.insert(std::make_pair(name, std::unique_ptr<PerfStatsRecord>(ptr)));
        return ptr;
      }
      else
      {
        return it->second.get();
      }
    }

  private:
    tContainer m_stats;
};

} // End of namespace ocl

} // End of namespace utils


// These operators have to be declared outside of the namespace to work properly
inline std::ostream & operator<<(std::ostream & os, const cl_float4 & v)
{
  return os << "[" << v.s[0] << "," << v.s[1] << "," << v.s[2] << "," << v.s[3] << "]";
}

inline std::ostream & operator<<(std::ostream & os, const cl_uint4 & v)
{
  return os << "[" << v.s[0] << "," << v.s[1] << "," << v.s[2] << "," << v.s[3] << "]";
}

#endif
