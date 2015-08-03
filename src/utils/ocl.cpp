#include "utils/ocl.h"
#include "utils/macros.h"
#include "utils/debug.h"

#include <boost/compute/platform.hpp>
#include <boost/compute/system.hpp>
#include <boost/compute/device.hpp>
#include <iomanip>


namespace utils {

namespace ocl {

///////////////////////////////////////////////////////////////////////////////
// Error handling

const char *errorToStr(cl_int err)
{
  switch (err)
  {
    case CL_SUCCESS: return "Success!";
    case CL_DEVICE_NOT_FOUND: return "Device not found.";
    case CL_DEVICE_NOT_AVAILABLE: return "Device not available";
    case CL_COMPILER_NOT_AVAILABLE: return "Compiler not available";
    case CL_MEM_OBJECT_ALLOCATION_FAILURE: return "Memory object allocation failure";
    case CL_OUT_OF_RESOURCES: return "Out of resources";
    case CL_OUT_OF_HOST_MEMORY: return "Out of host memory";
    case CL_PROFILING_INFO_NOT_AVAILABLE: return "Profiling information not available";
    case CL_MEM_COPY_OVERLAP: return "Memory copy overlap";
    case CL_IMAGE_FORMAT_MISMATCH: return "Image format mismatch";
    case CL_IMAGE_FORMAT_NOT_SUPPORTED: return "Image format not supported";
    case CL_BUILD_PROGRAM_FAILURE: return "Program build failure";
    case CL_MAP_FAILURE: return "Map failure";
    case CL_INVALID_VALUE: return "Invalid value";
    case CL_INVALID_DEVICE_TYPE: return "Invalid device type";
    case CL_INVALID_PLATFORM: return "Invalid platform";
    case CL_INVALID_DEVICE: return "Invalid device";
    case CL_INVALID_CONTEXT: return "Invalid context";
    case CL_INVALID_QUEUE_PROPERTIES: return "Invalid queue properties";
    case CL_INVALID_COMMAND_QUEUE: return "Invalid command queue";
    case CL_INVALID_HOST_PTR: return "Invalid host pointer";
    case CL_INVALID_MEM_OBJECT: return "Invalid memory object";
    case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR: return "Invalid image format descriptor";
    case CL_INVALID_IMAGE_SIZE: return "Invalid image size";
    case CL_INVALID_SAMPLER: return "Invalid sampler";
    case CL_INVALID_BINARY: return "Invalid binary";
    case CL_INVALID_BUILD_OPTIONS: return "Invalid build options";
    case CL_INVALID_PROGRAM: return "Invalid program";
    case CL_INVALID_PROGRAM_EXECUTABLE: return "Invalid program executable";
    case CL_INVALID_KERNEL_NAME: return "Invalid kernel name";
    case CL_INVALID_KERNEL_DEFINITION: return "Invalid kernel definition";
    case CL_INVALID_KERNEL: return "Invalid kernel";
    case CL_INVALID_ARG_INDEX: return "Invalid argument index";
    case CL_INVALID_ARG_VALUE: return "Invalid argument value";
    case CL_INVALID_ARG_SIZE: return "Invalid argument size";
    case CL_INVALID_KERNEL_ARGS: return "Invalid kernel arguments";
    case CL_INVALID_WORK_DIMENSION: return "Invalid work dimension";
    case CL_INVALID_WORK_GROUP_SIZE: return "Invalid work group size";
    case CL_INVALID_WORK_ITEM_SIZE: return "Invalid work item size";
    case CL_INVALID_GLOBAL_OFFSET: return "Invalid global offset";
    case CL_INVALID_EVENT_WAIT_LIST: return "Invalid event wait list";
    case CL_INVALID_EVENT: return "Invalid event";
    case CL_INVALID_OPERATION: return "Invalid operation";
    case CL_INVALID_GL_OBJECT: return "Invalid OpenGL object";
    case CL_INVALID_BUFFER_SIZE: return "Invalid buffer size";
    case CL_INVALID_MIP_LEVEL: return "Invalid mip-map level";
  }

  return "Unknown OpenCL error";
}


///////////////////////////////////////////////////////////////////////////////
// Information display

#define CASEL(x) case x: return #x

static const char *cacheTypeToString(cl_device_mem_cache_type type)
{
  switch (type)
  {
    CASEL(CL_NONE);
    CASEL(CL_READ_ONLY_CACHE);
    CASEL(CL_READ_WRITE_CACHE);
  }

  return "Unknown";
}


static const char *localMemTypeToString(cl_device_local_mem_type type)
{
  switch (type)
  {
    CASEL(CL_NONE);
    CASEL(CL_LOCAL);
    CASEL(CL_GLOBAL);
  }

  return "Unknown";
}

#undef CASEL

namespace {

template<typename T>
std::ostream & operator<<(std::ostream & os, const std::vector<T> & v)
{
  os << "[";
  if (!v.empty()) os << v[0];

  for (unsigned int i = 1; i < v.size(); ++i)
  {
    os << ", " << v[i];
  }

  return os << "]";
}

}


void printDeviceInfo(const boost::compute::device & dev)
{
  INFOM("====================================================================");

  INFOM("name                                      : " << dev.name());

  // misc information
  INFOM("--- misc info -----------------------------");
  INFOM("CL_DEVICE_ERROR_CORRECTION_SUPPORT        : " << dev.get_info<cl_bool>(CL_DEVICE_ERROR_CORRECTION_SUPPORT));
  INFOM("CL_DEVICE_HOST_UNIFIED_MEMORY             : " << dev.get_info<cl_bool>(CL_DEVICE_HOST_UNIFIED_MEMORY));
  INFOM("CL_DEVICE_MAX_COMPUTE_UNITS               : " << dev.get_info<cl_uint>(CL_DEVICE_MAX_COMPUTE_UNITS));
  INFOM("CL_DEVICE_PROFILING_TIMER_RESOLUTION (ns) : " << dev.get_info<size_t>(CL_DEVICE_PROFILING_TIMER_RESOLUTION));

  // global memory information
  INFOM("--- global memory info --------------------");
  INFOM("CL_DEVICE_GLOBAL_MEM_CACHE_SIZE (KB)      : " << dev.get_info<cl_ulong>(CL_DEVICE_GLOBAL_MEM_CACHE_SIZE) / 1024);
  INFOM("CL_DEVICE_GLOBAL_MEM_CACHE_TYPE           : " << cacheTypeToString(dev.get_info<cl_device_mem_cache_type>(CL_DEVICE_GLOBAL_MEM_CACHE_TYPE)));
  INFOM("CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE (B)   : " << dev.get_info<cl_uint>(CL_DEVICE_GLOBAL_MEM_CACHELINE_SIZE));
  INFOM("CL_DEVICE_GLOBAL_MEM_SIZE (MB)            : " << dev.get_info<cl_ulong>(CL_DEVICE_GLOBAL_MEM_SIZE) / 1024 / 1024);

  // local memory information
  INFOM("--- local memory info ---------------------");
  INFOM("CL_DEVICE_LOCAL_MEM_SIZE (KB)             : " << dev.get_info<cl_ulong>(CL_DEVICE_LOCAL_MEM_SIZE) / 1024);
  INFOM("CL_DEVICE_LOCAL_MEM_TYPE                  : " << localMemTypeToString(dev.get_info<cl_device_local_mem_type>(CL_DEVICE_LOCAL_MEM_TYPE)));

  // constant memory information
  INFOM("--- constant memory info ------------------");
  INFOM("CL_DEVICE_MAX_CONSTANT_ARGS               : " << dev.get_info<cl_uint>(CL_DEVICE_MAX_CONSTANT_ARGS));
  INFOM("CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE (KB)   : " << dev.get_info<cl_ulong>(CL_DEVICE_MAX_CONSTANT_BUFFER_SIZE) / 1024);

  // memory object information
  INFOM("--- memory object info --------------------");
  INFOM("CL_DEVICE_MAX_MEM_ALLOC_SIZE (MB)         : " << dev.get_info<cl_ulong>(CL_DEVICE_MAX_MEM_ALLOC_SIZE) / 1024 / 1024);

  // work group information
  INFOM("--- work group info -----------------------");
  INFOM("CL_DEVICE_MAX_WORK_GROUP_SIZE             : " << dev.get_info<size_t>(CL_DEVICE_MAX_WORK_GROUP_SIZE));
  INFOM("CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS        : " << dev.get_info<size_t>(CL_DEVICE_MAX_WORK_ITEM_DIMENSIONS));
  INFOM("CL_DEVICE_MAX_WORK_ITEM_SIZES             : " << dev.get_info<std::vector<size_t>>(CL_DEVICE_MAX_WORK_ITEM_SIZES));

  INFOM("====================================================================");
}


void printDeviceInfo(const std::vector<boost::compute::device> & devs)
{
  for (const boost::compute::device & dev : devs)
  {
    printDeviceInfo(dev);
  }
}


void printKernelInfo(boost::compute::kernel & kern)
{
  std::vector<boost::compute::device> devs(kern.get_program().get_devices());

  INFOM("====================================================================");
  INFOM("kernel : " << kern.name());

  for (boost::compute::device & dev : devs)
  {
    INFOM("device : " << dev.name());

    // This info is only supported for built-in kernels
    //INFOM("CL_KERNEL_GLOBAL_WORK_SIZE                   : " << kern.get_work_group_info<std::vector<size_t>>(dev, CL_KERNEL_GLOBAL_WORK_SIZE));

    // this info will return reasonable values only in case  __attribute__((reqd_work_group_size(X, Y, Z)))
    // has been defined for this kernel in OpenCL source
    //INFOM("CL_KERNEL_COMPILE_WORK_GROUP_SIZE            : " << kern.get_work_group_info<std::vector<size_t>>(dev, CL_KERNEL_COMPILE_WORK_GROUP_SIZE));

    INFOM("CL_KERNEL_WORK_GROUP_SIZE                    : " << kern.get_work_group_info<size_t>  (dev, CL_KERNEL_WORK_GROUP_SIZE));
    INFOM("CL_KERNEL_LOCAL_MEM_SIZE                     : " << kern.get_work_group_info<cl_ulong>(dev, CL_KERNEL_LOCAL_MEM_SIZE));
    INFOM("CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE : " << kern.get_work_group_info<size_t>  (dev, CL_KERNEL_PREFERRED_WORK_GROUP_SIZE_MULTIPLE));
    INFOM("CL_KERNEL_PRIVATE_MEM_SIZE                   : " << kern.get_work_group_info<cl_ulong>(dev, CL_KERNEL_PRIVATE_MEM_SIZE));
  }

  INFOM("====================================================================");
}


///////////////////////////////////////////////////////////////////////////////
// Device and platform management

bool selectGLDeviceAndPlatform(cl_device_id *device, cl_platform_id *platform)
{
  assert(device != nullptr);
  assert(platform != nullptr);

  /* first try to get the necessary extension */
  clGetGLContextInfoKHR_fn clGetGLContextInfoKHR = (clGetGLContextInfoKHR_fn) clGetExtensionFunctionAddress("clGetGLContextInfoKHR");
  if (clGetGLContextInfoKHR == nullptr)
  {
    ERRORM("clGetGLContextInfoKHR extension function not supproted.");
    return false;
  }

#if defined(FLUIDSIM_OS_MAC)
  cl_context_properties props[] = {
    CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE,
    (cl_context_properties) CGLGetShareGroup(CGLGetCurrentContext()),
    0
  };

  cl_int err = clGetGLContextInfoKHR(props,                                  // the OpenGL context
                                     CL_CURRENT_DEVICE_FOR_GL_CONTEXT_KHR,   // get the id of the device currently executing OpenGL
                                     sizeof(*device),
                                     device,
                                     nullptr);
  if (err != CL_SUCCESS)
  {
    ERRORM("Failed to retrieve the OpenCL id of the device executing OpenGL: " << errorToStr(err));
    return false;
  }

  /* get the platform associated with the device */
  err = clGetDeviceInfo(*device, CL_DEVICE_PLATFORM, sizeof(*platform), platform, nullptr);
  if (err != CL_SUCCESS)
  {
    ERRORM("Failed to retirieve platform id for the device executing OpenGL: " << errorToStr(err));
    return false;
  }
#else
  cl_context_properties props[] = {
#if defined(FLUIDSIM_OS_UNIX)
    CL_GL_CONTEXT_KHR, (cl_context_properties) glXGetCurrentContext(),
    CL_GLX_DISPLAY_KHR, (cl_context_properties) glXGetCurrentDisplay(),
    CL_CONTEXT_PLATFORM, (cl_context_properties) nullptr,
#elif defined(FLUIDSIM_OS_WIN)
    CL_GL_CONTEXT_KHR, (cl_context_properties) wglGetCurrentContext(),
    CL_WGL_HDC_KHR, (cl_context_properties) wglGetCurrentDC(),
    CL_CONTEXT_PLATFORM, (cl_context_properties) nullptr,
#else
# error "Unsupported OS platform"
#endif
    0
  };

  std::vector<boost::compute::platform> platform_list = boost::compute::system::platforms();

  for (const boost::compute::platform & p : platform_list)
  {
    WARNM("platform: " << p.name());
    props[5] = (cl_context_properties) p.id();
    cl_int err = clGetGLContextInfoKHR(props,                                 // the OpenGL context
                                       CL_CURRENT_DEVICE_FOR_GL_CONTEXT_KHR,  // get the id of the device currently executing OpenGL
                                       sizeof(*device),
                                       device,
                                       nullptr);
    if ((err == CL_SUCCESS) && (boost::compute::device(*device).type() == CL_DEVICE_TYPE_GPU))
    {
      *platform = (cl_platform_id) props[5];
      return true;
    }
    else
    {
      WARNM("clGetGLContextInfoKHR: " << errorToStr(err));
    }
  }
#endif

  return false;
}


bool selectPlatformAndDevice(cl_device_id *device, cl_platform_id *platform,
                             cl_device_type dev_type)
{
  assert(device != nullptr);
  assert(platform != nullptr);

  /* enumerate all available devices */
  std::vector<boost::compute::device> device_list = boost::compute::system::devices();

  /* select the one that is the most suitable according to given preferences */
  for (const boost::compute::device & d : device_list)
  {
    if (d.type() & dev_type)
    {
      *device = d.id();
      *platform = d.platform().id();
      return true;
    }
  }

  // no suitable device found
  return false;
}


bool findGPUDevice(boost::compute::device & dev)
{
  std::vector<boost::compute::device> devs(boost::compute::system::devices());

  for (const boost::compute::device & dev_ : devs)
  {
    if (dev_.type() == CL_DEVICE_TYPE_GPU)
    {
      INFOM("Using device: " << dev_.name());
      dev = dev_;
      return true;
    }
  }

  dev = boost::compute::system::default_device();

  WARNM("Could not find any GPU device. Using device: " << dev.name());

  return false;
}


///////////////////////////////////////////////////////////////////////////////
// OpenGL interoperability

bool initCLGLContext(boost::compute::context & ctx)
{
  INFOM("Initializing OpenCL/OpenGL shared context");

#if 1
  /* select appropriate device and platform */
  cl_platform_id platform = nullptr;
  cl_device_id device = nullptr;

  if ((!utils::ocl::selectGLDeviceAndPlatform(&device, &platform)) &&
      (!utils::ocl::selectPlatformAndDevice(&device, &platform)))
  {
    ERRORM("Failed to select an appropriate device or platform");
    return false;
  }

  /* setup context */
  cl_context_properties props[] = {
#if defined(FLUIDSIM_OS_MAC)
    CL_CONTEXT_PROPERTY_USE_CGL_SHAREGROUP_APPLE,
    (cl_context_properties) CGLGetShareGroup(CGLGetCurrentContext()),
#elif defined(FLUIDSIM_OS_UNIX)
    CL_GL_CONTEXT_KHR, (cl_context_properties) glXGetCurrentContext(),
    CL_GLX_DISPLAY_KHR, (cl_context_properties) glXGetCurrentDisplay(),
    CL_CONTEXT_PLATFORM, (cl_context_properties) platform,
#elif defined(FLUIDSIM_OS_WIN)
    CL_GL_CONTEXT_KHR, (cl_context_properties) wglGetCurrentContext(),
    CL_WGL_HDC_KHR, (cl_context_properties) wglGetCurrentDC(),
    CL_CONTEXT_PLATFORM, (cl_context_properties) platform,
#else
# error "Unsupported OS platform"
#endif
    0
  };

  /* create opencl context */
  cl_int err = CL_SUCCESS;
  cl_context ctx_ = clCreateContext(props, 1, &device, nullptr, nullptr, &err);
  if (err != CL_SUCCESS)
  {
    ERRORM("Failed to create OpenCL context: " << utils::ocl::errorToStr(err));
    return false;
  }

  ctx = boost::compute::context(ctx_, false);
  utils::ocl::printDeviceInfo(ctx.get_devices());
#else
  try
  {
    /* create opencl context */
    ctx = boost::compute::opengl_create_shared_context();
    utils::ocl::printDeviceInfo(ctx.get_devices());
  }
  catch (const boost::compute::opencl_error & e)
  {
    ERRORM(e.what());
    return false;
  }
  catch (const boost::compute::unsupported_extension_error & e)
  {
    ERRORM(e.what());
    return false;
  }
  catch (const std::exception & e)
  {
    ERRORM("An unexpected error occured during opencl context initialization: " << e.what());
    return false;
  }
#endif

  INFOM("Successfully initialized OpenCL context");
  INFOM("Using device   : " << ctx.get_device().name());
  INFOM("Using platform : " << ctx.get_device().platform().name());

  return true;
}


bool GLBuffer::bufferData(const GLvoid *data, GLsizeiptr size, AccessType at, GLenum usage)
{
  assert(m_ctx != nullptr);

  OGLF->glGetError();  // clear any previous errors

  OGLF->glBindBuffer(GL_ARRAY_BUFFER, m_vbo);
  OGLF->glBufferData(GL_ARRAY_BUFFER, size, data, usage);

  {
    GLenum err = OGLF->glGetError();
    if (err != GL_NO_ERROR)
    {
      ERRORM("Failed to buffer data: " << ogl::errorToStr(err));
      return false;
    }
  }

  OGLF->glBindBuffer(GL_ARRAY_BUFFER, 0);

  cl_int err = CL_SUCCESS;
  cl_mem mem = clCreateFromGLBuffer(m_ctx, at, m_vbo, &err);
  if (err != CL_SUCCESS)
  {
    ERRORM("Failed to create OpenCL buffer: " << ocl::errorToStr(err));
    return false;
  }

  clReleaseMemObject(m_mem);

  m_mem = mem;

  return true;
}


///////////////////////////////////////////////////////////////////////////////
// Performance counters and profiling

void CL_CALLBACK Timer::cl_callback(cl_event ev, cl_int /* status */, void *timer)
{
  Timer *t = (Timer *) timer;

  if (ev == t->m_start)
  {
    cl_ulong e = 0;
    cl_int err = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_END, sizeof(e), &e, nullptr);
    if (err != CL_SUCCESS)
    {
      WARNM("Timer::cl_callback: Failed to query start time: "
            << boost::compute::opencl_error::to_string(err)
            << " (" << err << ")");
    }
    t->m_start_time = e;
  }
  else if (ev == t->m_end)
  {
    cl_ulong s = 0;
    cl_int err = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_START, sizeof(s), &s, nullptr);
    if (err != CL_SUCCESS)
    {
      WARNM("Timer::cl_callback: Failed to query end time: "
            << boost::compute::opencl_error::to_string(err)
            << " (" << err << ")");
    }
    t->m_end_time = s;
  }
  else
  {
    WARNM("Timer::cl_callback: The event trapped is not the start neither the end event");
  }
}


void PerfStatsRecord::print(const std::string & name, std::ostream & os)
{
  os << "+----------------------------------------------------------------------+" << std::endl;
  os << "| " << std::setw(50) << std::left << name << " | " << std::setw(8) << std::right << m_count << " events |" << std::endl;
  os << "+----------------+-----------------+-----------------+-----------------+" << std::endl;
  os << "|    Type        |     min (ms)    |     max (ms)    |     avg (ms)    |" << std::endl;
  os << "+----------------+-----------------+-----------------+-----------------+" << std::endl;

  os << std::left;
  os << "| Queue time     | " << std::setw(15) << (m_queue_time.min() * 1.0e-6)  << " | "
                              << std::setw(15) << (m_queue_time.max() * 1.0e-6)  << " | "
                              << std::setw(15) << (m_queue_time.avg() * 1.0e-6)  << " |"
                              << std::endl;
  os << "| Submit time    | " << std::setw(15) << (m_submit_time.min() * 1.0e-6) << " | "
                              << std::setw(15) << (m_submit_time.max() * 1.0e-6) << " | "
                              << std::setw(15) << (m_submit_time.avg() * 1.0e-6) << " |"
                              << std::endl;
  os << "| Execution time | " << std::setw(15) << (m_exec_time.min() * 1.0e-6)   << " | "
                              << std::setw(15) << (m_exec_time.max() * 1.0e-6)   << " | "
                              << std::setw(15) << (m_exec_time.avg() * 1.0e-6)   << " |"
                              << std::endl;
  os << "+----------------+-----------------+-----------------+-----------------+" << std::endl;

  os << std::right;

  return;
}


void CL_CALLBACK PerfStats::cl_callback(cl_event ev, cl_int /* status */, void *stats_rec)
{
  cl_ulong time_queued = 0;
  cl_ulong time_submited = 0;
  cl_ulong time_started = 0;
  cl_ulong time_finished = 0;
  cl_int err = CL_SUCCESS;

  err = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_QUEUED, sizeof(time_queued), &time_queued, nullptr);
  if (err != CL_SUCCESS) WARNM("Failed to query queued time: " << ocl::errorToStr(err) << " (" << err << ")");

  err = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_SUBMIT, sizeof(time_submited), &time_submited, nullptr);
  if (err != CL_SUCCESS) WARNM("Failed to query submit time: " << ocl::errorToStr(err) << " (" << err << ")");

  err = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_START, sizeof(time_started), &time_started, nullptr);
  if (err != CL_SUCCESS) WARNM("Failed to query start time: " << ocl::errorToStr(err) << " (" << err << ")");

  err = clGetEventProfilingInfo(ev, CL_PROFILING_COMMAND_END, sizeof(time_finished), &time_finished, nullptr);
  if (err != CL_SUCCESS) WARNM("Failed to query end time: " << ocl::errorToStr(err) << " (" << err << ")");

  //LIBOCL_INFO("time_queued   : " << time_queued);
  //LIBOCL_INFO("time_submited : " << time_submited);
  //LIBOCL_INFO("time_started  : " << time_started);
  //LIBOCL_INFO("time_finished : " << time_finished);

  //LIBOCL_INFO("queued   : " << (time_submited - time_queued));
  //LIBOCL_INFO("submited : " << (time_started - time_submited));
  //LIBOCL_INFO("finished : " << (time_finished - time_started));

  static_cast<PerfStatsRecord *>(stats_rec)->add(time_submited - time_queued,
                                                 time_started - time_submited,
                                                 time_finished - time_started);

  //clReleaseEvent(ev);

  return;
}

}

} // End of namespace utilsd
