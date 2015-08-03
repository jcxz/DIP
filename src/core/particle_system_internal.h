/**
 * This header file is only meant to be included in *.cpp files
 * and contains some common macros and possibly in future other
 * functionality that can be shared across various particle system
 * classes, but that should not be exposed in public header files.
 */

//#define BOOST_COMPUTE_DEBUG_KERNEL_COMPILATION

#define DEBUG_SOURCE

#define deviceCall(kern, name, gws, lws, ...) \
  runKernel(gws, lws, (utils::ocl::KernelArgs(kern, name), __VA_ARGS__))

#define deviceCall2(kern, name, gws, ...) \
  runKernel(gws, (utils::ocl::KernelArgs(kern, name), __VA_ARGS__))
