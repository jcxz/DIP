/**
 * This header defines stuff that is shared between OpenCL and C++
 */

#ifndef SPH_OCL_COMMON_H
#define SPH_OCL_COMMON_H

#ifdef __OPENCL_VERSION__
// definition inside of OpenCL programs
typedef float4 tFloat4;
typedef float tFloat;
typedef ulong tULong;
typedef uint tUInt;
typedef uint4 tUInt4;
#else
#include <boost/compute/cl.hpp>
typedef cl_float4 tFloat4;
typedef cl_float tFloat;
typedef cl_ulong tULong;
typedef cl_uint tUInt;
typedef cl_uint4 tUInt4;
#endif

// A structure of simulation parameters
struct SimParams
{
  tFloat4 volumemin;
  tFloat4 volumemax;
  //tFloat4 gravitation;
  tULong seed;
  tFloat simscale;
  tFloat radius2;
  //tFloat polykern;
  tFloat mass_polykern;
  tFloat restdensity;
  tFloat intstiffness;
  tUInt numparticles;
  tFloat smoothradius;
  //tFloat viscosity;
  //tFloat lapkern;
  tFloat vterm;
  tFloat spikykern_half;
  tFloat slope;               // TODO: <unused>
  tFloat leftwave;            // TODO: <unused>
  tFloat rightwave;           // TODO: <unused>
  tFloat deltatime;
  tFloat limit;
  tFloat extstiffness;
  tFloat extdamping;
  tFloat radius;
  tFloat mass;
  tFloat time;                // TODO: <unused>
  tUInt flags;
  tFloat4 cell_size;          // TODO: if only uniformly sized cells were considered, then a single float is sufficient
  tUInt4 grid_size;
  //tFloat cell_size;
  //tUInt grid_size;
  tFloat4 grid_scale;         // TODO: <unused>
  tFloat4 top_face;
  tFloat4 bottom_face;
  tFloat4 front_face;
  tFloat4 back_face;
  tFloat4 left_face;
  tFloat4 right_face;
};

typedef struct SimParams tSimParams;

#define INVALID_UNIFORM_GRID_CELL_VALUE ((tUInt) (-1))

#endif // SPH_OCL_COMMON_H
