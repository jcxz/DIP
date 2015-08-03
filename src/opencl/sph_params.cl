#ifdef DEBUG_SOURCE
#include "sph_ocl_common.h"
#endif

#pragma OPENCL EXTENSION cl_amd_printf : enable


/**
 * A kernel to print simulation parameters for debugging purposes
 */
__kernel void sph_params_print(__constant const tSimParams *params)
{
  int i = get_global_id(0);
  if (i == 0)
  {
    printf((__constant char *) "=====================================================\n"
                               "volumemin      : %v4f\n", params->volumemin);
    printf((__constant char *) "volumemax      : %v4f\n", params->volumemax);
    printf((__constant char *) "seed           : %lu\n", params->seed);
    //printf((__constant char *) "gravitation    : %v4f\n", params->gravitation);
    printf((__constant char *) "simscale       : %f\n", params->simscale);
    printf((__constant char *) "radius2        : %f\n", params->radius2);
    //printf((__constant char *) "polykern       : %f\n", params->polykern);
    printf((__constant char *) "mass_polykern  : %f\n", params->mass_polykern);
    printf((__constant char *) "restdensity    : %f\n", params->restdensity);
    printf((__constant char *) "intstiffness   : %f\n", params->intstiffness);
    printf((__constant char *) "numparticles   : %u\n", params->numparticles);
    printf((__constant char *) "smoothradius   : %f\n", params->smoothradius);
    //printf((__constant char *) "viscosity      : %f\n", params->viscosity);
    //printf((__constant char *) "lapkern        : %f\n", params->lapkern);
    printf((__constant char *) "vterm          : %f\n", params->vterm);
    printf((__constant char *) "spikykern_half : %f\n", params->spikykern_half);
    printf((__constant char *) "slope          : %f\n", params->slope);
    printf((__constant char *) "leftwave       : %f\n", params->leftwave);
    printf((__constant char *) "rightwave      : %f\n", params->rightwave);
    printf((__constant char *) "deltatime      : %f\n", params->deltatime);
    printf((__constant char *) "limit          : %f\n", params->limit);
    printf((__constant char *) "extstiffness   : %f\n", params->extstiffness);
    printf((__constant char *) "extdamping     : %f\n", params->extdamping);
    printf((__constant char *) "radius         : %f\n", params->radius);
    printf((__constant char *) "mass           : %f\n", params->mass);
    printf((__constant char *) "time           : %f\n", params->time);
    printf((__constant char *) "flags          : %u\n", params->flags);
    printf((__constant char *) "cell_size      : %v4f\n", params->cell_size);
    printf((__constant char *) "grid_size      : %v4u\n", params->grid_size);
    printf((__constant char *) "grid_scale     : %v4u\n", params->grid_scale);
  }
}
