#pragma OPENCL EXTENSION cl_amd_printf : enable

#ifdef DEBUG_SOURCE
#include "sph_ocl_common.h"
#endif

inline float random(ulong seed, float min, float max)
{
  seed = seed * (seed * (seed * (seed * (seed + 1u) + 1u) + 1u) +1u);
  float t = (float) (seed) / (float) ((ulong) (-1));
  return mix(min, max, t);
}


/**
 * A kernel to initialize particle buffers with reasonable default values
 */
__kernel void sph_optimized_reset(__global         float4 *position,
                                  __global         float4 *velocity,
                                  __global         float4 *prev_velocity,
                                  __global         float *pressure,
                                  __global         float *density,
                                  __global         float4 *force,
                                  __constant const tSimParams *params)
{
  int gid = get_global_id(0);

#if 0
  position[gid] = (float4) (random(params->seed + gid + 1ul, params->volumemin.x, params->volumemax.x),
                            random(params->seed + gid + 2ul, params->volumemin.y, params->volumemax.y),
                            random(params->seed + gid + 3ul, params->volumemin.z, params->volumemax.z),
                            1.0f);
#else
  int dim_x = (params->volumemax.x - params->volumemin.x) * 0.75f;
  int dim_y = (params->volumemax.y - params->volumemin.y) * 0.75f;
  //int dim_x = 15;
  //int dim_y = 15;
  //position[gid] = (float4) (gid % dim_x, (gid / dim_x) % dim_y, gid / (dim_x * dim_y), 1.0f);
  position[gid] = (float4) (gid % dim_x           + params->volumemin.x * 0.75f,
                            (gid / dim_x) % dim_y + params->volumemin.y * 0.75f,
                            gid / (dim_x * dim_y),
                            1.0f);
#endif

  velocity[gid] = (float4) (0.0f);
  prev_velocity[gid] = (float4) (0.0f);
  pressure[gid] = 0;
  density[gid] = 0;
  force[gid] = (float4) (0.0f);

  //printf("sph_reset: seed == %u\n", params->seed);
  //printf("sph_reset: position[%u] == [%v4f]\n", gid, position[gid]);
  //printf("sph_reset: position[%u] == [%v4f]\n", gid, velocity[gid]);
  //printf("sph_reset: position[%u] == [%v4f]\n", gid, prev_velocity[gid]);
  //printf("sph_reset: position[%u] == %f\n", gid, pressure[gid]);
  //printf("sph_reset: position[%u] == %f\n", gid, density[gid]);
  //printf("sph_reset: position[%u] == [%v4f]\n", gid, force[gid]);
}


__kernel void sph_optimized_compute_pressure(__global   const float4* pos,
                                             __global         float*  density,
                                             __global         float*  pressure,
                                             __constant const tSimParams *params)
{
  int i = get_global_id(0);

  float sum = 0.0f;
  float4 pos_reg = pos[i];

  for (int j = 0; j < params->numparticles; ++j)
  {
    float4 d = (pos_reg - pos[j]) * params->simscale;

    // note for myself:
    // the dot product of float4 is defined as x*x + y*y + z*z + w*w,
    // but since we set w to 1.0 in reset kernel and by subtracting
    // pos[i] - pos[j] we get 0,
    // this should not be a problem, in case
    // this is not true, then this has to be changed
    // to dsq = dot(d.xyz, d.xyz) or something similar
    float sqr = dot(d, d);
    float c = (params->radius2 - sqr) * ((params->radius2 > sqr) & (i != j));
    sum += c * c * c;
  }

  // polykern konstanta
  float ro = sum * params->mass_polykern; //mass * polykern;

  //                          kludova hustota vody pri 20 C     plynova konstanta (ci sa to bude chovat skor ako plyn)
  pressure[i] = (ro - params->restdensity) * params->intstiffness;   // vzorec
  density[i] = 1.0f / ro;
}


__kernel void sph_optimized_compute_force(__global   const float4* pos,
                                          __global   const float*  density,
                                          __global   const float*  pressure,
                                          __global         float4* forces,
                                          __global   const float4* vel,
                                          __constant const tSimParams *params)
{
#if 0
  // VERZIA KDE JE ODSTRANENY if (i != j)
  int i = get_global_id(0);

  float4 force = (float4) (0.0f, 0.0f, 0.0f, 0.0f);

  //float vterm = params->lapkern * params->viscosity;

  for (int j = 0; j < params->numparticles; ++j)
  {
    float4 d = (pos[i] - pos[j]) * params->simscale;
    float sqr = dot(d, d);

    if ((params->radius2 > sqr) & (i != j))
    {
      float r = sqrt(sqr);
      float c = (params->smoothradius - r);
      //float pterm = -0.5f * c * params->spikeykern * (pressure[i] + pressure[j]) / r;
      float pterm = c * params->spikykern_half * (pressure[i] + pressure[j]) / r;
      float dterm = c * density[i] * density[j];

      force += (pterm * d + params->vterm * (vel[j] - vel[i])) * dterm;
    }
  }

  forces[i] = force;
#elif 1
  // VERZIA KDE JE ODSTRANENY if (i != j) A pos[i] JE ULOZENE V REGISTRI
  int i = get_global_id(0);

  float4 force = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
  float4 pos_reg = pos[i];

  //float vterm = params->lapkern * params->viscosity;

  for (int j = 0; j < params->numparticles; ++j)
  {
    float4 d = (pos_reg - pos[j]) * params->simscale;
    float sqr = dot(d, d);

    if ((params->radius2 > sqr) & (i != j))
    {
      float r = sqrt(sqr);
      float c = (params->smoothradius - r);
      //float pterm = -0.5f * c * spikeykern * (pressure[i] + pressure[j]) / r;
      float pterm = c * params->spikykern_half * (pressure[i] + pressure[j]) / r;
      float dterm = c * density[i] * density[j];

      force += (pterm * d + params->vterm * (vel[j] - vel[i])) * dterm;
    }
  }

  forces[i] = force;
#elif 0
  // VERZIA KDE SU UPLNE ODSTRANENE VSETKY PODMIENKY, ALE HODNOTY PRE i-tu CASTICU SU VZDY ZNOVU-NACITANE Z PAMATE
  int i = get_global_id(0);

  float4 force = (float4) (0.0f, 0.0f, 0.0f, 0.0f);

  //float vterm = params->lapkern * params->viscosity;

  for (int j = 0; j < params->numparticles; ++j)
  {
    float4 d = (pos[i] - pos[j]) * params->simscale;
    float sqr = dot(d, d);

    int pred = (params->radius2 > sqr) & (i != j);

    {
      float r = sqrt(sqr);
      float c = (params->smoothradius - r);
      //float pterm = -0.5f * c * spikeykern * (pressure[i] + pressure[j]) / r;
      float pterm = c * params->spikykern_half * (pressure[i] + pressure[j]) / r;
      float dterm = c * density[i] * density[j];

      float dterm2 = select(0.0f, dterm, pred);
      float pterm2 = select(0.0f, pterm, pred);

      force += (pterm2 * d + params->vterm * (vel[j] - vel[i])) * dterm2;
    }
  }

  forces[i] = force;
#elif 0
  // VERZIA KDE SU UPLNE ODSTRANENE VSETKY PODMIENKY A HODNOTY Z PAMATE PRE i-tu CASTICU SU V REGISTROCH
  int i = get_global_id(0);

  float4 force = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
  float4 pos_reg = pos[i];
  float4 vel_reg = vel[i];
  float pressure_reg = pressure[i];
  float density_reg = density[i];

  //float vterm = params->lapkern * params->viscosity;

  for (int j = 0; j < params->numparticles; ++j)
  {
    float4 d = (pos_reg - pos[j]) * params->simscale;
    float sqr = dot(d, d);

    int pred = (params->radius2 > sqr) & (i != j);

    float r = sqrt(sqr);
    float c = (params->smoothradius - r);
    //float pterm = -0.5f * c * spikeykern * (pressure[i] + pressure[j]) / r;
    float pterm = c * params->spikykern_half * (pressure_reg + pressure[j]) / r;
    float dterm = c * density_reg * density[j];

    float dterm2 = select(0.0f, dterm, pred);
    float pterm2 = select(0.0f, pterm, pred);

    force += (pterm2 * d + params->vterm * (vel[j] - vel_reg)) * dterm2;
  }

  forces[i] = force;
#else
  // VERZIA KDE SU UPLNE ODSTRANENE VSETKY PODMIENKY, HODNOTY Z PAMATE PRE i-tu CASTICU SU V REGISTROCH, JE OPTIMALIZOVANA ARITMETIKA
  int i = get_global_id(0);

  float4 force = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
  float4 pos_reg = pos[i];
  float4 vel_reg = vel[i];
  float pressure_reg = pressure[i];
  float density_reg = density[i];

  //float vterm = params->lapkern * params->viscosity;

  for (int j = 0; j < params->numparticles; ++j)
  {
    float4 d = (pos_reg - pos[j]) * params->simscale;
    float sqr = dot(d, d);

    int pred = (params->radius2 > sqr) & (i != j);

    float r = native_rsqrt(sqr);
    float c = (params->smoothradius - sqr * r);
    //float pterm = -0.5f * c * spikeykern * (pressure[i] + pressure[j]) / r;
    float pterm = c * params->spikykern_half * (pressure_reg + pressure[j]) * r;
    float dterm = c * density_reg * density[j];

    float dterm2 = select(0.0f, dterm, pred);
    float pterm2 = select(0.0f, pterm, pred);

    force += (pterm2 * d + params->vterm * (vel[j] - vel_reg)) * dterm2;
  }

  forces[i] = force;
#endif
}


/**
 * The code has been written according to http://www.rchoetzlein.com/eng/graphics/fluids.htm
 * and http://joeyfladderak.com/portfolio-items/sph-fluid-simulation/
 */

#define DRAIN_MASK    (1 << 0)
#define WAVE_MASK     (1 << 1)
#define FOUNTAIN_MASK (1 << 2)

#if 0
__kernel void sph_optimized_compute_step(__global         float4* position,
                                         __global   const float4* forces,
                                         __global         float4* velocity,
                                         __global         float4* prevvelocity,
                                         __constant const tSimParams *params)
{
  int i = get_global_id(0);

  float4 norm = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
  float diff;

  float4 pos = position[i];
  float4 prevvel = prevvelocity[i];
  float4 accel = forces[i] * params->mass;

  float speed = accel.x * accel.x + accel.y * accel.y + accel.z * accel.z;
  if (speed > (params->limit * params->limit))
  {
    accel *= params->limit / sqrt(speed);
  }

  /* Y-axis wall */
  diff = 2.0f * params->radius - (pos.y - params->volumemin.y - (pos.x - params->volumemin.x) * params->slope) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (-params->slope, 1.0f - params->slope, 0.0f, 0.0f);
    //norm.x = -params->slope;
    //norm.z = 0.0f;
    //norm.y = 1.0f - params->slope;  // = float4(-0, 0, 1.0 - 0, 0);
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * params->radius - (params->volumemax.y - pos.y) * params->simscale;
  if (diff > 0.0001f)
  {
    //float4 norm = (float4) (0.0f, 0.0f, -1.0f, 0.0f);
    float4 norm = (float4) (0.0f, -1.0f, 0.0f, 0.0f);
    //norm.x = 0.0f;
    //norm.z = -1.0f;
    //norm.y = 0.0f;  //float4( 0, 0, -1, 0 );
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* X-axis walls */
  diff = 2.0f * params->radius - (pos.x - params->volumemin.x + (sin(params->time * 10.0f) - 1.0f + (pos.y * 0.025f) * 0.25f) * params->leftwave) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (1.0f, 0.0f, 0.0f, 0.0f);
    //norm.x = 1.0f;
    //norm.y = 0.0f;
    //norm.z = 0.0f;  //float4( 1.0, 0, 0, 0 );
    float adj = (params->leftwave + 1.0f) * params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * params->radius - (params->volumemax.x - pos.x + (sin(params->time * 10.0f) - 1.0f) * params->rightwave) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (-1.0f, 0.0f, 0.0f, 0.0f);
    //norm.x = -1.0f;
    //norm.y = 0.0f;
    //norm.z = 0.0f;  //float4( -1, 0, 0, 0 );
    float adj = (params->rightwave + 1.0f) * params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* Z-axis walls */
  diff = 2.0f * params->radius - (pos.z - params->volumemin.z) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (0.0f, 0.0f, 1.0f, 0.0f);
    //norm.x = 0.0f;
    //norm.z = 1.0f;
    //norm.y = 0.0f;  //float4( 0, 1, 0, 0 );
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * params->radius - (params->volumemax.z - pos.z) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (0.0f, 0.0f, -1.0f, 0.0f);
    //norm.x = 0.0f;
    //norm.z = -1.0f;
    //norm.y = 0.0f;  //float4( 0, -1, 0, 0 );
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* drain */
#if 1
  if (params->flags & DRAIN_MASK)
  {
    float dx = (0 - pos.x) * params->simscale;            // dist in cm
    float dy = (params->volumemin.y - pos.y) * params->simscale;
    float dz = (0 - pos.z) * params->simscale;
    float dsq = (dx * dx + dy * dy + dz * dz);
    //if (0.0001f > dsq)
    if (0.005f > dsq)
    {
      pos.x = params->volumemin.x * 0.5f;
      pos.y = params->volumemax.y * 0.5f;
      pos.z = params->volumemin.z;
      accel.z += 50;
      accel.x += 1;
      accel.y += 20;
    }
  }
#else
  if (params->flags & DRAIN_MASK)
  {
    diff = 2 * params->radius - (pos.z - params->volumemin.z - 15) * params->simscale;
    if (diff < 2 * params->radius && diff > 0.0001f && (fabs(pos.x) > 3 || fabs(pos.y) > 3))
    {
      norm = (float4) (0, 0, 1, 0);
      //adj = stiff * diff - damp * dot(norm, prevvelocity[i]);
      float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvelocity[i]);
      accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
    }
  }
#endif

  /* wave */
  if (params->flags & WAVE_MASK)
  {
    if (pos.x < 0) accel.x += 20;
  }

  /* fountain */
  if (params->flags & FOUNTAIN_MASK)
  {
    float dx = (0 - pos.x) * params->simscale;            // dist in cm
    float dy = (params->volumemin.y - pos.y) * params->simscale;
    float dz = (0 - pos.z) * params->simscale;
    float dsq = (dx * dx + dy * dy + dz * dz);
    //if (0.0005f > dsq)
    if (0.005f > dsq)
    {
      accel.y += 140;
    }
  }

  accel.y += -9.8f;  //float4(0, -9.8, 0, 0);
  //accel += gravitation;

  // Leapfrog Integration ----------------------------
  float4 vel = velocity[i];
  float4 vnext = accel * params->deltatime + vel;  // v(t+1/2) = v(t-1/2) + a(t) dt
  prevvel = (vel + vnext) * 0.5f;                  // v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5     used to compute forces later
  vel = vnext;
  vnext *= params->deltatime / params->simscale;
  pos += vnext;                                    // p(t+1) = p(t) + v(t+1/2) dt

  velocity[i] = vel;
  prevvelocity[i] = prevvel;
  position[i] = pos;
}
#endif







///// TOTO JE ZJEDNODUSENA VERZIA, KTORA PODPORUJE PRELIEVANIE TEKUTINY (SLOSHING)
#if 1
__kernel void sph_optimized_compute_step(__global         float4 *position,
                                         __global   const float4 *forces,
                                         __global         float4 *velocity,
                                         __global         float4 *prevvelocity,
                                         __constant const tSimParams *params)
{
  int i = get_global_id(0);

  float4 norm = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
  float diff;

  float4 pos = position[i];
  float4 prevvel = prevvelocity[i];
  float4 accel = forces[i] * params->mass;

  float speed = accel.x * accel.x + accel.y * accel.y + accel.z * accel.z;
  if (speed > (params->limit * params->limit))
  {
    accel *= params->limit / sqrt(speed);
  }

  // Vypocty kolizii so stenami obaloveho telesa tekutiny vychadzaju
  // zo vzorca na vypocet vzdialenosti bodu od roviny, ktory ma tvar:
  //
  // d = (A * x + B * y + C * z - D) / sqrt(A * A + B * B + C * C)
  //
  // kazdu stenu obaloveho objemu mam reprezentovanu normalou danej steny
  // (roviny v ktorej lezi) a vzdialenostou D, co odpoveda premennym A, B, C, D
  // vo vyssie uvedenom vzorci.
  // V pamati mam kazdu roviny reprezentovanu vektorom 4 cisel.
  // Zlozky vektora x, y, z reprezentuju normalu a zlozka w reprezentuje posunutie D.
  // Normala je normalizovana na dlzku 1, takze v kode mi staci vzdialenost pocitat
  // len takto:
  //
  // d = A * x + B * y + C * z - D
  //
  // co je vlastne vektorovy sucin bodu, ktory prave vysetrujem a normaly danej steny
  // minus posunutie D

  const float lala = 1.0f; //0.1f; //1.0f; //0.01f;
  const float collision_treshold = 0.0001f;
  //const float diff_max = 10.0f; //0.01f; //2.0f * params->radius; //0.03f;

  /* Y-axis wall */
  // dolna stena (min_y)
  diff = dot(pos.xyz, params->bottom_face.xyz) - params->bottom_face.w;
  diff = 2.0f * params->radius - diff * params->simscale;
  if (diff > collision_treshold)
  {
    //if (diff > diff_max) diff = diff_max;
    //float4 norm = (float4) (0.0f, 1.0f, 0.0f, 0.0f);
    float4 norm = params->bottom_face; norm.w = 0.0f;
    float adj = params->extstiffness * diff * lala - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  // horna stena (max_y)
  diff = dot(pos.xyz, params->top_face.xyz) - params->top_face.w;
  diff = 2.0f * params->radius - diff * params->simscale;
  if (diff > collision_treshold)
  {
    //if (diff > diff_max) diff = diff_max;
    //float4 norm = (float4) (0.0f, -1.0f, 0.0f, 0.0f);
    float4 norm = params->top_face; norm.w = 0.0f;
    float adj = params->extstiffness * diff * lala - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* X-axis walls */
  // lava stena (min_x)
  diff = dot(pos.xyz, params->left_face.xyz) - params->left_face.w;
  diff = 2.0f * params->radius - diff * params->simscale;
  if (diff > collision_treshold)
  {
    //if (diff > diff_max) diff = diff_max;
    //float4 norm = (float4) (1.0f, 0.0f, 0.0f, 0.0f);
    float4 norm = params->left_face; norm.w = 0.0f;
    float adj = params->extstiffness * diff * lala - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  // prava stena (max_x)
  diff = dot(pos.xyz, params->right_face.xyz) - params->right_face.w;
  diff = 2.0f * params->radius - diff * params->simscale;
  if (diff > collision_treshold)
  {
    //if (diff > diff_max) diff = diff_max;
    //float4 norm = (float4) (-1.0f, 0.0f, 0.0f, 0.0f);
    float4 norm = params->right_face; norm.w = 0.0f;
    float adj = params->extstiffness * diff * lala - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* Z-axis walls */
  // zadna stena (min_z)
  diff = dot(pos.xyz, params->back_face.xyz) - params->back_face.w;
  diff = 2.0f * params->radius - diff * params->simscale;
  if (diff > collision_treshold)
  {
    //if (diff > diff_max) diff = diff_max;
    //float4 norm = (float4) (0.0f, 0.0f, 1.0f, 0.0f);
    float4 norm = params->back_face; norm.w = 0.0f;
    float adj = params->extstiffness * diff * lala - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  // predna stena (max_z)
  diff = dot(pos.xyz, params->front_face.xyz) - params->front_face.w;
  diff = 2.0f * params->radius - diff * params->simscale;
  if (diff > collision_treshold)
  {
    //if (diff > diff_max) diff = diff_max;
    //float4 norm = (float4) (0.0f, 0.0f, -1.0f, 0.0f);
    float4 norm = params->front_face; norm.w = 0.0f;
    float adj = params->extstiffness * diff * lala - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* drain */
#if 1
  if (params->flags & DRAIN_MASK)
  {
    float dx = (0.0f - pos.x) * params->simscale;            // dist in cm
    float dy = (params->volumemin.y - pos.y) * params->simscale;
    float dz = (0.0f - pos.z) * params->simscale;
    float dsq = (dx * dx + dy * dy + dz * dz);
    //if (0.0001f > dsq)
    if (0.005f > dsq)
    {
      diff = 2.0f * params->radius - (params->volumemax.z - pos.z) * params->simscale;

      pos.x = params->volumemax.x; //pos.x = 0.0f; //volumemin.x * 0.5f;
      pos.y = 0.0f; //volumemax.y * 0.5f;
      pos.z = 0.0f; //volumemin.z;

      //float4 norm = (float4) (0.0f, 0.5f, 0.5f, 0.0f);
      //float4 norm = (float4) (0.0f, sin(3.14f / 0.25f), cos(3.14f / 0.25f), 0.0f);
      float4 norm = (float4) (0.5f, -0.25f, 0.0f, 0.0f);
      norm = normalize(norm);
      float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
      accel += adj * norm;// * 0.75f;

      //accel.z =  0.4f;
      //accel.x =  0.4f;
      //accel.y = -0.4f;
    }
  }
#elif 0
  if (params->flags & DRAIN_MASK)
  {
    float dx = (0.0f - pos.x) * params->simscale;            // dist in cm
    float dy = (params->volumemin.y - pos.y) * params->simscale;
    float dz = (0.0f - pos.z) * params->simscale;
    float dsq = (dx * dx + dy * dy + dz * dz);
    if (0.0001f > dsq)
    {
      pos.x = params->volumemin.x * 0.5f;
      pos.y = params->volumemax.y * 0.5f;
      pos.z = params->volumemin.z;
      accel.z += 50;
      accel.x += 1;
      accel.y += 20;
    }
  }
#else
  if (params->flags & DRAIN_MASK)
  {
    diff = 2 * params->radius - (pos.z - params->volumemin.z - 15.0f) * params->simscale;
    if (diff < 2 * params->radius && diff > 0.0001f && (fabs(pos.x) > 3 || fabs(pos.y) > 3))
    {
      norm = (float4) (0, 0, 1, 0);
      //adj = stiff * diff - damp * dot(norm, prevvelocity[i]);
      float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvelocity[i]);
      accel.x += adj * norm.x;
      accel.y += adj * norm.y;
      accel.z += adj * norm.z;
    }
  }
#endif

  /* wave */
  if (params->flags & WAVE_MASK)
  {
    if (pos.x < 0.0f) accel.x += 20.0f;
  }

  /* fountain */
  if (params->flags & FOUNTAIN_MASK)
  {
    float dx = (0.0f - pos.x) * params->simscale;            // dist in cm
    float dy = (params->volumemin.y - pos.y) * params->simscale;
    float dz = (0.0f - pos.z) * params->simscale;
    float dsq = (dx * dx + dy * dy + dz * dz);
    //if (0.0005f > dsq)
    if (0.005f > dsq)
    {
      accel.y += 140.0f;
    }
  }

  accel.y += -9.8f;
  //accel += gravitation;

  //float speed = accel.x * accel.x + accel.y * accel.y + accel.z * accel.z;
  //if (speed > (params->limit * params->limit))
  //{
  //  accel *= params->limit / sqrt(speed);
  //}

  // Leapfrog Integration ----------------------------
#if 0
  float4 vel = velocity[i];
  float4 vnext = accel * params->deltatime + vel;  // v(t+1/2) = v(t-1/2) + a(t) dt
  prevvel = (vel + vnext) * 0.5f;          // v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5     used to compute forces later
  vel = vnext;
  vnext *= params->deltatime / params->simscale;
  pos += vnext;                       // p(t+1) = p(t) + v(t+1/2) dt

  velocity[i] = vel;
  prevvelocity[i] = prevvel;
  position[i] = pos;
#else
  float4 vel = velocity[i];
  float4 vnext = accel * params->deltatime + vel;                      // v(t+1/2) = v(t-1/2) + a(t) dt
  velocity[i] = vnext;
  prevvelocity[i] = (vel + vnext) * 0.5f;                              // v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5     used to compute forces later
  //prevvelocity[i] = vel;   // toto funguje tiez, ale ma to horsie matematicke vlastnosti
  position[i] = pos + vnext * (params->deltatime / params->simscale);  // p(t+1) = p(t) + v(t+1/2) dt
#endif
}
#endif
























///// TOTO JE NOVA ZJEDNODUSENA VERZIA
#if 0
__kernel void sph_optimized_compute_step(__global       float4* position,
                                         __global const float4* forces,
                                         __global       float4* velocity,
                                         __global       float4* prevvelocity,
                                         float slope,
                                         float leftwave,
                                         float rightwave,
                                         float deltatime,
                                         float limit,
                                         float extstiffness,
                                         float extdamping,
                                         float radius,
                                         float4 volumemin,
                                         float4 volumemax,
                                         float simscale,
                                         float mass,
                                         float time,
                                         uint flags)
                                         //float4 gravitation)
{
  int i = get_global_id(0);

  float4 norm = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
  float diff;

  float4 pos = position[i];
  float4 prevvel = prevvelocity[i];
  float4 accel = forces[i] * mass;

  float speed = accel.x * accel.x + accel.y * accel.y + accel.z * accel.z;
  if (speed > (limit * limit))
  {
    accel *= limit / sqrt(speed);
  }

  /* Y-axis wall */
  //diff = 2.0f * radius - (pos.y - volumemin.y - (pos.x - volumemin.x) * slope) * simscale;
  diff = 2.0f * radius - (pos.y - volumemin.y) * simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (-slope, 1.0f - slope, 0.0f, 0.0f);
    float adj = extstiffness * diff - extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * radius - (volumemax.y - pos.y) * simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (0.0f, 0.0f, -1.0f, 0.0f);
    float adj = extstiffness * diff - extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* X-axis walls */
  //diff = 2.0f * radius - (pos.x - volumemin.x + (sin(time * 10.0f) - 1.0f + (pos.y * 0.025f) * 0.25f) * leftwave) * simscale;
  diff = 2.0f * radius - (pos.x - volumemin.x) * simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (1.0f, 0.0f, 0.0f, 0.0f);
    float adj = (leftwave + 1.0f) * extstiffness * diff - extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  //diff = 2.0f * radius - (volumemax.x - pos.x + (sin(time * 10.0f) - 1.0f) * rightwave) * simscale;
  diff = 2.0f * radius - (volumemax.x - pos.x) * simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (-1.0f, 0.0f, 0.0f, 0.0f);
    float adj = (rightwave + 1.0f) * extstiffness * diff - extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* Z-axis walls */
  diff = 2.0f * radius - (pos.z - volumemin.z) * simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (0.0f, 0.0f, 1.0f, 0.0f);
    float adj = extstiffness * diff - extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * radius - (volumemax.z - pos.z) * simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (0.0f, 0.0f, -1.0f, 0.0f);
    float adj = extstiffness * diff - extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* drain */
#if 1
  if (flags & DRAIN_MASK)
  {
    float dx = (0 - pos.x) * simscale;            // dist in cm
    float dy = (volumemin.y - pos.y) * simscale;
    float dz = (0 - pos.z) * simscale;
    float dsq = (dx * dx + dy * dy + dz * dz);
    if (0.0001f > dsq)
    {
      diff = 2.0f * radius - (volumemax.z - pos.z) * simscale;

      pos.x = volumemax.x; //pos.x = 0.0f; //volumemin.x * 0.5f;
      pos.y = 0.0f; //volumemax.y * 0.5f;
      pos.z = 0.0f; //volumemin.z;

      //float4 norm = (float4) (0.0f, 0.5f, 0.5f, 0.0f);
      //float4 norm = (float4) (0.0f, sin(3.14f / 0.25f), cos(3.14f / 0.25f), 0.0f);
      float4 norm = (float4) (0.5f, -0.25f, 0.0f, 0.0f);
      norm = normalize(norm);
      float adj = extstiffness * diff - extdamping * dot(norm, prevvel);
      accel += adj * norm * 0.75f;

      //accel.z =  0.4f;
      //accel.x =  0.4f;
      //accel.y = -0.4f;
    }
  }
#elif 0
  if (flags & DRAIN_MASK)
  {
    float dx = (0 - pos.x) * simscale;            // dist in cm
    float dy = (volumemin.y - pos.y) * simscale;
    float dz = (0 - pos.z) * simscale;
    float dsq = (dx * dx + dy * dy + dz * dz);
    if (0.0001f > dsq)
    {
      pos.x = volumemin.x * 0.5f;
      pos.y = volumemax.y * 0.5f;
      pos.z = volumemin.z;
      accel.z += 50;
      accel.x += 1;
      accel.y += 20;
    }
  }
#else
  if (flags & DRAIN_MASK)
  {
    diff = 2 * radius - (pos.z - volumemin.z - 15) * simscale;
    if (diff < 2 * radius && diff > 0.0001f && (fabs(pos.x) > 3 || fabs(pos.y) > 3))
    {
      norm = (float4) (0, 0, 1, 0);
      //adj = stiff * diff - damp * dot(norm, prevvelocity[i]);
      float adj = extstiffness * diff - extdamping * dot(norm, prevvelocity[i]);
      accel.x += adj * norm.x;
      accel.y += adj * norm.y;
      accel.z += adj * norm.z;
    }
  }
#endif

  /* wave */
  if (flags & WAVE_MASK)
  {
    if (pos.x < 0) accel.x += 20;
  }

  /* fountain */
  if (flags & FOUNTAIN_MASK)
  {
    float dx = (0 - pos.x) * simscale;            // dist in cm
    float dy = (volumemin.y - pos.y) * simscale;
    float dz = (0 - pos.z) * simscale;
    float dsq = (dx * dx + dy * dy + dz * dz);
    if (0.0005f > dsq)
    {
      accel.y += 140;
    }
  }

  accel.y += -9.8f;  //float4(0, -9.8, 0, 0);
  //accel += gravitation;

  // Leapfrog Integration ----------------------------
  float4 vel = velocity[i];
  float4 vnext = accel * deltatime + vel;  // v(t+1/2) = v(t-1/2) + a(t) dt
  prevvel = (vel + vnext) * 0.5f;          // v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5     used to compute forces later
  vel = vnext;
  vnext *= deltatime / simscale;
  pos += vnext;                       // p(t+1) = p(t) + v(t+1/2) dt

  velocity[i] = vel;
  prevvelocity[i] = prevvel;
  position[i] = pos;
}
#endif
