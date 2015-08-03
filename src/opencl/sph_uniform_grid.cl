#ifdef DEBUG_SOURCE
#include "sph_ocl_common.h"
#endif

#pragma OPENCL EXTENSION cl_amd_printf : enable


//=================================================================================================================
// MISC HELPER FUNCTIONS

// A Function to compute a random value in range [min, max] given a seed
inline float random(ulong seed, float min, float max)
{
  seed = seed * (seed * (seed * (seed * (seed + 1u) + 1u) + 1u) +1u);
  float t = (float) (seed) / (float) ((ulong) (-1));
  return mix(min, max, t);
}


//=================================================================================================================
// UNIFORM GRID HELPER FUNCTIONS

// A function to test whether a given cell of uniform grid contains any particles
inline int is_cell_start_valid(uint start)
{
  return start != INVALID_UNIFORM_GRID_CELL_VALUE;
}

// A function to compute 3D cell id given a particle space position
inline int4 cell_id_3D(float4 pos, __constant const tSimParams *params)
{
  //return convert_int4(((pos - params->volumemin) * params->simscale) / params->cell_size);
  //return convert_int4(((pos - params->volumemin) / params->cell_size) * params->simscale);
  //return convert_int4((pos - params->volumemin) * (params->simscale / params->cell_size));
  return convert_int4((pos - params->volumemin) / params->cell_size);
  //return convert_int4(floor((pos - params->volumemin) / params->cell_size));
  //return convert_int4(trunc((pos - params->volumemin) / params->cell_size));
  //return convert_int4((pos - params->volumemin) * params->grid_scale / params->cell_size);
}

// A function to compute 1D cell id given a 3D cell id
inline uint cell_id_1D(int4 cell_id, __constant const tSimParams *params)
{
#if 1
  cell_id.x = cell_id.x & (params->grid_size.x - 1);  // wrap grid, assumes size is power of 2
  cell_id.y = cell_id.y & (params->grid_size.y - 1);
  cell_id.z = cell_id.z & (params->grid_size.z - 1);

  //cell_id.x = cell_id.x % params->grid_size.x;
  //cell_id.y = cell_id.y % params->grid_size.y;
  //cell_id.z = cell_id.z % params->grid_size.z;

  return cell_id.x + cell_id.y * params->grid_size.x +
         cell_id.z * params->grid_size.x * params->grid_size.y;
#else
  cell_id.x = cell_id.x & (params->grid_size - 1);  // wrap grid, assumes size is power of 2
  cell_id.y = cell_id.y & (params->grid_size - 1);
  cell_id.z = cell_id.z & (params->grid_size - 1);
  return cell_id.x + cell_id.y * params->grid_size +
         cell_id.z * params->grid_size * params->grid_size;
#endif
}


//=================================================================================================================
// DEBUGGING FUNCTIONS

/**
 * A kernel to count for each particle the number of neighbouring particles that
 * are within its smoothing radius.
 * This kernel uses the brute force method.
 */
__kernel void sph_uniform_grid_count_interactions(__global   const float4* pos,
                                                  __global         uint *radius_cnts,
                                                  __constant const tSimParams *params)
{
  int gid = get_global_id(0);

  float4 pos_reg = pos[gid];
  uint radius_cnt = 0;

  for (int part = 0; part < params->numparticles; ++part)
  {
    float4 d = (pos_reg - pos[part]) * params->simscale;
    float sqr = d.x * d.x + d.y * d.y + d.z * d.z;
    if (((params->radius2 > sqr) && (gid != part))) ++radius_cnt;
  }

  radius_cnts[gid] = radius_cnt;
}


/**
 * A kernel to count for each particle the number of neighbouring particles
 * that are within its smoothing radius as well as the total number of particles
 * that have to be checked to get the answer.
 * This kernel uses the uniform grid method.
 */
__kernel void sph_uniform_grid_count_grid_interactions(__global   const uint *cell_start,
                                                       __global   const uint *cell_end,
                                                       __global   const float4* pos,
                                                       __global         uint *radius_cnts,
                                                       __global         float *ratios,
                                                       __constant const tSimParams *params)
{
#if 0
  int gid = get_global_id(0);

  float4 pos_reg = pos[gid];
  int4 grid_pos = cell_id_3D(pos_reg, params);
  uint radius_cnt = 0;   // the count of particles that within the radius
  uint total_cnt = 0;    // the count of all particles that had to be checked but were not within the radius

  for (int k = -1; k <= 1; ++k)
  {
    for (int j = -1; j <= 1; ++j)
    {
      for (int i = -1; i <= 1; ++i)
      {
        uint cell_id = cell_id_1D(grid_pos + (int4) (i, j, k, 0), params);
        uint start = cell_start[cell_id];
        if (is_cell_start_valid(start))
        {
          uint end = cell_end[cell_id];
          for (int part = start; part < end; ++part)
          {
            float4 d = (pos_reg - pos[part]) * params->simscale;
            float sqr = d.x * d.x + d.y * d.y + d.z * d.z;
            if (((params->radius2 > sqr) && (gid != part))) ++radius_cnt;
            ++total_cnt;
          }
        }
      }
    }
  }

  radius_cnts[gid] = radius_cnt;
  // compute the percentage of particles that are really within
  // the smoothing radius to the total number of all neighbouring particles
  ratios[gid] = (((float) radius_cnt) / ((float) total_cnt)) * 100.0f;
#else
  int gid = get_global_id(0);

  float4 pos_reg = pos[gid];
  //int4 grid_pos_min = cell_id_3D(pos_reg - params->cell_size * 4.0f, params);
  //int4 grid_pos_max = cell_id_3D(pos_reg + params->cell_size * 4.0f, params);
  int4 grid_pos_min = cell_id_3D(pos_reg - params->cell_size, params);
  int4 grid_pos_max = cell_id_3D(pos_reg + params->cell_size, params);
  uint radius_cnt = 0;   // the count of particles that within the radius
  uint total_cnt = 0;    // the count of all particles that had to be checked but were not within the radius
  int4 grid_pos;

  for (grid_pos.z = grid_pos_min.z; grid_pos.z <= grid_pos_max.z; ++grid_pos.z)
  {
    for (grid_pos.y = grid_pos_min.y; grid_pos.y <= grid_pos_max.y; ++grid_pos.y)
    {
      for (grid_pos.x = grid_pos_min.x; grid_pos.x <= grid_pos_max.x; ++grid_pos.x)
      {
        uint cell_id = cell_id_1D(grid_pos, params);
        uint start = cell_start[cell_id];
        if (is_cell_start_valid(start))
        {
          uint end = cell_end[cell_id];
          for (int part = start; part < end; ++part)
          {
            float4 d = (pos_reg - pos[part]) * params->simscale;
            float sqr = d.x * d.x + d.y * d.y + d.z * d.z;
            if (((params->radius2 > sqr) && (gid != part))) ++radius_cnt;
            ++total_cnt;
          }
        }
      }
    }
  }

  radius_cnts[gid] = radius_cnt;
  // compute the percentage of particles that are really within
  // the smoothing radius to the total number of all neighbouring particles
  ratios[gid] = (((float) radius_cnt) / ((float) total_cnt)) * 100.0f;
#endif
}


__kernel void sph_uniform_grid_calc_interaction_stats(__global const uint *grid_radius_cnts,
                                                      __global const uint *brute_force_radius_cnts,
                                                      __global       uint *diffs)
{
  int gid = get_global_id(0);
  diffs[gid] = abs((int) grid_radius_cnts[gid] - (int) brute_force_radius_cnts[gid]);
}


/**
 * A kernel to help collect statistics on uniform grid occupancy
 */
__kernel void sph_uniform_grid_calc_occupancy(__global   const uint *cell_start,
                                              __global   const uint *cell_end,
                                              __global         uint *part_count)
{
  int gid = get_global_id(0);
  uint start = cell_start[gid];
  uint end = cell_end[gid];
  part_count[gid] = (is_cell_start_valid(start)) ? end - start : 0;
}


//=================================================================================================================
// The ACTUAL SIMULATION

/**
 * A kernel to initialize particle buffers with reasonable default values
 */
__kernel void sph_uniform_grid_reset(__global         float4 *position,
                                     __global         float4 *velocity,
                                     __global         float4 *prev_velocity,
                                     __global         float *pressure,
                                     __global         float *density,
                                     __global         float4 *force,

                                     __global         float4 *position2,
                                     __global         float4 *velocity2,
                                     __global         float4 *prev_velocity2,
                                     __global         float *pressure2,
                                     __global         float *density2,
                                     __global         float4 *force2,

                                     __global         uint *hashes,
                                     __global         uint *indices,

                                     __constant const tSimParams *params)
{
  int gid = get_global_id(0);

#if 0
  float4 pos = (float4) (random(params->seed + gid + 1ul, params->volumemin.x, params->volumemax.x),
                         random(params->seed + gid + 2ul, params->volumemin.y, params->volumemax.y),
                         random(params->seed + gid + 3ul, params->volumemin.z, params->volumemax.z),
                         1.0f);
#else
  int dim_x = (params->volumemax.x - params->volumemin.x) * 0.75f;
  int dim_y = (params->volumemax.y - params->volumemin.y) * 0.75f;
  //int dim_x = 15;
  //int dim_y = 15;
  //float4 pos = (float4) (gid % dim_x, (gid / dim_x) % dim_y, gid / (dim_x * dim_y), 1.0f);
  float4 pos = (float4) (gid % dim_x           + params->volumemin.x * 0.75f,
                         (gid / dim_x) % dim_y + params->volumemin.y * 0.75f,
                         gid / (dim_x * dim_y),
                         1.0f);
#endif

  velocity[gid]       = (float4) (0.0f);
  prev_velocity[gid]  = (float4) (0.0f);
  pressure[gid]       = 0.0f;
  density[gid]        = 0.0f;
  force[gid]          = (float4) (0.0f);
  position[gid]       = pos;

  velocity2[gid]      = (float4) (0.0f);
  prev_velocity2[gid] = (float4) (0.0f);
  pressure2[gid]      = 0.0f;
  density2[gid]       = 0.0f;
  force2[gid]         = (float4) (0.0f);
  position2[gid]      = pos;

  hashes[gid]         = cell_id_1D(cell_id_3D(pos, params), params);
  indices[gid]        = gid;

  //printf("sph_reset: seed == %u\n", params->seed);
  //printf("sph_reset: position[%u] == [%v4f]\n", gid, position[gid]);
  //printf("sph_reset: position[%u] == [%v4f]\n", gid, velocity[gid]);
  //printf("sph_reset: position[%u] == [%v4f]\n", gid, prev_velocity[gid]);
  //printf("sph_reset: position[%u] == %f\n", gid, pressure[gid]);
  //printf("sph_reset: position[%u] == %f\n", gid, density[gid]);
  //printf("sph_reset: position[%u] == [%v4f]\n", gid, force[gid]);
}


__kernel void sph_uniform_grid_reset_grid_cells(__global uint *cell_start)
{
  int i = get_global_id(0);
  cell_start[i] = INVALID_UNIFORM_GRID_CELL_VALUE;
}


__kernel void sph_uniform_grid_calc_hash(__global   const float4 *pos,
                                         __global         uint *hashes,
                                         __global         uint *indices,
                                         __constant const tSimParams *params)
{
#if 0
#if 1
  int i = get_global_id(0);

  // zisti poziciu castice v uniformnej mriezke
  uint4 grid_pos = convert_uint4(pos[i] / params->cell_size);

  // vypocet hashu (v podstate len prevod pozicie v uniformnej mriezke do 1D)
  uint hash = grid_pos.x + grid_pos.y * params->grid_size.x +
              grid_pos.z * params->grid_size.x * params->grid_size.y;

  hashes[i] = hash;
  indices[i] = i;
#else
  int i = get_global_id(0);

  // zisti poziciu castice v uniformnej mriezke
  uint4 grid_pos = convert_uint4(pos[i] / params->cell_size);

  // vypocet hashu (v podstate len prevod pozicie v uniformnej mriezke do 1D)
  uint hash = grid_pos.x + grid_pos.y * params->grid_size +
              grid_pos.z * params->grid_size * params->grid_size;

  hashes[i] = hash;
  indices[i] = i;
#endif
#else
  int i = get_global_id(0);
  hashes[i] = cell_id_1D(cell_id_3D(pos[i], params), params);
  indices[i] = i;
#endif
}


// tento kernel inicializuje uniformnu mriezku a preusporiada castice, tak aby isli podla poradia
__kernel void sph_uniform_grid_reorder(__global   const uint *sorted_hashes,
                                       __global   const uint *sorted_indices,
                                       __global         uint *cell_start,
                                       __global         uint *cell_end,

                                       // stare, nezoradene polozky
                                       __global   const float4* pos,
                                       //__global   const float*  density,
                                       //__global   const float*  pressure,
                                       //__global   const float4* forces,
                                       __global   const float4* vel,
                                       __global   const float4* prevvel,

                                       // nove zoradene polozky
                                       __global         float4* sorted_pos,
                                       //__global         float*  sorted_density,
                                       //__global         float*  sorted_pressure,
                                       //__global         float4* sorted_forces,
                                       __global         float4* sorted_vel,
                                       __global         float4* sorted_prevvel,

                                       // parametre simulacie
                                       __constant const tSimParams *params,

                                       // cache pre hash kluce (musi mat velkost (get_local_size(0) + 1) * sizeof(cl_uint))
                                       __local          uint *hash_cache)
{
  int i = get_global_id(0);
  int li = get_local_id(0);

  // nacitanie hashov z globalnej pamaci do lokalnej pamaci
  if (i < params->numparticles)
  {
    hash_cache[li + 1] = sorted_hashes[i];
    if ((li == 0) && (i > 0)) hash_cache[0] = sorted_hashes[i - 1];
    //if (li == 0) hash_cache[0] = (i > 0) ? sorted_hashes[i - 1] : (uint) (-1);
  }

  barrier(CLK_LOCAL_MEM_FENCE);

  // urcenie cell start a cell end na zaklade hashu, ktory je vlastne len linearizovany index cellu v gride
  if (i < params->numparticles)
  {
    //      prev       !=      current
    //if (hash_cache[li] != hash_cache[li + 1])
    if ((i == 0) || (hash_cache[li] != hash_cache[li + 1]))
    {
      if (i > 0) cell_end[hash_cache[li]] = i;
      cell_start[hash_cache[li + 1]] = i;
    }

    if (i == (params->numparticles - 1)) cell_end[hash_cache[li + 1]] = i + 1;

    // prezoradenie castic podla pola indices
    int old_i = sorted_indices[i];
    sorted_pos[i]      = pos[old_i];
    //sorted_density[i]  = density[old_i];
    //sorted_pressure[i] = pressure[old_i];
    //sorted_forces[i]   = forces[old_i];
    sorted_vel[i]      = vel[old_i];
    sorted_prevvel[i]  = prevvel[old_i];
  }
}


// input: sorted pos, density and pressure arrays and uniform grid
__kernel void sph_uniform_grid_compute_pressure(__global   const uint *cell_start,
                                                __global   const uint *cell_end,
                                                __global   const float4* pos,
                                                __global         float*  density,
                                                __global         float*  pressure,
                                                __constant const tSimParams *params)
{
#if 1
  int gid = get_global_id(0);

  float sum = 0.0f;
  float4 pos_reg = pos[gid];
  int4 grid_pos = cell_id_3D(pos_reg, params);

  // prejdenie okolitymi bunkami mriezky
  for (int k = -1; k <= 1; ++k)
  {
    for (int j = -1; j <= 1; ++j)
    {
      for (int i = -1; i <= 1; ++i)
      {
        uint cell_id = cell_id_1D(grid_pos + (int4) (i, j, k, 0), params);
        uint start = cell_start[cell_id];
        if (is_cell_start_valid(start))
        {
          uint end = cell_end[cell_id];

          // prejdenie bunkou mriezky
          for (int part = start; part < end; ++part)
          {
            float4 d = (pos_reg - pos[part]) * params->simscale;
            float sqr = d.x * d.x + d.y * d.y + d.z * d.z;
            float c = (params->radius2 - sqr) * ((params->radius2 > sqr) & (gid != part));
            sum += c * c * c;
          }
        }
      }
    }
  }

  // polykern konstanta
  float ro = sum * params->mass_polykern; //mass * polykern;

  //                          kludova hustota vody pri 20 C     plynova konstanta (ci sa to bude chovat skor ako plyn)
  pressure[gid] = (ro - params->restdensity) * params->intstiffness;   // vzorec
  density[gid] = 1.0f / ro;
#else
  int gid = get_global_id(0);

  float sum = 0.0f;
  float4 pos_reg = pos[gid];
  //int4 grid_pos_min = cell_id_3D(pos_reg - params->cell_size, params);
  //int4 grid_pos_max = cell_id_3D(pos_reg + params->cell_size, params);
  int4 grid_pos_min = cell_id_3D(pos_reg - params->cell_size * 2.0f, params);
  int4 grid_pos_max = cell_id_3D(pos_reg + params->cell_size * 2.0f, params);
  int4 grid_pos;

  // prejdenie okolitymi bunkami mriezky
  for (grid_pos.z = grid_pos_min.z; grid_pos.z <= grid_pos_max.z; ++grid_pos.z)
  {
    for (grid_pos.y = grid_pos_min.y; grid_pos.y <= grid_pos_max.y; ++grid_pos.y)
    {
      for (grid_pos.x = grid_pos_min.x; grid_pos.x <= grid_pos_max.x; ++grid_pos.x)
      {
        uint cell_id = cell_id_1D(grid_pos, params);
        uint start = cell_start[cell_id];
        if (is_cell_start_valid(start))
        {
          uint end = cell_end[cell_id];

          // prejdenie bunkou mriezky
          for (int part = start; part < end; ++part)
          {
            float4 d = (pos_reg - pos[part]) * params->simscale;
            float sqr = d.x * d.x + d.y * d.y + d.z * d.z;
            float c = (params->radius2 - sqr) * ((params->radius2 > sqr) & (gid != part));
            sum += c * c * c;
          }
        }
      }
    }
  }

  // polykern konstanta
  float ro = sum * params->mass_polykern; //mass * polykern;

  //                          kludova hustota vody pri 20 C     plynova konstanta (ci sa to bude chovat skor ako plyn)
  pressure[gid] = (ro - params->restdensity) * params->intstiffness;   // vzorec
  density[gid] = 1.0f / ro;
#endif
}


__constant int offsets[27][3] = {
   { -1, -1, -1 },
   { -1, -1,  0 },
   { -1, -1,  1 },
   { -1,  0, -1 },
   { -1,  0,  0 },
   { -1,  0,  1 },
   { -1,  1, -1 },
   { -1,  1,  0 },
   { -1,  1,  1 },

   {  0, -1, -1 },
   {  0, -1,  0 },
   {  0, -1,  1 },
   {  0,  0, -1 },
   {  0,  0,  0 },
   {  0,  0,  1 },
   {  0,  1, -1 },
   {  0,  1,  0 },
   {  0,  1,  1 },

   {  1, -1, -1 },
   {  1, -1,  0 },
   {  1, -1,  1 },
   {  1,  0, -1 },
   {  1,  0,  0 },
   {  1,  0,  1 },
   {  1,  1, -1 },
   {  1,  1,  0 },
   {  1,  1,  1 },
};


__kernel //__attribute__((reqd_work_group_size(256, 1, 1)))
         //__attribute__((work_group_size_hint(256, 1, 1)))
         void sph_uniform_grid_compute_force(__global   const uint *cell_start,
                                             __global   const uint *cell_end,
                                             __global   const float4* pos,
                                             __global   const float*  density,
                                             __global   const float*  pressure,
                                             __global         float4* forces,
                                             __global   const float4* vel,
                                             __constant const tSimParams *params)
{
#if 1
  //__local int4 grid_pos[256];
  //__local float4 vel_gid[256];

  int gid = get_global_id(0);
  //int lid = get_local_id(0);

  float4 force = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
  float4 pos_reg = pos[gid];
  int4 grid_pos = cell_id_3D(pos_reg, params);
  //grid_pos[lid] = cell_id_3D(pos_reg, params);
  float pressure_gid = pressure[gid];
  float4 vel_gid = vel[gid];
  //vel_gid[lid] = vel[gid];

  // prejdenie okolitymi bunkami mriezky
  for (int k = -1; k <= 1; ++k)
  {
    for (int j = -1; j <= 1; ++j)
    {
      for (int i = -1; i <= 1; ++i)
      {
        uint cell_id = cell_id_1D(grid_pos + (int4) (i, j, k, 0), params);
        //uint cell_id = cell_id_1D(grid_pos[lid] + (int4) (i, j, k, 0), params);
        uint start = cell_start[cell_id];
        if (is_cell_start_valid(start))
        {
          uint end = cell_end[cell_id];

          // prejdenie bunkou mriezky
          //#pragma unroll
          for (int part = start; part < end; ++part)
          {
            float4 d = (pos_reg - pos[part]) * params->simscale;
            float sqr = d.x * d.x + d.y * d.y + d.z * d.z;

            if ((params->radius2 > sqr) & (gid != part))
            {
              float r = sqrt(sqr);
              float c = (params->smoothradius - r);
              float pterm = c * params->spikykern_half * (pressure_gid + pressure[part]) / r;
              float dterm = c * density[part];

              force += (pterm * d + params->vterm * (vel[part] - vel_gid)) * dterm;
              //force += (pterm * d + params->vterm * (vel[part] - vel_gid[lid])) * dterm;
            }
          }
        }
      }
    }
  }

  forces[gid] = force * density[gid];
#endif

#if 0
  int gid = get_global_id(0);

  float4 force = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
  float4 pos_reg = pos[gid];
  //int4 grid_pos = cell_id_3D(pos_reg, params);
  int4 grid_pos = convert_int4((pos_reg - params->volumemin) / params->cell_size);
  float pressure_gid = pressure[gid];
  float4 vel_gid = vel[gid];

  // prejdenie okolitymi bunkami mriezky
  //for (int k = -1; k <= 1; ++k)
  //{
  //  for (int j = -1; j <= 1; ++j)
  //  {
  //    for (int i = -1; i <= 1; ++i)
  //    {
  //      int cell_id = (grid_pos.x + i) & (params->grid_size.x - 1);  // wrap grid, assumes size is power of 2
  //      int tmp_y = (grid_pos.y + j) & (params->grid_size.y - 1);
  //      int tmp_z = (grid_pos.z + k) & (params->grid_size.z - 1);

  //#pragma unroll
  for (int i = 0; i < 27; ++i)
  {
    int cell_id = (grid_pos.x + offsets[i][0]) & (params->grid_size.x - 1);  // wrap grid, assumes size is power of 2
    int tmp_y   = (grid_pos.y + offsets[i][1]) & (params->grid_size.y - 1);
    int tmp_z   = (grid_pos.z + offsets[i][2]) & (params->grid_size.z - 1);

        cell_id = cell_id + tmp_y * params->grid_size.x + tmp_z * params->grid_size.x * params->grid_size.y;

        uint start = cell_start[cell_id];
        //if (is_cell_start_valid(start))
        if (start != INVALID_UNIFORM_GRID_CELL_VALUE)
        {
          uint end = cell_end[cell_id];

          // prejdenie bunkou mriezky
          //#pragma unroll 10
          for (int part = start; part < end; ++part)
          {
            float4 d = (pos_reg - pos[part]) * params->simscale;
            float sqr = d.x * d.x + d.y * d.y + d.z * d.z;

            if ((params->radius2 > sqr) & (gid != part))
            {
              float r = sqrt(sqr);
              float c = (params->smoothradius - r);
              float pterm = c * params->spikykern_half * (pressure_gid + pressure[part]) / r;
              float dterm = c * density[part];
              force += (pterm * d + params->vterm * (vel[part] - vel_gid)) * dterm;
            }
          }
        }
  }

  //    }
  //  }
  //}

  forces[gid] = force * density[gid];
#endif

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///// EXPERIMENTALNE VERZIE, KDE SOM SA POKUSAL ZNIZIT REGISTER USAGE

#if 0
  int gid = get_global_id(0);

  float4 force = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
  float4 pos_reg = pos[gid];
  //int4 grid_pos = cell_id_3D(pos_reg, params);
  int4 grid_pos = convert_int4((pos_reg - params->volumemin) / params->cell_size);
  float pressure_gid = pressure[gid];
  float4 vel_gid = vel[gid];

  // prejdenie okolitymi bunkami mriezky
  //for (int k = -1; k <= 1; ++k)
  //{
  //  for (int j = -1; j <= 1; ++j)
  //  {
  //    for (int i = -1; i <= 1; ++i)
  //    {
  //      int cell_id = (grid_pos.x + i) & (params->grid_size.x - 1);  // wrap grid, assumes size is power of 2
  //      int tmp_y = (grid_pos.y + j) & (params->grid_size.y - 1);
  //      int tmp_z = (grid_pos.z + k) & (params->grid_size.z - 1);

  //#pragma unroll
  //#pragma unroll 27
  #pragma unroll 1
  for (int i = 0; i < 27; ++i)
  {
    uint cell_id = (grid_pos.x + offsets[i][0]) & (params->grid_size.x - 1);  // wrap grid, assumes size is power of 2
    uint tmp_y   = (grid_pos.y + offsets[i][1]) & (params->grid_size.y - 1);
    uint tmp_z   = (grid_pos.z + offsets[i][2]) & (params->grid_size.z - 1);

        //cell_id = cell_id + tmp_y * params->grid_size.x + tmp_z * params->grid_size.x * params->grid_size.y;
        cell_id = mad24(tmp_z, mul24(params->grid_size.x, params->grid_size.y), mad24(tmp_y, params->grid_size.x, cell_id));

        uint start = cell_start[cell_id];
        //if (is_cell_start_valid(start))
        if (start != INVALID_UNIFORM_GRID_CELL_VALUE)
        {
          uint end = cell_end[cell_id];

          // prejdenie bunkou mriezky
          //#pragma unroll 10
          #pragma unroll 1
          for (int part = start; part < end; ++part)
          {
            float4 d = (pos_reg - pos[part]) * params->simscale;
            //float sqr = d.x * d.x + d.y * d.y + d.z * d.z;
            float r = d.x * d.x + d.y * d.y + d.z * d.z;

            //if ((params->radius2 > sqr) & (gid != part))
            if ((params->radius2 > r) & (gid != part))
            {
              //float r = sqrt(sqr);
              //float c = (params->smoothradius - r);
              //float pterm = c * params->spikykern_half * (pressure_gid + pressure[part]) / r;
              //float dterm = c * density[part];
              //force += (pterm * d + params->vterm * (vel[part] - vel_gid)) * dterm;

#if 0
              r = sqrt(r);
              float pterm = pressure_gid;
              pterm += pressure[part];
              pterm /= r;
              pterm *= params->spikykern_half;
              r = (params->smoothradius - r);
              pterm *= r;
              d *= pterm;
              float dterm = density[part];
              dterm *= r;
              float4 vterm = vel[part];
              vterm -= vel_gid;
              vterm *= params->vterm;
              vterm += d;
              vterm *= dterm;
              force += vterm;
#else
              r = native_sqrt(r);
              float pterm = pressure_gid;
              pterm += pressure[part];
              pterm = native_divide(pterm, r);
              pterm *= params->spikykern_half;
              r = params->smoothradius - r;
              pterm *= r;
              d *= pterm;
              float dterm = density[part];
              dterm *= r;
              float4 vterm = vel[part];
              vterm -= vel_gid;
              //vterm *= params->vterm;
              //vterm += d;
              vterm = mad(vterm, params->vterm, d);
              //vterm *= dterm;
              //force += vterm;
              force = mad(vterm, dterm, force);
#endif
            }
          }
        }
  }

  //    }
  //  }
  //}

  forces[gid] = force * density[gid];
#endif

#if 0
  __local float4 vel_gid[256];
  __local float4 dist[256];
  __local float4 *d;
  __local float4 *p_vel_gid;

  int gid = get_global_id(0);

  {
    int lid = get_local_id(0);
    d         = dist    + lid;
    p_vel_gid = vel_gid + lid;
  }

  float4 force = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
  float4 pos_reg = pos[gid];
  int4 grid_pos = convert_int4((pos_reg - params->volumemin) / params->cell_size);
  float pressure_gid = pressure[gid];
  //float4 vel_gid = vel[gid];
  //vel_gid[lid] = vel[gid];
  *p_vel_gid = vel[gid];

  // prejdenie okolitymi bunkami mriezky
  #pragma unroll 1
  for (int i = 0; i < 27; ++i)
  {
    uint cell_id = (grid_pos.x + offsets[i][0]) & (params->grid_size.x - 1);  // wrap grid, assumes size is power of 2
    uint tmp_y   = (grid_pos.y + offsets[i][1]) & (params->grid_size.y - 1);
    uint tmp_z   = (grid_pos.z + offsets[i][2]) & (params->grid_size.z - 1);
    cell_id = mad24(tmp_z, mul24(params->grid_size.x, params->grid_size.y), mad24(tmp_y, params->grid_size.x, cell_id));

    uint start = cell_start[cell_id];
    if (start != INVALID_UNIFORM_GRID_CELL_VALUE)
    {
      uint end = cell_end[cell_id];

      // prejdenie bunkou mriezky
      #pragma unroll 1
      for (int part = start; part < end; ++part)
      {
        //  float4 d = (pos_reg - pos[part]) * params->simscale;
        //  //float sqr = d.x * d.x + d.y * d.y + d.z * d.z;
        //  float r = d.x * d.x + d.y * d.y + d.z * d.z;

        *d = (pos_reg - pos[part]) * params->simscale;
        float r = (*d).x * (*d).x + (*d).y * (*d).y + (*d).z * (*d).z;

        //if ((params->radius2 > sqr) & (gid != part))
        if ((params->radius2 > r) & (gid != part))
        {
          //float r = sqrt(sqr);
          //float c = (params->smoothradius - r);
          //float pterm = c * params->spikykern_half * (pressure_gid + pressure[part]) / r;
          //float dterm = c * density[part];
          //force += (pterm * d + params->vterm * (vel[part] - vel_gid)) * dterm;

          r = native_sqrt(r);
          float pterm = pressure_gid;
          pterm += pressure[part];
          pterm = native_divide(pterm, r);
          pterm *= params->spikykern_half;
          r = params->smoothradius - r;
          pterm *= r;
          //d *= pterm;
          *d *= pterm;
          float dterm = density[part];
          dterm *= r;
          float4 vterm = vel[part];
          //vterm -= vel_gid;
          //vterm -= vel_gid[lid];
          vterm -= *p_vel_gid;
          //vterm = mad(vterm, params->vterm, d);
          vterm = mad(vterm, params->vterm, *d);
          force = mad(vterm, dterm, force);
        }
      }
    }
  }

  forces[gid] = force * density[gid];
#endif

#if 0
  int gid = get_global_id(0);

  float4 force = (float4) (0.0f, 0.0f, 0.0f, 0.0f);
  float8 reg = (float8) (pos[gid], vel[gid]);
  //float4 pos_reg = pos[gid];
  //int4 grid_pos = convert_int4((pos_reg - params->volumemin) / params->cell_size);
  int4 grid_pos = convert_int4((reg.s0123 - params->volumemin) / params->cell_size);
  float pressure_gid = pressure[gid];
  //float4 vel_gid = vel[gid];

  // prejdenie okolitymi bunkami mriezky
  #pragma unroll 1
  for (int i = 0; i < 27; ++i)
  {
    //uint cell_id = (grid_pos.x + offsets[i][0]) & (params->grid_size.x - 1);  // wrap grid, assumes size is power of 2
    //uint tmp_y   = (grid_pos.y + offsets[i][1]) & (params->grid_size.y - 1);
    //uint tmp_z   = (grid_pos.z + offsets[i][2]) & (params->grid_size.z - 1);
    //cell_id = mad24(tmp_z, mul24(params->grid_size.x, params->grid_size.y), mad24(tmp_y, params->grid_size.x, cell_id));

    //uint cell_id = mad24((grid_pos.z + offsets[i][2]) & (params->grid_size.z - 1),
    //                     mul24(params->grid_size.x, params->grid_size.y),
    //                     mad24((grid_pos.y + offsets[i][1]) & (params->grid_size.y - 1), params->grid_size.x,
    //                           (grid_pos.x + offsets[i][0]) & (params->grid_size.x - 1)));

    uint cell_id = (offsets[i][0]) & (params->grid_size.x - 1);  // wrap grid, assumes size is power of 2
    uint tmp_y   = (offsets[i][1]) & (params->grid_size.y - 1);
    uint tmp_z   = (offsets[i][2]) & (params->grid_size.z - 1);
    cell_id += tmp_y + tmp_z;

    //uint cell_id = (offsets[i][0]) & (params->grid_size.x - 1);  // wrap grid, assumes size is power of 2
    //uint tmp_y   = (offsets[i][1]) & (params->grid_size.y - 1);
    //uint tmp_z   = (offsets[i][2]) & (params->grid_size.z - 1);
    //cell_id = mad24(tmp_z, mul24(params->grid_size.x, params->grid_size.y), mad24(tmp_y, params->grid_size.x, cell_id));

    uint start = cell_start[cell_id];
    if (start != INVALID_UNIFORM_GRID_CELL_VALUE)
    {
      uint end = cell_end[cell_id];

      // prejdenie bunkou mriezky
      #pragma unroll 1
      for (int part = start; part < end; ++part)
      {
        //float4 d = (pos_reg - pos[part]) * params->simscale;
        float4 d = (reg.s0123 - pos[part]) * params->simscale;
        //float sqr = d.x * d.x + d.y * d.y + d.z * d.z;
        float r = d.x * d.x + d.y * d.y + d.z * d.z;

        //if ((params->radius2 > sqr) & (gid != part))
        if ((params->radius2 > r) & (gid != part))
        //if ((params->radius2 > r) && (gid != part))
        //if (fabs(params->radius2 - r) < params->radius2)
        {
          //float r = sqrt(sqr);
          //float c = (params->smoothradius - r);
          //float pterm = c * params->spikykern_half * (pressure_gid + pressure[part]) / r;
          //float dterm = c * density[part];
          //force += (pterm * d + params->vterm * (vel[part] - vel_gid)) * dterm;

          r = native_sqrt(r);
          float pterm = pressure_gid;
          pterm += pressure[part];
          pterm = native_divide(pterm, r);
          pterm *= params->spikykern_half;
          r = params->smoothradius - r;
          pterm *= r;
          d *= pterm;

          float4 vterm = vel[part];
          //vterm -= vel_gid;
          vterm -= reg.s4567;
          //float4 vterm = mad(-1.0f, vel_gid, vel[part]);
          vterm = mad(vterm, params->vterm, d);

          float dterm = density[part];
          dterm *= r;
          force = mad(vterm, dterm, force);
        }
      }
    }
  }

  forces[gid] = force * density[gid];
#endif
}


#define DRAIN_MASK    (1 << 0)
#define WAVE_MASK     (1 << 1)
#define FOUNTAIN_MASK (1 << 2)

#if 0
__kernel void sph_uniform_grid_compute_step(__global         float4* position,
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
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * params->radius - (params->volumemax.y - pos.y) * params->simscale;
  if (diff > 0.0001f)
  {
    //float4 norm = (float4) (0.0f, 0.0f, -1.0f, 0.0f);
    float4 norm = (float4) (0.0f, -1.0f, 0.0f, 0.0f);
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* X-axis walls */
  diff = 2.0f * params->radius - (pos.x - params->volumemin.x + (sin(params->time * 10.0f) - 1.0f + (pos.y * 0.025f) * 0.25f) * params->leftwave) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (1.0f, 0.0f, 0.0f, 0.0f);
    float adj = (params->leftwave + 1.0f) * params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * params->radius - (params->volumemax.x - pos.x + (sin(params->time * 10.0f) - 1.0f) * params->rightwave) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (-1.0f, 0.0f, 0.0f, 0.0f);
    float adj = (params->rightwave + 1.0f) * params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* Z-axis walls */
  diff = 2.0f * params->radius - (pos.z - params->volumemin.z) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (0.0f, 0.0f, 1.0f, 0.0f);
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * params->radius - (params->volumemax.z - pos.z) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (0.0f, 0.0f, -1.0f, 0.0f);
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
    if (0.0005f > dsq)
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


///// TOTO JE ZJEDNODUSENA VERZIA, KTORA PODPORUJE PRELIEVANIE
#if 1
__kernel void sph_uniform_grid_compute_step(__global         float4 *position,
                                            __global   const float4 *forces,
                                            __global         float4 *velocity,
                                            __global         float4 *prevvelocity,
                                            __global         uint *hashes,
                                            __global         uint *indices,
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
#if 0
  float4 vel = velocity[i];
  float4 vnext = accel * params->deltatime + vel;                      // v(t+1/2) = v(t-1/2) + a(t) dt
  velocity[i] = vnext;
  prevvelocity[i] = (vel + vnext) * 0.5f;                              // v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5     used to compute forces later
  //prevvelocity[i] = vel;   // toto funguje tiez, ale ma to horsie matematicke vlastnosti
  position[i] = pos + vnext * (params->deltatime / params->simscale);  // p(t+1) = p(t) + v(t+1/2) dt
#else
  float4 vel = velocity[i];
  float4 vnext = accel * params->deltatime + vel;                      // v(t+1/2) = v(t-1/2) + a(t) dt
  velocity[i] = vnext;
  prevvelocity[i] = (vel + vnext) * 0.5f;                              // v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5     used to compute forces later
  //prevvelocity[i] = vel;   // toto funguje tiez, ale ma to horsie matematicke vlastnosti
  float4 pos_new = pos + vnext * (params->deltatime / params->simscale);  // p(t+1) = p(t) + v(t+1/2) dt
  position[i] = pos_new;
  hashes[i] = cell_id_1D(cell_id_3D(pos_new, params), params);
  indices[i] = i;
#endif
#endif
}
#endif























///// TOTO JE ZJEDNODUSENA VERZIA, KTORA FUNGUJE S tSimParams
#if 0
__kernel void sph_uniform_grid_compute_step(__global         float4* position,
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
  diff = 2.0f * params->radius - (pos.y - params->volumemin.y) * params->simscale;
  if (diff > 0.0001f)
  {
    //float4 norm = (float4) (-params->slope, 1.0f - params->slope, 0.0f, 0.0f);
    float4 norm = (float4) (0.0f, 1.0f, 0.0f, 0.0f);
    //float4 norm = (float4) (-0.2425f, 0.9714f, 0.0f, 0.0f);
    //float4 norm = (float4) (-0.062378f, 0.99805f, 0.0f, 0.0f);
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * params->radius - (params->volumemax.y - pos.y) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (0.0f, -1.0f, 0.0f, 0.0f);
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* X-axis walls */
  diff = 2.0f * params->radius - (pos.x - params->volumemin.x) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (1.0f, 0.0f, 0.0f, 0.0f);
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * params->radius - (params->volumemax.x - pos.x) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (-1.0f, 0.0f, 0.0f, 0.0f);
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* Z-axis walls */
  diff = 2.0f * params->radius - (pos.z - params->volumemin.z) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (0.0f, 0.0f, 1.0f, 0.0f);
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * params->radius - (params->volumemax.z - pos.z) * params->simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (0.0f, 0.0f, -1.0f, 0.0f);
    float adj = params->extstiffness * diff - params->extdamping * dot(norm, prevvel);
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
