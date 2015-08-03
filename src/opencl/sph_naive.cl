//#pragma OPENCL EXTENSION cl_amd_printf : enable


inline float random(ulong seed, float min, float max)
{
  seed = seed * (seed * (seed * (seed * (seed + 1u) + 1u) + 1u) +1u);
  float t = (float) (seed) / (float) ((ulong) (-1));
  return mix(min, max, t);
}


/**
 * A kernel to initialize particle buffers with reasonable default values
 */
__kernel void sph_naive_reset(__global float4 *position,
                              __global float4 *velocity,
                              __global float4 *prev_velocity,
                              __global float *pressure,
                              __global float *density,
                              __global float4 *force,
                              float4 volume_min,
                              float4 volume_max,
                              ulong seed)
{
  ulong gid = get_global_id(0);

#if 0
  position[gid] = (float4) (random(seed + gid + 1ul, volume_min.x, volume_max.x),
                            random(seed + gid + 2ul, volume_min.y, volume_max.y),
                            random(seed + gid + 3ul, volume_min.z, volume_max.z),
                            1.0f);
#else
  int dim_x = (volume_max.x - volume_min.x) * 0.75f;
  int dim_y = (volume_max.y - volume_min.y) * 0.75f;
  //int dim_x = 15;
  //int dim_y = 15;
  //position[gid] = (float4) (gid % dim_x, (gid / dim_x) % dim_y, gid / (dim_x * dim_y), 1.0f);
  position[gid] = (float4) (gid % dim_x           + volume_min.x * 0.75f,
                            (gid / dim_x) % dim_y + volume_min.y * 0.75f,
                            gid / (dim_x * dim_y),
                            1.0f);
#endif

  velocity[gid] = (float4)(0.0f);
  prev_velocity[gid] = (float4)(0.0f);
  pressure[gid] = 0;
  density[gid] = 0;
  force[gid] = (float4)(0.0f);

  //printf("sph_reset: seed == %u\n", seed);
  //printf("sph_reset: position[%u] == [%v4f]\n", gid, position[gid]);
  //printf("sph_reset: position[%u] == [%v4f]\n", gid, velocity[gid]);
  //printf("sph_reset: position[%u] == [%v4f]\n", gid, prev_velocity[gid]);
  //printf("sph_reset: position[%u] == %f\n", gid, pressure[gid]);
  //printf("sph_reset: position[%u] == %f\n", gid, density[gid]);
  //printf("sph_reset: position[%u] == [%v4f]\n", gid, force[gid]);
}


/**
 * The code has been written according to http://www.rchoetzlein.com/eng/graphics/fluids.htm
 * and http://joeyfladderak.com/portfolio-items/sph-fluid-simulation/
 */
__kernel void sph_naive_compute_pressure(__global const float4* pos,
                                         __global       float*  density,
                                         __global       float*  pressure,
                                         float simscale,
                                         //float smoothradius,
                                         float radius2,
                                         //float mass,
                                         //float polykern,
                                         float mass_polykern,
                                         float restdensity,
                                         float intstiffness,
                                         uint numparticles)
{
  int i = get_global_id(0);

  // The calculation is split in two loops to avoid
  // computing with self and to avoid an unnecessary
  // condition insede the for loop

  float sum = 0.0f;

  for (int j = 0; j < i; ++j)
  {
    float4 d = (pos[i] - pos[j]) * simscale;

    // note for myself:
    // the dot product of float4 is defined as x*x + y*y + z*z + w*w,
    // but since we set w to 1.0 in reset kernel and by subtracting
    // pos[i] - pos[j] we get 0,
    // this should not be a problem, in case
    // this is not true, then this has to be changed
    // to dsq = dot(d.xyz, d.xyz) or something similar
    float sqr = dot(d, d);

    if (radius2 > sqr)
    {
      float c = radius2 - sqr;
      sum += c * c * c;
    }
  }

  for (int j = i + 1; j < numparticles; ++j)
  {
    float4 d = (pos[i] - pos[j]) * simscale;
    float sqr = dot(d, d);

    if (radius2 > sqr)
    {
      float c = radius2 - sqr;
      sum += c * c * c;
    }
  }

  // polykern konstanta
  float ro = sum * mass_polykern; //mass * polykern;

  //                          kludova hustota vody pri 20 C     plynova konstanta (ci sa to bude chovat skor ako plyn)
  pressure[i] = (ro - restdensity) * intstiffness;   // vzorec
  density[i] = 1.0f / ro;
}


/**
 * The code has been written according to http://www.rchoetzlein.com/eng/graphics/fluids.htm
 * and http://joeyfladderak.com/portfolio-items/sph-fluid-simulation/
 */
__kernel void sph_naive_compute_force(__global const float4* pos,
                                      __global const float*  density,
                                      __global const float*  pressure,
                                      __global       float4* forces,
                                      __global const float4* vel,
                                      float simscale,
                                      float smoothradius,
                                      float radius2,
                                      //float viscosity,
                                      //float lapkern,
                                      float vterm,
                                      float spikykern_half,
                                      unsigned int numparticles)
{
  int i = get_global_id(0);

  float4 force = (float4) (0.0f, 0.0f, 0.0f, 0.0f);

  //float vterm = lapkern * viscosity;

  for (int j = 0; j < i; ++j)
  {
    float4 d = (pos[i] - pos[j]) * simscale;
    float sqr = dot(d, d);

    if (radius2 > sqr)
    {
      float r = sqrt(sqr);
      float c = (smoothradius - r);
      //float pterm = -0.5f * c * spikeykern * (pressure[i] + pressure[j]) / r;
      float pterm = c * spikykern_half * (pressure[i] + pressure[j]) / r;
      float dterm = c * density[i] * density[j];

      force += (pterm * d + vterm * (vel[j] - vel[i])) * dterm;
    }
  }

  for (int j = i + 1; j < numparticles; ++j)
  {
    float4 d = (pos[i] - pos[j]) * simscale;
    float sqr = dot(d, d);

    if (radius2 > sqr)
    {
      float r = sqrt(sqr);
      float c = (smoothradius - r);
      //float pterm = -0.5f * c * spikeykern * (pressure[i] + pressure[j]) / r;
      float pterm = c * spikykern_half * (pressure[i] + pressure[j]) / r;
      float dterm = c * density[i] * density[j];

      force += (pterm * d + vterm * (vel[j] - vel[i])) * dterm;
    }
  }

  forces[i] = force;
}


/**
 * The code has been written according to http://www.rchoetzlein.com/eng/graphics/fluids.htm
 * and http://joeyfladderak.com/portfolio-items/sph-fluid-simulation/
 */

#define DRAIN_MASK    (1 << 0)
#define WAVE_MASK     (1 << 1)
#define FOUNTAIN_MASK (1 << 2)

__kernel void sph_naive_compute_step(__global       float4* position,
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
  diff = 2.0f * radius - (pos.y - volumemin.y - (pos.x - volumemin.x) * slope) * simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (-slope, 1.0f - slope, 0.0f, 0.0f);
    //norm.x = -slope;
    //norm.z = 0.0f;
    //norm.y = 1.0f - slope;  // = float4(-0, 0, 1.0 - 0, 0);
    float adj = extstiffness * diff - extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * radius - (volumemax.y - pos.y) * simscale;
  if (diff > 0.0001f)
  {
    //float4 norm = (float4) (0.0f, 0.0f, -1.0f, 0.0f);
    float4 norm = (float4) (0.0f, -1.0f, 0.0f, 0.0f);
    //norm.x = 0.0f;
    //norm.z = -1.0f;
    //norm.y = 0.0f;  //float4( 0, 0, -1, 0 );
    float adj = extstiffness * diff - extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* X-axis walls */
  diff = 2.0f * radius - (pos.x - volumemin.x + (sin(time * 10.0f) - 1.0f + (pos.y * 0.025f) * 0.25f) * leftwave) * simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (1.0f, 0.0f, 0.0f, 0.0f);
    //norm.x = 1.0f;
    //norm.y = 0.0f;
    //norm.z = 0.0f;  //float4( 1.0, 0, 0, 0 );
    float adj = (leftwave + 1.0f) * extstiffness * diff - extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * radius - (volumemax.x - pos.x + (sin(time * 10.0f) - 1.0f) * rightwave) * simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (-1.0f, 0.0f, 0.0f, 0.0f);
    //norm.x = -1.0f;
    //norm.y = 0.0f;
    //norm.z = 0.0f;  //float4( -1, 0, 0, 0 );
    float adj = (rightwave + 1.0f) * extstiffness * diff - extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  /* Z-axis walls */
  diff = 2.0f * radius - (pos.z - volumemin.z) * simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (0.0f, 0.0f, 1.0f, 0.0f);
    //norm.x = 0.0f;
    //norm.z = 1.0f;
    //norm.y = 0.0f;  //float4( 0, 1, 0, 0 );
    float adj = extstiffness * diff - extdamping * dot(norm, prevvel);
    accel += adj * norm;
  }

  diff = 2.0f * radius - (volumemax.z - pos.z) * simscale;
  if (diff > 0.0001f)
  {
    float4 norm = (float4) (0.0f, 0.0f, -1.0f, 0.0f);
    //norm.x = 0.0f;
    //norm.z = -1.0f;
    //norm.y = 0.0f;  //float4( 0, -1, 0, 0 );
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
    //if (0.0001f > dsq)
    if (0.005f > dsq)
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
      accel.x += adj * norm.x; accel.y += adj * norm.y; accel.z += adj * norm.z;
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
  float4 vnext = accel * deltatime + vel;  // v(t+1/2) = v(t-1/2) + a(t) dt
  prevvel = (vel + vnext) * 0.5f;          // v(t+1) = [v(t-1/2) + v(t+1/2)] * 0.5     used to compute forces later
  vel = vnext;
  vnext *= deltatime / simscale;
  pos += vnext;                       // p(t+1) = p(t) + v(t+1/2) dt

  velocity[i] = vel;
  prevvelocity[i] = prevvel;
  position[i] = pos;
}
