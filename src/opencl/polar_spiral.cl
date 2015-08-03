
inline ulong srand(ulong seed, ulong gid)
{
  return seed * ((gid << 1) + 1) * 69069L + 1;
}

inline float random(ulong *seed, float min, float max)
{
  *seed = *seed * 69069L + 1;               // modulo 2^32 je implicitne
  return ((float) (*seed)) / (((float) ULONG_MAX) + 1); // konverzia do <0, 1)
}



__kernel void polar_spiral(__global float4 *positions,
                           __global float4 *colors,
                           float3 position,
                           ulong seed)
{
  uint gid = get_global_id(0);
  float gidf = gid * 0.1f;

  //seed = srand(seed, gid);
  srand(&seed, gid << 1);

  positions[gid] = (float4) (native_sin(gidf) * (gidf * 1.5f),
                             native_cos(gidf) * gidf,
                             -2.0f * gidf,
                             1.0f);

  colors[gid] = (float4) (random(&seed, 0.1f, 1.0f),
                          random(&seed, 0.1f, 1.0f),
                          random(&seed, 0.1f, 1.0f),
                          1.0f);
}
