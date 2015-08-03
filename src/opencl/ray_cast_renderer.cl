
#pragma OPENCL EXTENSION cl_khr_3d_image_writes : enable


inline float4 normalize_3d_helper(float4 v)
{
  float len = sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
  if (len > 0.0f)
  {
    v.x /= len;
    v.y /= len;
    v.z /= len;
  }
  return v;
}


inline uint get_cell_cnt(uint i, uint j, uint k,
                         __global const uint *p,
                         uint w, uint h, uint d)
{
  if ((i > 0) && (i < w) && (j > 0) && (j < h) && (k > 0) && (k < d))
  {
    return p[k * w * h + j * w + i];
  }
  else
  {
    return 0;
  }
}


inline float get_grad(uint i, uint j, uint k,
                      __global const float *p,
                      uint w, uint h, uint d)
{
  if ((i > 0) && (i < w) && (j > 0) && (j < h) && (k > 0) && (k < d))
  {
    return p[k * w * h + j * w + i];
  }
  else
  {
    return 0.0f;
  }
}


__kernel void ray_cast_renderer_calculate_particle_cnts(__global const uint *cell_starts,
                                                        __global const uint *cell_ends,
                                                        __global       uint *cell_cnts)
{
  int gid = get_global_id(0);
  //cell_cnts[gid] = (cell_ends[gid] - cell_starts[gid]);// * 100;
  uint start = cell_starts[gid];
  uint end = cell_ends[gid];
  uint cnt = (start == ((uint) -1)) ? 0 : end - start;
  cell_cnts[gid] = cnt; //* 40; //100; //25; //50;
}


// vypocet gradientov pomocou metody centralnych diferencii
__kernel void ray_cast_renderer_calculate_gradients(__global const uint *cell_cnts,
                                                    __global       float *grads_x,
                                                    __global       float *grads_y,
                                                    __global       float *grads_z,
                                                             const uint width,
                                                             const uint height,
                                                             const uint depth)
{
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);

#if 0
  float grad_x = get_cell_cnt(i - 1, j, k, cell_cnts, width, height, depth) -
                 get_cell_cnt(i + 1, j, k, cell_cnts, width, height, depth);
  float grad_y = get_cell_cnt(i, j - 1, k, cell_cnts, width, height, depth) -
                 get_cell_cnt(i, j + 1, k, cell_cnts, width, height, depth);
  float grad_z = get_cell_cnt(i, j, k - 1, cell_cnts, width, height, depth) -
                 get_cell_cnt(i, j, k + 1, cell_cnts, width, height, depth);
#else
  float grad_x = get_cell_cnt(i + 1, j, k, cell_cnts, width, height, depth) -
                 get_cell_cnt(i - 1, j, k, cell_cnts, width, height, depth);
  float grad_y = get_cell_cnt(i, j + 1, k, cell_cnts, width, height, depth) -
                 get_cell_cnt(i, j - 1, k, cell_cnts, width, height, depth);
  float grad_z = get_cell_cnt(i, j, k + 1, cell_cnts, width, height, depth) -
                 get_cell_cnt(i, j, k - 1, cell_cnts, width, height, depth);
#endif

  int idx = k * width * height + j * width + i;
  grads_x[idx] = grad_x;
  grads_y[idx] = grad_y;
  grads_z[idx] = grad_z;
}


__kernel void ray_cast_renderer_calculate_normals(__global     const uint *cell_cnts,
                                                  __global     const float *grads_x,
                                                  __global     const float *grads_y,
                                                  __global     const float *grads_z,
                                                  __write_only       image3d_t tex,
                                                               const uint width,
                                                               const uint height,
                                                               const uint depth)
{
#if 0
  // vyhladenie (spriemerovanie) vypocitanych gradientov
  const int radius = 1;              // velkost okolia pre filtraciu
  const int sample_cnt = 3 * 3 * 3;  // celkovy pocet pixlov vo filtrovanom okoli

  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);
  int idx = k * width * height + j * width + i;

  float4 sum = (float4) (0.0f, 0.0f, 0.0f, 0.0f);

  for (int kk = -radius; kk <= radius; ++kk)
  {
    for (int jj = -radius; jj < radius; ++jj)
    {
      for (int ii = -radius; ii < radius; ++ii)
      {
        sum.x += get_grad(i + ii, j + jj, k + kk, grads_x, width, height, depth);
        sum.y += get_grad(i + ii, j + jj, k + kk, grads_y, width, height, depth);
        sum.z += get_grad(i + ii, j + jj, k + kk, grads_z, width, height, depth);
      }
    }
  }

  sum.x = sum.x / ((float) sample_cnt);
  sum.y = sum.y / ((float) sample_cnt);
  sum.z = sum.z / ((float) sample_cnt);

  float grad_size = sqrt(sum.x * sum.x + sum.y * sum.y + sum.z * sum.z);

  // toto je kvoli tomu, ze v pixely v OpenGL texture maju rozsah <0, 1>
  sum = normalize(sum);
  //sum = normalize_3d_helper(sum);

  // * 0.5f + 0.5f je kvoli tomu, ze v pixely v OpenGL texture maju rozsah <0, 1>
  sum.x = sum.x * 0.5f + 0.5f;
  sum.y = sum.y * 0.5f + 0.5f;
  sum.z = sum.z * 0.5f + 0.5f;
  /*
  //sum.w = cell_cnts[idx] / 255.0f;  // 255.0f je len odhad najvacsej hodnoty
  sum.w = (cell_cnts[idx] * 20.0f) / 255.0f;  // 255.0f je len odhad najvacsej hodnoty
  //sum.w = (cell_cnts[idx] * 25.0f) / 255.0f;  // 255.0f je len odhad najvacsej hodnoty

  //sum.w = cell_cnts[idx];
  */

  sum.w = grad_size / 255.0f;   // 255.0f je len odhad najvacsej hodnoty

  write_imagef(tex, (int4) (i, j, k, 0), sum);
#else
  int i = get_global_id(0);
  int j = get_global_id(1);
  int k = get_global_id(2);

  float4 grad;
  grad.x = get_grad(i, j, k, grads_x, width, height, depth);
  grad.y = get_grad(i, j, k, grads_y, width, height, depth);
  grad.z = get_grad(i, j, k, grads_z, width, height, depth);
  grad.w = 0.0f;

  float grad_size = sqrt(grad.x * grad.x + grad.y * grad.y + grad.z * grad.z);

  grad = normalize(grad);

  // * 0.5f + 0.5f je kvoli tomu, ze v pixely v OpenGL texture maju rozsah <0, 1>
  grad.x = grad.x * 0.5f + 0.5f;
  grad.y = grad.y * 0.5f + 0.5f;
  grad.z = grad.z * 0.5f + 0.5f;
  grad.w = grad_size / 255.0f;   // 255.0f je len odhad najvacsej hodnoty

  //int idx = k * width * height + j * width + i;
  //grad.w = (cell_cnts[idx] * 20.0f) / 255.0f;

  write_imagef(tex, (int4) (i, j, k, 0), grad);
#endif
}
