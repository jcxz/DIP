#version 330

uniform sampler2D tex_particle_depth;
uniform vec2 dx;
uniform vec2 dy;
uniform float cx;
uniform float cy;
uniform float cx2;
uniform float cy2;
uniform float cx2cy2;

in vec2 tc;           // input texture coordinates of the screen space quad

out float frag_depth;



float readDepth(vec2 pos)
{
  return texture2D(tex_particle_depth, pos).r;// * 2.0f - 1.0f;
}

// Mean curvature. From "Screen Space Fluid Rendering with Curvature Flow"
vec3 meanCurvature(vec2 pos, float zc)
{
  // First and second order derivatives of depth from depth buffer
  // Derivatives on edges are forced to be zero to prevent smoothing.
#if 0
  float zdxp = readDepth(pos + dx);
  float zdxn = readDepth(pos - dx);
  float zdx = (zdxp == 0.0f || zdxn == 0.0f) ? 0.0f : 0.5f * (zdxp - zdxn);

  float zdyp = readDepth(pos + dy);
  float zdyn = readDepth(pos - dy);
  float zdy = (zdyp == 0.0f || zdyn == 0.0f) ? 0.0f : 0.5f * (zdyp - zdyn);

  float zdx2 = zdxp + zdxn - 2.0f * zc;
  float zdy2 = zdyp + zdyn - 2.0f * zc;
#else
  vec2 zx = vec2(readDepth(pos + dx), readDepth(pos - dx));
  float zdx = 0.5f * (zx.x - zx.y);
  //zdx = any(equal(zx, vec2(0.0f, 0.0f))) ? 0.0f : zdx;
  zdx = int(all(notEqual(zx, vec2(0.0f, 0.0f)))) * zdx;

  vec2 zy = vec2(readDepth(pos + dy), readDepth(pos - dy));
  float zdy = 0.5f * (zy.x - zy.y);
  //zdy = any(equal(zy, vec2(0.0f, 0.0f))) ? 0.0f : zdy;
  zdy = int(all(notEqual(zy, vec2(0.0f, 0.0f)))) * zdy;

  float zdx2 = zx.x + zx.y - 2.0f * zc;
  float zdy2 = zy.x + zy.y - 2.0f * zc;
#endif

  // Second order finite differences, alternating variables
  float zdxpyp = readDepth(pos + dx + dy);
  float zdxnyn = readDepth(pos - dx - dy);
  float zdxpyn = readDepth(pos + dx - dy);
  float zdxnyp = readDepth(pos - dx + dy);
  float zdxy = (zdxpyp + zdxnyn - zdxpyn - zdxnyp) / 4.0f;

  // Normalization term
  float d = cy2 * zdx * zdx + cx2 * zdy * zdy + cx2cy2 * zc * zc;

  // Derivative of the normalization term d in x and y directotions
  float ddx = cy2 * zdx * zdx2 + cx2 * zdy * zdxy + cx2cy2 * zc * zdx;
  float ddy = cy2 * zdx * zdxy + cx2 * zdy * zdy2 + cx2cy2 * zc * zdy;

  // Mean curvature
  float ex = zdx * ddx - zdx2 * d;
  float ey = zdy * ddy - zdy2 * d;

  float h = 0.5f * ((cy * ex + cx * ey) / pow(d, 1.5f));
  //float h = 0.5f * (cy * ex + cx * ey) * inversesqrt(d) * (1.0f / d);
  //float h = 0.5f * (cy * ex + cx * ey) / (sqrt(d) * d);

  // these all seem to have some numerical issues
  //float h = 0.5f * (cy * ex + cx * ey) / sqrt(d * d * d);
  //float h = 0.5f * (cy * ex + cx * ey) * inversesqrt(d * d * d);
  //float h = 0.5f * (cy * ex + cx * ey) * inversesqrt(d * d * d) * sign(d);
  //float h = 0.5f * (cy * ex + cx * ey) * inversesqrt(abs(d * d * d));
  //float h = 0.5f * (cy * ex + cx * ey) * inversesqrt(abs(d * d * d)) * sign(d);

  return vec3(zdx, zdy, h);
}


void main(void)
{
  float depth = readDepth(tc);
  if (depth == 0.0f)
  {
    frag_depth = 0.0f;
  }
  else
  {
    const float dt = 0.00055f;
    const float dzt = 1000.0f;
    vec3 dxyz = meanCurvature(tc, depth);

    // Vary contribution with absolute depth differential - trick from pySPH
    frag_depth = depth + dxyz.z * dt * (1.0f + (abs(dxyz.x) + abs(dxyz.y)) * dzt);
    //frag_depth = 0.0f;
    //frag_depth = dxyz.z * 0.5f + 0.5f;
    //frag_depth = depth + dxyz.z * dt;
  }
}
