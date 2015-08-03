#version 330

// This pass disables depth testing and uses additive blending
// so there is no need to bother with all the depth calculations
// and stuff

in float particle_radius;

out float thickness;


void main(void)
{
  // determine whether the point should be drawn (i.e. is on the sphere)
  vec2 xy = vec2(gl_PointCoord.s * 2.0f - 1.0f, 1.0f - gl_PointCoord.t * 2.0f);
  float xy_sqr = dot(xy, xy);
  if (xy_sqr > 1.0f) discard;

  // write the resulting thickness
  //thickness = 1.0f;
  //thickness = 1.0f - length(xy);
  thickness = 1.0f - sqrt(xy_sqr);
  //thickness = 0.1f;

  //thickness = sqrt(1.0f - xy_sqr) * particle_radius;
}
