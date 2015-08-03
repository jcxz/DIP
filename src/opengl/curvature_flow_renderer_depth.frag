#version 330

uniform mat4 proj;         // projection matrix
//uniform vec2 screen_size;

in vec3 particle_position;
in float particle_radius;

out float frag_depth;


void main(void)
{
  // determine whether the point should be drawn (i.e. is on the sphere)
  vec2 xy = vec2(gl_PointCoord.s * 2.0f - 1.0f, 1.0f - gl_PointCoord.t * 2.0f);
  float xy_sqr = dot(xy, xy);
  if (xy_sqr > 1.0f) discard;

  // calculate normal
  // Note that apart from numerical errors this normal should already be normalized,
  // since I am considering the sphere radius as 1
  vec3 N = vec3(xy, sqrt(1.0f - xy_sqr));

  // calculate point on the sphere
  //vec4 P = vec4(particle_position + N * particle_radius / screen_size.y, 1.0f);
  vec4 P = vec4(particle_position + N * particle_radius, 1.0f);

  // calculate depth
  float clip_depth = (proj[0][2] * P.x + proj[1][2] * P.y + proj[2][2] * P.z + proj[3][2]);
  float depth = clip_depth / (proj[0][3] * P.x + proj[1][3] * P.y + proj[2][3] * P.z + proj[3][3]);
  //float depth = (proj[2][2] * P.z + proj[3][2]) / -P.z;   // if I assumed perspective projection
  gl_FragDepth = 0.5f * depth + 0.5f;

  //frag_depth = 0.5f * depth + 0.5f;
  //frag_depth = (normalize(vec4(proj * P))).z * 0.5f + 0.5f;
  //frag_depth = P.z;
  frag_depth = clip_depth;
}
