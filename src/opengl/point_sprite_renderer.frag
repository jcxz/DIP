#version 330

uniform mat4 proj;         // projection matrix
uniform vec2 screen_size;

uniform vec3 light_pos;    // camera space position of the light
uniform vec3 light_col_a;
uniform vec3 light_col_d;
uniform vec3 light_col_s;

in vec3 particle_color;
in vec3 particle_position;
in float particle_radius;

out vec4 frag_color;


void main(void)
{
  // determine whether the point should be drawn (i.e. is on the sphere)
  vec2 xy = vec2(gl_PointCoord.s * 2.0f - 1.0f, 1.0f - gl_PointCoord.t * 2.0f);
  float xy_sqr = dot(xy, xy);
  if (xy_sqr > 1.0f) discard;

  // calculate normal
  vec3 N = vec3(xy, sqrt(1.0f - xy_sqr));
  //N = normalize(N);

  // calculate point on the sphere
  //vec4 P = vec4(particle_position + N * particle_radius / screen_size.y, 1.0f);
  vec4 P = vec4(particle_position + N * particle_radius, 1.0f);

  // calculate depth
  //vec4 clip_pos = proj * P;
  //float depth = clip_pos.z / clip_pos.w;
  //float depth = dot(proj[2], P) / dot(proj[3], P);
  float depth = (proj[0][2] * P.x + proj[1][2] * P.y + proj[2][2] * P.z + proj[3][2]) /
                (proj[0][3] * P.x + proj[1][3] * P.y + proj[2][3] * P.z + proj[3][3]);
  gl_FragDepth = 0.5f * depth + 0.5f;
  //gl_FragDepth = ((depth + 1.0f) / 2.0) * gl_DepthRange.diff;
  //gl_FragDepth = ((gl_DepthRange.diff * depth) + gl_DepthRange.near + gl_DepthRange.far) / 2.0;

  // calculate the light direction and normalize it
  vec3 L = normalize(light_pos - P.xyz);

  // calculate the cosine between N and L vectors
  float cos_NL = max(dot(N, L), 0.0f);

  // write the resulting color
  frag_color = vec4(light_col_a * particle_color + light_col_d * particle_color * cos_NL, 1.0f);
}
