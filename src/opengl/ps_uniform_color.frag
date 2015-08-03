#version 330

uniform vec3 light_pos;    // camera space position of the light
uniform vec3 light_col_a;
uniform vec3 light_col_d;
uniform vec3 light_col_s;

uniform vec3 particle_col;  // particle's color

in vec3 o_normal;        // normal to the point on model's surface
in vec3 o_surf_pos;      // the point for which we are calculating the lighting model


void main(void)
{
  //  gl_FragColor = vec4(1.0f, 1.0f, 1.0f, 1.0f);

  /* normalize the normal */
  vec3 N = normalize(o_normal);

  /* calculate the light direction and normalize it */
  vec3 L = normalize(light_pos - o_surf_pos);

  /* calculate the cosine between N and L vectors */
  float cos_NL = dot(N, L);
  if (cos_NL < 0.0f) cos_NL = 0.0f;

  /* write the resulting color */
  gl_FragColor = vec4(light_col_a * particle_col + light_col_d * particle_col * cos_NL, 1.0f);
}
