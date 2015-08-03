#version 330

uniform vec3 light_pos;
uniform vec4 light_col_a;
uniform vec4 light_col_k;
uniform vec4 light_col_s;

in vec3 o_normal;
in vec4 o_particle_col;

void main(void)
{
  //gl_FragColor = vec4(1.0f, 1.0f, 1.0f, 1.0f);
  float dp = dot(o_normal, light_pos);
  if (dp < 0.0f) dp = 0.0f;
  gl_FragColor = light_col_a * o_particle_col + light_col_k * o_particle_col * dp;
}
