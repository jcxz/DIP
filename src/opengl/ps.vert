#version 330

uniform mat4 mvp;

layout(location = 0) in vec3 pos;           // model's vertex position
layout(location = 1) in vec3 normal;        // model's normal
layout(location = 4) in vec3 particle_pos;  // particle's position
layout(location = 5) in vec4 particle_col;  // particle's color


out vec3 o_normal;
out vec4 o_particle_col;

void main(void)
{
  gl_Position = mvp * vec4(pos + particle_pos, 1.0f);
  o_normal = normal;
  o_particle_col = particle_col;
  //gl_Position = 0.005 * vec4(pos + particle_pos, 1.0f);
  //gl_Position = vec4((pos + particle_pos) * 0.02, 1.0f);
}
