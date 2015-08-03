#version 330

uniform mat4 mv;           // model-view matrix for vertices
uniform mat4 proj;         // projection matrix
//uniform mat4 mvp;          // model-view-projection matrix for vertices

uniform bool use_uniform_color;       // whether to use the color from attribute or the color provided by uniform
uniform vec3 particle_uniform_color;  // particle's color

uniform vec2 screen_size;

layout(location = 4) in vec3 particle_pos;  // particle's position
layout(location = 5) in vec3 particle_col;  // particle's color

out vec3 particle_position;
out float particle_radius;
out vec3 particle_color;


void main(void)
{
  // send particle's color to the fragment shader
  particle_color = int(use_uniform_color)  * particle_uniform_color +
                   int(!use_uniform_color) * particle_col;

#if 1
  // transform particle's position and calculate its radius
  vec4 P = mv * vec4(particle_pos, 1.0f);
  particle_position = P.xyz;
  //particle_radius = screen_size.y / (-4.0f * P.z);
  particle_radius = 1.0f / (-4.0f * P.z);

  //particle_radius = 1.0f / 50.0f; //1.0f / (-4.0f * P.z);
  //particle_radius = 1.0f / (-2.0f * P.z);
  //particle_radius = 1.0f / (-0.5f * P.z);

  // set up variables for rasterizer
  gl_Position = proj * P;
  gl_PointSize = screen_size.y * particle_radius;
  //gl_PointSize = particle_radius;
#else
  // transform particle's position and calculate its radius
  vec4 P = mv * vec4(particle_pos, 1.0f);
  particle_position = P.xyz;
  particle_radius = 15.0f;

  // set up variables for rasterizer
  gl_Position = proj * P;
  gl_PointSize = 30.0f;
#endif
}
