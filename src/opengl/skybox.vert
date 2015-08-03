#version 330

uniform mat4 mvp;                  // model-view-projection matrix
uniform float sky_box_cube_size;   // skybox scale

in vec3 pos;

out vec3 tc;

void main(void)
{
  //tc = pos * sky_box_cube_size;
  tc = pos;
  gl_Position = mvp * vec4(pos * sky_box_cube_size, 1.0f);
  //gl_Position = mvp * vec4(pos, 1.0f);
}
