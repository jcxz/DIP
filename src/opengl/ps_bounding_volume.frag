#version 330

uniform vec3 col;  // color of the bounding volume


void main(void)
{
  gl_FragColor = vec4(col, 1.0f);
}
