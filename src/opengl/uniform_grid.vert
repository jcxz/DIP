#version 330

uniform mat4 mv;          // model-view matrix for vertices
uniform mat4 proj;        // projection matrix
uniform float cell_size;  // uniform grid cell size
uniform float len;        // uniform grid cell size
uniform int w;            // uniform grid dimensions


void main(void)
{
  int id1 = gl_VertexID / 2;
  int id2 = gl_VertexID % 2;
  vec2 pos = vec2(id1 % w, id1 / w) * cell_size;
  vec4 off = vec4(len, len, len, 0.0f) * 4;
  gl_Position = proj * mv * vec4(vec4(pos, id2 * len * 8, 1.0f) - off);
}
