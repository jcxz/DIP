#version 330

const vec2 quad[4] = {
  vec2(-1.0f, -1.0f),
  vec2( 1.0f, -1.0f),
  vec2(-1.0f,  1.0f),
  vec2( 1.0f,  1.0f)
};


out vec2 tc;

void main(void)
{
  vec2 p = quad[gl_VertexID];
  tc = p * 0.5f + 0.5f;
  gl_Position = vec4(p, 0.0f, 1.0f);
}
