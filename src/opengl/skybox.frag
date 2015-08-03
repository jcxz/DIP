#version 330

uniform samplerCube tex_sky_box_cubemap;

in vec3 tc;

out vec4 frag_color;

void main(void)
{
  frag_color = texture(tex_sky_box_cubemap, tc);
}
