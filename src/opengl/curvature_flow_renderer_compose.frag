#version 330

#define DISPLAY_BY_MODE
//#define DISPLAY_COMPOSED
//#define DISPLAY_COPUTED_NORMALS
//#define DISPLAY_WORLD_SPACE_REFLECTIONS
//#define DISPLAY_VIEW_SPACE_REFLECTIONS

uniform sampler2D tex_particle_depth;
uniform sampler2D tex_particle_thickness;
uniform samplerCube tex_sky_box_cubemap;
uniform vec2 screen_size;
uniform mat4 proj;
uniform mat4 mv;
uniform mat3 camera_mv;
uniform mat3 mv_inv;
uniform float cx;
uniform float cy;
uniform vec2 dx;
uniform vec2 dy;
uniform float proj_00_inv;
uniform float proj_11_inv;
uniform int display_mode;

in vec2 tc;

out vec4 frag_color;



vec3 eyespacePos(vec2 pos, float depth)
{
  pos = (pos - vec2(0.5f)) * 2.0f;
  return depth * vec3(-pos.x * proj[0][0], -pos.y * proj[1][1], 1.0f);
}


// Compute eye-space normal. Adapted from PySPH.
vec3 eyespaceNormal(vec2 pos, float zc)
{
  // give borders a special treatment - in case that
  // one of the derivatives is not defined, try the other one
  float zdxp = texture2D(tex_particle_depth, pos + dx).r;
  float zdxn = texture2D(tex_particle_depth, pos - dx).r;
  float zdx  = (zdxp == 0.0f) ? (zdxn == 0.0f ? 0.0f : (zc - zdxn)) : (zdxp - zc);

  float zdyp = texture2D(tex_particle_depth, pos + dy).r;
  float zdyn = texture2D(tex_particle_depth, pos - dy).r;
  float zdy  = (zdyp == 0.0f) ? (zdyn == 0.0f ? 0.0f : (zc - zdyn)) : (zdyp - zc);

  float wx = proj_00_inv + cx * floor(pos.x * (screen_size.x - 1.0f));
  float wy = proj_11_inv + cy * floor(pos.y * (screen_size.y - 1.0f));

  vec3 pdx = normalize(vec3(cx * zc + wx * zdx, wy * zdx, zdx));
  vec3 pdy = normalize(vec3(wx * zdy, cy * zc + wy * zdy, zdy));

  return normalize(cross(pdx, pdy));
}


// Calculate fresnel coefficient
// Schlicks approximation is for lamers
float fresnel(float rr1, float rr2, vec3 n, vec3 d)
{
  float r = rr1 / rr2;
  float theta1 = dot(n, -d);
  float theta2 = sqrt(1.0f - r * r * (1.0f - theta1 * theta1));

  // Figure out what the Fresnel equations say about what happens next
  float rs = (rr1 * theta1 - rr2 * theta2) / (rr1 * theta1 + rr2 * theta2);
  rs = rs * rs;
  float rp = (rr1 * theta2 - rr2 * theta1) / (rr1 * theta2 + rr2 * theta1);
  rp = rp * rp;

  return (rs + rp) / 2.0f;
}


vec4 compose(float depth)
{
#if 0
  const vec3 light_dir = vec3(1.0f, 1.0f, 1.0f);
  vec3 normal = normalize(eyespaceNormal(tc, depth));
  vec3 eye_pos = normalize(eyespacePos(tc, depth));
  eye_pos.xz = -eye_pos.xz;
#else
  vec3 normal = normalize(eyespaceNormal(tc, depth));
  vec3 eye_pos = eyespacePos(tc, depth);
  eye_pos.xz = -eye_pos.xz;
  //vec3 light_dir = normalize(vec3(10.0f, 10.0f, 10.0f)) - eye_pos;
  vec3 light_dir = vec3(-1000.0f, -1000.0f, -1000.0f) - eye_pos;     // looks fucking awesome
  //vec3 light_dir = vec3(-1000.0f, 1000.0f, -1000.0f) - eye_pos;
  eye_pos = normalize(eye_pos);
#endif

  float lambert = max(0.0f, dot(normalize(light_dir), normal));
  vec3 R = normalize(reflect(eye_pos, normal));
  //vec3 R = normalize(refract(eye_pos, normal, 1.0f / 1.3333f));
  float specular = clamp(fresnel(1.0f, 1.5f, normal, eye_pos), 0.0f, 0.4f);

  //vec4 environment_color = vec4(0.8f, 0.8f, 0.8f, 1.0f);
  //vec4 environment_color = texture(tex_sky_box_cubemap, inverse(mat3(mv)) * R);
  vec4 environment_color = texture(tex_sky_box_cubemap, mv_inv * R);

  float thickness = texture2D(tex_particle_thickness, tc).r / 10.0f;
  //float thickness = texture2D(tex_particle_thickness, tc).r;
  //float thickness = texture2D(tex_particle_thickness, tc).r / 2.0f;

  vec4 particle_color = exp(-vec4(0.6f, 0.2f, 0.05f, 3.0f) * thickness);    // toto je original
  //vec4 particle_color = exp(-vec4(1.0f, 1.0f, 1.0f, 3.0f) * thickness);
  //vec4 particle_color = exp(-vec4(0.05f, 0.05f, 0.05f, 3.0f) * thickness);
  //vec4 particle_color = exp(-vec4(0.5f, 0.5f, 0.05f, 0.5f) * thickness);
  //vec4 particle_color = exp(vec4(0.5f, 0.2f, 0.05f, -3.0f) * thickness);        // fire
  //vec4 particle_color = exp(vec4(-0.6f, 0.2f, 0.05f, -3.0f) * thickness);         // chemical
  //vec4 particle_color = exp(vec4(-0.6f, -0.4f, -0.2f, -1.0f) * thickness);

  particle_color.w    = clamp(1.0f - particle_color.w, 0.0f, 1.0f);
  particle_color.rgb  = (lambert + 0.4f) * particle_color.rgb * (1.0f - specular) +
                        specular * environment_color.rgb;

  return particle_color;
}


// Function to calculate refractions
vec4 pureRefractions(float depth)
{
  vec3 normal  = normalize(eyespaceNormal(tc, depth));
  vec3 I = normalize(eyespacePos(tc, depth));
  //I.xz = -I.xz;
  I = -I;   // takto to nemam flipnute

  //vec3 R = normalize(reflect(I, normal));
  vec3 R = normalize(refract(I, normal, 1.0f / 1.3333f));
  //vec3 R = normalize(refract(I, normal, 1.0f / 1.1f));

  //return texture(tex_sky_box_cubemap, mv_inv * R);
  return texture(tex_sky_box_cubemap, inverse(camera_mv) * R);
  //return texture(tex_sky_box_cubemap, R);
  //return texture(tex_sky_box_cubemap, R * mv_inv);    // toto je zle, tu sa zobrazuju vsetky mozne strany skyboxu
}

//////////////////////////////////////////////////////////////////////////////////
////  The actual composition shader
#ifdef DISPLAY_COMPOSED
void main(void)
{
  float depth = texture2D(tex_particle_depth, tc).r;
  if (depth == 0.0f) discard;
  frag_color = compose(depth);
}
#endif


//////////////////////////////////////////////////////////////////////////////////
////  An alternative composition shader, that can switch between normal rendering
////  vs. displaying depth, thickness or smoothed normals
#ifdef DISPLAY_BY_MODE
void main(void)
{
  float depth = texture2D(tex_particle_depth, tc).r;
  if (depth == 0.0f) discard;

  if (display_mode == 0)
  { // regular display
    frag_color = compose(depth);
  }
  if ((display_mode == 1) || (display_mode == 2))
  { // depth display
    float near = -5.0f;
    float far = -1000.0f;
    vec3 d = vec3(far * (depth + near) / (depth * (far - near))) * 0.5f + 0.5f;
    frag_color = vec4(vec3((d - 0.8f) * 2.0f), 1.0f);
  }
  else if (display_mode == 3)
  { // thickness display
    frag_color = vec4(vec3(texture2D(tex_particle_thickness, tc).r / 50.0f), 1.0f);
    //frag_color = vec4(vec3(texture2D(tex_particle_thickness, tc).r), 1.0f);
  }
  else if (display_mode == 4)
  { // smoothed normals display
    frag_color = vec4(normalize(eyespaceNormal(tc, depth)) * 0.5f + 0.5f, 1.0f);
  }
  else if (display_mode == 5)
  { // pure refractions
    frag_color = pureRefractions(depth);
  }
  else
  { // in any other case default to regular display
    frag_color = compose(depth);
  }
}
#endif


//////////////////////////////////////////////////////////////////////////////////
////  Version to display calculated normal vectors
#ifdef DISPLAY_COPUTED_NORMALS
void main(void)
{
  float depth = texture2D(tex_particle_depth, tc).r;
  if (depth == 0.0f) discard;
  frag_color = vec4(normalize(eyespaceNormal(tc, depth)) * 0.5f + 0.5f, 1.0f);
}
#endif


//////////////////////////////////////////////////////////////////////////////////
////  Experimental pure reflections/refractions version - first transformed
////  to world space and then computed - zda sa, ze toto nefunguje
#ifdef DISPLAY_WORLD_SPACE_REFLECTIONS
void main(void)
{
  float depth = texture2D(tex_particle_depth, tc).r;
  if (depth == 0.0f) discard;

  //vec3 normal = inverse(transpose(mat3(mv))) * eyespaceNormal(tc);
  vec3 normal = mv_inv * eyespaceNormal(tc, depth);
  //vec3 normal = eyespaceNormal(tc) * mv_inv;
  normal = normalize(normal);

  vec3 I = eyespacePos(tc, depth);
  // //I = (vec4(I, 1.0f) * inverse(mv)).xyz;
  // I = (inverse(mv) * vec4(I, 1.0f)).xyz;

  //I = I * mv_inv;
  I = mv_inv * I;

  I = normalize(I);
  I.xz = -I.xz;

  vec3 R = normalize(reflect(I, normal));
  //vec3 R = normalize(refract(I, normal, 1.0f / 1.3333f));
  frag_color = texture(tex_sky_box_cubemap, R);
}
#endif


//////////////////////////////////////////////////////////////////////////////////
////  Experimental pure reflections/refractions version - computed in view space
////  and then transformed to world space
#ifdef DISPLAY_VIEW_SPACE_REFLECTIONS
void main(void)
{
  float depth = texture2D(tex_particle_depth, tc).r;
  if (depth == 0.0f) discard;
  frag_color = pureRefractions(depth);
}
#endif
