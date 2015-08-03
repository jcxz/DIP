#include "core/fluid_particle_system.h"



cl_float4 FluidParticleSystem::calcGravitationVector(int rx, int ry)
{
  double dir_s, dir_c, mult_x, mult_y, mult_z;
  int corr_z = 1;

  double sinx = sin((double) rx * (3.141592f / 180.0f));
  if (sinx < 0)
  {
    dir_s = sinx;
    dir_c = 1+sinx;
  }
  else
  {
    dir_s = sinx;
    dir_c = 1-sinx;
  }

  if (abs(rx) % 360 > 90 && abs(rx) % 360 < 270)
  //if (fmod(abs(rx), 360) > 90 && fmod(abs(rx), 360) < 270)
  {
      corr_z = -1;
  }

  mult_y = cos((double) ry * (3.141592f / 180.0f));
  mult_x = sin((double) ry * (3.141592f / 180.0f)) * dir_s;
  mult_z = sin((double) ry * (3.141592f / 180.0f)) * dir_c;

  cl_float4 F;

  F.s[0] = mult_x * -9.81f;
  F.s[1] = mult_y * -9.81f;
  F.s[2] = corr_z * mult_z * 9.81f;
  F.s[3] = 0.0f;

  return F; /* vektor gravitace */
}
