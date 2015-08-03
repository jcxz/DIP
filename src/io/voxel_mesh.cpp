#include "io/voxel_mesh.h"
#include "utils/debug.h"

#include <fstream>
#include <QFile>



namespace io {

bool VoxelMesh::loadFromRaw(const char *filename, int w, int h, int d)
{
  QFile file(filename);
  if (!file.open(QFile::ReadOnly))
  {
    ERRORM("Failed to open file: " << filename);
    return false;
  }

  QByteArray data(file.readAll());
  const char * __restrict__ p_data = data.data();

  m_width = w;
  m_height = h;
  m_depth = d;
  m_data.clear();

  for (int k = 0; k < d; ++k)
  {
    for (int j = 0; j < h; ++j)
    {
      for (int i = 0; i < w; ++i)
      {
        int idx = i + j * w + k * w * h;
        if (p_data[idx] != 0)
        {
          //cl_float4 pos = { i, j, k, 1.0f };
          //cl_float4 pos = { k, j, i, 1.0f };
          cl_float4 pos = { j, i, d - k, 1.0f };
          //cl_float4 pos = { j, k, i, 1.0f };
          m_data.push_back(pos);
        }
      }
    }
  }

  return true;
}


void VoxelMesh::dump(std::ostream & os)
{
  int cnt = count();
  if (cnt <= 0)
  {
    WARNM("Not dumping the voxel mesh because it is empty");
    return;
  }

  os << "Voxel count: " << cnt << std::endl;

  for (int i = 0; i < cnt; ++i)
  {
    cl_float4 pos = m_data[i];
    os << "[" << i << "]: ("
       << pos.s[0] << "," << pos.s[1] << "," << pos.s[2] << ")"
       << std::endl;
  }
}

} // End of namspace io
