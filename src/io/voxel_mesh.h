#ifndef VOXEL_MESH_H
#define VOXEL_MESH_H

#include <vector>
#include <cstdint>
#include <iosfwd>
#include <boost/compute/cl.hpp>


namespace io {

class VoxelMesh
{
  public:
    VoxelMesh(void)
      : m_width(0)
      , m_height(0)
      , m_depth(0)
      , m_data()
    { }

    bool isValid(void) const { return !m_data.empty(); }

    const cl_float4 *data(void) const { return m_data.data(); }
    cl_float4 *data(void) { return m_data.data(); }

    // returns the actual number of voxels
    int count(void) const { return m_data.size(); }

    // these return the dimensions maximal dimensions
    // of the voxelized object
    int width(void) const { return m_width; }
    int height(void) const { return m_height; }
    int depth(void) const { return m_depth; }

    bool loadFromRaw(const char *filename, int w, int h, int d);

    void dump(std::ostream & os);

  private:
    int m_width;
    int m_height;
    int m_depth;
    std::vector<cl_float4> m_data;
};

} // End of namespace io

#endif // VOXEL_MESH_H
