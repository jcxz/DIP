#ifndef GEOM_H
#define GEOM_H

#include "utils/ogl.h"

#include <ostream>


namespace utils {

namespace geom {

// class encapsulating information on mesh geometry
struct Model
{
  GLenum mode;    /// what type of primitives to render
  GLsizei count;  /// number of vertices
  GLuint vbo;     /// vertex buffer object handle
  GLuint vao;     /// vertex array object handle

  Model(void)
    : mode(GL_TRIANGLES),
      count(0),
      vbo(0),
      vao(0)
  {
    OGLF->glGenVertexArrays(1, &vao);
    OGLF->glGenBuffers(1, &vbo);
  }

  ~Model(void)
  {
    OGLF->glDeleteBuffers(1, &vbo);
    OGLF->glDeleteVertexArrays(1, &vao);
  }

  friend std::ostream & operator<<(std::ostream & os, const Model & geom)
  {
    return os << "Model(" << ogl::primitiveToStr(geom.mode) << ", "
              << geom.count << ", "
              << geom.vbo << ", "
              << geom.vao << ")";
  }
};


/**
 * Generates sphere vertices and loads them to GPU
 *
 * @param mesh mesh object where the information about loaded geometry will be stored
 * @param r radius of the sphere
 *
 * @return true when the geometry has been successfully loaded to GPU, false otherwise
 */
bool genSphere(Model & model, float r = 1.0f);


/**
 * Generates prism vertices and loads them to GPU
 *
 * @param mesh mesh object where the information about loaded geometry will be stored
 * @param a first side of the prism
 * @param a third side of the prism
 * @param a second side of the prism
 *
 * @return true when the geometry has been successfully loaded to GPU, false otherwise
 */
bool genPrism(Model & model, float a = 2.0f, float b = 2.0f, float c = 2.0f);


/**
 */
bool gen2DTriangle(Model & model);

/**
 * Generates square vertices (position and texture coordinates) and loads them to GPU
 */
bool gen2DRectangle(Model & model, float w = 2.0f, float h = 2.0f);

} // End of namesapce geom

} // End of namespace utils

#endif
