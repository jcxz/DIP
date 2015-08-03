#include "utils/ogl.h"
#include "utils/debug.h"
#include "utils/misc.h"

#include <QOpenGLTexture>
#include <iomanip>



namespace {

// A helper function to load the texture from a given file
// A properly generated texture object has to be bound via
// a call to glBindTexture before calling this function
bool loadTextureHelper(const char *name,
                       GLenum target = GL_TEXTURE_2D,
                       GLint internal_format = GL_RGBA)
{
  QImage img(name);
  if (img.isNull())
  {
    ERRORM("Failed to load the image from file: " << name);
    return false;
  }

  if ((img.format() != QImage::Format_ARGB32) &&
      (img.format() != QImage::Format_RGB32))// &&
      //(img.format() != QImage::Format_RGBA8888))
  {
    img = img.convertToFormat(QImage::Format_ARGB32);
    if (img.isNull())
    {
      ERRORM("Failed to convert image to ARGB32 format");
      return false;
    }
  }

  OGLF->glTexImage2D(target, 0, internal_format, img.width(), img.height(),
                     0, GL_BGRA, GL_UNSIGNED_BYTE, img.constBits());

  return true;
}

} // End of private namespace



namespace utils {

namespace ogl {

///////////////////////////////////////////////////////////////////////////////
// Error handling

const char *errorToStr(GLenum err)
{
  #define CASE(e) case e: return #e

  switch (err)
  {
    CASE(GL_NO_ERROR);
    CASE(GL_INVALID_ENUM);
    CASE(GL_INVALID_VALUE);
    CASE(GL_INVALID_OPERATION);
    CASE(GL_INVALID_FRAMEBUFFER_OPERATION);
    CASE(GL_OUT_OF_MEMORY);
  }

  #undef CASE

  return "UNKNOWN_GL_ERROR";
}

///////////////////////////////////////////////////////////////////////////////
// Miscellaneous utility functions

const char *primitiveToStr(GLenum primitive)
{
  #define CASE(p) case p: return #p

  switch (primitive)
  {
    CASE(GL_POINTS);
    CASE(GL_LINE_STRIP);
    CASE(GL_LINE_LOOP);
    CASE(GL_LINES);
    CASE(GL_LINE_STRIP_ADJACENCY);
    CASE(GL_LINES_ADJACENCY);
    CASE(GL_TRIANGLE_STRIP);
    CASE(GL_TRIANGLE_FAN);
    CASE(GL_TRIANGLES);
    CASE(GL_TRIANGLE_STRIP_ADJACENCY);
    CASE(GL_TRIANGLES_ADJACENCY);
    CASE(GL_PATCHES);
  }

  #undef CASE

  return "Unknown primitive type";
}

///////////////////////////////////////////////////////////////////////////////
// Texture helper functions

bool load2DTexture(QOpenGLTexture & tex,
                   const char *name,
                   GLint internal_format,
                   GLint filter_mode, GLint clamp_mode)
{
  if (!tex.create())
  {
    ERRORM("Failed to create 2D texture's OpenGL object name");
    return false;
  }

  tex.bind();
  utils::misc::ScopedGuard guard([&tex] { tex.release(); });
  (void) guard;

  if (!loadTextureHelper(name, GL_TEXTURE_2D, internal_format)) return false;

  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, filter_mode);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, filter_mode);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, clamp_mode);
  OGLF->glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, clamp_mode);

  return true;
}


bool loadSkyBoxTexture(QOpenGLTexture & tex,
                       const char *posx, const char *negx,
                       const char *posy, const char *negy,
                       const char *posz, const char *negz)
{
  if (!tex.create())
  {
    ERRORM("Failed to create skybox texture's OpenGL object name");
    return false;
  }

  tex.bind();
  utils::misc::ScopedGuard guard([&tex] { tex.release(); });
  (void) guard;

  if (!loadTextureHelper(posx, GL_TEXTURE_CUBE_MAP_POSITIVE_X)) return false;
  if (!loadTextureHelper(negx, GL_TEXTURE_CUBE_MAP_NEGATIVE_X)) return false;
  if (!loadTextureHelper(posy, GL_TEXTURE_CUBE_MAP_POSITIVE_Y)) return false;
  if (!loadTextureHelper(negy, GL_TEXTURE_CUBE_MAP_NEGATIVE_Y)) return false;
  if (!loadTextureHelper(posz, GL_TEXTURE_CUBE_MAP_POSITIVE_Z)) return false;
  if (!loadTextureHelper(negz, GL_TEXTURE_CUBE_MAP_NEGATIVE_Z)) return false;

  OGLF->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  OGLF->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
  OGLF->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);
  OGLF->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
  OGLF->glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);

  return true;
}

///////////////////////////////////////////////////////////////////////////////
// Performance counters and profiling

void PerfStatsRecord::print(const std::string & name, std::ostream & os) const
{
  os << "+----------------------------------------------------------------+" << std::endl;
  os << "| " << std::setw(44) << std::left << name << " | " << std::setw(8) << std::right
             << m_stats.count() << " events |" << std::endl;
  os << "+----------+-----------------+-----------------+-----------------+" << std::endl;
  os << "|   Type   |     min (ms)    |     max (ms)    |     avg (ms)    |" << std::endl;
  os << "+----------+-----------------+-----------------+-----------------+" << std::endl;

  os << std::left;
  os << "| time     | " << std::setw(15) << (m_stats.min())  << " | "
                        << std::setw(15) << (m_stats.max())  << " | "
                        << std::setw(15) << (m_stats.avg())  << " |"
                        << std::endl;
  os << "+----------+-----------------+-----------------+-----------------+" << std::endl;

  os << std::right;

  return;
}


void PerfStats::begin(const std::string & name)
{
  tContainer::iterator it = m_stats.find(name);
  if (it == m_stats.end())
  {
    PerfStatsRecord *rec = new PerfStatsRecord;
    m_stats.insert(std::make_pair(name, std::unique_ptr<PerfStatsRecord>(rec)));
    rec->begin();
  }
  else
  {
    if (!it->second->restart())
    {
      WARNM("The previous OpenGL timer query for the record "
            << name
            << " has not yet finished");
    }
  }
}


void PerfStats::end(const std::string & name)
{
  tContainer::iterator it = m_stats.find(name);
  if (it == m_stats.end())
  {
    WARNM("No record named " << name);
  }
  else
  {
    it->second->end();
  }
}


void PerfStats::printAsTable(std::ostream & os) const
{
  for (auto & it : m_stats)
  {
    it.second->print(it.first, os);
    os << std::endl;
  }
}

} // End of namespace ogl

} // End of namespace utils
