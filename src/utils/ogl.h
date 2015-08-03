#ifndef OGL_H
#define OGL_H

#include "utils/stats.h"
#include "utils/debug.h"

#include <QOpenGLContext>
#include <QOpenGLFunctions>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>
#include <QOpenGLTimeMonitor>
#include <memory>
#include <unordered_map>
#include <cassert>


class QOpenGLTexture;

namespace utils {

namespace ogl {

///////////////////////////////////////////////////////////////////////////////
// OpenGL function pointers

inline QOpenGLFunctions_3_3_Core *get_ogl_functions(void)
{
  QOpenGLContext *ctx = QOpenGLContext::currentContext();
  assert(ctx != nullptr);

  QOpenGLFunctions_3_3_Core *f = ctx->versionFunctions<QOpenGLFunctions_3_3_Core>();
  assert(f != nullptr);

  f->initializeOpenGLFunctions();

  return f;
}

inline bool init(void) { return get_ogl_functions() != nullptr; }

///////////////////////////////////////////////////////////////////////////////
// Shader helper functions

inline bool buildShaderProgram(QOpenGLShaderProgram & prog,
                               const QString & vertex_shader,
                               const QString & fragment_shader,
                               const QString & geometry_shader,
                               const QString & tess_ctrl_shader,
                               const QString & tess_eval_shader)
{
  bool ok = true;

  if (!vertex_shader.isNull())    ok = ok && prog.addShaderFromSourceFile(QOpenGLShader::Vertex,   vertex_shader);
  if (!fragment_shader.isNull())  ok = ok && prog.addShaderFromSourceFile(QOpenGLShader::Fragment, fragment_shader);
  if (!geometry_shader.isNull())  ok = ok && prog.addShaderFromSourceFile(QOpenGLShader::Geometry, geometry_shader);
  if (!tess_ctrl_shader.isNull()) ok = ok && prog.addShaderFromSourceFile(QOpenGLShader::TessellationControl, tess_ctrl_shader);
  if (!tess_eval_shader.isNull()) ok = ok && prog.addShaderFromSourceFile(QOpenGLShader::TessellationEvaluation, tess_eval_shader);

  return ok && prog.link();
}


inline bool buildShaderProgram(QOpenGLShaderProgram & prog,
                               const QString & vertex_shader,
                               const QString & fragment_shader)
{
  return buildShaderProgram(prog, vertex_shader, fragment_shader, QString(), QString(), QString());
}

///////////////////////////////////////////////////////////////////////////////
// Error handling

/** A function to convert OpenGL error enum to a string */
const char *errorToStr(GLenum err);

/** A function to directly return an OpenGL error string */
inline const char *errorString(void)
{
  return errorToStr(glGetError());
}

///////////////////////////////////////////////////////////////////////////////
// Miscellaneous utility functions

/** returns a string name for a given OpenGL primtive kind */
const char *primitiveToStr(GLenum primitive);

///////////////////////////////////////////////////////////////////////////////
// Texture helper functions

bool load2DTexture(QOpenGLTexture & tex,
                   const char *name,
                   GLint internal_format = GL_RGBA,
                   GLint filter_mode = GL_LINEAR,
                   GLint clamp_mode = GL_CLAMP_TO_EDGE);

bool loadSkyBoxTexture(QOpenGLTexture & tex,
                       const char *posx, const char *negx,
                       const char *posy, const char *negy,
                       const char *posz, const char *negz);

///////////////////////////////////////////////////////////////////////////////
// Performance counters and profiling

class PerfStatsRecord
{
  private:
    //typedef utils::stats::Statistics<float> tStats;
    typedef utils::stats::Statistics<double> tStats;

  public:
    PerfStatsRecord(void)
      : m_monitor()
      , m_stats()
    {
      m_monitor.setSampleCount(2);
      if (!m_monitor.create()) WARNM("Failed to create QOpenGLTimeMonitor");
    }

    const QOpenGLTimeMonitor & monitorRef(void) const { return m_monitor; }
    const tStats & statsRef(void) const { return m_stats; }

    // a method to begin the timer (corresponds with glBeginQuery)
    void begin(void)
    {
      m_monitor.reset();
      m_monitor.recordSample();
    }

    // a method to end the timer (corresponds with glEndQuery)
    void end(void)
    {
      m_monitor.recordSample();
      //m_monitor.reset();
    }

    // a method to record the meassured milliseconds in statistics
    // this method will block the cpu if necessary
    void flush(void)
    {
      double duration = m_monitor.waitForIntervals()[0];
      m_stats.add(duration / 1000000.0);
    }

    // This method will record the time of previous query and start a new one.
    // In case that the previous query has not finsihed yet, this function is a no op.
    // The method returns true in case a new query has been started and false otherwise.
    bool restart(void)
    {
      if (m_monitor.isResultAvailable())
      {
        flush();
        begin();
        return true;
      }
      return false;
    }

    void print(const std::string & name, std::ostream & os) const;

  private:
    QOpenGLTimeMonitor m_monitor;  // the OpenGL query object
    tStats m_stats;                // statistics from the recorded times
};


class PerfStats
{
  private:
    typedef std::unordered_map<std::string, std::unique_ptr<PerfStatsRecord>> tContainer;

  public:
    PerfStats(void) : m_stats() { }

    ~PerfStats(void) { flush(); }  // wait for all the pending timers to finish

    void clear(void) { return m_stats.clear(); }

    void resetCounter(const std::string & name)
    { m_stats[name].reset(new PerfStatsRecord); }

    // This function will flush all the pending timers
    // and thus might block the cpu
    void flush(void) { for (auto & it : m_stats) it.second->flush(); }

    void begin(const std::string & name);
    void end(const std::string & name);

    void printAsTable(std::ostream & os) const;

    friend std::ostream & operator<<(std::ostream & os, const PerfStats & stats)
    {
      stats.printAsTable(os);
      return os;
    }

  private:
    tContainer m_stats;
};

} // End of namespace ogl

} // End of namespace utils

#define OGLF (utils::ogl::get_ogl_functions())

#endif // OGL_H
