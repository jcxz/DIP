#ifndef BASE_RENDERER_H
#define BASE_RENDERER_H

#include "utils/ogl.h"
#include "utils/maths.h"

#include <QMatrix4x4>
#include <QOpenGLTexture>
#include <QOpenGLShaderProgram>
#include <QOpenGLBuffer>
#include <QOpenGLVertexArrayObject>


class QQuaternion;
class QVector3D;

class BaseRenderer
{
  private:
    static constexpr float SKY_BOX_SIZE_HALF = 40.0f;
    static constexpr float SKY_BOX_SIZE = 2.0f * SKY_BOX_SIZE_HALF;

  public:
    BaseRenderer(void)
#if 0
      : m_vbo_positions(0)
      , m_vbo_colors(0)
      , m_num_particles(0)
#endif
      : m_proj()
      , m_use_lighting(false)
      , m_light_pos(QVector3D(0.0f, 1.0f, 1.0f))
      , m_light_ambient_col(QVector3D(0.0f, 0.0f, 0.0f))
      , m_light_diffuse_col(QVector3D(0.0f, 0.0f, 0.0f))
      , m_light_specular_col(QVector3D(0.0f, 0.0f, 0.0f))
      , m_use_uniform_color(true)
      , m_particle_col(QVector3D(0.5f, 0.5f, 1.0f))
      , m_default_fbo(0)
      , m_stats()
      , m_prog_sky_box()
      , m_vao_sky_box_cube()
      , m_vbo_sky_box_cube()
      , m_tex_sky_box(QOpenGLTexture::TargetCubeMap)
      , m_draw_sky_box(false)
      , m_prog_bounding_volume()
      , m_vao_bounding_volume()
      //, m_volume_size()
      , m_draw_bounding_volume(true)
      , m_prog_uniform_grid()
      , m_vao_uniform_grid()
      , m_mixing(false)
      , m_clear_col(0.0f, 0.0f, 0.0f, 0.0f)
      , m_print_stats(true)
    {
      if (!init())
        throw std::runtime_error("Failed to initialize BaseRenderer");
    }

    virtual ~BaseRenderer(void)
    {
      if (m_print_stats)
      {
        std::cerr << "Rendering performance statistics:\n" <<  m_stats << std::endl;
      }
    }

#if 0
    // getters/setters for vertex buffer objects shared with OpenCL
    GLuint vboPositionsID(void) const { return m_vbo_positions; }
    void setVboPositionsID(GLuint id) { m_vbo_positions = id; }

    GLuint vboColorsID(void) const { return m_vbo_colors; }
    void setVboColorsID(GLuint id) { m_vbo_colors = id; }

    // particle count settings
    void setParticleCount(size_t part_num) { m_num_particles = part_num; }

    // getters/setters for fluid container size
    void setBoundingVolumeSize(const QVector3D & bv_size) { m_volume_size = bv_size; }
#endif

    // bounding volume display
    void setDrawBoundingVolume(bool enabled) { m_draw_bounding_volume = enabled; }
    bool toggleDrawBoundingVolume(void) { return m_draw_bounding_volume = !m_draw_bounding_volume; }

    // coloring mode
    void setUseUniformColor(bool enabled) { m_use_uniform_color = enabled; }
    void setParticleColor(const QVector3D & particle_col) { m_particle_col = particle_col; }

    // mixing mode
    void setMixing(bool enabled) { m_mixing = enabled; }

    // lighting settings
    void setUseLighting(bool enabled) { m_use_lighting = enabled; }
    void setLightPosition(const QVector3D & light_pos) { m_light_pos = light_pos; }
    void setLightAmbientColor(const QVector3D & ambient_col) { m_light_ambient_col = ambient_col; }
    void setLightDiffuseColor(const QVector3D & diffuse_col) { m_light_diffuse_col = diffuse_col; }
    void setLightSpecularColor(const QVector3D & specular_col) { m_light_specular_col = specular_col; }

    // projection settings
    void setProjection(const QMatrix4x4 & proj) { m_proj = proj; }
    void setPerspectiveProjection(int width, int height)
    {
      m_proj.setToIdentity();
      //m_proj.perspective(30.0f, ((float) width) / ((float) height), 0.01f, 1000.0f);
      //m_proj.perspective(30.0f, ((float) width) / ((float) height), 0.01f, 100.0f);
      //m_proj.perspective(30.0f, ((float) width) / ((float) height), 1.0f, 100.0f);
      //m_proj.perspective(30.0f, ((float) width) / ((float) height), 5.0f, 100.0f);

      //m_proj.perspective(30.0f, ((float) width) / ((float) height), 5.0f,
      //                   std::ceil(utils::maths::cubeDiagonalLength(SKY_BOX_SIZE)) + 5.0f);

      m_proj.perspective(30.0f, ((float) width) / ((float) height), 5.0f, 1000.0f);
      //m_proj.perspective(30.0f, ((float) width) / ((float) height), 5.0f, 200.0f);
      //m_proj.perspective(30.0f, ((float) width) / ((float) height), 10.0f, 100.0f);
    }

    // default frame buffer (necessary due to Qt)
    void setDefaultFBO(int fbo) { m_default_fbo = fbo; }

    // OpenGL framebuffer clear color and Sky Box options
    void setDrawSkyBox(bool enabled) { m_draw_sky_box = enabled; }
    void setClearColor(const QVector4D & c) { m_clear_col = c; }
    bool setSkyBoxTextures(const char *posx, const char *negx,
                           const char *posy, const char *negy,
                           const char *posz, const char *negz);

    // performance counters
    void setPrintStatsOnExit(bool enabled = true) { m_print_stats = enabled; }
    void clearPerformanceCounters(void) { m_stats.clear(); }
    const utils::ogl::PerfStats & performanceCounters(void) const { return m_stats; }

    // initialization and resizing
    bool reset(int w, int h) { return reset_impl(w, h); }
    bool resize(int w, int h);

    // rendering
    void render(const QQuaternion & rotation,
                const QVector3D & scale,
                const QVector3D & translation,
                const QQuaternion & camera_rotation,
                GLuint vbo_positions, GLuint vbo_colors,
                size_t part_cnt, const QVector3D & volume_size);

    //void renderPreview(const QQuaternion & rotation,
    //                   const QVector3D & scale,
    //                   const QVector3D & translation)
    //{ return render(rotation, scale, translation, m_data.maxSize() / 4); }

    void renderUniformGrid(const QQuaternion & rotation,
                           const QVector3D & scale,
                           const QVector3D & translation,
                           const size_t grid_w,
                           const size_t grid_h,
                           const size_t grid_d);

  protected:
    virtual bool resize_impl(int /* w */, int /* h */) { return true; }
    virtual bool reset_impl(int w, int h) = 0;
    virtual void render_impl(const QQuaternion & rotation,
                             const QVector3D & scale,
                             const QVector3D & translation,
                             const QQuaternion & camera_rotation,
                             GLuint vbo_positions,
                             GLuint vbo_colors,
                             size_t part_cnt,
                             const QVector3D & volume_size) = 0;

  protected:
    void renderBBox(const QQuaternion & rotation,
                    const QVector3D & scale,
                    const QVector3D & translation,
                    const QVector3D & size);

    void renderSkyBox(const QQuaternion & rotation,
                      const QVector3D & /* scale */,
                      const QVector3D &translation);

    inline void bindSkyBoxTexture(void) { m_tex_sky_box.bind(); }

  private:
    bool init(void);
    bool initBoundingBox(void);
    bool initUniformGrid(void);
    bool initSkyBox(void);

  protected:
    // simulation data shared with OpenCL
#if 0
    GLuint m_vbo_positions;
    GLuint m_vbo_colors;

    size_t m_num_particles;      // number of particles in simulation
#endif

    // projection
    QMatrix4x4 m_proj;

    // lighting parameters
    bool m_use_lighting;
    QVector3D m_light_pos;
    QVector3D m_light_ambient_col;
    QVector3D m_light_diffuse_col;
    QVector3D m_light_specular_col;

    // color settings
    bool m_use_uniform_color;  // whether to use the same color for all particles or per particle color
    QVector3D m_particle_col;  // color used for drawing the particles

    // id defaultneho frame buffer objektu (aby som vedel ako zresetovat frame buffer
    // po skonceni kreslenia do mojho frame bufer-u, pretoze QOpenGLWidget interne pouziva
    // na kreslenie svoj vlastny frame buffer)
    // default frame buffer id - this is necessary because of Qt,
    // which itself uses frambuffer for rendering
    GLuint m_default_fbo;

    // statistics for meassuring performance
    utils::ogl::PerfStats m_stats;

  private:
    // Stuff for sky box
    QOpenGLShaderProgram m_prog_sky_box;             // shader program to draw a Sky Box
    QOpenGLVertexArrayObject m_vao_sky_box_cube;     // vertex array object for the skybox cube
    QOpenGLBuffer m_vbo_sky_box_cube;                // vertex buffer object for the sky box cube
    QOpenGLTexture m_tex_sky_box;                    // cube map with the sky box texture
    bool m_draw_sky_box;                             // whether to draw sky box or a regular monochrome background

    // Stuff for fluid bounding volume
    QOpenGLShaderProgram m_prog_bounding_volume;     // volume rendering program
    QOpenGLVertexArrayObject m_vao_bounding_volume;  // vertex array object for the bounding volume
    //QVector3D m_volume_size;                       // the size of the fluid container
    bool m_draw_bounding_volume;                     // whether to display bounding volume or not

    // Stuff for uniform grid
    QOpenGLShaderProgram m_prog_uniform_grid;        // uniform grid rendering program
    QOpenGLVertexArrayObject m_vao_uniform_grid;     // vertex array object for uniform grid

    // other display options
    bool m_mixing;                                   // this option is used to sync with simulation,
                                                     // when the mixing option is turned on in the simulation
    QVector4D m_clear_col;                           // the color used for frame buffer clearing

    // statistics
    bool m_print_stats;                              // whether printing statistics in destructor is enabled

};

#endif // BASE_RENDERER_H
