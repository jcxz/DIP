#ifndef CURVATURE_FLOW_RENDERER_H
#define CURVATURE_FLOW_RENDERER_H

#include "core/base_renderer.h"


class CurvatureFlowRenderer : public BaseRenderer
{
  public:
    enum DisplayMode {
      MODE_Lighting       = 0,
      MODE_Depth          = 1,
      MODE_SmoothedDepth  = 2,
      MODE_Thickness      = 3,
      MODE_Normals        = 4,
      MODE_PureRefraction = 5
    };

  private:
    static constexpr int CURVATURE_FLOW_STEP_COUNT = int(MODE_Normals) + 1;
    static constexpr int DISPLAY_MODE_COUNT = int(MODE_PureRefraction) + 1;
    static constexpr int DEFAULT_SMOOTH_ITER_CNT = 49; //25; //10; //2; //1; //9; //49;  // use only odd numbers until this gets fixed
    static constexpr int DEFAULT_THICKNESS_SUBSAMPLING = 4; //1; //2;

  public:
    CurvatureFlowRenderer(void)
      : BaseRenderer()
      , m_w(0.0f)
      , m_h(0.0f)
      , m_vao()
      , m_prog_depth()
      , m_prog_thickness()
      , m_prog_smooth()
      , m_prog_compose()
      , m_fbo_depth(0)
      , m_rbo_depth(0)
      , m_tex_depth{ 0, 0 }
      , m_active_depth_tex(0)
      , m_fbo_thickness(0)
      , m_tex_thickness(QOpenGLTexture::Target2D)
      , m_smooth_iter_cnt(DEFAULT_SMOOTH_ITER_CNT)
      , m_thickness_scale(DEFAULT_THICKNESS_SUBSAMPLING)
      , m_display_mode(MODE_Lighting)
    {
      if (!init())
        throw std::runtime_error("Failed to initialize CurvatureFlowRenderer");
    }

    ~CurvatureFlowRenderer(void);

    int smoothingIterationCount(void) const { return m_smooth_iter_cnt; }
    void setSmoothingIterationCount(int cnt) { m_smooth_iter_cnt = cnt; }

    // note: this function has to be called right before reseting
    // this renderer, otherwise it will produce wrong results
    int thicknessSubsamplingFactor(void) const { return m_thickness_scale; }
    void setThicknessSubsamplingFactor(int factor) { m_thickness_scale = factor; }

    DisplayMode displayMode(void) const { return DisplayMode(m_display_mode); }
    void setDisplayMode(DisplayMode mode) { m_display_mode = mode; }
    void toggleDisplayMode(void)
    { m_display_mode = ((m_display_mode + 1) % DISPLAY_MODE_COUNT); }

    void toggleCurvatureFlowSteps(void)
    { m_display_mode = ((m_display_mode + 1) % CURVATURE_FLOW_STEP_COUNT); }

  protected:
    virtual bool resize_impl(int w, int h) override;
    virtual bool reset_impl(int w, int h) override;
    virtual void render_impl(const QQuaternion & rotation,
                             const QVector3D & scale,
                             const QVector3D & translation,
                             const QQuaternion & camera_rotation,
                             GLuint vbo_positions,
                             GLuint /* vbo_colors */,
                             size_t part_cnt,
                             const QVector3D & /* volume_size */) override;

  private:
    void renderThickness(const QMatrix4x4 & mv, size_t part_cnt);
    void renderDepth(const QMatrix4x4 & mv, size_t part_cnt);
    void smoothDepth(void);

    // combine together the thickness, smoothed depth and
    // lighting calculation
    void renderComposed(const QMatrix4x4 & mv, const QMatrix4x4 & camera_rot);

    bool resetDepthFBO(int w, int h);
    bool resetThicknessFBO(int w, int h);

    bool init(void);

    // depth texture management
    GLuint activeDepthTexture(void) const { return m_tex_depth[m_active_depth_tex]; }
    GLuint inactiveDepthTexture(void) const { return m_tex_depth[1 - m_active_depth_tex]; }
    void swapActiveDepthTextures(void) { m_active_depth_tex = 1 - m_active_depth_tex; }
    GLenum drawBufferAttachement(void) { return GL_COLOR_ATTACHMENT0 + 1 - m_active_depth_tex; }

  private:
    // window dimensions
    GLfloat m_w;
    GLfloat m_h;

    // vertex array object (in newer OpenGL
    // it is necessary to have one bound before drawing)
    QOpenGLVertexArrayObject m_vao;

    // OpenGL shaders
    QOpenGLShaderProgram m_prog_depth;
    QOpenGLShaderProgram m_prog_thickness;
    QOpenGLShaderProgram m_prog_smooth;
    QOpenGLShaderProgram m_prog_compose;

    // frame bufer for particle depths and depth smoothing
    GLuint m_fbo_depth;      // frame buffer for ray exit points
    GLuint m_rbo_depth;      // a render buffer object as a storage for depth buffer
    GLuint m_tex_depth[2];   // depth attachement for particle depth frame buffer
    int m_active_depth_tex;  // to keep track of which depth texture buffer is active

    // frame buffer for particle thickness
    GLuint m_fbo_thickness;
    QOpenGLTexture m_tex_thickness;  //

    int m_smooth_iter_cnt;   // number of smoothing iterations
    int m_thickness_scale;   // the factor by how much is the thickness subsampled
    int m_display_mode;      // whether to display normals, depth or thickness or do
                             // regular lighting calculations
};

#endif // CURVATURE_FLOW_RENDERER_H
