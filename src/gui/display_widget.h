#ifndef DISPLAY_WIDGET_H
#define DISPLAY_WIDGET_H

#include "utils/gllogger.h"
#include "utils/trackball.h"
#include "io/voxel_mesh.h"
#include "core/base_renderer.h"
#include "core/instancing_renderer.h"
#include "core/point_sprite_renderer.h"
#include "core/curvature_flow_renderer.h"
//#include "core/texture_renderer.h"
#include "core/ray_cast_renderer.h"
#include "core/sph_naive.h"
#include "core/sph_optimized.h"
#include "core/sph_uniform_grid.h"
#include "core/debug_particle_system.h"

#include <QOpenGLWidget>
#include <memory>



class QPainter;

class DisplayWidget : public QOpenGLWidget
{
    Q_OBJECT

  public:
    enum RendererType {
      REN_Debug,
      REN_Instancing,
      REN_PointSprite,
      REN_CurvatureFlow,
      REN_Texture,
      REN_RayCast,
      REN_Invalid
    };

    enum ParticleSystemType {
      PS_SPHUniformGrid,
      PS_SPHOptimized,
      PS_SPHNaive,
      PS_Debug,
      PS_Invalid
    };

    enum ResetMode {
      RM_Prism,
      RM_ArmadilloFilled,
      RM_ArmadilloHollow,
      RM_BunnyFilled,
      RM_BunnyHollow
    };

    enum SkyBoxType {
      SKY_Forest,
      SKY_PalmTrees,
      SKY_Clouds,
      SKY_Debug,
      SKY_DebugLarge
    };

  private:
    static constexpr int RENDERER_COUNT = int(REN_Invalid);
    static constexpr int PARTICLE_SYSTEM_COUNT = int(PS_Invalid);
    static constexpr int DEFAULT_GRID_W = 128;
    static constexpr int DEFAULT_GRID_H = 128;
    static constexpr int DEFAULT_GRID_D = 128;
    static constexpr int DEFAULT_PARTICLE_COUNT = 120000; //1920000; //960000; //480000; //240000; //120000; //60000; //30000; //15000; //10000; //2000;
    static constexpr int SMALL_FONT_SIZE = 8;
    static constexpr int MEDIUM_FONT_SIZE = 10;
    static constexpr int LARGE_FONT_SIZE = 12;
    static constexpr int LINE_SPACING_SIZE = 5;
    static constexpr float DEFAULT_SCALE = 0.05f; //0.15f;

  public:
    explicit DisplayWidget(QWidget *parent = 0)
      : QOpenGLWidget(parent)
      , m_cl_ctx()
      , m_sph_params(nullptr)
      , m_sph_naive(nullptr)
      , m_sph_optimized(nullptr)
      , m_sph_uniform_grid(nullptr)
      , m_debug_system(nullptr)
      , m_cur_ps(nullptr)
      , m_cur_ps_type(PS_Invalid)
      , m_ren_instancing(nullptr)
      , m_ren_point_sprite(nullptr)
      , m_ren_curvature_flow(nullptr)
      , m_ren_ray_cast(nullptr)
      , m_cur_ren(nullptr)
      , m_cur_ren_type(REN_Invalid)
      , m_cur_skybox_type(SKY_Forest)
      , m_initialized(false)
      , m_reset_mode(RM_ArmadilloHollow)
      , m_user_grid_w(DEFAULT_GRID_W)
      , m_user_grid_h(DEFAULT_GRID_H)
      , m_user_grid_d(DEFAULT_GRID_D)
      , m_user_particle_count(DEFAULT_PARTICLE_COUNT)
      , m_track_ball()
      , m_scale(DEFAULT_SCALE)
      , m_track_ball_camera()
      , m_pan_pos()
      , m_pan_dx_prev(0.0f)
      , m_pan_dy_prev(0.0f)
      , m_pan_dx(0.0f)
      , m_pan_dy(0.0f)
      , m_z_dist_pos(0)
      , m_z_dist(0.0f)
      , m_z_dist_prev(0.0f)
      , m_detail(0)
      , m_high_quality(true)
      , m_particle_count(DEFAULT_PARTICLE_COUNT)
      , m_auto_subsampling(true)
      , m_display_info(true)
      , m_display_help(false)
      , m_pause(false)
      , m_step_frame(false)
      , m_calc_normals(true)
      , m_mixing(false)
      , m_sim_frame_num(0)
      , m_frame_num(0)
      , m_gl_timer_id(0)
      , m_total_time_ms(0.0f)
      , m_sim_time_ms(0.0f)
      , m_ren_time_ms(0.0f)
      , m_freeze_timers(false)
      , m_mesh_bunny()
      , m_mesh_armadillo()
      , m_logger()
    {
      // to make keyPressEvent work with this widget
      setFocusPolicy(Qt::ClickFocus);
    }

    ~DisplayWidget(void);

    // particle system
    QString currentParticleSystemName(void) const
    { return particleSystemTypeToString(currentParticleSystemType()); }

    ParticleSystemType currentParticleSystemType(void) const
    { return m_cur_ps_type; }

    static QString particleSystemTypeToString(ParticleSystemType t);

    // renderer
    QString currentRendererName(void) const
    { return rendererTypeToString(m_cur_ren_type); }

    RendererType currentRendererType(void) const
    { return m_cur_ren_type; }

    static QString rendererTypeToString(RendererType t);

  signals:
    void error(const QString & msg);

  public slots:
    void toggleParticleSystem(void);
    void setParticleSystem(ParticleSystemType type);
    void toggleRenderer(void);
    void setRenderer(RendererType type);
    void setDetail(int level);
    void setAutoSubsampling(bool enabled);
    void setHighQuality(void);
    void setLowQuality(void);
    void setDisplayBBox(bool enabled);
    void setUseShading(bool enabled);
    void setLightPosition(const QVector3D & light_pos);
    void setLightAmbientColor(const QVector3D & ambient_col);
    void setLightDiffuseColor(const QVector3D & diffuse_col);

  protected:
    virtual void initializeGL(void) override;
    virtual void paintGL(void) override;
    virtual void resizeGL(int w, int h) override;

    virtual void keyPressEvent(QKeyEvent *event) override;
    virtual void keyReleaseEvent(QKeyEvent *event) override;

    virtual void mousePressEvent(QMouseEvent *event) override;
    virtual void mouseReleaseEvent(QMouseEvent *event) override;
    virtual void mouseMoveEvent(QMouseEvent *event) override;
    virtual void wheelEvent(QWheelEvent *event) override;

    virtual void timerEvent(QTimerEvent * /* event */) override;

  private:
    // helper functions to switch skybox
    void setSkyBox(SkyBoxType type);

    // helper functions to toggle and set particle system
    void setParticleSystem_impl(ParticleSystemType type);
    ParticleSystemType toggleParticleSystemType(void)
    { return ParticleSystemType((int(m_cur_ps_type) + 1) % PARTICLE_SYSTEM_COUNT); }

    // helper functions to toggle and set renderer
    void setRenderer_impl(RendererType type);
    RendererType toggleRendererType(void)
    { return RendererType((int(m_cur_ren_type) + 1) % RENDERER_COUNT); }

    // restart functions
    bool particleSystemReset(void);
    bool rendererReset(void);
    bool simulationReset(void);

    // helper functions to reset manipulators
    void resetXYPanning(void);
    void resetZTranslation(void);

    // helper functions to display information
    void drawHUD(QPainter & painter);
    int displayInfo(QPainter & painter, int height, int font_size);
    int displayHelp(QPainter & painter, int height, int font_size);

    // other helper functions
    io::VoxelMesh *loadVoxelMesh(io::VoxelMesh & mesh);

  private:
    // OpenCL context that is shared by all simulation and rendering classes
    boost::compute::context m_cl_ctx;

    // sph simulation parameters
    std::unique_ptr<SPHParams> m_sph_params;

    // particle system
    std::unique_ptr<SPHNaive> m_sph_naive;
    std::unique_ptr<SPHOptimized> m_sph_optimized;
    std::unique_ptr<SPHUniformGrid> m_sph_uniform_grid;
    std::unique_ptr<DebugParticleSystem> m_debug_system;
    ParticleSystem *m_cur_ps;
    ParticleSystemType m_cur_ps_type;

    // renderer
    std::unique_ptr<InstancingRenderer> m_ren_instancing;
    std::unique_ptr<PointSpriteRenderer> m_ren_point_sprite;
    std::unique_ptr<CurvatureFlowRenderer> m_ren_curvature_flow;
    std::unique_ptr<RayCastRenderer> m_ren_ray_cast;
    BaseRenderer *m_cur_ren;
    RendererType m_cur_ren_type;

    // SkyBox
    SkyBoxType m_cur_skybox_type;

    // a flag to check successfull initialization
    bool m_initialized;

    // how will the simulation be initialized after reset
    ResetMode m_reset_mode;
    size_t m_user_grid_w;
    size_t m_user_grid_h;
    size_t m_user_grid_d;
    unsigned int m_user_particle_count;

    // scene manipulation
    utils::TrackBall m_track_ball;
    float m_scale;

    // camera manipulation
    utils::TrackBall m_track_ball_camera;

    // panning in xy plane
    QPoint m_pan_pos;
    float m_pan_dx_prev;
    float m_pan_dy_prev;
    float m_pan_dx;
    float m_pan_dy;

    // z translation
    int m_z_dist_pos;
    float m_z_dist;
    float m_z_dist_prev;

    // other settings for renderer (in case of renderer switch, they must be set again)
    int m_detail;
    bool m_high_quality;    // whether the user manipulates with the scene (then render in high quality) or not (render in low quality)

    // settings to be used for the simulation and rendering
    size_t m_particle_count;

    // options
    bool m_auto_subsampling;
    bool m_display_info;
    bool m_display_help;
    bool m_pause;
    bool m_step_frame;
    bool m_calc_normals;
    bool m_mixing;          // whether the fluid is going to react to bounding volume rotation

    // frame counters
    unsigned int m_sim_frame_num;  // current frame number of the particle simulation
    unsigned int m_frame_num;      // current frame number (this counter is not stopped when simulation is paused)

    // timers and frame times
    GLuint m_gl_timer_id;          // OpenGL timer query object
    float m_total_time_ms;
    float m_sim_time_ms;
    float m_ren_time_ms;
    bool m_freeze_timers;

    // voxel data used to initialize the fluid
    io::VoxelMesh m_mesh_bunny;
    io::VoxelMesh m_mesh_armadillo;

    // OpenGL debug logging
    utils::ogl::Logger m_logger;
};

#endif // DISPLAYWIDGET_H
