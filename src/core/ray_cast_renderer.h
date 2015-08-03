#ifndef RAY_CAST_RENDERER_H
#define RAY_CAST_RENDERER_H

#include "core/base_renderer.h"

#include <QOpenGLVertexArrayObject>
#include <QOpenGLBuffer>
#include <QOpenGLFramebufferObject>
#include <boost/compute/buffer.hpp>
#include <boost/compute/command_queue.hpp>
#include <boost/compute/interop/opengl/opengl_texture.hpp>


class RayCastRenderer : public BaseRenderer
{
  public:
    RayCastRenderer(const boost::compute::context & ctx,
                    const boost::compute::command_queue & queue)
      : BaseRenderer()
      , m_prog_ray_cast()
      , m_vao()
      , m_vbo_cube(QOpenGLBuffer::VertexBuffer)
      , m_ibo_cube(QOpenGLBuffer::IndexBuffer)
      , m_fbo(0)
      , m_color_attach(0)
      , m_data(QOpenGLTexture::Target3D)
      , m_data_width(0)
      , m_data_height(0)
      , m_data_depth(0)
      , m_transfer_func(QOpenGLTexture::Target1D)
      , m_use_transfer_func(true)
      , m_cl_ctx()
      , m_queue()
      , m_prog()
      , m_calculate_particle_cnts_kernel()
      , m_calculate_gradients_kernel()
      , m_calculate_normals_kernel()
      , m_cell_cnts_buf()
      , m_grads_x_buf()
      , m_grads_y_buf()
      , m_grads_z_buf()
      , m_data_cl_img()
      , m_cell_starts_buf()
      , m_cell_ends_buf()
    {
      if (!initGL())
        throw std::runtime_error("RayCastRenderer: Failed to initialize OpenGL.");

      if (!initCL(ctx, queue))
        throw std::runtime_error("RayCastRenderer: Failed to initialize OpenCL.");
    }

    ~RayCastRenderer(void);

    void setGridSize(size_t width, size_t height, size_t depth)
    {
      m_data_width = width;
      m_data_height = height;
      m_data_depth = depth;
    }

    void setGrid(const boost::compute::buffer & cell_starts,
                 const boost::compute::buffer & cell_ends)
    {
      m_cell_starts_buf = cell_starts;
      m_cell_ends_buf = cell_ends;
    }

    void calcNormals(void);

  protected:
    virtual bool resize_impl(int w, int h) override;
    virtual bool reset_impl(int w, int h) override;
    virtual void render_impl(const QQuaternion & rotation,
                             const QVector3D & scale,
                             const QVector3D & translation,
                             const QQuaternion & /* camera_rotation */,
                             GLuint /* vbo_positions */,
                             GLuint /* vbo_colors */,
                             size_t /* part_cnt */,
                             const QVector3D & volume_size) override;
                             //float peel_depth,
                             //int detail

  private:
    bool isDataValid(void) const { return (m_data.isCreated()); }
    bool initGL(void);
    bool initCL(const boost::compute::context & ctx,
                const boost::compute::command_queue & queue);
    bool resetFramebuffer(int w, int h);
    bool resetDataTexture(void);
    bool resetTransferFunctionTexture(void);
    bool resetBuffers(void);

  private:
    // OpenGL rendering programs
    QOpenGLShaderProgram m_prog_gen_back_faces;
    QOpenGLShaderProgram m_prog_ray_cast;

    // OpenGL stuff needed for rendering
    QOpenGLVertexArrayObject m_vao;
    QOpenGLBuffer m_vbo_cube;
    QOpenGLBuffer m_ibo_cube;

    // OpenGL Frame buffer
    GLuint m_fbo;               // frame buffer for ray exit points
    GLuint m_color_attach;      // color attachement for frame buffer

    // The fluid data that has to be rendered
    QOpenGLTexture m_data;      // the actual fluid data
    size_t m_data_width;        // actual width of data texture
    size_t m_data_height;       // actual height of data texture
    size_t m_data_depth;        // actual depth of data texture

    // transfer function texture
    QOpenGLTexture m_transfer_func;
    bool m_use_transfer_func;

    // OpenCL context and command queue
    boost::compute::context m_cl_ctx;
    boost::compute::command_queue m_queue;

    // OpenCL program and kernels
    boost::compute::program m_prog;
    boost::compute::kernel m_calculate_particle_cnts_kernel;
    boost::compute::kernel m_calculate_gradients_kernel;
    boost::compute::kernel m_calculate_normals_kernel;

    // OpenCL buffers
    boost::compute::buffer m_cell_cnts_buf;
    boost::compute::buffer m_grads_x_buf;
    boost::compute::buffer m_grads_y_buf;
    boost::compute::buffer m_grads_z_buf;
    boost::compute::opengl_texture m_data_cl_img;

    // these two buffers are not created in this class, but should be supplied by user
    boost::compute::buffer m_cell_starts_buf;
    boost::compute::buffer m_cell_ends_buf;
};

#endif // RAY_CAST_RENDERER_H
