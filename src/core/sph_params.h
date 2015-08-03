#ifndef SPH_PARAMS_H
#define SPH_PARAMS_H

#include "core/sph_ocl_common.h"

#include <boost/compute/program.hpp>
#include <boost/compute/kernel.hpp>
#include <boost/compute/buffer.hpp>
#include <boost/compute/event.hpp>
#include <ostream>


class QMatrix4x4;
class QStringList;

/**
 * @brief The SPHParams class wraps the SimParams
 * structure that is shared with OpenCL a provides some
 * C++ candy on top of it.
 */
class SPHParams
{
  private:
    static const tFloat4 TOP_FACE_NORMAL;
    static const tFloat4 BOTTOM_FACE_NORMAL;
    static const tFloat4 FRONT_FACE_NORMAL;
    static const tFloat4 BACK_FACE_NORMAL;
    static const tFloat4 LEFT_FACE_NORMAL;
    static const tFloat4 RIGHT_FACE_NORMAL;

    enum Effects {
      EFFECT_DRAIN    = (1 << 0),
      EFFECT_WAVE     = (1 << 1),
      EFFECT_FOUNTAIN = (1 << 2),
      EFFECT_NONE     = (0 << 0)
    };

  public:
    explicit SPHParams(const boost::compute::context & ctx)
      : m_prog()
      , m_print_kernel()
      , m_ev()
      , m_params_buf()
      , m_params()
      , m_wave_start(0.0f)
    { init(ctx); }

    // getters
    size_t gridCellCount(void) const
    {
      return m_params.grid_size.s[0] * m_params.grid_size.s[1] * m_params.grid_size.s[2];
    }

    size_t gridWidth(void) const { return m_params.grid_size.s[0]; }
    size_t gridHeight(void) const { return m_params.grid_size.s[1]; }
    size_t gridDepth(void) const { return m_params.grid_size.s[2]; }

    tFloat smoothRadius(void) const { return m_params.smoothradius; }
    tFloat4 cellSize(void) const { return m_params.cell_size; }
    tUInt4 gridSize(void) const { return m_params.grid_size; }
    tFloat4 volumeMin(void) const { return m_params.volumemin; }
    tFloat4 volumeMax(void) const { return m_params.volumemax; }
    tFloat volumeWidth(void) const { return m_params.volumemax.s[0] - m_params.volumemin.s[0]; }
    tFloat volumeHeight(void) const { return m_params.volumemax.s[1] - m_params.volumemin.s[1]; }
    tFloat volumeDepth(void) const { return m_params.volumemax.s[2] - m_params.volumemin.s[2]; }

    const boost::compute::buffer & buffer(void) const { return m_params_buf; }

    // particle count
    tUInt particleCount(void) const { return m_params.numparticles; }
    void setParticleCount(const tUInt cnt) { m_params.numparticles = cnt; }

    // face rotation
    void rotateTopFace(const QMatrix4x4 & m)
    //{ return rotateFaceNormal(m, m_params.top_face); }
    {
      m_params.top_face = rotateFaceNormal(m, TOP_FACE_NORMAL, m_params.top_face.s[3]);
    }

    void rotateBottomFace(const QMatrix4x4 & m)
    {
      m_params.bottom_face = rotateFaceNormal(m, BOTTOM_FACE_NORMAL, m_params.bottom_face.s[3]);
    }

    void rotateFrontFace(const QMatrix4x4 & m)
    {
      m_params.front_face = rotateFaceNormal(m, FRONT_FACE_NORMAL, m_params.front_face.s[3]);
    }

    void rotateBackFace(const QMatrix4x4 & m)
    {
      m_params.back_face = rotateFaceNormal(m, BACK_FACE_NORMAL, m_params.back_face.s[3]);
    }

    void rotateLeftFace(const QMatrix4x4 & m)
    {
      m_params.left_face = rotateFaceNormal(m, LEFT_FACE_NORMAL, m_params.left_face.s[3]);
    }

    void rotateRightFace(const QMatrix4x4 & m)
    {
      m_params.right_face = rotateFaceNormal(m, RIGHT_FACE_NORMAL, m_params.right_face.s[3]);
    }

    // simulation time
    tFloat time(void) const { return m_params.time; }
    void setTime(const tFloat time) { m_params.time = time; }
    void advanceTime(float time_step) { m_params.time += time_step; }

    // a function to advance the simulation time and
    // all other dependant parameters
    void advanceSimulation(float time_step) { advanceTime(time_step); advanceWave(time_step); }

    // effect manipulation
    tUInt effectFlags(void) const { return m_params.flags; }
    //void setEffectFlags(const tUInt flags) { m_params.flags = flags; }
    void clearEffectFlags(void) { m_params.flags = 0; }
    void activateDrain(void) { deactivateFountain(); m_params.flags |= EFFECT_DRAIN; }
    void deactivateDrain(void) { m_params.flags &= ~(EFFECT_DRAIN); }
    bool toggleDrain(void) { deactivateFountain(); return (m_params.flags ^= EFFECT_DRAIN) != 0; }

    void emitWave(void) { m_params.flags |= EFFECT_WAVE; m_wave_start = m_params.time; }
    void advanceWave(float time_step)
    {
      // kill the wave after 50 frames
      if (m_params.flags & EFFECT_WAVE)
      {
        if ((m_params.time - m_wave_start) > (time_step * 50))
        {
          m_params.flags &= ~(EFFECT_WAVE);
        }
      }
    }

    void activateFountain(void) { deactivateDrain(); m_params.flags |= EFFECT_FOUNTAIN; }
    void deactivateFountain(void) { m_params.flags &= ~(EFFECT_FOUNTAIN); }
    bool toggleFountain(void) { deactivateDrain(); return (m_params.flags ^= EFFECT_FOUNTAIN) != 0; }

    // random seed
    tULong seed(void) const { return m_params.seed; }
    void setSeed(const tULong seed) { m_params.seed = seed; }
    void genSeed(void) { m_params.seed = (tULong) std::time(nullptr); }

    // reinitialization
    void resetToDefaults(void);

    // grid setters
    void setVolume(const tFloat4 & volume_min, const tFloat4 & volume_max);
    void setRoundOffCorrectGridSize(const tUInt4 & size, const tFloat correction = 1.05f);
    void setGridSize(const tUInt4 & size)
    { return setRoundOffCorrectGridSize(size, 1.0f); }

    void setGridSize(const tUInt width, const tUInt height, const tUInt depth)
    { return setRoundOffCorrectGridSize({ width, height, depth, 1 }, 1.0f); }

    void setGridSizeUniform(const tUInt size)
    { return setRoundOffCorrectGridSize({ size, size, size, 1 }, 1.0f); }

    cl_float3 boundingVolumeSize(void) const
    {
      return {
        (m_params.volumemax.s[0] - m_params.volumemin.s[0]) * 0.5f,
        (m_params.volumemax.s[1] - m_params.volumemin.s[1]) * 0.5f,
        (m_params.volumemax.s[2] - m_params.volumemin.s[2]) * 0.5f,
      };
    }

    // GPU handlers
    boost::compute::event & upload(boost::compute::command_queue & queue);
    void enqueuePrintGPUBuffer(boost::compute::command_queue & queue);
    void allocGPUBuffer(const boost::compute::context & ctx)
    { m_params_buf = boost::compute::buffer(ctx, sizeof(m_params)); }

    // debugging
    void dump(std::ostream & os) const;

    // loading from command line arguments
    bool parseCmdArgs(const QStringList & args);

  private:
    void recalcDParams(void);
    tFloat4 rotateFaceNormal(const QMatrix4x4 & mat, const tFloat4 & face, const tFloat d);
    //void rotateFaceNormal(const QMatrix4x4 & mat, tFloat4 & face);
    void init(const boost::compute::context & ctx);

  private:
    boost::compute::program m_prog;
    boost::compute::kernel m_print_kernel;
    boost::compute::event m_ev;             // used to ensure that parameter transfer
                                            // from previous frame is completed
    boost::compute::buffer m_params_buf;
    tSimParams m_params;

    // helper variables
    float m_wave_start;        // a helper variable to control the wave effect
};


inline std::ostream & operator<<(std::ostream & os, const SPHParams & p)
{
  p.dump(os);
  return os;
}

#endif // SPH_PARAMS_H
