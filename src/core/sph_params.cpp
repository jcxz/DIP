#include "core/sph_params.h"
#include "utils/ocl.h"

#include <QFile>
#include <ctime>
#include <cmath>


const tFloat4 SPHParams::TOP_FACE_NORMAL    = {  0.0f, -1.0f,  0.0f, 0.0f };
const tFloat4 SPHParams::BOTTOM_FACE_NORMAL = {  0.0f,  1.0f,  0.0f, 0.0f };
const tFloat4 SPHParams::FRONT_FACE_NORMAL  = {  0.0f,  0.0f, -1.0f, 0.0f };
const tFloat4 SPHParams::BACK_FACE_NORMAL   = {  0.0f,  0.0f,  1.0f, 0.0f };
const tFloat4 SPHParams::LEFT_FACE_NORMAL   = {  1.0f,  0.0f,  0.0f, 0.0f };
const tFloat4 SPHParams::RIGHT_FACE_NORMAL  = { -1.0f,  0.0f,  0.0f, 0.0f };



#if 0
void SPHParams::resetToDefaults(void)
{
#define SIM_SCALE ((tFloat) (0.004f))
#define SMOOTH_RADIUS ((tFloat) (0.01f))
#define RADIUS2 ((tFloat) ((SMOOTH_RADIUS) * (SMOOTH_RADIUS)))

// m = Ro * (V / n)
// Ro   ... hustota tekutiny
// V    ... objem
// n    ... pocet castic
#define MASS ((tFloat) (0.00020543f))
#define POLYKERN ((tFloat) (315.0f / (64.0f * 3.141592 * pow(SMOOTH_RADIUS, 9))))
#define RESTDENSITY ((tFloat) (600.0f))
#define INTSTIFFNESS ((tFloat) (1.0f))
#define MASS_POLYKERN ((tFloat) ((MASS) * (POLYKERN)))

#define VISCOSITY ((tFloat) (0.2f))
#define LAPKERN ((tFloat) (45.0f / (3.141592 * pow(SMOOTH_RADIUS, 6))))
#define VTERM ((tFloat) ((LAPKERN) * (VISCOSITY)))
#define SPIKEYKERN ((tFloat) (-45.0f / (3.141592 * pow(SMOOTH_RADIUS, 6))))
#define SPIKEYKERN_HALF ((tFloat) ((SPIKEYKERN) * (-0.5f)))

#define SLOPE ((tFloat) (0.0f))        // ???
#define LEFTWAVE ((tFloat) (0.0f))     // ???
#define RIGHTWAVE ((tFloat) (0.0f))   // ???
#define DELTATIME ((tFloat) (.003f))   // ???
#define LIMIT ((tFloat) (200.0f))
#define EXTSTIFFNESS ((tFloat) (10000.0f))
#define EXTDAMPING ((tFloat) (256.0f))

//#define RADIUS ((tFloat) (0.004f))
//#define RADIUS ((tFloat) (1.0f / 64.0f))
#define RADIUS ((tFloat) (1.0f / 96.0f))
//#define RADIUS ((tFloat) (1.0f / 128.0f))

  m_params.volumemin = { 0.0f, 0.0f, 0.0f, 0.0f };
  m_params.volumemax = { 0.0f, 0.0f, 0.0f, 0.0f };
  m_params.seed = (tULong) (std::time(nullptr));
  m_params.simscale = SIM_SCALE;
  m_params.radius2 = RADIUS2;
  //cl_float polykern = POLYKERN;
  m_params.mass_polykern = MASS_POLYKERN;
  m_params.restdensity = RESTDENSITY;
  m_params.intstiffness = INTSTIFFNESS;
  m_params.numparticles = 0;
  m_params.smoothradius = SMOOTH_RADIUS;
  //m_params.viscosity = VISCOSITY;
  //m_params.lapkern = LAPKERN;
  m_params.vterm = VTERM;
  m_params.spikykern_half = SPIKEYKERN_HALF;
  m_params.slope = SLOPE;
  m_params.leftwave = LEFTWAVE;
  m_params.rightwave = RIGHTWAVE;
  m_params.deltatime = DELTATIME;
  m_params.limit = LIMIT;
  m_params.extstiffness = EXTSTIFFNESS;
  m_params.extdamping = EXTDAMPING;
  m_params.radius = RADIUS;
  m_params.mass = MASS;
  m_params.time = 0.0f;
  m_params.flags = 0;

  //m_params.cell_size = { 2.0f * radius, 2.0f * radius, 2.0f * radius, 0.0f };
  m_params.cell_size = { 0.0f, 0.0f, 0.0f, 1.0f };
  m_params.grid_size = { 64, 64, 64, 1 };
  //m_params.cell_size = 0.0f;
  //m_params.grid_size = 64;
}
#else
void SPHParams::resetToDefaults(void)
{
#define SIM_SCALE ((tFloat) (0.005f))
#define SMOOTH_RADIUS ((tFloat) (0.01f))
#define RADIUS2 ((tFloat) ((SMOOTH_RADIUS) * (SMOOTH_RADIUS)))

// m = Ro * (V / n)
// Ro   ... hustota tekutiny
// V    ... objem
// n    ... pocet castic
#define MASS ((tFloat) (0.00020543f))
#define POLYKERN ((tFloat) (315.0f / (64.0f * 3.141592 * pow(SMOOTH_RADIUS, 9))))
#define RESTDENSITY ((tFloat) (600.0f))
#define INTSTIFFNESS ((tFloat) (1.5f))
#define MASS_POLYKERN ((tFloat) ((MASS) * (POLYKERN)))

#define VISCOSITY ((tFloat) (0.35f))
#define LAPKERN ((tFloat) (45.0f / (3.141592 * pow(SMOOTH_RADIUS, 6))))
#define VTERM ((tFloat) ((LAPKERN) * (VISCOSITY)))
#define SPIKEYKERN ((tFloat) (-45.0f / (3.141592 * pow(SMOOTH_RADIUS, 6))))
#define SPIKEYKERN_HALF ((tFloat) ((SPIKEYKERN) * (-0.5f)))

#define SLOPE ((tFloat) (0.0f))        // ???
#define LEFTWAVE ((tFloat) (0.0f))     // ???
#define RIGHTWAVE ((tFloat) (0.0f))   // ???
#define DELTATIME ((tFloat) (0.003f))   // ???
#define LIMIT ((tFloat) (150.0f))
#define EXTSTIFFNESS ((tFloat) (50000.0f))
#define EXTDAMPING ((tFloat) (100.0f))

#define RADIUS ((tFloat) (1.0f / 96.0f))
//#define RADIUS ((tFloat) (0.02f))

  m_params.volumemin = { 0.0f, 0.0f, 0.0f, 0.0f };
  m_params.volumemax = { 0.0f, 0.0f, 0.0f, 0.0f };
  m_params.seed = (tULong) (std::time(nullptr));
  m_params.simscale = SIM_SCALE;
  m_params.radius2 = RADIUS2;
  //cl_float polykern = POLYKERN;
  m_params.mass_polykern = MASS_POLYKERN;
  m_params.restdensity = RESTDENSITY;
  m_params.intstiffness = INTSTIFFNESS;
  m_params.numparticles = 0;
  m_params.smoothradius = SMOOTH_RADIUS;
  //m_params.viscosity = VISCOSITY;
  //m_params.lapkern = LAPKERN;
  m_params.vterm = VTERM;
  m_params.spikykern_half = SPIKEYKERN_HALF;
  m_params.slope = SLOPE;
  m_params.leftwave = LEFTWAVE;
  m_params.rightwave = RIGHTWAVE;
  m_params.deltatime = DELTATIME;
  m_params.limit = LIMIT;
  m_params.extstiffness = EXTSTIFFNESS;
  m_params.extdamping = EXTDAMPING;
  m_params.radius = RADIUS;
  m_params.mass = MASS;
  m_params.time = DELTATIME;
  m_params.flags = 0;

  m_params.cell_size = { 0.0f, 0.0f, 0.0f, 1.0f };
  m_params.grid_size = { 1, 1, 1, 1 };
  m_params.grid_scale = { 1.0f, 1.0f, 1.0f, 1.0f };

  // nastavenie stien bounding objemu kvapaliny
  m_params.top_face    = TOP_FACE_NORMAL;    //{  0.0f, -1.0f,  0.0f, 0.0f };
  m_params.bottom_face = BOTTOM_FACE_NORMAL; //{  0.0f,  1.0f,  0.0f, 0.0f };
  m_params.front_face  = FRONT_FACE_NORMAL;  //{  0.0f,  0.0f, -1.0f, 0.0f };
  m_params.back_face   = BACK_FACE_NORMAL;   //{  0.0f,  0.0f,  1.0f, 0.0f };
  m_params.left_face   = LEFT_FACE_NORMAL;   //{  1.0f,  0.0f,  0.0f, 0.0f };
  m_params.right_face  = RIGHT_FACE_NORMAL;  //{ -1.0f,  0.0f,  0.0f, 0.0f };
}
#endif


void SPHParams::setVolume(const tFloat4 & volume_min, const tFloat4 & volume_max)
{
  // THIS IS CORRECT
  //volume_min = { -32.0f, -32.0f, -32.0f, 1.0f };
  //volume_max = {  32.0f,  32.0f,  32.0f, 1.0f };

  //volume_min = { -64.0f, -64.0f, -64.0f, 1.0f };
  //volume_max = {  64.0f,  64.0f,  64.0f, 1.0f };

  //volume_min = { -128.0f, -128.0f, -128.0f, 1.0f };
  //volume_max = {  128.0f,  128.0f,  128.0f, 1.0f };

  //volume_min = { -256.0f, -256.0f, -256.0f, 1.0f };
  //volume_max = {  256.0f,  256.0f,  256.0f, 1.0f };

  //volume_min = { -512.0f, -512.0f, -512.0f, 1.0f };
  //volume_max = {  512.0f,  512.0f,  512.0f, 1.0f };

  //volume_min = { -128.0f, -64.0f, -128.0f, 1.0f };
  //volume_max = {  128.0f,  64.0f,  128.0f, 1.0f };

  //volume_min = { -256.0f, -64.0f, -256.0f, 1.0f };
  //volume_max = {  256.0f,  64.0f,  256.0f, 1.0f };

  m_params.volumemin = volume_min;
  m_params.volumemax = volume_max;

  // velkost jednej bunky gridu
  //m_params.cell_size.s[0] = (2.0f * m_params.smoothradius) / m_params.simscale;
  //m_params.cell_size.s[1] = (2.0f * m_params.smoothradius) / m_params.simscale;
  //m_params.cell_size.s[2] = (2.0f * m_params.smoothradius) / m_params.simscale;

  //m_params.cell_size.s[0] = (8.0f * m_params.smoothradius) / m_params.simscale;
  //m_params.cell_size.s[1] = (8.0f * m_params.smoothradius) / m_params.simscale;
  //m_params.cell_size.s[2] = (8.0f * m_params.smoothradius) / m_params.simscale;

  m_params.cell_size.s[0] = (m_params.smoothradius) / m_params.simscale;
  m_params.cell_size.s[1] = (m_params.smoothradius) / m_params.simscale;
  m_params.cell_size.s[2] = (m_params.smoothradius) / m_params.simscale;

  //m_params.cell_size.s[0] = (0.5f * m_params.smoothradius) / m_params.simscale;
  //m_params.cell_size.s[1] = (0.5f * m_params.smoothradius) / m_params.simscale;
  //m_params.cell_size.s[2] = (0.5f * m_params.smoothradius) / m_params.simscale;

  //m_params.cell_size.s[0] = (0.25f * m_params.smoothradius) / m_params.simscale;
  //m_params.cell_size.s[1] = (0.25f * m_params.smoothradius) / m_params.simscale;
  //m_params.cell_size.s[2] = (0.25f * m_params.smoothradius) / m_params.simscale;

  m_params.cell_size.s[3] = 1.0f;

  // floatovy pocet buniek v danom objeme
  cl_float4 tmp;
  tmp.s[0] = (volume_max.s[0] - volume_min.s[0]) / m_params.cell_size.s[0];
  tmp.s[1] = (volume_max.s[1] - volume_min.s[1]) / m_params.cell_size.s[1];
  tmp.s[2] = (volume_max.s[2] - volume_min.s[2]) / m_params.cell_size.s[2];

  // integerovy pocet buniek v danom objeme
  m_params.grid_size.s[0] = ceil(tmp.s[0]);
  m_params.grid_size.s[1] = ceil(tmp.s[1]);
  m_params.grid_size.s[2] = ceil(tmp.s[2]);
  m_params.grid_size.s[3] = 1;

  // velkost korekcie, aby sedela floatova a integerova velkost pri pocitani na GPU
  m_params.grid_scale.s[0] = tmp.s[0] / (m_params.grid_size.s[0] * m_params.cell_size.s[0]);
  m_params.grid_scale.s[1] = tmp.s[1] / (m_params.grid_size.s[1] * m_params.cell_size.s[1]);
  m_params.grid_scale.s[2] = tmp.s[2] / (m_params.grid_size.s[2] * m_params.cell_size.s[2]);
  m_params.grid_scale.s[3] = 1.0f;

  // nastavenie ostatnych zavislych parametrov
  recalcDParams();
}


void SPHParams::setRoundOffCorrectGridSize(const tUInt4 & size, const tFloat correction)
{
#if 0
  // This is also a correct version, only here I experimented with a little safe boundary
  // to correct for the round off errors
  m_volume_min = { -33.6f, -33.6f, -33.6f, 1.0f };
  m_volume_max = {  33.6f,  33.6f,  33.6f, 1.0f };

  m_sim_params.volumemin = m_volume_min;
  m_sim_params.volumemax = m_volume_max;

  // velkost jednej bunky gridu
  m_params.cell_size.s[0] = (1.05f * m_params.smoothradius) / m_params.simscale;
  m_params.cell_size.s[1] = (1.05f * m_params.smoothradius) / m_params.simscale;
  m_params.cell_size.s[2] = (1.05f * m_params.smoothradius) / m_params.simscale;
  m_params.cell_size.s[3] = 1.0f;

  // integerovy pocet buniek v danom objeme
  m_params.grid_size.s[0] = (m_volume_max.s[0] - m_volume_min.s[0]) / m_params.cell_size.s[0];
  m_params.grid_size.s[1] = (m_volume_max.s[1] - m_volume_min.s[1]) / m_params.cell_size.s[1];
  m_params.grid_size.s[2] = (m_volume_max.s[2] - m_volume_min.s[2]) / m_params.cell_size.s[2];
  m_params.grid_size.s[3] = 1;
#else
  // pevne dane : grid_size + cell_size
  // dopocitane : volume
  // popis      : v tomto pripade su pevne dane grid_size a cell_size a tomuto sa prisposobi maximalny objem kvapaliny
  m_params.grid_size.s[0] = size.s[0];
  m_params.grid_size.s[1] = size.s[1];
  m_params.grid_size.s[2] = size.s[2];
  m_params.grid_size.s[3] = 1; //size.s[3];

  m_params.cell_size.s[0] = (correction * m_params.smoothradius) / m_params.simscale;
  m_params.cell_size.s[1] = (correction * m_params.smoothradius) / m_params.simscale;
  m_params.cell_size.s[2] = (correction * m_params.smoothradius) / m_params.simscale;
  m_params.cell_size.s[3] = 1.0f;

  m_params.volumemin.s[0] = -0.5f * m_params.cell_size.s[0] * m_params.grid_size.s[0];
  m_params.volumemin.s[1] = -0.5f * m_params.cell_size.s[1] * m_params.grid_size.s[1];
  m_params.volumemin.s[2] = -0.5f * m_params.cell_size.s[2] * m_params.grid_size.s[2];
  m_params.volumemin.s[3] =  1.0f;

  m_params.volumemax.s[0] =  0.5f * m_params.cell_size.s[0] * m_params.grid_size.s[0];
  m_params.volumemax.s[1] =  0.5f * m_params.cell_size.s[1] * m_params.grid_size.s[1];
  m_params.volumemax.s[2] =  0.5f * m_params.cell_size.s[2] * m_params.grid_size.s[2];
  m_params.volumemax.s[3] =  1.0f;
#endif

  // nastavenie ostatnych zavislych parametrov
  recalcDParams();
}


/**
 * Returns a reference to internal event which can be used to wait for the upload to finish
 */
boost::compute::event & SPHParams::upload(boost::compute::command_queue & queue)
{
  // ensure that previous frame has already finished
  if (m_ev)
  {
    m_ev.wait();
    //m_ev = boost::compute::event(); // uvolni event
  }

  cl_event ev;
  cl_int err = clEnqueueWriteBuffer(queue, m_params_buf.get(), CL_FALSE,
                                    0, sizeof(m_params), &m_params,
                                    0, nullptr, &ev);
  if (err != CL_SUCCESS)
  {
    WARNM("Failed to update simulation parameters: " <<
          boost::compute::opencl_error::to_string(err) <<
          "(" << err << ")");
  }

  m_ev = boost::compute::event(ev, false);  // do not retain the event since it was just created

  return m_ev;
}


void SPHParams::enqueuePrintGPUBuffer(boost::compute::command_queue & queue)
{
  //m_print_kernel.set_arg(0, m_params_buf);

  size_t gws = 1;
  cl_int err = clEnqueueNDRangeKernel(queue, m_print_kernel, 1,
                                      nullptr, &gws, nullptr,
                                      0, nullptr, nullptr);
  if (err != CL_SUCCESS)
  {
    WARNM("Failed to enqueue test simulation kernel: " <<
          boost::compute::opencl_error::to_string(err) <<
          "(" << err << ")");
  }
}


void SPHParams::dump(std::ostream & os) const
{
  os << "SPHParams:" << std::endl;
  os << "  volumemin      : " << m_params.volumemin      << std::endl;
  os << "  volumemax      : " << m_params.volumemax      << std::endl;
  os << "  seed           : " << m_params.seed           << std::endl;
  os << "  simscale       : " << m_params.simscale       << std::endl;
  os << "  radius2        : " << m_params.radius2        << std::endl;
  os << "  mass_polykern  : " << m_params.mass_polykern  << std::endl;
  os << "  restdensity    : " << m_params.restdensity    << std::endl;
  os << "  intstiffness   : " << m_params.intstiffness   << std::endl;
  os << "  numparticles   : " << m_params.numparticles   << std::endl;
  os << "  smoothradius   : " << m_params.smoothradius   << std::endl;
  os << "  vterm          : " << m_params.vterm          << std::endl;
  os << "  spikykern_half : " << m_params.spikykern_half << std::endl;
  os << "  slope          : " << m_params.slope          << std::endl;
  os << "  leftwave       : " << m_params.leftwave       << std::endl;
  os << "  rightwave      : " << m_params.rightwave      << std::endl;
  os << "  deltatime      : " << m_params.deltatime      << std::endl;
  os << "  limit          : " << m_params.limit          << std::endl;
  os << "  extstiffness   : " << m_params.extstiffness   << std::endl;
  os << "  extdamping     : " << m_params.extdamping     << std::endl;
  os << "  radius         : " << m_params.radius         << std::endl;
  os << "  mass           : " << m_params.mass           << std::endl;
  os << "  time           : " << m_params.time           << std::endl;
  os << "  flags          : " << m_params.flags          << std::endl;
  os << "  cell_size      : " << m_params.cell_size      << std::endl;
  os << "  grid_size      : " << m_params.grid_size      << std::endl;
  os << "  grid_scale     : " << m_params.grid_scale     << std::endl;
  os << "  top_face       : " << m_params.top_face       << std::endl;
  os << "  bottom_face    : " << m_params.bottom_face    << std::endl;
  os << "  front_face     : " << m_params.front_face     << std::endl;
  os << "  back_face      : " << m_params.back_face      << std::endl;
  os << "  left_face      : " << m_params.left_face      << std::endl;
  os << "  right_face     : " << m_params.right_face     << std::endl;
}


bool SPHParams::parseCmdArgs(const QStringList & args)
{
  int i = 0;
  bool ok = true;

  while (i < args.size())
  {
    const QString & arg = args[i];
    if (arg == "-particle-count")
    {
      ++i;
      tUInt cnt = args[i].toUInt(&ok);
      if (!ok)
      {
        ERRORM("Given particle count is not a number");
        return false;
      }

      m_params.numparticles = cnt;
    }
    else if (arg == "-grid-size")
    {
      ++i;

      tUInt grid_w = 0;
      tUInt grid_h = 0;
      tUInt grid_d = 0;

      QVector<QStringRef> grid_sizes(args[i].splitRef(QChar('x')));
      if (grid_sizes.size() == 1)
      {
        grid_w = grid_sizes[0].toUInt(&ok);
        if (!ok)
        {
          ERRORM("Grid size is invalid");
          return false;
        }
        grid_h = grid_d = grid_w;
      }
      else if (grid_sizes.size() == 3)
      {
        ok = true;
        bool ok2 = true;
        grid_w = grid_sizes[0].toUInt(&ok2); ok = ok && ok2;
        grid_h = grid_sizes[1].toUInt(&ok2); ok = ok && ok2;
        grid_d = grid_sizes[2].toUInt(&ok2); ok = ok && ok2;
        if (!ok)
        {
          ERRORM("Some of the given grid sizes are not valid numbers");
          return false;
        }
      }
      else
      {
        ERRORM("Need grid dimensions in format: '-grid-size WIDTHxHEIGHTxDEPTH'"
               " or in format: '-grid-size SIZE'");
        return false;
      }

      setGridSize(grid_w, grid_h, grid_d);
    }

    ++i;
  }

  return true;
}


void SPHParams::recalcDParams(void)
{
  //
  // This function recomputes the D parameter of every face (the plane in which it lies)
  // of the bounding volume of the fluid.
  // The formula for the D parameter was derived from parametric formula of plane:
  //
  // D = A * x + B * y + C * z
  //

  m_params.top_face.s[3]    = m_params.volumemax.s[1] * m_params.top_face.s[1];
  m_params.bottom_face.s[3] = m_params.volumemin.s[1] * m_params.bottom_face.s[1];
  m_params.front_face.s[3]  = m_params.volumemax.s[2] * m_params.front_face.s[2];
  m_params.back_face.s[3]   = m_params.volumemin.s[2] * m_params.back_face.s[2];
  m_params.left_face.s[3]   = m_params.volumemin.s[0] * m_params.left_face.s[0];
  m_params.right_face.s[3]  = m_params.volumemax.s[0] * m_params.right_face.s[0];
}


#if 1
tFloat4 SPHParams::rotateFaceNormal(const QMatrix4x4 & mat, const tFloat4 & face, const tFloat d)
{
  QVector3D v(face.s[0], face.s[1], face.s[2]);
  v = mat * v;
  v.normalize();
  return { v.x(), v.y(), v.z(), d };
}
#else
void SPHParams::rotateFaceNormal(const QMatrix4x4 & mat, tFloat4 & face)
{
  QVector3D v(face.s[0], face.s[1], face.s[2]);
  v = mat * v;
  face.s[0] = v.x();
  face.s[1] = v.y();
  face.s[2] = v.z();
}
#endif


void SPHParams::init(const boost::compute::context & ctx)
{
  try
  {
    // build opencl program
#ifdef DEBUG_SOURCE
    QFile f(":/src/opencl/sph_params.cl");
    if (!f.open(QFile::ReadOnly)) return false;

    m_prog = boost::compute::program::create_with_source(f.readAll().toStdString(), ctx);
    m_prog.build("-I../DIP/src/core -DDEBUG_SOURCE");
#else
    QFile f_cl(":/src/opencl/sph_params.cl");
    QFile f_h(":/src/core/sph_ocl_common.h");
    if ((!f_cl.open(QFile::ReadOnly)) || (!f_h.open(QFile::ReadOnly)))
    {
      ERRORM("Failed to open OpenCL program source files");
      throw std::runtime_error("Failed to open either " +
                               f_cl.fileName().toStdString() + " or " +
                               f_h.fileName().toStdString());
    }

    QByteArray src(f_h.readAll());
    src.append(f_cl.readAll());

    m_prog = boost::compute::program::create_with_source(src.toStdString(), ctx);
    m_prog.build();
#endif

    // print build log in case there are any warnings
    std::string log(m_prog.build_log());
    if (!log.empty())
    {
      WARNM("---------------------- Build log ---------------------------------");
      WARNM(log);
      WARNM("------------------- End of Build log -----------------------------");
    }

    // create kernels
    m_print_kernel = m_prog.create_kernel("sph_params_print");
    utils::ocl::printKernelInfo(m_print_kernel);

    // allocate GPU buffer
    allocGPUBuffer(ctx);

    // reset to default parameter settings
    resetToDefaults();

    m_print_kernel.set_arg(0, m_params_buf);
  }
  catch (const boost::compute::opencl_error & e)
  {
    ERRORM(e.what());
    if (e.error_code() == CL_BUILD_PROGRAM_FAILURE) ERRORM(m_prog.build_log());
    throw;
  }
}
