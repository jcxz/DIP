#include "core/particle_system_internal.h"
#include "core/sph_uniform_grid.h"
#include "utils/macros.h"
#include "utils/debug.h"
#include "utils/ocl.h"
#include "io/voxel_mesh.h"

#include <QFile>
#include <ctime>
#include <boost/compute/algorithm/sort_by_key.hpp>
#include <boost/compute/algorithm/sort.hpp>
#include <boost/compute/algorithm/reduce.hpp>
#include <boost/compute/algorithm/count.hpp>

//#define DEBUG_PARAMS
#define BLOCK_SIZE 256




void SPHUniformGrid::rotate(const QQuaternion & quat)
{
  QMatrix4x4 m;
  m.rotate(quat);

  //m = m.inverted();

  //m.rotate(45, 0.0f, 0.0f, 1.0f);
  //m.rotate(45, 0.0f, 1.0f, 0.0f);
  //m.rotate(45, 1.0f, 0.0f, 0.0f);
  //m.rotate(180, 0.0f, 0.0f, 1.0f);

  //m.rotate(10, 0.0f, 0.0f, -1.0f);
  //m.rotate(10, 0.0f, -1.0f, 0.0f);
  //m.rotate(10, -1.0f, 0.0f, 0.0f);

  m_sim_params->rotateTopFace(m);
  m_sim_params->rotateBottomFace(m);
  m_sim_params->rotateFrontFace(m);
  m_sim_params->rotateBackFace(m);
  m_sim_params->rotateLeftFace(m);
  m_sim_params->rotateRightFace(m);
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Simulation initialization

bool SPHUniformGrid::reset_impl(const io::VoxelMesh *mesh)
{
  if (mesh) m_sim_params->setParticleCount(mesh->count());

  // reset simulation parameters
  m_sim_params->setTime(0.0f);
  m_sim_params->clearEffectFlags();

  // total number of cells in the uniform grid structure
  size_t grid_cell_count = m_sim_params->gridCellCount();
  size_t part_cnt = m_sim_params->particleCount();

  INFOM("particle radius       : " << m_sim_params->smoothRadius());
  INFOM("grid_size             : " << m_sim_params->gridSize());
  INFOM("cell_size             : " << m_sim_params->cellSize());
  INFOM("volume min            : " << m_sim_params->volumeMin());
  INFOM("volume max            : " << m_sim_params->volumeMax());
  INFOM("cell_count            : " << grid_cell_count);

  // calculate optimal global and local work sizes and buffer padding
  size_t block_size = BLOCK_SIZE;
  size_t grid_size = (part_cnt + block_size - 1) / block_size;
  size_t local_work_size = block_size;
  size_t global_work_size = grid_size * block_size;
  size_t padding = global_work_size - part_cnt;

  INFOM("Global Work Size      : " << global_work_size);
  INFOM("Local Work Size       : " << local_work_size);
  INFOM("Memory Object padding : " << padding);
  INFOM(*m_sim_params);

  // allocate buffers on GPU
  if (!m_particle_pos_buf.bufferData(nullptr, (part_cnt + padding) * sizeof(cl_float4),
                                     utils::ocl::GLBuffer::WRITE_ONLY))
  {
    ERRORM("SPHUniformGrid: Failed to initialize position GLBuffer");
    return false;
  }

  if (!m_particle_pos_sorted_buf.bufferData(nullptr, (part_cnt + padding) * sizeof(cl_float4),
                                            utils::ocl::GLBuffer::WRITE_ONLY))
  {
    ERRORM("SPHUniformGrid: Failed to initialize position GLBuffer");
    return false;
  }

  try
  {
    m_velocity_buf[0] = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float4));
    m_prev_velocity_buf[0] = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float4));
    m_pressure_buf[0] = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float));
    m_density_buf[0] = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float));
    m_force_buf[0] = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float4));

    m_velocity_buf[1] = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float4));
    m_prev_velocity_buf[1] = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float4));
    m_pressure_buf[1] = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float));
    m_density_buf[1] = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float));
    m_force_buf[1] = boost::compute::buffer(m_cl_ctx, (part_cnt + padding) * sizeof(cl_float4));

    //m_part_hash_buf.resize(part_cnt + padding);
    //m_part_index_buf.resize(part_cnt + padding);
    m_part_hash_buf = boost::compute::vector<cl_uint>(part_cnt + padding, m_cl_ctx);
    m_part_index_buf = boost::compute::vector<cl_uint>(part_cnt + padding, m_cl_ctx);

    // the uniform grid itself
    m_cell_starts_buf = boost::compute::buffer(m_cl_ctx, grid_cell_count * sizeof(cl_uint));
    m_cell_ends_buf = boost::compute::buffer(m_cl_ctx, grid_cell_count * sizeof(cl_uint));

    // buffers used for statistics and debugging
    m_cell_part_counts_buf = boost::compute::vector<cl_uint>(grid_cell_count, m_cl_ctx);
    m_grid_stats_ratios_buf = boost::compute::vector<cl_float>(part_cnt + padding, m_cl_ctx);
    m_grid_radius_cnts_buf = boost::compute::vector<cl_uint>(part_cnt + padding, m_cl_ctx);
#ifdef SIM_DEBUG
    m_brute_force_radius_cnts_buf = boost::compute::vector<cl_uint>(part_cnt + padding, m_cl_ctx);
    m_grid_stats_diffs_buf = boost::compute::vector<cl_uint>(part_cnt + padding, m_cl_ctx);
#endif
  }
  catch (const boost::compute::opencl_error & e)
  {
    ERRORM("SPHUniformGrid: An OpenCL error occured: " << e.what());
    return false;
  }
  catch (const std::exception & e)
  {
    ERRORM("SPHUniformGrid: An unexpected error occured during SPHUniformGrid initialization: "
           << e.what());
    return false;
  }

  // upload simulation parameters to gpu
  m_stats.addEvent(m_sim_params->upload(m_queue), "Update SPHUniformGrid params");

  // run the reset kernel to initialize particle data
#if 1
  cl_mem buffers[] = { m_particle_pos_buf.getCLID() };
  utils::ocl::GLSyncHandler sync(m_queue, sizeof(buffers) / sizeof(buffers[0]), buffers);
  if (!sync) return false;
#else
  utils::ocl::GLSyncHandler sync(m_queue, m_particle_pos_buf.getCLID());
  if (!sync) return false;
#endif

  if (mesh == nullptr)
  {
    // This can be optionally be solved by doing e.g.: m_velocity_buf[activeBuf()],
    // but this way it is more safe and the here code is called only once during the initialization,
    // so what is the big deal anyways
    return deviceCall(m_sph_reset_kernel, "sph_uniform_grid_reset",
                      global_work_size, local_work_size,
                      m_particle_pos_buf.getCLID(),
                      m_velocity_buf[0],
                      m_prev_velocity_buf[0],
                      m_pressure_buf[0],
                      m_density_buf[0],
                      m_force_buf[0],

                      m_particle_pos_sorted_buf.getCLID(),
                      m_velocity_buf[1],
                      m_prev_velocity_buf[1],
                      m_pressure_buf[1],
                      m_density_buf[1],
                      m_force_buf[1],

                      m_part_hash_buf.get_buffer(),
                      m_part_index_buf.get_buffer(),
                      m_sim_params->buffer());
  }
  else
  {
    // otherwise in case the initial particle positions are provided
    cl_float4 *pos        = utils::ocl::mapBuffer<cl_float4>(m_queue, m_particle_pos_buf, part_cnt, CL_MAP_WRITE);
    cl_float4 *pos_sorted = utils::ocl::mapBuffer<cl_float4>(m_queue, m_particle_pos_sorted_buf, part_cnt, CL_MAP_WRITE);
    const cl_float4 *mesh_data = mesh->data();

    if ((pos == nullptr) || (pos_sorted == nullptr))
    {
      ERRORM("Failed to map pos and pos_sorted buffers");
      utils::ocl::unmapBuffer(m_queue, m_particle_pos_sorted_buf, pos_sorted);
      utils::ocl::unmapBuffer(m_queue, m_particle_pos_buf, pos);
      return false;
    }

    cl_float4 bbox_min = m_sim_params->volumeMin();

    float transl_x = (m_sim_params->volumeWidth()  - mesh->width())  * 0.5f;
    float transl_y = (m_sim_params->volumeHeight() - mesh->height()) * 0.5f;
    float transl_z = (m_sim_params->volumeDepth()  - mesh->depth())  * 0.5f;

    //float transl_x = 0.0f;
    //float transl_y = 0.0f;
    //float transl_z = 0.0f;

    std::cout << "min_x=" << bbox_min.s[0] << ", min_y=" << bbox_min.s[1] << ", min_z=" << bbox_min.s[2] << std::endl;
    std::cout << "transl_x=" << transl_x << ", transl_y=" << transl_y << ", transl_z=" << transl_z << std::endl;

    for (int i = 0; i < mesh->count(); ++i)
    {
      cl_float4 p = mesh_data[i];
      //p.s[0] -= transl_x;
      //p.s[1] -= transl_y;
      //p.s[2] -= transl_z;
      //p.s[0] += transl_x;
      //p.s[1] += transl_y;
      //p.s[2] += transl_z;
      p.s[0] += bbox_min.s[0] + transl_x;
      p.s[1] += bbox_min.s[1] + transl_y;
      p.s[2] += bbox_min.s[2] + transl_z;
      pos[i] = p;
      pos_sorted[i] = p;
    }

    utils::ocl::unmapBuffer(m_queue, m_particle_pos_sorted_buf, pos_sorted);
    utils::ocl::unmapBuffer(m_queue, m_particle_pos_buf, pos);

    //utils::ocl::writeBuffer(m_queue, m_particle_pos_buf,        mesh->data(), part_cnt);
    //utils::ocl::writeBuffer(m_queue, m_particle_pos_sorted_buf, mesh->data(), part_cnt);

    utils::ocl::zeroBuffer<cl_float4>(m_queue, m_velocity_buf[0], part_cnt);
    utils::ocl::zeroBuffer<cl_float4>(m_queue, m_velocity_buf[1], part_cnt);
    utils::ocl::zeroBuffer<cl_float4>(m_queue, m_prev_velocity_buf[0], part_cnt);
    utils::ocl::zeroBuffer<cl_float4>(m_queue, m_prev_velocity_buf[1], part_cnt);
    utils::ocl::zeroBuffer<cl_float>(m_queue, m_pressure_buf[0], part_cnt);
    utils::ocl::zeroBuffer<cl_float>(m_queue, m_pressure_buf[1], part_cnt);
    utils::ocl::zeroBuffer<cl_float>(m_queue, m_density_buf[0], part_cnt);
    utils::ocl::zeroBuffer<cl_float>(m_queue, m_density_buf[1], part_cnt);
    utils::ocl::zeroBuffer<cl_float4>(m_queue, m_force_buf[0], part_cnt);
    utils::ocl::zeroBuffer<cl_float4>(m_queue, m_force_buf[1], part_cnt);
    return deviceCall(m_sph_calc_hash_kernel, "sph_uniform_grid_calc_hash",
                      global_work_size, local_work_size,
                      m_particle_pos_buf.getCLID(),
                      m_part_hash_buf.get_buffer(),
                      m_part_index_buf.get_buffer(),
                      m_sim_params->buffer());
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Main simulation code

void SPHUniformGrid::update_impl(float time_step)
{
  if (m_cl_timer.isReady()) m_cl_timer.start();

  // upload simulation parameters to gpu
  m_stats.addEvent(m_sim_params->upload(m_queue), "Update SPHUniformGrid params");

#ifdef DEBUG_PARAMS
  m_sim_params.enqueuePrintGPUBuffer(m_queue);
#endif

  // synchronise with OpenGL
  cl_mem buffers[] = {
    m_particle_pos_buf.getCLID(),
    m_particle_pos_sorted_buf.getCLID()
  };

  utils::ocl::GLSyncHandler sync(m_queue, sizeof(buffers) / sizeof(buffers[0]), buffers);
  if (!sync) return;

  //if (!acquireGL<2>({ m_particle_pos_buf, m_particle_pos_sorted_buf })) return;
  //
  //utils::misc::ScopedGuard guard([this] {
  //  releaseGL<2>({ m_particle_pos_buf, m_particle_pos_sorted_buf });
  //});

  //if (!acquireGL<2>({ &m_particle_pos_buf, &m_particle_pos_sorted_buf }))
  //  return;

  //std::array<const utils::ocl::GLBuffer *, 2> a { &m_particle_pos_buf, &m_particle_pos_sorted_buf };
  //if (!acquireGL<2>(a)) return;

  // calculate global work size based on local work size
  size_t part_cnt = m_sim_params->particleCount();
  size_t block_size = BLOCK_SIZE;
  size_t grid_size = (part_cnt + block_size - 1) / block_size;
  size_t local_work_size = block_size;
  size_t global_work_size = grid_size * block_size;

  // compute grid hashes
  //deviceCall(m_sph_calc_hash_kernel, "sph_uniform_grid_calc_hash",
  //           global_work_size, local_work_size,
  //           m_particle_pos_buf.getCLID(),
  //           m_part_hash_buf.get_buffer(),
  //           m_part_index_buf.get_buffer(),
  //           m_sim_params->buffer());

  // reset grid cells
  resetGrid();

  // sort
  boost::compute::sort_by_key(m_part_hash_buf.begin(),
                              m_part_hash_buf.begin() + part_cnt,
                              m_part_index_buf.begin(),
                              m_queue);

  // compute cell starts and cell ends and reorder data
  deviceCall(m_sph_reorder_kernel, "sph_uniform_grid_reorder",
             global_work_size, local_work_size,
             m_part_hash_buf.get_buffer(),
             m_part_index_buf.get_buffer(),
             m_cell_starts_buf,
             m_cell_ends_buf,
             m_particle_pos_buf.getCLID(),
             //m_density_buf[activeBuf()],
             //m_pressure_buf[activeBuf()],
             //m_force_buf[activeBuf()],
             m_velocity_buf[activeBuf()],
             m_prev_velocity_buf[activeBuf()],
             m_particle_pos_sorted_buf.getCLID(),
             //m_density_buf[inactiveBuf()],
             //m_pressure_buf[inactiveBuf()],
             //m_force_buf[inactiveBuf()],
             m_velocity_buf[inactiveBuf()],
             m_prev_velocity_buf[inactiveBuf()],
             m_sim_params->buffer(),
             utils::ocl::KernelArgs::LocalArg<cl_uint>(local_work_size + 1));

  swapActiveBuffers();

  // compute stats
  if (m_compute_sim_stats)
  {
    calcAndPrintGridOccupancy();
    calcAndPrintInteractions();
  }

  // compute pressure
  deviceCall(m_sph_compute_pressure_kernel, "sph_uniform_grid_compute_pressure",
             global_work_size, local_work_size,
             m_cell_starts_buf,
             m_cell_ends_buf,
             m_particle_pos_buf.getCLID(),
             m_density_buf[activeBuf()],
             m_pressure_buf[activeBuf()],
             m_sim_params->buffer());

  // compute force
  deviceCall(m_sph_compute_force_kernel, "sph_uniform_grid_compute_force",
             global_work_size, local_work_size,
             m_cell_starts_buf,
             m_cell_ends_buf,
             m_particle_pos_buf.getCLID(),
             m_density_buf[activeBuf()],
             m_pressure_buf[activeBuf()],
             m_force_buf[activeBuf()],
             m_velocity_buf[activeBuf()],
             m_sim_params->buffer());

  // integrate
  deviceCall(m_sph_compute_step_kernel, "sph_uniform_grid_compute_step",
             global_work_size, local_work_size,
             m_particle_pos_buf.getCLID(),
             m_force_buf[activeBuf()],
             m_velocity_buf[activeBuf()],
             m_prev_velocity_buf[activeBuf()],
             m_part_hash_buf.get_buffer(),
             m_part_index_buf.get_buffer(),
             m_sim_params->buffer());

  // advance simulation time
  m_sim_params->advanceSimulation(time_step);  // 3.0f;

  if (m_cl_timer.isReady()) m_cl_timer.stop();

  //std::cout << "=================================================================================" << std::endl;
  //std::cout << "Time           : " << m_sim_params->time() << std::endl;
  //std::cout << "Particle count : " << part_cnt             << std::endl;
  //dumpUniformGrid(std::cout);
  //std::cout << "=================================================================================" << std::endl;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Initialization (program compilation, kernel creation, ...)

bool SPHUniformGrid::init(void)
{
  try
  {
    // build opencl program
#ifdef DEBUG_SOURCE
    QFile f(":/src/opencl/sph_uniform_grid.cl");
    if (!f.open(QFile::ReadOnly)) return false;

    m_sph_prog = boost::compute::program::create_with_source(f.readAll().toStdString(), m_cl_ctx);
    m_sph_prog.build("-I../DIP/src/core -DDEBUG_SOURCE");
#else
    QFile f_cl(":/src/opencl/sph_uniform_grid.cl");
    QFile f_h(":/src/core/sph_ocl_common.h");
    if ((!f_cl.open(QFile::ReadOnly)) || (!f_h.open(QFile::ReadOnly)))
    {
      ERRORM("Failed to open OpenCL program source files");
      return false;
    }

    QByteArray src(f_h.readAll());
    src.append(f_cl.readAll());

    m_sph_prog = boost::compute::program::create_with_source(src.toStdString(), m_cl_ctx);
    m_sph_prog.build();
#endif

    // print build log in case there are any warnings
    std::string log(m_sph_prog.build_log());
    if (!log.empty())
    {
      WARNM("---------------------- Build log ---------------------------------");
      WARNM(log);
      WARNM("------------------- End of Build log -----------------------------");
    }

    // create kernels
    m_sph_reset_kernel = m_sph_prog.create_kernel("sph_uniform_grid_reset");
    utils::ocl::printKernelInfo(m_sph_reset_kernel);

    m_sph_compute_step_kernel = m_sph_prog.create_kernel("sph_uniform_grid_compute_step");
    utils::ocl::printKernelInfo(m_sph_compute_step_kernel);

    m_sph_compute_force_kernel = m_sph_prog.create_kernel("sph_uniform_grid_compute_force");
    utils::ocl::printKernelInfo(m_sph_compute_force_kernel);

    m_sph_compute_pressure_kernel = m_sph_prog.create_kernel("sph_uniform_grid_compute_pressure");
    utils::ocl::printKernelInfo(m_sph_compute_pressure_kernel);

    m_sph_calc_hash_kernel = m_sph_prog.create_kernel("sph_uniform_grid_calc_hash");
    utils::ocl::printKernelInfo(m_sph_calc_hash_kernel);

    m_sph_reorder_kernel = m_sph_prog.create_kernel("sph_uniform_grid_reorder");
    utils::ocl::printKernelInfo(m_sph_reorder_kernel);

    m_sph_calc_occupancy_kernel = m_sph_prog.create_kernel("sph_uniform_grid_calc_occupancy");
    utils::ocl::printKernelInfo(m_sph_calc_occupancy_kernel);

    m_sph_count_grid_interactions_kernel = m_sph_prog.create_kernel("sph_uniform_grid_count_grid_interactions");
    utils::ocl::printKernelInfo(m_sph_count_grid_interactions_kernel);

#ifdef SIM_DEBUG
    m_sph_count_interactions_kernel = m_sph_prog.create_kernel("sph_uniform_grid_count_interactions");
    utils::ocl::printKernelInfo(m_sph_count_interactions_kernel);

    m_sph_calc_interaction_stats_kernel = m_sph_prog.create_kernel("sph_uniform_grid_calc_interaction_stats");
    utils::ocl::printKernelInfo(m_sph_calc_interaction_stats_kernel);
#endif
  }
  catch (const boost::compute::opencl_error & e)
  {
    ERRORM("SPHUniformGrid: An OpenCL error occured: " << e.what());
    if (e.error_code() == CL_BUILD_PROGRAM_FAILURE) ERRORM(m_sph_prog.build_log());
    return false;
  }
  catch (const std::exception & e)
  {
    ERRORM("SPHUniformGrid: An unexpected error occured during SPHUniformGrid initialization: "
           << e.what());
    return false;
  }

  m_particle_pos_sorted_buf.setCLContext(m_cl_ctx);

  return true;
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Helper functions for the main simulation code (i.e update_impl function)

void SPHUniformGrid::resetGrid(void)
{
  const cl_uint pattern = INVALID_UNIFORM_GRID_CELL_VALUE;

  size_t grid_cell_count = m_sim_params->gridCellCount();

  cl_int err = clEnqueueFillBuffer(m_queue, m_cell_starts_buf.get(),
                                   &pattern, sizeof(pattern),
                                   0, (grid_cell_count) * sizeof(cl_uint),
                                   0, nullptr, m_stats.event("Clear cell_starts buffer"));
  if (err != CL_SUCCESS)
  {
    WARNM("SPHUniformGrid: Failed to reset cell_starts buffer: " <<
          boost::compute::opencl_error::to_string(err) <<
          "(" << err << ")");
  }
}


////////////////////////////////////////////////////////////////////////////////////////////////////
// Debugging and statistics

template <typename T>
void SPHUniformGrid::calcAndPrintStats(const boost::compute::vector<T> & buf,
                                       size_t n, const std::string & msg)
{
#if 1
  T sum = 0;
  boost::compute::reduce(buf.begin(), buf.begin() + n, &sum, m_queue);

  T min = std::numeric_limits<T>::max();
  boost::compute::reduce(buf.begin(), buf.begin() + n, &min,
                         boost::compute::min<T>(), m_queue);

  T max = std::numeric_limits<T>::min();
  boost::compute::reduce(buf.begin(), buf.begin() + n, &max,
                         boost::compute::max<T>(), m_queue);

  boost::compute::sort(buf.begin(), buf.begin() + n, m_queue);
  T median = buf[buf.size() / 2];

  cl_uint zero_cnt = 0;
  zero_cnt = boost::compute::count(buf.begin(), buf.begin() + n, zero_cnt, m_queue);
#else
  boost::compute::sort(buf.begin(), buf.begin() + n, m_queue);

  T sum = 0;
  boost::compute::reduce(buf.begin(), buf.begin() + n, &sum, m_queue);

  cl_uint zero_cnt = 0;
  zero_cnt = boost::compute::count(buf.begin(), buf.begin() + n, zero_cnt, m_queue);

  T min = buf[0];
  T max = buf[n - 1];
  T median = buf[n / 2];
#endif

  float avg = (float(sum) / float(n - zero_cnt));

  INFOM(msg << ": avg=" << avg << ", min=" << min << ", max=" <<
        max << ", median=" << median << ", zero_cnt=" << zero_cnt);
}


void SPHUniformGrid::calcAndPrintGridOccupancy(void)
{
  size_t grid_cell_count = m_sim_params->gridCellCount();

  deviceCall2(m_sph_calc_occupancy_kernel,
              "sph_uniform_grid_calc_occupancy",
              grid_cell_count,
              m_cell_starts_buf,
              m_cell_ends_buf,
              m_cell_part_counts_buf.get_buffer());

  calcAndPrintStats(m_cell_part_counts_buf,
                    grid_cell_count,
                    "Grid Occupancy stats");
}


void SPHUniformGrid::calcAndPrintInteractions(void)
{
  using utils::ocl::mapBuffer;
  using utils::ocl::unmapBuffer;

  size_t part_cnt = m_sim_params->particleCount();

  // count particle interactions when using the uniform grid method
  deviceCall2(m_sph_count_grid_interactions_kernel,
              "sph_uniform_grid_count_grid_interactions",
              part_cnt,
              m_cell_starts_buf,
              m_cell_ends_buf,
              m_particle_pos_buf.getCLID(),
              m_grid_radius_cnts_buf.get_buffer(),
              m_grid_stats_ratios_buf.get_buffer(),
              m_sim_params->buffer());

  calcAndPrintStats(m_grid_stats_ratios_buf,
                    part_cnt,
                    "Grid Interaction stats");

#ifdef SIM_DEBUG
  // count particle interactions when using the bruteforce method
  deviceCall2(m_sph_count_interactions_kernel,
              "sph_uniform_grid_count_interactions",
              part_cnt,
              m_particle_pos_buf.getCLID(),
              m_brute_force_radius_cnts_buf.get_buffer(),
              m_sim_params->buffer());

  //
  deviceCall2(m_sph_calc_interaction_stats_kernel,
              "sph_uniform_grid_calc_interaction_stats",
              part_cnt,
              m_grid_radius_cnts_buf.get_buffer(),
              m_brute_force_radius_cnts_buf.get_buffer(),
              m_grid_stats_diffs_buf.get_buffer());

  cl_uint cnt_non_zero = boost::compute::count_if(m_grid_stats_diffs_buf.begin(),
                                                  m_grid_stats_diffs_buf.begin() + part_cnt,
                                                  boost::compute::_1 != 0, m_queue);

  cl_uint sum = 0;
  boost::compute::reduce(m_grid_stats_diffs_buf.begin(),
                         m_grid_stats_diffs_buf.begin() + part_cnt,
                         &sum, m_queue);
  if (sum != 0)
  {
    INFOM("!!!!!!!!!!!!!!!!!!!! Difference between grid and brute force: "
          << (int) sum << ", nonzero count:" << cnt_non_zero);

    cl_uint *p = mapBuffer<cl_uint>(m_queue, m_grid_stats_diffs_buf.get_buffer(), part_cnt);
    cl_uint *grid_p = mapBuffer<cl_uint>(m_queue, m_grid_radius_cnts_buf.get_buffer(), part_cnt);
    cl_uint *brute_force_p = mapBuffer<cl_uint>(m_queue, m_brute_force_radius_cnts_buf.get_buffer(),
                                                part_cnt);

    for (size_t i = 0; i < part_cnt; ++i)
    {
      if (p[i] != 0)
      {
        std::cerr << i << ": " << grid_p[i] << "<->"
                  << brute_force_p[i] << ": " << (int) p[i]
                  << std::endl;
      }
    }

    unmapBuffer(m_queue, m_brute_force_radius_cnts_buf.get_buffer(), brute_force_p);
    unmapBuffer(m_queue, m_grid_radius_cnts_buf.get_buffer(), grid_p);
    unmapBuffer(m_queue, m_grid_stats_diffs_buf.get_buffer(), p);
  }
  else
  {
    INFOM("Difference between grid and brute force: " << (int) sum);
  }
#endif
}


void SPHUniformGrid::dumpUniformGrid(std::ostream & os)
{
  using utils::ocl::mapBuffer;

  size_t grid_w = m_sim_params->gridWidth();
  size_t grid_h = m_sim_params->gridHeight();
  size_t grid_d = m_sim_params->gridDepth();
  size_t grid_cell_count = grid_w * grid_h * grid_d;
  size_t part_cnt = m_sim_params->particleCount();

  cl_uint *p_cell_starts = mapBuffer<cl_uint>(m_queue, m_cell_starts_buf, grid_cell_count);
  cl_uint *p_cell_ends = mapBuffer<cl_uint>(m_queue, m_cell_ends_buf, grid_cell_count);
  //cl_float4 *p_pos = mapBuffer<cl_float4>(m_queue, m_particle_pos_buf, part_cnt);
  //cl_float *p_density = mapBuffer<cl_float>(m_queue, m_density_buf[activeBuf()], part_cnt);
  //cl_float *p_pressure = mapBuffer<cl_float>(m_queue, m_pressure_buf[activeBuf()], part_cnt);
  //cl_float4 *p_forces = mapBuffer<cl_float4>(m_queue, m_force_buf[activeBuf()], part_cnt);

  os << "=================================================================================" << std::endl;

  size_t cells_cnt = 0;
  size_t part_sum = 0;
  for (size_t k = 0; k < grid_d; ++k)
  {
    for (size_t j = 0; j < grid_h; ++j)
    {
      for (size_t i = 0; i < grid_w; ++i)
      {
        int idx = i + j * grid_w + k * grid_w * grid_h;
        cl_int start = p_cell_starts[idx];
        cl_int end = p_cell_ends[idx];
        if (start != (cl_int) (INVALID_UNIFORM_GRID_CELL_VALUE))
        {
          ++cells_cnt;

          cl_int diff = end - start;
          part_sum += diff;
          if (diff <= 0)
          {
            os << "!!!!!!!!!!!! [" << i << "," << j << "," << k << "]: ("
               << start << "," << end << ")\n";
            continue;
          }

          if (diff > 1)  os << "*** ";
#if 1
          os << "[" << i << "," << j << "," << k << "]: (" << start << "," << end << ")\n";
#else
          // print position
          os << "[" << i << "," << j << "," << k << "]: (" << start << "," << end << ") -> ";
          utils::debug::printArray1D(p_pos, diff, os);
          os << "\n";

          // print density
          os << "                  -> ";
          utils::debug::printArray1D(p_density, diff, os);
          os << "\n";

          // print pressure
          os << "                  -> ";
          utils::debug::printArray1D(p_pressure, diff, os);
          os << "\n";

          // print forces
          os << "                  -> ";
          utils::debug::printArray1D(p_forces, diff, os);
          os << "\n";
#endif
        }
      }
    }
  }

  os << "Number of occupied cells: " << cells_cnt << std::endl;
  if (part_sum != part_cnt) os << "@@@@@@@@@@@@@@";
  os << "Particle checksum: " << part_sum << std::endl;

  //utils::ocl::unmapBuffer(m_queue, m_force_buf[activeBuf()], p_forces);
  //utils::ocl::unmapBuffer(m_queue, m_pressure_buf[activeBuf()], p_pressure);
  //utils::ocl::unmapBuffer(m_queue, m_density_buf[activeBuf()], p_density);
  //utils::ocl::unmapBuffer(m_queue, m_particle_pos_buf, p_pos);
  utils::ocl::unmapBuffer(m_queue, m_cell_ends_buf, p_cell_ends);
  utils::ocl::unmapBuffer(m_queue, m_cell_starts_buf, p_cell_starts);
}
