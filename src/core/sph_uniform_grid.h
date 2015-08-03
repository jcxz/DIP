#ifndef SPH_UNIFORM_GRID_H
#define SPH_UNIFORM_GRID_H

#include "core/fluid_particle_system.h"
#include "core/sph_params.h"

#include <boost/compute/program.hpp>
#include <boost/compute/kernel.hpp>
#include <boost/compute/buffer.hpp>
#include <boost/compute/container/vector.hpp>

//#define SIM_DEBUG

namespace test { class SPHUniformGridTestAdapter; }


class SPHUniformGrid : public FluidParticleSystem
{
  public:
    SPHUniformGrid(const boost::compute::context & ctx, SPHParams *params)
      : FluidParticleSystem(ctx, params)
      , m_sph_prog()
      , m_sph_reset_kernel()
      , m_sph_compute_step_kernel()
      , m_sph_compute_force_kernel()
      , m_sph_compute_pressure_kernel()
      , m_sph_calc_hash_kernel()
      , m_sph_reorder_kernel()
      , m_sph_calc_occupancy_kernel()
      , m_sph_count_grid_interactions_kernel()
#ifdef SIM_DEBUG
      , m_sph_count_interactions_kernel()
      , m_sph_calc_interaction_stats_kernel()
#endif
      , m_active_buf(0)
      , m_particle_pos_sorted_buf()
      , m_velocity_buf()
      , m_pressure_buf()
      , m_density_buf()
      , m_force_buf()
      , m_prev_velocity_buf()
      , m_part_hash_buf()
      , m_part_index_buf()
      , m_cell_starts_buf()
      , m_cell_ends_buf()
      , m_cell_part_counts_buf()
      , m_grid_stats_ratios_buf()
      , m_grid_radius_cnts_buf()
#ifdef SIM_DEBUG
      , m_brute_force_radius_cnts_buf()
      , m_grid_stats_diffs_buf()
#endif
      //, m_compute_sim_stats(true)
      , m_compute_sim_stats(false)
    {
      // initialize OpenCL context, compile kernels
      if (!init())
      {
        throw std::runtime_error("Failed to construct SPHUniformGrid: OpenCL initialization failed");
      }
    }

    const boost::compute::buffer & cellStartsBuffer(void) const { return m_cell_starts_buf; }
    const boost::compute::buffer & cellEndsBuffer(void) const { return m_cell_ends_buf; }

    // This function rotates the simulation box, so that the fluid can
    // react to the bounding volume being rotated
    void rotate(const QQuaternion & quat);

  protected:
    virtual bool reset_impl(const io::VoxelMesh * /* mesh */) override;
    virtual void update_impl(float time_step) override;

  private:
    // initializes the OpenCL program and kernel for SPH simulation
    bool init(void);

    // simulation helpers
    void resetGrid(void);

    // buffer swapping
    int activeBuf(void) const { return m_active_buf; }
    int inactiveBuf(void) const { return 1 - m_active_buf; }

    void swapActiveBuffers(void)
    {
      m_particle_pos_buf.swap(m_particle_pos_sorted_buf);
      m_active_buf = 1 - m_active_buf;
    }

    // statistics and debugging
    template <typename T>
    void calcAndPrintStats(const boost::compute::vector<T> & buf, size_t n,
                           const std::string & msg = std::string());
    void calcAndPrintGridOccupancy(void);
    void calcAndPrintInteractions(void);

    // printing
    template <typename T, typename BufType>
    void dumpBuffer(const BufType & buf, std::ostream & os)
    {
      return utils::ocl::dumpBuffer<T, BufType>(m_queue, buf, os,
                                                m_sim_params->particleCount(),
                                                10);
    }

    void dumpUniformGrid(std::ostream & os);

  private:
    // OpenCL programs
    boost::compute::program m_sph_prog;

    // OpenCL kernels
    boost::compute::kernel m_sph_reset_kernel;                    // a reset kernel for the SPH simulation
    boost::compute::kernel m_sph_compute_step_kernel;             // a kernel to compute a single SPH step
    boost::compute::kernel m_sph_compute_force_kernel;            // kernel for computing forces
    boost::compute::kernel m_sph_compute_pressure_kernel;         // kernel for computing the pressure inside of the fluid
    boost::compute::kernel m_sph_calc_hash_kernel;                // a kernel for computing particle hash
    boost::compute::kernel m_sph_reorder_kernel;                  // a kernel to reorder particles sorted according to hash
    boost::compute::kernel m_sph_calc_occupancy_kernel;           // a helper kernel for collecting statistics on grid occupancy
    boost::compute::kernel m_sph_count_grid_interactions_kernel;  // a kernel to count the number of interacting particles using uniform grid
#ifdef SIM_DEBUG
    boost::compute::kernel m_sph_count_interactions_kernel;       // a kernel to count the number of interacting particles using brute force
    boost::compute::kernel m_sph_calc_interaction_stats_kernel;   // a kernel to process data gathered in count_interactions
                                                                  // and calc_interaction_stats
#endif

    // buffers for SPH simulation
    int m_active_buf;
    utils::ocl::GLBuffer m_particle_pos_sorted_buf;
    boost::compute::buffer m_velocity_buf[2];
    boost::compute::buffer m_pressure_buf[2];
    boost::compute::buffer m_density_buf[2];
    boost::compute::buffer m_force_buf[2];
    boost::compute::buffer m_prev_velocity_buf[2];

    // vectors for sorting the particles pointed to from uniform grid
    boost::compute::vector<cl_uint> m_part_hash_buf;   // this buffer stores particle hashes
    boost::compute::vector<cl_uint> m_part_index_buf;  // this buffer stores sorted particle indices

    // these buffers basically define the grid structure
    boost::compute::buffer m_cell_starts_buf;               // a buffer to store cell starting indices
    boost::compute::buffer m_cell_ends_buf;                 // a buffer to store cell ending indices

    // buffers used for statistics and debugging
    boost::compute::vector<cl_uint> m_cell_part_counts_buf;        // a buffer to store particle counts in cells
    boost::compute::vector<cl_float> m_grid_stats_ratios_buf;      //
    boost::compute::vector<cl_uint> m_grid_radius_cnts_buf;        // a buffer to store the count of particles interacting in grid
#ifdef SIM_DEBUG
    boost::compute::vector<cl_uint> m_brute_force_radius_cnts_buf; // a buffer to store the count of interacting particles
                                                                   // when using brute force method
    boost::compute::vector<cl_uint> m_grid_stats_diffs_buf;        //
#endif
    bool m_compute_sim_stats;                                      // whether gathering simulation statstics is enabled

    // give the test adapter class access to the internals
    friend class test::SPHUniformGridTestAdapter;
};

#endif // SPH_UNIFORM_GRID_H
