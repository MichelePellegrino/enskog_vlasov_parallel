/*! \file density.hpp
 *  \brief Header containing the class implementing density kernel
 */

#ifndef EV_DENSITIES_HPP
#define EV_DENSITIES_HPP

#include "matrix.hpp"
#include "types.hpp"
#include "utility.hpp"
#include "motherbase.hpp"

#include <cmath>

/*! \class DensityKernel
 *  \brief Class for density and reduced density computation
 *
 *  It stores the values of averaging weights, binned particles, number density
 *  (reduced and averaged). Moreover, it encapsulates the convolutioner for
 *  averages computation. Finally, defines map buffer to locate particles in and
 *  neighbouring a given cell.
 */
class DensityKernel : protected Motherbase
{

private:

  // Utilities (maybe reference would be better?)
  real_number reduce_factor;    /*!< Specific covolume (to reduce number density) */
  int ns_x, ns_y;               /*!< Half-number of weights in each direction     */

  ev_matrix::SlideMaskMatrix<real_number> stencil_x;  /*!< Stencil values along x to compute averaging weights  */
  ev_matrix::SlideMaskMatrix<real_number> stencil_y;  /*!< Stencil values along y to compute averaging weights  */
  ev_matrix::SlideMaskMatrix<real_number> weights;    /*!< Weights to average density values                    */

  real_number x_min, y_min;     /*!< Grid boundaries                              */
  int n_cutoff_x, n_cutoff_y;   /*!< Maximum cutoff (needed to reserve halo size) */

  // Problem: casting Eigen datatype -> solved!

  ev_matrix::MaskMatrix<int> n_part_cell;                         /*!< Number of particles for each cell            */
  ev_matrix::HaloMaskMatrix<real_number> num_dens_cell;           /*!< Number density value in each cell            */
  ev_matrix::HaloMaskMatrix<real_number> reduced_density;         /*!< Reduced density values in each cell          */
  ev_matrix::MaskMatrix<real_number> average_reduced_density;     /*!< Averaged reduced density values in each cell */

  ev_matrix::MatrixConvolutioner<real_number> avg_convolutioner;  /*!< Convolutioner computing averaged density     */

  // PARTICLES-CELL MAP BUFFERS
  /*!
   *  DensityKernel class defines maps to locate particles. They are obtained by
   *  computing cumulated density values (putting cells in lexico-graphic order).
   */
  std::vector<int> idx_map, cum_num, raw_num;
  void compute_ind_map_part(void);

public:

  // Init
  DensityKernel(DSMC*);
  ~DensityKernel() = default;

  // Each step of density kernel
  void binning (void);
  void fill_dummy_field (void);
  void compute_reduced_density (void);
  void compute_avg_density (void);

  // Density kernel in a packet
  void perform_density_kernel(void);

  // GETTERS
  inline int get_n_cutoff_x(void) const { return n_cutoff_x; }
  inline int get_n_cutoff_y(void) const { return n_cutoff_y; }
  inline ev_matrix::HaloMaskMatrix<real_number>& get_num_dens_cell(void) { return num_dens_cell; }
  inline const ev_matrix::HaloMaskMatrix<real_number>& get_num_dens_cell(void) const { return num_dens_cell; }

  // void compute_density (void);

  inline const int get_npc(int i, int j) const { return n_part_cell(i,j); }
  inline const real_number get_aveta(int i, int j) const { return average_reduced_density(i,j); }
  inline const real_number get_numdens(int i, int j) const { return num_dens_cell(i,j); }
  inline const ev_matrix::MaskMatrix<int>& get_npc(void) const { return n_part_cell; }
  inline const ev_matrix::MaskMatrix<real_number>& get_aveta(void) const { return average_reduced_density; }
  inline const int iof(int k) const { return cum_num[k]; }
  inline const int ind(int k) const { return idx_map[k]; }

  // DEBUG
  // # # # # #
  void print_binned_particles(void) const;
  void print_reduced_numdens(void) const;
  void print_reduced_aveta(void) const;
  // # # # # #

};

#endif /* EV_DENSITIES_HPP */
