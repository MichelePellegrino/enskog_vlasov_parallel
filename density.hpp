#ifndef EV_DENSITIES_HPP
#define EV_DENSITIES_HPP

#include "matrix.hpp"
#include "types.hpp"
#include "utility.hpp"
#include "motherbase.hpp"

#include <cmath>

// Contents:
//  1) references to: 1.1) grid, 1.2) ensemble
//  2) no. particles per cell (and function for binning)
//  3) number density (and num_dens_dummy ???)
//  4) eta and aveta...

/*! \class DensityCalculator
 *  \brief Class for density and reduced density computation
 */
class DensityKernel : protected Motherbase
{

private:

  // Utilities (maybe reference would be better?)
  real_number reduce_factor;
  int ns_x, ns_y;

  ev_matrix::SlideMaskMatrix<real_number> stencil_x;
  ev_matrix::SlideMaskMatrix<real_number> stencil_y;
  ev_matrix::SlideMaskMatrix<real_number> weights;

  real_number x_min, y_min;
  int n_cutoff_x, n_cutoff_y;

  // Problem: casting Eigen datatype -> solved!

  ev_matrix::MaskMatrix<int> n_part_cell;
  ev_matrix::HaloMaskMatrix<real_number> num_dens_cell;
  ev_matrix::HaloMaskMatrix<real_number> reduced_density;
  ev_matrix::MaskMatrix<real_number> average_reduced_density;

  ev_matrix::MatrixConvolutioner<real_number> avg_convolutioner;

  // Particle-cell map buffers
  std::vector<int> idx_map, cum_num, raw_num;
  void compute_ind_map_part(void);

public:

  // Init
  DensityKernel(DSMC*);
  ~DensityKernel() = default;

  void binning (void);
  void fill_dummy_field (void);
  void compute_reduced_density (void);
  void compute_avg_density (void);

  void perform_density_kernel(void);

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

};

#endif /* EV_DENSITIES_HPP */
