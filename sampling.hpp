#ifndef EV_SAMPLING_HPP
#define EV_SAMPLING_HPP

#include "motherbase.hpp"
#include "matrix.hpp"

class Sampler : protected Motherbase
{

private:

  int rank, size;

  int outer_counter = 0;

  ev_matrix::MaskMatrix<int> inner_counter;
  ev_matrix::MaskMatrix<real_number> inner_counter_cast;

  // LOCAL QUANTITIES
  ev_matrix::MaskMatrix<real_number> vx_avg;
  ev_matrix::MaskMatrix<real_number> vy_avg;
  ev_matrix::MaskMatrix<real_number> vz_avg;
  ev_matrix::MaskMatrix<real_number> temp_avg;
  ev_matrix::MaskMatrix<real_number> pxx_avg;
  ev_matrix::MaskMatrix<real_number> pyy_avg;
  ev_matrix::MaskMatrix<real_number> pzz_avg;
  ev_matrix::MaskMatrix<real_number> pxy_avg;
  ev_matrix::MaskMatrix<real_number> pxz_avg;
  ev_matrix::MaskMatrix<real_number> pyz_avg;
  ev_matrix::MaskMatrix<real_number> qx_avg;
  ev_matrix::MaskMatrix<real_number> qy_avg;
  ev_matrix::MaskMatrix<real_number> qz_avg;
  ev_matrix::MaskMatrix<real_number> numdens_avg;
  ev_matrix::MaskMatrix<real_number> fx_avg;
  ev_matrix::MaskMatrix<real_number> fy_avg;

  ev_matrix::MaskMatrix<real_number> ph_avg;

  // GLOBAL QUANTITIES

  ev_matrix::MaskMatrix<real_number> global_temp_avg;
  ev_matrix::MaskMatrix<real_number> global_ph_avg;
  ev_matrix::MaskMatrix<real_number> global_numdens_avg;

  int n_cells_x;
  int low_x, low_y;
  inline std::pair<real_number, real_number> lexico_inv_sub(int idx) const
  {
    int j = idx / n_cells_x;
    int i = idx - j * n_cells_x;
    return std::make_pair(i, j);
  }
  inline std::pair<real_number, real_number> map_to_domain(int i_loc, int j_loc) const
  {
    return std::make_pair(i_loc + low_x, j_loc + low_y);
  }

public:

  Sampler(DSMC*);
  ~Sampler() = default;

  // Sampling
  void reset(void);
  void sample(void);
  void average(void);

  // Communication
  void gather_samples(void);

  // Getters

  inline const ev_matrix::MaskMatrix<real_number>& get_temp_avg(void) const { return global_temp_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_temp_avg(void) { return global_temp_avg; }

  inline const ev_matrix::MaskMatrix<real_number>& get_ph_avg(void) const { return global_ph_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_ph_avg(void) { return global_ph_avg; }

  inline const ev_matrix::MaskMatrix<real_number>& get_numdens_avg(void) const { return global_numdens_avg; }
  inline ev_matrix::MaskMatrix<real_number>& get_numdens_avg(void) { return global_numdens_avg; }

};

#endif /* EV_SAMPLING_HPP */
