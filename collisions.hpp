#ifndef EV_COLLISIONS_HPP
#define EV_COLLISIONS_HPP

#include "motherbase.hpp"
#include "matrix.hpp"

#include <valarray>
#include <cmath>
#include <algorithm>
#include <vector>
#include <array>

#ifndef TEST_COEFF_MULT
#define TEST_COEFF_MULT 5
#endif

#ifndef N_COLLISION_SUBDOM
#define N_COLLISION_SUBDOM 4
#endif

class CollisionHandler : protected Motherbase
{

private:

  std::array<int, N_COLLISION_SUBDOM> subdom_order;

  ev_matrix::MaskMatrix<real_number> a11;
  ev_matrix::MaskMatrix<real_number> vrmax11;

  ev_matrix::MaskMatrix<real_number> anew;

  int n_coll;
  ev_matrix::MaskMatrix<int> n_coll_cell;

  std::valarray<real_number> scaled_k;

  // References
  const int& npart;
  const real_number& xmin, xmax, ymin, ymax;
  const real_number& rdx, rdy;
  const real_number& sigma, delta_t;

  // Utilities
  real_number vr = 0.0;
  const real_number lx_sub, ux_sub, ly_sub, uy_sub;

  inline void gen_scaled_k(void)
  {
    rng->sample_unit_sphere(scaled_k[0], scaled_k[1], scaled_k[2]);
    scaled_k *= sigma;
  }

public:

  CollisionHandler(DSMC*);
  ~CollisionHandler() = default;

  void shuffle_order(void);

  void compute_majorants(void);
  void compute_collision_number(void);

  // DEBUG
  void print_subdomain_order(void) const;

};

#endif /* EV_COLLISIONS_HPP */
