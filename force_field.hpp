#ifndef EV_FORCE_FIELD_HPP
#define EV_FORCE_FIELD_HPP

#include "motherbase.hpp"
#include "matrix.hpp"
#include "integration.hpp"
#include "parallel_environment.hpp"
#include "topology.hpp"

#include <functional>

using namespace ev_numeric;

#define DUMMY_A -1
#define DUMMY_B 1

class ForceField : protected Motherbase
{

private:

  int n_cutoff_x, n_cutoff_y;

  real_number diamol;
  real_number dx, dy;

  const real_number zero_threshold = 1.0e-4;
  real_number dist2 = 0;

  std::function<real_number(real_number)> kernel_function;
  std::function<real_number(real_number)> psi = [this](const real_number& z)
    -> real_number
  {
    return kernel_function( sqrt( dist2 + z*z ) );
  };

  NumericalIntegrator<Finite> finite_integrator;
  NumericalIntegrator<Infinite> infinite_integrator;

  ev_matrix::SlideMaskMatrix<real_number> kernel_matrix;

  void compute_kernel_matrix (void);
  real_number compute_integral (void);

  ev_matrix::MaskMatrix<real_number> force_x_matrix;
  ev_matrix::MaskMatrix<real_number> force_y_matrix;

  ev_matrix::MatrixConvolutioner<real_number> force_x_convolutioner;
  ev_matrix::MatrixConvolutioner<real_number> force_y_convolutioner;

public:

  ForceField(DSMC*);
  ~ForceField() = default;

  void compute_force_field(void);

  inline int get_n_cutoff_x(void) const { return n_cutoff_x; }
  inline int get_n_cutoff_y(void) const { return n_cutoff_y; }

  inline real_number get_force_x(int i, int j) { return force_x_matrix(i,j); }
  inline real_number get_force_y(int i, int j) { return force_y_matrix(i,j); }

  inline const ev_matrix::SlideMaskMatrix<real_number>& get_kernel_matrix(void) const { return kernel_matrix; }
  inline const ev_matrix::MaskMatrix<real_number>& get_force_x(void) { return force_x_matrix; }
  inline const ev_matrix::MaskMatrix<real_number>& get_force_y(void) { return force_y_matrix; }

};

#endif /* EV_FORCE_FIELD_HPP */
