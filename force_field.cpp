#include "force_field.hpp"
#include "configuration.hpp"
#include "density.hpp"
#include "grid.hpp"

ForceField::ForceField(DSMC* dsmc):

  Motherbase(dsmc),
  n_cutoff_x( density->get_n_cutoff_x() ),
  n_cutoff_y( density->get_n_cutoff_y() ),
  diamol( species->get_diam_fluid() ),
  dx( grid->get_dx() ),
  dy( grid->get_dy() ),
  kernel_function( potential->get_pot_kernel() ),
  finite_integrator( psi, DUMMY_A, DUMMY_B ),
  infinite_integrator( psi, DUMMY_A, DUMMY_B ),

  kernel_matrix(n_cutoff_x, n_cutoff_y, 0.0),
  kernel_matrix_x(n_cutoff_x, n_cutoff_y, 0.0),
  kernel_matrix_y(n_cutoff_x, n_cutoff_y, 0.0),

  force_x_matrix(topology->get_idx_lx(par_env->get_rank()), topology->get_idx_ux(par_env->get_rank()),
    topology->get_idx_ly(par_env->get_rank()), topology->get_idx_uy(par_env->get_rank()), 0.0),
  force_y_matrix(topology->get_idx_lx(par_env->get_rank()), topology->get_idx_ux(par_env->get_rank()),
    topology->get_idx_ly(par_env->get_rank()), topology->get_idx_uy(par_env->get_rank()), 0.0),

  force_x_convolutioner(force_x_matrix, kernel_matrix_x,
    density->get_num_dens_cell(), 0.0),
  force_y_convolutioner(force_y_matrix, kernel_matrix_y,
    density->get_num_dens_cell(), 0.0)

{

  if(par_env->get_rank()==MPI_MASTER)
    std::cout << "### COMPUTING POTENTIAL KERNEL MATRIX ###" << std::endl;
  compute_kernel_matrix();

}

void
ForceField::compute_kernel_matrix
(void)
{
  if (par_env->get_rank() == MPI_MASTER)
  {
    real_number pot_int(0.0), distx(0.0), disty(0.0);
    for (int i = -n_cutoff_x; i<=n_cutoff_x; ++i)
    {
      for (int j = -n_cutoff_y; j<=n_cutoff_y; ++j)
      {
        distx = i*dx;
        disty = j*dy;
        dist2 = distx*distx + disty*disty;
        pot_int = compute_integral();
        kernel_matrix(i,j) = pot_int;
        kernel_matrix_x(i,j) = pot_int*distx*dx*dy;
        kernel_matrix_y(i,j) = pot_int*disty*dx*dy;
      }
    }
  }
  par_env->broadcast(*kernel_matrix.data(), kernel_matrix.size());
  par_env->broadcast(*kernel_matrix_x.data(), kernel_matrix.size());
  par_env->broadcast(*kernel_matrix_y.data(), kernel_matrix.size());
}

real_number
ForceField::compute_integral
(void)
{
  real_number tmp = dist2 - diamol*diamol;
  if ( tmp > 0.0 )
  {
    return (
      infinite_integrator.integrate(psi, -CUTOFF_Z, -sqrt(tmp))
      + finite_integrator.integrate(psi, -sqrt(tmp), sqrt(tmp))
      + infinite_integrator.integrate(psi, sqrt(tmp), CUTOFF_Z)   );
  }
  else if ( tmp < 0.0 )
  {
    return (
      infinite_integrator.integrate(psi, -CUTOFF_Z, -sqrt(-tmp))
      + infinite_integrator.integrate(psi, sqrt(-tmp), CUTOFF_Z)  );
  }
  else
  {
    return (
      infinite_integrator.integrate(psi, -CUTOFF_Z, -ZERO_THRESHOLD)
      + infinite_integrator.integrate(psi, ZERO_THRESHOLD, CUTOFF_Z) );
  }
}

void
ForceField::compute_force_field
(void)
{
  force_x_convolutioner.convolute();
  // force_x_matrix /= (dx*dy);
  force_y_convolutioner.convolute();
  // force_y_matrix /= (dx*dy);
}
