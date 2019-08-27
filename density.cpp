#include "density.hpp"
#include "grid.hpp"
#include "species.hpp"
#include "particles.hpp"
#include "configuration.hpp"
#include "parallel_environment.hpp"
#include "topology.hpp"

// DEBUG
// # # # # #
#include <iostream>
#include <fstream>
// # # # # #

DensityKernel::DensityKernel
(DSMC* dsmc):

  Motherbase(dsmc),

  reduce_factor( ( ev_const::pi/6.0 ) * ev_utility::power<3>(species->get_diam_fluid()) ),

  ns_x( (int)( species->get_hdiam_fluid() / ( grid->get_dx() * ev_const::sqrt2) ) ),
  ns_y( (int)( species->get_hdiam_fluid() / ( grid->get_dy() * ev_const::sqrt2) ) ),

  stencil_x(ns_x, ns_y, 0.0),
  stencil_y(ns_x, ns_y, 0.0),
  weights(ns_x, ns_y, 0.0),

  x_min(grid->get_x_min()),
  y_min(grid->get_y_min()),

  n_cutoff_x( (int)(conf->get_x_extra()/grid->get_dx())+1 ),
  n_cutoff_y( (int)(conf->get_x_extra()/grid->get_dx())+1 ),

  n_part_cell( topology->get_idx_lx(par_env->get_rank()), topology->get_idx_ux(par_env->get_rank()),
    topology->get_idx_ly(par_env->get_rank()), topology->get_idx_uy(par_env->get_rank()), 0 ),

  num_dens_cell( ev_matrix::MaskMatrix<real_number> (
    -n_cutoff_x+topology->get_idx_lx(par_env->get_rank()),
    topology->get_idx_ux(par_env->get_rank())+n_cutoff_x,
    -n_cutoff_y+topology->get_idx_ly(par_env->get_rank()),
    topology->get_idx_uy(par_env->get_rank())+n_cutoff_y ),
    topology->get_topology_map(),
    n_cutoff_x, n_cutoff_y ),

  reduced_density( ev_matrix::MaskMatrix<real_number> (
    -n_cutoff_x+topology->get_idx_lx(par_env->get_rank()),
    topology->get_idx_ux(par_env->get_rank())+n_cutoff_x,
    -n_cutoff_y+topology->get_idx_ly(par_env->get_rank()),
    topology->get_idx_uy(par_env->get_rank())+n_cutoff_y ),
    topology->get_topology_map(),
    n_cutoff_x, n_cutoff_y ),

  average_reduced_density( topology->get_idx_lx(par_env->get_rank()), topology->get_idx_ux(par_env->get_rank()),
    topology->get_idx_ly(par_env->get_rank()), topology->get_idx_uy(par_env->get_rank()), 0.0 ),

  avg_convolutioner( average_reduced_density, weights, reduced_density, 0.0 ),

  idx_map( ensemble->get_n_particles(), 0 ),
  cum_num( grid->get_n_cells(), 0 ),
  raw_num( grid->get_n_cells(), 0 )

  {
    // Initialize weight
    if (par_env->get_rank()==MPI_MASTER)
    {
      real_number dx = grid->get_dx();
      real_number dy = grid->get_dy();
      real_number sigma = species->get_diam_fluid();
      real_number hsigma = sigma/2.0;
      real_number sum_w = 0.0;
      std::cout << "### COMPUTING AVERAGING WEIGHTS ###" << std::endl;
      for (int i =-ns_x; i<=ns_x; ++i)
      {
        for (int j =-ns_y; j<=ns_y; ++j)
        {
          stencil_x(i, j) = i * dx;
          stencil_y(i, j) = j * dy;
          weights(i,j) = 12.0 / ( ev_const::pi * ev_utility::power<3>(sigma) )
            * sqrt( hsigma * hsigma - stencil_x(i,j) * stencil_x(i,j)
            - stencil_y(i,j) * stencil_y(i,j) ) * dx * dy;
          sum_w = sum_w + weights(i,j);
        }
      }
      weights /= sum_w;
    }
    par_env->broadcast(*weights.data(), weights.size());
    compute_ind_map_part();
  }

void
DensityKernel::binning
(void)
{
  // Binning and setting each particle index
  n_part_cell = 0;
  int NP = ensemble->get_n_particles();
  int i,j;
  for (int k = 0; k<NP; ++k)
  {
    i = ensemble->get_cx(k);
    j = ensemble->get_cy(k);
    n_part_cell(i,j) += 1;
  }
  // Setting cumulate density
  int NC = grid->get_n_cells();
  cum_num.assign(NC, 0);
  // There is surely a more elegant way to proceed, but I don't have time!
  // (the 'if' condition is so inefficient...)
  for (int k = 1; k<NC; ++k)
  {
    i = grid->lexico_inv(k-1).first;
    j = grid->lexico_inv(k-1).second;
    if ( i >= n_part_cell.get_lx() && i < n_part_cell.get_ux()
      && j >= n_part_cell.get_ly() && j < n_part_cell.get_uy() )
        cum_num[k] = cum_num[k-1] + n_part_cell(i, j);
  }
}

void
DensityKernel::compute_ind_map_part
(void)
{
  int NP = ensemble->get_n_particles();
  int jc, k;
  idx_map.assign(NP, 0);
  raw_num.assign(grid->get_n_cells(), 0);
  for (int jp = 0; jp<NP; ++jp)
  {
    jc = grid->lexico(ensemble->get_cx(jp), ensemble->get_cy(jp));
    raw_num[jc] += 1;
    k = cum_num[jc] + raw_num[jc];
    idx_map[k] = jp;
  }
}

void
DensityKernel::fill_dummy_field
(void)
{

  // Copying inner data
  num_dens_cell.copy_patch<int>(n_part_cell, n_part_cell.get_lx(), n_part_cell.get_ly());
  num_dens_cell /= grid->get_cell_volume();

  // CASE: periodic B.C. on each and every edge...
  for(int r = 0; r < N_BUF; ++r )
  {
    if ( num_dens_cell.rank_from_idx(r) == par_env->get_rank() )
      num_dens_cell.set_outer_block( r, num_dens_cell.get_inner_block( ev_matrix::reflect_idx(r) ) );
  }

  int idx_recv = MPI_NO_RANK;
  int rank_send = MPI_NO_RANK;
  int n_recv = 0;
  int rank = par_env->get_rank(), size = par_env->get_size();
  par_env->barrier();
  for ( int k = 0; k<size; ++k )
  {
    if (rank == k)
    {
      for (int idx = 0; idx<8; ++idx)
      {
        rank_send = num_dens_cell.rank_from_idx(idx);
        if (rank_send != MPI_NO_RANK && rank_send!=rank )
        {
          par_env->send(idx, rank_send);
          num_dens_cell.send_block(idx);
        }
      }
    }
    else
    {
      n_recv = num_dens_cell.n_idx_from_rank(k);
      // std::cout << "rank " << rank << " should receive " << num_dens_cell.n_idx_from_rank(k) << " times" << std::endl;
      for (int i = 0; i<n_recv; ++i)
      {
        par_env->receive(idx_recv, k);
        idx_recv = num_dens_cell.reflect_idx(idx_recv);
        num_dens_cell.recv_block(idx_recv);
      }
    }
    par_env->barrier();
  }

}

void
DensityKernel::compute_reduced_density
(void)
{
  reduced_density = num_dens_cell;
  reduced_density *= reduce_factor;
}

void
DensityKernel::compute_avg_density
(void)
{
  avg_convolutioner.convolute();
}

void
DensityKernel::perform_density_kernel
(void)
{
  binning();
  fill_dummy_field();
  compute_reduced_density();
  compute_avg_density();
  // DEBUG
  // # # # # #
  // print_binned_particles();
  // print_reduced_numdens();
  // print_reduced_aveta();
  // # # # # #
}

// TESTING

void
DensityKernel::print_binned_particles
(void) const
{
  for (int r = 0; r<par_env->get_size(); ++r) {
    if ( r==par_env->get_rank() )
      std::cout << "rank " << r << " : binned particles = " << n_part_cell.sum() << std::endl;
  }
  par_env->barrier();
}

void
DensityKernel::print_reduced_numdens
(void) const
{
  for (int r = 0; r<par_env->get_size(); ++r) {
    if ( r==par_env->get_rank() )
      std::cout << "rank " << r << " : reduced numdens = " << num_dens_cell.sum() << std::endl;
  }
  par_env->barrier();
}

void
DensityKernel::print_reduced_aveta
(void) const
{
  for (int r = 0; r<par_env->get_size(); ++r) {
    if ( r==par_env->get_rank() )
      std::cout << "rank " << r << " : reduced aveta = " << average_reduced_density.sum() << std::endl;
  }
  par_env->barrier();
}
