#include "sampling.hpp"
#include "parallel_environment.hpp"
#include "grid.hpp"
#include "topology.hpp"
#include "particles.hpp"
#include "density.hpp"

Sampler::Sampler(DSMC* dsmc):

  Motherbase(dsmc),

  rank(par_env->get_rank()),
  size(par_env->get_size()),

  inner_counter      ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0 ),
  inner_counter_cast ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),

  vx_avg      ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  vy_avg      ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  vz_avg      ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  temp_avg    ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  pxx_avg     ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  pyy_avg     ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  pzz_avg     ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  pxy_avg     ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  pxz_avg     ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  pyz_avg     ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  qx_avg      ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  qy_avg      ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  qz_avg      ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),
  numdens_avg ( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0 ),

  global_vx_avg       ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_vy_avg       ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_vz_avg       ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_temp_avg     ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_pxx_avg      ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_pyy_avg      ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_pzz_avg      ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_pxy_avg      ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_pxz_avg      ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_pyz_avg      ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_qx_avg       ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_qy_avg       ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_qz_avg       ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),
  global_numdens_avg  ( 0, grid->get_n_cells_x(), 0, grid->get_n_cells_y(), 0.0 ),

  n_cells_x (topology->get_idx_ux_rank() - topology->get_idx_lx_rank()),
  low_x (topology->get_idx_lx_rank()),
  low_y (topology->get_idx_ly_rank())

  { }

void
Sampler::reset
(void)
{
  inner_counter = 0;
  inner_counter_cast = 0.0;
  vx_avg = 0.0;
  vy_avg = 0.0;
  vz_avg = 0.0;
  temp_avg = 0.0;
  pxx_avg = 0.0;
  pyy_avg = 0.0;
  pzz_avg = 0.0;
  pxy_avg = 0.0;
  pxz_avg = 0.0;
  pyz_avg = 0.0;
  qx_avg = 0.0;
  qy_avg = 0.0;
  qz_avg = 0.0;
}

void
Sampler::sample
(void)
{
  int nc = topology->get_n_cells();
  int i, j, i_map, j_map, idx_p, idx_c_map;
  real_number vx, vy, vz, e_kin;
  outer_counter++;
  // # # # # #
  // DEBUG !!!
  // # # # # #
  for ( int idx_c = 0; idx_c<nc; ++idx_c )
  {
    i = map_to_domain( lexico_inv_sub(idx_c).first, lexico_inv_sub(idx_c).second ).first;
    j = map_to_domain( lexico_inv_sub(idx_c).first, lexico_inv_sub(idx_c).second ).second;
    idx_c_map = grid->lexico(i,j);
    for ( int k = density->iof(idx_c_map); k < density->iof(idx_c_map+1); ++k )
    {
      idx_p = density->ind(k);
      vx = ensemble->get_vx(idx_p);
      vy = ensemble->get_vy(idx_p);
      vz = ensemble->get_vz(idx_p);
      inner_counter(i,j)++;
      vx_avg(i,j) += vx;
      vy_avg(i,j) += vy;
      vz_avg(i,j) += vz;
      pxx_avg(i,j) += vx*vx;
      pyy_avg(i,j) += vy*vy;
      pzz_avg(i,j) += vz*vz;
      pxy_avg(i,j) += vx*vy;
      pxz_avg(i,j) += vx*vz;
      pyz_avg(i,j) += vy*vz;
      e_kin = vx*vx + vy*vy + vz*vz;
      temp_avg(i,j) += e_kin;
      qz_avg(i,j) += vz*e_kin;
      qx_avg(i,j) += vx*e_kin;
      qy_avg(i,j) += vy*e_kin;
    }
  }
  // # # # # #
}

void
Sampler::average
(void)
{

  inner_counter_cast.copy_cast<int>(inner_counter);
  vx_avg /= inner_counter_cast;
  vy_avg /= inner_counter_cast;
  vz_avg /= inner_counter_cast;
  pxx_avg /= inner_counter_cast;  pxx_avg -= vx_avg*vx_avg;
  pyy_avg /= inner_counter_cast;  pyy_avg -= vy_avg*vy_avg;
  pzz_avg /= inner_counter_cast;  pzz_avg -= vz_avg*vz_avg;
  pxy_avg /= inner_counter_cast;  pxy_avg -= vx_avg*vy_avg;
  pxz_avg /= inner_counter_cast;  pxz_avg -= vx_avg*vz_avg;
  pyz_avg /= inner_counter_cast;  pyz_avg -= vy_avg*vz_avg;
  qz_avg /= 2.0*inner_counter_cast;
  qx_avg /= 2.0*inner_counter_cast;
  qy_avg /= 2.0*inner_counter_cast;
  temp_avg /= inner_counter_cast;
  temp_avg = ( temp_avg -
    vx_avg*vx_avg - vy_avg*vy_avg - vz_avg*vz_avg ) / 3.0;
  numdens_avg = inner_counter_cast / ( (double)outer_counter * grid->get_cell_volume() );

  outer_counter = 0;

}

void
Sampler::gather_samples
(void)
{

  int lx_send = topology->get_idx_lx(rank), ux_send = topology->get_idx_ux(rank);
  int ly_send = topology->get_idx_ly(rank), uy_send = topology->get_idx_uy(rank);

  int lx_recv, ux_recv, ly_recv, uy_recv;

  if (rank != MPI_MASTER)
  {
    // SEND
    numdens_avg.send_block(lx_send, ly_send, ux_send-lx_send, uy_send-ly_send, MPI_MASTER);
    // temp_avg.send_block(lx_send, ly_send, ux_send-lx_send, uy_send-ly_send, MPI_MASTER);
  }
  else
  {
    for (int r = 0; r<size; ++r)
    {
      // RECV
      if ( r==MPI_MASTER )
      {
        // Copy its own data
        global_numdens_avg.set_block(lx_send, ux_send, ly_send, uy_send, numdens_avg);
        // global_temp_avg.set_block(lx_send, ux_send, ly_send, uy_send, temp_avg);
      }
      else
      {
        // Receive data
        lx_recv = topology->get_idx_lx(r);
        ux_recv = topology->get_idx_ux(r);
        ly_recv = topology->get_idx_ly(r);
        uy_recv = topology->get_idx_uy(r);
        global_numdens_avg.recv_block(lx_recv, ly_recv, ux_recv-lx_recv, uy_recv-ly_recv, r);
        // global_temp_avg.recv_block(lx_recv, ly_recv, ux_recv-lx_recv, uy_recv-ly_recv, r);
      }
    }
  }

}
