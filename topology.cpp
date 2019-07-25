#include "topology.hpp"
#include "parallel_environment.hpp"
#include "grid.hpp"

// DEBUG
// # # # # #
#include <fstream>
// # # # # #

Topology::Topology
(DSMC* dsmc):
  Motherbase(dsmc),
  rank(par_env->get_rank()),
  size(par_env->get_size()),
  n_cutoff_x( (int)(conf->get_x_extra()/grid->get_dx())+1 ),
  n_cutoff_y( (int)(conf->get_x_extra()/grid->get_dx())+1 ),
  idx_lx(size,0),
  idx_ly(size,0),
  idx_ux(size,0),
  idx_uy(size,0),
  xmin_sub(size,0.0),
  ymin_sub(size,0.0),
  xmax_sub(size,0.0),
  ymax_sub(size,0.0),
  topology_map( -n_cutoff_x, grid->get_n_cells_x()+n_cutoff_x,
    -n_cutoff_y, grid->get_n_cells_y()+n_cutoff_y, MPI_NO_RANK )
{
  fill_topology_map();
  // DEBUG
  // # # # # #
  /*
  std::cout << "RANK " << rank << " HAS LIMITS : [" << idx_lx[rank] <<
    "," << idx_ux[rank] << ";" << idx_ly[rank] << "," << idx_uy[rank] << "]" << std::endl;
  */
  // # # # # #
  // DEBUG
  // # # # # #
  if (rank==MPI_MASTER)
  {
    std::ofstream file1("output_files/topology.txt");
    file1 << topology_map;
    file1.close();
  }
  // # # # # #
  for (int r = 0; r<size; ++r)
  {
    if (rank==r)
    {
      std::cout << "Rank " << r << " has limits: ["<<idx_lx[r]<<","<<idx_ux[r]<<";"<<idx_ly[r]<<","<<idx_uy[r]<<"]" << std::endl;
      // std::cout << "Rank " << r << " has limits: ["<<xmin_sub[r]<<","<<xmax_sub[r]<<";"<<ymin_sub[r]<<","<<ymax_sub[r]<<"]" << std::endl;
    }
  }
  // # # # # #
}

int
Topology::tag_subdom
(int i, int j)
{
  return topology_map(i,j);
}

void
Topology::fill_topology_map
(void)
{
  // Dummy versions: quarters with 4 cores
  assert(size==4 && "It works with 4, don't try anything else");
  int offset_x = grid->get_n_cells_x() / 2;
  int offset_y = grid->get_n_cells_y() / 2;
  int rem_x = (grid->get_n_cells_x()%2>0);
  int rem_y = (grid->get_n_cells_y()%2>0);
  idx_lx[rank] = (rank%2) * offset_x;
  idx_ly[rank] = (rank>=2) * offset_y;
  idx_ux[rank] = idx_lx[rank] + offset_x + (rank%2) * rem_x;
  idx_uy[rank] = idx_ly[rank] + offset_y + (rank>=2) * rem_y;
  par_env->all_gather_inp(idx_lx[0]);
  par_env->all_gather_inp(idx_ly[0]);
  par_env->all_gather_inp(idx_ux[0]);
  par_env->all_gather_inp(idx_uy[0]);
  /* PERIODIC B.C. */
  for (int r = 0; r<4; ++r)
  {
    topology_map.set_block(idx_lx[r], idx_ux[r], idx_ly[r], idx_uy[r], r);
    topology_map.set_block(idx_lx[r], idx_ux[r],
      (idx_ly[r]==0) * grid->get_n_cells_y() - (idx_ly[r]!=0) * n_cutoff_y,
      (idx_ly[r]==0) * (grid->get_n_cells_y()+n_cutoff_y) - (idx_ly[r]!=0) * 0,
      r);
    topology_map.set_block(
      (idx_lx[r]==0) * grid->get_n_cells_x() - (idx_lx[r]!=0) * n_cutoff_x,
      (idx_lx[r]==0) * (grid->get_n_cells_x()+n_cutoff_x) - (idx_lx[r]!=0) * 0,
      idx_ly[r], idx_uy[r],
      r);
    topology_map.set_block(
      (idx_lx[r]==0) * grid->get_n_cells_x() - (idx_lx[r]!=0) * n_cutoff_x,
      (idx_lx[r]==0) * (grid->get_n_cells_x()+n_cutoff_x) - (idx_lx[r]!=0) * 0,
      (idx_ly[r]==0) * grid->get_n_cells_y() - (idx_ly[r]!=0) * n_cutoff_y,
      (idx_ly[r]==0) * (grid->get_n_cells_y()+n_cutoff_y) - (idx_ly[r]!=0) * 0,
      r);

  }
  /*
  switch (topology_type)
  {
    case StripesX:
      // ...
      break;
    case StripesY:
      // ...
      break;
    case Quarters:
      // ...
      break;
    default:
      throw "Invalid MPI topology";
  }
  */
  fill_limits();
}

void
Topology::fill_limits
(void)
{
  real_number dx = grid->get_dx(), dy =grid->get_dy();
  for (int r = 0; r<size; ++r)
  {
    xmin_sub[r] = dx*idx_lx[r];
    ymin_sub[r] = dy*idx_ly[r];
    xmax_sub[r] = dx*idx_ux[r];
    ymax_sub[r] = dy*idx_uy[r];
  }
}
