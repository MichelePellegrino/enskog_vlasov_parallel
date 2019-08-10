#include "topology.hpp"
#include "parallel_environment.hpp"
#include "grid.hpp"

#include <cmath>

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
  idx_lx_rank(0),
  idx_ly_rank(0),
  idx_ux_rank(0),
  idx_uy_rank(0),
  topology_map( -n_cutoff_x, grid->get_n_cells_x()+n_cutoff_x,
    -n_cutoff_y, grid->get_n_cells_y()+n_cutoff_y, MPI_NO_RANK )
{

  fill_topology_map();
  setup_quarters();

  if (rank==MPI_MASTER)
  {
    std::ofstream file1("output_files/topology.txt");
    file1 << topology_map;
    file1.close();
  }

  // DEBUG
  // # # # # #
  /*
  for (int r = 0; r<size; ++r)
  {
    if (rank==r)
    {
      std::cout << " >> Rank " << r << " has limits: ["<<idx_lx[r]<<","<<idx_ux[r]<<";"<<idx_ly[r]<<","<<idx_uy[r]<<"]" << std::endl;
      std::cout << "   quarters: ["<<quarter_lx[0]<<","<<quarter_ux[0]<<";"<<quarter_ly[0]<<","<<quarter_uy[0]<<"]" <<
        " ["<<quarter_lx[1]<<","<<quarter_ux[1]<<";"<<quarter_ly[1]<<","<<quarter_uy[1]<<"]" <<
        " ["<<quarter_lx[2]<<","<<quarter_ux[2]<<";"<<quarter_ly[2]<<","<<quarter_uy[2]<<"]" <<
        " ["<<quarter_lx[3]<<","<<quarter_ux[3]<<";"<<quarter_ly[3]<<","<<quarter_uy[3]<<"]" << std::endl;
    }
    par_env->barrier();
  }
  */
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

  switch (topology_type)
  {
    case StripesX:
      if (rank == MPI_MASTER) std::cout << " >> Topology: stripes (along x dir.)" << std::endl;
      setup_stripes_x();
      break;
    case StripesY:
      if (rank == MPI_MASTER) std::cout << " >> Topology: stripes (along y dir.)" << std::endl;
      setup_stripes_y();
      break;
    case Quadrants:
      if (rank == MPI_MASTER) std::cout << " >> Topology: quadrants" << std::endl;
      setup_quadrants();
      break;
    default:
      if (rank == MPI_MASTER) throw "|!| WARNING: invalid MPI topology";
  }

  par_env->all_gather_inp(idx_lx[0]);
  par_env->all_gather_inp(idx_ly[0]);
  par_env->all_gather_inp(idx_ux[0]);
  par_env->all_gather_inp(idx_uy[0]);

  // PERIODIC B.C.
  for (int r = 0; r<size; ++r)
    topology_map.set_block(idx_lx[r], idx_ux[r], idx_ly[r], idx_uy[r], r);
  int nx = grid->get_n_cells_x(), ny = grid->get_n_cells_y();
  int i_remap, j_remap;
  for ( int i = topology_map.get_lx(); i<topology_map.get_ux(); ++i )
  {
    for ( int j = topology_map.get_ly(); j<topology_map.get_uy(); ++j )
    {
      i_remap = (i+nx)%nx;
      j_remap = (j+ny)%ny;
      topology_map(i,j) = topology_map(i_remap, j_remap);
    }
  }

  idx_lx_rank = idx_lx[par_env->get_rank()];
  idx_ly_rank = idx_ly[par_env->get_rank()];
  idx_ux_rank = idx_ux[par_env->get_rank()];
  idx_uy_rank = idx_uy[par_env->get_rank()];
}

void
Topology::setup_stripes_x
(void)
{
  int offset_x = grid->get_n_cells_x() / size;
  int rem_x = grid->get_n_cells_x()%size;
  idx_ly[rank] = 0;
  idx_lx[rank] = (rank % size) * offset_x + std::min( rem_x, rank % size );
  idx_uy[rank] = grid->get_n_cells_y();
  idx_ux[rank] = idx_lx[rank] + offset_x + (rank % size < rem_x);
}

void
Topology::setup_stripes_y
(void)
{
  int offset_y = grid->get_n_cells_y() / size;
  int rem_y = grid->get_n_cells_y()%size;
  idx_lx[rank] = 0;
  idx_ly[rank] = (rank % size) * offset_y + std::min( rem_y, rank % size );
  idx_ux[rank] = grid->get_n_cells_x();
  idx_uy[rank] = idx_ly[rank] + offset_y + (rank % size < rem_y);
}

void
Topology::setup_quadrants
(void)
{
  int size_rem = size;
  int n_proc_x, n_proc_y, int_sqrt;
  do {
    int_sqrt = (int)sqrt(size_rem);
    if ( size%int_sqrt == 0 )
    {
      n_proc_x = int_sqrt;
      n_proc_y = size/n_proc_x;
    }
    else
      size_rem--;
  } while (size%int_sqrt != 0);
  int offset_x = grid->get_n_cells_x() / n_proc_x;
  int offset_y = grid->get_n_cells_y() / n_proc_y;
  int rem_x = grid->get_n_cells_x()%n_proc_x;
  int rem_y = grid->get_n_cells_y()%n_proc_y;
  idx_lx[rank] = (rank % n_proc_x) * offset_x + std::min( rem_x, rank % n_proc_x );
  idx_ly[rank] = (rank / n_proc_x) * offset_y + std::min( rem_y, rank / n_proc_x );
  idx_ux[rank] = idx_lx[rank] + offset_x + (rank % n_proc_x < rem_x);
  idx_uy[rank] = idx_ly[rank] + offset_y + (rank / n_proc_x < rem_y);
}

void
Topology::setup_quarters
(void)
{
  // assert(N_COLLISION_SUBDOM>=4 && N_COLLISION_SUBDOM%2==0 && "Not yet implemented for a generic number of collisional subdomains");
  assert(N_COLLISION_SUBDOM==4 && "Not yet implemented for a no. subdomains != 4");
  int div_x = sqrt(N_COLLISION_SUBDOM);
  int div_y = N_COLLISION_SUBDOM/div_x;
  int delta_x = (idx_ux_rank-idx_lx_rank)/div_x;
  int delta_y = (idx_uy_rank-idx_ly_rank)/div_y;
  int reminder_x = ((idx_ux_rank-idx_lx_rank)%div_x>0);
  int reminder_y = ((idx_uy_rank-idx_ly_rank)%div_y>0);
  // *** NB: ACCOUNT FOR THE REMINDER !!! ***
  for (int i = 0; i<N_COLLISION_SUBDOM; ++i)
  {
    quarter_lx[i] = idx_lx_rank + delta_x*(i%2);
    quarter_ly[i] = idx_ly_rank + delta_y*(i>=2);
    quarter_ux[i] = quarter_lx[i] + delta_x + (rank%2) * reminder_x;
    quarter_uy[i] = quarter_ly[i] + delta_y + (rank>=2) * reminder_y;
    // quarter_ux[i] = idx_lx_rank + delta_x*(i%div_x+1);
    // quarter_uy[i] = idx_ly_rank + delta_y*(i%div_y+1);
    n_cells_quarter[i] = ( quarter_ux[i]-quarter_lx[i] ) * ( quarter_uy[i]-quarter_ly[i] );
  }
}
