#include "particles.hpp"
#include "parallel_environment.hpp"
#include "grid.hpp"
#include "boundary.hpp"
#include "configuration.hpp"
#include "topology.hpp"

Ensemble::Ensemble
(DSMC* dsmc):
  Motherbase(dsmc),
  n_particles(conf->get_n_part()),
  T_ini(conf->get_T_ini()),
  mass(species->get_mass_fluid()),
  domain_rank(par_env->get_rank()),
  n_subdom(par_env->get_size()),
  particles(n_particles)
  {
    commit_particle_type();
    /* Init. of the incidence matrix (for communication) */
    inc_matrix.resize(par_env->get_size2(), 0);
    /* Populate the ensemble */
    if(par_env->is_root())
      std::cout << "### POPULATING ENSEMBLE ###" << std::endl;
    populate();
  }

void
Ensemble::populate
(void)
{

  // Basic type
  for ( int i = 0; i<n_particles; ++i )
  {

    particles[i].xp = grid->get_x_min() + rng->sample_uniform() * ( grid->get_x_max() - grid->get_x_min() );
    particles[i].yp = grid->get_y_min() + rng->sample_uniform() * ( grid->get_y_max() - grid->get_y_min() );
    particles[i].cell_x = (int) ( (particles[i].xp-grid->get_x_min() ) / grid->get_dx() );
    particles[i].cell_y = (int) ( (particles[i].yp-grid->get_y_min() ) / grid->get_dy() );
    rng->sample_box_muller (
      mass, 0.0, 0.0, T_ini,
      particles[i].vx,
      particles[i].vy,
      particles[i].vz );
    particles[i].p_tag = i;
    particles[i].r_tag = topology->tag_subdom(particles[i].cell_x, particles[i].cell_y);

    add_to_inc_matrix(i);

  }

  par_env->barrier();
  exchange_particles();

}

void
Ensemble::exchange_particles
(void)
{
  
  fill_inc_matrix();

    for ( int i = 0; i<n_subdom; ++i )
    {
      // IF 'i' IS THE CURRENT PROCESS, SEND
      if ( i == domain_rank )
      {
        for ( int j = 0; j<n_subdom; ++j )
        {
          // Avoid 1) sending to itself; 2) sending 0 particles
          if ( i!=j && inc_matrix[i*n_subdom+j]>0 )
          {
            /* USE A SEND VERSION FROM PAR. ENV. !!! -> MPI_PARTICLE_TYPE should be in parallel_environment.hpp */
            MPI_Send ( temp_send_map[j].data(), temp_send_map[j].size(), MPI_PARTICLE_TYPE, j, 0, MPI_COMM_WORLD );
          }
        }
      }
      // If 'i' is not the curr. process, receive
      else
      {
        if ( inc_matrix[i*n_subdom+domain_rank]>0 )
        {
          /* USE A RECV VERSION FROM PAR. ENV. !!! -> MPI_PARTICLE_TYPE should be in parallel_environment.hpp */
          MPI_Recv ( &(temp_send_map[domain_rank][idx_r]), inc_matrix[i*n_subdom+domain_rank],
            MPI_PARTICLE_TYPE, i, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
          idx_r += inc_matrix[i*n_subdom+domain_rank];
        }
      }
    }
    // LET'S SWAP TO BE MORE EFFICIENT
    n_particles = n_particles + inc_particles - dep_particles;
    std::swap(particles, temp_send_map[domain_rank]);

}

void
Ensemble::fill_inc_matrix
(void)
{

  // Broadcast incidence matrix (so now it's filled)
  for ( int i = 0; i<n_subdom; ++i )
    par_env->broadcast(inc_matrix[i*n_subdom], n_subdom, i);
  // Sync. processes before send/recv routines
  par_env->barrier();
  // Make a balance (in/out)
  for ( int i = domain_rank*n_subdom; i<(domain_rank+1)*n_subdom; ++i )
    dep_particles += inc_matrix[i];
  for ( int i = domain_rank; i<par_env->get_size2(); i+=n_subdom )
    inc_particles += inc_matrix[i];
  // DEBUG
  // # # # # #
  assert( n_particles+inc_particles-dep_particles>=0 && "Can't balance in/out particles" );
  assert( dep_particles>=0 && "Negative outflow of departing particles" );
  assert( inc_particles>=0 && "Negative inflow of departing particles" );
  // # # # # #
  temp_send_map[domain_rank].resize(n_particles+inc_particles-dep_particles);
  idx_r = n_particles - dep_particles;

}

/*! \fn
 *  \brief [...]
 */
void
Ensemble::add_to_inc_matrix
(int i)
{

  temp_send_map[particles[i].r_tag].push_back(particles[i]);
  if ( particles[i].r_tag != domain_rank )
    inc_matrix[ domain_rank*n_subdom + particles[i].r_tag ]++;

}

/*! \fn clear_buffers
 *  \brief Clears send/receive buffers
 */
void
Ensemble::clear_buffers
(void)
{

  inc_particles = 0;
  dep_particles = 0;
  inc_matrix.assign(par_env->get_size2(), 0);

  temp_send_map.clear();

}

void
Ensemble::commit_particle_type
(void)
{
  /* Definition of custom data type */
  MPI_Datatype oldtypes[2];
  int blockcounts[2];
  MPI_Aint offsets[2], lb, extent;
  offsets[0] = 0;
  oldtypes[0] = MPI_DOUBLE;
  blockcounts[0] = 5;         // 5 doubles: xp, yp, vx, vy, vz
  MPI_Type_get_extent(MPI_DOUBLE, &lb, &extent);
  offsets[1] = 5 * extent;
  oldtypes[1] = MPI_INT;
  blockcounts[1] = 4;         // 4 integers: cell_x, cell_y, p_tag, r_tag
  MPI_Type_create_struct(2, blockcounts, offsets, oldtypes, &MPI_PARTICLE_TYPE);
  MPI_Type_commit(&MPI_PARTICLE_TYPE);
}

/*
void
Ensemble::test_stream
(real_number dt)
{
  // PERIODIC B.C. ON EVERY EDGE
  for ( int i = 0; i<n_particles; ++i )
  {
    particles[i].xp = particles[i].xp + particles[i].vx*dt;
    particles[i].yp = particles[i].yp + particles[i].vy*dt;
    if ( particles[i].xp > grid->get_x_max() )
      particles[i].xp = grid->get_x_min() + particles[i].xp - grid->get_x_max();
    else if ( particles[i].xp < grid->get_x_min() )
      particles[i].xp = grid->get_x_max() - grid->get_x_min() + particles[i].xp;
    if ( particles[i].yp > grid->get_y_max() )
      particles[i].yp = grid->get_y_min() + particles[i].yp - grid->get_y_max();
    else if ( particles[i].yp < grid->get_y_min() )
      particles[i].yp = grid->get_y_max() - grid->get_y_min() + particles[i].yp;
  }
}
*/
