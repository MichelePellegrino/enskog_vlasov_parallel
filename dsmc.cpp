#include "dsmc.hpp"
#include "parallel_environment.hpp"
#include "display.hpp"
#include "configuration.hpp"
#include "boundary.hpp"
#include "grid.hpp"
#include "topology.hpp"
#include "particles.hpp"
#include "density.hpp"
#include "force_field.hpp"
#include "advection.hpp"
#include "collisions.hpp"

DSMC::DSMC(const DefaultString& file_name):
par_env (
  new ParallelEnvironment(this)
),
io_hand (
  new IOHandler(this)
),
conf (
  new ConfigurationReader(this, file_name)
),
rng (
  new RandomEngine(par_env->seed_rank(conf->get_seed()))
),
species (
  new Species(conf->get_diam_fluid(), conf->get_diam_solid(),
    conf->get_mass_fluid(), conf->get_mass_solid())
),
times (
  new Times(conf->get_nc_out(), conf->get_restart(), conf->get_t_ini(),
    conf->get_t_max(), conf->get_t_im(), conf->get_delta_t())
),
boundary (
  new Boundary(this)
),
grid (
  new Grid(this)
),
topology (
  new Topology(this)
),
ensemble (
  new Ensemble(this)
),
density (
  new DensityKernel(this)
),
potential (
  new SutherlandMie(conf->get_phi11(), (conf->get_diam_fluid()/2.0),
    conf->get_gamma11())
),
mean_field (
  new ForceField(this)
),
time_marching (
  new TimeMarching<TM>(this)
),
collision_handler (
  new CollisionHandler(this)
),
correlation (),
stopwatch(4, "millisecond")
{
  // TESTING
  par_env->barrier();
  test_density();
  test_force_field();
  test_time_marching();
  // if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: collision in quarter-subdomain map ###" << std::endl;
  // collision_handler->test_map_to_domain_collision();
  test_collisions();
}


// Testing
void
DSMC::test_density
(void)
{
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: performing density kernel ###" << std::endl;
  stopwatch.local_start(DENSITY_TAG);
  density->perform_density_kernel();
  par_env->barrier();
  stopwatch.local_stop(DENSITY_TAG);
  if (par_env->get_rank()==MPI_MASTER) stopwatch.show_local_elapsed(DENSITY_TAG);
}

void
DSMC::test_force_field
(void)
{
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: computing force field ###" << std::endl;
  stopwatch.local_start(FORCES_TAG);
  mean_field->compute_force_field();
  par_env->barrier();
  stopwatch.local_stop(FORCES_TAG);
  if (par_env->get_rank()==MPI_MASTER) stopwatch.show_local_elapsed(FORCES_TAG);
}

void
DSMC::test_time_marching
(void)
{
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: time-marching ###" << std::endl;
  stopwatch.local_start(ADVECT_TAG);
  time_marching->update_ensemble();
  par_env->barrier();
  stopwatch.local_stop(ADVECT_TAG);
  if (par_env->get_rank()==MPI_MASTER) stopwatch.show_local_elapsed(ADVECT_TAG);
}

void
DSMC::test_collisions
(void)
{
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: performing collisions ###" << std::endl;
  stopwatch.local_start(COLLISION_TAG);
  collision_handler->compute_majorants();
  collision_handler->perform_collisions();
  par_env->barrier();
  stopwatch.local_stop(COLLISION_TAG);
  if (par_env->get_rank()==MPI_MASTER) stopwatch.show_local_elapsed(COLLISION_TAG);
}

// Initialize
void
DSMC::initialize_simulation
(void)
{

}

void
DSMC::dsmc_iteration
(void)
{

}

void
DSMC::dsmc_loop
(void)
{

}
