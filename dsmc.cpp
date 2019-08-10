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
correlation ()
{
  // TESTING
  /*
  par_env->barrier();
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: binning and number density ###" << std::endl;
  density->binning();
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: filling density field ###" << std::endl;
  density->fill_dummy_field();
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: computing reduced density ###" << std::endl;
  density->compute_reduced_density();
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: computing averaged density ###" << std::endl;
  density->compute_avg_density();
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: computing force field ###" << std::endl;
  mean_field->compute_force_field();
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: time-marching ###" << std::endl;
  time_marching->update_ensemble();
  // if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: shuffling collisional subdomains ###" << std::endl;
  // int n_tests = 5;
  // for (int i = 0; i<n_tests; ++i)
  // {
  //   collision_handler->shuffle_order();
  //   if (par_env->get_rank()==2)
  //     collision_handler->print_subdomain_order();
  // }
  // par_env->barrier();
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: computing majorants ###" << std::endl;
  collision_handler->compute_majorants();
  par_env->barrier();
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: collision in quarter-subdomain map ###" << std::endl;
  collision_handler->test_map_to_domain_collision();
  if (par_env->get_rank()==MPI_MASTER) std::cout << "### TEST: performing collisions ###" << std::endl;
  collision_handler->perform_collisions();
  */
}
