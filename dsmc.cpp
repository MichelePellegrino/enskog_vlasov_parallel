/*! \file dsmc.cpp
 *  \brief Source code for DSMC class
 */

#include "dsmc.hpp"
#include "parallel_environment.hpp"
#include "display.hpp"
#include "configuration.hpp"
#include "boundary.hpp"
#include "grid.hpp"
#include "topology.hpp"
#include "particles.hpp"
#include "thermostat.hpp"
#include "density.hpp"
#include "force_field.hpp"
#include "advection.hpp"
#include "collisions.hpp"
#include "sampling.hpp"
#include "output.hpp"

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
thermostat (
  new Thermostat(this)
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
sampler (
  new Sampler(this)
),
output (
  new Output(this)
),
correlation (),
stopwatch (N_TAGS, "millisecond"),
n_iter_thermo ( conf->get_niter_thermo() ),
n_iter_sample ( conf->get_niter_sampling() )
{

  /*!
   *  Establish whether the mean field has to be computed; mean field computation
   *  is very time consuming: better de-activate the routine if not needed
   */
  if ( conf->get_mean_f_gg() == 'y' || conf->get_mean_f_gg() == 'Y' )
  {
    mean_field_gg = true;
    if (par_env->is_root())
      std::cout << "### MEAN-FIELD ON ###" << std::endl;
  }
  else
  {
    mean_field_gg = false;
    if (par_env->is_root())
      std::cout << "### MEAN-FIELD OFF ###" << std::endl;
  }

  initialize_simulation();
  par_env->barrier();

  /*!
   *  Preliminary tests: comment the ones that are not needed
   */
  // test_sampling();
  // test_output();

  /*!
   *  DSMC loop is tested for analysis purposes (e.g. speed-up)
   */
  test_loop(DEFAULT_DUMMY_TEST_ITER);

}


// TESTING FUNCTIONALITIES

/*! \fn void DSMC::test_density (void)
    \brief Tests the density kernel and stores partial time
*/
void
DSMC::test_density
(void)
{
  if ( par_env->is_root() ) std::cout << "### TEST: performing density kernel ###" << std::endl;
  stopwatch.local_start(DENSITY_TAG);
  density->perform_density_kernel();
  par_env->barrier();
  stopwatch.local_stop(DENSITY_TAG);
  stored_elapsed_times[DENSITY_TAG].push_back(stopwatch.get_local_elapsed(DENSITY_TAG));
  if ( par_env->is_root() ) stopwatch.show_local_elapsed(DENSITY_TAG);
}

/*! \fn void DSMC::test_force_field (void)
    \brief Tests mean-field computation and stores partial time
*/
void
DSMC::test_force_field
(void)
{
  if ( par_env->is_root() ) std::cout << "### TEST: computing force field ###" << std::endl;
  stopwatch.local_start(FORCES_TAG);
  mean_field->compute_force_field();
  par_env->barrier();
  stopwatch.local_stop(FORCES_TAG);
  stored_elapsed_times[FORCES_TAG].push_back(stopwatch.get_local_elapsed(FORCES_TAG));
  if ( par_env->is_root() ) stopwatch.show_local_elapsed(FORCES_TAG);
}

/*! \fn void DSMC::test_time_marching (void)
    \brief Tests an advection iteration and stores partial time
*/
void
DSMC::test_time_marching
(void)
{
  if ( par_env->is_root() ) std::cout << "### TEST: time-marching ###" << std::endl;
  stopwatch.local_start(ADVECT_TAG);
  time_marching->update_ensemble();
  par_env->barrier();
  stopwatch.local_stop(ADVECT_TAG);
  stored_elapsed_times[ADVECT_TAG].push_back(stopwatch.get_local_elapsed(ADVECT_TAG));
  if ( par_env->is_root() ) stopwatch.show_local_elapsed(ADVECT_TAG);
}

/*! \fn void DSMC::test_collisions (void)
    \brief Tests collisions kernel and stores partial time
*/
void
DSMC::test_collisions
(void)
{
  if ( par_env->is_root() ) std::cout << "### TEST: performing collisions ###" << std::endl;
  stopwatch.local_start(COLLISION_TAG);
  collision_handler->perform_collisions();
  par_env->barrier();
  stopwatch.local_stop(COLLISION_TAG);
  stored_elapsed_times[COLLISION_TAG].push_back(stopwatch.get_local_elapsed(COLLISION_TAG));
  if ( par_env->is_root() ) stopwatch.show_local_elapsed(COLLISION_TAG);
}

/*! \fn void DSMC::test_thermostat (void)
    \brief Tests a thermostat iteration
*/
void
DSMC::test_thermostat
(void)
{
  if ( par_env->is_root() ) std::cout << "### TEST: applying thermostat ###" << std::endl;
  thermostat->rescale_velocity();
}

/*! \fn void DSMC::test_sampling (void)
    \brief Tests sampling, averaging and gathering, stores partial time
*/
void
DSMC::test_sampling
(void)
{
  if ( par_env->is_root() ) std::cout << "### TEST: sampling ###" << std::endl;
  stopwatch.local_start(SAMPLING_TAG);
  sampler->reset();
  sampler->sample();
  sampler->average();
  sampler->gather_samples();
  stopwatch.local_stop(SAMPLING_TAG);
  stored_elapsed_times[SAMPLING_TAG].push_back(stopwatch.get_local_elapsed(SAMPLING_TAG));
  if ( par_env->is_root() ) stopwatch.show_local_elapsed(SAMPLING_TAG);
}

/*! \fn void DSMC::test_output (void)
    \brief Tests output

    Outputs (1) samples, (2) collisional statistics, (3) elapsed times
*/
void
DSMC::test_output
(void)
{
  if ( par_env->is_root() ) std::cout << "### TEST: output ###" << std::endl;
  output_all_samples();
  output_collision_statistics();
  output_elapsed_times();
}

/*! \fn void DSMC::test_loop (int n_tests)
    \brief Tests n_tests iterations of the DSMC loop
*/
void
DSMC::test_loop
(int n_tests)
{
  for (int i = 0; i<n_tests; ++i)
  {
    if ( par_env->is_root() ) std::cout << " >> test iteration: " << i << std::endl;
    if ( mean_field_gg )
      test_force_field();
    test_time_marching();
    if ( !SHUT_THERMO )
      test_thermostat();
    test_density();
    if ( !SHUT_COLLISIONS )
      test_collisions();
    test_sampling();
  }
  collision_handler->gather_collisions();
  test_output();
}


// DSMC SIMULATION STEPS

/*! \fn void DSMC::initialize_simulation (void)
    \brief Initialize density fielda and computes first estimate of coll. majorants
*/
void
DSMC::initialize_simulation
(void)
{
  if ( par_env->is_root() ) std::cout << "### INITIALIZINING DENSITY FIELD ###" << std::endl;
  density->perform_density_kernel();
  if ( !SHUT_COLLISIONS )
  {
    if ( par_env->is_root() ) std::cout << "### INITIALIZINING COLLISIONS MAJORANTS ###" << std::endl;
    collision_handler->compute_majorants();
  }
}

void
DSMC::dsmc_iteration
(void)
{
  // TO BE CONTINUED ...
}

void
DSMC::dsmc_loop
(void)
{
  // TO BE CONTINUED ...
}


// OUTPUT FUNCTIONALITIES

/*! \fn void DSMC::output_all_samples (void)
    \brief Outputs temperature, pressure and density; no time tag
*/
void
DSMC::output_all_samples
(void)
{
  if ( par_env->is_root() )
  {
    std::cout << "### OUTPUT ALL SAMPLES ###" << std::endl;
    output->output_sample(sampler->get_temp_avg(), "output_files/samples/test_sample_temperature.txt");
    output->output_sample(sampler->get_ph_avg(), "output_files/samples/test_sample_hpressure.txt");
    output->output_sample(sampler->get_numdens_avg(), "output_files/samples/test_sample_numdens.txt");
  }
}

/*! \fn void DSMC::output_all_samples (real_number t)
    \brief Outputs temperature, pressure and density; with time tag
*/
void
DSMC::output_all_samples
(real_number t)
{
  if ( par_env->is_root() )
  {
    std::cout << "### OUTPUT ALL SAMPLES ###" << std::endl;
    output->output_sample(sampler->get_temp_avg(), "output_files/samples/test_sample_temperature",t);
    output->output_sample(sampler->get_ph_avg(), "output_files/samples/test_sample_hpressure",t);
    output->output_sample(sampler->get_numdens_avg(), "output_files/samples/test_sample_numdens", t);
  }
}

/*! \fn void DSMC::output_collision_statistics (void)
    \brief Outputs collisions statistics (fake, real, total, out of bound)
*/
void
DSMC::output_collision_statistics
(void)
{
  if ( par_env->is_root() )
  {
    output->output_vector(collision_handler->get_n_fake_store(), "output_files/collisions_fake.txt");
    output->output_vector(collision_handler->get_n_real_store(), "output_files/collisions_real.txt");
    output->output_vector(collision_handler->get_n_total_store(), "output_files/collisions_total.txt");
    output->output_vector(collision_handler->get_n_out_store(), "output_files/collisions_out.txt");
  }
}


/*! \fn void DSMC::output_elapsed_times (void)
    \brief Outputs partial CPU times for each sub-routine
*/
void
DSMC::output_elapsed_times
(void)
{
  if ( par_env->is_root() )
  {
    output->output_vector(stored_elapsed_times[DENSITY_TAG], "output_files/times_density.txt");
    output->output_vector(stored_elapsed_times[FORCES_TAG], "output_files/times_forces.txt");
    output->output_vector(stored_elapsed_times[ADVECT_TAG], "output_files/times_advection.txt");
    output->output_vector(stored_elapsed_times[COLLISION_TAG], "output_files/times_collision.txt");
    output->output_vector(stored_elapsed_times[SAMPLING_TAG], "output_files/times_sampling.txt");
  }
}
