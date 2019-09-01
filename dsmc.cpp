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

  initialize_simulation();
  par_env->barrier();

  // TESTING MODULES
  // test_sampling();
  // test_output();

  // TESTING LOOP
  test_loop(DEFAULT_DUMMY_TEST_ITER);

}


// Testing
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

void
DSMC::test_thermostat
(void)
{
  if ( par_env->is_root() ) std::cout << "### TEST: applying thermostat ###" << std::endl;
  thermostat->rescale_velocity();
}

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

void
DSMC::test_output
(void)
{
  if ( par_env->is_root() ) std::cout << "### TEST: output ###" << std::endl;
  output_all_samples();
  output_collision_statistics();
  output_elapsed_times();
}

void
DSMC::test_loop
(int n_tests)
{
  for (int i = 0; i<n_tests; ++i)
  {
    if ( par_env->is_root() ) std::cout << " >> test iteration: " << i << std::endl;
    // DEBUG
    // # # # # #
    test_force_field();
    // # # # # #
    test_time_marching();
    test_density();
    test_collisions();
    test_thermostat();
    test_sampling();
  }
  collision_handler->gather_collisions();
  test_output();
}

// Initialize
void
DSMC::initialize_simulation
(void)
{
  if ( par_env->is_root() ) std::cout << "### INITIALIZINING DENSITY FIELD ###" << std::endl;
  density->perform_density_kernel();
  if ( par_env->is_root() ) std::cout << "### INITIALIZINING COLLISIONS MAJORANTS ###" << std::endl;
  collision_handler->compute_majorants();
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

void
DSMC::output_all_samples
(void)
{
  if ( par_env->is_root() )
  {
    std::cout << "### OUTPUT ALL SAMPLES ###" << std::endl;
    output->output_sample(sampler->get_vx_avg(), "output_files/samples/test_sample_vx.txt");
    output->output_sample(sampler->get_vy_avg(), "output_files/samples/test_sample_vy.txt");
    output->output_sample(sampler->get_vz_avg(), "output_files/samples/test_sample_vz.txt");
    output->output_sample(sampler->get_temp_avg(), "output_files/samples/test_temp_avg.txt");
    output->output_sample(sampler->get_pxx_avg(), "output_files/samples/test_sample_pxx.txt");
    output->output_sample(sampler->get_pyy_avg(), "output_files/samples/test_sample_pyy.txt");
    output->output_sample(sampler->get_pzz_avg(), "output_files/samples/test_sample_pzz.txt");
    output->output_sample(sampler->get_pxy_avg(), "output_files/samples/test_sample_pxy.txt");
    output->output_sample(sampler->get_pxz_avg(), "output_files/samples/test_sample_pxz.txt");
    output->output_sample(sampler->get_pyz_avg(), "output_files/samples/test_sample_pyz.txt");
    output->output_sample(sampler->get_qx_avg(), "output_files/samples/test_sample_qx.txt");
    output->output_sample(sampler->get_qy_avg(), "output_files/samples/test_sample_qy.txt");
    output->output_sample(sampler->get_qz_avg(), "output_files/samples/test_sample_qz.txt");
    output->output_sample(sampler->get_numdens_avg(), "output_files/samples/test_sample_numdens.txt");
  }
}

void
DSMC::output_all_samples
(real_number t)
{
  if ( par_env->is_root() )
  {
    std::cout << "### OUTPUT ALL SAMPLES ###" << std::endl;
    output->output_sample(sampler->get_vx_avg(), "output_files/samples/test_sample_vx", t);
    output->output_sample(sampler->get_vy_avg(), "output_files/samples/test_sample_vy", t);
    output->output_sample(sampler->get_vz_avg(), "output_files/samples/test_sample_vz", t);
    output->output_sample(sampler->get_temp_avg(), "output_files/samples/test_temp_avg", t);
    output->output_sample(sampler->get_pxx_avg(), "output_files/samples/test_sample_pxx", t);
    output->output_sample(sampler->get_pyy_avg(), "output_files/samples/test_sample_pyy", t);
    output->output_sample(sampler->get_pzz_avg(), "output_files/samples/test_sample_pzz", t);
    output->output_sample(sampler->get_pxy_avg(), "output_files/samples/test_sample_pxy", t);
    output->output_sample(sampler->get_pxz_avg(), "output_files/samples/test_sample_pxz", t);
    output->output_sample(sampler->get_pyz_avg(), "output_files/samples/test_sample_pyz", t);
    output->output_sample(sampler->get_qx_avg(), "output_files/samples/test_sample_qx", t);
    output->output_sample(sampler->get_qy_avg(), "output_files/samples/test_sample_qy", t);
    output->output_sample(sampler->get_qz_avg(), "output_files/samples/test_sample_qz", t);
    output->output_sample(sampler->get_numdens_avg(), "output_files/samples/test_sample_numdens", t);
  }
}

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
