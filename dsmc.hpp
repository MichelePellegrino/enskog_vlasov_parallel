#ifndef EV_DSMC_HPP
#define EV_DSMC_HPP

#include "types.hpp"
#include "random.hpp"
#include "stopwatch.hpp"
#include "species.hpp"
#include "times.hpp"
#include "potential.hpp"
#include "correlations.hpp"

#include <cassert>
#include <map>
#include <vector>

#ifndef DEFAULT_ITER_THERMO
#define DEFAULT_ITER_THERMO 10
#endif

#ifndef DEFAULT_ITER_SAMPLE
#define DEFAULT_ITER_SAMPLE 100
#endif

#ifndef DEFAULT_DUMMY_ITER
#define DEFAULT_DUMMY_ITER 500
#endif

#ifndef DEFAULT_DUMMY_TEST_ITER
#define DEFAULT_DUMMY_TEST_ITER 10
#endif

#ifndef SHUT_COLLISIONS
#define SHUT_COLLISIONS false
#endif

#ifndef SHUT_THERMO
#define SHUT_THERMO false
#endif

#define DENSITY_TAG   0
#define FORCES_TAG    1
#define ADVECT_TAG    2
#define COLLISION_TAG 3
#define SAMPLING_TAG  4

#define N_TAGS        5

class ParallelEnvironment;
class IOHandler;
class ConfigurationReader;
class Boundary;
class Grid;
class Topology;
class Ensemble;
class Thermostat;
class DensityKernel;
class ForceField;
class CollisionHandler;
class Sampler;
class Output;

template<MarchingType tm_type> class TimeMarching;

class DSMC
{

public:
  typedef ev_random::CustomRngObject<RNG> RandomEngine;
  typedef ev_correlations::CorrelationFunction<CORR> CorrelationFun;

private:

  DefaultPointer<ParallelEnvironment> par_env;
  DefaultPointer<IOHandler> io_hand;
  DefaultPointer<ConfigurationReader> conf;

  DefaultPointer<RandomEngine> rng;

  DefaultPointer<Species> species;
  DefaultPointer<Times> times;
  DefaultPointer<Boundary> boundary;
  DefaultPointer<Grid> grid;

  DefaultPointer<Topology> topology;

  DefaultPointer<Ensemble> ensemble;
  DefaultPointer<Thermostat> thermostat;

  DefaultPointer<DensityKernel> density;
  DefaultPointer<NondirectionalPairPotential> potential;
  DefaultPointer<ForceField> mean_field;
  DefaultPointer<TimeMarching<TM>> time_marching;
  DefaultPointer<CollisionHandler> collision_handler;

  DefaultPointer<Sampler> sampler;
  DefaultPointer<Output> output;

  CorrelationFun correlation;

  Stopwatch<DefaultWatchPrecision> stopwatch;

  int n_iter_thermo = DEFAULT_ITER_THERMO;
  int n_iter_sample = DEFAULT_ITER_SAMPLE;
  bool mean_field_gg;

  std::map < int, std::vector<int> > stored_elapsed_times;

public:

  DSMC(const DefaultString&);
  ~DSMC() = default;

  // Getters
  inline DefaultPointer<ParallelEnvironment>& get_par_env() { return par_env; }
  inline DefaultPointer<IOHandler>& get_io_hand() { return io_hand; }
  inline DefaultPointer<ConfigurationReader>& get_conf() { return conf; }
  inline DefaultPointer<RandomEngine>& get_rng() { return rng; }
  inline DefaultPointer<Species>& get_species() { return species; }
  inline DefaultPointer<Times>& get_times() { return times; }
  inline DefaultPointer<Boundary>& get_boundary() { return boundary; }
  inline DefaultPointer<Grid>& get_grid() { return grid; }
  inline DefaultPointer<Topology>& get_topology() { return topology; }
  inline DefaultPointer<Ensemble>& get_ensemble() { return ensemble; }
  inline DefaultPointer<Thermostat>& get_thermostat() { return thermostat; }
  inline DefaultPointer<DensityKernel>& get_density() { return density; }
  inline DefaultPointer<NondirectionalPairPotential>& get_potential() { return potential; }
  inline DefaultPointer<ForceField>& get_mean_field() { return mean_field; }
  inline DefaultPointer<TimeMarching<TM>>& get_time_marching() { return time_marching; }
  inline DefaultPointer<CollisionHandler>& get_collision_handler() { return collision_handler; }
  inline DefaultPointer<Sampler>& get_sampler() { return sampler; }
  inline DefaultPointer<Output>& get_output() { return output; }

  inline CorrelationFun& get_correlation() { return correlation; }

  // Testing
  void test_density(void);
  void test_force_field(void);
  void test_time_marching(void);
  void test_collisions(void);
  void test_thermostat(void);
  void test_sampling(void);
  void test_output(void);
  void test_loop(int);

  // Initialize
  void initialize_simulation(void);
  void dsmc_iteration(void);
  void dsmc_loop(void);

  // Output
  void output_all_samples(void);
  void output_all_samples(real_number);
  void output_collision_statistics(void);
  void output_elapsed_times(void);

};

#endif /* DSMC_HPP */
