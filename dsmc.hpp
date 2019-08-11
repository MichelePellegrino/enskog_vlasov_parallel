#ifndef EV_DSMC_HPP
#define EV_DSMC_HPP

#include "types.hpp"
#include "random.hpp"
#include "species.hpp"
#include "times.hpp"
#include "potential.hpp"
#include "correlations.hpp"

#include <cassert>

class ParallelEnvironment;
class IOHandler;
class ConfigurationReader;
class Boundary;
class Grid;
class Topology;
class Ensemble;
class DensityKernel;
class ForceField;
class CollisionHandler;

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
  DefaultPointer<DensityKernel> density;
  DefaultPointer<NondirectionalPairPotential> potential;
  DefaultPointer<ForceField> mean_field;
  DefaultPointer<TimeMarching<TM>> time_marching;
  DefaultPointer<CollisionHandler> collision_handler;

  CorrelationFun correlation;

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
  inline DefaultPointer<DensityKernel>& get_density() { return density; }
  inline DefaultPointer<NondirectionalPairPotential>& get_potential() { return potential; }
  inline DefaultPointer<ForceField>& get_mean_field() { return mean_field; }
  inline DefaultPointer<TimeMarching<TM>>& get_time_marching() { return time_marching; }
  inline DefaultPointer<CollisionHandler>& get_collision_handler() { return collision_handler; }
  inline CorrelationFun& get_correlation() { return correlation; }

};

#endif /* DSMC_HPP */
