#ifndef EV_MOTHERBASE_HPP
#define EV_MOTHERBASE_HPP

#include "dsmc.hpp"

class Motherbase
{

public:
  typedef DSMC::RandomEngine RandomEngine;
  typedef DSMC::CorrelationFun CorrelationFun;

protected:

  DSMC* dsmc;

  DefaultPointer<ParallelEnvironment>& par_env;
  DefaultPointer<IOHandler>& io_hand;
  DefaultPointer<ConfigurationReader>& conf;

  DefaultPointer<RandomEngine>& rng;

  DefaultPointer<Species>& species;
  DefaultPointer<Times>& times;
  DefaultPointer<Boundary>& boundary;
  DefaultPointer<Grid>& grid;

  DefaultPointer<Topology>& topology;

  DefaultPointer<Ensemble>& ensemble;
  DefaultPointer<DensityKernel>& density;
  DefaultPointer<NondirectionalPairPotential>& potential;
  DefaultPointer<ForceField>& mean_field;
  DefaultPointer<TimeMarching<TM>>& time_marching;
  DefaultPointer<CollisionHandler>& collision_handler;

  CorrelationFun& correlation;

public:

  Motherbase(DSMC* _dsmc_):

    dsmc(_dsmc_),

    par_env           (dsmc->get_par_env()),
    io_hand           (dsmc->get_io_hand()),
    conf              (dsmc->get_conf()),

    rng               (dsmc->get_rng()),

    species           (dsmc->get_species()),
    times             (dsmc->get_times()),

    boundary          (dsmc->get_boundary()),
    grid              (dsmc->get_grid()),

    topology          (dsmc->get_topology()),

    ensemble          (dsmc->get_ensemble()),

    density           (dsmc->get_density()),
    potential         (dsmc->get_potential()),
    mean_field        (dsmc->get_mean_field()),
    time_marching     (dsmc->get_time_marching()),
    collision_handler (dsmc->get_collision_handler()),

    correlation       (dsmc->get_correlation())

    { }

  virtual ~Motherbase() = default;

};

#endif /* EV_MOTHERBASE_HPP */
