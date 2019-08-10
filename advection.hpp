#ifndef EV_ADVECTION_HPP
#define EV_ADVECTION_HPP

#include "types.hpp"
#include "motherbase.hpp"
#include "times.hpp"
#include "particles.hpp"
#include "grid.hpp"
#include "species.hpp"

class AbstractTimeMarching : protected Motherbase
{
protected:
  const real_number& delta_t;
  const int& n_particles;
  const real_number& xmin, xmax, ymin, ymax;
  const real_number& delta_x, delta_y;
  const real_number& mass;
public:
  AbstractTimeMarching(DSMC* _dsmc_):
    Motherbase(_dsmc_),
    delta_t( times->get_delta_t() ),
    n_particles( ensemble->get_n_particles() ),
    xmin( grid->get_x_min() ),
    xmax( grid->get_x_max() ),
    ymin( grid->get_y_min() ),
    ymax( grid->get_y_max() ),
    delta_x( grid->get_dx() ),
    delta_y( grid->get_dy() ),
    mass( species->get_mass_fluid() )
    { }
  virtual ~AbstractTimeMarching() = default;
  virtual void update_ensemble() = 0;
};

template <MarchingType dummy_rng_type>
class TimeMarching : public AbstractTimeMarching
{
public:
  TimeMarching(DSMC* _dsmc_):
    AbstractTimeMarching(_dsmc_)
    {
      throw "Invalid spacification for TimeMarching template";
    }
  ~TimeMarching() = default;
};

template <>
class TimeMarching<Standard> : public AbstractTimeMarching
{
public:
  TimeMarching<Standard>(DSMC* _dsmc_):
    AbstractTimeMarching(_dsmc_)
    {

    }
  ~TimeMarching<Standard>() = default;
  virtual void update_ensemble() override
  {
    real_number ax, ay;
    ensemble->clear_buffers();
    for ( int i = 0; i<n_particles; ++i )
    {
      /* GET FORCES */
      ax = mean_field->get_force_x( ensemble->data()[i].cell_x, ensemble->data()[i].cell_y ) / mass;
      ay = mean_field->get_force_y( ensemble->data()[i].cell_x, ensemble->data()[i].cell_y ) / mass;
      /* UPDATE POSITIONS */
      ensemble->data()[i].xp += ensemble->data()[i].vx * delta_t + 0.5 * ax * delta_t * delta_t;
      ensemble->data()[i].yp += ensemble->data()[i].vy * delta_t + 0.5 * ay * delta_t * delta_t;
      /* PERIODIC BOUNDARY CONDITIONS ... */
      while ( ensemble->data()[i].xp >= xmax )  ensemble->data()[i].xp += xmin - xmax;
      while ( ensemble->data()[i].xp < xmin )   ensemble->data()[i].xp += xmax - xmin;
      while ( ensemble->data()[i].yp >= ymax )  ensemble->data()[i].yp += ymin - ymax;
      while ( ensemble->data()[i].yp < ymin )   ensemble->data()[i].yp += ymax - ymin;
      /* UPDATE VELOCITIES */
      ensemble->data()[i].vx += ax * delta_t;
      ensemble->data()[i].vy += ay * delta_t;
      /* UPDATE CELL */
      ensemble->data()[i].cell_x = (int)( (ensemble->data()[i].xp - xmin) / delta_x );
      ensemble->data()[i].cell_y = (int)( (ensemble->data()[i].yp - ymin) / delta_y );
      /* UPDATE RANK */
      ensemble->data()[i].r_tag = topology->tag_subdom(ensemble->data()[i].cell_x, ensemble->data()[i].cell_y);
      /* INCIDENCE MATRIX */
      ensemble->add_to_inc_matrix(i);
    }
    par_env->barrier();
    ensemble->exchange_particles();
  }
};

#endif /* EV_ADVECTION_HPP */
