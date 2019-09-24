/*! \file advection.hpp
 *  \brief Header containing the classes for ensemble advection
 *
 *  Static-time polimorphism allows the definition of multiple advection schemes
 *  via template specialization.
 */

// NB:  no other time marching techniques have been defined apart from the
//      standard one; it would be better to refactor this class hierarchy.


#ifndef EV_ADVECTION_HPP
#define EV_ADVECTION_HPP

#include "types.hpp"
#include "motherbase.hpp"
#include "times.hpp"
#include "particles.hpp"
#include "grid.hpp"
#include "species.hpp"


/*! \class AbstractTimeMarching
 *  \brief Class containing all data needed by any advection scheme
 */
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


/*! \class TimeMarching
 *  \brief Base for advection classes hierarchy
 *
 *  This template class needs to be specialized in order to produce a valid time
 *  -marching scheme.
 */
template <MarchingType dummy_advection_type>
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


/*! \class TimeMarching<Standard>
 *  \brief Standard advection scheme
 *
 *  This is the standard explicit foreard advection scheme adopted by standard
 *  DSMC simulations.
 */
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
    // DEBUG
    // # # # # #
    int counter_out_of_bound;
    // # # # # #
    real_number ax, ay;
    ensemble->clear_buffers();
    for ( int i = 0; i<n_particles; ++i )
    {
      // DEBUG
      // # # # # #
      counter_out_of_bound = 0;
      // # # # # #
      /* GET FORCES */
      ax = mean_field->get_force_x( ensemble->data()[i].cell_x, ensemble->data()[i].cell_y ) / mass;
      ay = mean_field->get_force_y( ensemble->data()[i].cell_x, ensemble->data()[i].cell_y ) / mass;
      /* UPDATE POSITIONS */
      ensemble->data()[i].xp += ensemble->data()[i].vx * delta_t + 0.5 * ax * delta_t * delta_t;
      ensemble->data()[i].yp += ensemble->data()[i].vy * delta_t + 0.5 * ay * delta_t * delta_t;
      /* PERIODIC BOUNDARY CONDITIONS ... */
      // In the following it is tested whether a particle go through the domain
      // more than 10 times (the number is arbitrary; it is just a rough test to
      // ensure velocities do not 'explode').
      while ( ensemble->data()[i].xp >= xmax )
      {
        counter_out_of_bound++;
        ensemble->data()[i].xp += xmin - xmax;
        assert(counter_out_of_bound<10 && "WTF!");
      }
      while ( ensemble->data()[i].xp < xmin )
      {
        counter_out_of_bound++;
        ensemble->data()[i].xp += xmax - xmin;
        assert(counter_out_of_bound<10 && "WTF!");
      }
      while ( ensemble->data()[i].yp >= ymax )
      {
        counter_out_of_bound++;
        ensemble->data()[i].yp += ymin - ymax;
        assert(counter_out_of_bound<10 && "WTF!");
      }
      while ( ensemble->data()[i].yp < ymin )
      {
        counter_out_of_bound++;
        ensemble->data()[i].yp += ymax - ymin;
        assert(counter_out_of_bound<10 && "WTF!");
      }
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
    /*!
    *   Particles need to be exchanged between processes AFTER advection (this
    *   is one advantage of the replicated grid approach).
    */
    if ( par_env->is_root() )
      std::cout << " >> exchanging particles..." << std::endl;
    ensemble->exchange_particles();
  }
};

#endif /* EV_ADVECTION_HPP */
