/*! \file thermostat.hpp
 *  \brief Header containing a class implementing the thermostat
 */

#ifndef EV_THERMOSTAT_HPP
#define EV_THERMOSTAT_HPP

#include "motherbase.hpp"
#include <valarray>
#include <cmath>

class Thermostat : protected Motherbase
{
private:
  int niter_thermo;
  real_number T_ref;
  std::valarray<real_number> v_tmp;
  real_number t_tmp;
  const int& n_part;
public:
  Thermostat(DSMC*);
  ~Thermostat() = default;
  void rescale_velocity(void);
  int get_niter_thermo() const { return niter_thermo; }
  real_number get_T_ref() const { return T_ref; }
};

#endif /* EV_THERMOSTAT_HPP */
