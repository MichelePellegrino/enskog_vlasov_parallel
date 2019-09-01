#include "thermostat.hpp"
#include "configuration.hpp"
#include "particles.hpp"

Thermostat::Thermostat(DSMC* dsmc):
  Motherbase(dsmc),
  niter_thermo( conf->get_niter_thermo() ),
  T_ref( conf->get_T_ref() ),
  v_tmp(0.0, 3),
  t_tmp(0.0),
  n_part( ensemble->get_n_particles() )
  { }

void
Thermostat::rescale_velocity
(void)
{

  v_tmp = 0.0;
  t_tmp = 0.0;

  real_number vx, vy, vz;

  for (int i = 0; i<n_part; ++i)
  {
    vx = ensemble->get_vx(i);
    vy = ensemble->get_vy(i);
    vz = ensemble->get_vz(i);
    v_tmp[0] += vx;
    v_tmp[1] += vy;
    v_tmp[2] += vz;
    t_tmp += vx*vx + vy*vy + vz*vz;
  }

  v_tmp /= (double)n_part;
  t_tmp = ( t_tmp/(double)n_part -
    (v_tmp[0]*v_tmp[0] + v_tmp[1]*v_tmp[1] + v_tmp[2]*v_tmp[2]) ) / 3.0;

  t_tmp = sqrt(t_tmp/T_ref);

  for (int i = 0; i<n_part; ++i)
  {
    ensemble->get_vx(i) /= t_tmp;
    ensemble->get_vy(i) /= t_tmp;
    ensemble->get_vz(i) /= t_tmp;
  }

}
