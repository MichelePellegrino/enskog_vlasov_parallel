#ifndef EV_TIMES_HPP
#define EV_TIMES_HPP

#include <vector>

class Times
{

private:
  int nc_out, restart;
  real_number t_ini, t_max, t_im, delta_t;        // What are tempo and tc_end ?
  std::vector<real_number> tc;

public:

  Times( int _nc_out_, int _restart_,
    real_number _t_ini_, real_number _t_max_, real_number _t_im_, real_number _delta_t_ ):
    nc_out(_nc_out_),
    restart(_restart_),
    t_ini(_t_ini_),
    t_max(_t_max_),
    t_im(_t_im_),
    delta_t(_delta_t_),
    tc(nc_out+1, 0.0)
    {
      for (int i=0; i<nc_out+1; i++)
        tc[i] = t_im + (t_max - t_im) * (real_number)i / (real_number)nc_out;
    }

  ~Times() = default;

  const int get_nc_out(void) const { return nc_out; }
  const int get_restart(void) const { return restart; }

  const real_number get_t_ini(void) const { return t_ini; }
  const real_number get_t_max(void) const { return t_max; }
  const real_number get_t_im(void) const { return t_im; }

  const real_number& get_delta_t(void) const { return delta_t; }

  const real_number get_tc(int i) const { return tc[i]; }

};

#endif /* EV_TIMES_HPP */
