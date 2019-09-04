#include "grid.hpp"
#include "boundary.hpp"
#include "configuration.hpp"

// DEBUG
// # # # # #
#include "parallel_environment.hpp"
// # # # # #

#include <cassert>

Grid::Grid
(DSMC* dsmc):
  Motherbase(dsmc),
  n_cells_x(conf->get_n_cells_x()),
  n_cells_y(conf->get_n_cells_y()),
  n_cells(n_cells_x*n_cells_y),
  x_min(-conf->get_x_min()),
  x_max(conf->get_x_max()),
  y_min(-conf->get_y_min()),
  y_max(conf->get_y_max()),
  x_extra(conf->get_x_extra()),
  y_extra(conf->get_y_extra()),
  dx(conf->get_dx()),
  dy(conf->get_dy()),
  rdx(1.0/dx),
  rdy(1.0/dy),
  xc(n_cells_x),
  yc(n_cells_y)
  {

    // Coerence check x1
    if ( boundary->get_pct().find(boundary->get_wall_condition()[0]) == boundary->get_pct().end() )
    {
      assert(x_min==boundary->get_Lx1() && "In case of logical b.c. x_min and Lx1 have to be the same");
      x_lim_m = -boundary->get_Lx1();
    }
    else
    {
      assert(x_min!=boundary->get_Lx1() && "In case of physical b.c. x_min and Lx1 can't be equal");
      x_lim_m = -(boundary->get_Lx1() - species->get_diam_gw());
    }
    // Coerence check x2
    if ( boundary->get_pct().find(boundary->get_wall_condition()[2]) == boundary->get_pct().end() )
    {
      assert(x_max==boundary->get_Lx2() && "In case of logical b.c. x_min and Lx2 have to be the same");
      x_lim_p = boundary->get_Lx2();
    }
    else
    {
      assert(x_max!=boundary->get_Lx2() && "In case of physical b.c. x_min and Lx2 can't be equal");
      x_lim_p = boundary->get_Lx2() - species->get_diam_gw();
    }
    // Coerence check y1
    if ( boundary->get_pct().find(boundary->get_wall_condition()[1]) == boundary->get_pct().end() )
    {
      assert(y_min==boundary->get_Ly1() && "In case of logical b.c. x_min and Ly1 have to be the same");
      y_lim_m = -boundary->get_Ly1();
    }
    else
    {
      assert(y_min!=boundary->get_Ly1() && "In case of physical b.c. x_min and Ly1 can't be equal");
      y_lim_m = -(boundary->get_Ly1() - species->get_diam_gw());
    }
    // Coerence check y2
    if ( boundary->get_pct().find(boundary->get_wall_condition()[3]) == boundary->get_pct().end() )
    {
      assert(y_max==boundary->get_Ly2() && "In case of logical b.c. x_min and Ly2 have to be the same");
      y_lim_p = boundary->get_Ly2();
    }
    else
    {
      assert(y_max!=boundary->get_Ly2() && "In case of physical b.c. x_min and Ly2 can't be equal");
      y_lim_p = boundary->get_Ly2() - species->get_diam_gw();
    }
    // Set cell volume (basic case)
    /* TEMP: */
    // real_number channel_area = (x_lim_p-x_lim_m) * (y_lim_p-y_lim_m);
    real_number channel_area = (x_max-x_min) * (y_max-y_min);
    real_number dg = species->get_diam_fluid();
    /* TEMP: */
    real_number homogeneous_density = 6.0 * conf->get_eta_liq0() / (ev_const::pi * dg*dg*dg);
    // real_number homogeneous_density = 6.0 * conf->get_eta_liq1() / (ev_const::pi * dg*dg*dg);
    real_number channel_volume = conf->get_n_part() / homogeneous_density;
    channel_section = channel_volume / channel_area;
    cell_volume = dx*dy*channel_section;
    // Initialize centroids
    for ( int i = 0; i < n_cells_x; ++i)  xc[i] = x_min + ( i+0.5 ) * dx;
    for ( int j = 0; j < n_cells_y; ++j)  yc[j] = y_min + ( j+0.5 ) * dy;

    // DEBUG
    // # # # # #
    // print_info();
    if ( par_env->is_root() )
      std::cout << "### CHECK: n_part/(cell_volume*n_cells) = " << conf->get_n_part()/(cell_volume*n_cells) << std::endl;
    // # # # # #

  }


// DEBUG

void
Grid::print_info
(void) const
{
  for ( int r = 0; r<par_env->get_size(); ++r )
  {
    if ( r == par_env->get_rank() )
      std::cout << "rank " << r << " : cell volume = " << cell_volume << ";\t dx = " << dx << std::endl;
  }
  par_env->barrier();
}
