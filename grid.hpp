#ifndef EV_GRID_HPP
#define EV_GRID_HPP

#include <valarray>
#include <utility>

#include "motherbase.hpp"

class Grid : protected Motherbase
{

private:

  int n_cells_x, n_cells_y;                         /*!< Number of computational cells in each direction */
  int n_cells;
  real_number x_min, x_max, y_min, y_max;           /*!< Domain limits (whole 'logical' domain) */
  real_number x_extra, y_extra;                     /*!< Extra buffer between 'logical' and 'physical' domain (dimension) */
  real_number dx, dy;                               /*!< Meshwidth (and meshwidth^-1) */
  real_number rdx, rdy;
  real_number x_lim_p, x_lim_m, y_lim_p, y_lim_m;   /*!< Width of the buffer between walls and last position (dep. on boundary conditions) */
  real_number channel_section;
  real_number cell_volume;                          /*!< Volume of the computational cell */
  std::valarray<real_number> xc, yc;                /*!< Centroids */

public:

  Grid (DSMC*);
  ~Grid() = default;

  inline int lexico(int i, int j) const
  {
    return i + j * n_cells_x;
  }

  inline std::pair<real_number, real_number> lexico_inv(int idx) const
  {
    int j = idx / n_cells_x;
    int i = idx - j * n_cells_x;
    return std::make_pair(i, j);
  }

  // 'Getter' (and implicitly 'setter') methods ...

  inline const real_number& get_dx(void) const { return dx; }
  inline const real_number& get_dy(void) const { return dy; }
  inline const real_number& get_rdx(void) const { return rdx; }
  inline const real_number& get_rdy(void) const { return rdy; }

  inline const int& get_n_cells_x(void) const { return n_cells_x; }
  inline const int& get_n_cells_y(void) const { return n_cells_y; }
  inline int get_n_cells(void) const { return n_cells; }

  inline real_number get_x_extra(void) const { return x_extra; }
  inline real_number get_y_extra(void) const { return y_extra; }
  inline const real_number& get_x_min(void) const { return x_min; }
  inline const real_number& get_x_max(void) const { return x_max; }
  inline real_number get_y_min(void) const { return y_min; }
  inline real_number get_y_max(void) const { return y_max; }

  inline real_number get_cell_volume(void) const { return cell_volume; }

  inline real_number get_xlp(void) const { return x_lim_p; }
  inline real_number get_xlm(void) const { return x_lim_m; }
  inline real_number get_ylp(void) const { return y_lim_p; }
  inline real_number get_ylm(void) const { return y_lim_m; }

  inline const real_number get_xc(int i) const { return xc[i]; }
  inline const real_number get_yc(int j) const { return yc[j]; }

};


#endif /* EV_GRID_HPP */
