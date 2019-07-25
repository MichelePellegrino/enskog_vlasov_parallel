#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include "motherbase.hpp"
#include "matrix.hpp"
#include "configuration.hpp"

#include <vector>

using namespace ev_matrix;

enum TopologyType
{
  StripesX,
  StripesY,
  Quarters,
};

class Topology : protected Motherbase
{

private:

  const int rank;
  const int size;

  int n_cutoff_x;
  int n_cutoff_y;

  std::vector<int> idx_lx, idx_ly, idx_ux, idx_uy;

  std::vector<real_number> xmin_sub, ymin_sub, xmax_sub, ymax_sub;

  MaskMatrix<int> topology_map;
  TopologyType topology_type = Quarters;

  void fill_topology_map();
  void fill_limits();

public:

  Topology(DSMC*);
  ~Topology() = default;

  int tag_subdom(int, int);

  inline int get_idx_lx(int r) { return idx_lx[r]; }
  inline int get_idx_ly(int r) { return idx_ly[r]; }
  inline int get_idx_ux(int r) { return idx_ux[r]; }
  inline int get_idx_uy(int r) { return idx_uy[r]; }

  inline int get_xmin_sub(int r) { return xmin_sub[r]; }
  inline int get_ymin_sub(int r) { return ymin_sub[r]; }
  inline int get_xmax_sub(int r) { return xmax_sub[r]; }
  inline int get_ymax_sub(int r) { return ymax_sub[r]; }

  MaskMatrix<int>& get_topology_map(void) { return topology_map; }

  inline int get_n_cells(int r) { return (idx_ux[r]-idx_lx[r])*(idx_uy[r]-idx_ly[r]); }

};

#endif /* TOPOLOGY_HPP */
