#ifndef TOPOLOGY_HPP
#define TOPOLOGY_HPP

#include "motherbase.hpp"
#include "matrix.hpp"
#include "configuration.hpp"

#include <vector>
#include <array>

#ifndef N_COLLISION_SUBDOM
#define N_COLLISION_SUBDOM 4
#endif

#ifndef DEFAULT_TOPOLOGY
#define DEFAULT_TOPOLOGY Quadrants
#endif

using namespace ev_matrix;

enum TopologyType
{
  StripesX,
  StripesY,
  Quadrants,
};

class Topology : protected Motherbase
{

private:

  const int rank;
  const int size;

  int n_cutoff_x;
  int n_cutoff_y;

  std::vector<int> idx_lx, idx_ly, idx_ux, idx_uy;
  int idx_lx_rank, idx_ly_rank, idx_ux_rank, idx_uy_rank;

  std::array<int, N_COLLISION_SUBDOM> quarter_lx, quarter_ly, quarter_ux, quarter_uy;
  std::array<int, N_COLLISION_SUBDOM> n_cells_quarter;

  MaskMatrix<int> topology_map;
  TopologyType topology_type = DEFAULT_TOPOLOGY;

  void fill_topology_map(void);
  void setup_stripes_x(void);
  void setup_stripes_y(void);
  void setup_quadrants(void);

  void setup_quarters(void);

public:

  Topology(DSMC*);
  ~Topology() = default;

  int tag_subdom(int, int);

  inline int get_idx_lx(int r) { return idx_lx[r]; }
  inline int get_idx_ly(int r) { return idx_ly[r]; }
  inline int get_idx_ux(int r) { return idx_ux[r]; }
  inline int get_idx_uy(int r) { return idx_uy[r]; }

  inline int get_idx_lx_rank() { return idx_lx_rank; }
  inline int get_idx_ly_rank() { return idx_ly_rank; }
  inline int get_idx_ux_rank() { return idx_ux_rank; }
  inline int get_idx_uy_rank() { return idx_uy_rank; }

  inline int get_quarter_lx(int q) { return quarter_lx[q]; }
  inline int get_quarter_ly(int q) { return quarter_ly[q]; }
  inline int get_quarter_ux(int q) { return quarter_ux[q]; }
  inline int get_quarter_uy(int q) { return quarter_uy[q]; }
  inline int get_n_cells_quarter(int q) { return n_cells_quarter[q]; }

  inline MaskMatrix<int>& get_topology_map(void) { return topology_map; }
  inline int get_topology_map(int i, int j) { return topology_map(i,j); }

  inline int get_n_cells(int r) { return (idx_ux[r]-idx_lx[r])*(idx_uy[r]-idx_ly[r]); }
  inline int get_n_cells() { return (idx_ux[rank]-idx_lx[rank])*(idx_uy[rank]-idx_ly[rank]); }

};

#endif /* TOPOLOGY_HPP */
