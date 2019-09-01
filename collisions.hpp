#ifndef EV_COLLISIONS_HPP
#define EV_COLLISIONS_HPP

#include "motherbase.hpp"
#include "matrix.hpp"
#include "topology.hpp"

#include <valarray>
#include <cmath>
#include <algorithm>
#include <vector>
#include <array>
#include <map>

#ifndef TEST_COEFF_MULT
#define TEST_COEFF_MULT 5
#endif

#ifndef DEFAULT_ALPHA_1
#define DEFAULT_ALPHA_1 1e-2
#endif

#ifndef DEFAULT_ALPHA_2
#define DEFAULT_ALPHA_2 0.99
#endif

struct VelocityType
{
  real_number vx, vy, vz;
  VelocityType() = default;
  VelocityType(real_number _vx_, real_number _vy_, real_number _vz_):
    vx(_vx_), vy(_vy_), vz(_vz_) { }
  ~VelocityType() = default;
};

struct CollisionCounter
{
  int total, real, fake, out_bound;
  CollisionCounter() = default;
  CollisionCounter(int tt, int tr, int fl, int fk):
    total(tt), real(tr), fake(fl), out_bound(fk) { }
  ~CollisionCounter() = default;
};

class CollisionHandler : protected Motherbase
{

private:

  template <class data_type>
  using process_map = std::map< int, std::vector< data_type > >;

  // Parallel env
  const int size, rank;

  // Global number of collisions
  int n_fake = 0;
  int n_real = 0;
  int n_total = 0;

  // Vectors for storage
  std::vector<int> n_fake_store;
  std::vector<int> n_real_store;
  std::vector<int> n_total_store;
  std::vector<int> n_out_store;

  std::array<int, N_COLLISION_SUBDOM> subdom_order;

  ev_matrix::MaskMatrix<real_number> a11;
  ev_matrix::MaskMatrix<real_number> vrmax11;

  // Computed at next iteration
  ev_matrix::MaskMatrix<real_number> anew;
  ev_matrix::MaskMatrix<real_number> vrmaxnew;

  int n_coll;
  ev_matrix::MaskMatrix<int> n_coll_cell;

  std::valarray<real_number> scaled_k;
  std::valarray<real_number> rel_vel;
  std::valarray<real_number> delta;

  // Index for choosing cells within collisional subdomain
  std::map< int, std::vector<int> > cells_ind;
  std::array<int, N_COLLISION_SUBDOM> n_cells_x;
  std::array<int, N_COLLISION_SUBDOM> low_x_quad;
  std::array<int, N_COLLISION_SUBDOM> low_y_quad;

  // References
  const int& npart;
  const real_number& xmin, xmax, ymin, ymax;
  const real_number& rdx, rdy;
  const real_number& sigma, delta_t;

  // Utilities
  real_number vr = 0.0, scalar_prod = 0.0;
  const real_number lx_sub, ux_sub, ly_sub, uy_sub;

  inline void gen_scaled_k(void)
  {
    rng->sample_unit_sphere(scaled_k[0], scaled_k[1], scaled_k[2]);
    scaled_k *= sigma;
  }

  // Controlling number of collisions
  const real_number alpha_1 = DEFAULT_ALPHA_1;
  const real_number alpha_2 = DEFAULT_ALPHA_2;
  void update_majorants(void);

  inline std::pair<real_number, real_number> lexico_inv_quad(int idx, int q) const
  {
    int j = idx / n_cells_x[q];
    int i = idx - j * n_cells_x[q];
    return std::make_pair(i, j);
  }

  inline std::pair<real_number, real_number> map_to_domain(int i_loc, int j_loc, int q) const
  {
    return std::make_pair(i_loc + low_x_quad[q], j_loc + low_y_quad[q]);
  }

  void init_cell_ind_map (void);
  void setup_cell_ind_map (void);

  // PARALLEL COMMUNICATION
  /* IT MAY BE SIMPLIFIED, MOST PROBABLY */
  std::map< int, int > n_parallel_collisions_send;
  std::map< int, int > n_parallel_collisions_recv;
  process_map< int > cells1_process_map;
  process_map< int > particle1_process_map;
  process_map< std::valarray<real_number> >  scaled_k_process_map;
  process_map< int > cells2_process_map_send;
  process_map< int > cells2_process_map_recv;
  process_map< int > hcells2_process_map_send;
  process_map< int > hcells2_process_map_recv;
  process_map< int > numdens_process_map_send;
  process_map< int > numdens_process_map_recv;
  process_map< real_number > aveta_mid_process_map_send;
  process_map< real_number > aveta_mid_process_map_recv;
  process_map< VelocityType > velocity_process_map_send;
  process_map< VelocityType > velocity_process_map_recv;
  process_map< VelocityType > delta_process_map_send;
  process_map< VelocityType > delta_process_map_recv;
  process_map< int > perform_collision_send;
  process_map< int > perform_collision_recv;

  MPI_Datatype MPI_VELOCITY_TYPE;
  void commit_velocity_type(void);

  CollisionCounter perform_parallel_collisions (void);
  void reset_parallel_buffers (void);

public:

  CollisionHandler(DSMC*);
  ~CollisionHandler() = default;

  void shuffle_order(void);

  void compute_majorants(void);
  void compute_collision_number(void);
  void perform_collisions(void);

  void gather_collisions(void);

  // DEBUG
  void print_subdomain_order(void) const;
  void print_est_collisions(void) const;
  void print_reduced_constants(void) const;
  void test_map_to_domain_collision(void);

  inline std::vector<int>& get_n_fake_store(void) { return n_fake_store; }
  inline const std::vector<int>& get_n_fake_store(void) const { return n_fake_store; }
  inline std::vector<int>& get_n_real_store(void) { return n_real_store; }
  inline const std::vector<int>& get_n_real_store(void) const { return n_real_store; }
  inline std::vector<int>& get_n_total_store(void) { return n_total_store; }
  inline const std::vector<int>& get_n_total_store(void) const { return n_total_store; }
  inline std::vector<int>& get_n_out_store(void) { return n_out_store; }
  inline const std::vector<int>& get_n_out_store(void) const { return n_out_store; }

};

#endif /* EV_COLLISIONS_HPP */
