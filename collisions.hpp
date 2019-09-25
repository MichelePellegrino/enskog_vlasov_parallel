/*! \file collisions.hpp
 *  \brief Header containing class for collisions simulation
 */

#ifndef EV_COLLISIONS_HPP
#define EV_COLLISIONS_HPP

#include "motherbase.hpp"
#include "matrix.hpp"
#include "topology.hpp"

/* Valarray turns out to be handy when dealing to fixed-dimension array behaving
    like mathematical vectors */
#include <valarray>
#include <cmath>
#include <algorithm>
#include <vector>
#include <array>
#include <map>

/*! \def TEST_COEFF_MULT
    \brief Default value for the initial number of evaluation to estimate A_i and C_i
*/
#ifndef TEST_COEFF_MULT
#define TEST_COEFF_MULT 5
#endif

/*! \def DEFAULT_ALPHA_1
    \brief Parameter tuning out-of-bound/real ratio
*/
#ifndef DEFAULT_ALPHA_1
#define DEFAULT_ALPHA_1 1e-2
#endif

/*! \def DEFAULT_ALPHA_2
    \brief Parameter tuning the 'lowering' of estimated majorant
*/
#ifndef DEFAULT_ALPHA_2
#define DEFAULT_ALPHA_2 0.99
#endif

/*! \struct VelocityType
    \brief Struct storing velocity vector components (to be converted into MPI data type)
*/
struct VelocityType
{
  real_number vx, vy, vz;
  VelocityType() = default;
  VelocityType(real_number _vx_, real_number _vy_, real_number _vz_):
    vx(_vx_), vy(_vy_), vz(_vz_) { }
  ~VelocityType() = default;
};

/*! \struct CollisionCounter
    \brief Helper struct used to counter parallel collisions
*/
struct CollisionCounter
{
  int total, real, fake, out_bound;
  CollisionCounter() = default;
  CollisionCounter(int tt, int tr, int fl, int fk):
    total(tt), real(tr), fake(fl), out_bound(fk) { }
  ~CollisionCounter() = default;
};

// NB: the variable names used in this class are not at all intuitive (more
//     meaningful names may be employed).

/*! \class CollisionHandler
 *  \brief Class for collisions simulation
 *
 *  CollisionHandler wraps collisions counters and all methods to intialize majorants,
 *  to compute the expected no. of collisions and to simulate collisions
 */
class CollisionHandler : protected Motherbase
{

private:

  /* Definition of map < process, data >, used for parallel communication in the
     collisional phase */
  template <class data_type>
  using process_map = std::map< int, std::vector< data_type > >;

  // Parallel environment
  const int size, rank;

  // Global number of collisions
  int n_fake = 0;       /*!< Number of false collisions                 */
  int n_real = 0;       /*!< Number of true collisions                  */
  int n_total = 0;      /*!< Total number  collisions                   */
  int out_bound = 0;    /*!< True collisions for which probability > 1  */

  // Vectors for storage
  std::vector<int> n_fake_store;    /*!< Number of false collisions for each time-step        */
  std::vector<int> n_real_store;    /*!< Number of true collisions for each time-step         */
  std::vector<int> n_total_store;   /*!< Number of total collisions for each time-step        */
  std::vector<int> n_out_store;     /*!< Number of out-of-bound collisions for each time-step */

  /*!
   *  The order of collisional subdomains is selected at random, to avoid possible
   *  spurious correlations
   */
  std::array<int, N_COLLISION_SUBDOM> subdom_order;   /*!< Order of coll. subdomains */

  ev_matrix::MaskMatrix<real_number> a11;       /*!< Majorant A_i (density * correlation) */
  ev_matrix::MaskMatrix<real_number> vrmax11;   /*!< Majorant C_i (relative speed)        */

  // Computed at next iteration
  ev_matrix::MaskMatrix<real_number> anew;      /*!< Updated values for A_i */
  ev_matrix::MaskMatrix<real_number> vrmaxnew;  /*!< Updated values for C_i */

  int n_coll;                                   /*!< Total number of collisions         */
  ev_matrix::MaskMatrix<int> n_coll_cell;       /*!< Number of collisions for each cell */

  std::valarray<real_number> scaled_k;          /*!< Scaled vector pointing to target cell (sigma*k)  */
  std::valarray<real_number> rel_vel;           /*!< Relative velocity between particles              */
  std::valarray<real_number> delta;             /*!< Vector to be added/subtracted to update velocity */

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

  /*! \fn void CollisionHandler::gen_scaled_k(void)
      \brief Creates the vector poiting to the cell jc2 (scaled by sigma)
  */
  inline void gen_scaled_k(void)
  {
    rng->sample_unit_sphere(scaled_k[0], scaled_k[1], scaled_k[2]);
    scaled_k *= sigma;
  }

  // Controlling number of collisions
  const real_number alpha_1 = DEFAULT_ALPHA_1;  /*!< First coefficient for collision number control   */
  const real_number alpha_2 = DEFAULT_ALPHA_2;  /*!< Second coefficient for collision number control  */

  /*! \fn void CollisionHandler::update_majorants(void)
      \brief Updates majorants according to the predefined coefficients alpha_1, alpha_1
  */
  void update_majorants(void);

  /*! \fn std::pair<real_number, real_number> CollisionHandler::lexico_inv_quad(int idx, int q) const
      \brief Returns the local coordinates (x_q, y_q) inside the collisional quadrant
  */
  inline std::pair<real_number, real_number> lexico_inv_quad(int idx, int q) const
  {
    int j = idx / n_cells_x[q];
    int i = idx - j * n_cells_x[q];
    return std::make_pair(i, j);
  }

  /*! \fn std::pair<real_number, real_number> CollisionHandler::map_to_domain(int i_loc, int j_loc, int q) const
      \brief Maps the local coordinates in quandrant q into the global coordinates of the whole domain
  */
  inline std::pair<real_number, real_number> map_to_domain(int i_loc, int j_loc, int q) const
  {
    return std::make_pair(i_loc + low_x_quad[q], j_loc + low_y_quad[q]);
  }

  /*! \fn void CollisionHandler::init_cell_ind_map(void)
      \brief Initialize the vector containing the index for each cell
  */
  void init_cell_ind_map (void);

  /*! \fn void CollisionHandler::setup_cell_ind_map(void)
      \brief Set-up the vector containing the index for each cell

      Cells indices are randomly selected from this array; once the cell is selected
      the index is erased from the array
  */
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

  MPI_Datatype MPI_VELOCITY_TYPE;   /*!< MPI custom data type for velocity vectorial values */
  void commit_velocity_type(void);  /*!< Commitment of MPI custom data type for velocity    */

  CollisionCounter perform_parallel_collisions (void);
  void reset_parallel_buffers (void);

public:

  CollisionHandler(DSMC*);
  ~CollisionHandler() = default;

  // ...
  void shuffle_order(void);

  // Each step of collisional stage (including initialization)
  void compute_majorants(void);
  void compute_collision_number(void);
  void perform_collisions(void);

  // Collisional stage in a packet
  void gather_collisions(void);

  // DEBUG
  void print_subdomain_order(void) const;
  void print_est_collisions(void) const;
  void print_reduced_constants(void) const;
  void test_map_to_domain_collision(void);

  // GETTERS
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
