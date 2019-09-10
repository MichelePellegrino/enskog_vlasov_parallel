#ifndef EV_PARTICLES_HPP
#define EV_PARTICLES_HPP

#include "motherbase.hpp"

#include <vector>
#include <map>

/*! \struct Particle
 *  \brief A struct for a single particle
 */
struct Particle
{
  real_number xp, yp;
  real_number vx, vy, vz;
  int cell_x, cell_y;
  int p_tag;
  int r_tag;
};

// NB CLEAN THE CODE AND REMOVE UNUSED BUFFERS

/*! \class Ensemble
 *  \brief A class for an ensemble of particles
 */
class Ensemble : protected Motherbase
{

private:

  int n_particles;                    /*!< Total number of particles */
  real_number vx_ini, vy_ini, vz_ini;
  real_number T_ini;
  real_number mass;

  int domain_rank;                    /*!< Rank of this process */
  int n_subdom;                       /*!< Number of subdomains */

  std::vector<Particle> particles;    /*!< Local vector of particles */
  MPI_Datatype MPI_PARTICLE_TYPE;     /*!< MPI custom datatype for a Particle */

  void commit_particle_type(void);

  /*! \fn populate
   *  \brief Randomly populate the phase space
   */
  void populate(void);

  std::vector<int> inc_matrix;        /*!< Incidence matrix (may be sparse) */
  int inc_particles = 0;              /*!< ... */
  int dep_particles = 0;              /*!< ... */
  int idx_r = 0;                      /*!< ... */
  // MAP ( j, [a_1j, a_2j, ...] )
  // BEING 'j' THE INDEX OF THE PROCESS WHERE TO SEND PARTICLES
  std::map< int, std::vector<Particle> > temp_send_map;

  /*! \fn fill_inc_matrix
   *  \brief Fills the incidence matrix
   */
  void fill_inc_matrix(void);

public:

  Ensemble(DSMC*);
  ~Ensemble() = default;

  // void test_stream(real_number dt = 0.1);

  /* UTILITIES */
  // void save_to_file(const std::string&) const;

  /*! \fn clear_buffers
   *  \brief Clears send/receive buffers
   */
  void clear_buffers(void);

  /*! \fn
   *  \brief [...]
   */
  void add_to_inc_matrix(int);

  /*! \fn exchange_particles
   *  \brief Exchage particles between processes, based on their position
   */
  void exchange_particles(void);

  // Parameter getters
  inline const int& get_n_particles(void) const { return n_particles; }

  // Data getters
  inline std::vector<Particle>& data(void) { return particles; }
  inline const std::vector<Particle>& data(void) const { return particles; }

  inline const real_number get_xp(int k) const { return particles[k].xp; }
  inline real_number& get_xp(int k) { return particles[k].xp; }
  inline const real_number get_yp(int k) const { return particles[k].yp; }
  inline real_number& get_yp(int k) { return particles[k].yp; }

  inline const real_number get_vx(int k) const { return particles[k].vx; }
  inline real_number& get_vx(int k) { return particles[k].vx; }
  inline const real_number get_vy(int k) const { return particles[k].vy; }
  inline real_number& get_vy(int k) { return particles[k].vy; }
  inline const real_number get_vz(int k) const { return particles[k].vz; }
  inline real_number& get_vz(int k) { return particles[k].vz; }

  inline int get_cx(int k) const { return particles[k].cell_x; }
  inline int get_cy(int k) const { return particles[k].cell_y; }

};


#endif /* EV_PARTICLES_HPP */
