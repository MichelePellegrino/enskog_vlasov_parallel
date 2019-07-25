#ifndef EV_CONFIGURATION_HPP
#define EV_CONFIGURATION_HPP

#include "motherbase.hpp"

#include <string>
#include <array>

class DSMC;

class ConfigurationReader : protected Motherbase
{

private:

  // File name
  DefaultString conf_file_name;             // Name of the configuration file

  // Tags
  DefaultString numtest;                    // ORIGINAL: Numtest (numero sequenziale della prova)

  // Molecular properties
  real_number mass_fluid;                   // ORIGINAL: Massa1 (massa molecolare del liquido)
  real_number diam_fluid;                   // ORIGINAL: Diamol1 (diametro molecolare del liquido)
  char pot_gas;                             // ORIGINAL: pot_gas (tipo di potenziale per gas-gas)
  char mean_f_gg;                           // ORIGINAL: c_med_gg (campo medio gas-gas [yes/no])
  real_number phi11;                        // ORIGINAL: Phi11 (parametro potenziale PHI11)
  real_number gamma11;                      // ORIGINAL: Gamma11 (esponente potenziale PHI11 [p])
  // real_number alpha11;                   // ORIGINAL: alpha11 (    "     "   [e])
  // Numero punti per spline pot gas-gas
  real_number mass_solid;                   // ORIGINAL: Massa2 (massa molecolare solido)
  real_number diam_solid;                   // ORIGINAL: Diamol2 (diametro molecolare wall)
  char pot_wall;                            // ORIGINAL: pot_wall (tipo di potenziale gas-wall)
  char mean_f_gw;                           // ORIGINAL: c_med_gw (campo medio gas-wall [yes/no])
  real_number phi12;                        // ORIGINAL: Phi12 (parametro potenziale PHI12)
  real_number gamma12;                      // ORIGINAL: Gamma12 (esponente potenziale PHI12 [p])
  // real_number alpha12;                   // ORIGINAL: alpha12 (    "	    "		[e])
  // Numero punti per spline pot gas-sol

  // Phisical state of fluid / walls
  real_number T_ini;                        // ORIGINAL: T_ini (temperatura iniziale del liquido)
  real_number T_ini1;                       // ORIGINAL: T_ini1	---> (?????)
  real_number eta_w0;                       // ORIGINAL: Etaw0 (densità della parete, per campo medio)
  real_number eta_w1;                       // ORIGINAL: Etaw1 (densità della parete, per correlazione)
  real_number eta_liq0;                     // ORIGINAL: Etaliq0 (densità del gas)
  real_number eta_liq1;                     // ORIGINAL: Etaliq1 (densità della fase liquida)
  int liq_interf;                           // ORIGINAL: liq_interf	---> (?????)
  real_number x_liq_interf;		              // ORIGINAL: x_liq_interf ---> (?????)
  real_number y_liq_interf;			            // ORIGINAL: y_liq_interf ---> (?????)
  real_number r_liq_interf;				          // ORIGINAL: r_liq_interf ---> (?????)
  real_number liq_drop_offset;              // ORIGINAL: liq_drop_offset ---> (?????)
  int niter_thermo;			                    // ORIGINAL: niter_thermo (number of iterations for thermostat)
  real_number T_ref;					              // ORIGINAL: T_ref (reference temperature thermostat)
  bool fix_baricentre;				              // ORIGINAL: bari (keep baricentre fixed or not [1/0])
  std::array<char, 4> wall_cond;            // ORIGINAL: natura del bordo x1
  bool set_px;                              // ORIGINAL: set_px	(sets/ignore periodicity in mean f x [1/0])
  bool set_py;                              // ORIGINAL: set_py	(			"		 "		f y [1/0])
  bool set_eta_px;                          // ORIGINAL: set_eta_px	(sets/ignore periodicity in avdens x [1/0])
  bool set_eta_py;                          // ORIGINAL: set_eta_py (   "		"		avdens y [1/0])
  std::array<real_number, 4> T_w;           // ORIGINAL: temperatura parete
  real_number f_x_1;                        // ORIGINAL: f_x_1 ---> (?????)
  real_number f_x_2;                        // ORIGINAL: f_x_2 ---> (?????)
  real_number f_y_1;		                    // ORIGINAL: f_y_1 ---> (?????)
  real_number f_y_2;			                  // ORIGINAL: f_y_2 ---> (?????)
  std::array<real_number, 4> p_e;           // ORIGINAL: probabilità di perdita parete
  std::array<real_number, 4> U_wx;          // ORIGINAL: Uwx (velocità pareti x)
  std::array<real_number, 4> U_wy;          // ORIGINAL: Uwx (velocità pareti y)
  real_number L_x_1;	                      // ORIGINAL: Lx1 (semi-distanza pareti x1)
  real_number L_x_2;		                    // ORIGINAL: Lx2 (   "		"	 x2)
  real_number L_y_1;			                  // ORIGINAL: Ly1 (   "		"	 y1)
  real_number L_y_2;				                // ORIGINAL: Ly2 (	 "		"	 y2)
  real_number F_x_ext;                      // ORIGINAL: Fx ---> (?????)
  // real_number F_y_ext;                   // ...

  // Computational and numerical parameters
  // int n_Iw;                              // ORIGINAL: n_Iw (numero punti per spline potenziale gas-wall)
  // int n_If;                              // ORIGINAL: n_If (numero punti per spline potenziale gas-gas)
  real_number x_min;                        // ORIGINAL: xmin (...)
  real_number x_max;			                  // ORIGINAL: xmax (...)
  real_number y_min;		                    // ORIGINAL: ymin (...)
  real_number y_max;                        // ORIGINAL: ymax (...)
  real_number x_extra;                      // ORIGINAL: xextra (troncamento campo medio)
  real_number y_extra;                      // ORIGINAL: yextra (   "		 "	  )
  int n_cells_x;                            // ORIGINAL: ncellex (numero celle in dir. x)
  int n_cells_y;                            // ORIGINAL: ncelley (numero celle in dir. y)
  int n_part;                               // ORIGINAL: Npart (numero particelle liquido)
  bool mean_vel;                            // ORIGINAL: mean_vel ---> (?????)
  bool restart;                             // ORIGINAL: restart ---> (?????)
  real_number vx_ini;                       // ORIGINAL: vx_ini (...)
  real_number vy_ini;                       // ORIGINAL: vy_ini (...)
  real_number vz_ini;                       // ORIGINAL: vz_ini (...)
  real_number t_ini;		                    // ORIGINAL: Tini	(...)
  real_number t_max;	                      // ORIGINAL: Tmax	(...)
  real_number t_im;                         // ORIGINAL: Tim ---> (?????)
  int nc_out;                               // ORIGINAL: nc_out (numero intervalli tra Tim e Tmax)
  int iter_out_ini;                         // ORIGINAL: iter_out_ini	---> (?????)
  real_number delta_t;                      // ORIGINAL: Deltat (time-step)
  int qwrite;                               // ORIGINAL: qwrite	---> (?????)
  int seed;                                 // ORIGINAL: Seme (...)
  int Nv;                                   // ORIGINAL: Nv	---> (?????)
  int routine_choice;                       // ORIGINAL: routine_choice ---> (?????)
  int split_type;                           // ORIGINAL: split_type	---> (?????)
  char c_med_comp_type;                     // ORIGINAL: c_med_comp_type (choice for build_I_matrix routine)
  bool collstat;                            // ORIGINAL: collstat	(collisions statistics off/on [0/1])
  int ndom;                                 // ORIGINAL: ndom	(number of subdomains for statistics)

public:

  ConfigurationReader(DSMC*, const DefaultString&);
  ~ConfigurationReader() = default;

  // Getter mathods

  inline int get_seed() const { return seed; }
  inline int get_L_x_1() const { return L_x_1; }
  inline int get_L_x_2() const { return L_x_2; }
  inline int get_L_y_1() const { return L_y_1; }
  inline int get_L_y_2() const { return L_y_2; }
  inline const std::array<char, 4>& get_wall_cond() const { return wall_cond; }
  inline const std::array<real_number, 4>& get_p_e() const { return p_e; }

  inline real_number get_diam_fluid() { return diam_fluid; }
  inline real_number get_diam_solid() { return diam_solid; }
  inline real_number get_mass_fluid() { return mass_fluid; }
  inline real_number get_mass_solid() { return mass_solid; }

  inline real_number get_phi11() { return phi11; }
  inline real_number get_gamma11() { return gamma11; }
  inline int get_n_cells_x() const { return n_cells_x; }
  inline int get_n_cells_y() const { return n_cells_y; }
  inline real_number get_x_min() const { return x_min; }
  inline real_number get_x_max() const { return x_max; }
  inline real_number get_y_min() const { return y_min; }
  inline real_number get_y_max() const { return y_max; }
  inline real_number get_x_extra() const { return x_extra; }
  inline real_number get_y_extra() const { return y_extra; }

  inline real_number get_eta_liq0() const { return eta_liq0; }
  inline real_number get_eta_liq1() const { return eta_liq1; }

  inline int get_n_part() const { return n_part; }

  inline int get_nc_out() const { return nc_out; }
  inline int get_restart() const { return restart; }
  inline real_number get_t_ini() const { return t_ini; }
  inline real_number get_t_max() const { return t_max; }
  inline real_number get_t_im() const { return t_im; }
  inline real_number get_delta_t() const { return delta_t; }

  inline int get_niter_thermo() const { return niter_thermo; }
  inline real_number get_T_ref() const { return T_ref; }
  inline real_number get_T_ini() const { return T_ini; }

private:

  void read_conf_file( void );

};

#endif /* EV_CONFIGURATION_HPP */
