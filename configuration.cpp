#include "configuration.hpp"
#include "parallel_environment.hpp"

#include <fstream>
#include <sstream>
#include <iostream>

ConfigurationReader::ConfigurationReader
(DSMC* dsmc, const DefaultString& file_name):
Motherbase(dsmc), conf_file_name(file_name)
{
  par_env->barrier();
  read_conf_file();
  setup_initial_configuration();
}

void
ConfigurationReader::read_conf_file
(void)
{

  if( par_env->is_root() )
  {

  std::ifstream fs;
  std::stringstream ss;
  ss << std::scientific;
  fs.open(conf_file_name);
  DefaultString line_buffer;

  std::cout << "### READING CONFIGURATION PARAMETERS ###" << std::endl;
  std::cout << " >> Configuration file: " << conf_file_name << std::endl;

  getline (fs, line_buffer);                // ##  **DEFINIZIONE DELLA PROVA
  getline (fs, line_buffer);                // #
  getline (fs, line_buffer);                // # NUMERO SEQUENZIALE DELLA PROVA ======>.... droppa2

  ss << line_buffer;
  ss.seekg(45);
  ss >> numtest;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // #
  getline (fs, line_buffer);                // # ** PROPRIETA' MOLECOLARI
  getline (fs, line_buffer);                // #
  getline (fs, line_buffer);                // # MASSA MOLECOLARE LIQUIDO ============>.... 0.1000e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> mass_fluid;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # DIAMETRO MOLECOLARE LIQUIDO =========>.... 0.1000e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> diam_fluid;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Tipo di potenziale gas-gas          =>.... p

  ss << line_buffer;
  ss.seekg(45);
  ss >> pot_gas;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Campo medio gas-gas yes or no       =>.... y

  ss << line_buffer;
  ss.seekg(45);
  ss >> mean_f_gg;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # PARAMETRO POTENZIALE PHI11 ==========>.... 0.1000e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> phi11;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # ESPONENTE POTENZIALE PHI11 ==========>.... 0.6000e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> gamma11;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Numero punti per spline pot gas-gas =>....     1200
  getline (fs, line_buffer);                // # MASSA MOLECOLARE SOLIDO =============>.... 0.2000e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> mass_solid;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # DIAMETRO MOLECOLARE SOLIDO ==========>.... 0.1000e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> diam_solid;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Tipo di potenziale gas-solido       =>.... e

  ss << line_buffer;
  ss.seekg(45);
  ss >> pot_wall;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Campo medio gas-wall yes or no      =>.... y

  ss << line_buffer;
  ss.seekg(45);
  ss >> mean_f_gw;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # PARAMETRO POTENZIALE PHI12 ==========>.... 0.1500e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> phi12;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # ESPONENTE POTENZIALE PHI12 ==========>.... 0.5405e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> gamma12;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Numero punti per spline pot gas-sol =>....     1200
  getline (fs, line_buffer);                // #
  getline (fs, line_buffer);                // # ** STATO FISICO DELLE PARETI          :
  getline (fs, line_buffer);                // #
  getline (fs, line_buffer);                // # Temp. iniziale liquido ==============>.... 0.5000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> T_ini;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Temp. iniziale liquido ==============>.... 0.5000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> T_ini1;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # DENSITA' DELLA PARETE per campo medio>.... 0.7000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> eta_w0;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # DENSITA' DELLA PARETE per correlazio >.... 0.7000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> eta_w1;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # DENSITA' DEL GAS     ================>.... 0.6490e-02

  ss << line_buffer;
  ss.seekg(45);
  ss >> eta_liq0;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # DENSITA' DEL LIQUIDO ================>.... 0.3763e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> eta_liq1;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Presenza di strato liquido iniziale==>....       13

  ss << line_buffer;
  ss.seekg(45);
  ss >> liq_interf;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # spessore strato liquido iniziale x ==>.... 0.1800e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> x_liq_interf;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # spessore strato liquido iniziale y ==>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> y_liq_interf;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Raggio goccia iniziale liquida     ==>.... 0.3000e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> r_liq_interf;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Offset goccia liquida iniziale     ==>.... 0.1400e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> liq_drop_offset;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Number of iterations for thermostat =>....      50

  ss << line_buffer;
  ss.seekg(45);
  ss >> niter_thermo;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Reference temperature for thermostat=>.... 0.5000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> T_ref;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Keep baricentre fixed (1) or not (0)=>....        0

  ss << line_buffer;
  ss.seekg(45);
  ss >> fix_baricentre;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Natura del Bordo x1, parete o libero =>... h

  ss << line_buffer;
  ss.seekg(45);
  ss >> wall_cond[0];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Natura del Bordo x2, parete o libero =>... h

  ss << line_buffer;
  ss.seekg(45);
  ss >> wall_cond[2];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Natura del Bordo y1, parete o libero =>... w

  ss << line_buffer;
  ss.seekg(45);
  ss >> wall_cond[1];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Natura del Bordo y2, parete o libero =>... m

  ss << line_buffer;
  ss.seekg(45);
  ss >> wall_cond[3];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Sets(1),ignores(0) periodicity in mean f x        1

  ss << line_buffer;
  ss.seekg(45);
  ss >> set_px;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Sets(1),ignores(0) periodicity in mean f y        0

  ss << line_buffer;
  ss.seekg(45);
  ss >> set_py;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Sets(1),ignores(0) periodicity in avdens.x        1

  ss << line_buffer;
  ss.seekg(45);
  ss >> set_eta_px;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Sets(1),ignores(0) periodicity in avdens.y        0

  ss << line_buffer;
  ss.seekg(45);
  ss >> set_eta_py;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # TEMPERATURA PARETE X1 ===============>.... 0.5000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> T_w[0];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # TEMPERATURA PARETE X2 ===============>.... 0.5000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> T_w[2];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # TEMPERATURA PARETE Y1 ===============>.... 0.5000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> T_w[1];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # TEMPERATURA PARETE Y2 ===============>.... 0.5000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> T_w[3];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Fattore flusso parete x1  ===========>.... 0.1000e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> f_x_1;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Fattore flusso parete x2  ===========>.... 0.1000e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> f_x_2;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Fattore flusso parete y1  ===========>.... 0.1000e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> f_y_1;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Fattore flusso parete y2  ===========>.... 0.1000e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> f_y_2;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Probabilita perdita parete x1 =======>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> p_e[0];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Probabilita perdita parete x2 =======>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> p_e[2];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Probabilita perdita parete y1 =======>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> p_e[1];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Probabilita perdita parete y2 =======>.... 0.5000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> p_e[3];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # VELOCITA' PARETI  1 X ===============>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> U_wx[0];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # VELOCITA' PARETI  1 Y ===============>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> U_wy[0];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # VELOCITA' PARETI  2 X ===============>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> U_wx[1];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # VELOCITA' PARETI  2 Y ===============>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> U_wy[1];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # VELOCITA' PARETI  3 X ===============>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> U_wx[2];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # VELOCITA' PARETI  3 Y ===============>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> U_wy[2];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # VELOCITA' PARETI  4 X ===============>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> U_wx[3];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # VELOCITA' PARETI  4 Y ===============>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> U_wy[3];
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # SEMI-DISTANZA PARETI X1 =============>.... 0.2500e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> L_x_1;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # SEMI-DISTANZA PARETI X2 =============>.... 0.2500e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> L_x_2;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # SEMI-DISTANZA PARETI Y1 =============>.... 0.2500e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> L_y_1;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # SEMI-DISTANZA PARETI Y2 =============>.... 0.2900e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> L_y_2;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # INTENSITA' CAMPO ESTERNO(POISEUILLE)=>.... 0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> F_x_ext;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // #
  getline (fs, line_buffer);                // # ** PARAMETRI COMPUTAZIONALI
  getline (fs, line_buffer);                // #
  getline (fs, line_buffer);                // # XMIN ================================>.... 0.3500e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> x_min;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # XMAX ================================>.... 0.3500e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> x_max;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # YMIN ================================>.... 0.2900e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> y_min;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # YMAX ================================>.... 0.2900e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> y_max;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # XEXTRA (TRONCAMENTO CAMPO MEDIO)=====>     0.4000e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> x_extra;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # YEXTRA (TRONCAMENTO CAMPO MEDIO)=====>     0.4000e01

  ss << line_buffer;
  ss.seekg(45);
  ss >> y_extra;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # NUMERO CELLE DIR. X =================>....      700

  ss << line_buffer;
  ss.seekg(45);
  ss >> n_cells_x;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # NUMERO CELLE DIR. Y =================>....      580

  ss << line_buffer;
  ss.seekg(45);
  ss >> n_cells_y;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # NUMERO PARTICELLE LIQUIDO ===========>....  3000000

  ss << line_buffer;
  ss.seekg(45);
  ss >> n_part;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # mean_vel=0 (0) or /=0 (1)  ==========>            1

  ss << line_buffer;
  ss.seekg(45);
  ss >> mean_vel;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Restart 0 no 1 yes 2 yes mirror  ====>            0

  ss << line_buffer;
  ss.seekg(45);
  ss >> restart;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # vx_ini ==============================>     0.0000e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> vx_ini;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # vy_ini ==============================>     -0.600e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> vy_ini;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # vz_ini ==============================>     -0.600e00

  ss << line_buffer;
  ss.seekg(45);
  ss >> vz_ini;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # TINI ================================>     0.0000e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> t_ini;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # TMAX ================================>     1.0000e03

  ss << line_buffer;
  ss.seekg(45);
  ss >> t_max;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # TIM =================================>     2.0000e02

  ss << line_buffer;
  ss.seekg(45);
  ss >> t_im;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Numero intervalli tra tim e tmax ====>           1

  ss << line_buffer;
  ss.seekg(45);
  ss >> nc_out;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Partenza del contatore per intervalli>           0

  ss << line_buffer;
  ss.seekg(45);
  ss >> iter_out_ini;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # DELTAT ==============================>     0.5000e-02

  ss << line_buffer;
  ss.seekg(45);
  ss >> delta_t;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Qwrite                                            0

  ss << line_buffer;
  ss.seekg(45);
  ss >> qwrite;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # SEME  ===============================>           76

  ss << line_buffer;
  ss.seekg(45);
  ss >> seed;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # NUMERO NODI VELOCITA' FUNZ. DISTR. ==>          100

  ss << line_buffer;
  ss.seekg(45);
  ss >> Nv;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Routine to compute collisions choice=>     2

  ss << line_buffer;
  ss.seekg(45);
  ss >> routine_choice;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Split type for parallel collisions  =>     1

  ss << line_buffer;
  ss.seekg(45);
  ss >> split_type;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Choice for build_I_matrix routine   =>     o

  ss << line_buffer;
  ss.seekg(45);
  ss >> c_med_comp_type;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Collision statistics off (0) on (1) =>     0

  ss << line_buffer;
  ss.seekg(45);
  ss >> collstat;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Number of subdomains for statistics =>     6

  ss << line_buffer;
  ss.seekg(45);
  ss >> ndom;
  ss.clear(); ss.str(std::string());

  getline (fs, line_buffer);                // # Number of iterations for sampling

  ss << line_buffer;
  ss.seekg(45);
  ss >> niter_sampling;
  ss.clear(); ss.str(std::string());

  fs.close();

  }

  par_env->broadcast(mean_f_gg);

  par_env->broadcast(seed);
  par_env->broadcast(L_x_1);
  par_env->broadcast(L_x_2);
  par_env->broadcast(L_y_1);
  par_env->broadcast(L_y_2);

  for (int i = 0; i<4; ++i)
  {
    par_env->broadcast(p_e[i]);
    par_env->broadcast(wall_cond[i]);
  }

  par_env->broadcast(diam_fluid);
  par_env->broadcast(diam_solid);
  par_env->broadcast(mass_fluid);
  par_env->broadcast(mass_solid);

  par_env->broadcast(phi11);
  par_env->broadcast(gamma11);

  par_env->broadcast(n_cells_x);
  par_env->broadcast(n_cells_y);
  par_env->broadcast(x_min);
  par_env->broadcast(x_max);
  par_env->broadcast(y_min);
  par_env->broadcast(y_max);
  par_env->broadcast(x_extra);
  par_env->broadcast(y_extra);

  par_env->broadcast(liq_interf);
  par_env->broadcast(x_liq_interf);
  par_env->broadcast(y_liq_interf);

  par_env->broadcast(eta_liq0);
  par_env->broadcast(eta_liq1);

  if ( par_env->is_root() )
    std::cout << " >> npc = " << n_part/(double)(n_cells_x*n_cells_y) << std::endl;

  par_env->broadcast(n_part);
  npart_global = n_part;
  int rem = n_part % par_env->get_size();
  n_part = n_part / par_env->get_size();
  if ( par_env->get_rank()<rem ) n_part++;

  par_env->broadcast(nc_out);
  par_env->broadcast(restart);
  par_env->broadcast(t_ini);
  par_env->broadcast(t_max);
  par_env->broadcast(t_im);
  par_env->broadcast(delta_t);

  par_env->broadcast(niter_thermo);
  par_env->broadcast(niter_sampling);
  par_env->broadcast(T_ref);
  par_env->broadcast(T_ini);

  dx = (double)(x_max+x_min) / (double)n_cells_x;
  dy = (double)(y_max+y_min) / (double)n_cells_y;

  par_env->broadcast(vx_ini);
  par_env->broadcast(vy_ini);
  par_env->broadcast(vz_ini);

}

void
ConfigurationReader::setup_initial_configuration
(void)
{

  int rem;

  if ( par_env->is_root() )
    std::cout << "### SETTING-UP INITIAL CONFIGURATION ###" << std::endl;

  real_number tmp0, tmp1;
  switch( liq_interf )
  {
    case 0:     /* Uniform phase */
      if ( par_env->is_root() )
        std::cout << "### CASE: uniform single phase ###" << std::endl;
      // real_number channel_area = (x_lim_p-x_lim_m) * (y_lim_p-y_lim_m);
      channel_area = (x_max+x_min) * (y_max+y_min);
      homogeneous_density = 6.0 * eta_liq0 / (ev_const::pi * diam_fluid*diam_fluid*diam_fluid);
      channel_volume = npart_global / homogeneous_density;
      channel_section = channel_volume / channel_area;
      cell_volume = dx*dy*channel_section;
      // Distributing
      npart0 = n_part;  npart1 = 0;
      break;
    case 5:     /* Initial liquid layer in the middle, parallel to the x axis */
      if ( par_env->is_root() )
       std::cout << "### CASE: initial liquid layer in the middle (parallel to the x axis) ###" << std::endl;
      channel_area0 = (x_max+x_min) * (y_max+y_min-y_liq_interf);
      channel_area1 = (x_max+x_min) * y_liq_interf;
      homogeneous_density0 = 6.0 * eta_liq0 / (ev_const::pi * diam_fluid*diam_fluid*diam_fluid);
      homogeneous_density1 = 6.0 * eta_liq1 / (ev_const::pi * diam_fluid*diam_fluid*diam_fluid);
      tmp0 = 1.0/(homogeneous_density0*channel_area0);
      tmp1 = 1.0/(homogeneous_density1*channel_area1);
      npart1 = (int)(npart_global*tmp0/(tmp0+tmp1));
      npart0 = npart_global - npart1;
      channel_volume0 = npart0 / homogeneous_density0;
      channel_volume1 = npart1 / homogeneous_density1;
      channel_section = (channel_volume0+channel_volume1) / (channel_area0+channel_area1);
      cell_volume = dx*dy*channel_section;
      // Distributing
      rem = npart0 % par_env->get_size(); npart0 = npart0 / par_env->get_size();
      if ( par_env->get_rank()<rem ) npart0++;
      rem = npart1 % par_env->get_size(); npart1 = npart1 / par_env->get_size();
      if ( par_env->get_rank()<rem ) npart1++;
      break;
    case 6:     /* Initial liquid layer in the middle, parallel to the y axis */
      if ( par_env->is_root() )
        std::cout << "### CASE: initial liquid layer in the middle (parallel to the y axis) ###" << std::endl;
      channel_area0 = (x_max+x_min-x_liq_interf) * (y_max+y_min);
      channel_area1 = x_liq_interf * (y_max+y_min);
      homogeneous_density0 = 6.0 * eta_liq0 / (ev_const::pi * diam_fluid*diam_fluid*diam_fluid);
      homogeneous_density1 = 6.0 * eta_liq1 / (ev_const::pi * diam_fluid*diam_fluid*diam_fluid);
      tmp0 = 1.0/(homogeneous_density0*channel_area0);
      tmp1 = 1.0/(homogeneous_density1*channel_area1);
      npart1 = (int)(npart_global*tmp0/(tmp0+tmp1));
      npart0 = npart_global - npart1;
      channel_volume0 = npart0 / homogeneous_density0;
      channel_volume1 = npart1 / homogeneous_density1;
      channel_section = (channel_volume0+channel_volume1) / (channel_area0+channel_area1);
      cell_volume = dx*dy*channel_section;
      // Distributing
      rem = npart0 % par_env->get_size(); npart0 = npart0 / par_env->get_size();
      if ( par_env->get_rank()<rem ) npart0++;
      rem = npart1 % par_env->get_size(); npart1 = npart1 / par_env->get_size();
      if ( par_env->get_rank()<rem ) npart1++;
      break;
    default:
      if ( par_env->is_root() )
        std::cerr << "Unrecognized configuration" << std::endl;
  }
}
