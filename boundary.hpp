/*! \file boundary.hpp
 *  \brief Header containing the class representing boundary conditions
 */

#ifndef BOUNDARY_HPP
#define BOUNDARY_HPP

#include "motherbase.hpp"

// Include array (edges) and set (tags for b.c.)
#include <array>
#include <set>

// Redundant forward declaration (can be avoided).
class DSMC;


/*! \class Boundary
 *  \brief Class storing data regarding the position and the nature of b.c.
 *
 *  This class simply stores edges lengths and related boundary condistions.
 *  Solid walls behaviour is not included (separate class needed).
 */
class Boundary : protected Motherbase
{

private:

  /*              3 (x2)
                - - - - -
              |           |
      4 (y2)  |           |   2 (y1)
              |           |
                - - - - -
                  1 (x1)                */

  real_number Lx1, Lx2, Ly1, Ly2;         /*!< Distance from the origin of the walls (solid, holes or prescribed) */
  std::array<real_number, 4> p_escape;    /*!< Probability for Maxwellian wall with escaping molecules.           */
  std::array<char, 4> wall_condition;     /*!< Boundary conditions                                                */

  /*! \var cond_tags
   *  \brief Set containing all possible labels for b.c.
   *
   *  'i' = inlet
   *  'o' = outlet
   *  'r' = reflection
   *  'p' = periodic
   *  'e' = Maxwell (out)
   *  'm' = Maxwell (in/out)
   *  'w' = wall
   *  'h' = hole
   */
  static const std::set<char> cond_tags;      /*!< All possible b.c. tags           */
  static const std::set<char> log_cond_tags;  /*!< All possible logical b.c. tags   */
  static const std::set<char> pys_cond_tags;  /*!< All possible physical b.c. tags  */

public:

  Boundary(DSMC*);
  ~Boundary() = default;

  // GETTERS
  inline real_number get_Lx1(void) const { return Lx1; }
  inline real_number get_Ly1(void) { return Ly1; }
  inline real_number get_Lx2(void) const { return Lx2; }
  inline real_number get_Ly2(void) { return Ly2; }
  inline const std::array<real_number, 4>& get_p_escape(void) const { return p_escape; }
  inline const std::array<char, 4>& get_wall_condition() const { return wall_condition; }
  inline const std::set<char>& get_ct(void) const { return cond_tags; }
  inline const std::set<char>& get_lct(void) const { return log_cond_tags; }
  inline const std::set<char>& get_pct(void) const { return pys_cond_tags; }

};

#endif /* BOUNDARY_HPP */
