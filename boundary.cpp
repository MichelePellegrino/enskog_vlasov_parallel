#include "boundary.hpp"
#include "configuration.hpp"

const std::set<char> Boundary::cond_tags = {'i', 'o', 'r', 'p', 'e', 'm', 'w', 'h'};
const std::set<char> Boundary::log_cond_tags = {'i', 'o', 'r', 'p', 'e', 'm'};
const std::set<char> Boundary::pys_cond_tags = {'w', 'h'};

Boundary::Boundary
(DSMC* dsmc):
  Motherbase(dsmc),
  Lx1( conf->get_L_x_1() ),
  Lx2( conf->get_L_x_2() ),
  Ly1( conf->get_L_y_1() ),
  Ly2( conf->get_L_y_2() ),
  p_escape( conf->get_p_e() ),
  wall_condition( conf->get_wall_cond() )
  {
    for (auto it = wall_condition.cbegin(); it!=wall_condition.cend(); ++it)
    {
      assert(cond_tags.find(*it)!=cond_tags.end() && "Invalid boundary condition");
    }
    assert (
      !( ( wall_condition[0]=='p' && wall_condition[2]!='p' ) || ( wall_condition[0]!='p' && wall_condition[2]=='p' ) )
      && "Invalid periodic conditions on x-wall" );
    assert (
      !( ( wall_condition[1]=='p' && wall_condition[3]!='p' ) || ( wall_condition[1]!='p' && wall_condition[3]=='p' ) )
      && "Invalid periodic conditions on y-wall" );
  }
