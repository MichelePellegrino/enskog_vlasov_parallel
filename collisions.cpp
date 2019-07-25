#include "collisions.hpp"
#include "parallel_environment.hpp"
#include "topology.hpp"
#include "particles.hpp"
#include "density.hpp"
#include "grid.hpp"

CollisionHandler::CollisionHandler
(DSMC* dsmc):

  Motherbase(dsmc),

  a11( topology->get_idx_lx(par_env->get_rank()), topology->get_idx_ux(par_env->get_rank()),
    topology->get_idx_ly(par_env->get_rank()), topology->get_idx_uy(par_env->get_rank()), 0 ),

  vrmax11( topology->get_idx_lx(par_env->get_rank()), topology->get_idx_ux(par_env->get_rank()),
    topology->get_idx_ly(par_env->get_rank()), topology->get_idx_uy(par_env->get_rank()), 0 ),

  anew( topology->get_idx_lx(par_env->get_rank()), topology->get_idx_ux(par_env->get_rank()),
    topology->get_idx_ly(par_env->get_rank()), topology->get_idx_uy(par_env->get_rank()), 0 ),

  n_coll(0.0),
  n_coll_cell( topology->get_idx_lx(par_env->get_rank()), topology->get_idx_ux(par_env->get_rank()),
    topology->get_idx_ly(par_env->get_rank()), topology->get_idx_uy(par_env->get_rank()), 0 ),

  scaled_k( 0.0, 3 ),

  npart( ensemble->get_n_particles() ),

  xmin( grid->get_x_min() ),
  xmax( grid->get_x_max() ),
  ymin( grid->get_y_min() ),
  ymax( grid->get_y_max() ),

  rdx( grid->get_rdx() ),
  rdy( grid->get_rdy() ),

  sigma( species->get_diam_fluid() ),
  delta_t( times->get_delta_t() ),

  lx_sub( a11.get_lx() ),
  ux_sub( a11.get_ux()  ),
  ly_sub( a11.get_ly()  ),
  uy_sub( a11.get_uy()  )

  {

    for (int i = 0; i<N_COLLISION_SUBDOM; ++i)
      subdom_order[i] = i;

  }

void
CollisionHandler::shuffle_order
(void)
{

  if (par_env->get_rank() == MPI_MASTER)
  {
    std::random_shuffle(subdom_order.begin(), subdom_order.end());
  }

  par_env->broadcast(*subdom_order.data(), N_COLLISION_SUBDOM);
  par_env->barrier();

}

void
CollisionHandler::compute_majorants
(void)
{
  int idx_p1, ip1, jp1, idx_ck, ick, jck, ichk, jchk, jjp2, jp2;
  real_number xk, yk, xkh, ykh, chi11;
  int np2;
  int ntest = TEST_COEFF_MULT*npart;
  a11.fill(0.0);
  anew.fill(0.0);
  vrmax11.fill(0.0);
  for (int itest = 0; itest<ntest; itest++)
  {
    idx_p1 = (int)( rng->sample_uniform() * npart );
    ip1 = (int)( ( ensemble->get_xp(idx_p1) - xmin ) * rdx );
    jp1 = (int)( ( ensemble->get_yp(idx_p1) - ymin ) * rdy );
    gen_scaled_k();
    xk = ensemble->get_xp(idx_p1) - scaled_k[0];
    yk = ensemble->get_yp(idx_p1) - scaled_k[1];
    ick = (int)( ( xk - xmin ) * rdx );
    jck = (int)( ( yk - ymin ) * rdy );
    if ( ( ick >= lx_sub && ick < ux_sub ) && ( jck >= ly_sub && jck < uy_sub ) )
    {
      xkh = xk + scaled_k[0]/2.0;
      ykh = yk + scaled_k[1]/2.0;
      idx_ck = grid->lexico(ick, jck);
      if ( density->get_npc(ick, jck) >= 1 )
      {
        ichk = (int)( (xkh - xmin) * rdx );
        jchk = (int)( (ykh - ymin) * rdy );
        chi11 = correlation( density->get_aveta(ip1,jp1) );
        a11(ip1, jp1) = std::max(a11(ip1, jp1), density->get_numdens(ip1, jp1) * chi11);
        a11(ichk, jchk) = std::max(a11(ichk, jchk), density->get_numdens(ichk, jchk) * chi11);
        anew(ip1, jp1) = a11(ip1, jp1);
        np2 = density->iof(idx_ck);
        jjp2 = np2 + (int)( rng->sample_uniform() * density->get_npc(ick, jck) );
        jp2 = density->ind(jjp2);
        vr = sqrt (
            ev_utility::power<2>( ensemble->get_vx(jp2) - ensemble->get_vx(idx_p1) )
          + ev_utility::power<2>( ensemble->get_vy(jp2) - ensemble->get_vy(idx_p1) )
          + ev_utility::power<2>( ensemble->get_vz(jp2) - ensemble->get_vz(idx_p1) )
        );
        vrmax11(ip1, jp1) = std::max( vrmax11(ip1, jp1), vr );
        vrmax11(ick,jck) = std::max( vrmax11(ick,jck), vr );
      }
    }
  }
}

void
CollisionHandler::compute_collision_number
(void)
{
  real_number pfnc, cnc;
  n_coll_cell = 0;
  n_coll = 0;
  for (int i = lx_sub; i<ux_sub; ++i)
  {
    for (int j = ly_sub; j<uy_sub; ++j)
    {
      cnc = ev_const::pi2*sigma*sigma*a11(i,j)*vrmax11(i,j)*delta_t;
      n_coll_cell(i,j) = (int)cnc;
      pfnc = cnc - n_coll_cell(i,j);
      if ( rng->sample_uniform() < pfnc ) n_coll_cell(i,j)++;
    }
    n_coll = n_coll_cell.sum();
  }
}


// DEBUG e TESTING

void
CollisionHandler::print_subdomain_order
(void) const
{
  for (int i = 0; i<N_COLLISION_SUBDOM; ++i)
    std::cout << subdom_order[i] << "\t";
  std::cout << std::endl;
}
