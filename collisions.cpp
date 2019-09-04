#include "collisions.hpp"
#include "parallel_environment.hpp"
#include "topology.hpp"
#include "particles.hpp"
#include "density.hpp"
#include "grid.hpp"

#include <set>

CollisionHandler::CollisionHandler
(DSMC* dsmc):

  Motherbase(dsmc),

  size(par_env->get_size()),
  rank(par_env->get_rank()),

  n_fake_store(),
  n_real_store(),
  n_total_store(),
  n_out_store(),

  a11( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),

  vrmax11( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),

  anew( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),

  vrmaxnew( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),

  n_coll(0.0),
  n_coll_cell( topology->get_idx_lx(rank), topology->get_idx_ux(rank),
    topology->get_idx_ly(rank), topology->get_idx_uy(rank), 0.0 ),

  scaled_k( 0.0, 3 ),
  rel_vel( 0.0, 3 ),
  delta( 0.0, 3 ),

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
  ux_sub( a11.get_ux() ),
  ly_sub( a11.get_ly() ),
  uy_sub( a11.get_uy() )

  {

    commit_velocity_type();
    for (int i = 0; i<N_COLLISION_SUBDOM; ++i)
      subdom_order[i] = i;
    init_cell_ind_map();

  }

void
CollisionHandler::shuffle_order
(void)
{

  if (rank == MPI_MASTER)
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
  vrmaxnew.fill(0.0);
  // DEBUG
  // # # # # #
  // real_number correlation_sum(0.0);
  // # # # # #
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
        // DEBUG
        // # # # # #
        // correlation_sum += chi11;
        // # # # # #
        a11(ip1, jp1) = std::max( a11(ip1, jp1), density->get_numdens(ip1, jp1) * chi11 );
        a11(ick, jck) = std::max( a11(ick, jck), density->get_numdens(ick, jck) * chi11 );
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
        vrmaxnew(ip1, jp1) = vrmax11(ip1, jp1);
        vrmaxnew(ick, jck) = vrmax11(ick, jck);
      }
    }
  }
  // DEBUG
  // # # # # #
  // print_reduced_constants();
  // std::cout << "rank " << rank << " : correlation_sum = " << correlation_sum << std::endl;
  // # # # # #
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
  }
  n_coll = n_coll_cell.sum();
}

void
CollisionHandler::init_cell_ind_map
(void)
{
  int nc;
  for (int q = 0; q<N_COLLISION_SUBDOM; ++q)
  {
    nc = topology->get_n_cells_quarter(q);
    n_cells_x[q] = topology->get_quarter_ux(q) - topology->get_quarter_lx(q);
    low_x_quad[q] = topology->get_quarter_lx(q);
    low_y_quad[q] = topology->get_quarter_ly(q);
    cells_ind[q].resize(nc);
    for (int i = 0; i<nc; ++i)
      cells_ind[q][i] = i;
  }
}

void
CollisionHandler::setup_cell_ind_map
(void)
{
  int nc;
  for (int q = 0; q<N_COLLISION_SUBDOM; ++q)
  {
    nc = cells_ind[q].size();
    for (int i = 0; i<nc; ++i)
      cells_ind[q][i] = i;
  }
}

void
CollisionHandler::perform_collisions
(void)
{

  int nc;
  int idx, idx_cell1, i_cell1, j_cell1, idx_p1;
  int idx_cell2, idx_hcell2, i_cell2, j_cell2, i_hcell2, j_hcell2, idx_p2;
  int i_cell1_quad, j_cell1_quad;
  int rank_coll_part;
  real_number xk, yk, xkh, ykh, aa, fk;
  CollisionCounter parallel_collision_counter;
  // Set-up
  setup_cell_ind_map();
  compute_collision_number();
  shuffle_order();
  // Reset number of collisions
  int out_bound = 0;
  n_fake = 0; n_real = 0; n_total = 0;
  // DEBUG
  // # # # # #
  print_est_collisions();
  int parallel_counter = 0, communication_counter = 0;
  int null_density_counter = 0, scalar_prod_counter = 0;
  // # # # # #
  for (int q : subdom_order)
  {
    reset_parallel_buffers();
    nc = topology->get_n_cells_quarter(q);
    while ( nc > 0 )
    {
      // (1) Select a cell at random with with equiprobability
      idx = (int)( rng->sample_uniform() * nc );
      idx_cell1 = cells_ind[q][idx];
      cells_ind[q][idx] = cells_ind[q][nc-1];
      nc -= 1;
      i_cell1_quad = lexico_inv_quad(idx_cell1, q).first;
      j_cell1_quad = lexico_inv_quad(idx_cell1, q).second;
      i_cell1 = map_to_domain(i_cell1_quad, j_cell1_quad, q).first;
      j_cell1 = map_to_domain(i_cell1_quad, j_cell1_quad, q).second;
      // (2) Select a particle belonging to cell at random with equiprobability
      for (int i1 = 0; i1 < n_coll_cell(i_cell1, j_cell1); i1++)
      {
        idx_p1 = density->iof(idx_cell1) + (int)( rng->sample_uniform() *
          density->get_npc(i_cell1, j_cell1) );
        idx_p1 = density->ind(idx_p1);
        // (3) Select, at random, the unit k-vector
        gen_scaled_k();
        xk = ensemble->get_xp(idx_p1) - scaled_k[0];
        yk = ensemble->get_yp(idx_p1) - scaled_k[1];
        xkh = xk + scaled_k[0]/2.0;
        ykh = yk + scaled_k[1]/2.0;
        // BOUNDARY CONDITIONS: we suppose periodic b.c.
        if ( xk <= xmin || xk >= xmax )
        {
          xk = xk - round( (xk-0.5*(xmax+xmin))/(xmax-xmin) ) * (xmax-xmin);
          xkh = xkh - round( (xkh-0.5*(xmax+xmin))/(xmax-xmin) ) * (xmax-xmin);
        }
        if ( yk <= ymin || yk >= ymax )
        {
          yk = yk - round( (yk-0.5*(ymax+ymin))/(ymax-ymin) ) * (ymax-ymin);
          ykh = ykh - round( (ykh-0.5*(ymax+ymin))/(ymax-ymin) ) * (ymax-ymin);
        }
        i_cell2 = (int)( (xk-xmin)*rdx );
        j_cell2 = (int)( (yk-ymin)*rdy );
        idx_cell2 = grid->lexico(i_cell2, j_cell2);
        // For DEBUG reasons: it doesn't consider minimum-image convention:
        // # # # # #
        // i_hcell2 = (int)( (xkh-xmin)*rdx );
        // j_hcell2 = (int)( (ykh-ymin)*rdy );
        i_hcell2 = ( i_cell2 + i_cell1 ) / 2;
        j_hcell2 = ( j_cell2 + j_cell1 ) / 2;
        // # # # # #
        idx_hcell2 = grid->lexico(i_hcell2, j_hcell2);
        rank_coll_part = topology->get_topology_map(i_cell2, j_cell2);
        if ( rank_coll_part == rank )
        {
          // DEBUG
          // # # # # #
          parallel_counter++;
          // # # # # #
          // 'SERIAL' PART
          if( ( density->get_npc(i_cell2, j_cell2)>0 ) )
          {
            // DEBUG
            // # # # # #
            null_density_counter++;
            // # # # # #
            // (4) Select at random with equiprobability a particle in the cell where the k-vector points
            idx = density->iof(idx_cell2) + (int)( rng->sample_uniform() * density->get_npc(i_cell2, j_cell2) );
            idx_p2 = density->ind(idx);
            rel_vel[0] = ensemble->get_vx(idx_p2) - ensemble->get_vx(idx_p1);
            rel_vel[1] = ensemble->get_vy(idx_p2) - ensemble->get_vy(idx_p1);
            rel_vel[2] = ensemble->get_vz(idx_p2) - ensemble->get_vz(idx_p1);
            vr = sqrt( rel_vel[0]*rel_vel[0]+rel_vel[1]*rel_vel[1]+rel_vel[2]*rel_vel[2] );
            vrmaxnew(i_cell1, j_cell1) = std::max( vrmaxnew(i_cell1, j_cell1), vr );
            vrmaxnew(i_cell2, j_cell2) = std::max( vrmaxnew(i_cell2, j_cell2), vr );
            scalar_prod = rel_vel[0]*scaled_k[0] + rel_vel[1]*scaled_k[1] + rel_vel[2]*scaled_k[2];
            scalar_prod *= sigma;
            // DEBUG -> PROBLEM TO BE SOLVED !!!!!
            // # # # # #
            if ( topology->get_topology_map(i_hcell2, j_hcell2) != rank )
            {
              std::cout << "(" << i_cell1 << "," << j_cell1 << ")" << std::endl;
              std::cout << "(" << i_hcell2 << "," << j_hcell2 << ")" << std::endl;
              std::cout << "(" << i_cell2 << "," << j_cell2 << ")" << std::endl;
            }
            // # # # # #
            aa = density->get_numdens(i_cell2, j_cell2) * correlation(
              density->get_aveta(grid->lexico_inv(idx_hcell2).first, grid->lexico_inv(idx_hcell2).second) );
            anew(i_cell1, j_cell1) = std::max( anew(i_cell1, j_cell1), aa );
            anew(i_cell2, j_cell2) = std::max( anew(i_cell2, j_cell2),
              density->get_numdens(i_cell1, j_cell1) * aa / density->get_numdens(i_cell2, j_cell2) );
            if( a11(i_cell2, j_cell2)==0.0 )
              a11(i_cell2, j_cell2) = anew(i_cell2, j_cell2);
            if( scalar_prod > 0.0 )
            {
              // DEBUG
              // # # # # #
              scalar_prod_counter++;
              // # # # # #
              fk = scalar_prod * aa / ( a11(i_cell1, j_cell1) * vrmax11(i_cell1, j_cell1) );
              if (fk > 1.0)
                out_bound++;
              if ( rng->sample_uniform() < fk )
              {
                n_real++;
                delta = scaled_k*scalar_prod*species->get_diam_fluid();
                ensemble->get_vx(idx_p1) += delta[0];
                ensemble->get_vy(idx_p1) += delta[1];
                ensemble->get_vz(idx_p1) += delta[2];
                ensemble->get_vx(idx_p2) -= delta[0];
                ensemble->get_vy(idx_p2) -= delta[1];
                ensemble->get_vz(idx_p2) -= delta[2];
              }
              else
              {
                n_fake++;
              }
              n_total++;
            }
          }
        }
        else
        {
          // PARALLEL COMMUNICATION - PREPARATION
          // DEBUG
          // # # # # #
          communication_counter++;
          // # # # # #
          n_parallel_collisions_send[rank_coll_part]++;
          cells1_process_map[rank_coll_part].push_back(grid->lexico(i_cell1, j_cell1));
          particle1_process_map[rank_coll_part].push_back(idx_p1);
          cells2_process_map_send[rank_coll_part].push_back(idx_cell2);
          hcells2_process_map_send[rank_coll_part].push_back(idx_hcell2);
          scaled_k_process_map[rank_coll_part].emplace_back(scaled_k);
        }
      }
    }
    parallel_collision_counter = perform_parallel_collisions();
    n_total += parallel_collision_counter.total;
    n_real += parallel_collision_counter.real;
    n_fake += parallel_collision_counter.fake;
    out_bound += parallel_collision_counter.out_bound;
    par_env->barrier();
  }

  n_fake_store.push_back(n_fake);
  n_real_store.push_back(n_real);
  n_total_store.push_back(n_total);
  n_out_store.push_back(out_bound);

  // Display statistics
  for (int r = 0; r<size; ++r)
  {
    if (rank==r)
    {
      std::cout << "collisions performed by rank " << r << std::endl;
      // DEBUG
      // # # # # #
      std::cout << "local / communication : " << parallel_counter << " / " << communication_counter << std::endl;
      // std::cout << "density counter = " << null_density_counter << ";\t scalar counter = " << scalar_prod_counter << std::endl;
      // # # # # #
      std::cout << "total = " << n_total << "\t real = " << n_real << "\t fake = " << n_fake << "\t out-range = " << out_bound << std::endl;
    }
    par_env->barrier();
  }

  update_majorants();

}

CollisionCounter
CollisionHandler::perform_parallel_collisions
(void)
{
  int idx_cell2, i_cell2, j_cell2, i_hcell2, j_hcell2, idx, idx_p2, i_cell1, j_cell1;
  int numdens_cell2, aa, fk;
  CollisionCounter collision_counter(0, 0, 0, 0);
  for (int rs = 0; rs<size; ++rs)
  {
    if ( rs == rank )
    {
      // Send buffer with cell-2 index
      for ( int rr = 0; rr<size; ++rr )
      {
        if ( rr!=rs )
        {
          MPI_Send(&n_parallel_collisions_send[rr], 1, MPI_INT, rr, 0, MPI_COMM_WORLD);
          MPI_Send(&cells2_process_map_send[rr][0], n_parallel_collisions_send[rr], MPI_INT, rr, 0, MPI_COMM_WORLD);
          MPI_Send(&hcells2_process_map_send[rr][0], n_parallel_collisions_send[rr], MPI_INT, rr, 0, MPI_COMM_WORLD);
        }
      }
    }
    else
    {
      // Receive buffer with cell-2 index
      MPI_Recv(&n_parallel_collisions_recv[rs], 1, MPI_INT, rs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      cells2_process_map_recv[rs].resize(n_parallel_collisions_recv[rs]);
      MPI_Recv(&cells2_process_map_recv[rs][0], n_parallel_collisions_recv[rs], MPI_INT, rs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      hcells2_process_map_recv[rs].resize(n_parallel_collisions_recv[rs]);
      MPI_Recv(&hcells2_process_map_recv[rs][0], n_parallel_collisions_recv[rs], MPI_INT, rs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // Select a particle at random, for each cell
      for (int i = 0; i<n_parallel_collisions_recv[rs]; ++i)
      {
        idx_cell2 = cells2_process_map_recv[rs][i];
        i_cell2 = grid->lexico_inv(idx_cell2).first;
        j_cell2 = grid->lexico_inv(idx_cell2).second;
        i_hcell2 = grid->lexico_inv(hcells2_process_map_recv[rs][i]).first;
        j_hcell2 = grid->lexico_inv(hcells2_process_map_recv[rs][i]).second;
        numdens_cell2 = density->get_numdens(i_cell2, j_cell2);
        numdens_process_map_send[rs].push_back( numdens_cell2 );
        if ( topology->get_topology_map(i_hcell2, j_hcell2) == rank)
          aveta_mid_process_map_send[rs].push_back( density->get_aveta(i_hcell2, j_hcell2) );
        else
          aveta_mid_process_map_send[rs].push_back(0.0);
        if ( numdens_cell2>0 )
        {
          idx = density->iof(idx_cell2) + (int)( rng->sample_uniform() * density->get_npc(i_cell2, j_cell2) );
          idx_p2 = density->ind(idx);
          // Store its velocity
          velocity_process_map_send[rs].emplace_back(
            ensemble->get_vx(idx_p2),
            ensemble->get_vy(idx_p2),
            ensemble->get_vz(idx_p2) );
        }
        else
        {
          velocity_process_map_send[rs].emplace_back(0.0, 0.0, 0.0);
        }
      }
    }
  }
  for (int rs = 0; rs<size; ++rs)
  {
    if ( rs == rank )
    {
      // Send buffer with velocities
      for ( int rr = 0; rr<size; ++rr )
      {
        if ( rr!=rs )
        {
          MPI_Send(&velocity_process_map_send[rr][0], n_parallel_collisions_recv[rr], MPI_VELOCITY_TYPE, rr, 0, MPI_COMM_WORLD);
          MPI_Send(&numdens_process_map_send[rr][0], n_parallel_collisions_recv[rr], MPI_INT, rr, 0, MPI_COMM_WORLD);
          MPI_Send(&aveta_mid_process_map_send[rr][0], n_parallel_collisions_recv[rr], MPI_DOUBLE, rr, 0, MPI_COMM_WORLD);
        }
      }
    }
    else
    {
      velocity_process_map_recv[rs].resize(n_parallel_collisions_send[rs]);
      MPI_Recv(&velocity_process_map_recv[rs][0], n_parallel_collisions_send[rs], MPI_VELOCITY_TYPE, rs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      numdens_process_map_recv[rs].resize(n_parallel_collisions_send[rs]);
      MPI_Recv(&numdens_process_map_recv[rs][0], n_parallel_collisions_send[rs], MPI_INT, rs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      aveta_mid_process_map_recv[rs].resize(n_parallel_collisions_send[rs]);
      MPI_Recv(&aveta_mid_process_map_recv[rs][0], n_parallel_collisions_send[rs], MPI_DOUBLE, rs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      // Compute relative speed and scalar product with k
      for ( int i = 0; i<n_parallel_collisions_send[rs]; ++i  )
      {
        rel_vel[0] = ensemble->get_vx( particle1_process_map[rs][i] ) - velocity_process_map_recv[rs][i].vx;
        rel_vel[1] = ensemble->get_vy( particle1_process_map[rs][i] ) - velocity_process_map_recv[rs][i].vy;
        rel_vel[2] = ensemble->get_vz( particle1_process_map[rs][i] ) - velocity_process_map_recv[rs][i].vz;
        scalar_prod = rel_vel[0]*scaled_k_process_map[rs][i][0]
                    + rel_vel[1]*scaled_k_process_map[rs][i][1]
                    + rel_vel[2]*scaled_k_process_map[rs][i][2];
        i_hcell2 = grid->lexico_inv(hcells2_process_map_send[rs][i]).first;
        j_hcell2 = grid->lexico_inv(hcells2_process_map_send[rs][i]).second;
        if ( topology->get_topology_map(i_hcell2, j_hcell2) == rank)
          aa = numdens_process_map_recv[rs][i] * correlation ( density->get_aveta(i_hcell2, j_hcell2) );
        else
          aa = numdens_process_map_recv[rs][i] * correlation ( aveta_mid_process_map_recv[rs][i] );
        i_cell1 = grid->lexico_inv( cells1_process_map[rs][i] ).first;
        j_cell1 = grid->lexico_inv( cells1_process_map[rs][i] ).second;
        fk = scalar_prod * aa / ( a11(i_cell1, j_cell1) * vrmax11(i_cell1, j_cell1) );

        perform_collision_send[rs].push_back(
          (scalar_prod > 0) * (numdens_process_map_recv[rs][i] > 0) * (rng->sample_uniform() < fk) );
        delta = scaled_k_process_map[rs][i]*scalar_prod*species->get_diam_fluid();
        delta_process_map_send[rs].emplace_back(delta[0], delta[1], delta[2]);
        if ( perform_collision_send[rs][i] )
        {
          ensemble->get_vx( particle1_process_map[rs][i] ) += delta[0];
          ensemble->get_vy( particle1_process_map[rs][i] ) += delta[1];
          ensemble->get_vz( particle1_process_map[rs][i] ) += delta[2];
        }
        if ( numdens_process_map_recv[rs][i] > 0 )
        {
          collision_counter.total++;
          if ( perform_collision_send[rs][i] )
          {
            collision_counter.real++;
            if ( fk > 1.0 )
              collision_counter.out_bound++;
          }
          else
          {
            collision_counter.fake++;
          }
        }
      }
    }
  }
  for (int rs = 0; rs<size; ++rs)
  {
    if ( rs == rank )
    {
      for ( int rr = 0; rr<size; ++rr )
      {
        if ( rr!=rs )
        {
          MPI_Send(&perform_collision_send[rr][0], n_parallel_collisions_send[rr], MPI_INT, rr, 0, MPI_COMM_WORLD);
          MPI_Send(&delta_process_map_send[rr][0], n_parallel_collisions_send[rr], MPI_VELOCITY_TYPE, rr, 0, MPI_COMM_WORLD);
        }
      }
    }
    else
    {
      perform_collision_recv[rs].resize(n_parallel_collisions_recv[rs]);
      MPI_Recv(&perform_collision_recv[rs][0], n_parallel_collisions_recv[rs], MPI_INT, rs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      delta_process_map_recv[rs].resize(n_parallel_collisions_recv[rs]);
      MPI_Recv(&delta_process_map_recv[rs][0], n_parallel_collisions_recv[rs], MPI_VELOCITY_TYPE, rs, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
      for (int i = 0; i<n_parallel_collisions_recv[rs]; ++i)
      {
        if (perform_collision_recv[rs][i])
        {
          ensemble->get_vx(idx_p2) -= delta_process_map_recv[rs][i].vx;
          ensemble->get_vy(idx_p2) -= delta_process_map_recv[rs][i].vy;
          ensemble->get_vz(idx_p2) -= delta_process_map_recv[rs][i].vz;
        }
      }
    }
  }
  return collision_counter;
}

void
CollisionHandler::update_majorants
(void)
{

  if ( (double)n_fake > alpha_1*(double)n_real )
  {
    a11 = anew;
    vrmax11 = vrmaxnew;
  }
  else
  {
    a11 = alpha_2*a11;
    vrmax11 = alpha_2*vrmax11;
  }

}

void
CollisionHandler::reset_parallel_buffers
(void)
{
  for (int r = 0; r<size; ++r)
  {
    n_parallel_collisions_send[r] = 0;
    n_parallel_collisions_recv[r] = 0;
    cells1_process_map[r].clear();
    particle1_process_map[r].clear();
    cells2_process_map_send[r].clear();
    hcells2_process_map_send[r].clear();
    numdens_process_map_send[r].clear();
    aveta_mid_process_map_send[r].clear();
    velocity_process_map_send[r].clear();
    scaled_k_process_map[r].clear();
    delta_process_map_send[r].clear();
    perform_collision_send[r].clear();
  }
}

void
CollisionHandler::commit_velocity_type
(void)
{
  /* Definition of custom data type */
  MPI_Datatype oldtypes[1];
  int blockcounts[1];
  MPI_Aint offsets[1];
  offsets[0] = 0;
  oldtypes[0] = MPI_DOUBLE;
  blockcounts[0] = 3;         // 3 doubles: vx, vy, vz
  MPI_Type_create_struct(1, blockcounts, offsets, oldtypes, &MPI_VELOCITY_TYPE);
  MPI_Type_commit(&MPI_VELOCITY_TYPE);
}

// #######################
// ### DEBUG e TESTING ###
// #######################

void
CollisionHandler::print_subdomain_order
(void) const
{
  for (int i = 0; i<N_COLLISION_SUBDOM; ++i)
    std::cout << subdom_order[i] << "\t";
  std::cout << std::endl;
}

void
CollisionHandler::print_est_collisions
(void) const
{
  for (int r = 0; r<size; ++r) {
    if (r==rank)
      std::cout << "rank " << r << " simulates " << n_coll << " collisions" << std::endl;
  }
  par_env->barrier();
}

void
CollisionHandler::print_reduced_constants
(void) const
{
  for (int r = 0; r<size; ++r) {
    if (r==rank)
      std::cout << "rank " << r << " : sum(a11) = " << a11.sum() << ";\t sum(vrm11) = " << vrmax11.sum() << std::endl;
  }
  par_env->barrier();
}

void
CollisionHandler::test_map_to_domain_collision
(void)
{
  setup_cell_ind_map();
  // PRELIMINARY TEST
  int sum_cells_quarters = 0, sum_ind_quarters = 0;
  for (int q : subdom_order)
  {
    sum_cells_quarters += topology->get_n_cells_quarter(q);
    sum_ind_quarters += cells_ind[q].size();
  }
  std::cout << "rank " << rank << " has " << sum_cells_quarters << "/" << topology->get_n_cells() << "/" << sum_ind_quarters << " cells" << std::endl;
  // COLLISION TEST
  int i, j;
  std::set<std::pair<int,int>> cells_coordinates;
  for (int q : subdom_order)
  {
    for (int k : cells_ind[q])
    {
      i = lexico_inv_quad(k, q).first;
      j = lexico_inv_quad(k, q).second;
      auto p = map_to_domain(i, j, q);
      if ( cells_coordinates.find(p) != cells_coordinates.end() )
      {
        std::cout << "Collision for cell ("<<p.first<<","<<p.second<<"), quarter "<< q <<" of rank "<< rank << std::endl;
      }
      else
        cells_coordinates.insert(p);
    }
  }
}

void
CollisionHandler::gather_collisions
(void)
{
  par_env->all_reduce_inp(n_fake_store[0], MPI_SUM, n_fake_store.size());
  par_env->all_reduce_inp(n_real_store[0], MPI_SUM, n_real_store.size());
  par_env->all_reduce_inp(n_total_store[0], MPI_SUM, n_total_store.size());
  par_env->all_reduce_inp(n_out_store[0], MPI_SUM, n_out_store.size());
}
