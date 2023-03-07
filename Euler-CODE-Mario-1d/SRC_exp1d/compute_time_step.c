/***************************************************************************
                             compute_time_step.c
                             -------------------
    This is compute_time_step: it computes the time steps associated to
     first and second time level such that the past shield condition
                is satisfied for a linear scalar equation
                             -------------------
    begin                : Mon May 6 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern int NE, size, model, NN ;

extern struct element_struct *element ;
extern struct node_struct *node ;

extern int initial_state, time_step ;

extern double *temp_vector, *normal, *Lambda, CFL ;

extern double dt, shield_factor, time, t_max, q_ratio, gm, zero_r ;

extern void Waves( double * ) ; /* normals, variables, WAVES */
extern void local_average( int ) ;

void compute_time_step()
{
     int e, v, j, w, n ;
     double ratio, delta_t0, volume, max_v, max_c, max_h, speed ;
     

     for ( n = 0 ; n < NN ; n ++ ) node[n].dtau = 0. ;

     q_ratio = dt ;

     if ( time_step == 1 )  q_ratio = 1. ;
     
        delta_t0 = infinity ;

        for ( e = 0 ; e < NE ; e ++ )
            {
              max_v = max_c = 0. ;

              for ( v = 0 ; v < 2 ; v ++ )
                  {
                          n = element[e].node[v] ;

                          volume = fabs( node[n].Z[2][1] ) ;
                          if ( volume > max_v ) max_v = volume ;

                          volume = node[n].Z[2][0] ; if ( volume < zero_r ) volume = zero_r ;
                          volume = sqrt( gm*node[n].Z[2][2]/volume ) ; 
 
                          if ( volume > max_c ) max_c = volume ;
                  }

               speed = ( max_v + max_c ) ;    
 
               node[element[e].node[0]].dtau += speed ;
               node[element[e].node[1]].dtau += speed ;
            }

         for ( n = 0 ; n < NN ; n ++ )
             {
               ratio = node[n].vol/( node[n].dtau + 1.e-20 ) ;
               if ( delta_t0 > ratio ) delta_t0 = ratio ;
             }

         delta_t0 *= CFL ;

        if ( time + delta_t0 > t_max ) delta_t0 = t_max - time ;

        dt = delta_t0 ;

        time += dt ;

	q_ratio /= dt ;
}

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
