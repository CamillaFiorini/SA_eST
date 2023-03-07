/***************************************************************************
                              upwind_matrices.c
                              -----------------
 This is upwind_matrices: the model-independent computation of the positive
   upwind perameters K_PLUS fully using the eigenspectrum of the system
                              is performed here
                             -------------------
    begin                : Tue May 14 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern struct element_struct *element ;
extern struct node_struct *node ;

extern double zero_r, dt, **Right, **Left, *Lambda, **temp_mat, alpha, speed_of_sound0, total_enthalpy0 ;
extern double ***K_i_p0, ***K_i_p1, *temp_normal, gm, c_tau, ***K1, ref_vel, ubar0, vbar0 ;

extern int size, scheme ;

extern void Eigenvectors() ;       /* normals, variables, RIGHT EIGENVECTORS */
extern void Waves() ;        /* normals, variables, WAVES */
extern void A_times_B( double **, double **, double ** ) ;

void compute_upwind_Keys( int e )
{
     int i, j, vv, kk, n ;
     double  abss, volume, max_v, max_c, max_h, hh, prod ;
     double H, un, u, v, k, nx, ny, gm1, threegm, a ;


    Eigenvectors() ;
    Waves() ;

    hh = 0.01 ; 
    
    H       = total_enthalpy0 ;
    u       = ubar0 ;
    k       = 0.5*u*u ;        
    gm1     = gm - 1. ;
    threegm = 3. - gm ;
    a       = speed_of_sound0 ;

    max_h = max_v = max_c = 0.0 ;

/*
    K1[1][0][0] = 0. ;
    K1[1][0][1] = 1. ;
    K1[1][0][2] = 0. ;
     
    K1[1][1][0] = -threegm*k ;
    K1[1][1][1] = threegm*u ;
    K1[1][1][2] = gm1 ;
     
    K1[1][2][0] = u*( gm1*k - H ) ;
    K1[1][2][1] = H - 2.*gm1*k ;
    K1[1][2][2] = gm*u ;
*/

    for ( i = 0 ; i < size ; i ++ )
        {
          for ( j = 0 ; j < size ; j ++ ) 
              {
                K1[0][i][j] = 0.0 ;
                K1[1][i][j] = 0.0 ;
               
                K_i_p1[1][i][j] = 0.0 ; // 0.5*K[1][i][j] ;
                K_i_p1[0][i][j] = 0.0 ; // -K_i_p1[1][i][j] ;

                for ( kk = 0 ; kk < size ; kk ++ )
                    {
                       abss = fabs( Lambda[kk] ) ;
                       if ( abss < 2.*hh ) abss = hh + 0.25*( abss*abss )/hh ;

                       
                       prod = 0.5*Right[i][kk]*Lambda[kk]*Left[kk][j]/abss ;
                       K1[0][i][j] -= prod ;
                       K1[1][i][j] += prod ;

                       prod = Right[i][kk]*Left[kk][j]*0.5 ;
                       K_i_p1[1][i][j] += prod*( 1. + Lambda[kk]/abss ) ;
                       K_i_p1[0][i][j] += prod*( 1. - Lambda[kk]/abss ) ;
                    }
              }
        }
  

      for ( vv = 0 ; vv < 2 ; vv ++ )
          {
            n = element[e].node[vv] ;
   
            volume = fabs( node[n].Z[2][1] ) ;
            if ( volume > max_v ) max_v = volume ;

            volume = node[n].Z[2][0] ; 
            if ( volume < zero_r ) volume = zero_r ;
 
            volume = sqrt( gm*node[n].Z[2][2]/volume ) ; 
            if ( volume > max_c ) max_c = volume ;
          }


     alpha = max_v + max_c ;
}

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
