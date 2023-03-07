/***************************************************************************
                                  model_0.c
                                  ---------
  This is model_0: it contains all the model dependent functions used for
        the solution of the 2D EULER equations with PERFECT GAS EOS
                             -------------------
    begin                : Tue May 7 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern struct node_struct *node ;
extern struct element_struct *element ;
extern struct boundary_struct  *boundary ;
extern struct numerical_int_struct numerical_int ;
extern struct boundary_nodes_struct *b_nodes ; 

extern int initial_state, size, face_q_pts, volume_q_pts, NBN, NE, NN ;

extern double *vel, **W, time, dt, **Right, **Left, *Lambda, ***K_i_p1 ;
extern double *phi_a_w0, *phi_a_w1, *phi_node, gm, *phi_w, *temp_vector, *work_vector ;
extern double *work_vector0, *work_vector1, *P_bar, **temp_mat, *FLUX, *temp_normal ;
extern double pressure0, speed_of_sound0, total_enthalpy0, *FLUX, ubar0 ;
extern double pressure1, zero_r, speed_of_sound1, total_enthalpy1 ;

extern void residual_update( int ) ;
/* node global index, local residual distributed to node, time level to update (1/2) */

extern void A_times_v( double **, double *, double * ) ;
extern void A_times_B( double **, double **, double ** ) ;
extern void initvec( double * ) ;
extern void initmat( double ** ) ;
extern void invertmat( double **, double ** ) ;
extern void Eigenvectors() ;
extern void Waves() ;
extern void P_2_Z( double *, double * ) ;
extern void Z_2_P( double *, double * ) ;
extern void first_stremline( int ) ;


void ad_flux( double *P )
{
     FLUX[0] = P[0]*P[1] ;
     FLUX[1] = FLUX[0]*P[1] + P[2] ;
     FLUX[2] = FLUX[0]*( gm*P[2]/( P[0]*( gm - 1.0 ) ) + 0.5*P[1]*P[1] ) ;
}


void P_2_Z( double *P, double *Z )
{
    if ( P[0] < zero_r )
        Z[0] = zero_r ;
    else
    Z[0] = P[0] ;
    if ( P[0] < zero_r )
        { Z[1] = 0. ; }
    else {
    Z[1] = P[1]/P[0] ;/* x - velocity */
    Z[2] = ( gm - 1.0 )*( P[2] - 0.5*P[0]*Z[1]*Z[1] ) ;/*pressure */
    if ( Z[2] < zero_r ) Z[2] = zero_r ;
     }

}

void Z_2_P( double *Z, double *P )
{
    if ( Z[0] < zero_r )
         P[0] = zero_r ;
     else 
    P[0] = Z[0] ;
    if ( Z[0] < zero_r )
 	{
         P[1] = 0. ;
	}
     else {
    P[1] = Z[0]*Z[1] ; }
  
    if ( Z[2] > zero_r )
    P[2] = Z[2]/( gm - 1.0 ) + 0.5*Z[0]*Z[1]*Z[1] ;
    else 
    P[2] = zero_r/( gm - 1.0 ) + 0.5*Z[0]*Z[1]*Z[1] ;
}

/****************************************/
/**      Eigenvalues for 2D Euler      **/
/****************************************/

void Waves()
{
     double edge ;

     Lambda[1] = ubar0 ;
     Lambda[2] = Lambda[1] + speed_of_sound0 ;
     Lambda[0] = Lambda[1] - speed_of_sound0 ;
}


/****************************************/
/**     Eigenvectors for 2D Euler      **/
/****************************************/

void  Eigenvectors()
{
	int m,l ;
      double edge, nx, ny ;
      double rho, u, v, a, p, H, k ;

      p   = pressure0 ;
      rho = P_bar[0] ;
      a   = speed_of_sound0 ;
      u   = ubar0 ;
      k   = u*u*0.5 ;
      H   = total_enthalpy0 ;

	  // u
      Right[0][1] = 1.0 ;
      Right[1][1] = u ;
      Right[2][1] = k ;
   
	  // u + c
      Right[0][2] = rho/a ;
      Right[1][2] = rho*( 1. + u/a ) ;
      Right[2][2] = rho*H/a + rho*u  ;

	  // u - c
      Right[0][0] = rho/a ;
      Right[1][0] = rho*( -1. + u/a ) ;
      Right[2][0] = rho*H/a - rho*u  ;

	  // u
      Left[1][0] = 1.0 - ( gm -		1. )*k/( a*a ) ;
      Left[1][1] = ( gm - 1.)*u/( a*a ) ;
      Left[1][2] = -( gm - 1.)/( a*a ) ;

	  // u + c
      Left[2][0] = 0.5*( gm - 1.)*k/( rho*a ) - 0.5*u/rho ;
      Left[2][1] = -( gm - 1.)*0.5*u/( rho*a ) + 0.5/rho ;
      Left[2][2] = 0.5*( gm - 1.)/( rho*a ) ;

	  // u - c
      Left[0][0] = 0.5*( gm - 1.)*k/( rho*a ) + 0.5*u/rho ;
      Left[0][1] = -( gm - 1.)*0.5*u/( rho*a ) - 0.5/rho ;
      Left[0][2] = 0.5*( gm - 1.)/( rho*a ) ;
}



/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/



