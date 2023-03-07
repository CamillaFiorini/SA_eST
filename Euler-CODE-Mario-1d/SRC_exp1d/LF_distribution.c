/***************************************************************************
                               N_distribution.c
                               ----------------
  This is N_distribution: here the distribution of the cell fluctuation is
     performed through the N scheme distribution function KP(U_i-U_c)
                             -------------------
    begin                : Tue May 14 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern struct element_struct *element ;
extern struct node_struct *node ;

extern double *work_vector, *work_vector0, *work_vector1, *work_vector2 ;
extern double **PHI_d, **PHIN_d, *phi_a_w0, *phi_a_w1, *phi_w, **PHIN_scal ;
extern double ***K_i_p1, ***K_i_p0, **sum_K, **sum_K_1, ***K1 ;
extern double *temp_vector, *normal, *vel, *phi_node, **Left, **Right, **temp_mat ;
extern double **W,  **Z, *U_C, dt, **phi_t_w, alpha, **dU  ;

extern double diag0, diag1, q_ratio, ubar0, vbar0, speed_of_sound0, shield_factor ;

extern int size, lump_type ;

extern void residual_update( int ) ; /* node global index, local residual distributed to node, time level to update (1/2) */

void LF_unsteady( int e )
{
     int j ,v, n1, k, l ;
     double betaP, arg, r_max, u_max ;
     double length, theta, a, b, c, tau ;
     double diss[3], lump, mariuz ;
     double p_max, p_min, coeff ;

     for ( j = 0 ; j < size ; j ++ )
         {
// LF dissipation
          diss[j]  = W[1][j] - W[0][j] ; diss[j] *= alpha/2. ;

//LF scheme
           PHIN_d[0][j] = dt*phi_a_w1[j]/2. - dt*diss[j] ;
           PHIN_d[1][j] = dt*phi_a_w1[j]/2. + dt*diss[j] ;
         }
// Update residual

    for ( v = 0 ; v < 2 ; v ++ )
        {
          for ( k = 0 ; k < size ; k ++ )
                phi_node[k] = PHIN_d[v][k] ;

            residual_update( element[e].node[v] ) ;
         }
}


/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
