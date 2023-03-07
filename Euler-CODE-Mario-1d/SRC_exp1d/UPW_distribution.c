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

extern void velocity() ;
extern void initmat( double ** ) ;
extern void add2mat( double **, double ** ) ;
extern void initvec( double * ) ;
extern void add2vec( double *, double * ) ;
extern void A_times_v( double **, double *, double * ) ;
extern void residual_update( int ) ; /* node global index, local residual distributed to node, time level to update (1/2) */
extern void invertmat( double **, double ** ) ;
extern void decompose( double *, double *, double * ) ;
/* direction for the decomposition, nodal N scheme residual, scalar residuals */
extern void limiter( double *, double **, double *, int, int ) ;
/* scalar components of the cell fluctuation, N scheme scalar residuals, limited scalar nodal residuals, target node, nodes involved in the limiting */
extern void recompose( double *, double *, double * ) ;
/* direction for the decomposition, scalar P scheme residuals, P scheme residual */

extern void Eigenvectors( double * ) ;

void UPW_unsteady( int e )
{
     int j ,v, n1, k, l ;
     double betaP, arg, r_max, u_max ;
     double length, theta, a, b, c, tau ;
     double diss[3], lump, mariuz ;
     double p_max, p_min, coeff ;

     for ( j = 0 ; j < size ; j ++ )
         {
//UPW scheme
           PHIN_d[0][j] = 0. ;
           PHIN_d[1][j] = 0. ;
            
           for ( k = 0 ; k < size ; k ++)
               {
                 PHIN_d[0][j] += dt*K_i_p1[0][j][k]*phi_a_w1[k] ;
                 PHIN_d[1][j] += dt*K_i_p1[1][j][k]*phi_a_w1[k] ;
               }
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