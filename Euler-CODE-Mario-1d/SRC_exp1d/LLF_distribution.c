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

extern int size, lump_type, stab ;

extern void residual_update( int ) ; /* node global index, local residual distributed to node, time level to update (1/2) */
extern void Eigenvectors() ;

void LLF_unsteady( int e )
{
     int j ,v, n1, k, l ;
     double volume, betaP, arg, r_max, u_max ;
     double length, theta, a, b, c, tau ;
     double DDV[2][3], diss[3], lump, mariuz ;
     double p_max, p_min, coeff ;

     lump = 1.0*lump_type ;     

     volume = element[e].volume/2.0 ;

     length = dt*element[e].volume ;
     length = sqrt( length ) ;

     for ( j = 0 ; j < size ; j ++ )
         {
// LF dissipation
          diss[j]  = W[1][j] - W[0][j] ; diss[j] *= alpha/2. ;

//LF scheme
           PHIN_d[0][j] = dU[0][j] + dt*phi_a_w1[j]/2. - dt*diss[j] ;
           PHIN_d[1][j] = dU[1][j] + dt*phi_a_w1[j]/2. + dt*diss[j] ;
//
// Mass matrix corrections
//
          DDV[0][j] = -( 1. - lump )*( 2.*dU[0][j] + dU[1][j] )/3. - lump*dU[0][j] ;
          DDV[1][j] = -( 1. - lump )*( 2.*dU[1][j] + dU[0][j] )/3. - lump*dU[1][j] ;
         }

     Eigenvectors() ;

/* limiting and redistribution of the scalar residuals */


      for ( k = 0 ; k < size ; k ++ )
          {
            work_vector0[k] = 0. ;
            PHIN_scal[0][k] = 0. ;
            PHIN_scal[1][k] = 0. ;

            for ( l = 0; l < size ; l ++ )
                {
                  work_vector0[k] += Left[k][l]*phi_w[l] ;
                  PHIN_scal[0][k] += Left[k][l]*PHIN_d[0][l] ;
                  PHIN_scal[1][k] += Left[k][l]*PHIN_d[1][l] ;
                }
          }

//
// Entropy component: work_vector[1] ;
//
          theta = length*length/( fabs( work_vector0[1] )  + 1.e-20 ) ;
          if ( theta >= 1. )           theta = 1.0 ;
          theta *= stab ;

           mariuz = 1.e-15 ;

     for ( v = 0 ; v < 2 ; v ++ ) 
          node[element[e].node[v]].tt += (1./2.)*theta*element[e].volume/node[element[e].node[v]].vol ;

       for ( v = 0 ; v < 2 ; v ++ )
           {
             for ( j = 0 ; j < size ; j ++ )
                 {
                   if ( fabs( work_vector0[j] ) > 1.e-20 )
                      {
                        arg = PHIN_scal[0][j]/work_vector0[j] ;
                        if ( arg > 0. ) betaP  = arg + mariuz ;
                        else betaP = mariuz ;
                        //betaP = ( arg*arg + arg )/( 1. + arg*arg ) ; //Van Albada

                       arg = PHIN_scal[1][j]/work_vector0[j] ;
                       if ( arg > 0. ) betaP  += arg + mariuz ;
                       else betaP += mariuz ;
                       // betaP += ( arg*arg + arg )/( 1. + arg*arg ) ; //Van Albada


                       arg    = PHIN_scal[v][j]/work_vector0[j] ;
                       if ( arg > 0. ) betaP  = ( arg + mariuz )/betaP ;
                       else betaP = mariuz/betaP ;
//                        betaP = ( ( arg*arg + arg )/( 1. + arg*arg ) )/betaP  ; //Van Albada


                       work_vector1[j] = betaP*work_vector0[j] ;
                     }
                  else work_vector1[j] = 0. ;
                }

            for ( k = 0 ; k < size ; k ++ )
                {
                 //phi_node[k] = DDV[v][k]  + theta*phi_w[k]/2. ; // for Bc
                 phi_node[k] = DDV[v][k] ; // for LLF

		 // Projecting back + dissipation......
                  for ( l = 0 ; l < size ; l ++ )
                  //    phi_node[k] += ( 1. - theta )*Right[k][l]*work_vector1[l] + theta*K1[v][k][l]*phi_w[l] ; // for Bc
                       phi_node[k] += Right[k][l]*work_vector1[l] + theta*K1[v][k][l]*phi_w[l] ; // For LLF
                }

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
