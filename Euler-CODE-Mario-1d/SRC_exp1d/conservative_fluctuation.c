/***************************************************************************
                          conservative_fluctuation.c
                          --------------------------
       This is conservative_fluctuation: here the contour integral of
                    the (space-time) fluxes is performed
                             -------------------
    begin                : Fri May 10 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern struct element_struct *element ;
extern struct node_struct *node ;
extern struct numerical_int_struct numerical_int ;

extern int size, face_q_pts, volume_q_pts, iteration ;

extern double *temp_vector, *temp_normal, *F_loc, *PHI_loc, *X_temp, ***K0, ***K1, *work_vector ;
extern double **dU, **W, dt, **phi_t_w, *phi_w, *phi_a_w0, *phi_a_w1, *FLUX  ;
extern double  *work_vector1, *work_vector0, q_ratio ;
extern double a1[3], a2[3], a3[3] ;
extern double b1[3], b2[3], b3[3] ;


extern void ( *fluctuation2 )( int ) ;
extern void initvec( double * ) ;
extern void A_times_v( double **, double *, double * ) ;
extern void add2vec( double *, double * ) ;
extern void ad_flux( double * ) ;
extern void Z_2_P( double *, double * ) ;
extern void P_2_Z( double *, double * ) ;

/*********************************************/
/*********************************************/
/**      unsteady homogeneous systems       **/
/*********************************************/
/*********************************************/

void STfluctuation2( int e )
{
     int i, v, itr, i1, i2, i3, vv ;
     double V_T ;

     itr = iteration - 1 ;

     V_T = element[e].volume/2.0 ;


     

     for ( v = 0 ; v < 2 ; v ++ )
         {
           vv = element[e].node[v] ;
           for ( i = 0 ; i < size ; i ++ )
                 dU[v][i] = V_T*( b3[itr]*node[vv].P[2][i] + b2[itr]*node[vv].P[1][i] + b1[itr]*node[vv].P[0][i] ) ;
         }

     fluctuation2( e ) ;

     for ( i = 0 ; i < size ; i ++ )
          phi_w[i] = dU[0][i] + dU[1][i] + dt*phi_a_w1[i] ;
}

void fluctuationRK1( int e ) // only P1
{
     int f, i, i1, i2, i0, itr ;

     i0 = element[e].node[0] ;
     i1 = element[e].node[1] ;

     P_2_Z( node[i0].P[0], work_vector0 ) ;
     ad_flux( work_vector0 ) ;
     for( i = 0 ; i < size ; i ++ )
          phi_a_w1[i] = -FLUX[i] ;
     P_2_Z( node[i1].P[0], work_vector0 ) ;
     ad_flux( work_vector0 ) ;
     for( i = 0 ; i < size ; i ++ )
          phi_a_w1[i] += FLUX[i] ; 
}


void fluctuationRK2( int e )//only P1
{
     int f, i, i1, i2, i0, itr ;


     itr = iteration - 1 ;

     switch ( itr )
        {    
         case 0 : 
         fluctuationRK1( e ) ;
         break ;
         case 1 :
         i0 = element[e].node[0] ;
         i1 = element[e].node[1] ;

// time tn

         P_2_Z( node[i0].P[0], work_vector0 ) ;
         ad_flux( work_vector0 ) ;
         for( i = 0 ; i < size ; i ++ )
         phi_a_w1[i] = -0.5*FLUX[i] ;
         P_2_Z( node[i1].P[0], work_vector0 ) ;
         ad_flux( work_vector0 ) ;
         for( i = 0 ; i < size ; i ++ )
         phi_a_w1[i] += 0.5*FLUX[i] ;

// time t1
//       
         P_2_Z( node[i0].P[2], work_vector0 ) ;
         ad_flux( work_vector0 ) ;
         for( i = 0 ; i < size ; i ++ )
         phi_a_w1[i] -= 0.5*FLUX[i] ;
         P_2_Z( node[i1].P[2], work_vector0 ) ;
         ad_flux( work_vector0 ) ;
         for( i = 0 ; i < size ; i ++ )
         phi_a_w1[i] += 0.5*FLUX[i] ;
         break ;
// end switch  
 
       }
}

void fluctuationRK3( int e )// only P1
{
     int f, i, i1, i2, i0, itr ;


     itr = iteration - 1 ;

     switch ( itr )
        {    
         case 0 : 
         fluctuationRK1( e ) ;
         break ;
         case 1 :
         fluctuationRK2( e ) ;
         for ( i = 0 ; i < size ; i ++ ) 
               phi_a_w1[i] *= 0.5 ;
         break ;
         case 2 :
         i0 = element[e].node[0] ;
         i1 = element[e].node[1] ;

// time tn
// 
         
         P_2_Z( node[i0].P[0], work_vector0 ) ;
         ad_flux( work_vector0 ) ;
         for( i = 0 ; i < size ; i ++ )
         phi_a_w1[i] = -FLUX[i]/6. ;
         P_2_Z( node[i1].P[0], work_vector0 ) ;
         ad_flux( work_vector0 ) ;
         for( i = 0 ; i < size ; i ++ )
         phi_a_w1[i] += FLUX[i]/6. ;

// time t2
//
         P_2_Z( node[i0].P[2], work_vector0 ) ;
         ad_flux( work_vector0 ) ;
         for( i = 0 ; i < size ; i ++ )
         phi_a_w1[i] -= 2.*FLUX[i]/3. ;
         P_2_Z( node[i1].P[2], work_vector0 ) ;
         ad_flux( work_vector0 ) ;
         for( i = 0 ; i < size ; i ++ )
         phi_a_w1[i] += 2.*FLUX[i]/3. ;         
 
// time t1
//
         P_2_Z( node[i0].P[1], work_vector0 ) ;
         ad_flux( work_vector0 ) ;
         for( i = 0 ; i < size ; i ++ )
         phi_a_w1[i] -= FLUX[i]/6. ;
         P_2_Z( node[i1].P[1], work_vector0 ) ;
         ad_flux( work_vector0 ) ;
         for( i = 0 ; i < size ; i ++ )
         phi_a_w1[i] += FLUX[i]/6. ;

// end switch  
 
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
