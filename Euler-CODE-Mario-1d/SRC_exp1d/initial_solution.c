/***************************************************************************
                              initial_solution.c
                              ------------------
This is initial_solution: depending on the value of initial_state, it reads
          or just sets up the initial solution for the computation
                             --------------------
    begin                : Mon May 6 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be,


 ***************************************************************************/
#include "common.h"

extern int initial_state ;

extern int NN ;

extern struct node_struct *node ;

extern double gm ;

extern void read_solution() ;
extern void P_2_Z( double *, double * ) ;
extern void Z_2_P( double *, double * ) ;

void initial_solution()
{
     int n ;
     double Ms, Xs, rb, xc, yc, r ;
     double rho1, p1, u1, h1, v1 ;
     double rho2, p2, u2, h2, v2, omega, dp, pm, p4 ;
     double U0, U1, U2, U3 ;
     double U0b, U1b, U2b, U3b ;
     double V0, V1, V2, V3  ;

     switch ( initial_state )
            {
              case 1:/* Riemann Problem */
                    /* rho u v p for Lax*/
                    /*   PRB%VL = (/ 0.445, 0.698, 0.0, 0.0, 3.528 /) */
                      /* PRB%VR = (/   0.5,   0.0, 0.0, 0.0, 0.571 /) */
                     for ( n = 0 ; n < NN ; n ++ )
                         {
                           if ( node[n].coordinate < 0.5 )
                              {
                                  rho1  = 0.445 ;
                                  p1    = 3.528 ;
                                  u1    = 0.698 ;
                                node[n].P[2][0] = rho1 ;
                                node[n].P[2][1] = rho1*u1 ;
                                node[n].P[2][2] = p1/( gm - 1. ) + rho1*0.5*u1*u1 ;
                              }
                           else
                              {
                                  rho1  = 0.5 ;
                                  p1    = 0.571 ;
                                  u1    = 0.0 ;
                                node[n].P[2][0] = rho1 ;
                                node[n].P[2][1] = rho1*u1 ;
                                node[n].P[2][2] = p1/( gm - 1. ) + rho1*0.5*u1*u1 ;
                              }

				 
                           P_2_Z( node[n].P[2], node[n].Z[2] ) ;
                         }
                     break ;

               case 2:/* Q-1D Oscillating Riemann Problem (Shu & Osher) */
                     for ( n = 0 ; n < NN ; n ++ )
                         {
                           if ( node[n].coordinate > -4.0 )
                              {
		                p1 = 1.0 ;
		                u1 = 0.0 ;
		                h1 = 3.5/( 1.0 + 0.2*sin( 5.0*node[n].coordinate ) ) ;
		
		                U0 = p1*gm/( h1*( gm - 1.0 ) ) ;

                                node[n].P[2][0] = U0  ;
                                node[n].P[2][1] = U0*u1 ;
                                node[n].P[2][2] = p1/( gm - 1.0 ) + 0.5*u1*u1*U0 ;
                              }
                           else
                              {
		                p1 = 10.333330 ;
		                u1 = 2.629367 ;
		                h1 = 9.3765398 ;
		
		                U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
		                U3 = h1/gm ;

                                node[n].P[2][0] = U0 ;
                                node[n].P[2][1] = U0*u1 ;
                                node[n].P[2][2] = U3*U0 + 0.5*u1*u1*U0 ;
                              }

                         P_2_Z( node[n].P[2], node[n].Z[2] ) ;
                         }
                     break ;
	            case 3:/*Right Moving Planar Shock: Defined by Ms and Xs*/

                                // Right state

		                rho1 = 1.4 ;
		                p1   = 1.0 ;
		                u1   = 0.0 ;
		                h1   = gm*p1/( ( gm - 1.0 )*rho1 ) ;
		
		                U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
		                U3 = h1/gm ;

		                U1 = U0*u1 ;
		                U2 = U0*v1 ;
		                U3 = U0*U3 + 0.5*u1*u1*U0 ;

                                // Initial position and shock Mach number		    
 
		                Xs   = .50 ;
		                Ms   = 10. ; 

                                // Left state from analytical jump conditions
                                // for right-moving shocks

		                rho2 = rho1*( gm + 1.0 )*Ms*Ms/( 2.0 + ( gm - 1.0 )*Ms*Ms ) ;
		                p2   = p1*( 2.0*gm*Ms*Ms/( gm + 1.0 ) - ( gm -1.0 )/( gm + 1.0 ) ) ;
		                u2   = u1 + 2.0/( gm + 1.0 )*( Ms - 1.0/Ms )*sqrt( gm*p1/rho1 ) ;
		                h2   = gm*p2/( ( gm - 1.0 )*rho2 ) ;

		                V0 = gm*p2/( ( gm - 1.0 )*h2 ) ;
		                V3 = h2/gm ;

		                V1 = V0*u2 ;
		                V2 = V0*v2 ;
		                V3 = V0*V3 + 0.5*u2*u2*V0 ;

				
		                for ( n = 0 ; n < NN ; n ++ )
			                  {
		                      if ( node[n].coordinate < Xs )
		                         {
			                         node[n].P[2][0] = V0 ; 	
			                         node[n].P[2][1] = V1 ;	
			                         node[n].P[2][2] = V3 ;	
		                         }
                          else			  
		                         {
  
					        node[n].P[2][0] = U0 ; 	
			                        node[n].P[2][1] = U1 ;	
			                        node[n].P[2][2] = U3 ;	
		                         }
		
		                         P_2_Z( node[n].P[2], node[n].Z[2] ) ;
			                  }                   	
	                  break ;
              case 4:/*double expansion*/

                                // Right state

		                rho1 = 1. ;
		                p1   = 0.4 ;
		                u1   = 2.0 ;
		                h1   = gm*p1/( ( gm - 1.0 )*rho1 ) ;
		
		                U0 = gm*p1/( ( gm - 1.0 )*h1 ) ;
		                U3 = h1/gm ;

		                U1 = U0*u1 ;
		                U3 = U0*U3 + 0.5*u1*u1*U0 ;

                                // Initial position and shock Mach number		    
 
		                Xs   = 0.5 ;
		                Ms   = 10. ; 

                                // Left state from analytical jump conditions
                                // for right-moving shocks

		                rho2 = 1. ;
		                p2   = 0.4 ;
		                u2   = -2. ;
		                h2   = gm*p2/( ( gm - 1.0 )*rho2 ) ;

		                V0 = gm*p2/( ( gm - 1.0 )*h2 ) ;
		                V3 = h2/gm ;

		                V1 = V0*u2 ;
		                V3 = V0*V3 + 0.5*u2*u2*V0 ;

				
		                for ( n = 0 ; n < NN ; n ++ )
			                  {
		                      if ( node[n].coordinate < Xs )
		                         {
			                         node[n].P[2][0] = V0 ; 	
			                         node[n].P[2][1] = V1 ;	
			                         node[n].P[2][2] = V3 ;	
		                         }
                          else			  
		                         {
  
					        node[n].P[2][0] = U0 ; 	
			                        node[n].P[2][1] = U1 ;	
			                        node[n].P[2][2] = U3 ;	
		                         }
		
		                         P_2_Z( node[n].P[2], node[n].Z[2] ) ;
			                  }                   	
	                  break ;
			case 5:/*Colella&Woodward's blast*/
                    
		                for ( n = 0 ; n < NN ; n ++ )
			                  {
		                      if ( node[n].coordinate < 0.1 )
		                         {
			                         node[n].Z[2][0] = 1. ; 	
			                         node[n].Z[2][1] = 0. ;	
			                         node[n].Z[2][2] = 1000. ;	
		                         }
                                     else {
                                                 if ( node[n].coordinate < 0.9 )			  
		                         		{
					        		node[n].Z[2][0] = 1. ; 	
			                        		node[n].Z[2][1] = 0. ;	
			                        		node[n].Z[2][2] = 0.1 ;	
		                         		}
                                     		else
                                        		{
					        		node[n].Z[2][0] = 1. ; 	
			                        		node[n].Z[2][1] = 0. ;	
			                        		node[n].Z[2][2] = 100. ;	
							}
					  }
		
		                         Z_2_P( node[n].Z[2], node[n].P[2] ) ;
			                  }                   	
	                  break ;
            case 6:/*Shock-contact interaction*/
                    rho1 = 1.4 ;
                    p1   = 1.0 ;
                    u1   = 1.0 ;
                   
    

                            // Initial position and shock Mach number

                    Ms   = 2.95 ;

                            // Left state from analytical jump conditions
                            // for right-moving shocks

                    rho2 = rho1*( gm + 1.0 )*Ms*Ms/( 2.0 + ( gm - 1.0 )*Ms*Ms ) ;
                    p2   = p1*( 2.0*gm*Ms*Ms/( gm + 1.0 ) - ( gm -1.0 )/( gm + 1.0 ) ) ;
                    u2   = u1 + 2.0/( gm + 1.0 )*( Ms - 1.0/Ms )*sqrt( gm*p1/rho1 ) ;
                    
                    printf("rho2=%le, rho1=%le\n",rho2,rho1) ;
                    printf("p2=%le, p1=%le\n",p2,p1) ;
                    printf("u2=%le, u1=%le\n",u2,u1) ;
 
                            for ( n = 0 ; n < NN ; n ++ )
                                  {
                                  if ( node[n].coordinate < 0.35 )
                                     {
                                         node[n].Z[2][0] = rho2 ;
                                         node[n].Z[2][1] = u2 ;
                                         node[n].Z[2][2] = p2 ;
                                     }
                                         else {
                                                     if ( node[n].coordinate < 1.75*0.35 )
                                             {
                                                node[n].Z[2][0] = rho1 ;
                                                node[n].Z[2][1] = u1 ;
                                                node[n].Z[2][2] = p1 ;
                                             }
                                                 else
                                                    {
                                                        node[n].Z[2][0] = rho1/10. ;
                                                                node[n].Z[2][1] = u1 ;
                                                                node[n].Z[2][2] = p1 ;
                                }
                          }
            
                                     Z_2_P( node[n].Z[2], node[n].P[2] ) ;
                                  }
                          break ;
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

