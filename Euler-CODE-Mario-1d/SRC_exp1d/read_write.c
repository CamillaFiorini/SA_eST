/***************************************************************************
                                 read_write.c
                                 ------------
            This is read_write: it contains every function reading
                           or writing the solution
                             -------------------
    begin                : Mon May 6 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern int NN, NE, movie, time_step ;
extern int NBF ;

extern double zero_r, gm, time, CFL ;

extern struct element_struct   *element ;
extern struct node_struct      *node ;
extern struct boundary_struct *boundary ;
extern struct boundary_nodes_struct *b_nodes ;


extern char out_file[MAX_CHAR] ;
extern char oneD_file[MAX_CHAR] ;
extern char initial_solution_file[MAX_CHAR] ;

FILE *out ;
FILE *initial ;
FILE *oneD ;

void write_solution()
{
           char paraview_file_name[MAX_CHAR], vigie_file[MAX_CHAR], oneD_file_name[MAX_CHAR] ;	
	   int e, n, TOT_N, TOT_E, lev, i ;
	
	   double rho, u, v, p, Ma, s, T, en, x ;
	    double Ms, Xs, rb, xc, yc, r ;
          double rho1, p1, u1, h1, v1, mm ;
           double rho2, p2, u2, h2, v2, omega, dp, pm, p4, h0  ;

	
	   lev   = 1 ;
	   TOT_N = NN ;
	   TOT_E = NE ;

	    printf( "          **********************************\n" ) ;
            printf( "          **     Writing 1d  file......   **\n" ) ;
            printf( "          **********************************\n" ) ;
		        printf( "\n"                                                ) ;

	 if(time_step == 0) sprintf(paraview_file_name, "../output/initial_state.dat", out_file);
        else  { if ( movie ) sprintf(paraview_file_name, "../output/%s%d.dat", out_file, time_step);
	 else 	sprintf(paraview_file_name, "../output/%s.dat", out_file ); }

	 out = fopen( paraview_file_name, "w" ) ;

        fprintf( out, "# 1d pts  CFL  time -->  %d   %le   %le\n",NN,CFL,time ) ;
        fprintf( out, "# Var. names: x, rho, u, p, Ma, s, tt\n" ) ;
	 for( i = 0 ; i < NN ; i++ )
            {
              x   = node[i].coordinate ;
              rho = node[i].P[2][0] ;
              u   = node[i].P[2][1]/rho ;
              p   = ( gm - 1. )*( node[i].P[2][2] - 0.5*rho*u*u ) ;
              Ma  = fabs( u )/sqrt( gm*fabs( p/rho ) ) ;
              s   = p/pow( rho, gm ) ;
              fprintf( out, "%le  %le  %le %le %le %le %le\n", x,rho,u,p,Ma,s , node[i].tt);
            }

	          fclose( out ) ;
}




/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
