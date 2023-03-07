/***************************************************************************
                              input_parameters.c
                               ----------------
This is input_parameters: it reads and saves all the necessary informations
         needed to set up and run the computation  from the inputfile
                             --------------------
    begin                : Wed Apr 24 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern int model, initial_state, scheme, NE, lump_type ;
extern int iteration_max, convergence_variable, time_step_max, stop, movie, info ;
extern int save_solution, out_format,  start_up, steady_check, out_lev, stab ;

extern double shield_factor, CFL, residual_treshold, t_max, gm, lim_norm, ref_vel, xleft,Lmax ;

extern char initial_solution_file[MAX_CHAR] ;
extern char grid_file[MAX_CHAR] ;
extern char out_file[MAX_CHAR] ;

void input_parameters()
{
      int done ;
      FILE *inputfile ;

    inputfile = fopen( "../textinput/inputfile-exp-1d.txt", "r" ) ;
      if ( !inputfile ) printf( "ERROR: inputfile-exp-1d.txt not found!\n" ) ;

      printf( "          *************************************\n" ) ;
      printf( "          **    Reading the inputfile.....   **\n" ) ;
      printf( "          *************************************\n" ) ;
      printf( "\n" ) ;

      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                                               **\n" ) ;
      fscanf( inputfile, "**           NEO's basic parameters              **\n" ) ;
      fscanf( inputfile, "**                                               **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                  PHYSICS                      **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "   Initial state                             %d    \n", &initial_state ) ;
      fscanf( inputfile, "   Gamma                                     %le    \n", &gm ) ;
      fscanf( inputfile, "   Reference speed                           %le    \n", &ref_vel ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                  GEOMETRY                     **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "   Left end of domain                       %le    \n", &xleft );
      fscanf( inputfile, "   Length of domain                         %le    \n", &Lmax ) ;
      fscanf( inputfile, "   Number of cells                          %d     \n", &NE ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                  NUMERICS                     **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "   Scheme:   %d                                     \n", &scheme ) ;
      fscanf( inputfile, "   Lumping type (0) selective, (1) global   %d    \n", &lump_type );
      fscanf( inputfile, "   Stabilization                            %d      \n",&stab );
      fscanf( inputfile, "   CFL                                      %le    \n", &CFL ) ;
      fscanf( inputfile, "   RK iterations                            %d    \n", &iteration_max ) ;
      fscanf( inputfile, "   Convergence variable (0->Neq.s-1)         %d    \n", &convergence_variable ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                   STOP                        **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "   Final time                               %le    \n", &t_max ) ;
      fscanf( inputfile, "   Maximum number of time steps              %d    \n", &time_step_max ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "**                  OUTPUT                       **\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "   Movie (0/1)                               %d    \n", &movie ) ;
      fscanf( inputfile, "   Steps to save                             %d    \n", &save_solution ) ;
      fscanf( inputfile, "   Information Level (0/1)                   %d    \n", &info ) ;
      fscanf( inputfile, "   Output file                               %s    \n", out_file ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "                      %d                           \n", &done ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;
      fscanf( inputfile, "***************************************************\n" ) ;

      if ( done != 8576 ) printf( "ERROR: inputfile was not read correctly!!" ) ;
      fclose( inputfile ) ;
}

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
