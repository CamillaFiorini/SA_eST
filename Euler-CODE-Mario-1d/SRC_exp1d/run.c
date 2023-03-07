/***************************************************************************
                                    run.c
                                    -----
                 This is run: it drives the real computation
                             -------------------
    begin                : Mon May 6 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern int initial_state, time_step, time_step_max, save_solution, info, counter_info  ;

extern double time, t_max ;

extern void initialize_solution() ;
extern void compute_time_step() ;
extern void pseudo_time_stepper() ;
extern void write_solution() ;

void run()
{
     int counter ;

	       printf( "\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "          **                                 **\n" ) ;
         printf( "          **             End of the          **\n" ) ;
         printf( "          **                                 **\n" ) ;
         printf( "          **        pre-processing phase     **\n" ) ;
         printf( "          **                                 **\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "          **                                 **\n" ) ;
         printf( "          ** The real computation starts.... **\n" ) ;
         printf( "          **                                 **\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "          *************************************\n" ) ;
         printf( "\n" ) ;

         counter_info = 0 ;

         time_step = 0 ;

         counter = 0 ;

         time = 0.0 ;

         while ( ( t_max > time + small  )&&( time_step < time_step_max ) )
               {
	
                 time_step++ ;

                 counter ++ ;

                 counter_info ++ ;

                 initialize_solution() ;
                 compute_time_step() ;
                 pseudo_time_stepper() ;

                if ( counter_info == info )		
                   {
                     printf( "                   Time Iteration = %d            \n", time_step ) ;
 		     printf( "                     Time Reached = %le            \n", time ) ;
		     printf( "\n" ) ;

                     counter_info = 0 ;
                   }

                if ( counter == save_solution )
                   {
                     counter = 0 ;

                     write_solution() ;
                   }
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
