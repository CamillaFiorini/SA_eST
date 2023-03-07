/***************************************************************************
                             boundary_residual.c
                             -------------------
  This is boundary_residual: it drives the computation of the additional
            residual coming from the boundary conditions.
               It is largely based on a weak approach
                             -------------------
    begin                : Wed May 15 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern struct node_struct *node ;
extern struct boundary_struct  *boundary ;
extern struct boundary_nodes_struct *b_nodes ;
extern char grid_file[MAX_CHAR] ;

extern double *temp_normal ;

extern int NN, NBN, initial_state, size, NBF ;

void boundary_conditions()
{
     int f, n ;
      int e1, e2, e3, e4;


     for ( n = 0 ; n < NN ; n ++ ) node[n].flag = 0 ;
                            
           switch ( initial_state )
                  {
                    case 1: // all to zero !!!
           for ( f = 0 ; f < size ; f ++ ) 
               {
                 node[0].Res[f] = 0. ; // Left end
                 node[NN-1].Res[f] = 0. ; // right end
               }
                            break ;
                    case 2: // all to zero !!!
           for ( f = 0 ; f < size ; f ++ ) 
                 node[0].Res[f] = 0. ; // Left end
                            break ;
                    case 4: // all to zero !!!
           for ( f = 0 ; f < size ; f ++ ) 
               {
                 node[0].Res[f] = 0. ; // Left end
                 node[NN-1].Res[f] = 0. ; // right end
               }
                            break ;
                    case 5: // all to zero !!!
                 node[0].Res[1] = 0. ; // Left end
                 node[NN-1].Res[1] = 0. ; // right end
                            break  ;
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
