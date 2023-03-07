/***************************************************************************
                           element_preprocessing.c
                           -----------------------
This is element_preprocessing: the geometry of each element is computed and
 processed here to have all the necessary informations for the computation
                             -------------------
    begin                : Thu May 2 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern int NN, NE, NBF, NBN, initial_state ;
extern double xleft, Lmax ;

extern struct element_struct *element ;
extern struct node_struct    *node ;
extern struct boundary_struct *boundary ;
extern struct boundary_nodes_struct   *b_nodes ; 
extern char grid_file[MAX_CHAR] ;

/************************************************/
/**  Preprocessing for linear P1 interpolation **/
/************************************************/

void element_preprocessing()
{
     int e, n, n1, n2, f, m, m1, m2 ;
     double x1, y1, nx1, ny1, nx2, ny2  ;
     double x2, y2, dummy;
      int  ver, e1, e2, e3, e4 ;

// Generation of the nodal coordinates  
//
     for ( n = 0 ; n < NN ; n ++ )
         {
           node[n].coordinate = xleft + Lmax*n/( 1.*NE ) ;
           node[n].vol        = 0. ;
           node[n].flag       = 1 ;
         }

/*******************************************/
/**    Computation of the nodal nomals    **/
/*******************************************/
 
     
     for ( e = 0 ; e < NE ; e++ )
	 {
           element[e].node[0] = e ;
           element[e].node[1] = e + 1 ;

    	   element[e].volume = node[e+1].coordinate - node[e].coordinate ;

           node[e].vol   += element[e].volume/2. ;
           node[e+1].vol += element[e].volume/2. ;
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
