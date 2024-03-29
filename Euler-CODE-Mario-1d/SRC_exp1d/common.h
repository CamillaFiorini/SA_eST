/***************************************************************************
                                  common.h
                                  -------
This is the common: it defines the data structures used in the computations
                             and some constants.
                             -------------------
    begin                : Wed Apr 24 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <malloc/malloc.h>
#include "f2c.h"
#include "clapack.h"

#define MAX(a,b) ((a)>(b) ? (a) : (b))
#define MIN(a,b) ((a)<(b) ? (a) : (b))

#define pi  3.14159265358979
#define MAX_CHAR 150
#define infinity ( double )( 1.e+50 )
#define epsilon ( double )( 1.e-14 )
#define small ( double )( 1.e-10 )
#define SIGN(a)  ((a)<=0 ? -1 : 1)
#define MINMOD(a) ( 0.5*( 1.0 + SIGN((a)) )*MIN(1,(a)) )



struct element_struct
{
       int *node ;
       double volume ;
} ;

struct node_struct
{
       int flag ;
       double **P ;
       double **Z ;
       double *Res ;
       double coordinate ;
       double dtau ;
       double vol ;
       double tt ;
} ;

struct numerical_int_struct
{
       double **face_coordinate ;
       double *face_weight ;
} ;



/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
