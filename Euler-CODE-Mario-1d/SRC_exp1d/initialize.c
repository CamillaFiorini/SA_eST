/***************************************************************************
                                 initialize.c
                                 ------------
This is initialize: it initializes some model, geometry, scheme and case
                      dependent constant and parameters
                             -------------------
    begin                : Wed Apr 24 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern int model, size, time_step_max, scheme ;
extern int vertex ;
extern int NBN, NN, NE, face_q_pts, volume_q_pts ;

extern double c_tau, zero_r ;

extern void read_spatial_grid_header() ;

void initialize()
{

     zero_r = 1.e-10 ;

     NN = NE + 1 ;

     size = 3 ;

     vertex = 2 ;

     face_q_pts = 3 ;
}

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
