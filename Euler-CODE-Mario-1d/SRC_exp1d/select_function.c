/***************************************************************************
                              select_function.c
                              -----------------
          This is select_function: function pointers are assigned here
                             -------------------
    begin                : Thu Apr 25 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern int problem ;
extern int scheme, iteration_max ;

/**************************************************************************/
/**                   SCHEME dependent functions                         **/
/**************************************************************************/

extern void ( *distribute_fluctuation )( int ) ;

extern void UPW_unsteady( int ) ;
extern void LF_unsteady( int ) ;
extern void LUPW_unsteady( int ) ;
extern void B_UPW_unsteady( int ) ;
extern void LLF_unsteady( int ) ;
extern void LW_LF_unsteady( int ) ;


extern void ( *fluctuation2 )( int ) ;

extern void fluctuationRK1( int ) ;
extern void fluctuationRK2( int ) ;
extern void fluctuationRK3( int ) ;

/**************************************************************************/
/**************************************************************************/
/**                  here begins SELECT_FUNCTION                         **/
/**************************************************************************/
/**************************************************************************/

void select_function()
{
/**************************************************************************/
/**                  PROBLEM DEPENDENT SETTINGS !!!                      **/
/**************************************************************************/


     switch ( scheme )
            {
              case 0:
                     distribute_fluctuation = UPW_unsteady ;
                     break ;
              case 1:
                     distribute_fluctuation = LF_unsteady ;
                     break ;
	      case 2:
                     distribute_fluctuation = LUPW_unsteady ;
                     break ;	     
	      case 3:
                     distribute_fluctuation = LLF_unsteady ;
                     break ;	     
	      case 4:
                     distribute_fluctuation = B_UPW_unsteady ;
                     break ;	     
              case 5:
                     distribute_fluctuation = LW_LF_unsteady ;
                     break ;	     
            }

      switch ( iteration_max )
             {
               case 1 :
                        fluctuation2 = fluctuationRK1 ;
                        break ; 
               case 2 :
                        fluctuation2 = fluctuationRK2 ;
                        break ; 
               case 3 :
                        fluctuation2 = fluctuationRK3 ;
                        break ; 
             }

/**************************************************************************/
/**************************************************************************/
/**                    here ends SELECT_FUNCTION                         **/
/**************************************************************************/
/**************************************************************************/
}

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/

