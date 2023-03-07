/***************************************************************************
                              memory_allocation.c
                              ------------------
                          This is memory_allocation:
               any vector or matrix used in the code was born here
                             --------------------
    begin                : Tue Apr 30 2002
    copyright            : (C) 2002 by Mario Ricchiuto
    email                : ricchiut@vki.ac.be
 ***************************************************************************/
#include "common.h"

extern int size, NE, NBF, thermal_vars, NN ;
extern int NBN, face_q_pts, volume_q_pts ;

extern struct element_struct   *element ;
extern struct node_struct      *node ;
extern struct numerical_int_struct numerical_int ;

extern double *vel, *work_vector, *work_vector0, *work_vector1, *work_vector2, **dU ;
extern double *P_bar, **Lambda, *temp_normal, *normal, **temp_mat, *temp_vector ;
extern double *phi_a_w0, *phi_a_w1, **phi_t_w, *F_loc, *U_C, *work, *PHI_loc, dt ;
extern double **W, ***Z, ***U, **PHI_d, **PHIN_d, **PHIN_scal, **Right, **Left ;
extern double ***K_i, ***K_i_p0, ***K_i_p1 ;
extern double **sum_K, **sum_K_1, *X_temp, *phi_node, *FLUX, *phi_w, *phi1_w ;
extern double ***K0, ***K1, ***K_i_n0, ***K_i_n1, **temp_mat1 ;

extern long int *ipiv ;

extern int    *MA_int_vector( long ) ;
extern double *MA_double_vector( long ) ;
extern double **MA_double_matrix( long , long ) ;
extern double ***MA_double_tensor( long , long, long ) ;

/*******************************************/
/**    Beginning of memory allocation     **/
/*******************************************/

void memory_allocation()
{
	   int e, i, j ;

/**********************************************/
/** memory allocation for structure: element **/
/**********************************************/

	       element = ( struct element_struct * ) calloc( NE + 1, sizeof( struct element_struct ) ) ;
	       for ( e = 0 ; e < NE + 1 ; e++ )
		           element[e].node            = ( int * ) calloc( 2, sizeof( int ) ) ;

/**********************************************/
/**  memory allocation for structure: node   **/
/**********************************************/		

	       node = ( struct node_struct * ) calloc( NN + 1, sizeof( struct node_struct ) ) ;
	       for ( i = 0 ; i < NN + 1 ; i++ )
	           {
	             node[i].P          = ( double ** ) calloc( 3, sizeof( double * ) ) ;
		
		           node[i].P[0]       = ( double * ) calloc( ( ( 3 )*size ), sizeof( double ) ) ;
		           for ( j = 1 ; j < 3 ; j++ )
		                 node[i].P[j] = node[i].P[j-1] + size  ;
		
	             node[i].Z          = ( double ** ) calloc( 3, sizeof( double * ) ) ;

		           node[i].Z[0]       = ( double * ) calloc( ( ( 3 )*size ), sizeof( double ) ) ;
		           for ( j = 1 ; j < 3 ; j++ )
		                 node[i].Z[j] = node[i].Z[j-1] + size  ;
		
		           node[i].Res        = ( double * ) calloc( size, sizeof( double ) ) ;
	           }


/**********************************************/
/**     numerical_integration structure      **/
/**********************************************/

         numerical_int.face_coordinate          = ( double ** ) calloc( face_q_pts, sizeof( double * ) ) ;

         numerical_int.face_coordinate[0]       = ( double * ) calloc( ( face_q_pts*2 ), sizeof( double ) ) ;
         for ( i = 1 ; i < face_q_pts ; i++ )
               numerical_int.face_coordinate[i] = numerical_int.face_coordinate[i-1] + 2 ;

         numerical_int.face_weight              = ( double * ) calloc( face_q_pts, sizeof( double ) ) ;

/**********************************************/
/**          vectors and matrices            **/
/**********************************************/

         temp_vector  = MA_double_vector( size ) ;
         work_vector  = MA_double_vector( size ) ;
         work_vector0 = MA_double_vector( size ) ;
         work_vector1 = MA_double_vector( size ) ;
         work_vector2 = MA_double_vector( size ) ;
         P_bar        = MA_double_vector( size ) ;
         Lambda       = MA_double_matrix( 2, size ) ;
         phi_node     = MA_double_vector( size ) ;

         phi_a_w0    = MA_double_vector( size ) ;
         phi_a_w1    = MA_double_vector( size ) ;

         FLUX   =  MA_double_vector( size ) ;
		
		     PHI_loc      = ( double * ) calloc( size, sizeof( double ) ) ;
		
		     F_loc        = MA_double_vector( size ) ;
         work         = MA_double_vector( size ) ;
         phi_w        = MA_double_vector( size ) ;
         phi1_w       = MA_double_vector( size ) ;
         U_C          = MA_double_vector( size ) ;

         ipiv = ( long int * ) calloc( size, sizeof( long int ) ) ;

         temp_normal = MA_double_vector( 2 + 1 ) ;
         normal      = MA_double_vector( 2 + 1 ) ;
         vel         = MA_double_vector( 2 ) ;
         X_temp      = MA_double_vector( 2 ) ;

         PHI_d       = MA_double_matrix( 3, size ) ;
         PHIN_d      = MA_double_matrix( 3, size ) ;
         PHIN_scal   = MA_double_matrix( 3, size ) ;
         phi_t_w     = MA_double_matrix( 3, size ) ;
         Right       = MA_double_matrix( size, size ) ;
         Left        = MA_double_matrix( size, size ) ;
         sum_K       = MA_double_matrix( size, size ) ;
         sum_K_1     = MA_double_matrix( size, size ) ;
         temp_mat    = MA_double_matrix( size, size ) ;
         temp_mat1   = MA_double_matrix( size, size ) ;
         W           = MA_double_matrix( 3, size ) ;
         Z           = MA_double_tensor( 2, 3, size ) ;
         U           = MA_double_tensor( 2, 3, size ) ;
         dU          = MA_double_matrix( 2, size ) ;


         K_i         = MA_double_tensor( 3, size, size ) ;
         K_i_p0      = MA_double_tensor( 3, size, size ) ;
         K_i_p1      = MA_double_tensor( 3, size, size ) ;
         K_i_n0      = MA_double_tensor( 3, size, size ) ;
         K_i_n1      = MA_double_tensor( 3, size, size ) ;
         K0          = MA_double_tensor( 3, size, size ) ;
         K1          = MA_double_tensor( 3, size, size ) ;
	                                                                                                   	
/*******************************************/
/**       End of memory allocation        **/
/*******************************************/	
}

/*****************************************************/
/**  Memory allocation function for integer vector  **/
/*****************************************************/

int *MA_int_vector( long length )
{
	int *v ;

	v = ( int * ) calloc( length, sizeof( int ) ) ;
	if ( !v ) printf("Memory allocation error \n") ;

	return v ;
}

/*****************************************************/
/**  Memory allocation function for double vector   **/
/*****************************************************/

double *MA_double_vector( long length )
{
	double *v ;

	v = ( double * ) calloc( length, sizeof( double ) ) ;
	if ( !v ) printf("Memory allocation error \n") ;

	return v ;
}

/*****************************************************/
/**  Memory allocation function for double matrix   **/
/*****************************************************/

double **MA_double_matrix( long row, long column )
{
	long i ;

	double **m ;

	m = ( double ** ) calloc( row, sizeof( double * ) ) ;
	if ( !m ) printf("Memory allocation error \n") ;

	m[0] = ( double * ) calloc( ( row*column ), sizeof( double ) ) ;
	if ( !m[0] ) printf("Memory allocation error \n") ;

	for ( i = 1 ; i < row ; i++ ) m[i] = m[i-1] + column ;

	return m ;
}

/*****************************************************/
/**  Memory allocation function for double tensor   **/
/*****************************************************/

double ***MA_double_tensor( long nx, long ny, long nz )
{
	int i ;

	double ***m ;

	m = ( double *** ) calloc( nx, sizeof( double ** ) ) ;
	for ( i = 0 ; i <= nx ; i++ ) m[i] = MA_double_matrix( ny, nz ) ;

	return m ;
}

/***************************************************************************
 *                                                                         *
 *   This program has been developed at the von Karman Institute for Fluid *
 *   Dynamics by Mario Ricchiuto. Any change or update must be done in a   *
 *   clean way and it must be clearly commented.                           *
 *   Mario Ricchiuto                                                       *
 *                                              23/04/2002                 *
 ***************************************************************************/
