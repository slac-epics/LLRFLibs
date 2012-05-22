/****************************************************
 *  MATH_complex_computation.h                                         
 * 
 *  The functions to perform complex number calculation
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: April 26, 2009 
 *
 *  Modified by: Zheqiao Geng @ DESY, zheqiao.geng@desy.de
 *  Modified on: The first development period until May 10, 2009
 *  Description: see the .c file    
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 04, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications          
 ****************************************************/
#ifndef MATH_COMPLEX_COMPUTATION_H
#define MATH_COMPLEX_COMPUTATION_H

#include <math.h>
#include <stdio.h>
#include "MATH_required_interface.h"

/*------------------------------------------------------------
 * the defintion of the complex vector
 *------------------------------------------------------------*/
typedef enum
{
    MODE_IQ,
    MODE_AP
} MATH_ENUM_COMPLEX_MODE;               /* enum for complex coordinator mode */

typedef struct
{
	int pointNum;                       /* the point number of the complex vector */
	MATH_ENUM_COMPLEX_MODE mode;        /* 0 - I/Q mode; 1 - A/P mode */
	double *real_amp;                   /* the buffer for real part or amplitude */
	double *imag_pha;                   /* the buffer for imag part or phase */
} MATH_DATA_COMPLEX;

/*------------------------------------------------------------
 * the functions to create and destroy the complex vector
 *------------------------------------------------------------*/
int MATH_complex_create_new(MATH_DATA_COMPLEX *data, int pointNum,        /* creat new buffers */
                            MATH_ENUM_COMPLEX_MODE mode);
                            
int MATH_complex_delete_new(MATH_DATA_COMPLEX *data);                     /* delete the newly created buffers */
                               
int MATH_complex_create_ext(MATH_DATA_COMPLEX *data,                      /* use existing buffers */
                            double *real_amp_Buf,
                            double *imag_pha_Buf,
                            int pointNum,
                            MATH_ENUM_COMPLEX_MODE mode);

/*------------------------------------------------------------
 * the functions for complex vector data access
 *------------------------------------------------------------*/
int MATH_complex_set_element(MATH_DATA_COMPLEX *data,                     /* set one value in the complex vector */
                            int id,
                            MATH_ENUM_COMPLEX_MODE mode,
                            double real_amp_Data,
                            double imag_pha_Data);
 
/*------------------------------------------------------------
 * the functions for complex vector computation
 *------------------------------------------------------------*/
int MATH_complex_polar(MATH_DATA_COMPLEX *inputData);                                              /* I/Q to A/P */
int MATH_complex_orthogonal(MATH_DATA_COMPLEX *inputData);                                         /* A/P to I/Q */
int MATH_complex_scale(MATH_DATA_COMPLEX *arg1, double arg2, MATH_DATA_COMPLEX *result);           /* scalar scale */
int MATH_complex_conj(MATH_DATA_COMPLEX *inputData, MATH_DATA_COMPLEX *result);                    /* conjugate */   
int MATH_complex_add(MATH_DATA_COMPLEX *arg1, MATH_DATA_COMPLEX *arg2, MATH_DATA_COMPLEX *result); /* + */                           
int MATH_complex_sub(MATH_DATA_COMPLEX *arg1, MATH_DATA_COMPLEX *arg2, MATH_DATA_COMPLEX *result); /* - */                      
int MATH_complex_mul(MATH_DATA_COMPLEX *arg1, MATH_DATA_COMPLEX *arg2, MATH_DATA_COMPLEX *result); /* * */                       
int MATH_complex_div(MATH_DATA_COMPLEX *arg1, MATH_DATA_COMPLEX *arg2, MATH_DATA_COMPLEX *result); /* / */                                                                

/*------------------------------------------------------------
 * the functions for common parts
 *------------------------------------------------------------*/
double MATH_norm_phase(double inI, double inQ);                /* normalize the phase to [-pi, pi] range */ 

/*------------------------------------------------------------
 * the functions for debugging
 *------------------------------------------------------------*/
void MATH_complex_print(MATH_DATA_COMPLEX *inputData);

#endif

