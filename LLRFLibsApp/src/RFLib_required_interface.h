/****************************************************
 *  RFLib_required_interface.h                                         
 * 
 *  The required interface for the library of LLRF Applications
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: Dec. 05, 2009
 *  NOTE: All functions here concern to mathematics can be replaced by the functions
 *        in the GNU Scientific Library (GSL), in case of killing bugs
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: June 07, 2010
 *  Description: Optimize the circle fitting, remove the points with large errors
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications 
 ****************************************************/
#ifndef RFLIB_REQUIRED_INTERFACE_H
#define RFLIB_REQUIRED_INTERFACE_H

#include "MEMM_buffer_operation.h"
#include "MEMM_circular_buf_operation.h"
#include "MATH_complex_computation.h"
#include "MATH_linear_equations.h"
#include "MATH_matrix_computation.h"
#include "MATH_polynomial.h"
#include "MATH_statistics.h"
#include "DSP_filtering.h"

#define RFLIB_MACRO_MATRIX_DATA_ID(A,j,k) MATH_MACRO_MATRIX_DATA_ID(A,j,k)

/*------------------------------------------------------------
 * the data type used for interfacing
 *------------------------------------------------------------*/
typedef MATH_DATA_MATRIX  RFLIB_DATA_MATRIX;
typedef MATH_DATA_COMPLEX RFLIB_DATA_COMPLEX;

/*------------------------------------------------------------
 * the functions for memory management
 *------------------------------------------------------------*/
int RFLib_common_buffer_create(double **buf, int pointNum);                       /* create a floating buffer and clear to zero */
int RFLib_common_buffer_delete(double **buf);                                     /* delete the newly created buffer */
int RFLib_common_buffer_copy(double *destBuf, double *srcBuf, int pointNum);      /* copy buffers */

int RFLib_common_buffer_init_linear(double *destBuf,                              /* init the buffer with linear function */
                                    double  initValue,
                                    double  stepValue,
                                    int pointNum);

int RFLib_common_buffer_shift(double *buf, int point_num, int shift_num);         /* shift the buffer */

/*------------------------------------------------------------
 * the functions for mathematic calculation
 * 1. Here use the local library
 * 2. Can also use GSL
 *------------------------------------------------------------*/
/* complex number calculation */
double RFLib_common_complex_abs(double data_I, double data_Q);                         /* amplitude of the single complex value */
int RFLib_common_complex_abs_vector_long(long *buf_I, long *buf_Q, 
                                         double *buf_A, int pointNum);                 /* amplitude of the complex vector, long */
                                             
int RFLib_common_complex_abs_vector_float(double *buf_I, double *buf_Q, 
                                          double *buf_A, int pointNum);                /* amplitude of the complex vector, float */
                                             
double RFLib_common_complex_angle(double data_I, double data_Q);                       /* phase of the single complex value, rad */

int RFLib_common_complex_polar(double *buf_I, double *buf_Q, int pointNum);            /* from IQ to AP, save to the same buffer */
int RFLib_common_complex_polar_d(double *buf_I, double *buf_Q, 
                                 double *buf_A, double *buf_P, int pointNum);          /* from IQ to AP, save to different buffers */

int RFLib_common_complex_orthogonal(double *buf_A, double *buf_P, int pointNum);       /* from AP to IQ, save to the same buffer */

int RFLib_common_complex_add(double *buf_I1, double *buf_Q1,                           /* complex addition */
                             double *buf_I2, double *buf_Q2,
                             double *buf_I3, double *buf_Q3, 
                             int pointNum);

int RFLib_common_complex_div(double *buf_I1, double *buf_Q1,                           /* complex division */
                             double *buf_I2, double *buf_Q2,
                             double *buf_I3, double *buf_Q3, 
                             int pointNum);

int RFLib_common_complex_mul_const(double *buf_I1, double *buf_Q1,                     /* complex multiplification by a constant */
                             double data_I2, double data_Q2,
                             double *buf_I3, double *buf_Q3, 
                             int pointNum);

int RFLib_common_complex_div_const(double *buf_I1, double *buf_Q1,                     /* complex divide by a constant */
                             double data_I2, double data_Q2,
                             double *buf_I3, double *buf_Q3, 
                             int pointNum);

int RFLib_common_complex_sub_scale(double *buf_I1, double *buf_Q1,                     /* complex sub with scale */
                                 double *buf_I2, double *buf_Q2,
                                 double data_scale, int pointNum);

/* matrix calculation */ 
int RFLib_common_matrix_create(RFLIB_DATA_MATRIX *matrix, int rowNum, int colNum);             /* create a new matrix */
int RFLib_common_matrix_delete(RFLIB_DATA_MATRIX *matrix);                                     /* delete the new created matrix */

int RFLib_common_matrix_create_ext(RFLIB_DATA_MATRIX *matrix, int rowNum, int colNum, double *extBuf);  /* create a matrix, with existing buffer */
int RFLib_common_matrix_set_element(MATH_DATA_MATRIX *matrixIn, int rowId, int colId, double elemData); /* set one element of the matrix */
                                 
int RFLib_common_matrix_add(RFLIB_DATA_MATRIX *matrix1, RFLIB_DATA_MATRIX *matrix2, RFLIB_DATA_MATRIX *matrixOut);
int RFLib_common_matrix_mul(RFLIB_DATA_MATRIX *matrix1, RFLIB_DATA_MATRIX *matrix2, RFLIB_DATA_MATRIX *matrixOut);

/* statistics */
double RFLib_common_mea_avg(const double *dataIn, int pointNum);                            /* measure the average */

/* multi-linear fitting and interpolation */
double RFLib_common_linear_interp(const double *dataIn, int pointNum, double dest_x);       /* linear interpolation */
int RFLib_common_linear_fit(RFLIB_DATA_MATRIX *A, double *b, double *x);                    /* multi-variable linear fitting */
int RFLib_common_poly_fit(const double *data, int pointNum, int poly_order, double *coef);  /* polynomial fitting, x is unit vector */

/* digital signal processing */
int RFLib_common_mea_derivative(double *dataIn, double *dataDer, double dt, int pointNum);  /* measure the derivative */

/* circle fitting */
int  RFLib_common_fit_circle(double *x, double *y, int point_num, double *result);
void RFLib_common_gen_circle(double x0, double y0, double a);

int  RFLib_common_fit_circle_opt(double *x, 
                                 double *y, 
                                 int point_num, 
                                 double *result,
                                 double xerr_threshold, 
                                 double yerr_threshold,
                                 double aerr_threshold, 
                                 int max_iteration);

/* cos fitting */
int RFLib_common_fitCos(double *pha_deg, double *value, int pno, double *phaOff_deg, double *amp, double *ampOff);


#endif













