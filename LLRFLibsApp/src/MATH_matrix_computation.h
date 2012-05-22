/****************************************************
 *  MATH_matrix_computation.h                                         
 * 
 *  The functions to perform matrix calculation
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: April 26, 2009    
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: The first development period until May 10, 2009
 *  Description: see the .c file                  
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 04, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library   
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Jan. 28, 2010
 *  Description: Change the matrix element read/write function to inline function (failed)
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications
 ****************************************************/
#ifndef MATH_MATRIX_COMPUTATION_H
#define MATH_MATRIX_COMPUTATION_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "MATH_required_interface.h"

/**
 * Macro for get the data ID in the matrix
 * A: pointer to the matrix struct
 * j: row id (start from 0)
 * k: column id (start from 0)
 */
#define MATH_MACRO_MATRIX_DATA_ID(A,j,k) (((A)->majorDim == MAJOR_COL)?(k * (A)->rowNum + j):(j * (A)->colNum + k))

/*------------------------------------------------------------
 * the defintion of the matrix
 *------------------------------------------------------------*/
typedef enum
{
	MAJOR_COL,
	MAJOR_ROW
} MATH_ENUM_MATRIX_MDIM;                        

typedef struct
{
	int rowNum;                             /* the row number, the row ID starts from zero, 0..rowNum-1 */
	int colNum;                             /* the column number, the column ID starts from zero, 0..colNum-1 */
	MATH_ENUM_MATRIX_MDIM majorDim;         /* the major dimention, if 1, data is stored row by row, Note: this is not visible to outside */
	double *data;       
} MATH_DATA_MATRIX;

/*------------------------------------------------------------
 * the functions to create and destroy the matrix
 *------------------------------------------------------------*/
int MATH_matrix_create_new(MATH_DATA_MATRIX *matrix, int rowNum, int colNum);                              /* create a new buffer */
int MATH_matrix_delete_new(MATH_DATA_MATRIX *matrix);                                                      /* delete the new buffer */
int MATH_matrix_create_ext(MATH_DATA_MATRIX *matrix, int rowNum, int colNum, double *extBuf);              /* use existing buffer */
int MATH_matrix_init(MATH_DATA_MATRIX *matrix);                                                            /* init the matrix */

/*------------------------------------------------------------
 * the functions for matrix data access
 *------------------------------------------------------------*/
double MATH_matrix_get_element(MATH_DATA_MATRIX *matrixIn, int rowId, int colId);                          /* read data from matrix */
void   MATH_matrix_set_element(MATH_DATA_MATRIX *matrixIn, int rowId, int colId, double elemData);         /* write data to matrix */
int    MATH_matrix_copy(MATH_DATA_MATRIX *matrixDest, MATH_DATA_MATRIX *matrixSrc);                        /* copy a matrix */

/*------------------------------------------------------------
 * the functions for matrix computation
 *------------------------------------------------------------*/
int MATH_matrix_transpose(MATH_DATA_MATRIX *matrixIn, MATH_DATA_MATRIX *matrixOut);                        /* transpose */
int MATH_matrix_add(MATH_DATA_MATRIX *matrix1, MATH_DATA_MATRIX *matrix2, MATH_DATA_MATRIX *matrixOut);    /* matrix add */
int MATH_matrix_mul(MATH_DATA_MATRIX *matrix1, MATH_DATA_MATRIX *matrix2, MATH_DATA_MATRIX *matrixOut);    /* matrix multiply */
int MATH_matrix_swap_row(MATH_DATA_MATRIX *matrixIn, int rowId1, int rowId2);                              /* exchange two rows */
int MATH_matrix_LU_decomp(MATH_DATA_MATRIX *A, int *permutation);                                          /* LU decomposition */
                   
/*------------------------------------------------------------
 * the functions for debugging
 *------------------------------------------------------------*/
void MATH_matrix_print(MATH_DATA_MATRIX *matrixIn);                                                                                         
                        
#endif
 

