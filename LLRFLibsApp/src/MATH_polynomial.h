/****************************************************
 *  MATH_polynomial.h                                         
 *                                      
 *  The functions to perform polynomial calculation
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
 *  Modified on: Jan. 22, 2010
 *  Description: Remove the scaling of the x and y when doing the polynomial fitting
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications
 ****************************************************/
#ifndef MATH_POLYNOMIAL_H
#define MATH_POLYNOMIAL_H

#include "MATH_required_interface.h"
#include "MATH_matrix_computation.h"
#include "MATH_linear_equations.h"
#include "MATH_statistics.h"

/*------------------------------------------------------------
 * the functions for polynomial fitting
 *------------------------------------------------------------*/
int MATH_poly_fit(const double *X0,
                  const double *Y0,
                  int pointNum,
                  int polyOrder,
                  double *P);
                  
int MATH_poly_val_vec(const double *X,
                  const double *P,
                  int pointNum,
                  int polyOrder,
                  int derOrder,
                  double *valY);                      
                  
double MATH_poly_val(double x,
                  const double *P,
                  int polyOrder,
                  int derOrder);                                                                 

/*------------------------------------------------------------
 * the functions for polynomial interpolation
 *------------------------------------------------------------*/                         
int MATH_poly_interp_val_vec(const double *sourceIndex,
                  const double *sourceData,
                  int sourcePointNum,
                  int sourcePolyOrder,
                  int sourceDerOrder,
                  int interpPointNum,
                  const double *interpIndex,
                  double *interpVal);
                                                          
double MATH_poly_interp_val(const double *sourceIndex,
                  const double *sourceData,
                  int sourcePointNum,
                  int sourcePolyOrder,                                         
                  int sourceDerOrder,
                  double interpIndex);                                                        

#endif

