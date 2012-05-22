/****************************************************
 *  MATH_linear_equations.h   
 *                                      
 *  The functions to solve linear equations
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: April 26, 2009   
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: The first development period until May 10, 2009
 *  Description: Solve the linear equation 
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 04, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications                   
 ****************************************************/
#ifndef MATH_LINEAR_EQUATIONS_H
#define MATH_LINEAR_EQUATIONS_H

#include <stdlib.h>
#include "MATH_required_interface.h"
#include "MATH_matrix_computation.h"

/**
 * Solve the linear equation A0*x = b0
 * x = A0 \ b0
 */

int MATH_linear_equations_solver(MATH_DATA_MATRIX *A0, double *b0, double *x);

#endif
 

