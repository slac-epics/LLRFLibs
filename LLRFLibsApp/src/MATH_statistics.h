/****************************************************
 *  MATH_statistics.h                                         
 * 
 *  The functions to perform statistics calculation
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
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications            
 ****************************************************/
#ifndef MATH_STATISTICS_H
#define MATH_STATISTICS_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "MATH_required_interface.h"
#include "MATH_polynomial.h"

/**
 * The operation to statistic the data in a buffer, be aware that the data are
 * viewed as a single colume vector
 */

double MATH_statistic_avg(const double *dataIn, int pointNum);
double MATH_statistic_rms(const double *dataIn, int pointNum);
double MATH_statistic_std(const double *dataIn, int pointNum);
double MATH_statistic_max(const double *dataIn, int pointNum);
double MATH_statistic_min(const double *dataIn, int pointNum);
double MATH_statistic_max_abs(const double *dataIn, int pointNum);
double MATH_statistic_rms_detrend(const double *dataIn, int pointNum, int polyOrder);

#endif
 

