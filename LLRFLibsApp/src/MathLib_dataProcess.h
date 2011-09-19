/****************************************************
 * MathLib_dataProcess.h
 * 
 * Header file for the mathematic data and routines
 *
 * Created by: Zheqiao Geng, gengzq@slac.stanford.edu
 * Created on: 2011.05.18
 * Description: Initial creation
 ****************************************************/
#ifndef MATHLIB_DATA_PROCESS_H
#define MATHLIB_DATA_PROCESS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/*======================================
 * Constant definitions
 *======================================*/

/*======================================
 * Data models
 *======================================*/

/*======================================
 * Math algorithms
 *======================================*/
/* Common */
#define MATHLIB_max(a, b)   ((a>b)?a:b)
#define MATHLIB_min(a, b)   ((a<b)?a:b)

/* Bits calculation */
#define MATHLIB_u32ToShortHi(data_ptr)  (*((short *)(data_ptr) + 1))        /* get the higher 16 bits of a unsigned int data as a signed short */
#define MATHLIB_u32ToShortLo(data_ptr)  (*((short *)(data_ptr)))            /* get the lower 16 bits of a unsigned int data as a signed short */

int MATHLIB_u32ToShortArray(unsigned int *data, int pointNum, short *arrayHi, short *arrayLo);
int MATHLIB_getEveryNSubArray(short *data, int pointNum, short *subArray, int idOffset, int everyN);

/* Statistics */
double MATHLIB_avg_double(double *dataIn, int pointNum);                    /* double version */
double MATHLIB_rms_double(double *dataIn, int pointNum);
double MATHLIB_std_double(double *dataIn, int pointNum);
double MATHLIB_max_double(double *dataIn, int pointNum);
double MATHLIB_min_double(double *dataIn, int pointNum);
double MATHLIB_max_abs_double(double *dataIn, int pointNum);

short MATHLIB_avg_short(short *dataIn, int pointNum);                       /* short version */
short MATHLIB_rms_short(short *dataIn, int pointNum);
short MATHLIB_std_short(short *dataIn, int pointNum);
short MATHLIB_max_short(short *dataIn, int pointNum);
short MATHLIB_min_short(short *dataIn, int pointNum);
short MATHLIB_max_abs_short(short *dataIn, int pointNum);

#ifdef __cplusplus
}
#endif

#endif




