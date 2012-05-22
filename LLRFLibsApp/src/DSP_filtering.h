/****************************************************
 *  DSP_filtering.h                                         
 * 
 *  The functions to perform the filtering
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: April 27, 2009        
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: The first development period until May 10, 2009
 *  Description: see the .c file       
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 03, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library   
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications    
 ****************************************************/
#ifndef DSP_FILTERING_H
#define DSP_FILTERING_H

#include "DSP_required_interface.h"

int DSP_FIR_filter(const double *rawData,                 /* the FIR filter */
                   const double *coef,
                   double *filteredData,
                   int pointNum,
                   int tapNum);

int DSP_FIR_filter_time_inverse(const double *rawData,   /* the time inversed FIR filter */
                   const double *coef,
                   double *filteredData,
                   int pointNum,
                   int tapNum);

int DSP_data_shift(double *data,                         /* shift the data waveform in time */
                   int pointNum, 
                   int shiftNum);

#endif

