/****************************************************
 *  MATH_required_interface.h                                         
 * 
 *  The required interface for the library of Mathematics
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: May 05, 2009            
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
#ifndef MATH_REQUIRED_INTERFACE_H
#define MATH_REQUIRED_INTERFACE_H

#include "LLRF_algorithm_lib_config.h"         /* the common configurations */
#include "MEMM_buffer_operation.h"             /* the services provided by the Memory Management library */

/*------------------------------------------------------------
 * the functions for memory management
 *------------------------------------------------------------*/
int MATH_buffer_create(double **buf, int pointNum);                       /* create a floating buffer and clear to zero */
int MATH_buffer_delete(double **buf);                                     /* delete the newly created buffer */
int MATH_buffer_copy(double *destBuf, double *srcBuf, int pointNum);      /* copy buffers */

int MATH_buffer_init_linear(double *destBuf,                              /* init the buffer with linear function */
                            double  initValue,
                            double  stepValue,
                            int pointNum);

#endif

