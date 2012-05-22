/****************************************************
 *  DSP_required_interface.h                                         
 * 
 *  The required services by the library of DSP
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: April 27, 2009        
 *
 *  Modified by: Zheqiao Geng @ DESY, zheqiao.geng@desy.de
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
#ifndef DSP_REQUIRED_INTERFACE_H
#define DSP_REQUIRED_INTERFACE_H

#include "LLRF_algorithm_lib_config.h"         /* the common configurations */
#include "MEMM_buffer_operation.h"             /* call the services of Memory Management */

/*------------------------------------------------------------
 * the functions for memory management
 *------------------------------------------------------------*/
int DSP_buffer_create(double **buf, int pointNum);                       /* create a floating buffer and clear to zero */
int DSP_buffer_delete(double **buf);                                     /* delete the newly created buffer */
int DSP_buffer_copy(double *destBuf, double *srcBuf, int pointNum);  /* copy buffers */

#endif

