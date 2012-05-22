/****************************************************
 *  MEMM_buffer_operation.h                                                                              
 * 
 *  The functions that handle the buffer creation, delete and init
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: May 05, 2009   
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
#ifndef MEMM_BUFFER_OPERATION_H
#define MEMM_BUFFER_OPERATION_H

#include <stdlib.h>
#include <string.h>

void *MEMM_buffer_create(unsigned int pointNum, unsigned int dataSize);   /* create a new buffer dynamically */
void *MEMM_buffer_delete(void *ptr);                                      /* delete an existing buffer */

int   MEMM_buffer_init(double *destBuf,                               /* init the buffer in several ways */
                       double  initValue,
                       double  stepValue,
                       double *sourceBuf,
                       int pointNum,
                       int sel);

#endif




