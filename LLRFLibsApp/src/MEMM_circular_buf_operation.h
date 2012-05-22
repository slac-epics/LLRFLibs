/****************************************************
 *  MEMM_circular_buf_operation.h                                         
 * 
 *  The functions that handle the circular buffer creation, delete and init
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: May 02, 2009     
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
#ifndef MEMM_CIRCULAR_BUF_OPERATION_H
#define MEMM_CIRCULAR_BUF_OPERATION_H

#include <stdio.h>
#include "MEMM_buffer_operation.h"

/*------------------------------------------------------------
 * Definition of the circular buffer
 *------------------------------------------------------------*/  
typedef struct
{
   int pointNum;                                                                      /* the point number in the circular buffer */
   int currentId;                                                                     /* the current Id of the buffer, always point to the position for next data */
   double *data;                                                                  /* the pointer for data buffer */
} MEMM_DATA_CIRCULAR_BUF;

/*------------------------------------------------------------
 * the functions for circular buffer operation
 *------------------------------------------------------------*/  
int MEMM_circular_buffer_create_new(MEMM_DATA_CIRCULAR_BUF *cirBuf, int pointNum);    /* create a new buffer */
int MEMM_circular_buffer_delete_new(MEMM_DATA_CIRCULAR_BUF *cirBuf);                  /* delete a newly defined buffer */

int MEMM_circular_buffer_create_ext(MEMM_DATA_CIRCULAR_BUF *cirBuf, int pointNum,     /* create a circular buffer with external defined buffer */
                                    double *data);                                   

int MEMM_circular_buffer_init(MEMM_DATA_CIRCULAR_BUF *cirBuf,                         /* init the circular buffer */
                              double  initValue, 
                              double  stepValue,
                              double *sourceBuf,
                              int sel);

int MEMM_circular_buffer_clear(MEMM_DATA_CIRCULAR_BUF *cirBuf);                       /* clear the buffer, set to zero */
int MEMM_circular_buffer_push(MEMM_DATA_CIRCULAR_BUF *cirBuf, double newData);    /* push one value to the buffer */
int MEMM_circular_buffer_reshape(MEMM_DATA_CIRCULAR_BUF *cirBuf);                     /* reshape the buffer to be linear buffer */
                             
/*------------------------------------------------------------
 * the functions for debugging
 *------------------------------------------------------------*/                                                     
void MEMM_circular_buffer_print(MEMM_DATA_CIRCULAR_BUF *cirBuf);

#endif

