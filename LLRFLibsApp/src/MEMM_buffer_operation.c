/****************************************************
 *  MEMM_buffer_operation.c                                                                             
 * 
 *  Realize the functions that handle the buffer creation, delete and init
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: May 05, 2009   
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: The first development period until May 10, 2009
 *  Description: The buffer operation
 *               1. create and delete
 *               2. buffer data init   
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 03, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library, realize
 *               the function to create, delete and init the buffer     
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications               
 ****************************************************/
#include "MEMM_buffer_operation.h"

/**
 * create a one-dimensional buffer, realize the function of:
 * #include <stdlib.h>
 * void *calloc(size_t num, size_t size);
 *
 * This is a basic function, do not need to collect run time message internally
 */
void *MEMM_buffer_create(unsigned int pointNum, unsigned int dataSize)
{
   return calloc(pointNum, dataSize);     
}

/**
 * delete a one-dimensional buffer, realize the function of:
 * #include <stdlib.h>
 * void free(void *ptr);
 *
 * This is a basic function, do not need to collect run time message internally
 */
void *MEMM_buffer_delete(void *ptr)
{
   free(ptr);
   return 0;      
}

/**
 * init the one-dimensional floating point buffer
 * sel = 0: init the buffer elements with "initValue"
 * sel = 1: init the buffer with a linear function of "initValue + i*stepValue",
 *          where i is the index of elements
 * sel = 2: init the buffer by copying the data from the "sourceBuf"
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal destination buffer
 *    2 - illegal point number, should be larger than 0
 *    3 - illegal source buffer if it is used to init the destination
 *    4 - illegal user selection, do nothing
 */
int MEMM_buffer_init(double *destBuf,                             
                     double  initValue,
                     double  stepValue,
                     double *sourceBuf,
                     int pointNum,
                     int sel)
{
   int i; 
  
   /* check the input parameters */
   if(!destBuf)                  return 1;                                             /* illegal destination buffer */
   if(pointNum <= 0)             return 2;                                             /* illegal point number, should be larger than 0 */
   if(sel == 2 && !sourceBuf)    return 3;                                             /* illegal source buffer if it is used to init the destination */
   
   /* init the buffer based on the selection */
   switch(sel)
   {
       case 0:                                                                         /* init the buffer with initValue */
           for(i = 0; i < pointNum; i ++) *(destBuf + i) = initValue;
           break;
           
       case 1:                                                                         /* init the buffer with linear function */
           for(i = 0; i < pointNum; i ++) *(destBuf + i) = initValue + i * stepValue;
           break;
           
       case 2:                                                                         /* init the buffer by copying another buffer */
           memcpy(destBuf, sourceBuf, pointNum * sizeof(double));
           break;
           	
       default:
           return 4;    	                                                           /* illegal user selection, do nothing */
   }
   
   return 0;                                                                           /* succeed */   
} 


