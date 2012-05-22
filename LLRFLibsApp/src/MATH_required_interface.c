/****************************************************
 *  MATH_required_interface.c                                         
 * 
 *  Realize of the required interface for the library of Mathematics
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: May 05, 2009        
 *
 *  Modified by: Zheqiao Geng @ DESY, zheqiao.geng@desy.de
 *  Modified on: The first development period until May 10, 2009
 *  Description: Define the required interface concern to the run
 *               time message generation and memory management  
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 04, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications             
 ****************************************************/
#include "MATH_required_interface.h"

/*------------------------------------------------------------
 * the functions for memory management
 *------------------------------------------------------------*/ 
/**
 * create a buffer in the heap and return the buffer address,
 * when using it, the pointer should be converted to the destination
 * type. simply call the function of library of MEMM
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal point number
 *    2 - failed to create the buffer
 */
int MATH_buffer_create(double **buf, int pointNum)
{
    int status;
    
    /* check the input parametes */  
    if(pointNum <= 0) return 1;
    
    /* create the buffer */
    *buf   = (double *)MEMM_buffer_create(pointNum, sizeof(double));  

    /* init the buffer to zero */
    status = MEMM_buffer_init(*buf, 0, 0, 0, pointNum, 0);
    
    if(!status) return 0;
    else        return 2;
} 

/**
 * delete a buffer, simply call the function of library of MEMM
 * the returned pointer is always 0 (NULL), transfer it to the 
 * pointer ptr for safety
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffer
 */
int MATH_buffer_delete(double **buf)
{
    /* check the input parameters */
    if(!*buf) return 1;
    
    /* delete the buffer */
    *buf  = (double *)MEMM_buffer_delete(*buf);     
    
    return 0;
}           

/**
 * copy the buffers, simply call the function of library of MEMM
 * the two buffer should have the same size or the destination should be
 * larger than the source
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffer
 *    2 - illegal point number
 *    3 - failed to copy
 */
int MATH_buffer_copy(double *destBuf, double *srcBuf, int pointNum)
{
    int status;
    
    /* check the input parameters */
    if(!destBuf || !srcBuf) return 1;
    if(pointNum <= 0)       return 2;

    status = MEMM_buffer_init(destBuf, 0, 0, srcBuf, pointNum, 2);
    
    if(!status) return 0;
    else        return 3;                   
}

/**
 * init one buffer with linear function data, simply call the function of library of MEMM
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffer
 *    2 - illegal point number
 *    3 - failed to init
 */
int MATH_buffer_init_linear(double *destBuf, double initValue, double stepValue, int pointNum)
{  
    int status;
    
    /* check the input parameters */
    if(!destBuf)       return 1;
    if(pointNum <= 0)  return 2;    
                            
    status = MEMM_buffer_init(destBuf, initValue, stepValue, 0, pointNum, 1);                                                                                     
    
    if(!status) return 0;
    else        return 3;     
}                            



