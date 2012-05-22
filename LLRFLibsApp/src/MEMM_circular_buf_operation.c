/****************************************************
 *  MEMM_circular_buf_operation.c                                         
 * 
 *  Realize the functions that handle the circular buffer creation, delete and init
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: May 02, 2009            
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: The first development period until May 10, 2009
 *  Description: The circular buffer operation
 *               1. create and delete            
 *               2. init and clear
 *               3. push data to the circular buffer
 *               4. read data
 *               5. reshape the circular buffer to linear
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 03, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications
 ****************************************************/
#include "MEMM_circular_buf_operation.h"

/**
 * create a circular buffer, the data buffer should be created newly
 * DO NOT forget to delete the buffer after using it!
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal destination buffer
 *    2 - illegal point number, should be larger than 0
 *    3 - failed to create new buffer
 */
int MEMM_circular_buffer_create_new(MEMM_DATA_CIRCULAR_BUF *cirBuf, int pointNum)
{
   /* check the input parameters */   
   if(!cirBuf)       return 1;                                                           /* illegal destination buffer */
   if(pointNum <= 0) return 2;                                                           /* illegal point number, should be larger than 0 */
 
   /* create the new buffer */
   cirBuf -> pointNum  = pointNum;
   cirBuf -> currentId = 0;                                                              /* the initial ID is set to zero */
   cirBuf -> data      = (double *)MEMM_buffer_create(pointNum, sizeof(double));         /* create a new buffer */  
   
   if(!cirBuf -> data) return 3;                                                         /* failed to create the new buffer */

   return 0;
}

/**
 * delete the circular buffer, only used when the data buffer is newly created
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal destination buffer
 */
int MEMM_circular_buffer_delete_new(MEMM_DATA_CIRCULAR_BUF *cirBuf)
{
   /* check the input parameters */   
   if(!cirBuf || !cirBuf -> data) return 1;                                              /* illegal destination buffer */
 
   cirBuf -> pointNum  = 0;
   cirBuf -> currentId = 0;
   cirBuf -> data      = (double *)MEMM_buffer_delete(cirBuf -> data);
        
   return 0;
}

/**
 * create a circular buffer, the data buffer is provided from outside
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal destination buffer
 *    2 - illegal point number, should be larger than 0
 *    3 - illegal external buffer
 */ 
int MEMM_circular_buffer_create_ext(MEMM_DATA_CIRCULAR_BUF *cirBuf, int pointNum, double *data)
{
   /* check the input parameters */   
   if(!cirBuf)       return 1;                      /* illegal destination buffer */
   if(pointNum <= 0) return 2;                      /* illegal point number, should be larger than 0 */
   if(!data)         return 3;                      /* illegal external buffer */
   
   cirBuf -> pointNum  = pointNum;
   cirBuf -> currentId = 0;
   cirBuf -> data      = data;   
   
   return 0;
}

/**
 * init the circular buffer 
 * sel = 0: init the buffer elements with "initValue"
 * sel = 1: init the buffer with a linear function of "initValue + i*stepValue",
 *          where i is the index of elements
 * sel = 2: init the buffer by copying the data from the "sourceBuf"
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal destination buffer
 *    2 - failed to init the buffer
 */
int MEMM_circular_buffer_init(MEMM_DATA_CIRCULAR_BUF *cirBuf,   
                              double  initValue, 
                              double  stepValue,
                              double *sourceBuf,
                              int sel)
{
   int status;

   /* check the input parameters */   
   if(!cirBuf) return 1;              /* illegal destination buffer */ 
    
   /* init the current ID */
   cirBuf -> currentId = 0;

   /* init the buffer */
   status = MEMM_buffer_init(cirBuf -> data, initValue, stepValue, sourceBuf, cirBuf -> pointNum, sel);

   return status;
}

/**
 * clear the circular buffer, set each element zero, without touching the currentId
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal destination buffer
 */
int MEMM_circular_buffer_clear(MEMM_DATA_CIRCULAR_BUF *cirBuf)
{
	 int i;
	 
   /* check the input parameters */   
   if(!cirBuf || !cirBuf -> data) return 1;              /* illegal destination buffer */ 
   
   /* clear the elements to zero */
   for(i = 0; i < cirBuf -> pointNum; i ++) *(cirBuf -> data + i) = 0.0;	
   
   return 0;
} 

/**
 * push new data to the circular buffer
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal destination buffer
 */
int MEMM_circular_buffer_push(MEMM_DATA_CIRCULAR_BUF *cirBuf, double newData)
{      
   /* check the input parameters */   
   if(!cirBuf || !cirBuf -> data) return 1;              /* illegal destination buffer */ 
   
   /* save the data to the current ID */
   *(cirBuf -> data + cirBuf -> currentId) = newData;                          
   cirBuf -> currentId ++;
   
   /* exceed the limitation of the buffer, wrap */
   if(cirBuf -> currentId == cirBuf -> pointNum)                
      cirBuf -> currentId = 0;
}

/**
 * reshape the buffer to linear index, the newest data is always at the end of the buffer
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal destination buffer
 *    2 - illegal point number, should be larger than 0
 *    3 - failed to perform reshaping (errors for creating temp buffers)
 */
int MEMM_circular_buffer_reshape(MEMM_DATA_CIRCULAR_BUF *cirBuf)
{
   double *tmpBuf1, *tmpBuf2;
   
   /* check the input parameters */   
   if(!cirBuf || !cirBuf -> data) return 1;              /* illegal destination buffer */ 
   if(cirBuf -> pointNum <= 0)    return 2;              /* illegal point number, should be larger than 0 */

   /* if the shape is happened to be linear now, do nothing */
   if(cirBuf -> currentId == 0)   return 0;
   
   /* init the variables */
   tmpBuf1 = (double *)MEMM_buffer_create(cirBuf -> currentId, sizeof(double));  
   tmpBuf2 = (double *)MEMM_buffer_create(cirBuf -> pointNum - cirBuf -> currentId, sizeof(double));  

   if(!tmpBuf1 || !tmpBuf2)       return 3;              /* failed to create temp buffers */
   
   /* save the data of the two part int the circular buffer to the two temp buffers */
   memcpy(tmpBuf1,
          cirBuf -> data,
          cirBuf -> currentId * sizeof(double));

   memcpy(tmpBuf2,
          cirBuf -> data + cirBuf -> currentId,
          (cirBuf -> pointNum - cirBuf -> currentId) * sizeof(double));
   
   /* copy the data back for reshape the buffer */
   memcpy(cirBuf -> data,
          tmpBuf2,
          (cirBuf -> pointNum - cirBuf -> currentId) * sizeof(double));
          
   memcpy(cirBuf -> data + cirBuf -> pointNum - cirBuf -> currentId,
          tmpBuf1,
          cirBuf -> currentId * sizeof(double));
          
   cirBuf -> currentId = 0;
   
   /* free the buffer */
   tmpBuf1 = (double *)MEMM_buffer_delete(tmpBuf1);
   tmpBuf2 = (double *)MEMM_buffer_delete(tmpBuf2);   

   return 0;
}

/*------------------------------------------------------------
 * the functions for debugging
 *------------------------------------------------------------*/     

/**
 * print the circular buffer
 */                                                
void MEMM_circular_buffer_print(MEMM_DATA_CIRCULAR_BUF *cirBuf)
{
   int i;
   
   printf("The circular buffer current ID: %d\n", cirBuf -> currentId);
   
   for(i = 0; i < cirBuf -> pointNum; i ++)
      printf("%f ", *(cirBuf -> data + i));
      
   printf("\n");     
}











