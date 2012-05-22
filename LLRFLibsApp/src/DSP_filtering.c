/****************************************************
 *  DSP_filtering.c                                         
 * 
 *  Realize the functions to perform the filtering
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: April 27, 2009    
 * 
 *  Modified by: Zheqiao Geng
 *  Modified on: April 27, 2009
 *  Description: Realize the basic path FIR filter
 *  To be done:  add the exception handling  
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: The first development period until May 10, 2009
 *  Description: Realize the time inversed FIR filter and data series time shift
 *               Implement the basic exception handling and run time message generation
 *  To be done:  1. detailed description of each function
 *               2. put "const" keyword before each input buffer pointer in order to avoid
 *                  unexpected change of the buffer content 
 *                  -- THIS SHOULD BE VALID FOR ALL FUNCTIONS  
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 03, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications                               
 ****************************************************/
#include "DSP_filtering.h"

/**
 * the direct implementation of the FIR filter
 *
 * Parameters:
 *       rawData:      the floating data series to be filtered
 *       coef:         the FIR coefficients
 *       filteredData: the floating data buffer to store the results
 *       pointNum:     the point number of the input data
 *       tapNum:       the tap number of the FIR filter
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffers
 *    2 - illegal point number or FIR tap number
 *
 * Note:
 *    a. The group delay has been compensated here
 */
int DSP_FIR_filter(const double *rawData,   
                   const double *coef,
                   double *filteredData,
                   int pointNum,
                   int tapNum)
{
   int i, j, half_tapNum;

   /* check the input parameters */
   if(!rawData || !coef || !filteredData)                return 1;
   if(pointNum <= 0 || tapNum <= 0 || tapNum > pointNum) return 2;
   
   /* calculate the shift number for group delay compensation */
   half_tapNum = tapNum / 2;
   
   /* calculate FIR filtering, the group delay is compensated */
   for(i = tapNum - 1; i < pointNum; i ++)
   {
      *(filteredData + i - half_tapNum) = 0.0;      
      for(j = 0; j < tapNum; j ++)
         *(filteredData + i - half_tapNum) += *(coef + j) * (*(rawData + i - j));
   }
   
   return 0;	
} 

/**
 * the FIR filter with the raw data inversed in time
 * 1. inverse the raw data
 * 2. FIR filter (note: for derivative, the sign may change by the time inverse)
 * 3. Inverse the result
 *
 * But in the implementation, the steps above is not followed, because it can be
 * done by arrange the address in a clever way :)
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffers
 *    2 - illegal point number or FIR tap number
 *
 * Note:
 *    a. The group delay has been compensated here
 */
int DSP_FIR_filter_time_inverse(const double *rawData,
                   const double *coef,
                   double *filteredData,
                   int pointNum,
                   int tapNum)
{                            
   int i, j, half_tapNum;

   /* check the input parameters */
   if(!rawData || !coef || !filteredData)                return 1;
   if(pointNum <= 0 || tapNum <= 0 || tapNum > pointNum) return 2;

   /* calculate the shift number for group delay compensation */
   half_tapNum = tapNum / 2;
      
   /* time inversed FIR filter, the group delay is compensated */                     
   for(i = pointNum - tapNum; i >= 0; i --)
   {
      *(filteredData + i + half_tapNum) = 0.0;      
      for(j = 0; j < tapNum; j ++)
         *(filteredData + i + half_tapNum) += *(coef + j) * (*(rawData + i + j));
   } 
   
   return 0;                        
}                       

/**
 * shift the data buffer, negative number is for left, positive number for right
 * after shifting, the empty buffer is filled by zero
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffers
 *    2 - illegal point number
 *    3 - failed to perform data shift
 *
 * To be done: more complete exception handling
 */                            
int DSP_data_shift(double *data, int pointNum, int shiftNum)
{
   int status; 
   double *tmpBuf;

   /* check the input parameters */
   if(!data)         return 1;
   if(pointNum <= 0) return 2;
   
   /* if the shift number is zero, do nothing */
   if(shiftNum == 0) return 0;
   
   /* init the temp buffer, and clear to zero */
   status = DSP_buffer_create(&tmpBuf, pointNum);

   /* copy the data to the temp buffer according to the shift setting */
   if(shiftNum > 0)                                                        /* right shifting */
      status = DSP_buffer_copy(tmpBuf + shiftNum, data, pointNum - shiftNum);
   else                                                                    /* left shifting */
      status = DSP_buffer_copy(tmpBuf, data - shiftNum, pointNum + shiftNum);

   /* copy the temp buffer to the destination */
   status = DSP_buffer_copy(data, tmpBuf, pointNum); 
     
   /* delete the buffer */
   status = DSP_buffer_delete(&tmpBuf);      
   
   if(!status) return 0;
   else        return 3;
}                           
                            
/*===========================================================
 * Note: there is no difference between the FIR filter and
 *       time inversed FIR filter if we consider the group delay
 *===========================================================*/ 
                            

