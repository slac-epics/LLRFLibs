/****************************************************
 *  MATH_statistics.c                                         
 * 
 *  Realize the functions to perform statistics calculation
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: April 26, 2009          
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: The first development period until May 10, 2009
 *  Description: Statistics operation
 *               1. average
 *               2. rms, rms with polynomial detrend
 *               3. standard deviation
 *               4. maximum, absolute maximum
 *               5. min
 *  To be done:  add exception handling and run time message  
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 04, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications               
 ****************************************************/
#include "MATH_statistics.h"

/**
 * measure the average of the input data
 */
double MATH_statistic_avg(const double *dataIn, int pointNum)
{
    /* Later, add exception handling here!!! */
	int i;
    double avgResult;
    
    avgResult = 0;
    for(i = 0; i < pointNum; i ++) 
       avgResult += *(dataIn + i);
    
    avgResult /= pointNum; 
    
    return avgResult;
} 

/**
 * measure root mean square value of the input data
 */
double MATH_statistic_rms(const double *dataIn, int pointNum)
{
   /* Later, add exception handling here!!! */
   int i;
   double rmsResult;
   
   rmsResult = 0;	
   for(i = 0; i < pointNum; i ++) 
      rmsResult += *(dataIn + i) * (*(dataIn + i));
    
   rmsResult /= pointNum;
   rmsResult = sqrt(rmsResult); 
    
   return rmsResult;
} 

/**
 * measure the standard deviation of the input data, taking the definition of
 * std = sqrt(1/n * sum((x-mean(x)).^2))
 */
double MATH_statistic_std(const double *dataIn, int pointNum)
{
   /* Later, add exception handling here!!! */
   int i;
   double stdResult, tmpAvg;
   
   /* get the average */
   tmpAvg = MATH_statistic_avg(dataIn, pointNum);
   
   /* get the rms with no average */
   stdResult = 0;
   for(i = 0; i < pointNum; i ++)
      stdResult += (*(dataIn + i) - tmpAvg) * (*(dataIn + i) - tmpAvg);
    
   /* stdResult /= (pointNum - 1); */
   stdResult /= pointNum;            /* There are two definitions of standard deviation, which one to use ? */
   stdResult = sqrt(stdResult); 
   
   return stdResult;
} 

/**
 * calculate the maximum value of the input data
 */
double MATH_statistic_max(const double *dataIn, int pointNum)
{
   int i;
   double result;
   
   result = *dataIn;
   
   for(i = 1; i < pointNum; i ++)
   {   
      if(*(dataIn + i) > result)
         result = *(dataIn + i);
   }            
   
   return result;      
}

/**
 * calculate the minimum value of the input data
 */
double MATH_statistic_min(const double *dataIn, int pointNum)
{
   int i;
   double result;
   
   result = *dataIn;
   
   for(i = 1; i < pointNum; i ++)
   {   
      if(*(dataIn + i) < result)
         result = *(dataIn + i);
   } 
   
   return result;                 
}

/**
 * calculate the maximum absolute value of the input data
 */
double MATH_statistic_max_abs(const double *dataIn, int pointNum)
{
   int i;
   double tmpData, result;
   
   result = fabs(*dataIn);
   
   for(i = 1; i < pointNum; i ++)
   {   
      tmpData = fabs(*(dataIn + i));
       
      if(tmpData > result)
         result = tmpData;
   }            
   
   return result;      
}

/**
 * the rms error with polynomial detrend
 */
double MATH_statistic_rms_detrend(const double *dataIn, int pointNum, int polyOrder)
{
   int i;
   int status;
   double *tmpXData, *tmpYData;                  /* the x, y data for fitting */
   double result;
   
   /* init the variables */	  
   status = MATH_buffer_create(&tmpXData, pointNum);
   if(status != 0) return 0.0;
   status = MATH_buffer_create(&tmpYData, pointNum);
   if(status != 0) return 0.0;

   status = MATH_buffer_init_linear(tmpXData, 0, 1, pointNum);
   
   /* here use the polynomial fitting for interpolation calculation */
   status = MATH_poly_interp_val_vec(tmpXData, dataIn, pointNum, polyOrder, 0, pointNum, tmpXData, tmpYData);   
     
   /* detrend with fitted values */
   for(i = 0; i < pointNum; i ++)
      *(tmpYData + i) = *(dataIn + i) - *(tmpYData + i);
      
   /* get the rms */
   result = MATH_statistic_rms(tmpYData, pointNum);
   
   /* free the memory */ 
   status = MATH_buffer_delete(&tmpXData);
   status = MATH_buffer_delete(&tmpYData);
   
   return result;            
}




