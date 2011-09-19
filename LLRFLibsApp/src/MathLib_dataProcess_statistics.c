/****************************************************
 * MathLib_dataProcess_statistics.c
 * 
 * Source file for the mathematic data and routines. For statistics.
 *
 * Created by: Zheqiao Geng, gengzq@slac.stanford.edu
 * Created on: 2011.05.18
 * Description: Initial creation
 ****************************************************/
#include "MathLib_dataProcess.h"

/**
 * Measure the average of the input data
 */
double MATHLIB_avg_double(double *dataIn, int pointNum)
{
    int i;
    double avgResult;
     
    /* Check the input */
    if(!dataIn || pointNum <= 0) return 0;
     
    /* Calculation */
    avgResult = 0;
    for(i = 0; i < pointNum; i ++) 
        avgResult += *(dataIn + i);
     
    avgResult /= pointNum; 
     
    return avgResult;
} 

short MATHLIB_avg_short(short *dataIn, int pointNum)
{
    int i;
    double avgResult;
     
    /* Check the input */
    if(!dataIn || pointNum <= 0) return 0;
     
    /* Calculation */
    avgResult = 0;
    for(i = 0; i < pointNum; i ++) 
        avgResult += (double)(*(dataIn + i));
     
    avgResult /= pointNum; 
     
    return (short)avgResult;
} 

/**
 * Measure root mean square value of the input data
 */
double MATHLIB_rms_double(double *dataIn, int pointNum)
{
    int i;
    double rmsResult;
    
    /* Check the input */
    if(!dataIn || pointNum <= 0) return 0;
     
     /* Calculation */
    rmsResult = 0;	
    for(i = 0; i < pointNum; i ++) 
        rmsResult += *(dataIn + i) * (*(dataIn + i));
     
    rmsResult /= pointNum;
    rmsResult = sqrt(rmsResult); 
     
    return rmsResult;
} 

short MATHLIB_rms_short(short *dataIn, int pointNum)
{
    int i;
    double rmsResult;
    
    /* Check the input */
    if(!dataIn || pointNum <= 0) return 0;
     
     /* Calculation */
    rmsResult = 0;	
    for(i = 0; i < pointNum; i ++) 
        rmsResult += (double)(*(dataIn + i) * (*(dataIn + i)));
     
    rmsResult /= pointNum;
    rmsResult = sqrt(rmsResult); 
     
    return (short)rmsResult;
} 

/**
 * Measure the standard deviation of the input data, taking the definition of
 * std = sqrt(1/n * sum((x-mean(x)).^2))
 */
double MATHLIB_std_double(double *dataIn, int pointNum)
{
    int i;
    double stdResult, tmpAvg;

    /* Check the input */
    if(!dataIn || pointNum <= 0) return 0;
         
    /* get the average */
    tmpAvg = MATHLIB_avg_double(dataIn, pointNum);
    
    /* get the rms with no average */
    stdResult = 0;
    for(i = 0; i < pointNum; i ++)
        stdResult += (*(dataIn + i) - tmpAvg) * (*(dataIn + i) - tmpAvg);
     
    /* stdResult /= (pointNum - 1); */
    stdResult /= pointNum;                /* There are two definitions of standard deviation, which one to use ? */
    stdResult = sqrt(stdResult); 

    return stdResult;
} 

short MATHLIB_std_short(short *dataIn, int pointNum)
{
    int i;
    double stdResult;
    short  tmpAvg;

    /* Check the input */
    if(!dataIn || pointNum <= 0) return 0;
         
    /* get the average */
    tmpAvg = MATHLIB_avg_short(dataIn, pointNum);
    
    /* get the rms with no average */
    stdResult = 0;
    for(i = 0; i < pointNum; i ++)
        stdResult += (double)((*(dataIn + i) - tmpAvg) * (*(dataIn + i) - tmpAvg));
     
    /* stdResult /= (pointNum - 1); */
    stdResult /= pointNum;                /* There are two definitions of standard deviation, which one to use ? */
    stdResult = sqrt(stdResult); 
    
    return (short)stdResult;
} 

/**
 * Calculate the maximum value of the input data
 */
double MATHLIB_max_double(double *dataIn, int pointNum)
{
    int i;
    double result;
    
    /* Check the input */
    if(!dataIn || pointNum <= 0) return 0;
     
    /* Calculation */
    result = *dataIn;
    
    for(i = 1; i < pointNum; i ++)
    {    
        if(*(dataIn + i) > result)
            result = *(dataIn + i);
    }                
    
    return result;        
}

short MATHLIB_max_short(short *dataIn, int pointNum)
{
    int i;
    short result;
    
    /* Check the input */
    if(!dataIn || pointNum <= 0) return 0;
     
    /* Calculation */
    result = *dataIn;
    
    for(i = 1; i < pointNum; i ++)
    {    
        if(*(dataIn + i) > result)
            result = *(dataIn + i);
    }                
    
    return result;
}

/**
 * Calculate the minimum value of the input data
 */
double MATHLIB_min_double(double *dataIn, int pointNum)
{
    int i;
    double result;

    /* Check the input */
    if(!dataIn || pointNum <= 0) return 0;
     
    /* Calculation */   
    result = *dataIn;
    
    for(i = 1; i < pointNum; i ++)
    {    
        if(*(dataIn + i) < result)
            result = *(dataIn + i);
    } 
    
    return result;                      
}

short MATHLIB_min_short(short *dataIn, int pointNum)
{
    int i;
    short result;

    /* Check the input */
    if(!dataIn || pointNum <= 0) return 0;
     
    /* Calculation */   
    result = *dataIn;
    
    for(i = 1; i < pointNum; i ++)
    {    
        if(*(dataIn + i) < result)
            result = *(dataIn + i);
    } 
    
    return result;                      
}

/**
 * Calculate the maximum absolute value of the input data
 */
double MATHLIB_max_abs_double(double *dataIn, int pointNum)
{
    int i;
    double tmpData, result;

    /* Check the input */
    if(!dataIn || pointNum <= 0) return 0;
     
    /* Calculation */   
    result = fabs(*dataIn);
    
    for(i = 1; i < pointNum; i ++)
    {    
        tmpData = fabs(*(dataIn + i));
         
        if(tmpData > result)
            result = tmpData;
    }                
    
    return result;        
}

short MATHLIB_max_abs_short(short *dataIn, int pointNum)
{
    int i;
    short tmpData, result;

    /* Check the input */
    if(!dataIn || pointNum <= 0) return 0;
     
    /* Calculation */   
    result = abs(*dataIn);
    
    for(i = 1; i < pointNum; i ++)
    {    
        tmpData = abs(*(dataIn + i));
         
        if(tmpData > result)
            result = tmpData;
    }                
    
    return result;        
}

