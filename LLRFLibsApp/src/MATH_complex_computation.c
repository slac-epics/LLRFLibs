/****************************************************
 *  MATH_complex_computation.c                                         
 * 
 *  Realization of the functions to perform complex number calculation
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: April 26, 2009       
 *  
 *  Modified by: Zheqiao Geng @ DESY, zheqiao.geng@desy.de
 *  Modified on: April 28, 2009
 *  Description: Realize the basic path of the complex calculation   
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: The first development period until May 10, 2009
 *  Description: Define the basic operation of the complex vector:
 *               1. create and destroy
 *               2. set element
 *               3. arithmetic (+ - * / scale)     
 *               4. polar and orthogonal
 *               5. conjugate
 *  To be done : add exception handling and run time message         
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 04, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library   
 *  Note: 
 *        a. All the phase data are in the unit of rad!
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications
 ****************************************************/
#include "MATH_complex_computation.h"

/*------------------------------------------------------------
 * the functions to create and destroy the complex numbers
 *------------------------------------------------------------*/
/**
 * create a complex vector, the buffer is created newly
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffer
 *    2 - illegal point number
 *    3 - failed to create buffer
 */
int MATH_complex_create_new(MATH_DATA_COMPLEX *data, int pointNum, MATH_ENUM_COMPLEX_MODE mode)
{
   int status;
   double *real_amp_Buf, *imag_pha_Buf;                      /* temp buffer address */
   
   /* check the input parameters */
   if(!data)          return 1;
   if(pointNum <= 0)  return 2; 
   
   /* create the buffers in the heap */
   status = MATH_buffer_create(&real_amp_Buf, pointNum);
   status = MATH_buffer_create(&imag_pha_Buf, pointNum);
   
   if(status != 0)    return 3;
   
   /* set to the complex data */
   data -> pointNum = pointNum;
   data -> mode     = mode;
   data -> real_amp = real_amp_Buf;
   data -> imag_pha = imag_pha_Buf;
   
   return 0; 
}

/**
 * delete the buffer of a complex vector, only works for the newly created case
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal complex vector
 *    2 - failed to delete buffer
 */ 
int MATH_complex_delete_new(MATH_DATA_COMPLEX *data)
{
   int status;
   
   /* check the input parameters */
   if(!data)   return 1;

   /* delete the buffer */   
   status = MATH_buffer_delete(&(data -> real_amp));
   status = MATH_buffer_delete(&(data -> imag_pha)); 
   data -> pointNum = 0;

   if(!status) return 0;
   else        return 2;
}
                               
/**
 * create a complex vector with external buffers, NEVER delete this one
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffer
 *    2 - illegal point number
 */
int MATH_complex_create_ext(MATH_DATA_COMPLEX *data,       
                            double *real_amp_Buf,
                            double *imag_pha_Buf,
                            int pointNum,
                            MATH_ENUM_COMPLEX_MODE mode)
{
   /* check the input parameters */
   if(!data || !real_amp_Buf || !imag_pha_Buf) return 1;
   if(pointNum <= 0)                           return 2; 
   
   /* init the complex data */	
   data -> pointNum = pointNum;
   data -> mode     = mode;
   data -> real_amp = real_amp_Buf;
   data -> imag_pha = imag_pha_Buf;         
   
   return 0;                      
}

/*------------------------------------------------------------
 * the functions for complex vector data access
 *------------------------------------------------------------*/

/**
 * set elements of the complex vector
 * if the input mode if different with the buffer mode, do the conversion
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal complex vector
 *    2 - illegal id number
 */
int MATH_complex_set_element(MATH_DATA_COMPLEX *data, 
                            int id,
                            MATH_ENUM_COMPLEX_MODE mode,
                            double real_amp_Data,
                            double imag_pha_Data)
{
   /* check the input parameters */
   if(!data)                            return 1;
   if(id < 0 || id >= data -> pointNum) return 2; 
   
   /* set the data */
   if(mode == data -> mode)
   {
      /* if the same, directly set the values */
      *(data -> real_amp + id) = real_amp_Data;
      *(data -> imag_pha + id) = imag_pha_Data;        
   }
   else if(mode == MODE_IQ)
   {
      /* means that the input is in I/Q, while the buffer in A/P */ 
      *(data -> real_amp + id) = sqrt(real_amp_Data * real_amp_Data + imag_pha_Data * imag_pha_Data);
      *(data -> imag_pha + id) = MATH_norm_phase(real_amp_Data, imag_pha_Data);
   }
   else if(mode == MODE_AP)
   {
      /* means that the input is in A/P, while the buffer in I/Q */
      *(data -> real_amp + id) = real_amp_Data * cos(imag_pha_Data);    
      *(data -> imag_pha + id) = real_amp_Data * sin(imag_pha_Data);  
   }
   
   return 0;
}
                               
/*------------------------------------------------------------
 * the functions for complex number computation
 *------------------------------------------------------------*/                               

/**
 * translate the rectagular coordinator to polar coordinator
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal complex vector
 */
int MATH_complex_polar(MATH_DATA_COMPLEX *inputData)
{
   int i;
   double tmpAbs, tmpAngle;
   
   /* check the input parameters */
   if(!inputData || !inputData -> real_amp || !inputData -> imag_pha) return 1;   
   
   /* if it is already polar, do nothing */
   if(inputData -> mode == MODE_AP) return 0;
     
   /* change to polar */
   inputData -> mode = MODE_AP;
     
   for(i = 0; i < inputData -> pointNum; i ++)
   {
      /* the amplitude and phase */
      tmpAbs   = sqrt(*(inputData -> real_amp + i) * (*(inputData -> real_amp + i)) + *(inputData -> imag_pha + i) * (*(inputData -> imag_pha + i)));
      tmpAngle = MATH_norm_phase(*(inputData -> real_amp + i), *(inputData -> imag_pha + i));
         
      /* save to buffer */
      *(inputData -> real_amp + i) = tmpAbs;
      *(inputData -> imag_pha + i) = tmpAngle;
   }    
   
   return 0;
}

/**
 * translate the coordinator polar to rectagular coordinator
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal complex vector
 */
int MATH_complex_orthogonal(MATH_DATA_COMPLEX *inputData)
{
   int i;
   double tmpReal, tmpImag;

   /* check the input parameters */
   if(!inputData || !inputData -> real_amp || !inputData -> imag_pha) return 1; 
     
   /* if it is already orthogonal, do nothing */
   if(inputData -> mode == MODE_IQ) return 0; 
  
   /* change to orthogonal */
   inputData -> mode = MODE_IQ;
   
   for(i = 0; i < inputData -> pointNum; i ++)
   {
      tmpReal = *(inputData -> real_amp + i) * cos(*(inputData -> imag_pha + i));
      tmpImag = *(inputData -> real_amp + i) * sin(*(inputData -> imag_pha + i));   
         
      *(inputData -> real_amp + i) = tmpReal;
      *(inputData -> imag_pha + i) = tmpImag;
   }    
   
   return 0;
}

/**
 * scale the complex vector by a real constant number, the mode of the complex
 * vector does not change (will not change the I/Q or A/P format)
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal complex vector
 */
int MATH_complex_scale(MATH_DATA_COMPLEX *arg1, double arg2, MATH_DATA_COMPLEX *result)
{
   int i;
   MATH_DATA_COMPLEX *tmpResult;

   /* check the input parameters */
   if(!arg1 || !arg1 -> real_amp || !arg1 -> imag_pha) return 1; 
   
   /* select the result storage area, if the buffer of result is empty, use the inital one */
   if(!result) tmpResult = arg1;
   else
   {
      tmpResult          = result;
      result -> pointNum = arg1 -> pointNum;
      result -> mode     = arg1 -> mode;
   }
      
   /* scale them */
   for(i = 0; i < arg1 -> pointNum; i ++)
   {
      if(arg1 -> mode == MODE_IQ)
      {
         /* the I/Q mode */
         *(tmpResult -> real_amp + i) = *(arg1 -> real_amp + i) * arg2;
         *(tmpResult -> imag_pha + i) = *(arg1 -> imag_pha + i) * arg2;
      }
      else
      {
         /* the A/P mode */
         *(tmpResult -> real_amp + i) = *(arg1 -> real_amp + i) * arg2;   
         *(tmpResult -> imag_pha + i) = *(arg1 -> imag_pha + i);
      }
   }
   
   return 0;
}

/**
 * calculate the conjugate of the complex vector, will not change the mode
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal complex vector
 */
int MATH_complex_conj(MATH_DATA_COMPLEX *inputData, MATH_DATA_COMPLEX *result)
{
   int i;
   MATH_DATA_COMPLEX *tmpResult;

   /* check the input parameters */
   if(!inputData || !inputData -> real_amp || !inputData -> imag_pha) return 1;
   
   /* select the result storage area */
   if(!result) tmpResult = inputData;
   else
   {
      tmpResult          = result;
      result -> pointNum = inputData -> pointNum;
      result -> mode     = inputData -> mode;
   } 	
   
   /* calculate the conjugate, same for both IQ and AP modes */
   for(i = 0; i < inputData -> pointNum; i ++)
   {
      *(tmpResult -> real_amp + i) =   (*(inputData -> real_amp + i));
      *(tmpResult -> imag_pha + i) = - (*(inputData -> imag_pha + i));
   }
   
   return 0;
} 

/**
 * add two complex vector, the result will be always in I/Q format, and the two
 * input will also be forced to I/Q mode
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal complex vector
 *    2 - illegal point number (only works when two vector have the same point number or the second one only have 1 point)
 * 
 * Note:
 *    a. the output are always in I/Q format!
 */                      
int MATH_complex_add(MATH_DATA_COMPLEX *arg1, MATH_DATA_COMPLEX *arg2, MATH_DATA_COMPLEX *result)
{
   int i;
   MATH_DATA_COMPLEX *tmpResult;

   /* check the input parameters */
   if(!arg1 || !arg1 -> real_amp || !arg1 -> imag_pha || !arg2 || !arg2 -> real_amp || !arg2 -> imag_pha) return 1;
   
   /* force the two input in I/Q mode */
   MATH_complex_orthogonal(arg1);
   MATH_complex_orthogonal(arg2);
   
   /* select the result storage area, can be stored in arg1 or result buffers */
   if(!result) tmpResult = arg1;
   else
   {
      tmpResult          = result;
      result -> pointNum = arg1 -> pointNum;
      result -> mode     = arg1 -> mode;
   } 

   /* sum them up */
   if(arg1 -> pointNum == arg2 -> pointNum)
   {
      /* if the two input complex buffer have the same length, add them one by one */
      for(i = 0; i < arg1 -> pointNum; i ++)
      {
         *(tmpResult -> real_amp + i) = *(arg1 -> real_amp + i) + *(arg2 -> real_amp + i);
         *(tmpResult -> imag_pha + i) = *(arg1 -> imag_pha + i) + *(arg2 -> imag_pha + i);
      }
   }
   else if(arg2 -> pointNum == 1)
   {
      /* if the second input complex buffer have only one point, also legal */
      for(i = 0; i < arg1 -> pointNum; i ++)
      {
         *(tmpResult -> real_amp + i) = *(arg1 -> real_amp + i) + *(arg2 -> real_amp);
         *(tmpResult -> imag_pha + i) = *(arg1 -> imag_pha + i) + *(arg2 -> imag_pha);
      }
   }
   else 
      return 2;
      
   return 0;
} 

/**
 * subtract two complex vector, the result will be always in I/Q format, and the two
 * input will also be forced to I/Q mode
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal complex vector
 *    2 - illegal point number (only works when two vector have the same point number or the second one only have 1 point)
 * 
 * Note:
 *    a. the output are always in I/Q format!
 */ 
int MATH_complex_sub(MATH_DATA_COMPLEX *arg1, MATH_DATA_COMPLEX *arg2, MATH_DATA_COMPLEX *result)
{
   int i;
   MATH_DATA_COMPLEX *tmpResult;

   /* check the input parameters */
   if(!arg1 || !arg1 -> real_amp || !arg1 -> imag_pha || !arg2 || !arg2 -> real_amp || !arg2 -> imag_pha) return 1;
    
   /* force the two input in I/Q mode */
   MATH_complex_orthogonal(arg1);
   MATH_complex_orthogonal(arg2);
   
   /* select the result storage area */
   if(!result) tmpResult = arg1;
   else
   {
      tmpResult          = result;
      result -> pointNum = arg1 -> pointNum;
      result -> mode     = arg1 -> mode;
   }

   /* subtract */
   if(arg1 -> pointNum == arg2 -> pointNum)
   {
      /* if the two input complex buffer have the same length, sub them one by one */
      for(i = 0; i < arg1 -> pointNum; i ++)
      {
         *(tmpResult -> real_amp + i) = *(arg1 -> real_amp + i) - *(arg2 -> real_amp + i);
         *(tmpResult -> imag_pha + i) = *(arg1 -> imag_pha + i) - *(arg2 -> imag_pha + i);
      }
   }
   else if(arg2 -> pointNum == 1)
   {
      /* if the second input complex buffer have only one point, also legal */
      for(i = 0; i < arg1 -> pointNum; i ++)
      {
         *(tmpResult -> real_amp + i) = *(arg1 -> real_amp + i) - *(arg2 -> real_amp);
         *(tmpResult -> imag_pha + i) = *(arg1 -> imag_pha + i) - *(arg2 -> imag_pha);
      }
   }
   else 
      return 2;
      
   return 0;	
} 

/**
 * multiply two complex vector, the result will be always in I/Q format, and the two
 * input will also be forced to I/Q mode
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal complex vector
 *    2 - illegal point number (only works when two vector have the same point number or the second one only have 1 point)
 * 
 * Note:
 *    a. the output are always in I/Q format!
 */ 
int MATH_complex_mul(MATH_DATA_COMPLEX *arg1, MATH_DATA_COMPLEX *arg2, MATH_DATA_COMPLEX *result)
{
   int i;
   MATH_DATA_COMPLEX *tmpResult;
   double tmpReal;     /* for multiplying, the temp variables are necessary */
   double tmpImag;
   
   /* check the input parameters */
   if(!arg1 || !arg1 -> real_amp || !arg1 -> imag_pha || !arg2 || !arg2 -> real_amp || !arg2 -> imag_pha) return 1;
    
   /* force the two input in I/Q mode */
   MATH_complex_orthogonal(arg1);
   MATH_complex_orthogonal(arg2);
   
   /* select the result storage area */
   if(!result) tmpResult = arg1;
   else
   {
      tmpResult          = result;
      result -> pointNum = arg1 -> pointNum;
      result -> mode     = arg1 -> mode;
   }

   /* multiply */
   if(arg1 -> pointNum == arg2 -> pointNum)
   {
      /* if the two input complex buffer have the same length, multiply them one by one */
      for(i = 0; i < arg1 -> pointNum; i ++)
      {
         tmpReal = *(arg1 -> real_amp + i) * (*(arg2 -> real_amp + i)) - *(arg1 -> imag_pha + i) * (*(arg2 -> imag_pha + i));
         tmpImag = *(arg1 -> real_amp + i) * (*(arg2 -> imag_pha + i)) + *(arg1 -> imag_pha + i) * (*(arg2 -> real_amp + i));      
         
         *(tmpResult -> real_amp + i) = tmpReal;
         *(tmpResult -> imag_pha + i) = tmpImag;
      }
   }
   else if(arg2 -> pointNum == 1)
   {
      /* if the second input complex buffer have only one point, also legal */
      for(i = 0; i < arg1 -> pointNum; i ++)
      {
         tmpReal = *(arg1 -> real_amp + i) * (*(arg2 -> real_amp)) - *(arg1 -> imag_pha + i) * (*(arg2 -> imag_pha));
         tmpImag = *(arg1 -> real_amp + i) * (*(arg2 -> imag_pha)) + *(arg1 -> imag_pha + i) * (*(arg2 -> real_amp));
         
         *(tmpResult -> real_amp + i) = tmpReal;
         *(tmpResult -> imag_pha + i) = tmpImag;
      }
   }	
   else 
      return 2;
      
   return 0;
}                                     

/**
 * divide two complex vector, the result will be always in I/Q format, and the two
 * input will also be forced to I/Q mode
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal complex vector
 *    2 - illegal point number (only works when two vector have the same point number or the second one only have 1 point)
 * 
 * Note:
 *    a. the output are always in I/Q format!
 */  
int MATH_complex_div(MATH_DATA_COMPLEX *arg1, MATH_DATA_COMPLEX *arg2, MATH_DATA_COMPLEX *result)
{
   int i;
   MATH_DATA_COMPLEX *tmpResult;
   double tmpReal;              /* for dividing, the temp variables are necessary */
   double tmpImag;
   double tmpArg2Amp_Square;    /* the amplitude square of the arg2 element, avoid repeating calculation */
   
   /* check the input parameters */
   if(!arg1 || !arg1 -> real_amp || !arg1 -> imag_pha || !arg2 || !arg2 -> real_amp || !arg2 -> imag_pha) return 1;
    
   /* force the two input in I/Q mode */
   MATH_complex_orthogonal(arg1);
   MATH_complex_orthogonal(arg2);
   
   /* select the result storage area */
   if(!result) tmpResult = arg1;
   else
   {
      tmpResult          = result;
      result -> pointNum = arg1 -> pointNum;
      result -> mode     = arg1 -> mode;
   }
   
   /* divide */
   if(arg1 -> pointNum == arg2 -> pointNum)
   {
      /* if the two input complex buffer have the same length, divide them one by one */
      for(i = 0; i < arg1 -> pointNum; i ++)
      {
         tmpArg2Amp_Square = *(arg2 -> real_amp + i) * (*(arg2 -> real_amp + i)) + *(arg2 -> imag_pha + i) * (*(arg2 -> imag_pha + i));
         
         if(tmpArg2Amp_Square != 0)
         {
            tmpReal = (*(arg1 -> real_amp + i) * (*(arg2 -> real_amp + i)) + *(arg1 -> imag_pha + i) * (*(arg2 -> imag_pha + i))) / tmpArg2Amp_Square;
            tmpImag = (*(arg1 -> imag_pha + i) * (*(arg2 -> real_amp + i)) - *(arg1 -> real_amp + i) * (*(arg2 -> imag_pha + i))) / tmpArg2Amp_Square;
         }
         else
         {
            tmpReal = 0;
            tmpImag = 0;    
         }      
         
         *(tmpResult -> real_amp + i) = tmpReal;
         *(tmpResult -> imag_pha + i) = tmpImag;
      }
   }
   else if(arg2 -> pointNum == 1)
   {
      /* if the second input complex buffer have only one point, also legal */
      for(i = 0; i < arg1 -> pointNum; i ++)
      {
         tmpArg2Amp_Square = *(arg2 -> real_amp) * (*(arg2 -> real_amp)) + *(arg2 -> imag_pha) * (*(arg2 -> imag_pha));
         
         if(tmpArg2Amp_Square != 0)
         {
            tmpReal = (*(arg1 -> real_amp + i) * (*(arg2 -> real_amp)) + *(arg1 -> imag_pha + i) * (*(arg2 -> imag_pha))) / tmpArg2Amp_Square;
            tmpImag = (*(arg1 -> imag_pha + i) * (*(arg2 -> real_amp)) - *(arg1 -> real_amp + i) * (*(arg2 -> imag_pha))) / tmpArg2Amp_Square;
         } 
         else
         {
            tmpReal = 0;
            tmpImag = 0;    
         } 
         
         *(tmpResult -> real_amp + i) = tmpReal;
         *(tmpResult -> imag_pha + i) = tmpImag;
      }
   }	
   else 
      return 2;
      
   return 0;
}                 

/*------------------------------------------------------------
 * the functions for common parts
 *------------------------------------------------------------*/

/**
 * measure the normalized phase of input I,Q pare, the phase is normalized to
 * the range of [-pi, pi]
 */
double MATH_norm_phase(double inI, double inQ)
{
   double tmpPhase;
   
   if(inI == 0 && inQ == 0)
      tmpPhase = 0;
   else if(inI == 0 && inQ > 0)
      tmpPhase =   LLRF_pi2;
   else if(inI == 0 && inQ < 0)
      tmpPhase = - LLRF_pi2;
   else if(inI > 0)
      tmpPhase = atan(inQ / inI);      
   else if(inQ > 0)   
      tmpPhase = atan(inQ / inI) + LLRF_pi;
   else
      tmpPhase = atan(inQ / inI) - LLRF_pi;   
      
   return tmpPhase;                      
}

/*------------------------------------------------------------
 * the functions for debugging
 *------------------------------------------------------------*/

/**
 * print the complex vector
 */
void MATH_complex_print(MATH_DATA_COMPLEX *inputData)
{
   int i;
   
   printf("Complex Vector: point number = %d, mode = ", inputData -> pointNum);
   
   if(inputData -> mode == MODE_IQ)     
      printf("I/Q\n");
   else
      printf("A/P\n");
      
   for(i = 0; i < inputData -> pointNum; i ++)
      printf("data%d = (%f, %f)\n", i, *(inputData -> real_amp + i), *(inputData -> imag_pha + i));
   
   printf("\n");
}














