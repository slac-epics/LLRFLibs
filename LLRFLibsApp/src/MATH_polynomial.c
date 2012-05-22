/****************************************************
 *  MATH_polynomial.c                                         
 *                                      
 *  Realize the functions to perform polynomial calculation
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: April 26, 2009               
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: The first development period until May 10, 2009
 *  Description: Polynomial operation
 *               1. polynomial fitting (polyfit)
 *               2. polynomial value and derivative (work together with the polyfit)
 *               3. polynomial value and derivative interpolation (based on 2)      
 *  To be done:  add exception handling and run time message  
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 04, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library   
 *   
 *  Modified by: Zheqiao Geng
 *  Modified on: Jan. 22, 2010
 *  Description: Remove the scaling of the x and y when doing the polynomial fitting
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications       
 ****************************************************/
#include "MATH_polynomial.h"

/*------------------------------------------------------------
 * the functions for polynomial fitting
 *------------------------------------------------------------*/

/**
 * P0 + P1x + P2x^2 + .. = y, x and y vectors are know, get the coefficients of P0, P1, P2...
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffer
 *    2 - illegal numbers
 *    3 - other failures
 */
int MATH_poly_fit(const double *X0,
                  const double *Y0,
                  int pointNum,
                  int polyOrder,
                  double *P)
{   
   int status;
   int i, j;
   int rowNum, colNum;                                /* row and column numbers for the coefficient matrix */
   double tmpData;	
   MATH_DATA_MATRIX A;                                /* the matrix for x^j, j = 0..polyOrder */  
   
   /*
   printf("point num = %d\n", pointNum);
   printf("poly order = %d\n", polyOrder);
   */
   
   /* check the input parameters */
   if(!X0 || !Y0 || !P)                                        return 1;
   if(pointNum <= 0 || polyOrder <=0 || polyOrder >= pointNum) return 2;   
   
   /* init the variables */
   rowNum = pointNum;
   colNum = polyOrder + 1;

   /* calculate the coefficient matrix */
   MATH_matrix_init(&A);                             /* to avoid warning from the matrix create function!!-Modify the  MATH_matrixCreate_NewBuf function !*/
   MATH_matrix_create_new(&A, rowNum, colNum);       /* create a new matrix, do not forget to delete it */
   
   for(i = 0; i < rowNum; i ++)
   {
      tmpData = 1;
      MATH_matrix_set_element(&A, i, 0, tmpData);    /* the first column is 1 */

      for(j = 1; j < colNum; j ++)                   /* other column is x^j */
      {
         tmpData *= *(X0 + i);
         MATH_matrix_set_element(&A, i, j, tmpData);       
      }         
   }
   
   /* solve the linear equations */
   MATH_linear_equations_solver(&A, (double *)Y0, P);
   
   /* free the buffer */ 
   status = MATH_matrix_delete_new(&A);             /* delete the newly designed matrix */
   
   if(!status) return 0;
   else        return 3;
} 

/**
 * calculate the function value of the polinomial at X, the output is a vector (data group) 
 * for the derOrder, 0 means get the polynomial values, 1 means get the first order derivative,
 * currently, only 0 and 1 are valid
 *
 * Reference: book of <<Numerical Recipes - The art of scientific computing>>, 3rd edition, P202
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffer
 *    2 - illegal numbers
 */                  
int MATH_poly_val_vec(const double *X,
                  const double *P,
                  int pointNum,
                  int polyOrder,
                  int derOrder,
                  double *valY)
{
   int i, j;
   double tmpVal, tmpDer;

   /* check the input parameters */
   if(!X || !P || !valY)                                              return 1;
   if(pointNum <= 0 || polyOrder <=0 || derOrder < 0 || derOrder > 1) return 2;  
   
   /* calculate the values */
   for(i = 0; i < pointNum; i ++)
   {
      /* get the polynomial value and derivative for the ith point */
      tmpVal = *(P + polyOrder);
      tmpDer = 0.0;
      
      for(j = polyOrder - 1; j >= 0; j --)
      {
         tmpDer = tmpDer * (*(X + i)) + tmpVal;
         tmpVal = tmpVal * (*(X + i)) + *(P + j);     
      }
      
      /* output based on the selection */
      if(derOrder == 0)
         *(valY + i) = tmpVal;
      else
         *(valY + i) = tmpDer;
   }
   
   return 0;
}          

/**
 * calculate the function value of the polinomial at x, the output is a single value 
 * for the derOrder, 0 means get the polynomial values, 1 means get the first order derivative,
 * currently, only 0 and 1 are valid
 *
 * Reference: book of <<Numerical Recipes - The art of scientific computing>>, 3rd edition, P202
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffer
 *    2 - illegal numbers
 */                          
double MATH_poly_val(double x,
                  const double *P,
                  int polyOrder,
                  int derOrder)                         
{  
   int j;
   double tmpVal, tmpDer;

   /* check the input parameters */
   if(!P)                                            return 1;
   if(polyOrder <=0 || derOrder < 0 || derOrder > 1) return 2;  
   
   /* get the polynomial value and derivative for the point */
   tmpVal = *(P + polyOrder);
   tmpDer = 0.0;
      
   for(j = polyOrder - 1; j >= 0; j --)
   {
      tmpDer = tmpDer * x + tmpVal;
      tmpVal = tmpVal * x + *(P + j);     
   }
      
   /* output based on the selection */
   if(derOrder == 0)
      return tmpVal;
   else
      return tmpDer;  
      
   return 0;
}                                
                                                  
/*------------------------------------------------------------
 * the functions for polynomial interpolation
 *------------------------------------------------------------*/                          

/**
 * interpolation of the polynomial value by fitting to the interpIndex,
 * the output will be a vector (data group)
 * for the derOrder, 0 means get the polynomial values, 1 means get the first order derivative,
 * currently, only 0 and 1 are valid
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffer
 *    2 - illegal numbers
 *    3 - failure
 */
int MATH_poly_interp_val_vec(const double *sourceIndex,
                  const double *sourceData,
                  int sourcePointNum,
                  int sourcePolyOrder,
                  int sourceDerOrder,
                  int interpPointNum,
                  const double *interpIndex,
                  double *interpVal)
{
   int status;
   double *P;                       /* the polynomial coefficients */
   
   /* check the input parameters */
   if(!sourceIndex || !sourceData || !interpIndex || !interpVal) return 1;
   if(sourcePointNum <=0 || sourcePolyOrder <=0 || interpPointNum <=0 || sourcePolyOrder > sourcePointNum || sourceDerOrder < 0 || sourceDerOrder > 1) return 2;  
   
   /* init the variables */
   status = MATH_buffer_create(&P, sourcePolyOrder + 1);
   if(status != 0) return 3;

   /* poly fitting */
   status = MATH_poly_fit(sourceIndex, sourceData, sourcePointNum, sourcePolyOrder, P);
   
   /* get the poly values */  
   status = MATH_poly_val_vec(interpIndex, P, interpPointNum, sourcePolyOrder, sourceDerOrder, interpVal);

   /* free the memory */
   status = MATH_buffer_delete(&P);
   
   if(!status) return 0;
   else        return 3;
} 

/**
 * interpolation of the polynomial value by fitting to the interpIndex,
 * the output will be a single value
 * for the derOrder, 0 means get the polynomial values, 1 means get the first order derivative,
 * currently, only 0 and 1 are valid
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffer
 *    2 - illegal numbers
 *    3 - failure
 */
double MATH_poly_interp_val(const double *sourceIndex,
                  const double *sourceData,
                  int sourcePointNum,
                  int sourcePolyOrder,                                         
                  int sourceDerOrder,
                  double interpIndex)
{
   int status;
   double *P;                       /* the polynomial coefficients */
   double tempVal;

   /* check the input parameters */
   if(!sourceIndex || !sourceData) return 1;
   if(sourcePointNum <=0 || sourcePolyOrder <=0 || sourcePolyOrder > sourcePointNum || sourceDerOrder < 0 || sourceDerOrder > 1) return 2;  
    
   /* init the variables */
   status = MATH_buffer_create(&P, sourcePolyOrder + 1);
   if(status != 0) return 3;
   
   /* poly fitting */
   status = MATH_poly_fit(sourceIndex, sourceData, sourcePointNum, sourcePolyOrder, P);
   
   /*
   printf("status of polyfit = %d\n", status);   
   printf("P1 = %f\n", *(P));
   printf("P2 = %f\n", *(P+1));
   */   
      
   /* get the poly values */
   tempVal = MATH_poly_val(interpIndex, P, sourcePolyOrder, sourceDerOrder);

   /* free the memory */
   status = MATH_buffer_delete(&P);
   
   return tempVal;
}


