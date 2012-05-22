/****************************************************
 *  MATH_linear_equations.c                                         
 *                                      
 *  Realize the functions to solve linear equations
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: April 26, 2009      
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: The first development period until May 10, 2009
 *  Description: Solve the linear equation          
 *  To be done:  add exception handling and run time message 
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 04, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library     
 *  To be done:  add exception handling and run time message
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications              
 ****************************************************/
#include "MATH_linear_equations.h"

/**
 * Solve the linear equation A0*x = b0
 * x = A0 \ b0
 * Reference: book of <<Numerical Recipes - The art of scientific computing>>, 3rd edition, P53
 *
 * To be done: check the singularity of the matrix, and if the matrix is singular, generate alarms
 */
int MATH_linear_equations_solver(MATH_DATA_MATRIX *A0, double *b0, double *x)
{                  
   int i, ii = 0, ip, j, m, n, *permu;                                     
   double sum, tmpData1;
   
   MATH_DATA_MATRIX A0T, A0TA0, matrixb, A0Tb, *A;              /* the variables for matrix calculation */
   double *b;
   
   /* init variables */
   m     = A0 -> rowNum;
   n     = A0 -> colNum;
   permu = (int *)calloc(n, sizeof(int));
   
   /* normalize the matrix of A0 if the matrix is not square */
   if(m != n)
   {
      MATH_matrix_create_new(&A0TA0, n, n);                  /* the new matrix for A' * A */
      MATH_matrix_create_new(&A0Tb, n, 1);                   /* the new matrix for A' * b */
      MATH_matrix_create_ext(&matrixb, m, 1, b0);            /* form a matrix of vector b0 for matrix calculation */
      
      MATH_matrix_transpose(A0, &A0T);                       /* the transpose uses the same buffer as input matrix */
      MATH_matrix_mul(&A0T, A0, &A0TA0);                     /* calculate A' * A */
      MATH_matrix_mul(&A0T, &matrixb, &A0Tb);                /* calculate A' * b */ 
      
      A = &A0TA0;
      b = A0Tb.data;
   }
   else
   {
      A = A0;
      b = b0;
   }
   
   /* get the LU decomposition of A */
   MATH_matrix_LU_decomp(A, permu);
   
   /* solve the equation */
   for(i = 0; i < n; i ++)
      *(x + i) = *(b + i);
      
   for(i = 0; i < n; i ++)
   {
      ip        = *(permu + i);
      sum       = *(x + ip);
      *(x + ip) = *(x + i);  
      
      if(ii != 0)
      {
         for(j = ii - 1; j < i; j ++)
            sum -= MATH_matrix_get_element(A, i, j) * (*(x + j));           
      }
      else if(sum != 0.0)
         ii = i + 1;
         
      *(x + i) = sum;    
   }
   
   for(i = n - 1; i >= 0; i --)
   {
      sum = *(x + i);
      for(j = i + 1; j < n; j ++)
         sum -= MATH_matrix_get_element(A, i, j) * (*(x + j));     
      
      *(x + i) = sum / MATH_matrix_get_element(A, i, i);
   }

   /* free the memory */
   free(permu);
   permu = 0;
   
   if(m != n)
   {
      MATH_matrix_delete_new(&A0TA0); 
      MATH_matrix_delete_new(&A0Tb);
   }
   
   return 0;
}


