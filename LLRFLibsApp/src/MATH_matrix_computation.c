/****************************************************
 *  MATH_matrix_computation.c                                         
 * 
 *  Realize the functions to perform matrix calculation
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: April 26, 2009                     
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: The first development period until May 10, 2009
 *  Description: Define the basic operation of the complex vector:
 *               1. create and destroy
 *               2. get and set element, matrix copy
 *               3. transpose
 *               4. multiply
 *               5. row swap
 *               6. LU decomposition
 *  To be done : add exception handling and run time message  
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Dec. 04, 2009
 *  Description: Reorganize the code for the LLRF Algorithm C Library  
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Jan. 28, 2010
 *  Description: Change the matrix element read/write function to inline function
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications   
 ****************************************************/
#include "MATH_matrix_computation.h"

/*------------------------------------------------------------
 * the functions to create and destroy the matrix
 *------------------------------------------------------------*/
 
/**
 *  create a matrix, with newly created buffer
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal matrix
 *    2 - illegal row number or column number
 *    3 - failed to create buffer
 * 
 * Note:
 *    1. To be added: first check if the matrix is empty or still occupied
 */
int MATH_matrix_create_new(MATH_DATA_MATRIX *matrix, int rowNum, int colNum)
{
   int status;
   double *data_Buf;

   /* check the input parameters */
   if(!matrix)                    return 1;
   if(rowNum <= 0 || colNum <= 0) return 2; 
   
   /* create buffer */
   status = MATH_buffer_create(&data_Buf, rowNum * colNum);
   
   if(status != 0)    return 3;
   
   /* init the matrix */
   matrix -> rowNum   = rowNum;
   matrix -> colNum   = colNum;
   matrix -> majorDim = MAJOR_COL;         
   matrix -> data     = data_Buf;      
   
   return 0;                    
}

/**
 * delete a matrix, only works on the matrix newly created buffer
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal complex vector
 *    2 - failed to delete buffer
 */
int MATH_matrix_delete_new(MATH_DATA_MATRIX *matrix)
{
   int status;
   
   /* check the input parameters */
   if(!matrix)   return 1;

   /* delete the buffer */
   matrix -> rowNum = 0;
   matrix -> colNum = 0;   
   status = MATH_buffer_delete(&(matrix -> data));

   if(!status) return 0;
   else        return 2;
}

/**
 * create a matrix, use the external buffer, NEVER delete it
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffer
 *    2 - illegal row number or column number
 */
int MATH_matrix_create_ext(MATH_DATA_MATRIX *matrix, int rowNum, int colNum, double *extBuf)
{
   /* check the input parameters */
   if(!matrix || !extBuf)         return 1;
   if(rowNum <= 0 || colNum <= 0) return 2; 
   
   /* init the matrix */
   matrix -> rowNum   = rowNum;
   matrix -> colNum   = colNum;
   matrix -> majorDim = MAJOR_COL;
   matrix -> data     = extBuf;  
   
   return 0;                             
}

/**
 * init the matrix to empty state
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal matrix
 */
int MATH_matrix_init(MATH_DATA_MATRIX *matrix)
{
   /* check the input parameters */
   if(!matrix) return 1;
   
   matrix -> rowNum = 0;
   matrix -> colNum = 0;
   matrix -> data   = 0;
   
   return 0;    
}

/*------------------------------------------------------------
 * the functions for matrix data access
 *------------------------------------------------------------*/
/**
 * get an element from the matrix, rowId shoud be within [0, rowNum-1], colId should be within [0, colNum]
 */
double MATH_matrix_get_element(MATH_DATA_MATRIX *matrixIn, int rowId, int colId)
{
   return (matrixIn -> majorDim == MAJOR_COL)?*(matrixIn -> data + colId * matrixIn -> rowNum + rowId):*(matrixIn -> data + rowId * matrixIn -> colNum + colId);
}

/**
 * set an element to the matrix
 */
void MATH_matrix_set_element(MATH_DATA_MATRIX *matrixIn, int rowId, int colId, double elemData)
{
   (matrixIn -> majorDim == MAJOR_COL)?(*(matrixIn -> data + colId * matrixIn -> rowNum + rowId) = elemData):(*(matrixIn -> data + rowId * matrixIn -> colNum + colId) = elemData);
} 

/**
 * init a matrix by copying an existing one, they have to have the same dimensions
 * the major dimension settings will also be copyed to the destination matrix
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal matrix
 *    2 - illegal dimensions
 *    3 - copy failure
 */
int MATH_matrix_copy(MATH_DATA_MATRIX *matrixDest, MATH_DATA_MATRIX *matrixSrc)
{
   int status;
   
   /* check the input parameters */
   if(!matrixDest || !matrixDest -> data || !matrixSrc || !matrixSrc -> data)                      return 1;
   if(matrixDest -> rowNum != matrixSrc -> rowNum ||  matrixDest -> colNum != matrixSrc -> colNum) return 2;   
                             
   /* copy the matrix data buffer, and copy the major dimension settings */
   matrixDest -> majorDim = matrixSrc -> majorDim;
   status = MATH_buffer_copy(matrixDest -> data, matrixSrc -> data, matrixDest -> rowNum * matrixDest -> colNum);
   
   if(!status) return 0;
   else        return 3;
}

/*------------------------------------------------------------
 * the functions for matrix computation
 *------------------------------------------------------------*/

/**
 * transpose a matrix, two cases
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal matrix
 */
int MATH_matrix_transpose(MATH_DATA_MATRIX *matrixIn, MATH_DATA_MATRIX *matrixOut)
{
   int tmpSwap;

   /* check the input parameters */
   if(!matrixIn || !matrixIn -> data) return 1;
   
   /* if the matrixOut is set to zero, replace the input matrix */
   if(!matrixOut)
   {
      /* first case, replace the initial matrix with the transposed one, only change the configuration (this is why we provide the major dimension settings!) */
      tmpSwap = matrixIn -> rowNum;
   
      matrixIn -> rowNum   = matrixIn -> colNum;
      matrixIn -> colNum   = tmpSwap;	
      matrixIn -> majorDim = (MATH_ENUM_MATRIX_MDIM)((matrixIn -> majorDim) ^ 1);
   }
   else
   {
      /* second case, build the matrixOut, but share the real data buffer with matrixIn */
      /* 此时matrixOut必须是空的，否则将失去其数据的指针，在free的时候将出错 */
      /* 注意：在所有可能改变非空数据buffer指针的时候，都要注意这一点. 比如对某指针初始化的时候使用了动态内存分配函数，在释放这部分内存之前
              如果把这个指针指向别的buffer了，那么释放这个指针的时候就相当于释放了错误的buffer! */
      /* 解决方案：在矩阵声明的时候，初始化rowNum，colNum和data指针都为0，在后面对其指针进行赋值时，检查指针的值，
              如果不是0，就有上面问题的可能，需要提供alarm或warning, see the function of matrixInit*/ 
       
      matrixOut -> rowNum   = matrixIn -> colNum;
      matrixOut -> colNum   = matrixIn -> rowNum;
      matrixOut -> majorDim = (MATH_ENUM_MATRIX_MDIM)((matrixIn -> majorDim) ^ 1);
      matrixOut -> data     = matrixIn -> data;
   } 
   
   return 0;  
} 

/**
 * add two input matrix, and store the result in matrixOut
 * the dimension of the two matrices have to be the same
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal matrix
 *    2 - illegal dimensions
 */
int MATH_matrix_add(MATH_DATA_MATRIX *matrix1, MATH_DATA_MATRIX *matrix2, MATH_DATA_MATRIX *matrixOut)
{
   int i, j;
   
   /* check the input parameters */
   if(!matrix1 || !matrix1 -> data || !matrix2 || !matrix2 -> data || !matrixOut || !matrixOut -> data) return 1;
   if(matrix1 -> rowNum != matrix2 -> rowNum || matrix1 -> colNum != matrix2 -> colNum) return 2;
      
   /* init the configuration of the output matrix */
   matrixOut -> rowNum   = matrix1 -> rowNum;
   matrixOut -> colNum   = matrix2 -> colNum;
   matrixOut -> majorDim = matrix1 -> majorDim;   
   
   /* calculate the matrix elements */
   for(i = 0; i < matrixOut -> rowNum; i ++)
      for(j = 0; j < matrixOut -> colNum; j ++)
          *(matrixOut -> data + MATH_MACRO_MATRIX_DATA_ID(matrixOut,i,j)) = *(matrix1 -> data + MATH_MACRO_MATRIX_DATA_ID(matrix1,i,j)) + 
                                                                            *(matrix2 -> data + MATH_MACRO_MATRIX_DATA_ID(matrix2,i,j));
   
   return 0;
} 

/**
 * multiply two input matrix, and store the result in matrixOut
 * the dimension of the two matrices have to be consistent for multiplication
 * n*m and m*k, the matrix output will be n*m
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal matrix
 *    2 - illegal dimensions
 */
int MATH_matrix_mul(MATH_DATA_MATRIX *matrix1, MATH_DATA_MATRIX *matrix2, MATH_DATA_MATRIX *matrixOut)
{
   int i, j, k;
   double tmpAccum;
   
   /* check the input parameters */
   if(!matrix1 || !matrix1 -> data || !matrix2 || !matrix2 -> data || !matrixOut || !matrixOut -> data) return 1;
   if(matrix1 -> colNum != matrix2 -> rowNum) return 2;
      
   /* init the configuration of the output matrix */
   matrixOut -> rowNum   = matrix1 -> rowNum;
   matrixOut -> colNum   = matrix2 -> colNum;
   matrixOut -> majorDim = matrix1 -> majorDim;     
   
   /* calculate the matrix elements */
   for(i = 0; i < matrixOut -> rowNum; i ++)
   {
      for(j = 0; j < matrixOut -> colNum; j ++)
      {
         tmpAccum = 0;
         for(k = 0; k < matrix1 -> colNum; k ++)
            /*tmpAccum += MATH_matrix_get_element(matrix1, i, k) * MATH_matrix_get_element(matrix2, k, j);*/
            tmpAccum += *(matrix1->data + MATH_MACRO_MATRIX_DATA_ID(matrix1,i,k)) * (*(matrix2->data + MATH_MACRO_MATRIX_DATA_ID(matrix2,k,j)));
         
         /*MATH_matrix_set_element(matrixOut, i, j, tmpAccum);*/
         *(matrixOut -> data + MATH_MACRO_MATRIX_DATA_ID(matrixOut,i,j)) = tmpAccum;
      }      
   }	
   
   return 0;
} 

/**
 * Swap the rows of the matrix, mainly used for LU decomposition
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal matrix
 *    2 - illegal dimensions
 *    3 - buffer error
 */
int MATH_matrix_swap_row(MATH_DATA_MATRIX *matrixIn, int rowId1, int rowId2)
{
   int status; 
   int j;
   double tmpData1;
   double tmpData2;
   double *tmpRowBuf;

   /* check the input parameters */
   if(!matrixIn || !matrixIn -> data)                                                           return 1;
   if(rowId1 < 0 || rowId1 >= matrixIn -> rowNum || rowId2 < 0 || rowId2 >= matrixIn -> rowNum) return 2;
   
   /* exchange the rows */
   if(matrixIn -> majorDim == MAJOR_COL)                                       /* the storage is by column, have to swap element by element */
   { 
      for(j = 0; j < matrixIn -> colNum; j ++)
      {
         tmpData1 = MATH_matrix_get_element(matrixIn, rowId1, j);
         tmpData2 = MATH_matrix_get_element(matrixIn, rowId2, j);
      
         MATH_matrix_set_element(matrixIn, rowId1, j, tmpData2);
         MATH_matrix_set_element(matrixIn, rowId2, j, tmpData1);
      }
   }
   else                                                                       /* the storage is by row, block copy can be used */
   {
      status = MATH_buffer_create(&tmpRowBuf, matrixIn -> colNum);   
      if(status != 0) return 3;
      status = MATH_buffer_copy(tmpRowBuf, matrixIn -> data + rowId1 * matrixIn -> colNum, matrixIn -> colNum);
      if(status != 0) return 3;
      status = MATH_buffer_copy(matrixIn -> data + rowId1 * matrixIn -> colNum, matrixIn -> data + rowId2 * matrixIn -> colNum, matrixIn -> colNum);
      if(status != 0) return 3;
      status = MATH_buffer_copy(matrixIn -> data + rowId2 * matrixIn -> colNum, tmpRowBuf, matrixIn -> colNum);
      if(status != 0) return 3;
      status = MATH_buffer_delete(&tmpRowBuf);                                                                    
      if(status != 0) return 3;
   }
   
   return 0;
}

/**
 * LU decomposition of the matrix
 * Reference: book of <<Numerical Recipes - The art of scientific computing>>, 3rd edition, P52
 *            Please note that there is problem of the code on the text book, use the code provided
 *            by the book CD
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal matrix
 *    2 - buffer error
 *    3 - matrix singular
 */
int MATH_matrix_LU_decomp(MATH_DATA_MATRIX *A, int *permutation)
{
    int status;
	int n = A -> colNum;
    double d;
    
    const double TINY = 1.0e-20;             /* a small number */
	int i, imax, j, k;
	double big, dum, sum, temp, tmpData1;
    double *vv;                              /* vv stors the implicit scaling of each row */

    /* check the input */
    if(!A || !A -> data || !permutation) return 1;

    /* create buffer */
    status = MATH_buffer_create(&vv, n); 
    if(status != 0) return 2;

	d = 1.0;
	for(i = 0; i < n; i ++) 
    {
	   	big = 0.0;
		for(j = 0; j < n; j ++)
		{
            if((temp = fabs(MATH_matrix_get_element(A, i, j))) > big)
               big = temp;
        }
		if(big == 0.0) 
		{
           /*printf("Singular matrix in routine ludcmp\n");*/
           return 3;
        }   
		vv[i] = 1.0 / big;
	}
	
	for(j = 0; j < n; j ++) 
    {
		for(i = 0; i < j; i ++) 
        {
            sum = MATH_matrix_get_element(A, i, j);
			for(k = 0; k < i; k ++)
               sum -= MATH_matrix_get_element(A, i, k) * MATH_matrix_get_element(A, k, j);

			MATH_matrix_set_element(A, i, j, sum);
		}
		
		big = 0.0;
		for(i = j; i < n; i ++) 
        {
			sum = MATH_matrix_get_element(A, i, j);
			for(k = 0; k < j; k ++)
               sum -= MATH_matrix_get_element(A, i, k) * MATH_matrix_get_element(A, k, j);
            MATH_matrix_set_element(A, i, j, sum);
            
			if((dum = vv[i] * fabs(sum)) >= big) 
            {
				big = dum;
				imax = i;
			}
		}
		
		if(j != imax) 
        {
            MATH_matrix_swap_row(A, imax, j);
			d = -d;
			vv[imax] = vv[j];
		}
		
		permutation[j] = imax;
		
		tmpData1 = MATH_matrix_get_element(A, j, j);
		if(tmpData1 == 0.0) 
           MATH_matrix_set_element(A, j, j, TINY);
           
		if(j != n-1) 
        {
			dum = 1.0 / MATH_matrix_get_element(A, j, j);
			for(i = j + 1; i < n; i ++)
            {
               tmpData1 = MATH_matrix_get_element(A, i, j);
               tmpData1 *= dum;
               MATH_matrix_set_element(A, i, j, tmpData1);
            }
		}
	}   
 
   MATH_buffer_delete(&vv);                   
}

/*------------------------------------------------------------
 * the functions for debugging
 *------------------------------------------------------------*/

/**
 * print the matrix
 */
void MATH_matrix_print(MATH_DATA_MATRIX *matrixIn)
{
   int i, j;
   
   printf("The matrix row and column numbers are: %d, %d\n", matrixIn -> rowNum, matrixIn -> colNum);
   
   for(i = 0; i < matrixIn -> rowNum; i ++)
   {
      for(j = 0; j < matrixIn -> colNum; j ++)
         printf("%f ", MATH_matrix_get_element(matrixIn, i, j));
      printf("\n");
   }
   
   printf("\n");
} 

