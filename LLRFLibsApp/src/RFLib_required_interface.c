/****************************************************
 *  RFLib_required_interface.c                                         
 * 
 *  Realize of the required interface for the library of LLRF Applications
 * 
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: Dec. 05, 2009       
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: June 07, 2010
 *  Description: Optimize the circle fitting, remove the points with large errors   
 *
 *  Modified by: Zheqiao Geng, gengzq@slac.stanford.edu
 *  Modified on: 2011.11.26
 *  Description: Immegrate to SLAC LLRF applications  
 * 
 *  Modified by: Zheqiao Geng
 *  Modified on: 2011.12.06
 *  Description: Add the function to do cos fitting
 ****************************************************/
#include "RFLib_required_interface.h"

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
int RFLib_common_buffer_create(double **buf, int pointNum)
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
int RFLib_common_buffer_delete(double **buf)
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
int RFLib_common_buffer_copy(double *destBuf, double *srcBuf, int pointNum)
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
int RFLib_common_buffer_init_linear(double *destBuf, double initValue, double stepValue, int pointNum)
{  
    int status;
    
    /* check the input parameters */
    if(!destBuf)       return 1;
    if(pointNum <= 0)  return 2;    
                            
    status = MEMM_buffer_init(destBuf, initValue, stepValue, 0, pointNum, 1);                                                                                     
    
    if(!status) return 0;
    else        return 3;     
}                            

/**
 * shift the buffer, negative number is for left, positive number for right
 * after shifting, the empty buffer is filled by zero
 *
 * Return value:
 *    0 - succeed
 *    1 - illegal buffers
 *    2 - illegal point number
 *    3 - failed to perform data shift
 *
 */
int RFLib_common_buffer_shift(double *buf, int point_num, int shift_num)
{
    return DSP_data_shift(buf, point_num, shift_num);
}

/*------------------------------------------------------------
 * the functions for mathematic calculation
 *------------------------------------------------------------*/
/**
 * amplitude of the single complex value
 */
double RFLib_common_complex_abs(double data_I, double data_Q)
{
    return sqrt(data_I * data_I + data_Q * data_Q);           
}

/**
 * amplitude of the complex vector, long
 */
int RFLib_common_complex_abs_vector_long(long *buf_I, long *buf_Q, 
                                             double *buf_A, int pointNum)
{
    int i;
    double data_I, data_Q;
            
    for(i = 0; i < pointNum; i ++) 
    {
        data_I       = (float)(*(buf_I + i));
        data_Q       = (float)(*(buf_Q + i));
        *(buf_A + i) = sqrt(data_I * data_I + data_Q * data_Q);
    }
    
    return 0;	
}                                        

/**
 * amplitude of the complex vector, float
 */
int RFLib_common_complex_abs_vector_float(double *buf_I, double *buf_Q, 
                                              double *buf_A, int pointNum)
{    
    int i;
    
    for(i = 0; i < pointNum; i ++)
    {
        *(buf_A + i) = sqrt(*(buf_I + i) * (*(buf_I + i)) + *(buf_Q + i) * (*(buf_Q + i)));	
    }                                              	
    
    return 0;
}


/**
 * phase of the single complex value, rad
 */
double RFLib_common_complex_angle(double data_I, double data_Q)
{
    return MATH_norm_phase(data_I, data_Q);           
}

/**
 * get the amplitude and phase (rad) from the I and Q, save to the same buffers
 */
int RFLib_common_complex_polar(double *buf_I, double *buf_Q, int pointNum)
{
    int status;
    RFLIB_DATA_COMPLEX c_data;
    
    /* create complex vectors */
    status = MATH_complex_create_ext(&c_data, buf_I, buf_Q, pointNum, MODE_IQ);

    /* ploar */
    status = MATH_complex_polar(&c_data);
    
    return status;        
}

/**
 * get the amplitude and phase (rad) from the I and Q, save to different buffers
 */
int RFLib_common_complex_polar_d(double *buf_I, double *buf_Q, 
                                     double *buf_A, double *buf_P, int pointNum)
{
    int status;
    RFLIB_DATA_COMPLEX c_data;
    
    /* first copy the buffers to the destination */
    status = MEMM_buffer_init(buf_A, 0, 0, buf_I, pointNum, 2);
    status = MEMM_buffer_init(buf_P, 0, 0, buf_Q, pointNum, 2);	
    
    /* create complex vectors */
    status = MATH_complex_create_ext(&c_data, buf_A, buf_P, pointNum, MODE_IQ);

    /* ploar */
    status = MATH_complex_polar(&c_data);
    
    return status;    
}                                     

/**
 * get the I and Q from amplitude and phase (rad), save to the same buffer
 */
int RFLib_common_complex_orthogonal(double *buf_A, double *buf_P, int pointNum)
{
    int status;
    RFLIB_DATA_COMPLEX c_data;
    
    /* create complex vectors */
    status = MATH_complex_create_ext(&c_data, buf_A, buf_P, pointNum, MODE_AP);

    /* ploar */
    status = MATH_complex_orthogonal(&c_data);
    
    return status; 
}

/**
 * complex addition
 * c1 = complex(buf_I1, buf_Q1);
 * c2 = complex(buf_I2, buf_Q2);
 * c3 = complex(buf_I3, buf_Q3);
 * c3 = c1 + c2;
 */
int RFLib_common_complex_add(double *buf_I1, double *buf_Q1,      
                                 double *buf_I2, double *buf_Q2,
                                 double *buf_I3, double *buf_Q3, 
                                 int pointNum)
{
    int status;
    RFLIB_DATA_COMPLEX c1, c2, c3;
    
    /* create complex vectors */
    status = MATH_complex_create_ext(&c1, buf_I1, buf_Q1, pointNum, MODE_IQ);
    status = MATH_complex_create_ext(&c2, buf_I2, buf_Q2, pointNum, MODE_IQ);
    status = MATH_complex_create_ext(&c3, buf_I3, buf_Q3, pointNum, MODE_IQ);
    
    /* division */        
    status = MATH_complex_add(&c1, &c2, &c3);
    
    return status;                
} 

/**
 * complex division
 * c1 = complex(buf_I1, buf_Q1);
 * c2 = complex(buf_I2, buf_Q2);
 * c3 = complex(buf_I3, buf_Q3);
 * c3 = c1 ./ c2;
 */
int RFLib_common_complex_div(double *buf_I1, double *buf_Q1,      
                                 double *buf_I2, double *buf_Q2,
                                 double *buf_I3, double *buf_Q3, 
                                 int pointNum)
{
    int status;
    RFLIB_DATA_COMPLEX c1, c2, c3;
    
    /* create complex vectors */
    status = MATH_complex_create_ext(&c1, buf_I1, buf_Q1, pointNum, MODE_IQ);
    status = MATH_complex_create_ext(&c2, buf_I2, buf_Q2, pointNum, MODE_IQ);
    status = MATH_complex_create_ext(&c3, buf_I3, buf_Q3, pointNum, MODE_IQ);
    
    /* division */        
    status = MATH_complex_div(&c1, &c2, &c3);
    
    return status;                
}                            

/**
 * complex multiplification by a constant
 * c1 = complex(buf_I1, buf_Q1);
 * c2 = complex(data_I2, data_Q2); is a constant
 * c3 = complex(buf_I3, buf_Q3);
 * c3 = c1 * c2;
 */
int RFLib_common_complex_mul_const(double *buf_I1, double *buf_Q1,                        
                                       double data_I2, double data_Q2,
                                       double *buf_I3, double *buf_Q3, 
                                       int pointNum)
{
    int status;
    RFLIB_DATA_COMPLEX c1, c2, c3;
    
    /* create complex vectors */
    status = MATH_complex_create_ext(&c1, buf_I1, buf_Q1, pointNum, MODE_IQ);
    status = MATH_complex_create_ext(&c2, &data_I2, &data_Q2, 1, MODE_IQ);
    status = MATH_complex_create_ext(&c3, buf_I3, buf_Q3, pointNum, MODE_IQ);
    
    /* multiplification */        
    status = MATH_complex_mul(&c1, &c2, &c3);
    
    return status;                                            
}

/**
 * complex division by a constant
 * c1 = complex(buf_I1, buf_Q1);
 * c2 = complex(data_I2, data_Q2); is a constant
 * c3 = complex(buf_I3, buf_Q3);
 * c3 = c1 / c2;
 */
int RFLib_common_complex_div_const(double *buf_I1, double *buf_Q1,                        
                                       double data_I2, double data_Q2,
                                       double *buf_I3, double *buf_Q3, 
                                       int pointNum)
{
    int status;
    RFLIB_DATA_COMPLEX c1, c2, c3;
    
    /* create complex vectors */
    status = MATH_complex_create_ext(&c1, buf_I1, buf_Q1, pointNum, MODE_IQ);
    status = MATH_complex_create_ext(&c2, &data_I2, &data_Q2, 1, MODE_IQ);
    status = MATH_complex_create_ext(&c3, buf_I3, buf_Q3, pointNum, MODE_IQ);
    
    /* multiplification */        
    status = MATH_complex_div(&c1, &c2, &c3);
    
    return status;                                            
}

/**
 * Measure the vector sub: vector1 - vector2*scale, results are stored in vector1
 * the second one should not be changed (here realize it from start, because if we
 * call the function in MATH, the second will also be changed)
 */
int RFLib_common_complex_sub_scale(double *buf_I1, double *buf_Q1,
                                       double *buf_I2, double *buf_Q2,
                                       double data_scale, int pointNum)
{
    int i;
   
    for(i = 0; i < pointNum; i ++)
    {
        *(buf_I1 + i) -= *(buf_I2 + i) * data_scale;
        *(buf_Q1 + i) -= *(buf_Q2 + i) * data_scale;       
    }   
    
    return 0;                                       
}
 
/**
 * Create a new matrix
 */
int RFLib_common_matrix_create(RFLIB_DATA_MATRIX *matrix, int rowNum, int colNum)
{
    int status;
   
    status = MATH_matrix_init(matrix); 
    status = MATH_matrix_create_new(matrix, rowNum, colNum);   
   
    return status;
}
 
/**
 * Delete a new created matrix
 */ 
int RFLib_common_matrix_delete(RFLIB_DATA_MATRIX *matrix)
{
    int status;
    
    status = MATH_matrix_delete_new(matrix);    
    
    return status;
}

/**
 * create a matrix with existing buffer
 */
int RFLib_common_matrix_create_ext(RFLIB_DATA_MATRIX *matrix, int rowNum, int colNum, double *extBuf)
{
    int status;
   
    status = MATH_matrix_init(matrix); 
    status = MATH_matrix_create_ext(matrix, rowNum, colNum, extBuf);   
   
    return status;    	
}

/**
 * Set element to the matrix
 */
int RFLib_common_matrix_set_element(MATH_DATA_MATRIX *matrixIn, int rowId, int colId, double elemData)
{
   int status;
   
   MATH_matrix_set_element(matrixIn, rowId, colId, elemData); 
   
   return status;   
}

/**
 * matrix add
 */
int RFLib_common_matrix_add(RFLIB_DATA_MATRIX *matrix1, RFLIB_DATA_MATRIX *matrix2, RFLIB_DATA_MATRIX *matrixOut)
{
    return MATH_matrix_add(matrix1, matrix2, matrixOut);	
}

/**
 * matrix multiply
 */
int RFLib_common_matrix_mul(RFLIB_DATA_MATRIX *matrix1, RFLIB_DATA_MATRIX *matrix2, RFLIB_DATA_MATRIX *matrixOut)
{
    return MATH_matrix_mul(matrix1, matrix2, matrixOut);		
}


/**
 * Function to measure the average of the buffer
 * TBD: exception handlings
 */
double RFLib_common_mea_avg(const double *dataIn, int pointNum)
{
    return MATH_statistic_avg(dataIn, pointNum);           
}

/**
 * Function to perform linear interpolation calculation
 * 1. The input data is fitted to a linear function, with x as the point number, and the
 *    function will return the value of the fitted linear function at the defined x (dest_x)
 * 2. TBD: exception handlings
 */
double RFLib_common_linear_interp(const double *dataIn, int pointNum, double dest_x)
{
    int status;
    double tmp_result;                                 /* store the result */
    double *buf_x;                                     /* the x buffer for fitting */
    
    /* creat buffer and init it for x */
    status = RFLib_common_buffer_create(&buf_x, pointNum);
    status = RFLib_common_buffer_init_linear(buf_x, 0, 1, pointNum);
    
    /* calculate the interpolation */
    tmp_result = MATH_poly_interp_val(buf_x, dataIn, pointNum, 1, 0, dest_x);
    
    /* delete the buffer */
    status = RFLib_common_buffer_delete(&buf_x);
    
    return tmp_result;
}

/**
 * Function for multi-variable linear fitting
 * Solve the problem of Ax = b, so x = A\b
 */
int RFLib_common_linear_fit(RFLIB_DATA_MATRIX *A, double *b, double *x)
{
    int status;
    
    status = MATH_linear_equations_solver(A, b, x);    
    
    return status;
}

/**
 * Function to perform polynomial fitting, the x data is the unit vector (1,2,3...)
 */
int RFLib_common_poly_fit(const double *data, int pointNum, int poly_order, double *coef)
{
    int status;    
    double *buf_x;                                     /* the x buffer for fitting */	
    
    /* creat buffer and init it for x */
    status = RFLib_common_buffer_create(&buf_x, pointNum);
    status = RFLib_common_buffer_init_linear(buf_x, 0, 1, pointNum);
    
    /* poly fit */
    status = MATH_poly_fit(buf_x, data, pointNum, poly_order, coef);
    
    /* delete the buffer */
    status = RFLib_common_buffer_delete(&buf_x);
    
    return status;
}

/**
 * Function to measure the derative of the input waveform
 * 32 tap FIR filter is used
 * the FIR coefficients is generated by the matlab script of sgfilter:
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%--------------------------------------------------------
% FIR filter design via local moving window LS fitting  -
% A magic smooth and derivative formula generator.      -
% By Dr Yangquan Chen		08-07-1999                   -
% Email=<yqchen@ieee.org>; URL=http://www.crosswinds.net/~yqchen/
% ------------------------------------------------------- 
% Purpose: general FIR design via local LS fitting. 
% total taps = nL+nR+1
% Format: function [c]=sgfilter(nL,nR,M,id)
%          -nL      -nL+1           -1                nR
% FIR=c(1)z   +c(2)z     +...+c(nL)z  +...+c(nL+nR+1)z
%
%       nR
%     ------
%     \                   j
%      >       c(nL+1+j) z
%     /
%     ------
%      j=-nL
% M: the order of LS fit at the moving window of [-nL, ... , nR]
% id: index for the derivative order
%		0: smooth filter, 
%		1: 1st order differentiator,
% 		2: 2nd order differentiator, ... 
% NOTE: M>=id, set M=(2~4)*(id+1) for reliably results.
%		  to do LS fit, M<nL+nR+1.
%---------------------------------------------------------------
function [c]=sgfilter(nL,nR,M,id)
% Savitzky-Golay smoothing filter.
if (id>M)
   disp('Error in id! (id<M)');return;
end
if (M>(nL+nR))
   disp('Error in M! (M<=nL+nR)');return;
end

A=zeros(nL+nR+1,M+1);
for i=-nL:nR;
   for j=0:M;
      A(i+nL+1,j+1)=i^j;
   end
end
h=zeros(M+1,1);
%h(1)=1;
h(id+1)=1;
b=(A'*A)\h;
c=zeros(nL+nR+1,1);
for n=-nL:nR
   nm=n.^[0:M];
   c(n+nL+1)=nm*b;
end
% coefficient for smoothing
return
 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
 */
int RFLib_common_mea_derivative(double *dataIn, double *dataDer, double dt, int pointNum)
{
    int i;
    int status;
   
    /* the FIR coefficients for derivative calculation */
    double fir_coef[32] = {5.681818181818181e-003,
     5.315249266862169e-003,
     4.948680351906157e-003,
     4.582111436950146e-003,
     4.215542521994134e-003,
     3.848973607038122e-003,
     3.482404692082111e-003,
     3.115835777126099e-003,
     2.749266862170087e-003,
     2.382697947214076e-003,
     2.016129032258064e-003,
     1.649560117302052e-003,
     1.282991202346041e-003,
     9.164222873900292e-004,
     5.498533724340174e-004,
     1.832844574780058e-004,
    -1.832844574780059e-004,
    -5.498533724340175e-004,
    -9.164222873900293e-004,
    -1.282991202346041e-003,
    -1.649560117302053e-003,
    -2.016129032258065e-003,
    -2.382697947214076e-003,
    -2.749266862170088e-003,
    -3.115835777126100e-003,
    -3.482404692082111e-003,
    -3.848973607038122e-003,
    -4.215542521994135e-003,
    -4.582111436950146e-003,
    -4.948680351906157e-003,
    -5.315249266862170e-003,
    -5.681818181818181e-003};
   
    /* check the input */
    if(!dataIn || !dataDer) return 1;
    if(pointNum <= 0)       return 2;
    if(dt <= 0)             return 3;

    /* scale the coefficients for time step */
    for(i = 0; i < 32; i ++) 
        fir_coef[i] /= dt; 
      
    /* FIR filtering */
    status = DSP_FIR_filter(dataIn, fir_coef, dataDer, pointNum, 32);
    
    return status;
}

/**
 * function to perform the circle fitting
 * x: x values
 * y: y values
 * result: [x0, y0, a]
 *
 * Matlab code for reference:
 ---------------------------------------
 function [x0, y0, a] = fun_fit_circle(x, y)
%============================================
% fit circle of x,y
% (x - x0)^2 + (y - y0)^2 = a^2
% Zheqiao Geng
% Feb 06, 2010
%============================================
n = length(x);

A = zeros(n, 3);
B = zeros(n, 1);

for ii = 1 : n
    A(ii, 1) = 2 * x(ii);
    A(ii, 2) = 2 * y(ii);
    A(ii, 3) = 1;
    
    B(ii)    = x(ii) * x(ii) + y(ii) * y(ii);
end;

X = A\B;

x0 = X(1);
y0 = X(2);
a  = sqrt(X(3) + X(1)^2 + X(2)^2);
 ---------------------------------------
 */
int RFLib_common_fit_circle(double *x, double *y, int point_num, double *result)
{
    int i;
    int status;
    
    RFLIB_DATA_MATRIX A;                    /* matrix for linear fitting */
    double *B;                              /* vector for linear fitting */
    double X[3];
    
    /* create the matrix and buffer */
    status   = RFLib_common_matrix_create(&A, point_num, 3);
    status   = RFLib_common_buffer_create(&B, point_num);
    
    /* make up the matrix */
    for(i = 0; i < point_num; i ++)
    {
        RFLib_common_matrix_set_element(&A, i, 0, *(x + i) * 2);	
        RFLib_common_matrix_set_element(&A, i, 1, *(y + i) * 2);
        RFLib_common_matrix_set_element(&A, i, 2, 1);
        
        *(B + i) = *(x + i) * (*(x + i)) + *(y + i) * (*(y + i));
    }
    
    /* linear fitting */
    status   = RFLib_common_linear_fit(&A, B, X);  
    
    *(result)     = X[0];
    *(result + 1) = X[1];
    *(result + 2) = sqrt(X[2] + X[0] * X[0] + X[1] * X[1]);
    
    /* delete the matrix and buffer */
    status   = RFLib_common_matrix_delete(&A);
    status   = RFLib_common_buffer_delete(&B); 	
    
    return status;
}

/**
 * generate a circle table x, y
 */
void RFLib_common_gen_circle(double x0, double y0, double a)
{
    /* TBD later */
}

/**
 * function to remove the points with large error during circle fitting
 * Input:
 *      x, y, point_num, result: same as the function for RFLib_common_fit_circle
 *      xerr_threshold: threshold for x error, as times of the rms value
 *      yerr_threshold: threshold for y error, as times of the rms value
 *      aerr_threshold: threshold for radius error, as ratio to the fitted radius of the circle
 *      max_iteration: the max iteration
 *
 * Matlab code as reference:
 ---------------------------------------
 function [x0, y0, a] = fun_fit_circle_opt(x, y, xerr_threshold, yerr_threshold, aerr_threshold)
%============================================
% fit circle of x,y
% (x - x0)^2 + (y - y0)^2 = a^2
% Zheqiao Geng
% Feb 06, 2010
%
% Zheqiao Geng
% June 4, 2010
% optimization, remove the wrong points
%============================================

n    = length(x);
flag = ones(n, 1);
x1   = zeros(n, 1);
y1   = zeros(n, 1);

% first we need to remove the points far from the average, which will
% destroy the fitting quickly
rms_x = std(x);
rms_y = std(y);
avg_x = mean(x);
avg_y = mean(y);

for ii = 1 : n
   if abs(x(ii) - avg_x) > xerr_threshold * rms_x || abs(y(ii) - avg_y) > yerr_threshold * rms_y
       flag(ii) = 0;
   end;
end;

% do the optimization iteratively
for ii = 1:100
   ii
   
   % select data
   pno  = 0;
   for jj = 1 : n
      if flag(jj) == 1
         pno     = pno + 1;
         x1(pno) = x(jj);
         y1(pno) = y(jj);
      end;
   end;

   x1 = x1(1:pno);
   y1 = y1(1:pno);

   % do a fit
   [x0, y0, a] = fun_fit_circle(x1, y1);

   % search the point with large errors
   err_pno = 0;
   for jj = 1 : n
       if flag(jj) == 1
           if abs(sqrt((x(jj) - x0)^2 + (y(jj) - y0)^2) - a) / a > aerr_threshold
               flag(jj) = 0;
               err_pno = err_pno + 1;
           end;
       end;
   end;
   
   err_pno
   
   if err_pno == 0
       break;
   end;
end;
 ---------------------------------------
 */
int  RFLib_common_fit_circle_opt(double *x, 
                                     double *y, 
                                     int point_num, 
                                     double *result,
                                     double xerr_threshold, 
                                     double yerr_threshold,
                                     double aerr_threshold, 
                                     int max_iteration)
{   
   int i, j;
   int status;
   int cur_pno, err_pno;             /* point number of this iteration and the error point number */
   
   double *x1, *y1;              /* store the valid date for fitting */
   double *flag;                 /* flags to show the status of the points, 1 means valid (should use int or bool in the future) */
   
   double rms_x, rms_y;          /* rms values (standard deviation) of x and y */
   double avg_x, avg_y;          /* average values of x and y */
   
   /* create new buffers */
   status   = RFLib_common_buffer_create(&x1, point_num);
   status   = RFLib_common_buffer_create(&y1, point_num);
   status   = RFLib_common_buffer_create(&flag, point_num);
   
   /* init the buffers */
   for(i = 0; i < point_num; i ++)
   {
       *(x1 + i)   = 0.0;
       *(y1 + i)   = 0.0;
       *(flag + i) = 1.0;	
   }
   
   /* remove the point large from the average */
   rms_x = MATH_statistic_std((const double *)x, point_num);
   rms_y = MATH_statistic_std((const double *)y, point_num);
   avg_x = MATH_statistic_avg((const double *)x, point_num);
   avg_y = MATH_statistic_avg((const double *)y, point_num);
   
   for(i = 0; i < point_num; i ++)
       if(fabs(*(x + i) - avg_x) > xerr_threshold * rms_x || fabs(*(y + i) - avg_y) > yerr_threshold * rms_y) *(flag + i) = 0.0;

   /* fit the circle iteratively */
   for(j = 0; j < max_iteration; j ++)
   {
       /* select data for the current iteration */
       cur_pno = 0;
       
       for(i = 0; i < point_num; i ++)
       {
           if(*(flag + i) == 1.0)
           {
               *(x1 + cur_pno) = *(x + i);
               *(y1 + cur_pno) = *(y + i);
               cur_pno ++;               	
           }     	
       } 
       
       /* make a fit */
       status = RFLib_common_fit_circle(x1, y1, cur_pno, result);
       
       /* search the point with large errors */
       err_pno = 0;
       
       for(i = 0; i < point_num; i ++)
       {
          if(*(flag + i) == 1.0)
          {
              if(sqrt((*(x + i) - *(result)) * (*(x + i) - *(result)) + (*(y + i) - *(result + 1)) * (*(y + i) - *(result + 1))) / *(result + 2) > aerr_threshold)
              {
                  *(flag + i) = 0.0;
                  err_pno ++;	
              }	
          }	
       }
   
       /* if all points OK, break the iteration */
       if(err_pno == 0) break;
   }
   
   /* delete the buffers */
   status   = RFLib_common_buffer_delete(&x1); 		
   status   = RFLib_common_buffer_delete(&y1); 	
   status   = RFLib_common_buffer_delete(&flag); 	
   
   return status;
}

/**
 * Do the fitting of a cos function, this is useful for phasing algorithm and I/Q calibration
 *     value = ampOff + amp * cos(pha_deg - phaOff_deg)
 * return:
 *     0 - successful
 *    -1 - failed
 */
int RFLib_common_fitCos(double *pha_deg, double *value, int pno, double *phaOff_deg, double *amp, double *ampOff)
{
    int i;
    int status;
    
    RFLIB_DATA_MATRIX A;        /* for multi linear fitting */
    double *B;
    double X[3] = {0};

    /* check the input */
    if(!pha_deg || !value || pno <= 0) return -1;    

    /* init the data structure */
    status = RFLib_common_matrix_create(&A, pno, 3);
    status = RFLib_common_buffer_create(&B, pno);

    if(status != 0) return -1;

    /* make up the matrix and solve the linear equation */
    for(i = 0; i < pno; i ++) {
        RFLib_common_matrix_set_element(&A, i, 0, sin(*(pha_deg + i) * LLRF_pi / 180.0));
        RFLib_common_matrix_set_element(&A, i, 1, cos(*(pha_deg + i) * LLRF_pi / 180.0));
        RFLib_common_matrix_set_element(&A, i, 2, 1);

        *(B + i) = *(value + i);
    }

    RFLib_common_linear_fit(&A, B, X);    

    /* get the output */
    if(phaOff_deg) *phaOff_deg  = RFLib_common_complex_angle(X[1], X[0]) * 180.0 / LLRF_pi;
    if(amp)        *amp         = RFLib_common_complex_abs(X[1], X[0]);
    if(ampOff)     *ampOff      = X[2];

    /* delete the matrix and buffers */
    RFLib_common_matrix_delete(&A);
    RFLib_common_buffer_delete(&B);

    return 0;
}



