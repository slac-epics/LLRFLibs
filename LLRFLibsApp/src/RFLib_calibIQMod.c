/****************************************************
 * RFLib_calibIQMod.c
 * 
 * Application to calibration the IQ modulator
 *
 * Created by: Zheqiao Geng, gengzq@slac.stanford.edu
 * Created on: 2011.11.26
 * Description: Initial creation
 ****************************************************/
#include "RFLib_calibIQMod.h"

#include "RFLib_required_interface.h"

/**
 * Routine to calculate the I/Q modulator distortions
 * Return:
 *     0        : successful
 *    -1        : failed
 */
int RFLIB_func_calibIQModDistortion(double *dataIn_I, double *dataIn_Q, double *dataOut_I, double *dataOut_Q, int pno,
                                    double *ki, double *kq, double *offset_I, double *offset_Q, double *phaImbalance_deg, double *ampImbalance)
{
    int i;
    int status;
    
    RFLIB_DATA_MATRIX A;        /* for multi linear fitting */
    double *B;
    double X[6] = {0};

    double var_ki;
    double var_kq;
    double var_offset_I;
    double var_offset_Q;
    double var_phaImbalance_deg;
    double var_ampImbalance;

    double a, b, c, d, k, l, gcosf, gsinf;      /* temporary variables */

    /*printf("----------------------\n");
    for(i = 0; i < pno; i ++) {
          printf("%lf, %lf, %lf, %lf\n", dataIn_I[i], dataIn_Q[i], dataOut_I[i], dataOut_Q[i]);      
    }
    printf("----------------------\n");*/

    /* check the input */
    if(!dataIn_I || !dataIn_Q || !dataOut_I || !dataOut_Q || pno <= 0) return -1;

    /* init the data structure */
    status = RFLib_common_matrix_create(&A, 2 * pno, 6);
    status = RFLib_common_buffer_create(&B, 2 * pno);

    if(status != 0) return -1;

    /* make up the matrix and solve the linear equation */
    for(i = 0; i < pno; i ++) {
        RFLib_common_matrix_set_element(&A, 2 * i, 0, *(dataIn_I + i));
        RFLib_common_matrix_set_element(&A, 2 * i, 1, *(dataIn_Q + i));
        RFLib_common_matrix_set_element(&A, 2 * i, 2, 0);
        RFLib_common_matrix_set_element(&A, 2 * i, 3, 0);
        RFLib_common_matrix_set_element(&A, 2 * i, 4, 1);
        RFLib_common_matrix_set_element(&A, 2 * i, 5, 0);

        RFLib_common_matrix_set_element(&A, 2 * i + 1, 0, 0);
        RFLib_common_matrix_set_element(&A, 2 * i + 1, 1, 0);
        RFLib_common_matrix_set_element(&A, 2 * i + 1, 2, *(dataIn_I + i));
        RFLib_common_matrix_set_element(&A, 2 * i + 1, 3, *(dataIn_Q + i));
        RFLib_common_matrix_set_element(&A, 2 * i + 1, 4, 0);
        RFLib_common_matrix_set_element(&A, 2 * i + 1, 5, 1);

        *(B + 2 * i)        = *(dataOut_I + i);
        *(B + 2 * i + 1)    = *(dataOut_Q + i);
    }

    /*printf("----------------------\n");
    MATH_matrix_print(&A);
    printf("----------------------\n");
    for(i = 0; i < 2* pno; i ++)
        printf("%lf\n", B[i]);
    printf("----------------------\n");*/

    RFLib_common_linear_fit(&A, B, X);    

    a = X[0];
    b = X[1];
    c = X[2];
    d = X[3];
    k = X[4];
    l = X[5];

    /* get the results */
    var_offset_I = (k * d - b * l) / (a * d - b * c);
    var_offset_Q = (a * l - k * c) / (a * d - b * c);

    var_ki       = (a - c * b / d) / (1 + (c / d) * (c / d));
    var_kq       = var_ki * c / d;

    gcosf        = c / var_kq;
    gsinf        = (b + var_kq) / var_ki;

    var_phaImbalance_deg    = RFLib_common_complex_angle(gcosf, gsinf) * 180.0 / LLRF_pi;
    var_ampImbalance        = RFLib_common_complex_abs(gcosf, gsinf);

    /* generate the output */
    if(ki)                  *ki                 = var_ki;
    if(kq)                  *kq                 = var_kq;
    if(offset_I)            *offset_I           = var_offset_I;
    if(offset_Q)            *offset_Q           = var_offset_Q;
    if(phaImbalance_deg)    *phaImbalance_deg   = var_phaImbalance_deg;
    if(ampImbalance)        *ampImbalance       = var_ampImbalance;

    /* delete the matrix and buffers */
    RFLib_common_matrix_delete(&A);
    RFLib_common_buffer_delete(&B);

    return 0;
}

/**
 * Routine to do the I/Q calibration
 * Input:
 *     arg      : address of the data structure
 * Return:
 *     0        : successful
 *    -1        : failed
 */
int RFLIB_appl_calibIQMod(RFLIB_struc_calibIQMod *arg)
{
    int status;

    double applied_phaImbalance_deg;        /* applied values */
    double applied_ampImbalance;
    double inc_distMatrix_a11;              /* increment distortion matrix */
    double inc_distMatrix_a12;
    double inc_distMatrix_a21;
    double inc_distMatrix_a22;

    /* check the input parameters */
    if(!arg) return -1;
    if(arg -> pno <= 0) return -1;

    /* calculate the distortion */
    status = RFLIB_func_calibIQModDistortion(arg -> dataIn_I,
                                             arg -> dataIn_Q, 
                                             arg -> dataOut_I, 
                                             arg -> dataOut_Q, 
                                             arg -> pno,
                                            &arg -> int_ki,
                                            &arg -> int_kq,
                                            &arg -> int_offset_I,
                                            &arg -> int_offset_Q,
                                            &arg -> int_phaImbalance_deg,
                                            &arg -> int_ampImbalance);

    /* calculate the new offset and pre-distortion matrix */
    applied_phaImbalance_deg    = arg -> int_phaImbalance_deg * arg -> phaImbalanceScale;
    applied_ampImbalance        = sqrt(arg -> int_ampImbalance);

    inc_distMatrix_a11          = 1;
    inc_distMatrix_a12          = -tan(applied_phaImbalance_deg * LLRF_pi / 180.0);
    inc_distMatrix_a21          = 0;
    inc_distMatrix_a22          = 1 / (applied_ampImbalance  * cos(applied_phaImbalance_deg * LLRF_pi / 180.0));

    arg -> corr_a11 = arg -> corr_a11_old * inc_distMatrix_a11 + arg -> corr_a12_old * inc_distMatrix_a21;
    arg -> corr_a12 = arg -> corr_a11_old * inc_distMatrix_a12 + arg -> corr_a12_old * inc_distMatrix_a22;
    arg -> corr_a21 = arg -> corr_a21_old * inc_distMatrix_a11 + arg -> corr_a22_old * inc_distMatrix_a21;
    arg -> corr_a22 = arg -> corr_a21_old * inc_distMatrix_a12 + arg -> corr_a22_old * inc_distMatrix_a22;

    arg -> offset_I = arg -> offset_I_old + arg -> offsetScale * arg -> int_offset_I;
    arg -> offset_Q = arg -> offset_Q_old + arg -> offsetScale * arg -> int_offset_Q;

    return status;
}







