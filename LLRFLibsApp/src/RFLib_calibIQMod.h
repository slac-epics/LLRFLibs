/****************************************************
 * RFLib_calibIQMod.h
 * 
 * Application to calibration the IQ modulator
 *
 * Created by: Zheqiao Geng, gengzq@slac.stanford.edu
 * Created on: 2011.11.26
 * Description: Initial creation
 ****************************************************/
#ifndef RFLIB_CALIB_IQMOD_H
#define RFLIB_CALIB_IQMOD_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#define RFLIB_APPL_IQCAL_MAX_PNO    360

#ifdef __cplusplus
extern "C" {
#endif

/*======================================
 * Data models
 *======================================*/
/* Data model for IQ modulator calibration */
typedef struct {    

    /* Input data */
    int     pno;                                    /* number of the points */

    double  dataIn_I[RFLIB_APPL_IQCAL_MAX_PNO];     /* input to the RF system, I component */
    double  dataIn_Q[RFLIB_APPL_IQCAL_MAX_PNO];     /* Q component */
    double  dataOut_I[RFLIB_APPL_IQCAL_MAX_PNO];    /* output from the RF system (IQ modulator output), I component */
    double  dataOut_Q[RFLIB_APPL_IQCAL_MAX_PNO];    /* Q component */

    /* Input parameters */
    double  offset_I_old;                       /* offset value from last step */
    double  offset_Q_old;   
    double  corr_a11_old;                       /* correction matrix of last step */
    double  corr_a12_old;
    double  corr_a21_old;
    double  corr_a22_old;

    double  phaImbalanceScale;                  /* scale factor for the phase imbalance (for smooth convergence) */
    double  offsetScale;                        /* scale factor for the offset (for smooth convergence) */

    /* Intemediate results */
    double  int_ki;                             /* the intermediate value from current calculation */
    double  int_kq;
    double  int_offset_I;                       /* offset from current calculation */
    double  int_offset_Q;
    double  int_phaImbalance_deg;               /* imbalance */
    double  int_ampImbalance;   

    /* Output results from current calculation */
    double  offset_I;                           /* offset value */
    double  offset_Q;   
    double  corr_a11;                           /* correction matrix */
    double  corr_a12;
    double  corr_a21;
    double  corr_a22;        

} RFLIB_struc_calibIQMod;

/*======================================
 * Routines
 *======================================*/
int RFLIB_appl_calibIQMod(RFLIB_struc_calibIQMod *arg);
int RFLIB_func_calibIQModDistortion(double *dataIn_I, double *dataIn_Q, double *dataOut_I, double *dataOut_Q, int pno,
                                    double *ki, double *kq, double *offset_I, double *offset_Q, double *phaImbalance_deg, double *ampImbalance);

#ifdef __cplusplus
}
#endif

#endif





