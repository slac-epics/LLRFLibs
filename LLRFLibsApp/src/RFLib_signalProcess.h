/****************************************************
 * RFLib_signalProcess.h
 * 
 * Header file for the RF Signal processing data and routines
 *
 * Created by: Zheqiao Geng, gengzq@slac.stanford.edu
 * Created on: 2011.05.18
 * Description: Initial creation
 *
 * Modified by: Zheqiao Geng
 * Modified on: 2011.05.24
 * Description: Simplify the definition of the RFWaveform, keep only the short buffer for I/Q
 *
 * Modified by: Zheqiao Geng
 * Modified on: 2011.06.28
 * Description: Add scaling to the amplitude calculation (to convert from digits to physical units),
 *              New algorithm for demodulation
 ****************************************************/
#ifndef RFLIB_SIGNAL_PROCESS_H
#define RFLIB_SIGNAL_PROCESS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/*======================================
 * Constant definitions
 *======================================*/
#define RFLIB_CONST_PI          3.141592653589793   /* pi */
#define RFLIB_CONST_2PI         6.283185307179586   /* 2 * pi */
#define RFLIB_CONST_PI2         1.570796326794897   /* pi / 2 */

#define RFLIB_CONST_WF_SIZE     65536               /* buffer size of the waveforms */

#define RFLIB_CONST_SAMPLE_NUM  14                  /* sample number */
#define RFLIB_CONST_CYCLE_NUM   3                   /* RF cycles fully sampled by the point number of RFLIB_CONST_SAMPLE_NUM */

/*======================================
 * Data models
 *======================================*/
/* Data model for the RF waveform */
typedef struct {    
    volatile unsigned long chId;                    /* (note: all chId is for 16 bits, can be any number) raw: any number; I/Q: 0,2,4,... */
    volatile unsigned short valid;                  /* 1 to indicate that this waveform is valid */

    long pointNum;                                  /* number of the effective points */
    long demodCoefIdCur;                            /* ID offset of the coefficient */
    
    double sampleFreq_MHz;                          /* sampling frequency of the waveform */
    double sampleDelay_ns;                          /* delay of the waveform respect to the main trigger */
    volatile double avgStartTime_ns;                /* average start time */
    volatile double avgTime_ns;                     /* average time (point index will be calculated from them) */
    
    short   wfRaw[RFLIB_CONST_WF_SIZE];             /* raw data from the ADC */

    short   wfI[RFLIB_CONST_WF_SIZE];               /* I/Q */
    short   wfQ[RFLIB_CONST_WF_SIZE];

    double  wfAmp[RFLIB_CONST_WF_SIZE];             /* amplitude waveform */
    double  wfPha_deg[RFLIB_CONST_WF_SIZE];         /* phase waveform in degree */

    volatile double ampScale;                       /* scale factor for amplitude calculation */
    volatile double phaOffset_deg;                  /* phase offset in degree */

    volatile double avgDataI;                       /* averaged I */
    volatile double avgDataQ;                       /* averaged Q */

    volatile double avgDataAmp;                     /* averaged ampligude */
    volatile double avgDataPha_deg;                 /* averaged phase in degree */
} RFLIB_struc_RFWaveform;

/* Data model for the analog waveform */
typedef struct {
    volatile unsigned long chId;                    /* (note: all chId is for 16 bits, can be any number) raw: any number; I/Q: 0,2,4,... */    

    long pointNum;                                  /* number of the effective points */
    
    double sampleFreq_MHz;                          /* sampling frequency of the waveform */
    double sampleDelay_ns;                          /* delay of the waveform respect to the main trigger */
    volatile double avgStartTime_ns;                /* average start time */
    volatile double avgTime_ns;                     /* average time (point index will be calculated from them) */

    short   wfRaw[RFLIB_CONST_WF_SIZE];             /* raw data from the ADC */
    double  wfAmp[RFLIB_CONST_WF_SIZE];             /* amplitude in physical unit */
    
    volatile double ampScale;                       /* scale factor for amplitude calculation */
    volatile double avgData;                        /* averaged data */
} RFLIB_struc_analogWaveform;

/*======================================
 * RF signal processing routines
 *======================================*/
/* Common routines */
#define RFLIB_degToRad(angle_deg)   ((angle_deg) * RFLIB_CONST_PI / 180.0)
#define RFLIB_radToDeg(angle_rad)   ((angle_rad) * 180.0 / RFLIB_CONST_PI)

#define RFLIB_vectorAmp(I, Q)       (sqrt((I) * (I) + (Q) * (Q)))                       /* calculate the amplitude of the vector I + jQ */
#define RFLIB_vectorPha_rad(I, Q)   (atan((Q) / (I)))                                   /* calculate the phase, radian */
#define RFLIB_vectorPha_deg(I, Q)   (atan((Q) / (I)) * 180.0 / RFLIB_CONST_PI)          /* calculate the phase, degree */

double  RFLIB_vectorPha_rad_norm(double I, double Q);                                   /* cacluate the phase (normalized to [-pi, pi] radian) */ 

#define RFLIB_vectorPha_deg_norm(I, Q) (RFLIB_vectorPha_rad_norm((I), (Q)) * 180.0 / RFLIB_CONST_PI)

double  RFLIB_normPha_deg(double pha_deg, double hiLimit, double loLimit);

/* Initialization */
int     RFLIB_initRFWaveform(RFLIB_struc_RFWaveform *wf, long pointNum);
int     RFLIB_initAnalogWaveform(RFLIB_struc_analogWaveform *wf, long pointNum);

/* RF signal processing */
int     RFLIB_rfDemod(short *rawDataBuf, int pointNum, int demodIdCur, short *IOutBuf, short *QOutBuf);

/* LLRF procedures */
int     RFLIB_RFWaveformDemod(RFLIB_struc_RFWaveform *wf);                              /* demodulate the RF waveform */
int     RFLIB_RFWaveformIQ2AP(RFLIB_struc_RFWaveform *wf);                              /* convert the I/Q to A/P waveform */
int     RFLIB_RFWaveformAvg(RFLIB_struc_RFWaveform *wf);                                /* average the RF waveform */

int     RFLIB_analogWaveformScale(RFLIB_struc_analogWaveform *wf);                      /* scale to physical unit */
int     RFLIB_analogWaveformAvg(RFLIB_struc_analogWaveform *wf);                        /* average */

#ifdef __cplusplus
}
#endif

#endif





