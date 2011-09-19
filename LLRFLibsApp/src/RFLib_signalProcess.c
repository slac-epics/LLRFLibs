/****************************************************
 * RFLib_signalProcess.c
 * 
 * Realize the RF Signal processing routines
 *
 * Created by: Zheqiao Geng, gengzq@slac.stanford.edu
 * Created on: 2011.05.18
 * Description: Itial creation
 *
 * Modified by: Zheqiao Geng
 * Modified on: 2011.05.24
 * Description: 1. Simplify the definition of the RFWaveform, keep only the short buffer for I/Q
 *              2. Simplify the routine implementation to handle only the short buffers
 *
 * Modified by: Zheqiao Geng
 * Modified on: 2011.08.27
 * Description: Add the function to limit the phase measurement results (-360, 360)
 ****************************************************/
#include "MathLib_dataProcess.h"
#include "RFLib_signalProcess.h"

/**
 * Init the RF waveform
 * Return:
 *   0          : Successful
 *  -1          : Failed
 */
int RFLIB_initRFWaveform(RFLIB_struc_RFWaveform *wf, long pointNum)
{
    /* check the input */
    if(!wf) return -1;

    /* clear */
    memset((void *)wf, 0, sizeof(RFLIB_struc_RFWaveform));
    
    /* init values */        
    wf -> pointNum  = pointNum;   
    wf -> ampScale  = 1.0;
    
    return 0;
}

/**
 * Init the analog waveform
 * Return:
 *   0          : Successful
 *  -1          : Failed
 */
int RFLIB_initAnalogWaveform(RFLIB_struc_analogWaveform *wf, long pointNum)
{
    /* check the input */
    if(!wf) return -1;

    /* clear */
    memset((void *)wf, 0, sizeof(RFLIB_struc_analogWaveform));
    
    /* init values */        
    wf -> pointNum  = pointNum;   
    wf -> ampScale  = 1.0;
    
    return 0;
}

/**
 * Calculate the phase (radian) of the input vector. The value is normalized to the range of [-pi, pi]
 */
double RFLIB_vectorPha_rad_norm(double I, double Q)
{
   double var_pha_rad;
   
   if(I == 0 && Q == 0)
      var_pha_rad = 0;
   else if(I == 0 && Q > 0)
      var_pha_rad =   RFLIB_CONST_PI2;
   else if(I == 0 && Q < 0)
      var_pha_rad = - RFLIB_CONST_PI2;
   else if(I > 0)
      var_pha_rad = atan(Q / I);      
   else if(Q > 0)   
      var_pha_rad = atan(Q / I) + RFLIB_CONST_PI;
   else
      var_pha_rad = atan(Q / I) - RFLIB_CONST_PI;   
      
   return var_pha_rad;                      
}

/**
 * Normalize the phase in degree to the limited range
 */
double RFLIB_normPha_deg(double pha_deg, double hiLimit, double loLimit)
{
    double result = pha_deg;

    while(result > hiLimit || result < loLimit) {
        result = (result > hiLimit) ? result - 360 : result;
        result = (result < loLimit) ? result + 360 : result;
    }

    return result;
}

/**
 * Demodulate the sample of the RF signal into I and Q components
 * Input:
 *   rawDataBuf         : The raw sampling
 *   sampleNum          : Sample number for integer RF cycles
 * Output:
 *   IOutBuf            : I component
 *   QOutBuf            : Q component
 * Return:
 *   0                  : Successful
 *  -1                  : Failed
 */
int RFLIB_rfDemod(short *rawDataBuf, int pointNum, int demodIdCur, short *IOutBuf, short *QOutBuf)
{
    int i;   
    int j;
    
    double var_cs[RFLIB_CONST_SAMPLE_NUM];              /* cos coefficient for discrete Fourior transform */
    double var_sn[RFLIB_CONST_SAMPLE_NUM];              /* sin coefficient for discrete Fourior transform */
    
    double var_pipeI[RFLIB_CONST_SAMPLE_NUM + 1] = {0.0};   /* shift registers */
    double var_pipeQ[RFLIB_CONST_SAMPLE_NUM + 1] = {0.0};
    
    double var_phaStep = RFLIB_CONST_CYCLE_NUM * 2 * RFLIB_CONST_PI / RFLIB_CONST_SAMPLE_NUM;
    
    double var_accI     = 0;
    double var_accQ     = 0;
    int    var_index    = demodIdCur;
    
    /* Check the input */
    if(!rawDataBuf || !IOutBuf || !QOutBuf || pointNum <= 0 || demodIdCur < 0) return -1;       /* we still can not guarantee the size of the buffers */
    
    /* Calculate the coefficients */
    for(i = 0; i < RFLIB_CONST_SAMPLE_NUM; i ++) {
        var_cs[i] = cos(i * var_phaStep);
        var_sn[i] = sin(i * var_phaStep);   
    }
    
    /* Demodulation */
    for(i = 0; i < pointNum; i ++) {        
        /* accumulation */
        var_accI += *(rawDataBuf + i) * var_cs[var_index] / RFLIB_CONST_SAMPLE_NUM;
        var_accQ += *(rawDataBuf + i) * var_sn[var_index] / RFLIB_CONST_SAMPLE_NUM;
        
        /* shift registers */
        for(j = RFLIB_CONST_SAMPLE_NUM; j > 0; j--) {
            var_pipeI[j] = var_pipeI[j - 1];
            var_pipeQ[j] = var_pipeQ[j - 1];
        }        
        
        /* append to the end */
        var_pipeI[0] = var_accI;
        var_pipeQ[0] = var_accQ;
        
        var_index ++;
        if(var_index >= RFLIB_CONST_SAMPLE_NUM) var_index = 0;
        
        /* get the demodulation result */        
        *(IOutBuf + i) = (short)(var_accI - var_pipeI[RFLIB_CONST_SAMPLE_NUM]);
        *(QOutBuf + i) = (short)(var_accQ - var_pipeQ[RFLIB_CONST_SAMPLE_NUM]);
    }
    
    return 0;
}

/**
 * Demodulate the RF waveform
 * Input:
 *   wf         : RF waveform object
 * Return:
 *   0          : Successful
 *  -1          : Failed
 */
int RFLIB_RFWaveformDemod(RFLIB_struc_RFWaveform *wf)
{
    int status = 0;
    int pno;

    /* Check the input */
    if(!wf) return -1;
    
    /* Demodulation */
    pno = MATHLIB_min(wf->pointNum, RFLIB_CONST_WF_SIZE);               /* to avoid overflow */

    status = RFLIB_rfDemod(wf->wfRaw, pno, (int)wf->demodCoefIdCur, wf->wfI, wf->wfQ);
    
    return status;
}

/**
 * Convert the I/Q waveform to A/P waveform
 * Input:
 *   wf         : RF waveform object
 * Return:
 *   0          : Successful
 *  -1          : Failed
 */
int RFLIB_RFWaveformIQ2AP(RFLIB_struc_RFWaveform *wf)
{
    int i;
    int pno;
    
    /* Check the input */
    if(!wf) return -1;
    
    /* Convert */
    pno = MATHLIB_min(wf->pointNum, RFLIB_CONST_WF_SIZE);                   /* to avoid overflow */

    /* calculate from the short buffer for I/Q */
    for(i = 0; i < pno; i ++) {
        wf->wfAmp[i]     = RFLIB_vectorAmp((double)wf->wfI[i], (double)wf->wfQ[i]) * wf->ampScale;
        wf->wfPha_deg[i] = RFLIB_normPha_deg(RFLIB_vectorPha_deg_norm((double)wf->wfI[i], (double)wf->wfQ[i]) + wf->phaOffset_deg, 180, -180);
    }
    
    return 0;
}

/**
 * Get the average I/Q/A/P
 * Input:
 *   wf         : RF waveform object
 * Return:
 *   0          : Successful
 *  -1          : Failed
 */
int RFLIB_RFWaveformAvg(RFLIB_struc_RFWaveform *wf)
{
    int avgStartPoint;
    int avgLength;

    /* Check the input */
    if(!wf) return -1;    

    /* get the average start point and length */
    avgStartPoint = (int)((wf->avgStartTime_ns - wf->sampleDelay_ns) * wf->sampleFreq_MHz / 1000.0);
    avgLength     = (int)(wf->avgTime_ns * wf->sampleFreq_MHz / 1000.0);      
    
    if(avgStartPoint + avgLength > RFLIB_CONST_WF_SIZE) {
        avgLength = RFLIB_CONST_WF_SIZE - avgStartPoint - 1;
    }

    /* Average from the short I/Q */
    wf->avgDataI = (double)MATHLIB_avg_short(wf->wfI + avgStartPoint, avgLength);
    wf->avgDataQ = (double)MATHLIB_avg_short(wf->wfQ + avgStartPoint, avgLength);

    /* Get the amplitude and phase from the averaged I/Q */    
    wf->avgDataAmp      = RFLIB_vectorAmp(wf->avgDataI, wf->avgDataQ) * wf->ampScale;
    wf->avgDataPha_deg  = RFLIB_normPha_deg(RFLIB_vectorPha_deg_norm(wf->avgDataI, wf->avgDataQ) + wf->phaOffset_deg, 180, -180);

    return 0;
}

/**
 * Scale the analog waveform to physical unit 
 * Input:
 *   wf         : analog waveform object
 * Return:
 *   0          : Successful
 *  -1          : Failed
 */
int RFLIB_analogWaveformScale(RFLIB_struc_analogWaveform *wf)
{
    int i;
    
    /* check the input */
    if(!wf) return -1;
    
    /* scale */
    for(i = 0; i < wf->pointNum; i ++) {
        wf->wfAmp[i] = (double)(wf-> wfRaw[i] * wf->ampScale);
    }    

    return 0;
}

/**
 * Average the analog waveform
 * Input:
 *   wf         : analog waveform object
 * Return:
 *   0          : Successful
 *  -1          : Failed
 */
int RFLIB_analogWaveformAvg(RFLIB_struc_analogWaveform *wf)
{
    int avgStartPoint;
    int avgLength;

    /* Check the input */
    if(!wf) return -1;    

    /* get the average start point and length */
    avgStartPoint = (int)((wf->avgStartTime_ns - wf->sampleDelay_ns) * wf->sampleFreq_MHz / 1000.0);
    avgLength     = (int)(wf->avgTime_ns * wf->sampleFreq_MHz / 1000.0);      
    
    if(avgStartPoint + avgLength > RFLIB_CONST_WF_SIZE) {
        avgLength = RFLIB_CONST_WF_SIZE - avgStartPoint - 1;
    }

    /* Average from raw data */
    wf->avgData = (double)MATHLIB_avg_short(wf->wfRaw + avgStartPoint, avgLength) * wf->ampScale;
    
    return 0;
}













