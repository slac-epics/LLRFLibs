/****************************************************
 * MathLib_dataProcess_bitCalc.c
 * 
 * Source file for the mathematic data and routines. For bit calculation.
 *
 * Created by: Zheqiao Geng, gengzq@slac.stanford.edu
 * Created on: 2011.05.18
 * Description: Initial creation
 ****************************************************/
#include "MathLib_dataProcess.h"

/**
 * Split the unsigned int array to two arraies, one from the higher 16 bits, another from the lower 16 bits
 * Input:
 *   data           : Source data of unsigned int
 *   pointNum       : Number of the points in data
 * Output:
 *   arrayHi        : Short array from the higher 16 bits
 *   arrayLo        : Short array from the lower 16 bits
 * Return:
 *   0              : Successful
 *  -1              : Failed
 */
int MATHLIB_u32ToShortArray(unsigned int *data, int pointNum, short *arrayHi, short *arrayLo)
{
    int i;

    /* Check the input */
    if(!data || pointNum <= 0) return -1;

    /* Split the data */
    if(arrayHi && arrayLo) {
        for(i = 0; i < pointNum; i ++) {
            *(arrayHi + i) = MATHLIB_u32ToShortHi(data + i);
            *(arrayLo + i) = MATHLIB_u32ToShortLo(data + i);
        }
    } else if(arrayHi) {
        for(i = 0; i < pointNum; i ++) {
            *(arrayHi + i) = MATHLIB_u32ToShortHi(data + i);
        }
    } else if(arrayLo) {
        for(i = 0; i < pointNum; i ++) {
            *(arrayLo + i) = MATHLIB_u32ToShortLo(data + i);
        }
    } else {
        return -1;
    }

    return 0;
}

/**
 * Get the subarray for every N point
 * Input:
 *   data           : Source data of short
 *   pointNum       : Number of the points in source data
 *   idOffset       : 0 count from the first point, 1 count from the second point
 *   everyN         : get data every N point (minimum 1)
 * Output:
 *   subArray       : subarray result (the buffer size should be large enough)
 * Return:
 *   0              : Successful
 *  -1              : Failed
 */
int MATHLIB_getEveryNSubArray(short *data, int pointNum, short *subArray, int idOffset, int everyN)
{
    int    i;
    int    dataId    = 0;

    short *dataStart = data + idOffset;
    int    pno       = pointNum - idOffset;

    /* check the input */
    if(!data || !subArray || idOffset < 0 || everyN < 1) return -1;

    /* get the data */
    for(i = 0; i < pno; i += everyN) {
        *(subArray + dataId) = *(dataStart + i);
        dataId ++;
    }

    return 0;
}








