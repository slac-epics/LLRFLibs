/****************************************************
 * LLRFLibs_misc.h
 * 
 * Some small and useful routines for the LLRF libs
 *
 * Created by: Zheqiao Geng, gengzq@slac.stanford.edu
 * Created on: 2011.10.13
 * Description: Initial creation
 ****************************************************/
#ifndef LLRFLIBS_MISC_H
#define LLRFLIBS_MISC_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#ifdef __cplusplus
extern "C" {
#endif

/**
 * String functions
 */
void stringReplace(char *str, char src, char des);
char *getSystemTime();

#ifdef __cplusplus
}
#endif

#endif





