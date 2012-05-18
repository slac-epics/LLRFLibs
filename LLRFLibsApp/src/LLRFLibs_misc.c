/****************************************************
 * LLRFLibs_misc.c
 * 
 * Some small and useful routines for the LLRF libs
 *
 * Created by: Zheqiao Geng, gengzq@slac.stanford.edu
 * Created on: 2011.10.13
 * Description: Initial creation
 ****************************************************/
#include "LLRFLibs_misc.h"

/**
 * replate the char of "src" to "des"
 */
void stringReplace(char *str, char src, char des)
{    
    char *p;
    while((p = strchr(str, src))) *p = des;
}

/**
 * Another way to get the system time
 */
char *getSystemTime()
{
    time_t timer;
    timer = time(0);
    return asctime(localtime(&timer));
}



