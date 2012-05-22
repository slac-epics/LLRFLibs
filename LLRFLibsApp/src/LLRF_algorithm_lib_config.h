/****************************************************
 *  LLRF_algorithm_lib_config.h                                       
 *   
 *  Common definitions for LLRF algorithm library  
 *
 *  Created by: Zheqiao Geng, zheqiao.geng@desy.de
 *  Created on: Nov. 28, 2009
 *  Description: Define the local floating point type  
 *
 *  Modified by: Zheqiao Geng
 *  Modified on: Feb 08, 2010
 *  Description: Add the definition for error handling system              
 ****************************************************/
#ifndef LLRF_ALGORITHM_LIB_CONFIG_H
#define LLRF_ALGORITHM_LIB_CONFIG_H

/**
 * Define the common constants may be used in the library
 */
#define LLRF_pi     3.141592653589793     /* pi */
#define LLRF_2pi    6.283185307179586     /* 2 * pi */
#define LLRF_pi2    1.570796326794897     /* pi / 2 */
#define LLRF_roQ13  1036.0                /* the roQ of the 1.3GHz cavity, Ohm */
#define LLRF_Z0     50.0                  /* the impedence of the transmission line, Ohm */
#define LLRF_f013   1300                  /* the frequency of the 1.3GHz cavity, MHz */
#define LLRF_f039   3900                  /* the frequency of the 3.9GHz cavity, MHz */

/**
 * Define the error code (refer to GNU Scientific Library)
 */ 
enum {
    LLRF_SUCCESS     = 0,                 /* success to perform functions */
    LLRF_FAILURE     = -1,                /* failure, the reason is unclear */
    LLRF_ERR_DOM     = 1,                 /* input domain error, e.g sqrt(-1) */
    LLRF_ERR_RANGE   = 2,                 /* output range error, e.g. exp(1e100) */
    LLRF_ERR_NOMEM   = 3,                 /* No memory available */
    LLRF_ERR_PTR     = 4,                 /* invalide pointer */   
    LLRF_ERR_TIMING  = 5,                 /* timing relation error (especially for RF pulse) */ 
    LLRF_ERR_BUFOF   = 6                  /* buffer overflow (index exceed the max size) */ 
};

#endif

