##########################################################
# Makefile
# 
# Makefile for the LLRFLibs
# 
# Created by: Zheqiao Geng, gengzq@slac.stanford.edu
# Created on: 2011.07.05
# Description: Initial creation 
##########################################################
TOP=../..

include $(TOP)/configure/CONFIG

#------------------------------------------------
# build the library for LLRFLibs
#------------------------------------------------
# ---- library name ----
LIBRARY_IOC = LLRFLibs

# ---- ported include files ----
INC += RFLib_signalProcess.h
INC += MathLib_dataProcess.h
INC += EPICSLib_wrapper.h
INC += LLRFLibs_misc.h
INC += LLRF_algorithm_lib_config.h
INC += MEMM_buffer_operation.h
INC += MEMM_circular_buf_operation.h
INC += MATH_required_interface.h
INC += MATH_complex_computation.h
INC += MATH_linear_equations.h
INC += MATH_matrix_computation.h
INC += MATH_polynomial.h
INC += MATH_statistics.h
INC += DSP_required_interface.h
INC += DSP_filtering.h
INC += RFLib_required_interface.h
INC += RFLib_calibIQMod.h

# ---- private source files ----
LLRFLibs_SRCS += RFLib_signalProcess.c
LLRFLibs_SRCS += MathLib_dataProcess_bitCalc.c
LLRFLibs_SRCS += MathLib_dataProcess_statistics.c
LLRFLibs_SRCS += LLRFLibs_misc.c
LLRFLibs_SRCS += MEMM_buffer_operation.c
LLRFLibs_SRCS += MEMM_circular_buf_operation.c
LLRFLibs_SRCS += MATH_required_interface.c
LLRFLibs_SRCS += MATH_complex_computation.c
LLRFLibs_SRCS += MATH_linear_equations.c
LLRFLibs_SRCS += MATH_matrix_computation.c
LLRFLibs_SRCS += MATH_polynomial.c
LLRFLibs_SRCS += MATH_statistics.c
LLRFLibs_SRCS += DSP_required_interface.c
LLRFLibs_SRCS += DSP_filtering.c
LLRFLibs_SRCS += RFLib_required_interface.c
LLRFLibs_SRCS += RFLib_calibIQMod.c

# ---- finally link to the EPICS Base libraries ----
LLRFLibs_LIBS += $(EPICS_BASE_IOC_LIBS)

include $(TOP)/configure/RULES
















