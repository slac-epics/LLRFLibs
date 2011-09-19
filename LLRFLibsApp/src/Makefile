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

# ---- private source files ----
LLRFLibs_SRCS += RFLib_signalProcess.c
LLRFLibs_SRCS += MathLib_dataProcess_bitCalc.c
LLRFLibs_SRCS += MathLib_dataProcess_statistics.c

# ---- finally link to the EPICS Base libraries ----
LLRFLibs_LIBS += $(EPICS_BASE_IOC_LIBS)

include $(TOP)/configure/RULES
















