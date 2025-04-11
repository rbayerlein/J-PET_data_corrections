/*********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*               (C) Copyright 2005-2012 Department of Radiology			  		*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/*********************************************************************************
*
*     Module Name:        LbTiming.h
*     Revision Number:    1.2
*     Date last revised:  19 December 2012
*     Programmer:         Steven Gillispie
*     Date Originated:    13 July 2005
*
*     Module Overview:    This module declares all constants, types, variables
*							and prototypes for the Timing library.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*     Global variables defined:   
*
*********************************************************************************/

#ifndef LIB_TIMING
#define LIB_TIMING


#ifdef GEN_UNIX
	#include <time.h>
#endif
#ifdef __MWERKS__
	#include <time.h>
#endif
#ifdef DARWIN
	#include <time.h>
#endif
#ifdef OSX
	#include <time.h>
#endif
#ifdef COCOA
	#include <time.h>
#endif
#ifdef WINNT
	#include <time.h>
#endif

#include "LbTypes.h"


/* CONSTANTS */

/* FLAGS */

/* TYPES */
typedef	struct {
	#undef				LBTM_TYPE_DEFINED
	
	/* The following types are only defined for these listed systems (as far as we know) */
	#ifdef GEN_UNIX
		clock_t				wallClockTime;		/* Real elapsed time */
		struct tms			cpuClockTime;		/* CPU elapsed time */
		#define			LBTM_TYPE_DEFINED
	#endif
	
	#ifdef __MWERKS__
		clock_t				wallClockTime;		/* Real elapsed time */
		clock_t				cpuClockTime;		/* CPU elapsed time */
		#define			LBTM_TYPE_DEFINED
	#endif
	
	#ifdef DARWIN
		clock_t				wallClockTime;		/* Real elapsed time */
		struct tms			cpuClockTime;		/* CPU elapsed time */
		#define			LBTM_TYPE_DEFINED
	#endif
	
	#ifdef OSX
		clock_t				wallClockTime;		/* Real elapsed time */
		struct tms			cpuClockTime;		/* CPU elapsed time */
		#define			LBTM_TYPE_DEFINED
	#endif
	
	#ifdef COCOA
		clock_t				wallClockTime;		/* Real elapsed time */
		struct tms			cpuClockTime;		/* CPU elapsed time */
		#define			LBTM_TYPE_DEFINED
	#endif
	
	/* Windows NT allows for tracking of CPU time spent on a task - in theory*/
	#ifdef WINNT
		clock_t				wallClockTime;		/* Real elapsed time */
		clock_t				cpuClockTime;		/* CPU elapsed time */
		#define			LBTM_TYPE_DEFINED
	#endif
	
	/* Other systems do not support timing capabilities (as far as we know) */
	#ifndef LBTM_TYPE_DEFINED
		int					wallClockTime;		/* Real elapsed time */
		int					cpuClockTime;		/* CPU elapsed time */
	#endif
} LbTmTimingType;		/* System-independent timing variable type */


/* PROTOTYPES */
void LbTmStartTiming(LbTmTimingType *startTime);
Boolean LbTmStopTiming(LbTmTimingType *startTime, 
						double *realTimeIntervalSecs, double *cpuTimeIntervalSecs);

#endif /* LIB_TIMING */
