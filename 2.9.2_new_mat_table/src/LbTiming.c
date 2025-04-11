/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 2005-2012 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
********************************************************************************/

/********************************************************************************
*
*     Module Name:        LbTiming.c
*     Revision Number:    1.3
*     Date last revised:  19 December 2012
*     Programmer:         Steven Gillispie
*     Date Originated:    13 July 2005
*
*     Module Overview:    This module provides routines for performing timing
*							tasks.
*						  MACHINE DEPENDENT						
*
*     References:         
*
*********************************************************************************
*
*     Global functions defined:
*		LbTmStartTiming
*		LbTmStopTiming
*
*     Global variables defined:   
*
**********************************************************************************
*
*	  Revision Section (Also update version number, if relevant)
*
*	  Programmer(s):		
*
*	  Revision date:		
*
*	  Revision description:
*
**********************************************************************************
*
*	  Revision Section (Also update version number, if relevant)
*
*	  Programmer(s):		Steven Gillispie
*
*	  Revision date:		19 December 2012
*
*	  Revision description:	Added back more timing options by system
*
**********************************************************************************
*
*	  Revision Section (Also update version number, if relevant)
*
*	  Programmer(s):		Steven Gillispie
*
*	  Revision date:		10 October 2012
*
*	  Revision description:	Reduced timing options by system to just two
*
********************************************************************************/


#include "SystemDependent.h"

#include "LbTiming.h"


/*	CONSTANTS */

/*	LOCAL TYPES */

/*  LOCAL GLOBALS */

/*	LOCAL MACROS */

/*	LOCAL FUNCTIONS */

/*	FUNCTIONS */

/*****************************************
*		LbTmStartTiming
*
*	Arguments:
*				LbTmTimingType		*startTime	- Pointer to start time (initialized here).
*
*	Purpose:	Start the timing and return the start time in startTime.
*				MACHINE DEPENDENT
*
*	Returns:	
*
******************************************/
void LbTmStartTiming(LbTmTimingType *startTime)

{
	/* Get the start times, if on a system that supports timing */
	#ifdef GEN_UNIX
		startTime->wallClockTime = times(&(startTime->cpuClockTime));
	#elif DARWIN
		startTime->wallClockTime = times(&(startTime->cpuClockTime));
	#elif OSX
		startTime->wallClockTime = times(&(startTime->cpuClockTime));
	#elif COCOA
		startTime->wallClockTime = times(&(startTime->cpuClockTime));
	#else
		startTime->wallClockTime = clock();
		startTime->cpuClockTime = startTime->wallClockTime;
	#endif
	
	/* Otherwise, just set start times to zero */
	#ifndef LBTM_TYPE_DEFINED
		startTime->wallClockTime = 0;
		startTime->cpuClockTime = 0;
	#endif
}


/*****************************************
*		LbTmStopTiming
*
*	Arguments:
*				LbTmTimingType	*startTime				- Pointer to start time.
*				double			*realTimeIntervalSecs	- Returned real time interval.
*				double			*cpuTimeIntervalSecs	- Returned CPU time interval.
*
*	Purpose:	Return the real and CPU elapsed time intervals since startTime.
*				MACHINE DEPENDENT
*
*	Returns:	True if valid times; false if timing not supported on this system.
*
******************************************/
Boolean LbTmStopTiming(LbTmTimingType *startTime, 
						double *realTimeIntervalSecs, double *cpuTimeIntervalSecs)

{
	LbTmTimingType		stopTime;			/* Local stop time */
	LbUsFourByte		clockTicks;			/* Time interval converter value */
	Boolean				result;				/* Result of function */
	
	
	/* Get the stop times, if on a system that supports timing */
	{
		#ifdef GEN_UNIX
			stopTime.wallClockTime = times(&(stopTime.cpuClockTime));
		#elif DARWIN
			stopTime.wallClockTime = times(&(stopTime.cpuClockTime));
		#elif OSX
			stopTime.wallClockTime = times(&(stopTime.cpuClockTime));
		#elif COCOA
			stopTime.wallClockTime = times(&(stopTime.cpuClockTime));
		#else
			stopTime.wallClockTime = clock();
			stopTime.cpuClockTime = stopTime.wallClockTime;
		#endif
		
		/* Otherwise, just set stop times to zero */
		#ifndef LBTM_TYPE_DEFINED
			stopTime.wallClockTime = 0;
			stopTime.cpuClockTime = 0;
		#endif
	}
	
	
	/* Get time interval converter value (timing granularity) */
	{
		#ifdef GEN_UNIX
			clockTicks = sysconf(_SC_CLK_TCK);
		#elif DARWIN
			clockTicks = sysconf(_SC_CLK_TCK);
		#elif OSX
			clockTicks = sysconf(_SC_CLK_TCK);
		#elif COCOA
			clockTicks = sysconf(_SC_CLK_TCK);
		#else
			clockTicks = CLOCKS_PER_SEC;
		#endif
		
		#ifndef LBTM_TYPE_DEFINED
			clockTicks = 1;
		#endif
	}
	
	
	/* Compute the real and CPU elapsed times to return */
	{
		#ifdef GEN_UNIX
			*realTimeIntervalSecs = (double)(stopTime.wallClockTime - startTime->wallClockTime) 
										/ clockTicks;
			*cpuTimeIntervalSecs = (double)
				((stopTime.cpuClockTime.tms_utime + stopTime.cpuClockTime.tms_stime) -
					(startTime->cpuClockTime.tms_utime + startTime->cpuClockTime.tms_stime)) 
				/ clockTicks;
		#elif DARWIN
			*realTimeIntervalSecs = (double)(stopTime.wallClockTime - startTime->wallClockTime) 
										/ clockTicks;
			*cpuTimeIntervalSecs = (double)
				((stopTime.cpuClockTime.tms_utime + stopTime.cpuClockTime.tms_stime) -
					(startTime->cpuClockTime.tms_utime + startTime->cpuClockTime.tms_stime)) 
				/ clockTicks;
		#elif OSX
			*realTimeIntervalSecs = (double)(stopTime.wallClockTime - startTime->wallClockTime) 
										/ clockTicks;
			*cpuTimeIntervalSecs = (double)
				((stopTime.cpuClockTime.tms_utime + stopTime.cpuClockTime.tms_stime) -
					(startTime->cpuClockTime.tms_utime + startTime->cpuClockTime.tms_stime)) 
				/ clockTicks;
		#elif COCOA
			*realTimeIntervalSecs = (double)(stopTime.wallClockTime - startTime->wallClockTime) 
										/ clockTicks;
			*cpuTimeIntervalSecs = (double)
				((stopTime.cpuClockTime.tms_utime + stopTime.cpuClockTime.tms_stime) -
					(startTime->cpuClockTime.tms_utime + startTime->cpuClockTime.tms_stime)) 
				/ clockTicks;
		#else
			*realTimeIntervalSecs = (double)(stopTime.wallClockTime - startTime->wallClockTime) 
										/ clockTicks;
			*cpuTimeIntervalSecs = (double)(stopTime.cpuClockTime - startTime->cpuClockTime) 
										/ clockTicks;
		#endif
		
		#ifndef LBTM_TYPE_DEFINED
			*realTimeIntervalSecs = 0.0;
			*cpuTimeIntervalSecs = 0.0;
		#endif
	}
	
	
	/* Determine the function return value */
	#ifndef LBTM_TYPE_DEFINED
		result = false;
	#else
		result = true;
	#endif
	
	
	return (result);
}
