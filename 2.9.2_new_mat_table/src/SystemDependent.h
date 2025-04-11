/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1994-2012 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		SystemDependent.h
*			Revision Number:	1.4
*			Date last revised:	19 December 2012
*			Programmer:			Robert Harrison
*			Date Originated:	8 Dec 2001
*
*			Module Overview:	An attempt to consolidate the OS-specific
*								library includes
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:
*
*			Global variables defined:		none
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		
*
*			Revision date:		
*
*			Revision description:
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		20 November 2003
*
*			Revision description:	Added changes for Mac CW 9 compatibility:
*										Removed obsolete #includes
*
*********************************************************************************/

#include <errno.h>
#include <stdlib.h>

#ifdef __MWERKS__
	#include	<fp.h>
	#include 	<Types.h>
	#include 	<ctype.h>
	#include	<Errors.h>
	#include	<float.h>
	#include	<limits.h>
	#include 	<String.h>
	#include 	<math.h>
	#define		MAXFLOAT	FLT_MAX
	#include 	<SIOUX.h>
	#include	<console.h>
#endif

/* old metroworks definition
#ifdef MWC
	#include 	<Math.h>
	#include 	<fp.h>
	#define		MAXFLOAT	FLT_MAX
#endif
*/

#ifdef DGUX 
	#include	<math.h>
	#include 	<values.h>
    #include 	<sys/m88kbcs.h>
#elif MPW
	#include	<Types.h>
	#include 	<ctype.h>
	#include	<ErrMgr.h>
	#include	<CursorCtl.h>
	#include	<Errors.h>
#elif __USE_TWRP__
	#include 	"TWrpDoc.h"
	#include	"TWrapper.h"
	#include	"TWrpLib.h"
	#define 	APP_DEBUG
#elif defined GEN_UNIX 
	#include 	<values.h>
	#include 	<sys/types.h>
	#include 	<sys/stat.h>
	#include 	<sys/times.h>
	#include 	<math.h>
	#include	<float.h>
	#include 	<unistd.h>

	#ifdef  LINUX
		#define		MAXFLOAT	FLT_MAX
	#endif

#elif defined AOS_VS
	#include 	<math.h>
#elif defined WINNT
	#include 	<string.h>
	#include	<stdlib.h>
	#include 	<math.h>
	#include	<assert.h>
//	#include 	<fp.h>
	#include	<ctype.h> // need include this inorder to use tolower()
	#define		MAXFLOAT	FLT_MAX
	/* windows doesn't have a built-in strtoll, copysign, so we #define them here... */
	# define strtoll _strtoi64
	#define copysign(val, signval) ( ((signval) < 0) ? (-fabs(val)) : (fasbs(val)) )
#elif defined DARWIN
	#include 	<unistd.h>
	#include 	<limits.h>
	#include 	<math.h>
	#include 	<float.h>
	#include	<sys/times.h>
#elif defined OSX
	#include 	<unistd.h>
	#include 	<limits.h>
	#include 	<math.h>
	#include 	<float.h>
	#include	<sys/times.h>
#elif defined COCOA
	#include 	<unistd.h>
	#include 	<limits.h>
	#include 	<math.h>
	#include 	<float.h>
	#include	<sys/times.h>
#endif


/* old sun definition
#ifdef SUN_OS
	#include 	<unistd.h>
#endif
*/

