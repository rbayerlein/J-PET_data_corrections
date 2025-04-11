/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 2006-2012 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		addrandUsr.c
*			Revision Number:	2.2
*			Date last revised:	1 February 2012
*			Programmer:			Robert Harrison
*			Date Originated:	17 February 2006
*
*			Module Overview:	This module is provided for individuals to
*								alter the processing of the decays in a time
*								window during addrandoms.c.  For instance,
*								this would be a place where the user could
*								implement changes to the way triples are
*								handled (SimSET rejects them by default) or
*								system specific detector/electronics 
*								algorithms.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:
*
*			Global variables defined:
*					addrandUsrParamsFPtr
*					addrandUsrModDecays1FPtr
*					addrandUsrModDecays2FPtr
*					addrandUsrTerminateFPtr
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
*			Revision date:		1 February 2012
*
*			Revision description:	Changed form of user functions to pointers
*
*********************************************************************************/

#define ADD_RANDOMS_USER

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"
#include "LbHeader.h"
#include "LbTiming.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhgMath.h"
#include "PhoHFile.h"
#include "PhgHdr.h"
#include "ProdTbl.h"
#include "PhoTrk.h"
#include "SubObj.h"
#include "EmisList.h"
#include "phg.h"
#include "Collimator.h"
#include "Detector.h"
#include "PhgBin.h"
#include "addrandoms.h"
#include "addrandUsr.h"


/* #include files with your user function software below.  We include the file for containing
 the triples processing user function example. */
#include "addrandUsr_procTriples.c"

/* Point to the user functions you want to use here.  Each user function is pointed to by 
	one of the variables below.  To *not* use a user function set the variable 
	to NULL, as the  commented-out line for each pointer shows below.  To use one, set
	the function pointer to the desired function.  Below we set the pointers to functions in
	the file addrandUsrSample.c, which is #included above. */
addrandUsrParamsFType		*addrandUsrInitializeFPtr = procTriplesInitialize;

addrandUsrDets1FType		*addrandUsrModDecays1FPtr = procTriplesProcess;

addrandUsrDets2FType		*addrandUsrModDecays2FPtr = NULL;

addrandUsrParamsFType		*addrandUsrTerminateFPtr = procTriplesTerminate;


#undef ADD_RANDOMS_USER
