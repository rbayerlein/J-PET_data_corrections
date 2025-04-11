/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 2006-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		addrandUsr.c
*			Revision Number:	2.4
*			Date last revised:	10 June 2013
*			Programmer:			Robert Harrison
*			Date Originated:	17 February 2006
*
*			Module Overview:	This module is provided for individuals to
*								alter the processing of the decays in a time
*								window during addrandoms.c.  For instance,
*								this would be a place where the user could
*								implement changes to the way triples are
*								handled (SimSET rejects them by default) or
*								add system-specific detector/electronics 
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
*			Revision date:		7 June 2013
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


/* Global variables */
/*
##
##
	Declare your globals within the addrandUsrVars struct below.  This will prevent any
		problems with naming and scope.
##
##
*/
struct {
	/* Your fields of choice go here */
	int		dummyVar;
} addrandUsrVars;


/* Local Prototypes
 ##
 ##
   Declare your local prototypes here as well - in particular, if you change the
    function pointers below (e.g., addrandUsrInitializeFPtr) to point to your own 
    functions, they would need prototypes here.
 
   We give the following examples of prototypes for the unused example
    functions shown at the end of the file - these show the parameters that would 
    be passed to the user functions if the pointers below were enabled:
*/

void addrandUsrInitialize(DetRunTimeParamsTy	*detParams);
void addrandUsrModDecays1(	timeWindowDetectionsTy	*timeWindowDetectionsPtr,
                         PhoHFileHkTy			*adrandRandomsFileHk,
                         LbUsFourByte			*numDecaysWritten,
                         LbUsFourByte			*numDecaysUnchanged,
                         LbUsFourByte			*numDecaysRandom,
                         LbUsFourByte			*numLostCorrectWindow,
                         LbUsFourByte			*numDecaysLostTriples);
Boolean addrandUsrModDecays2(	timeWindowDetectionsTy	*timeWindowDetectionsPtr,
                         PhoHFileHkTy			*adrandRandomsFileHk,
                         LbUsFourByte			*numDecaysWritten,
                         LbUsFourByte			*numDecaysUnchanged,
                         LbUsFourByte			*numDecaysRandom,
                         LbUsFourByte			*numLostCorrectWindow,
                         LbUsFourByte			*numDecaysLostTriples);
void addrandUsrTerminate(DetRunTimeParamsTy	*detParams);

/*
 ##
 ##
*/


/* Point to the user functions you want to use here.  Each user function is pointed to by
	one of the variables below.  To *not* use a user function (the default) set the variable 
	to NULL (as below).  To use one, change the NULL to the name of your desired user function, 
	as in the commented examples below each variable.  Prototypes for the sample functions in  
	the comments are supplied in the comment above as templates for your functions.  
   For more information on how to write user functions go to the user guide section 
    'User Functions' on SimSET's web page.  There are also examples of user functions
    given in the SimSET subdirectory samples/userFuncExamples.  These are described
    on the web page.
*/
addrandUsrParamsFType		*addrandUsrInitializeFPtr = NULL;
/*addrandUsrParamsFType		*addrandUsrInitializeFPtr = addrandUsrInitialize;*/

addrandUsrDets1FType		*addrandUsrModDecays1FPtr = NULL;
/*addrandUsrDets1FType		*addrandUsrModDecays1FPtr = addrandUsrModDecays1;*/

addrandUsrDets2FType		*addrandUsrModDecays2FPtr = NULL;
/*addrandUsrDets2FType		*addrandUsrModDecays2FPtr = addrandUsrModDecays2;*/

addrandUsrParamsFType		*addrandUsrTerminateFPtr = NULL;
/*addrandUsrParamsFType		*addrandUsrTerminateFPtr = addrandUsrTerminate;*/

/*
 ##
 ##
 
  The rest of the file gives the sample functions corresponding to the above local prototypes
    and commented-out alternate function pointers above.

  The functions below can be ignored - they, and their prototypes above, are here only to give 
  	users an idea of how to go about implementing user functions.  There is much more 
  	information on these functions on the web pages.
 ##
 ##
*/


/*********************************************************************************
 *
 *		Name:			addrandUsrInitialize
 *
 *		Summary:	Initialize any special user-created structures.
 *				For instance, the user could create an extra list-mode
 *				file here.
 *
 *		Arguments:
 *			DetRunTimeParamsTy	*detParams
 *						This is a ptr to the detector parameters. It
 *						is assumed that you will treat this as
 *						a read-only structure. Modifications should
 *						only be made after a careful study of
 *						the parameters involved.
 *
 *		Function return: None.
 *
 *********************************************************************************/
void addrandUsrInitialize(DetRunTimeParamsTy	*detParams)

{
	/* Do your initialization here */
	
	if (detParams) {};		/* Removes unused parameter compiler warning */
	
}


/*********************************************************************************
 *
 *		Name:		addrandUsrModDecays1
 *
 *		Summary:	This routine is called before the photons/decays are
 *				time-windowed.  All the decays within a time-window of each
 *				other (or actually longer: the window is extended each
 *				time another photon arrives) are passed in the
 *				timeWindowDetectionsTy structure.  The photons and decays
 *				may be altered/deleted as desired, e.g. to account for
 *				deadtime, different triples handling than SimSET's
 *				(currently all triples are deleted, even if they occur
 *				in separate rings of the tomograph).  Make sure you
 *				understand the data structure first!
 *
 *
 *		Arguments:
 *			timeWindowDetectionsTy	*timeWindowDetections -
 *					All the decays in the current time window
 *
 *		Function return: None
 *
 *********************************************************************************/
void addrandUsrModDecays1(	timeWindowDetectionsTy	*timeWindowDetectionsPtr,
							PhoHFileHkTy			*adrandRandomsFileHk,
							LbUsFourByte			*numDecaysWritten,
							LbUsFourByte			*numDecaysUnchanged,
							LbUsFourByte			*numDecaysRandom,
							LbUsFourByte			*numLostCorrectWindow,
							LbUsFourByte			*numDecaysLostTriples)

{
	if (timeWindowDetectionsPtr) {};	/* Removes unused parameter compiler warning */
	if (adrandRandomsFileHk) {};		/* Removes unused parameter compiler warning */
	if (numDecaysWritten) {};			/* Removes unused parameter compiler warning */
	if (numDecaysUnchanged) {};			/* Removes unused parameter compiler warning */
	if (numDecaysRandom) {};			/* Removes unused parameter compiler warning */
	if (numLostCorrectWindow) {};		/* Removes unused parameter compiler warning */
	if (numDecaysLostTriples) {};		/* Removes unused parameter compiler warning */
	
	
	/* Make your alterations here */
	do {
		
	} while (false);
	
	return;
}


/*********************************************************************************
 *
 *			Name:			addrandUsrModDecays2
 *
 *		Summary:	This routine is called after the photons/decays are
 *				time-windowed.  Only the decays that are actually going
 *				to be written out will be passed to this routine - many
 *				of the decays/time windows passed to addrandUsrModDecays1
 *				will not generate a call to this function.  This
 *				function gives the user a chance to change/reject the
 *				coincidences that are getting written out.
 *
 *		Arguments:
 *			timeWindowDetectionsTy	*timeWindowDetections -
 *					All the decays in the current time window
 *
 *		Function return: Boolean
 *			True to accept the coincidence, False to reject it.
 *
 *********************************************************************************/
Boolean addrandUsrModDecays2(	timeWindowDetectionsTy	*timeWindowDetectionsPtr,
								PhoHFileHkTy			*adrandRandomsFileHk,
								LbUsFourByte			*numDecaysWritten,
								LbUsFourByte			*numDecaysUnchanged,
								LbUsFourByte			*numDecaysRandom,
								LbUsFourByte			*numLostCorrectWindow,
								LbUsFourByte			*numDecaysLostTriples)

{
	Boolean	acceptCoincidence = false;	/* Should coincidence be added to history list? */
	
	if (timeWindowDetectionsPtr) {};	/* Removes unused parameter compiler warning */
	if (adrandRandomsFileHk) {};		/* Removes unused parameter compiler warning */
	if (numDecaysWritten) {};			/* Removes unused parameter compiler warning */
	if (numDecaysUnchanged) {};			/* Removes unused parameter compiler warning */
	if (numDecaysRandom) {};			/* Removes unused parameter compiler warning */
	if (numLostCorrectWindow) {};		/* Removes unused parameter compiler warning */
	if (numDecaysLostTriples) {};		/* Removes unused parameter compiler warning */
	
	
	/* Make your alterations here */
	do {
        
		/* Set acceptance to true if we made it to here */
		acceptCoincidence = true;
        
	} while (false);
	
	return (acceptCoincidence);
}


/*********************************************************************************
 *
 *		Name:		addrandUsrTerminate
 *
 *		Summary:	Do any termination required by the user functions (e.g.
 *				write to/close any user files).
 *
 *		Arguments:
 *			DetRunTimeParamsTy	*detParams
 *						This is a ptr to the detector parameters. It
 *						is assumed that you will treat this as
 *						a read-only structure. Modifications should
 *						only be made after a careful study of
 *						the parameters involved.
 *
 *		Function return: None.
 *
 *********************************************************************************/
void addrandUsrTerminate(DetRunTimeParamsTy	*detParams)

{
	/* Do your finishing here */
	
	if (detParams) {};		/* Removes unused parameter compiler warning */
	
}


#undef ADD_RANDOMS_USER
