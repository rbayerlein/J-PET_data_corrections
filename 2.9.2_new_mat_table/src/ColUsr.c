/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1994-2013 Department of Radiology                	 *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		ColUsr.c
*			Revision Number:	2.2
*			Date last revised:	26 June 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	11 July 1994
*
*			Module Overview:	This module is provided for individuals to
*								perform specific collimation operations on photons
*								that have reached the collimator.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:	none
*
*			Global variables defined:		
*					ColUsrInitializeFPtr;
*					ColUsrModPETPhotonsFPtr;
*					ColUsrModSPECTPhotonsFPtr;
*					ColUsrTerminateFPtr;
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
*			Revision date:		26 June 2013
*
*			Revision description:	Changed form of user functions to pointers
*
*********************************************************************************/

#define COL_USER


#include <stdio.h>


#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "PhgMath.h"
#include "ColUsr.h"
#include "CylPos.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "PhoHFile.h"
#include "phg.h"
#include "PhgBin.h"


/* Global variables */
/*
##
##
	Declare your globals within the ColUsrVars struct.  This will prevent any
		problems with naming and scope.
##
##
*/
struct {
	/* Your fields of choice go here */
	int		dummyVar;
} ColUsrVars;


/* Local Prototypes
 ##
 ##
   Declare your local prototypes here as well - in particular, if you change the
    function pointers below (e.g., ColUsrInitializeFPtr) to point to your own 
    functions, they would need prototypes here.
 
   We give the following examples of prototypes for the unused example
    functions shown at the end of the file - these show the parameters that would 
    be passed to the user functions if the pointers below were enabled:
*/

void			ColUsrInitialize(ColRunTimeParamsTy *ColRunTimeParamsPtr);
Boolean			ColUsrPETPhotons(ColRunTimeParamsTy *ColRunTimeParamsPtr, PHG_Decay *decayPtr,
					PHG_TrackingPhoton *photonPtr);
Boolean			ColUsrSPECTPhotons(ColRunTimeParamsTy *ColRunTimeParamsPtr,
					PHG_Decay *decayPtr,
					PHG_TrackingPhoton *photonPtr);
void			ColUsrTerminate(ColRunTimeParamsTy *ColRunTimeParamsPtr);

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
ColUsrParamsFType			*ColUsrInitializeFPtr = NULL;
/*ColUsrParamsFType			*ColUsrInitializeFPtr = ColUsrInitialize;*/

ColUsrTrackingFType			*ColUsrModPETPhotonsFPtr = NULL;
/*ColUsrTrackingFType		*ColUsrModPETPhotonsFPtr = ColUsrPETPhotons;*/

ColUsrTrackingFType			*ColUsrModSPECTPhotonsFPtr = NULL;
/*ColUsrTrackingFType		*ColUsrModSPECTPhotonsFPtr = ColUsrSPECTPhotons;*/

ColUsrParamsFType			*ColUsrTerminateFPtr = NULL;
/*ColUsrParamsFType			*ColUsrTerminateFPtr = ColUsrTerminate;*/

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
 *			Name:			ColUsrInitialize
 *
 *			Summary:		Initialize the binning module.
 *
 *			Arguments:
 *				ColRunTimeParamsTy	*ColRunTimeParamsPtr - This is a ptr to the actual
 *											data structure, think very
 *											carefully before you modify
 *											any of the fields. It is
 *											assumed you will use it to
 *											set your variables up for
 *											later processing.
 *
 *
 *			Function return: None.
 *
 *********************************************************************************/

void ColUsrInitialize(ColRunTimeParamsTy	*ColRunTimeParamsPtr)

{
	/* Do your initialization here */
	
	if (ColRunTimeParamsPtr) {};	/* Removes unused parameter compiler warning */
	
}


/*********************************************************************************
 *
 *			Name:			ColUsrPETPhotons
 *
 *			Summary:		Update the binning images with the current batch of
 *							detected photons.
 *
 *			Arguments:
 *				ColRunTimeParamsTy	*ColRunTimeParamsPtr	- The collimator data.
 *				PHG_Decay			*decayPtr				- The decay that started the process.
 *				PHG_TrackingPhoton 	*photonPtr				- The photon detected.
 *
 *			Function return: True to accept the coincidence, False to reject it.
 *
 *********************************************************************************/

Boolean ColUsrPETPhotons(	ColRunTimeParamsTy *ColRunTimeParamsPtr,
                         PHG_Decay *decayPtr,
                         PHG_TrackingPhoton *photonPtr)

{
	Boolean	acceptCoincidence = false;		/* Should we accept the coincidence? */
	
	/* Remove unused parameter compiler warnings */
	if (ColRunTimeParamsPtr) {};
	if (decayPtr) {};
	if (photonPtr) {};
	
	
	/* Do your PET binning here */
	do {
        
		/* Set acceptance to true if we made it to here */
		acceptCoincidence = true;
	} while (false);
	
	return (acceptCoincidence);
}


/*********************************************************************************
 *
 *			Name:			ColUsrSPECTPhotons
 *
 *			Summary:		Update the binning images with the current batch of
 *							detected photons.
 *
 *			Arguments:
 *				ColRunTimeParamsTy	*ColRunTimeParamsPtr	- The collimator data.
 *				PHG_Decay			*decayPtr				- The decay that started the process.
 *				PHG_TrackingPhoton	*photonPtr				- The photon detected.
 *
 *			Function return: True to accept the coincidence, False to reject it.
 *
 *********************************************************************************/

Boolean ColUsrSPECTPhotons(	ColRunTimeParamsTy *ColRunTimeParamsPtr,
                           PHG_Decay *decayPtr,
                           PHG_TrackingPhoton *photonPtr)

{
	Boolean	acceptPhoton = false;		/* Should we accept the photon? */
	
	/* Remove unused parameter compiler warnings */
	if (ColRunTimeParamsPtr) {};
	if (decayPtr) {};
	if (photonPtr) {};
	
	
	/* Do your SPECT collimation here */
	do {
        
		/* Set acceptance to true if we made it to here */
		acceptPhoton = true;
	} while (false);
	
	return (acceptPhoton);
}


/*********************************************************************************
 *
 *			Name:			ColUsrTerminate
 *
 *			Summary:		Process the bin buffers and terminate the module.
 *
 *			Arguments:
 *				ColRunTimeParamsTy	*ColRunTimeParamsPtr - This is a ptr to the actual
 *											data structure, think very
 *											carefully before you modify
 *											any of the fields. It is
 *											assumed you will want to
 *											change the image data,
 *											and it may be ok to change
 *											other fields - but it is
 *											not advised.
 *
 *			Function return: None.
 *
 *********************************************************************************/

void ColUsrTerminate(ColRunTimeParamsTy	*ColRunTimeParamsPtr)

{
	if (ColRunTimeParamsPtr) {};	/* Removes unused parameter compiler warning */
	
	
	/* Do your finishing here */
	
	/* NOTE that this routine is called before PhgBinTerminate so 
     any modification to the image data here will be reflected
     in the output files.
     */
}


#undef COL_USR
