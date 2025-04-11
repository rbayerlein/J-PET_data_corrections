/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1993-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		PhgUsrBin.c
*			Revision Number:	2.1
*			Date last revised:	26 June 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	28 July 1993
*
*			Module Overview:	This module is provided for individuals to
*								perform specific binning operations on PHG
*								data independent of the "built-in" functions
*								of the PHG.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:		none
*
*			Global variables defined:		
*					BinUsrInitializeFPtr;
*					BinUsrModPETPhotonsFPtr;
*					BinUsrModPETPhotons2FPtr;
*					BinUsrModSPECTPhotonsFPtr;
*					BinUsrModSPECTPhotons2FPtr;
*					BinUsrTerminateFPtr;
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
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*			Revision description:
*						- binning by random-state and crystal number.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:	Added support for time-of-flight
*
*********************************************************************************/

#define PHG_USR_BIN


#include <stdio.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbParamFile.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "PhgMath.h"
#include "ColTypes.h"
#include "DetTypes.h"
#include "ColParams.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhoHFile.h"
#include "phg.h"
#include "PhgUsrBin.h"
#include "PhgBin.h"


/* Global variables */
/*
##
##
	Declare your globals within the PhgUsrBinVars struct.  This will prevent any
		problems with naming and scope.
##
##
*/
struct {
	/* Your fields of choice go here */
	int		dummyVar;
} PhgUsrBinVars;


/* Local Prototypes
 ##
 ##
   Declare your local prototypes here as well - in particular, if you change the
    function pointers below (e.g., BinUsrInitializeFPtr) to point to your own 
    functions, they would need prototypes here.
 
   We give the following examples of prototypes for the unused example
    functions shown at the end of the file - these show the parameters that would 
    be passed to the user functions if the pointers below were enabled:
*/

void	PhgUsrBinInitialize(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData);
Boolean	PhgUsrBinPETPhotons(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData,
					PHG_Decay *decay,
					PHG_TrackingPhoton *bluePhoton,
					PHG_TrackingPhoton *pinkPhoton);
Boolean	PhgUsrBinSPECTPhotons(PHG_BinParamsTy *binParams ,PHG_BinDataTy *binData,
					PHG_Decay *decay,
					PHG_TrackingPhoton *photon);
void	PhgUsrBinTerminate(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData);
Boolean	PhgUsrBinSPECTPhotons2(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData,
					PHG_Decay *decay,
					PHG_TrackingPhoton *photon,
					LbUsFourByte *angleIndex,
					LbUsFourByte *distIndex,
					LbUsFourByte *scatterIndex,
					LbUsFourByte *energyIndex,
					LbUsFourByte *crystalIndex,
					LbUsFourByte *zIndex,
					LbUsFourByte *imageIndex);
Boolean	PhgUsrBinPETPhotons2(
			PHG_BinParamsTy *binParams,
			PHG_BinDataTy *binData,
			PHG_Decay *decay,
			PHG_TrackingPhoton *bluePhoton,
			PHG_TrackingPhoton *pinkPhoton,
			LbUsFourByte *angleIndex,
			LbUsFourByte *distIndex,
			LbUsFourByte *tofIndex,
			LbUsFourByte *scatter1Index,
			LbUsFourByte *scatter2Index,
			LbUsFourByte *crystal1Index,
			LbUsFourByte *crystal2Index,
			LbUsFourByte *energy1Index,
			LbUsFourByte *energy2Index,
			LbUsFourByte *zDownIndex,
			LbUsFourByte *zUpIndex,
			LbUsFourByte *thetaIndex,
			LbUsFourByte *phiIndex,
			LbUsFourByte *xrIndex,
			LbUsFourByte *yrIndex,
			LbUsFourByte *imageIndex,
			double		 *coincidenceWt,
			double		 *coincidenceSqWt);

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
BinUsrParamsFType			*BinUsrInitializeFPtr = NULL;
/*BinUsrParamsFType			*BinUsrInitializeFPtr = PhgUsrBinInitialize;*/

BinUsrPETTrackingFType		*BinUsrModPETPhotonsFPtr = NULL;
/*BinUsrPETTrackingFType	*BinUsrModPETPhotonsFPtr = PhgUsrBinPETPhotons;*/

BinUsrPETTrackingFType2		*BinUsrModPETPhotonsF2Ptr = NULL;
/*BinUsrPETTrackingFType2	*BinUsrModPETPhotonsF2Ptr = PhgUsrBinPETPhotons2;*/

BinUsrSPECTTrackingFType	*BinUsrModSPECTPhotonsFPtr = NULL;
/*BinUsrSPECTTrackingFType	*BinUsrModSPECTPhotonsFPtr = PhgUsrBinSPECTPhotons;*/

BinUsrSPECTTrackingFType2	*BinUsrModSPECTPhotonsF2Ptr = NULL;
/*BinUsrSPECTTrackingFType2	*BinUsrModSPECTPhotonsF2Ptr = PhgUsrBinSPECTPhotons2;*/

BinUsrParamsFType			*BinUsrTerminateFPtr = NULL;
/*BinUsrParamsFType			*BinUsrTerminateFPtr = PhgUsrBinTerminate;*/

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
 *			Name:			PhgUsrBinInitialize
 *
 *			Summary:		Initialize the binning module.
 *
 *			Arguments:
 *				PHG_BinParamsTy	*binParams - This is a ptr to the parameters
 *												used by the binning module. It
 *												is primarily assumed that you
 *												will treat this as a read-only
 *												structure. Modifying any fields
 *												should only be done after a
 *												careful study of the binning
 *												process.
 *				PHG_BiDataTy	*binData	 - This is a ptr to structure that
 *												contains the histogram data.
 *
 *
 *			Function return: None.
 *
 *********************************************************************************/

void PhgUsrBinInitialize(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData)

{
	/* Do your initialization here */
	
	/* Remove unused parameter compiler warnings */
	if (binParams) {};
	if (binData) {};
	
}


/*********************************************************************************
 *
 *			Name:			PhgUsrBinPETPhotons
 *
 *			Summary:		This routine is called before the binning module processes
 *							the coincidence.  You may alter the photon's parameters to alter
 *							the way the binning module handles it.  A nice trick to
 *							circumvent the binning module but use it's data files is
 *							to use the binData information to update the images and
 *							then return 'false' so the binning module goes onto the next
 *							coincidence pair.  The updated image will be written to disc
 *							even though the binning module "thinks" the pair was skipped.
 *
 *			Arguments:
 *				PHG_BinParamsTy		*binParams		- The binning parameters.
 *				PHG_BinDataTy		*binData			- The binning histograms.
 *				PHG_Decay			*decay			- The decay that started the process.
 *				PHG_TrackingPhoton *bluePhoton		- The blue photon detected.
 *				PHG_TrackingPhoton *pinkPhoton		- The  pink photon detected.
 *
 *			Function return: True to accept the coincidence, False to reject it.
 *
 *********************************************************************************/

Boolean PhgUsrBinPETPhotons(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData,
                            PHG_Decay *decay,
                            PHG_TrackingPhoton *bluePhoton,
                            PHG_TrackingPhoton *pinkPhoton)

{
	Boolean	acceptCoincidence = false;		/* Should the coincidence be processed by the binning module? */
	
	/* Remove unused parameter compiler warnings */
	if (binParams) {};
	if (binData) {};
	if (decay) {};
	if (bluePhoton) {};
	if (pinkPhoton) {};
	
	
	/* Do your PET binning here */
	do {
        
		/* Set acceptance to true if we made it to here */
		acceptCoincidence = true;
	} while (false);
	
	return (acceptCoincidence);
}


/*********************************************************************************
 *
 *			Name:			PhgUsrBinPETPhotons2
 *
 *			Summary:		This routine is called after the coincidence has been processed
 *							by the binning module but just before the images are updated.
 *							The argument list is very long but everything is included to
 *							allow users to modify any aspect of the coincidence about to
 *							be binned at this point.
 *			Arguments:
 *				PHG_BinParamsTy		*binParams		- The binning parameters.
 *				PHG_BinDataTy		*binData			- The binning histograms.
 *				PHG_Decay			*decay			- The decay that started the process.
 *				PHG_TrackingPhoton *bluePhoton		- The blue photon detected.
 *				PHG_TrackingPhoton *pinkPhoton		- The  pink photon detected.
 *				LbUsFourByte		*angleIndex			- Index for angle bin
 *				LbUsFourByte		*distIndex			- Index for dist bin
 *				LbUsFourByte		*tofIndex			- Index for time-of-flight bin
 *				LbUsFourByte		*scatter1Index		- Index for scatter 1 bin
 *				LbUsFourByte		*scatter2Index		- Index for scatter 2 bin
 *				LbUsFourByte		*energy1Index		- Index for energy 1 bin
 *				LbUsFourByte		*energy2Index		- Index for energy 2 bin
 *				LbUsFourByte		*zDownIndex			- Z index for photon having min(blue.y,pink.y)
 *				LbUsFourByte		*zUpIndex			- Z index for photon having max(blue.y,pink.y)
 *				LbUsFourByte		*thetaIndex			- Index for theta in 3DRP
 *				LbUsFourByte		*phiIndex			- Index for PHI in 3DRP
 *				LbUsFourByte		*xrIndex			- Xr index in 3DRP
 *				LbUsFourByte		*yrIndex			- Yr index in 3DRP
 *				LbUsFourByte		*imageIndex			- Index for image
 *				double				*coincidenceWt		- The weight
 *				double				coincidenceSqWt		- The weight squared
 *
 *			Function return: True to accept the coincidence, False to reject it.
 *
 *********************************************************************************/

Boolean PhgUsrBinPETPhotons2(
                             PHG_BinParamsTy 	*binParams,
                             PHG_BinDataTy 		*binData,
                             PHG_Decay 			*decay,
                             PHG_TrackingPhoton	*bluePhoton,
                             PHG_TrackingPhoton	*pinkPhoton,
                             LbUsFourByte		*angleIndex,
                             LbUsFourByte		*distIndex,
                             LbUsFourByte		*tofIndex,
                             LbUsFourByte		*scatter1Index,
                             LbUsFourByte		*scatter2Index,
                             LbUsFourByte		*crystal1Index,
                             LbUsFourByte		*crystal2Index,
                             LbUsFourByte		*energy1Index,
                             LbUsFourByte		*energy2Index,
                             LbUsFourByte		*zDownIndex,
                             LbUsFourByte		*zUpIndex,
                             LbUsFourByte		*thetaIndex,
                             LbUsFourByte		*phiIndex,
                             LbUsFourByte		*xrIndex,
                             LbUsFourByte		*yrIndex,
                             LbUsFourByte		*imageIndex,
                             double				*coincidenceWt,
                             double				*coincidenceSqWt)

{
	Boolean	acceptCoincidence = false;		/* Should the coincidence be processed by the binning module? */
	
	/* Remove unused parameter compiler warnings */
	if (binParams) {};
	if (binData) {};
	if (decay) {};
	if (bluePhoton) {};
	if (pinkPhoton) {};
	if (angleIndex) {};
	if (distIndex) {};
	if (tofIndex) {};
	if (scatter1Index) {};
	if (scatter2Index) {};
	if (crystal1Index) {};
	if (crystal2Index) {};
	if (energy1Index) {};
	if (energy2Index) {};
	if (zDownIndex) {};
	if (zUpIndex) {};
	if (thetaIndex) {};
	if (phiIndex) {};
	if (xrIndex) {};
	if (yrIndex) {};
	if (imageIndex) {};
	if (coincidenceWt) {};
	if (coincidenceSqWt) {};
	
	
	/* Do your PET binning here */
	do {
        
		/* Set acceptance to true if we made it to here */
		acceptCoincidence = true;
	} while (false);
	
	return (acceptCoincidence);
}


/*********************************************************************************
 *
 *			Name:			PhgUsrBinSPECTPhotons
 *
 *			Summary:		This routine is called before the binning module processes
 *							the photon.  You may alter the photon's parameters to alter
 *							the way the binning module handles it.  A nice trick to
 *							circumvent the binning module but use it's data files is
 *							to use the binData information to update the images and
 *							then return 'false' so the binning module goes onto the next
 *							coincidence pair.  The updated image will be written to disc
 *							even though the binning module "thinks" the pair was skipped.
 *
 *			Arguments:
 *				PHG_BinParamsTy		*binParams		- The binning parameters.
 *				PHG_BinDataTy		*binData			- The binning histograms.
 *				PHG_Decay			*decay			- The decay that started the process.
 *				PHG_TrackingPhoton *photon			- The photon detected.
 *			Function return: True to accept the coincidence, False to reject it.
 *
 *********************************************************************************/

Boolean PhgUsrBinSPECTPhotons(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData,
                              PHG_Decay *decay, PHG_TrackingPhoton *photon)

{
	Boolean	acceptPhoton = false;		/* Should the photon be processed by the binning module? */
	
	/* Remove unused parameter compiler warnings */
	if (binParams) {};
	if (binData) {};
	if (decay) {};
	if (photon) {};
	
	
	/* Do your SPECT binning here */
	do {
        
		/* Set acceptance to true if we made it to here */
		acceptPhoton = true;
	} while (false);
	
	return (acceptPhoton);
}


/*********************************************************************************
 *
 *			Name:			PhgUsrBinSPECTPhotons2
 *
 *			Summary:		This routine is called after the photon has been processed
 *							by the binning module but just before the images are updated.
 *							The argument list is very long but everything is included to
 *							allow users to modify any aspect of the photon about to
 *							be binned at this point.
 *			Arguments:
 *				PHG_BinParamsTy		*binParams			- The binning parameters.
 *				PHG_BinDataTy		*binData			- The binning histograms.
 *				PHG_Decay			*decay				- The decay that started the process.
 *				PHG_TrackingPhoton 	*photon				- The blue photon detected.
 *				PHG_TrackingPhoton	 *pinkPhoton		- The  pink photon detected.
 *				LbUsFourByte		*angleIndex			- Index for angle bin
 *				LbUsFourByte		*distIndex			- Index for dist bin
 *				LbUsFourByte		*scatterIndex		- Index for scatter 1 bin
 *				LbUsFourByte		*energyIndex		- Index for energy 1 bin
 *				LbUsFourByte		*zIndex				- Z index for photon having min(blue.y,pink.y)
 *				LbUsFourByte		*imageIndex			- Index for image
 *
 *			Function return: True to accept the coincidence, False to reject it.
 *
 *********************************************************************************/

Boolean PhgUsrBinSPECTPhotons2(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData,
                               PHG_Decay *decay,
                               PHG_TrackingPhoton *photon,
                               LbUsFourByte *angleIndex,
                               LbUsFourByte *distIndex,
                               LbUsFourByte *scatterIndex,
                               LbUsFourByte *energyIndex,
                               LbUsFourByte *crystalIndex,
                               LbUsFourByte *zIndex,
                               LbUsFourByte *imageIndex)

{
	Boolean	acceptPhoton = false;		/* Should the photon be processed by the binning module? */
	
	/* Remove unused parameter compiler warnings */
	if (binParams) {};
	if (binData) {};
	if (decay) {};
	if (photon) {};
	if (angleIndex) {};
	if (distIndex) {};
	if (scatterIndex) {};
	if (energyIndex) {};
	if (crystalIndex) {};
	if (zIndex) {};
	if (imageIndex) {};
	
	
	/* Do your SPECT binning here */
	do {
        
		/* Set acceptance to true if we made it to here */
		acceptPhoton  = true;
	} while (false);
	
	return (acceptPhoton);
}


/*********************************************************************************
 *
 *			Name:			PhgUsrBinTerminate
 *
 *			Summary:		Process the bin buffers and terminate the module.
 *
 *			Arguments:
 *				PHG_BinParamsTy	*binParams - This is a ptr to the parameters
 *												used by the binning module. It
 *												is primarily assumed that you
 *												will treat this as a read-only
 *												structure. Modifying any fields
 *												should only be done after a 
 *												careful study of the binning
 *												process.
 *				PHG_BinDataTy	*binData	 - This is a ptr to structure that
 *													contains the histogram data.
 *													
 *			Function return: None.
 *
 *********************************************************************************/

void PhgUsrBinTerminate(PHG_BinParamsTy	*binParams, PHG_BinDataTy *binData)

{
	/* Do your finishing here */
	
	/* Remove unused parameter compiler warnings */
	if (binParams) {};
	if (binData) {};
	
	/* NOTE that this routine is called before PhgBinTerminate so 
     any modification to the image data here will be reflected
     in the output files.
     */
}


#undef PHG_USR_BIN
