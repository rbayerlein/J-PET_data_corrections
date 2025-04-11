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
*			Module Name:		DetUsr.h
*			Revision Number:	2.0
*			Date last revised:	26 January 2012
*			Programmer:			Steven Vannoy
*			Date Originated:	11 July 1994
*
*			Module Overview:	Definitions for DetUsr.c.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:	none
*
*			Global variables defined:		
*					DetUsrInitializeFPtr;
*					DetUsrModPETPhotonsFPtr;
*					DetUsrModSPECTPhotonsFPtr;
*					DetUsrTerminateFPtr;
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
*			Revision date:		26 January 2012
*
*			Revision description:	Changed form of user functions to pointers
*
*********************************************************************************/

#ifndef DET_USR_HDR
#define DET_USR_HDR

#ifdef DET_USR
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */

/* DEFINITIONS */
typedef 	void DetUsrParamsFType(DetRunTimeParamsTy	*detParamsPtr);
typedef 	Boolean DetUsrTrackingFType(DetRunTimeParamsTy	*detParamsPtr,
										PHG_Decay			*decayPtr,
										PHG_TrackingPhoton	*photonPtr);

/* GLOBALS */
extern DetUsrParamsFType		*DetUsrInitializeFPtr;
extern DetUsrTrackingFType		*DetUsrModPETPhotonsFPtr;
extern DetUsrTrackingFType		*DetUsrModSPECTPhotonsFPtr;
extern DetUsrParamsFType		*DetUsrTerminateFPtr;

/* PROTOTYPES */

#undef LOCALE
#endif /* DET_USR_HDR */
