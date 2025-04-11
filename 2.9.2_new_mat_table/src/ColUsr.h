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
*			Module Name:		ColUsr.h
*			Revision Number:	2.0
*			Date last revised:	23 January 2012
*			Programmer:			Steven Vannoy
*			Date Originated:	11 July 1994
*
*			Module Overview:	Definitions for ColUsr.c.
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
*			Revision date:		23 January 2012
*
*			Revision description:	Changed form of user functions to pointers
*
*********************************************************************************/

#ifndef COL_USR_HDR
#define COL_USR_HDR

#ifdef COL_USR
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */

/* DEFINITIONS */
typedef 	void ColUsrParamsFType(ColRunTimeParamsTy	*colParamsPtr);
typedef 	Boolean ColUsrTrackingFType(ColRunTimeParamsTy	*colParamsPtr,
										PHG_Decay			*decayPtr,
										PHG_TrackingPhoton	*photonPtr);

/* GLOBALS */
extern ColUsrParamsFType		*ColUsrInitializeFPtr;
extern ColUsrTrackingFType		*ColUsrModPETPhotonsFPtr;
extern ColUsrTrackingFType		*ColUsrModSPECTPhotonsFPtr;
extern ColUsrParamsFType		*ColUsrTerminateFPtr;

/* PROTOTYPES */

#undef LOCALE
#endif /* COL_USR_HDR */
