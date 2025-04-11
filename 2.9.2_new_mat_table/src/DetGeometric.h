/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1998-2008 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		DetGeometric.h
*			Revision Number:	1.1
*			Date last revised:	8 January 2008
*			Programmer:			Steven Gillispie, Steven Vannoy
*			Date Originated:	21 November 2005
*
*			Module Overview:	Definitions for DetGeometric.c.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:	
*
*			Global variables defined:		None
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
*********************************************************************************/

#ifndef DETECTOR_GEOMETRIC_HDR
#define DETECTOR_GEOMETRIC_HDR


#ifdef DETECTOR_GEOMETRIC
	#define	LOCALE
#else
	#define LOCALE	extern
#endif


/* CONSTANTS */

/* TYPES */

/* GLOBALS */

/* PROTOTYPES */
Boolean DetGeomTrack(DetEn_DetectorTypeTy detectorType, 
					PHG_Decay *decayPtr, PHG_TrackingPhoton *photonPtr);
void DetGeometric(DetEn_DetectorTypeTy detectorType, 
					PHG_Decay *decayPtr,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					DetectedPhotonsTy *detPhotonsPtr);
void DetGeomCompCentroid(PHG_TrackingPhoton *photonPtr, 
					PHG_Position *centroidPos, double *depositedEnergy, 
					LbUsFourByte interactionIndex, Boolean *interactsSet);

#undef LOCALE
#endif /* DETECTOR_GEOMETRIC_HDR */
