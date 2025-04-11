/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1998-1999 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		DetCylinder.h
*			Revision Number:	1.0
*			Date last revised:	April 21, 1999
*			Programmer:			Steven Vannoy
*			Date Originated:	January 20, 1998
*
*			Module Overview:	Definitions for DetCylinder.c.
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
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*
*			Revision description:
*						- support for DetGeometric
*
*********************************************************************************/
#ifndef DETECTOR_CYLINDER_HDR
#define DETECTOR_CYLINDER_HDR


#include "Detector.h"


#ifdef DETECTOR_CYLINDER
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */

/* TYPES */
/* GLOBALS */

/* PROTOTYPES */
Boolean DetCylTrack(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *photonPtr);

void DetCylinder(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		DetectedPhotonsTy *detPhotonsPtr);


void DetCylInitCylinders(void);
double	DetCylGtTruncatedFreePaths(PHG_TrackingPhoton *photonPtr, Boolean firstTime,
			LbUsFourByte curRingIndex, LbUsFourByte curLayerIndex, double *weight);
Boolean DetCylProjectToDetector(PHG_TrackingPhoton *photonPtr, 
			LbFourByte *ringNum, LbFourByte *detNum);
void DetCylFindNextInteraction(PHG_TrackingPhoton *photonPtr, Boolean firstTime, 
								LbFourByte *curRingNum, LbFourByte *curDetNum, 
								detEn_ActionTy *actionType, double *fpToGo, 
								LbUsFourByte *detMaterial, Boolean *isActive);

#undef LOCALE
#endif /* DETECTOR_CYLINDER_HDR */
