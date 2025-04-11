/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1996-2006 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		DetPlanar.h
*			Revision Number:	1.5
*			Date last revised:	28 February 2006
*			Programmer:			Steven Vannoy
*			Date Originated:	5	December 1996
*
*			Module Overview:	Definitions for DetPlanar.c.
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
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*
*			Revision description:
*						- support for DetGeometric
*
*********************************************************************************/
#ifndef DETECTOR_PLANAR_HDR
#define DETECTOR_PLANAR_HDR

#ifdef DETECTOR_PLANAR
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */

/* TYPES */
/* GLOBALS */

/* PROTOTYPES */
void	DetPlanarSPECT(PHG_Decay *decayPtr,
			PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
			DetectedPhotonsTy *detPhotonsPtr);
			
void	DetDualHeaded(PHG_Decay *decayPtr,
			PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
			PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
			DetectedPhotonsTy *detPhotonsPtr);

void	DetGetRandomDetPosition(double *detPos);

void	DetGetDetAngle(PHG_Position *location, PHG_Direction *angle,
					double dPos, double *detPos);


void DetPlnrInitPhotons(PHG_Decay *decayPtr, PHG_TrackingPhoton *photonPtr);
double DetPlnrGtTruncatedFreePaths(PHG_TrackingPhoton *photonPtr, double *weight);
Boolean DetPlnrProjectToDetector(PHG_TrackingPhoton *photonPtr, 
			LbFourByte *ringNum, LbFourByte *detNum);
void DetPlnrFindNextInteraction(PHG_TrackingPhoton *photonPtr, 
								LbFourByte interactionNum, LbFourByte *curDetNum, 
								detEn_ActionTy *actionType, double *fpToGo, 
								LbUsFourByte *detMaterial, Boolean *isActive);

#undef LOCALE
#endif /* DETECTOR_PLANAR_HEADER */
