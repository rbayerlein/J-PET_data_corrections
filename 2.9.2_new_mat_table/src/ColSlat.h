/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1999 -1999 Department of Radiology                 *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		ColSlat.h
*			Revision Number:	1.0
*			Date last revised:	1 July 1999
*			Programmer:			Steven Vannoy
*			Date Originated:	1 July 19996
*
*			Module Overview:	Definitions for ColSlat.c.
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
*********************************************************************************/
#ifndef COLLIMATOR_SLAT_HDR
#define COLLIMATOR_SLAT_HDR

#ifdef COLLIMATOR_SLAT
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */

/* TYPES */
/* GLOBALS */

/* PROTOTYPES */
void	ColSlatSPECT(PHG_Decay *decayPtr,
			PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
			CollimatedPhotonsTy *colPhotonsPtr);
			
void ColSlatDualHeaded(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		CollimatedPhotonsTy *colPhotonsPtr);
void	ColSlatSetParamsFromDetector(double minAngle, double minZ, double maxZ);

#undef LOCALE
#endif /* COLLIMATOR_SLAT_HDR */
