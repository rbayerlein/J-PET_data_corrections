/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 2003-2006 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		addrandoms.h
*			Revision Number:	1.0
*			Date last revised:	20 February 2006
*			Programmer:			Robert Harrison
*			Date Originated:	17 February 2006
*
*			Module Overview:	Definitions for addrandoms.c.
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
#ifndef ADD_RAND_HDR
#define ADD_RAND_HDR

#ifdef ADD_RANDOMS_MAIN
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */
#define	ADRAND_MaxTWDecays 10					/* Maximum number of decays allowed in a time window -
												 as triples are rare, 10 should be enough */


/* TYPES */
typedef struct {
	PHG_Decay			decay;					/* decay info */
	LbUsFourByte		numBluePhotons;			/* number of blue photons in this decay */
	LbUsFourByte		numPinkPhotons;			/* number of pink photons in this decay */
	PHG_DetectedPhoton	bluePhotons[PHG_MAX_DETECTED_PHOTONS];
												/* the blue photons for this decay */
	double				blueDetectionTime[PHG_MAX_DETECTED_PHOTONS];
												/* time detected for each blue photon */
	PHG_DetectedPhoton	pinkPhotons[PHG_MAX_DETECTED_PHOTONS];
												/* the pink photons for this decay */
	double				pinkDetectionTime[PHG_MAX_DETECTED_PHOTONS];
												/* time detected for each pink photon */
} timeWindowDecayTy;

typedef struct {
	LbUsFourByte		numDecays;				/* number decays in this time window */
	double				lastDetectionTime;		/* last time at which a photon was detected */
	timeWindowDecayTy	decays[ADRAND_MaxTWDecays];
												/* the decays in this window */
} timeWindowDetectionsTy;



/* GLOBALS */


/* PROTOTYPES */

#undef LOCALE
#endif /* ADD_RAND_HDR */
