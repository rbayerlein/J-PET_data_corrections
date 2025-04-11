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
*			Module Name:		Detector.h
*			Revision Number:	1.4
*			Date last revised:	9 October 2012
*			Programmer:			Steven Vannoy
*			Date Originated:	11 July 1994
*
*			Module Overview:	Definitions for Detector.c.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:	
*
*			Global variables defined:		none
*
*********************************************************************************
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
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		1 December 2009
*
*			Revision description:	Added detHistFileHk to detDataTy.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		August 2006
*
*			Revision description:	Added DetIsBlock prototype.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		20 December 2005
*
*			Revision description:	Added new types to detEn_ActionTy
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:	Added support for time-of-flight blurring.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		February 2004
*
*			Revision description:	Added prototypes for the DetIs... functions.
*
*********************************************************************************/

#ifndef DETECTOR_HDR
#define DETECTOR_HDR

#ifdef DETECTOR
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */

/* TYPES */

/*	The following enumerate describes the types of actions to be
	taken as a result of projecting a photon to a new location.
*/
typedef enum { detEnAc_Null,
				detEnAc_Detect,
				detEnAc_Discard,
				detEnAc_Interact,
				detEnAc_Absorb,
				detEnAc_AxialCross,
				detEnAc_LayerCross, 
				detEnAc_ComptonScatter, 
				detEnAc_CohScatter
} detEn_ActionTy;

/* GLOBALS */

#define DET_STAT_BINS	32
#define DET_NUM_DEPTH_BINS		10					/*	Number of depth bins for local statistics */
#define	DET_RANGE_OF_DET_ANGLES	(PHGMATH_2PI)		/*  Detector goes all the way around */

LOCALE	char					detErrStr[1024];					/* Storage for creating error strings */

/* NOTE: These are globals only so they may be shared among detector modules */
typedef struct {
 double					detPlnDetectorDepth;					/*	Depth of detector in radial direction */
 double					detPlnTransLimit;						/* 	Limit in the transaxial direction of the planar detector */
 CylPosCylinderTy		detInBoundCyl;							/*	The inner bounding cylinder for the detector */

 double					detWeightAbsorbedBins[MAX_DET_INTERACTIONS];
 double					detWeightEscapedBins[MAX_DET_INTERACTIONS+1];
 double					detWeightAdjusted;
 LbUsEightByte			detNumReachedMaxInteractions;	/* Count of photons that were forced to absorb to to maximum interactions */

 PhoHFileHkTy			detHistFileHk;						/* The history file */
 PHG_TrackingPhoton		detDetectedBluePhotons[PHG_MAX_DETECTED_PHOTONS];
 PHG_TrackingPhoton		detDetectedPinkPhotons[PHG_MAX_DETECTED_PHOTONS];
 LbUsFourByte			detDetectedBluePhotonIndex;
 LbUsFourByte			detDetectedPinkPhotonIndex;
 CylPosCylinderTy		detOutBoundCyl;							/*	The outer bounding cylinder for the detector */
 double					detViewSize;							/*	Range of detector view, used for planar and dual headed planar */
 LbUsEightByte			detTotBluePhotons;					/* Counter for output report */
 LbUsEightByte			detTotPinkPhotons;					/* Counter for output report */
 LbUsEightByte			detTotAccBluePhotons;				/* Counter for output report */
 LbUsEightByte			detTotAccPinkPhotons;				/* Counter for output report */
 LbUsEightByte			detTotPhotonsDepositingEnergy;		/* Counter for output report */
 LbUsEightByte			detTotPhotonsAbsorbed;				/* Counter for output report */
 LbUsEightByte			detTotForcedAbsorptions;			/* Counter for output report */
 LbUsEightByte			detTotPhotonsPassingThrough;		/* Counter for output report */
 LbUsEightByte			detTotFirstTimeAbsorptions;			/* Counter for output report */
 LbUsEightByte			detTotReachingCrystal;				/* Counter for output report */
 double					detTotWtAbsorbed;
 double					detTotWtForcedAbsorbed;
 double					detTotWtFirstTimeAbsorbed;
 CylPosCylinderTy		detPlnrBigCylinder;
 double					detPlnrAngularCoverage;
 double					detPlnrDelta;						/* Size of detector positions */


#ifdef PHG_DEBUG
 LbUsEightByte 			detCylCountInteractions;
 LbUsEightByte 			detCylCountCohInteractions;
	Boolean				DetDiscardBlue;
	Boolean				DetDiscardPink;
#endif
} detDataTy;

LOCALE	detDataTy	detData[PHG_MAX_PARAM_FILES];

/* MACROS */
#define DET_DoEnergyBlur()	(DetRunTimeParams[DetCurParams].EnergyResolutionPercentage != -1.0)

#define DET_DoTofBlur()	(DetRunTimeParams[DetCurParams].PhotonTimeFWHM != 0.0)

/* PROTOTYPES */
double			DetGaussEnergyBlur(double energy);
double			DetGaussTimeBlur(double travelDistance);
void 			DetUpdateDetectedPhotonBlock(PHG_TrackingPhoton	*trackingPhotonPtr);
Boolean			DetInitialize(Boolean doHistory);
void			DetPrintParams(void);
void			DetPrintReport(void);
Boolean			DetIsSimplePet(void);
Boolean			DetIsSimpleSpect(void);
Boolean			DetIsUNCSpect(void);
Boolean			DetIsPlanar(void);
Boolean			DetIsCylindrical(void);
Boolean			DetIsBlock(void);
Boolean			DetIsDualHeaded(void);
double			DetGtOutsideRadius(void);
double			DetGtInsideRadius(void);
LbFourByte		DetGtNumViews(void);
void			DetPETPhotons(PHG_Decay *decay,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					DetectedPhotonsTy *colPhotonsPtr);
void			DetSPECTPhotons(PHG_Decay *decay,
					PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
					DetectedPhotonsTy *colPhotonsPtr);
void			DetTerminate(void);

#undef LOCALE
#endif /* DETECTOR_HEADER */
