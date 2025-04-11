/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1994-2011 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/

/*********************************************************************************
*
*			Module Name:		Collimator.h
*			Revision Number:	1.1
*			Date last revised:	8 December 2011
*			Programmer:			Steven Vannoy
*			Date Originated:	11 July 1994
*
*			Module Overview:	Definitions for Collimator.c.
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
*							- Eight-byte integer for number of decays support.
*							- changes that look forward to the possibility of
*							collimators with different radii for each layer and
*							non-circular collimators.
*
*********************************************************************************/
#ifndef COL_HDR
#define COL_HDR

#ifdef COLLIMATOR
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */

/* TYPES */
/* GLOBALS */
typedef	struct {
	Col_Slat_Slat_Ty			**colSlatSegs;						/* Segments for slat collimator */
	LbFourByte					*colSlatNumSegs;					/* How many segments per layer */
	 PhoHFileHkTy				colHistFileHk;						/* The history file */
	 PHG_TrackingPhoton			colDetectedBluePhotons[PHG_MAX_DETECTED_PHOTONS];
	 PHG_TrackingPhoton			colDetectedPinkPhotons[PHG_MAX_DETECTED_PHOTONS];
	 LbUsFourByte				colDetectedBluePhotonIndex;
	 LbUsFourByte				colDetectedPinkPhotonIndex;
	 LbFourByte					colCurLayer;						/* Current collimator layer */
	 LbUsFourByte				colCurSeg;							/* Current segment (if we are doing segmented collimators */
	 Col_Monte_Carlo_PET_SegTy	*colCurMCPETSeg;					/* The current MCPET collimator segment (if we are doing that type of collimator) */
	 double						colCurInnerRadius;					/* The inner radius of current layer of MCPET collimator */
	 double						colCurOuterRadius;					/* The outer radius of current layer of MCPET collimator */
	 double						colInnermostRadius;					/* The inner most radius of MCPET collimator */
	 double						colOutermostRadius;					/* The outer most radius of MCPET collimator */
	 CylPosCylinderTy			colInBoundCyl;						/* The inner bounding cylinder for the collimator */
	 CylPosCylinderTy			colOutBoundCyl;						/* The outer bounding cylinder for the collimator */
					/* Output report variables */
	 LbUsEightByte				colTotBluePhotons;					/* Counter for total incoming blue photons */
	 LbUsEightByte				colTotPinkPhotons;					/* Counter for total incoming pink photons */
	 LbUsEightByte				colTotReachingCollimator;			/* For slat collimators, the number of photons
	 																	that project to the collimator */
	 LbUsEightByte				colTotAccBluePhotons;				/* Counter for total accepted blue photons */
	 LbUsEightByte				colTotAccPinkPhotons;				/* Counter for total accepted pink photons */
	double						colTotScatWtPassedThroughCollimator;	/* Total weight of object-scattered photons
													that passed through the collimator without scattering again */
	double						colTotPrimWtPassedThroughCollimator;	/* Total weight of primary photons
													that passed through the collimator without scattering */
	 double						colTotBluePhotonWt;					/* Total blue photon weight entering collimator */
	 double						colTotPinkPhotonWt;					/* Total pink photon weight entering collimator */
	 double						colTotInCoincWt;					/* Total coincidence weight entering collimator */
	 double						colTotAccBluePhotonWt;				/* Total blue photon weight accepted by collimator */
	 double						colTotAccPinkPhotonWt;				/* Total pink photon weight accepted by collimator */
	 double						colTotAccCoincWt;					/* Total coincidence weight accepted by collimator */
	 double						colTotRejBluePhotonWt;				/* Total blue photon weight rejected by collimator */
	 double						colTotRejPinkPhotonWt;				/* Total pink photon weight rejected by collimator */
	 double						colAccPrimWeightSum;				/* Total weight primary photons accepted by the UNC Collimator */
	 double						colAccScatWeightSum;				/* Total weight scatter photons accepted by the UNC Collimator */
	 double						k1y;								/* calculated UNCCollimator constant */
	 double						k2y;								/* calculated UNCCollimator constant */
	 double						k3y;								/* calculated UNCCollimator constant */
	 double						k1z;								/* calculated UNCCollimator constant */
	 double						k2z;								/* calculated UNCCollimator constant */
	 double						k3z;								/* calculated UNCCollimator constant */
	 double						colDistOriginToColBack;				/* calculated distance to back of UNCCollimator */
	 double						colCellUnitArea;					/*	Area of collimator unit for circular holes in a hexagonal 
																		UNCCollimator closed pack arrangement
																	*/
	double						colAccAngle;						/*	acc. angle of UNCCollimator, used to calculate range over
																	which we need to have detector loop for given photon angle
																	--there might be more than one detector position at which the
																	photon can be detected
																	(Note that this might not be the traditional "acceptance angle" 
																	for converging beam collimators)
																	*/
	double						colRangeOfDetAngles;				/*	UNCCollimator: stop_angle - start_angle (of detector) */
} colDataTy;



LOCALE	colDataTy	colData[PHG_MAX_PARAM_FILES];

LOCALE Boolean				ColDiscardBlue[PHG_MAX_PARAM_FILES];
LOCALE Boolean				ColDiscardPink[PHG_MAX_PARAM_FILES];


/* MACROS */

/*********************************************************************************
*
*			Name:		COLColIsCircular
*
*			Summary:	Returns true if the collimator cylinder is
*				circular, false otherwise.  If there is no collimator,
*				returns true if the target cylinder is circular, false
*				otherwise.
*				(Currently the return is always true, but this will
*				change when the elliptical object/target cylinders are
*				introduced.)
*			Arguments:
*				
*			Function return: boolean, true if the collimator is circular.
*
*********************************************************************************/
#define COLColIsCircular() true

/* PROTOTYPES */
Boolean			ColInitialize(Boolean doHistory);
Boolean			ColIsSlat(void);
Boolean			ColIsUNC(void);
void			ColPrintParams(void);
void			ColPrintReport(void);
double			ColGtOutsideRadius(void);
double			ColGtInsideRadius(void);
void			ColPETPhotons(PHG_Decay *decay,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					CollimatedPhotonsTy *colPhotonsPtr);
void			ColSPECTPhotons(PHG_Decay *decay,
					PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
					CollimatedPhotonsTy *colPhotonsPtr);
Boolean			ColGetMontCarloPETGeom(LbPfHkTy paramFlHk,
					LbUsFourByte numParams);
void			ColUpdateCollimatedPhotonBlock(PHG_TrackingPhoton	*trackingPhotonPtr);
void			ColTerminate(void);
LbUsFourByte	ColGtNumViews(void);
void			ColGtZLimits(double *zMin, double *zMax);

#undef LOCALE
#endif /* COL_HDR */
