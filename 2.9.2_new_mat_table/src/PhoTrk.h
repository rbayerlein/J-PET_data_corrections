/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992-2013 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			PhoTrk.h
*     Revision Number:		1.1
*     Date last revised:	20 June 2013
*     Programmer:			Steven Vannoy
*     Date Originated:		Tuesday September 29, 1992
*
*     Module Overview:	Definitions for photon tracking module.
*
*     References:       'Photon Tracking Processes' PHG design.  
*
**********************************************************************************
*
*     Global functions defined:
*
*     Global variables defined:   
*
**********************************************************************************
*
*		Revision Section (Also update version number, if relevant)
*
*		Programmer(s);		
*
*		Revision date:		
*
*		Revision description:
*
**********************************************************************************
*
*		Revision Section (Also update version number, if relevant)
*
*		Programmer(s);		Steven Gillispie
*
*		Revision date:		20 June 2013
*
*		Revision description:	Moved detGeomDecidePhotonAction from DetGeometric.c
*									to here as PhoTrkDecidePhotonAction.
*
*********************************************************************************/

#ifndef PHOTRACK
#define PHOTRACK

#include <stdio.h>
#include <math.h>

#include "LbTypes.h"
#include "LbMacros.h"

#include "Photon.h"

/* CONSTANTS */
/* TYPES */

/* Action returned from calculation of new position */
typedef enum { PhoTrkDetect, PhoTrkDiscard, PhoTrkInteract} PhoTrkActionTy;

/*	The following enumerate describes the types of actions to be
	taken as a result of projecting a photon to a new location. */
typedef enum { phtrkEnAc_Null,
				phtrkEnAc_Discard,
				phtrkEnAc_Interact,
				phtrkEnAc_Absorb,
				phtrkEnAc_ComptonScatter, 
				phtrkEnAc_CohScatter
} phtrkEn_ActionTy;

/* Define types and constants for cone beam forced detection */
#define CBFD_NUM_RADIAL_BINS		6
#define CBFD_NUM_AXIAL_BINS			10
#define CBFD_NUM_E_IN_BINS			10
#define CBFD_NUM_MU_IN_BINS			10
#define CBFD_NUM_DELTA_PHI_BINS		20
#define CBFD_NUM_DELTA_MU_OUT_BINS	8

typedef double PhoTrkCBFDDeltaMuTy[CBFD_NUM_RADIAL_BINS][CBFD_NUM_AXIAL_BINS];
typedef double PhoTrkCBFDProbTy[CBFD_NUM_E_IN_BINS][CBFD_NUM_MU_IN_BINS][CBFD_NUM_AXIAL_BINS][CBFD_NUM_RADIAL_BINS][CBFD_NUM_DELTA_PHI_BINS][CBFD_NUM_DELTA_MU_OUT_BINS];
typedef double PhoTrkCBFDTotProbTy[CBFD_NUM_E_IN_BINS][CBFD_NUM_MU_IN_BINS][CBFD_NUM_AXIAL_BINS][CBFD_NUM_RADIAL_BINS];

/* Forced detection types */
typedef struct {
	LbUsFourByte	num_ipo;			/* Number of ipo values */
	double			min_ipo;			/* Minimum ipo  value */
	double			delta_ipo;			/* Change per bin */

	LbUsFourByte	num_iwo;			/* Number of iwo values */
	double			min_iwo;			/* Minimum iwo  value */
	double			max_iwo;			/* Maximum iwo value */
	double			range_normalizer;	/* Used to raise negative values to positive bin indexes */
	double			delta_iwo;			/* Change per bin */

	LbUsFourByte	num_iwi;			/* Number of iwi values */
	double			min_iwi;			/* Minimum iwi  value */
	double			delta_iwi;			/* Change per bin */

	LbUsFourByte	num_iei;			/* Number of iei values */
	double			min_iei;			/* Minimum iei  value */
	double			delta_iei;			/* Change per bin */
} PhoTrkFDTInfoTy;	
	
typedef struct {
	double	*iwoTable;					/* Table of iwo values */
	double	*ipoTable;					/* Table of ipo values */
	double	*ipoCumTable;				/* Table of ipo cumulative values */
} PhoTrkFDTIwiTy;

typedef PhoTrkFDTIwiTy	**PhoTrkFDTIeiTy;	/* IEI table */

/* PROTOTYPES */
phtrkEn_ActionTy PhoTrkDecidePhotonAction(LbUsFourByte material, double photonEnergy, 
							Boolean modelingAbsorption, Boolean modelingCohScatter);
void			PhoTrkCalcAttenuation(PHG_TrackingPhoton *trackingPhotonPtr,
					double *attenuationPtr);
void			PhoTrkAttInitForcedDetection(PHG_TrackingPhoton *trackingPhotonPtr);
Boolean			phoTrkInitializeCBFDTbl(void);
void 			PhoTrkAttemptForcedDetection(PHG_TrackingPhoton *trackingPhotonPtr);
PhoTrkActionTy 	PhoTrkCalcNewPosition(PHG_TrackingPhoton *trackingPhotonPtr,
					PHG_Position startingPosition, PHG_Direction startingDirection,
					double photonEnergy, double freePath_Length, PHG_Position *newPosition,
					double *distTraveled, double *freePathsUsed);
Boolean			PhoTrkInitialize(void);
void			PhoTrkTerminate(void);
void			 PhoTrkUpdatePhotonPosition(PHG_TrackingPhoton *trackingPhotonPtr,
					PHG_Position *newPositionPtr);
void	PhoTrkCalcRange(double freePaths,
				PHG_Position startingPos,
				PHG_Direction direction, 
				PHG_Position *finalPosPtr,
				Boolean	*discard,
				LbFourByte	*sliceIndex,
				LbFourByte	*xIndex,
				LbFourByte	*yIndex);
void			PhoTrkProject(PHG_Position *startPos, PHG_Direction *direction, 
					double distance, PHG_Position *newPos);
void			PhoTrkEnterNewSlice(PHG_Position position, PHG_Direction direction,
						LbUsFourByte sliceIndex, double distanceTraveled,
						double *nextXPtr, double *nextYPtr, double *nextZPtr,
			   			double *generalXPtr, double *generalYPtr, double *generalZPtr,
						LbFourByte *slicePtr, LbFourByte *xIndexPtr,
						LbFourByte *yIndexPtr);


/* MACROS */

/*********************************************************************************
*
*			Name:		PHOTRKIsConeBeamForcedDetection
*
*			Summary:	Test to see if cone beam forced detection is being performed.
*
*			Arguments:
*
*			Function return: True or False.
*
*********************************************************************************/
#define PHOTRKIsConeBeamForcedDetection()	((PHG_IsForcedDetection() == true) && PHG_IsSPECT() && (PHG_IsCollimateOnTheFly() == true)) 

#endif /* PHOTRACK */
