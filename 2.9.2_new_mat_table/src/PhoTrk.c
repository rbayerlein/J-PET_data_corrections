/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1992-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		PhoTrk.c
*			Revision Number:	1.4
*			Date last revised:	23 July 2013
*
*			Programmer:			Steven Vannoy
*			Date Originated:	Tuesday, September 29, 1992
*
*			Module Overview:	Photon Tracking Processes.
*
*			References:			'Photon Tracking' PHG design.
*
**********************************************************************************
*
*			Global functions defined:
*				PhoTrkDecidePhotonAction
*				PhoTrkAttemptForcedDetection
*				PhoTrkAttInitForcedDetection
*				PhoTrkCalcNewPosition
*				PhoTrkInitialize
*				PhoTrkTerminate
*				PhoTrkUpdatePhotonPosition
*			Global variables defined:		none
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s);		
*
*			Revision date:		
*
*			Revision description:
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s);		Steven Gillispie
*
*			Revision date:		20 June 2013
*
*			Revision description:	Moved detGeomDecidePhotonAction from DetGeometric.c
*									to here as PhoTrkDecidePhotonAction.
*
*********************************************************************************/

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ProdTbl.h"
#include "SubObj.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "PhgMath.h"
#include "ColUsr.h"
#include "CylPos.h"
#include "UNCCollimator.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "PhoHStat.h"
#include "EmisList.h"
#include "PhoTrk.h"
#include "PhoHFile.h"
#include "phg.h"
#include "PhgBin.h"

/* LOCAL CONSTANTS */
#define MAX_INDEX	256		/* Maximum voxel index (temporary) */
#define PHOTRK_MAX_CELLS	(MAX_INDEX * 2 * 6)

/* LOCAL TYPES */
typedef struct {
	double		distance;		/* Distance traveled through cell */
	double		attenuation;	/* Cell's attenuation */
	double		freePaths;		/* Cell's free path utilization */
	LbFourByte	xIndex;			/* X index of voxel represented by this cell */
	LbFourByte	yIndex;			/* Y index of voxel represented by this cell */
	LbFourByte	sliceIndex;		/* Z index of voxel represented by this cell */
} phoTrkCellInfoTy;

/* LOCAL GLOBALS */
Boolean						phoTrkIsInitialized = false;			/* Initialization flag */
char						phoTrkErrString[1024];					/* Handy error string */
phoTrkCellInfoTy			phoTrkCellInfoTable[PHOTRK_MAX_CELLS];	/* Should be dynamically allocated */
LbUsTwoByte					phoTrkCellsInUse;						/* Number of cells used */
LbUsTwoByte					phoTrkMaxCellCount = PHOTRK_MAX_CELLS;	/* Maximum cell count */
PhoTrkFDTInfoTy				phoTrkFDTInfo;							/* Forced detection information */
PhoTrkFDTIeiTy				phoTrkFDTable = 0;							/* Our phoTrkFDTable */
double						phoTrkTotalKN[1000];					/* Table = Integral of KN by energy */

PhoTrkCBFDDeltaMuTy			*phoTrkCBFDMaxDeltaMu;
PhoTrkCBFDDeltaMuTy			*phoTrkCBFDMinDeltaMu;
PhoTrkCBFDProbTy			*phoTrkCBFDProbDeltaPhiMu;
PhoTrkCBFDProbTy			*phoTrkCBFDCumProbDeltaPhiMu;
PhoTrkCBFDTotProbTy			*phoTrkCBFDTotalProbAccept;

PhoTrkFDTIeiTy				phoTrkFDTable2;							/* Our phoTrkFDTable */
double						phoTrkRadiusOfFocalCircle;				/* Used for cone beam FD */
double						phoTrkFocalLength;						/* Used for cone beam FD */	
double						phoTrkCollimatorZMin;					/* Used for cone beam FD */
double						phoTrkCollimatorZMax;					/* Used for cone beam FD */


/* LOCAL MACROS */
/*********************************************************************************
*
*			Name:		PHOTRKGetMinDeltaT
*
*			Summary:	Compare three delta T values, and return minimum.
*
*			Arguments:
*				dtx		- Delta t with respect to x axis.
*				dty		- Delta t with respect to y axis.
*				dtz		- Delta t with respect to z axis.
*
*			Function return: Minimum of 3 values.
*
*********************************************************************************/
#define	PHOTRKGetMinDeltaT(dtx, dty, dtz) ((((dtx) <= (dty)) && ((dtx) <= (dtz))) ? (dtx) : \
 (((dty) <= (dtx)) && ((dty) <= (dtz))) ? (dty) : (dtz))

/* LOCAL PROTOTYPES */
void			phoTrkCalcFreePaths(PHG_TrackingPhoton *trackingPhotonPtr,
					PHG_Position *startingPosPtr,
					double distanceToTravel,
					PHG_Direction direction,
					double photonEnergy,
					double *freePathsPtr);
Boolean			phoTrkCalcScatterAngleSPECT(PHG_TrackingPhoton	*trackingPhotonPtr);
Boolean			phoTrkCalcScatterAnglePET(PHG_TrackingPhoton	*trackingPhotonPtr);
Boolean				phoTrkCalcAcceptanceRange(PHG_TrackingPhoton	*trackingPhotonPtr,
					double *minSinePtr, double *maxSinePtr);
void 			phoTrkCalcAzimuthToDir(double inclination, double azimuth,
					PHG_Direction *dirPtr);
void			phoTrkCalcDirToAzimuth(PHG_Direction *dirPtr,
						double *inclinationPtr, double *azimuthPtr);
LbUsFourByte	phoTrkLookupProb(double *phiArray, LbUsFourByte numElements, double phi);
void			phoTrkCalcCritZoneFreePaths(PHG_TrackingPhoton	*trackingPhotonPtr,
						PHG_Intersection *intersectionPtr,
						double *freePathsToEnterPtr,
						double *freePathsToExitPtr);
void			phoTrkDoWeightWindow(PHG_TrackingPhoton *trackingPhotonPtr,
						LbTwoByte *numSplitsPtr, LbTwoByte *numRoulettesPtr, Boolean *discardItPtr);
Boolean			phoTrkInitializeFDTbl(PhoTrkFDTInfoTy	*fdInfoPtr);
double			phoTrkCalcCBFDOmega(double deltaRphi, double phoTrkRadiusOfFocalCircle, 
					double radialPos, double zPos);		
Boolean			phoTrkPositionIsAcceptable(PHG_TrackingPhoton	*photonPtr,
					double *minSinePtr, double *maxSinePtr);

/* Functions */
/*********************************************************************************
*
*			Name:		PhoTrkCalcAttenuation
*
*			Summary:	Calculate the attenuation a photon would undergo while
*						traveling a specified distance.
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*				double				distanceToTravel	- How far to track.
*				double				*attenuationPtr		- The attenuation
*
*			Function return: None.
*
*********************************************************************************/
void PhoTrkCalcAttenuation(PHG_TrackingPhoton *trackingPhotonPtr,
		double *attenuationPtr)
{
	double			distToObject;		/* Distance to object cylinder */
	double			distToTarget;		/* Distance to target cylinder */
	double			freePaths;			/* Free paths necessary to exit */
	PHG_Position	posOnTarCylinder;	/* Projected position on target cylinder */
	LbFourByte		cellIndex;			/* Index for looping through cells */
	
	do { /* Process Loop */
	
		/* Clear cells in use variable */
		phoTrkCellsInUse = 0;
		
		/* Copy current location */
		posOnTarCylinder = trackingPhotonPtr->location;
		
		/* See if we don't project to target cylinder within bounds */
		if (!CylPosProjectToTargetCylinder(&posOnTarCylinder,
				&(trackingPhotonPtr->angle), &distToTarget)) {
			PhgAbort("Photon does not reach target cylinder within bounds,"
				" something is wrong;", false);
		}
		
		/* Compute the distance to the object cylinder */
		CylPosCalcDistanceToObjectSurface(&(trackingPhotonPtr->location),
			&(trackingPhotonPtr->angle), &distToObject);
			
		/*	Compute attenuation (free paths) through object
			NOTE: phoTrkCalcFreePaths routine sets the local global 
			"phoTrkCellsInUse" for later use.
		*/
		phoTrkCalcFreePaths(trackingPhotonPtr, &(trackingPhotonPtr->location),
			distToObject, trackingPhotonPtr->angle,
			trackingPhotonPtr->energy, &freePaths);

		/* Loop through and accumulate the attenuation */
		*attenuationPtr = 0.0;
		for (cellIndex = 0; cellIndex < phoTrkCellsInUse; cellIndex++){
		
			*attenuationPtr += (phoTrkCellInfoTable[cellIndex].attenuation
				* phoTrkCellInfoTable[cellIndex].distance);
		}
		
	} while (false);
	
	/* Clear cells in use variable */
	phoTrkCellsInUse = 0;
}


/*********************************************************************************
*
*		Name:			PhoTrkDecidePhotonAction
*
*		Summary:		Decide, probabilistically, what action the photon will undergo.
*							
*		Arguments:
*			LbUsFourByte	material			- The interaction material.
*			double			photonEnergy		- The energy of the photon.
*			Boolean			modelingAbsorption	- If absorption is simulated.
*			Boolean			modelingCohScatter	- If coherent scatter is simulated.
*
*		Function return: A phtrkEn_ActionTy type of photon action.
*
*********************************************************************************/

phtrkEn_ActionTy PhoTrkDecidePhotonAction(LbUsFourByte material, double photonEnergy, 
							Boolean modelingAbsorption, Boolean modelingCohScatter)

{
	phtrkEn_ActionTy	action = phtrkEnAc_Interact;	/* Returned result */
	double				scatterProbability;				/* Probability of a scatter */
	double				comptonToScatterProbability;	/* Ratio of prob of Compton to 
															prob of scatter */
	double				interactionProbability;			/* Interaction-type probability */
	double				comptonScatterProb;				/* Probability of a Compton scatter */
	
	
	/* Look up the fixed probabilities that define the type of interaction; 
		note that these vary whether coherent scatter is being modeled or not */
	
	/* Get probability of any scatter */
	scatterProbability = 
		SubObjGetProbScatter(material, photonEnergy, modelingCohScatter);
	
	/* Get conditional probability of a Compton scatter */
	comptonToScatterProbability = 
		SubObjGetProbComptonCondnl(material, photonEnergy, modelingCohScatter);
	
	
	/* Get the random probability that determines the type of interaction */
	interactionProbability = PhgMathGetRandomNumber();
	
	
	/* Decide the photon's fate */
	if (modelingAbsorption) {
		#ifdef PHG_DEBUG
			if (PHGDEBUG_DetAbsorbOnly())
				interactionProbability = 1.0;
		#endif
		
		if (interactionProbability > scatterProbability) {
			/* The photon will be absorbed:  
				change the action status from interaction to absorption */
			action = phtrkEnAc_Absorb;
		}
		else {
			/* The photon may be scattered; compute the Compton scatter probability */
			comptonScatterProb = scatterProbability * comptonToScatterProbability;
		}
	}
	else {
		/* The photon may be scattered; compute the Compton scatter probability */
		/* (Note:  non-absorption changes the Compton-coherent scatter probability 
			from a conditional probability to an absolute probability) */
		comptonScatterProb = comptonToScatterProbability;
	}
	
	/* If the photon wasn't absorbed, check for a scatter */
	if (action == phtrkEnAc_Interact) {
		if (interactionProbability > comptonScatterProb) {
			/* Coherent scatter */
			action = phtrkEnAc_CohScatter;
		}
		else {
			/* A Compton scatter occurs */
			action = phtrkEnAc_ComptonScatter;
		}
	}
	
	return(action);
}


/*********************************************************************************
*
*			Name:		PhoTrkAttInitForcedDetection
*
*			Summary:	Attempt to extend the photon out through the acceptable
*						target surface area.
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*
*			Function return: None.
*
*********************************************************************************/
void PhoTrkAttInitForcedDetection(PHG_TrackingPhoton *trackingPhotonPtr)
{
	Boolean			photonDetected;		/* Did we detect the photon */
	double			distToObject;		/* Distance to object cylinder */
	double			distToTarget;		/* Distance to target cylinder */
	double			freePaths;			/* Free paths necessary to exit */
	PHG_Position	posOnTarCylinder;	/* Projected position on target cylinder */
	
	do { /* Process Loop */
	
		/* Clear cells in use variable */
		phoTrkCellsInUse = 0;

		/* Clear our detection flag */
		photonDetected = false;
		
		/* Increment statistics for forced detection attempts */
		PhoHStatUpdateForDetAttemptWeight(trackingPhotonPtr->decay_weight *
			trackingPhotonPtr->photon_primary_weight);
		
		/* See if we are outside acceptance angle */
		if (fabs(trackingPhotonPtr->angle.cosine_z) > PHGGetSineOfAccAngle()) {
			break;
		}
		
		/* Copy current location */
		posOnTarCylinder = trackingPhotonPtr->location;
		
		/* See if we don't project to target cylinder within bounds */
		if (!CylPosProjectToTargetCylinder(&posOnTarCylinder,
				&(trackingPhotonPtr->angle), &distToTarget)) {
			break;
		}
		
		/* Compute the distance to the object cylinder */
		CylPosCalcDistanceToObjectSurface(&(trackingPhotonPtr->location),
			&(trackingPhotonPtr->angle), &distToObject);
			
		/*	Compute attenuation (free paths) through object
			NOTE: phoTrkCalcFreePaths routine sets the local global 
			"phoTrkCellsInUse" for later use.
		*/
		phoTrkCalcFreePaths(trackingPhotonPtr, &(trackingPhotonPtr->location),
			distToObject, trackingPhotonPtr->angle,
			trackingPhotonPtr->energy, &freePaths);
		
		/* Adjust the photon's weight */
		trackingPhotonPtr->photon_primary_weight *=  exp(-freePaths);
		
		
		/* Set our success flag */
		photonDetected = true;
	} while (false);
	
	/* Handle result */
	if (photonDetected == true) {
		/* Update forced detection statistics */
		PhoHStatUpdateForDetHitWeight(trackingPhotonPtr->decay_weight *
			trackingPhotonPtr->photon_primary_weight);

		/* Update for the distance to the target cylinder */
		trackingPhotonPtr->travel_distance = distToTarget;
			
		/* Set our position on the cylinder */
		trackingPhotonPtr->location = posOnTarCylinder;
		
		/* Record the detection to emis list */
		EmisListDoDetection(trackingPhotonPtr);

		/* If this is not a track as scatter photon we 
			need to clear the cellsInUse flag, because this photon
			is going to be discarded and we don't want the list
			hanging around for the next photon (which might
			not be a track-as-primary photon, hence it would
			miss this routine.
		*/
		if (PHG_IsTrackAsScatter(trackingPhotonPtr) == false) {
			/* Clear cells in use variable */
			phoTrkCellsInUse = 0;
			
		}

	}		
}


/*********************************************************************************
 *
 *			Name:		PhoTrkAttemptForcedDetection
 *
 *			Summary:	Attempt to extend the photon out through the acceptable target
 *			  			surface area.
 *			Arguments:
 *				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
 *
 *			Function return: None.
 *
 *********************************************************************************/
void PhoTrkAttemptForcedDetection(PHG_TrackingPhoton *trackingPhotonPtr)
{
    Boolean				photonDetected;			/* Did the photon get detected */
    Boolean				discardIt;				/* Should photon be discarded (Russian Roulette) */
	double				freePathsToEnter;		/* Free Paths to enter */
    double				freePathsToExit;		/* Free paths to exit */
    double				freePathsToTravel;		/* Free paths to travel until interaction */
    double				freePathsToInteraction;	/* Free paths to travel until interaction */
    double				freePathsToDetection;	/* Free paths to get to detection */
    double				distanceTraveled;		/* Distance traveled through critical zone */
	double				fdWeight;				/* Weight of photon for forced detection */
	double				prob1;					/* Temp probability var */
	double				prob2;					/* Temp probability var 2 */
	double				randFromExp;			/* Random number from exponential districution. */
	double				fpToMoveCoeff;			/* Computed value used to determine number of free paths to move through the critical zone. */
    PHG_Position		interactionPosition;	/* Location of interaction */
    PHG_Position		detectionPosition;		/* Location of detection */
    PHG_Intersection	intersection;			/* The entry location of the critical zone */
    PhoTrkActionTy		action;					/* Action to be taken based on new position */
	LbTwoByte			numSplits;				/* Number of times to split this photon */
	LbTwoByte			numRoulettes;			/* Number of times to roulette */
	LbTwoByte			numFDAttempts;			/* Number of times to perform forced detection */
    LbTwoByte			attemptNum;   			/* Current FD attempt */
	LbUsFourByte		wholeNum;				/* Whole number temp variable */
	LbUsFourByte		cellIndex;				/* Current cell for tracking */
	LbUsFourByte		cellsInUseCopy;			/* Storage for a copy of this local global */
	PHG_TrackingPhoton	incomingPhoton;			/* Storage for restoring photon */
						
	do { /* Process Loop */
		
		/* Clear detection flag */
		photonDetected = false;

		/* Get critical zone intersection, IF there is not one then break out of the
			process loop causing the FD attempt to fail
		*/
		if (CylPosWillIntersectCritZone(trackingPhotonPtr->location,
				trackingPhotonPtr->angle, &intersection) == false) {
			
			break;
		}
		
		/* Set number of FD attempts to 1. If splitting was done, it will be adjusted below */
		numFDAttempts = 1;
		
		/* Save a copy of the incoming photon */
		incomingPhoton = *trackingPhotonPtr;

		/* Set FD weight to photon's weight. It may be adjusted as above */
		fdWeight = trackingPhotonPtr->photon_scatter_weight;
		
		/* Do weight windowing if we have a computed table and a valid ratio */
		if (PHG_IsProductivityComputed() && (PhgRunTimeParams.PhgMinWWRatio != 1.0)) {
			
			/* Now do weight windowing */
			phoTrkDoWeightWindow(trackingPhotonPtr,
				&numSplits, &numRoulettes, &discardIt);
			
			/* See if Russian Roulette was attempted */
			if (numRoulettes != -1) {
				
				/* Update roulette statistics */
				PhoHStatUpdateRoulettedPhotons(numRoulettes);
				
				/* See if we should discard the photon */
				if (discardIt == true)
				break;
			}

			/* See if Splitting was done */
			if (numSplits > 0) {
				LbInPrintf("\nnumSplits == %d", numSplits);
				if (numSplits >= (PHG_MAX_DETECTED_PHOTONS - 10))
					numSplits = PHG_MAX_DETECTED_PHOTONS - 10;
					
				/* Update statistics */
				PhoHStatUpdateSplitPhotons(numSplits);
				
				/* Adjust number of FD attempts (note that numSplits could be zero, its ok) */
				numFDAttempts = numFDAttempts + numSplits;
				
   				/* Compute the split weight */
				fdWeight = trackingPhotonPtr->photon_scatter_weight/numFDAttempts;
			}
		}
		
		
		/* Calculate free paths to enter critical zone and free paths to 
		   exit the critical zone.
		   
		   NOTE that it is possible for the freePathsToEnter to equal freePathsToExit,
		   in this case we bail out of the FD attempt.
		*/
		{
			/* Get CZ free paths */
			phoTrkCalcCritZoneFreePaths(trackingPhotonPtr,
				&intersection, &freePathsToEnter, &freePathsToExit);
			
			/* Check for realistic conditions  */
			if (PhgMathRealNumAreEqual(freePathsToEnter, freePathsToExit, 
		        -7, 0, 0, 0) == true)
				break;
		}
		
		/*	Adjust photon weight by probability that it will interact 
			in the critical zone
		*/
		{
			/* Compute first temporary probability */
			prob1 = exp(-freePathsToEnter);
			
			/* Compute second temprorary probability */
			prob2 = exp(-freePathsToExit);
			
			/* Adjust the weight */
			fdWeight = fdWeight * (prob1 - prob2);
		}

		/* Perform Forced detection */
		for (attemptNum = 1; attemptNum <= numFDAttempts; attemptNum++) {
			
			/* Clear photon detected flag each time */
			photonDetected = false;
			
			/* "Reset" photon weight */
			trackingPhotonPtr->photon_scatter_weight = fdWeight;
			
			/* Increase the photon's number of scatters by one */
			trackingPhotonPtr->num_of_scatters++;
			
			/* Increment statistics for forced detection attempts */
			PhoHStatUpdateForDetAttemptWeight(trackingPhotonPtr->decay_weight *
				trackingPhotonPtr->photon_scatter_weight);
			
			do { /* Process Loop */
				
			    /*	Pick a random location for interaction within the critical 
					zone.
				*/
				
				{
					/* We will compute a random number from a truncated 
					   exponential distribution with a range of 
					   freePathsToEnter -- freePathsToExit
					*/

					/* Start with a random from the exponential distribution */
					PhgMathGetTotalFreePaths(&randFromExp);
					
					/* Truncate to desired range */
					{
						/* Compute whole number part NOTE I AM IGNORING THE
						   POSSIBILITY OF INTEGER OVERFLOW
						*/
						wholeNum = (LbUsFourByte) (randFromExp/
							(freePathsToExit - freePathsToEnter));

						fpToMoveCoeff = ((randFromExp/
					   		(freePathsToExit - freePathsToEnter)) -
							wholeNum)
							* (freePathsToExit - freePathsToEnter);
						
					}

					/* Now determine how far to travel */
			    	freePathsToTravel = freePathsToEnter +
						fpToMoveCoeff;
				
					/* Clear looping variables */
					distanceTraveled  = 0;
					freePathsToInteraction = 0;
				
					/* Loop through cell table until freePathsToTravel have been used up */
					for (cellIndex = 0; cellIndex < phoTrkCellsInUse; cellIndex++) {
					
						/* If the free paths travelled so far + the free paths in this cell
							are greater than the selected free paths to interact, force the
							interaction in this cell
						*/ 
						 if ((freePathsToInteraction + phoTrkCellInfoTable[cellIndex].freePaths)
								>= freePathsToTravel) {
								
							/* Only use necessary amount of cell */
							distanceTraveled += (freePathsToTravel - freePathsToInteraction)/
								phoTrkCellInfoTable[cellIndex].attenuation;
								
							/* Set free paths used */
							freePathsToInteraction = freePathsToTravel;
				
							/* We are done */
							break;
								
						}
						else {
						
							/* Accumulate the distance traveled plus the free paths of this cell */
							distanceTraveled += phoTrkCellInfoTable[cellIndex].distance;
							freePathsToInteraction += phoTrkCellInfoTable[cellIndex].freePaths;
						}
					}
					
					/* Verify we did not fall out of the above loop without accumulating the
						freePathsToInteraction 
					*/
					#ifdef PHG_DEBUG
						if (cellIndex == phoTrkCellsInUse) {
							/* We should NOT have run out of cells before executing the if
								clause above
							*/
							PhgAbort("Error in FD calculation of scatter location (PhoTrkAttemptForcedDetection)",
								true);
						}
					#endif
					
					/* Update x/y/z index for interaction position */
					trackingPhotonPtr->xIndex = phoTrkCellInfoTable[cellIndex].xIndex;
					trackingPhotonPtr->yIndex = phoTrkCellInfoTable[cellIndex].yIndex;
					trackingPhotonPtr->sliceIndex = phoTrkCellInfoTable[cellIndex].sliceIndex;
					
					/* Update position */
					PhoTrkProject(&trackingPhotonPtr->location, &trackingPhotonPtr->angle,
						distanceTraveled, &interactionPosition);
						
				}
				
			    /* Save the new location */			
			    trackingPhotonPtr->location = 
					interactionPosition;
				
				/* Update the travel distance */
				trackingPhotonPtr->travel_distance += distanceTraveled;
/* 9/1/98 !!SV Changed from SubObjGetProbComptonScatter to just ProbScatter */
                /* Adjust the weight for probability of compton scatter */
                trackingPhotonPtr->photon_scatter_weight *=
                    SubObjGetProbScatterInObj(trackingPhotonPtr);
				
			    /* Try to update photon with scatter angle that assures 
				   detection.
				   Update weight accordingly.
				   If no detectable angle exits, break.
				*/
			    if (PHOTRKIsConeBeamForcedDetection()) {
			    	if (phoTrkCalcScatterAngleSPECT(trackingPhotonPtr) == false) {
			        	break;
			    	}
			    }
			    else {
			    	if (phoTrkCalcScatterAnglePET(trackingPhotonPtr) == false) {
			        	break;
			    	}
			    }
			    
					
			    /* Now track to detection */
			    {
					/* Select free path to get us out of
					   object cylinder without interactions
					*/
					freePathsToTravel = -1;
					
					/* Save the cells in use flag value, and set to zero.
						It will be restored after tracking the photon out of the cylinder.
					*/
					cellsInUseCopy = phoTrkCellsInUse;
					phoTrkCellsInUse = 0;
					/* Compute new location, result should be 
					   detection (but doesn't have to be).
					*/
					action = PhoTrkCalcNewPosition(
						trackingPhotonPtr,
						trackingPhotonPtr->location,
						trackingPhotonPtr->angle,
						trackingPhotonPtr->energy,
						freePathsToTravel,
						&detectionPosition,
						&distanceTraveled,
						&freePathsToDetection);
						
					/* Restore cells in use from copy */
					phoTrkCellsInUse = cellsInUseCopy;
			    }
				
			    /* If we missed the target, get out of here */
			    if (action != PhoTrkDetect) {
					break;
			    }

			    /* Adjust the weight based on probability
			       that interaction did not occur on the way 
			       out
				*/
				trackingPhotonPtr->photon_scatter_weight *= 
					exp(-freePathsToDetection);

				
				/* Update the travel distance to the object surface */
				trackingPhotonPtr->travel_distance += distanceTraveled;
				
				/* Update the location */
				trackingPhotonPtr->location = detectionPosition;
				
   			    /* If we made it here, we got detected */
			    photonDetected = true;
			} while (false);
			
			/* Record outcome */
			if (photonDetected == true) {
			
				/* Record detection */
				EmisListDoDetection(trackingPhotonPtr);
				
				/* Update forced detection statistics */
				{
					PhoHStatUpdateForDetHitWeight(
					   	trackingPhotonPtr->decay_weight *
						trackingPhotonPtr->photon_scatter_weight);
				}
			}

			/* Restore the incoming photon */
    	    *trackingPhotonPtr = incomingPhoton;

		} /* End of FD for loop */

	} while (false);
	/* Restore the incoming photon */
	*trackingPhotonPtr = incomingPhoton;

}

/*********************************************************************************
*
*			Name:		phoTrkCalcAcceptanceRange
*
*			Summary:		Calculate sine of minimum and maximum acceptance
*						angle for a photon to be detected from its current
*						position.
*
*
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr		- The tracking photon.
*				double				*minSinePtr				- Sine of minimum angle.
*				double				*maxSinePtr				- Size of maximum angle.
*
*			Function return: True if valid angle computed.
*
*********************************************************************************/
Boolean phoTrkCalcAcceptanceRange(PHG_TrackingPhoton	*trackingPhotonPtr,
			double *minSinePtr, double *maxSinePtr)
{
	Boolean	angleFound;			/* Did we find an acceptable angle? */
	double	rTar;				/* Radius of target cylinder */
	double	zMin;				/* Minimum z of target cylinder */
	double	zMax;				/* Maximum z of target cylinder */
	double	accSine;			/* Acceptance angle sine */
	double	rScat;				/* Distance of scatter point from center of cylinder */
	double	hAbove;				/* Height of scatter point above target cylinder (may be negative) */
	double	hBelow;				/* Height of scatter point below target cylinder (may be negative) */
	double	angleTopNearSine;	/* Sine of angle from scatter point to the nearest point at the top of the cylinder */
	double	angleTopFarSine;	/* Sine of angle from scatter point to the farthest point at the top of the cylinder */
	double	angleBotNearSine;	/* Sine of angle from scatter point to the nearest point at the bottom of the cylinder */
	double	angleBotFarSine;	/* Sine of angle from scatter point to the farthes point at the bottom of the cylinder */
	double	tempMax;			/* Temporary value for computing limits */
	double	tempMin;			/* Temporary value for computing limits */
	
	/* Assume no angle found */
	angleFound = false;
	
	do { /* Process Loop */
		/* Initialize computation variables */
		{
			/* Get limits on z range */
			CylPosGetLimitZRange(&zMin, &zMax);
	
			/* Get radius of target cylinder */
			rTar = CylPosGetTargetRadius();
			
			/* Get acceptance angle sine */
			accSine = PHGGetSineOfAccAngle();
		}
		
		/* Compute distance of scatter point to center of cylinder */
		rScat = PHGMATH_SquareRoot(PHGMATH_Square(trackingPhotonPtr->location.x_position) +
			PHGMATH_Square(trackingPhotonPtr->location.y_position));
			
		/* Compute height of scatter point above target cylinder */
		hAbove = trackingPhotonPtr->location.z_position - zMax;
		
		/* Compute height of scatter point below target cylinder */
		hBelow = trackingPhotonPtr->location.z_position - zMin;
		
		/* Compute sine of angle from scatter point to the nearest point at the top of the cylinder */
		angleTopNearSine = -(hAbove/
			PHGMATH_SquareRoot(PHGMATH_Square(rTar - rScat) + PHGMATH_Square(hAbove)));
		
		/* Compute sine of angle from scatter point to the farthest point at the top of the cylinder */
		angleTopFarSine = -(hAbove/
			PHGMATH_SquareRoot(PHGMATH_Square(rTar + rScat) + PHGMATH_Square(hAbove)));
		
		/* Compute sine of angle from scatter point to the nearst point at the bottom of the cylinder */
		angleBotNearSine = -(hBelow/
			PHGMATH_SquareRoot(PHGMATH_Square(rTar - rScat) + PHGMATH_Square(hBelow)));
		
		/* Compute sine of angle from scatter point to the farthest point at the bottom of the cylinder */
		angleBotFarSine = -(hBelow/
			PHGMATH_SquareRoot(PHGMATH_Square(rTar + rScat) + PHGMATH_Square(hBelow)));
	
	
		/* Select min angle sine limit */
		tempMin = ((angleBotNearSine < angleBotFarSine) ? angleBotNearSine: angleBotFarSine);
		
		/* Min acceptable gets max of target cylinder min and current location min */
		*minSinePtr = PHGMATH_Max(-accSine, tempMin);
		
		/* Select max angle sine limit */
		tempMax = PHGMATH_Max(angleTopNearSine, angleTopFarSine);
		
		/* Mfax acceptable gets min of target cylinder min and current location min */
		*maxSinePtr = ((tempMax < accSine) ? tempMax : accSine);
		
		/* Make sure we succeeded */
		if (*minSinePtr > *maxSinePtr)
			break;
			
		#ifdef PHG_DEBUG
			if ((*minSinePtr != -accSine) || (*maxSinePtr != accSine))
				angleFound = angleFound;
		#endif
		
		/* If we made it here, we got an acceptable angle */
		angleFound = true;
	} while (false);
	
	return (angleFound);
}
/*********************************************************************************
*
*			Name:		phoTrkPositionIsAcceptable	
*
*			Summary:	Determine if current photon position is acceptable
*						for performing forced detection.
*
*
*			Arguments:
*				PHG_TrackingPhoton	*photonPtr				- The tracking photon.
*				double				*minSinePtr				- Sine of minimum angle.
*				double				*maxSinePtr				- Size of maximum angle.
*
*			Function return: True if valid angle computed.
*
*********************************************************************************/
Boolean phoTrkPositionIsAcceptable(PHG_TrackingPhoton	*photonPtr,
			double *minSinePtr, double *maxSinePtr)
{
	Boolean posIsAcceptable;
	double r, zMin, zMax, safetyMargin;
	
	/* If not cone beam then use standard routine */
	if (PHOTRKIsConeBeamForcedDetection() == false) {
		return(phoTrkCalcAcceptanceRange(photonPtr, minSinePtr, maxSinePtr));
	}
	
	do {
	
		r = PHGMATH_RadialPos(photonPtr->location.x_position, photonPtr->location.y_position);
		safetyMargin = phoTrkFocalLength * ( (2.0*ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleRadius) /
				PHGMATH_SquareRoot( PHGMATH_Square(2.0*ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleRadius) +
				PHGMATH_Square(ColRunTimeParams[ColCurParams].UNCSPECTCol.Thickness) ) );	 /* This is approx focalLength * sin(thetaMaxDev) */
		zMin = phoTrkCollimatorZMin - safetyMargin;		
		zMax = phoTrkCollimatorZMax + safetyMargin;		
		
		if (((phoTrkFocalLength * (photonPtr->location.z_position/(r+phoTrkRadiusOfFocalCircle))) < zMax) &&
				((phoTrkFocalLength * (photonPtr->location.z_position/(r+phoTrkRadiusOfFocalCircle))) > zMin)) {	
			posIsAcceptable = true;
		}
		else {
			posIsAcceptable = false;
		}
			
	} while (false);
	
	return (posIsAcceptable);
}
/*********************************************************************************
*
*			Name:		phoTrkCalcDirToAzimuth
*
*			Summary:	Convert direction cosine to inclination and azimuth.
*
*
*			Arguments:
*				PHG_Direction	*dirPtr			- The incoming direction pointer.
*				double			*inclinationPtr	- The computed inclination.
*				double			*azimuthPtr		- The computed azimuth.
*			Function return: None.
*
*********************************************************************************/
void phoTrkCalcDirToAzimuth(PHG_Direction *dirPtr,
			double *inclinationPtr, double *azimuthPtr)
{
	/* Set inclination to z cosine */
	*inclinationPtr = dirPtr->cosine_z;
	
	/* Compute azimuth */
	*azimuthPtr = atan2(dirPtr->cosine_y, dirPtr->cosine_x);
}

/*********************************************************************************
*
*			Name:		phoTrkCalcAzimuthToDir
*
*			Summary:	Convert inclination and azimuth to direction cosine.
*
*
*			Arguments:
*				double			inclination		- The incoming inclination.
*				double			azimuth			- The  incoming azimuth.
*				PHG_Direction	*dirPtr			- The computed direction pointer.
*			Function return: None.
*
*********************************************************************************/
void phoTrkCalcAzimuthToDir(double inclination, double azimuth,
		PHG_Direction *dirPtr)
{
	double	zSineAbs;	/* Absolute value of sine */
	
	/* Set z cosine to inclination */
	dirPtr->cosine_z = inclination;
	
	/* Compute x/y cosines */
	zSineAbs = PHGMATH_SquareRoot(1 - PHGMATH_Square(dirPtr->cosine_z));
	
	dirPtr->cosine_x = PHGMATH_Cosine(azimuth) * zSineAbs;
	
	dirPtr->cosine_y = PHGMATH_Sine(azimuth) * zSineAbs;

}

/*********************************************************************************
*
*			Name:		phoTrkCalcFreePaths
*
*			Summary:	Determine the cumulative attenuation the photon will
*						encounter traveling a given distance through the object.
*
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The photon.
*				PHG_Position		*startingPosPtr		- The starting position.
*				double				distanceToTravel	- The distance to travel.
*				PHG_Direction		direction			- The direction of travel.
*				double				photonEnergy		- The photon's energy.
*				double				*freePathsPtr		- The accumulated free paths.
*
*			Function return: None.
*
*********************************************************************************/
void	phoTrkCalcFreePaths(PHG_TrackingPhoton	*trackingPhotonPtr,
				PHG_Position *startingPosPtr,
				double distanceToTravel,
				PHG_Direction direction, 
				double photonEnergy, double *freePathsPtr)
{
	double			cellAttenuation;		/* Attenuation of current cell */
	double			distanceTracked;		/* Total distance tracked */
	double			initialZDistance;		/* Initial distance from current z to z boundary */
	double			initialXDistance;		/* Initial distance from current x to x boundary */
	double			initialYDistance;		/* Initial distance from current y to y boundary */
	double			distToNextX;			/* Distance to next x crossing */
	double			distToNextY;			/* Distance to next y crossing */
	double			distToNextZ;			/* Distance to next z crossing */
	double			generalDistToX;			/* General distance to move relative to x axis */
	double			generalDistToY;			/* General distance to move relative to y axis */
	double			generalDistToZ;			/* General distance to move relative to z axis */
	double			nextDist;				/* Next incremental move */
	LbFourByte		newSliceIndex;			/* New slice index for entering new slice */
	PHG_Position	newPosition;			/* Used to update slice info */


	/* Clear "distance traveled" counter */
	distanceTracked = 0.0;
	
	/* Clear free paths counter */
	*freePathsPtr = 0.0;

	/* Clear distance marker */
	nextDist = 0.0;
	
	/* Get initial distance to voxel wall */
	SubObjGetInnerCellDistance(startingPosPtr, &direction, trackingPhotonPtr->sliceIndex,
		 trackingPhotonPtr->xIndex,  trackingPhotonPtr->yIndex,
		&initialXDistance, &initialYDistance, &initialZDistance);
	
	/* Avoid div by zero from direction cosines */
	{
		if ((direction.cosine_x <= 0.0000001) && 
				(direction.cosine_x >= -0.0000001)) {				
			direction.cosine_x = 
				(0.0000001 * ((direction.cosine_x < 0) ? -1 : 1));
		}
		if ((direction.cosine_y <= 0.0000001) && 
				(direction.cosine_y >= -0.0000001)) {				
			direction.cosine_y = 
				(0.0000001 * ((direction.cosine_y < 0) ? -1 : 1));
		}
		if ((direction.cosine_z <= 0.0000001) && 
				(direction.cosine_z >= -0.0000001)) {				
			direction.cosine_z = 
				(0.0000001 * ((direction.cosine_z < 0) ? -1 : 1));
		}
	}

	/* Calculate initial distances  (normalized to positive value) */
	{
		distToNextX = initialXDistance/direction.cosine_x;
		distToNextY = initialYDistance/direction.cosine_y;
		distToNextZ = initialZDistance/direction.cosine_z;
	}
	
	/* Calculate general distances */
	{
		generalDistToX = SUBOBJGetSliceAttVoxelWidth(trackingPhotonPtr->sliceIndex)/fabs(direction.cosine_x);
		generalDistToY = SUBOBJGetSliceAttVoxelHeight(trackingPhotonPtr->sliceIndex)/fabs(direction.cosine_y);
		generalDistToZ = SUBOBJGetSliceVoxelDepth(trackingPhotonPtr->sliceIndex)/fabs(direction.cosine_z);
	}
	
	do { /* Tracking Loop */

		/* Get shortest distance to next cell boundary */
		nextDist = PHOTRKGetMinDeltaT(distToNextX, distToNextY, distToNextZ);
		
		/* If the distance to the next voxel boundary is greater than our target distance
			(the entry/exit point of the critical zone), truncate the distance we
			travel through the voxel to the desired boundary.
		*/
		if (nextDist >= distanceToTravel)
			nextDist = distanceToTravel;
			
		/* Get attenutation of current cell, fails if we are out of the object */
		if (!SubObjGetCellAttenuation( trackingPhotonPtr->sliceIndex,  trackingPhotonPtr->xIndex, 
				 trackingPhotonPtr->yIndex, photonEnergy, &cellAttenuation)) {
			PhgAbort("Attempt to get attenuation for cell not in subobject (phoTrkCalcFreePaths)",
				true);
		}
		
		/* Save our cell information */
		{
			phoTrkCellInfoTable[phoTrkCellsInUse].distance = 
				nextDist - distanceTracked;
				
			phoTrkCellInfoTable[phoTrkCellsInUse].attenuation = 
				cellAttenuation;
				
			phoTrkCellInfoTable[phoTrkCellsInUse].freePaths = 
				(phoTrkCellInfoTable[phoTrkCellsInUse].distance
				* cellAttenuation);
				
			phoTrkCellInfoTable[phoTrkCellsInUse].xIndex = 
				trackingPhotonPtr->xIndex;
				
			phoTrkCellInfoTable[phoTrkCellsInUse].yIndex = 
				trackingPhotonPtr->yIndex;
	
			phoTrkCellInfoTable[phoTrkCellsInUse].sliceIndex = 
				trackingPhotonPtr->sliceIndex;
		}
		
		/* Add up our free paths */
		*freePathsPtr += phoTrkCellInfoTable[phoTrkCellsInUse].freePaths;
		
		/* Increment to next cell */
		phoTrkCellsInUse++;
		
		#ifdef PHG_DEBUG
			if (phoTrkCellsInUse == phoTrkMaxCellCount)
				PhgAbort("Attempt to index past cell info table size (phoTrkCalcFreePaths).",
					true);
		#endif
		
		/* Increment distance tracked, then see if we need to keep going */
		distanceTracked = nextDist;
		
		/* See if we have NOT traveled the total distance */
		if (!PhgMathRealNumAreEqual(distanceTracked, distanceToTravel, 
		        -7, 0, 0, 0)){
		
			/* See which value we used */
			if (nextDist == distToNextX) {
				/* From here on out each x crossing is the same distance */
				distToNextX += generalDistToX;
				
				/* Increment our x index */
				 trackingPhotonPtr->xIndex += ((direction.cosine_x >= 0) ? 1 : -1);
				
				#ifdef PHG_DEBUG
					if ((LbTwoByte)  trackingPhotonPtr->xIndex < 0)
						PhgAbort("Failed to realize photon reached object x boundary (phoTrkCalcFreePaths).",
							true);
				#endif
			}
			else if (nextDist == distToNextY) {
				/* From here on out each y crossing is the same distance */
				distToNextY += generalDistToY;
				
				/* Increment our y index */
				/* Notice that Y goes top to bottom so it is reveresed from X */
				 trackingPhotonPtr->yIndex += ((direction.cosine_y >= 0) ? -1 : 1);
				
				#ifdef PHG_DEBUG
					if ((LbTwoByte)  trackingPhotonPtr->yIndex < 0)
						PhgAbort("Failed to realize photon reached object y boundary (phoTrkCalcFreePaths).",
							true);
				#endif
	
			}
			else {
				
				/* We are going to a new slice so update the values */
				newSliceIndex =  trackingPhotonPtr->sliceIndex +
					((direction.cosine_z >= 0) ? 1 : -1);
	
				/* See if we went out the end of the object */
				if ((newSliceIndex < 0) || ((LbUsFourByte)newSliceIndex == SubObjNumSlices)) {
					/* We are done */
					break;
				}
	
				/* If we didn't break, we entered a new slice of the object; recalculate slice variables */						
				{
					/* Compute the new position */
					PhoTrkProject(startingPosPtr, &direction, distanceTracked, &newPosition);

					/* Update slice parameters */
					PhoTrkEnterNewSlice(newPosition, direction,
				   		 trackingPhotonPtr->sliceIndex, distanceTracked,
						&distToNextX, &distToNextY, &distToNextZ,
						&generalDistToX, &generalDistToY, &generalDistToZ,
						&(trackingPhotonPtr->sliceIndex), 
						&(trackingPhotonPtr->xIndex),
						&(trackingPhotonPtr->yIndex));
				}
			}
		}
		else {
			/* If we are equal within tolerance, we are equal */
			distanceTracked = distanceToTravel;
		}
	} while (distanceTracked < distanceToTravel);	
}


/*********************************************************************************
*
*			Name:		PhoTrkCalcRange
*
*			Summary:	Determine the distance a photon will travel due to 
*						positron range effects and compute the new decay location.
*
*			Arguments:
*				double				freePaths			- The free paths to travel.
*				PHG_Position		startingPos			- The starting position.
*				PHG_Direction		direction			- The direction of travel.
*				double				*finalPosPtr		- The photon's energy.
*				Boolean				*discard			- Flag for escaping.
*				LbFourByte			*sliceIndex			- Position index.
*				LbFourByte			*xIndex				- Position index.
*				LbFourByte			*yIndex				- Position index.
*
*			Function return: None.
*
*********************************************************************************/
void	PhoTrkCalcRange(double freePaths,
				PHG_Position startingPos,
				PHG_Direction direction, 
				PHG_Position *finalPosPtr,
				Boolean	*discard,
				LbFourByte	*sliceIndex,
				LbFourByte	*xIndex,
				LbFourByte	*yIndex)
{
	double			cellAttenuationR;		/* Attenuation of current cell */
	double			distanceTrackedR;		/* Total distance tracked */
	double			initialZDistanceR;		/* Initial distance from current z to z boundary */
	double			initialXDistanceR;		/* Initial distance from current x to x boundary */
	double			initialYDistanceR;		/* Initial distance from current y to y boundary */
	double			distToNextXR;			/* Distance to next x crossing */
	double			distToNextYR;			/* Distance to next y crossing */
	double			distToNextZR;			/* Distance to next z crossing */
	double			generalDistToXR;			/* General distance to move relative to x axis */
	double			generalDistToYR;			/* General distance to move relative to y axis */
	double			generalDistToZR;			/* General distance to move relative to z axis */
	double			nextDistR;				/* Next incremental move */
	double			distanceR;				/* Distance about to move */
	double			freePathsForMoveR;		/* Free paths used to make the next move */
	double			freePathsUsedR;			/* Free paths used to make the previous moves */
	LbFourByte		newSliceIndexR;			/* New slice index for entering new slice */
	PHG_Position	newPositionR;			/* Used to update slice info */
	double			distToObjectSurfaceR;	/* distance to object surface */
	

	/* Clear counters */
	distanceTrackedR = 0.0;
	distanceR = 0.0;
	freePathsUsedR = 0.0;
	freePathsForMoveR = 0.0;
	
	/* Clear distance marker */
	nextDistR = 0.0;
	
	/* Clear discard flag */
	*discard = false;
	
	/* Get initial distance to voxel wall */
	SubObjGetInnerCellDistance(&startingPos, &direction, *sliceIndex,
		 *xIndex,  *yIndex, &initialXDistanceR, &initialYDistanceR, &initialZDistanceR);
	
	/* Avoid div by zero from direction cosines */
	{
		if ((direction.cosine_x <= 0.0000001) && 
				(direction.cosine_x >= -0.0000001)) {				
			direction.cosine_x = 
				(0.0000001 * ((direction.cosine_x < 0) ? -1 : 1));
		}
		if ((direction.cosine_y <= 0.0000001) && 
				(direction.cosine_y >= -0.0000001)) {				
			direction.cosine_y = 
				(0.0000001 * ((direction.cosine_y < 0) ? -1 : 1));
		}
		if ((direction.cosine_z <= 0.0000001) && 
				(direction.cosine_z >= -0.0000001)) {				
			direction.cosine_z = 
				(0.0000001 * ((direction.cosine_z < 0) ? -1 : 1));
		}
	}

	/* Calculate initial distances  (normalized to positive value) */
	{
		distToNextXR = initialXDistanceR/direction.cosine_x;
		distToNextYR = initialYDistanceR/direction.cosine_y;
		distToNextZR = initialZDistanceR/direction.cosine_z;
	}
	
	/* Calculate general distances */
	{
		generalDistToXR = SUBOBJGetSliceAttVoxelWidth(*sliceIndex)/fabs(direction.cosine_x);
		generalDistToYR = SUBOBJGetSliceAttVoxelHeight(*sliceIndex)/fabs(direction.cosine_y);
		generalDistToZR = SUBOBJGetSliceVoxelDepth(*sliceIndex)/fabs(direction.cosine_z);
	}
	
	do { /* Tracking Loop */

		/* Get shortest distance to next cell boundary */
		nextDistR = PHOTRKGetMinDeltaT(distToNextXR, distToNextYR, distToNextZR);
		
		/* Calculate distance to object's cylindrical surface */
		CylPosCalcDistanceToObjectSurface(&startingPos, &direction, &distToObjectSurfaceR);

		/* See if this will take us out of the object.
			If so, adjust distance to go onto object cylinder
		*/
		if (nextDistR >= distToObjectSurfaceR) {

			/* Set distance to put us onto cylinder */
			nextDistR = distToObjectSurfaceR;		
		}

		/* Get attenutation of current cell, fails if we are out of the object */
		if (!SubObjGetCellAttenuation(*sliceIndex, *xIndex, *yIndex, 1000.0, &cellAttenuationR)) {
			PhgAbort("Attempt to get attenuation for cell not in subobject (PhoTrkCalcRange)",
				true);
		}
		
		/* Compute distance to move through this cell */
		distanceR = nextDistR - distanceTrackedR;
		
		/* Store free paths moved thus far */
		freePathsUsedR = freePathsForMoveR;

		/* Compute free paths used for this move */
		freePathsForMoveR += (distanceR * cellAttenuationR);
		
		/* See if this move will exhaust free paths available */
		if (freePathsForMoveR >= freePaths) {

			/* Truncate distance to move */
			distanceR = distanceTrackedR + (freePaths-freePathsUsedR)/cellAttenuationR;

			/* See if this will take us out of the object.
				If so, bolt
			*/
			if (distanceR >= distToObjectSurfaceR) {

				*discard = true;
				break;	
			}
			
			/* Compute the new position */
			PhoTrkProject(&startingPos, &direction, distanceR, finalPosPtr);
			
			/* Break out of this loop, we are done tracking */
			break;
		}

		/* Increment distance tracked, then see if we need to keep going */
		distanceTrackedR = nextDistR;

		/* See if this will take us out of the object.
			If so, bolt
		*/
		if (distanceTrackedR >= distToObjectSurfaceR) {

			*discard = true;
			break;	
		}
		
		/* See which value we used */
		if (nextDistR == distToNextXR) {
			/* From here on out each x crossing is the same distance */
			distToNextXR += generalDistToXR;
			
			/* Increment our x index */
			*xIndex += ((direction.cosine_x >= 0) ? 1 : -1);
			
			/* See if we leave the object */
			if ((*xIndex < 0) || ((LbUsFourByte)(*xIndex) >= SubObjObject[*sliceIndex].attNumXBins)) { 
				*discard = true;
				break;
			}
			
		}
		else if (nextDistR == distToNextYR) {
			/* From here on out each y crossing is the same distance */
			distToNextYR += generalDistToYR;
			
			/* Increment our y index */
			/* Notice that Y goes top to bottom so it is reveresed from X */
			 *yIndex += ((direction.cosine_y >= 0) ? -1 : 1);
			
			
			/* See if we leave the object */
			if ((*yIndex < 0) || ((LbUsFourByte)(*yIndex) >= SubObjObject[*sliceIndex].attNumYBins)) { 
				*discard = true;
				break;
			}
			

		}
		else {
			
			/* We are going to a new slice so update the values */
			newSliceIndexR =  *sliceIndex +
				((direction.cosine_z >= 0) ? 1 : -1);

			/* See if we went out the end of the object */
			if ((newSliceIndexR < 0) || ((LbUsFourByte)newSliceIndexR == SubObjNumSlices)) {
				*discard = true;
				break;
			}

			/* If we didn't break, we entered a new slice of the object; recalculate slice variables */						
			{
				/* Compute the new position */
				PhoTrkProject(&startingPos, &direction, distanceTrackedR, &newPositionR);

				/* Update slice parameters */
				PhoTrkEnterNewSlice(newPositionR, direction,
			   		 *sliceIndex, distanceTrackedR,
					&distToNextXR, &distToNextYR, &distToNextZR,
					&generalDistToXR, &generalDistToYR, &generalDistToZR,
					sliceIndex, 
					xIndex,
					yIndex);
			}
		}
	} while (true);		
}

/*********************************************************************************
*
*			Name:		phoTrkCalcCritZoneFreePaths
*
*			Summary:	Determine the cumulative attenuation the photon will
*						encounter on its trip through the critical zone. It is
*						assumed that the photon will enter the critical zone
*						before this function gets called.
*
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr		- The photon.
*				PHG_Intersection	*intersectionPtr		- The photon's critical zone intersection.
*				double				*freePathsToEnterPtr	- The free paths used getting into the zone.
*				double				*freePathsToExitPtr		- The free paths used getting to the outer edge.
*
*			Function return: None.
*
*********************************************************************************/
void phoTrkCalcCritZoneFreePaths(PHG_TrackingPhoton	*trackingPhotonPtr,
						PHG_Intersection *intersectionPtr,
						double *freePathsToEnterPtr,
						double *freePathsToExitPtr)
{
	LbUsFourByte	cellIndex;			/* Current cell for traversal */
	double			distanceTraveled;	/* Distance travelled through cell list */

	do { /* Process Loop */

		/* If we have a voxel (cell) list, then compute free paths from it */
		if (phoTrkCellsInUse != 0) {
			
			distanceTraveled = 0;
			*freePathsToExitPtr = 0;
			
			/* Loop through cell table until freePathsToTravel have been used up */
			for (cellIndex = 0; cellIndex < phoTrkCellsInUse; cellIndex++) {
				 if ((distanceTraveled + phoTrkCellInfoTable[cellIndex].distance)
						> intersectionPtr->distToExit) {
						
					/* Add freepaths used within this cell */
					*freePathsToExitPtr += 
						((intersectionPtr->distToExit -
						distanceTraveled) *
						phoTrkCellInfoTable[cellIndex].attenuation);

					/* We are done */
					break;
					
				}
				else {
					distanceTraveled += phoTrkCellInfoTable[cellIndex].distance;
					*freePathsToExitPtr += phoTrkCellInfoTable[cellIndex].freePaths;
				}
			}
		}
		else {
		
			/* We didn't have a cell list so track to exit */
			phoTrkCalcFreePaths(trackingPhotonPtr, &(trackingPhotonPtr->location),
					intersectionPtr->distToExit,
					trackingPhotonPtr->angle,
					trackingPhotonPtr->energy, freePathsToExitPtr);
		}
			
		/* If we are already in the critical zone then set freePathsToEnter = 0,
			otherwise track to the entrance of the zone
		*/
		if (intersectionPtr->distToEnter == 0.0) {
			
			/* Nothing used to get there */
			*freePathsToEnterPtr = 0.0;
		}
		else {

			distanceTraveled = 0;
			*freePathsToEnterPtr = 0;
			
			/* Loop through cell table until we enter the critical zone */
			for (cellIndex = 0; cellIndex < phoTrkCellsInUse; cellIndex++) {
				 if ((distanceTraveled + phoTrkCellInfoTable[cellIndex].distance)
						> intersectionPtr->distToEnter) {
						
					/* Add freepaths used within this cell */
					*freePathsToEnterPtr += 
						((intersectionPtr->distToEnter - 
						distanceTraveled) *
						phoTrkCellInfoTable[cellIndex].attenuation);
						
					/* We are done */
					break;
						
				}
				else {
					distanceTraveled += phoTrkCellInfoTable[cellIndex].distance;
					*freePathsToEnterPtr += phoTrkCellInfoTable[cellIndex].freePaths;
				}
			}
			#ifdef PHG_DEBUG
				/* If we got here because we ran out of cells, it is an error */
				if (cellIndex == phoTrkCellsInUse){
					PhgAbort("\nphoTrkCalcCritZoneFreePaths: Exhausted cell list while computing free paths to enter!\n",
						true);
				}
			#endif
		}
		
	} while (false);

	#ifdef PHG_DEBUG
		if (*freePathsToExitPtr < *freePathsToEnterPtr)
			PhgAbort("Invalid freePaths computation (phoTrkCalcCritZoneFreePaths).",
				true);
	#endif
}

/*********************************************************************************
*
*			Name:		PhoTrkEnterNewSlice
*
*			Summary:	Update tracking variables when entering a new slice.
*
*			Arguments:
*				PHG_Position		position			- The photon's position upon entry.
*				PHG_Direction		direction			- The photon's direction.
*				LbUsFourByte		sliceIndex			- The slice Index on entry.
*				LbUsFourByte		distanceTraveled	- The distance traveled so far.
*				double				*nextXPtr			- Next distance in x direction
*				double				*nextYPtr			- Next distance in y direction
*				double				*nextZPtr			- Next distance in z direction
*				double				*generalXPtr		- General distance to cell wall for new slice .
*				double				*generalYPtr		- General distance to cell wall for new slice .
*				double				*generalZPtr		- General distance to cell wall for new slice .
*				LbFourByte			*newSliceIndexPtr	- New slice Index.
*				LbFourByte			*newXIndexPtr		- New x Index.
*				LbFourByte			*newYIndexPtr		- New y index.
*
*			Function return: none
*
*********************************************************************************/
void PhoTrkEnterNewSlice(PHG_Position position, PHG_Direction direction,
				LbUsFourByte sliceIndex, double distanceTraveled,
				double *nextXPtr, double *nextYPtr, double *nextZPtr,
			    double *generalXPtr, double *generalYPtr, double *generalZPtr,
				LbFourByte *newSliceIndexPtr, LbFourByte *newXIndexPtr,
				LbFourByte *newYIndexPtr)
{
	double			initialXDist;		/* Initial distance to x crossing */
	double			initialYDist;		/* Initial distance to y crossing */
	double			initialZDist;		/* Initial distance to z crossing */
	double			newVoxelDepth;		/* Depth of new slice */
	double			newVoxelWidth;		/* Width of new slice voxels */
	double			newVoxelHeight;		/* Height of new voxel */
	double			currVoxelWidth;		/* Width of current slice voxels */
	double			currVoxelHeight;	/* Height of current slice voxels */
    double			currVoxelDepth;		/* Current voxel depth */

	/* Get parameters for current slice */
	currVoxelDepth = SUBOBJGetSliceVoxelDepth(sliceIndex);
	currVoxelWidth = SUBOBJGetSliceAttVoxelWidth(sliceIndex);
	currVoxelHeight = SUBOBJGetSliceAttVoxelHeight(sliceIndex);

	/* Compute new slice index */
	*newSliceIndexPtr = (sliceIndex + ((direction.cosine_z >= 0) ? 1 : -1));

	/* Get values for new slice */
	newVoxelDepth = SUBOBJGetSliceVoxelDepth(*newSliceIndexPtr);
	newVoxelWidth = SUBOBJGetSliceAttVoxelWidth(*newSliceIndexPtr);
	newVoxelHeight = SUBOBJGetSliceAttVoxelHeight(*newSliceIndexPtr);

	/* See if we need to update for z direction */
	if (currVoxelDepth != newVoxelDepth) {
		
		/* Calculate distance to edge of slice */
		if (direction.cosine_z >= 0) {
			initialZDist = fabs((SUBOBJGetSliceMaxZ(*newSliceIndexPtr) - 
				position.z_position)/direction.cosine_z);
		}
		else {
			initialZDist = fabs((SUBOBJGetSliceMinZ(*newSliceIndexPtr) - 
				position.z_position)/direction.cosine_z);
		}
		
		/* Store the new general Z distance */
		*generalZPtr = fabs(SUBOBJGetSliceVoxelDepth(*newSliceIndexPtr)/
			direction.cosine_z);

		/* Compute the next Z distance */
		*nextZPtr += initialZDist;
	}
	else {
		/* We entered because of change in slice so nextZPtr must get updated */
		*nextZPtr += *generalZPtr;
	}

	/* See if we need to update for x axis */
	if ((SUBOBJGetSliceMinX(sliceIndex) != SUBOBJGetSliceMinX(*newSliceIndexPtr))
		|| (currVoxelWidth != newVoxelWidth)) {

		/* Calculate the x index */
		*newXIndexPtr = (LbFourByte) floor(((position.x_position-
			SUBOBJGetSliceMinX(*newSliceIndexPtr)) *
			SUBOBJGetSliceAttNumXBins(*newSliceIndexPtr))/
			SUBOBJGetSliceWidth(*newSliceIndexPtr));
	
   		/* Calculate the x distance */
		if (direction.cosine_x >= 0) {
			initialXDist = fabs(((SUBOBJGetSliceMinX(*newSliceIndexPtr) + 
						 ((*newXIndexPtr + 1) * 
						  newVoxelWidth)) - 
						  position.x_position)/direction.cosine_x);
		}
		else {
			initialXDist = fabs(((SUBOBJGetSliceMinX(*newSliceIndexPtr) + 
				(*newXIndexPtr * newVoxelWidth)) 
				- position.x_position)/direction.cosine_x);
		}

		/* Compute the next X distance */
		*nextXPtr  = distanceTraveled + initialXDist;

		/* Compute the general x distance */
		*generalXPtr = fabs(newVoxelWidth/direction.cosine_x);
	}

	/* See if we need to update for the y direction */
	if (currVoxelHeight != newVoxelHeight){

		/* Calculate the y index */
		*newYIndexPtr = (LbFourByte) floor(((SUBOBJGetSliceMaxY(*newSliceIndexPtr) -
			position.y_position) * SUBOBJGetSliceAttNumYBins(*newSliceIndexPtr))
			/SUBOBJGetSliceSliceHeight(*newSliceIndexPtr));
		
		/* Calculate the y distance */
		if (direction.cosine_y >= 0) {
			initialYDist = ((SUBOBJGetSliceMaxY(*newSliceIndexPtr) -
				(*newYIndexPtr * newVoxelHeight)) - position.y_position);
				
			/* Do intermediate check */
			#ifdef PHG_DEBUG
				if ((initialYDist < 0) || (initialYDist > newVoxelHeight)) {
					PhgAbort("Invalid computation for y distance (PhoTrkEnterNewSlice)",
						true);
				}
			#endif
			
			/* Normalize to direction vector */
			initialYDist = 	initialYDist/direction.cosine_y;
			
		}
		else {
			initialYDist = ((SUBOBJGetSliceMaxY(*newSliceIndexPtr) -
				((*newYIndexPtr+1) * newVoxelHeight)) - position.y_position);

			#ifdef PHG_DEBUG
				if ((initialYDist > 0) || (initialYDist < -newVoxelHeight)) {
					PhgAbort("Invalid computation for y distance (PhoTrkEnterNewSlice)",
						true);
				}
			#endif
			
			initialYDist = 	initialYDist/direction.cosine_y;
		}

		/* Compute the next Y distance */
		*nextYPtr = distanceTraveled + initialYDist;

		/* Compute the general y distance */
		*generalYPtr = fabs(newVoxelHeight/direction.cosine_y);
	}

}

/*********************************************************************************
*
*			Name:		PhoTrkCalcNewPosition
*
*			Summary:	Determine the new_position for a photon at starting_position 
*						with photon_energy which must travel freepath_length 
*						through the attenuation subobject.  Clip the distance 
*						traveled to the limit cylinder, if necessary.
*
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The photon we are tracking.
*				PHG_Position		startingPosition	- The photon's starting position.
*				PHG_Direction		direction			- The photon's starting direction.
*				double				photonEnergy		- The current energy.
*				double				freePathLength		- Free paths to move.
*				PHG_Position		*newPosition		- The new position.
*				double				*distTraveled		- The distance the photon traveled.
*				double				*freePathsUsed		- The actual free paths used.
*
*			Function return: Action to take based on calculation.
*
*********************************************************************************/
PhoTrkActionTy PhoTrkCalcNewPosition(PHG_TrackingPhoton	*trackingPhotonPtr,
				PHG_Position startingPosition, PHG_Direction direction,
				double photonEnergy, double freePath_Length,  PHG_Position *newPosition,
				double *distTraveled, double *freePathsUsed)
{
	double			attenuation;			/* Attenuation for a given cell */
	double			distToObjectSurface;	/* Distance to object surface */
	double			distToTargetSurface=0;	/* Distance to target surface */
	double			initialZDistance;		/* Initial distance from current z to z boundary */
	double			initialXDistance;		/* Initial distance from current x to x boundary */
	double			initialYDistance;		/* Initial distance from current y to y boundary */
	double			limitZMin;				/* Minimum z value for limit cylinder */
	double			limitZMax;				/* Maximum z value for limit cylinder */
	double			distToNextX;			/* Distance to next x crossing */
	double			distToNextY;			/* Distance to next y crossing */
	double			distToNextZ;			/* Distance to next z crossing */
	double			generalDistToX;			/* General distance to move relative to x axis */
	double			generalDistToY;			/* General distance to move relative to x axis */
	double			generalDistToZ;			/* General distance to move relative to x axis */
	double			distanceTracked;		/* Total distance tracked */
	double			nextDist;				/* Next incremental move */
	double			nextFreePathsToUse;		/* Amount of free paths next move will use */
	LbFourByte		newSliceIndex;			/* Computed slice index for entering new slices */
	LbFourByte		cellIndex;				/* Computed slice index for entering new slices */
	PhoTrkActionTy	action;					/* Return value */
	

	/* Assume the photon will interact */
	action = PhoTrkInteract;
	
	/* Clear "distance traveled" counter */
	*distTraveled = 0.0;
	
	/* Clear free paths used counter */
	*freePathsUsed = 0.0;

	/* Set local limit vars (could be eliminated, but names would be long) */
	limitZMin = CylPosLimitCylinder.zMin;
	limitZMax = CylPosLimitCylinder.zMax;
		
	/* Set newPosition = startingPosition, note that this is important */
	*newPosition = startingPosition;

	do { /* Process Loop. Break occurs when photon interacts or leaves object */
		
		/* If we have a voxel (cell) list, then use it to update position */
		if (phoTrkCellsInUse != 0) {

			/* Clearloop variables */
			distanceTracked = 0;
			*freePathsUsed = 0.0;
			
			/* Loop through cell table until freePathsToTravel have been used up */
			for (cellIndex = 0; cellIndex < phoTrkCellsInUse; cellIndex++) {
				 
				 /* See if we will use available freepaths on this cell */
				 if ((*freePathsUsed + phoTrkCellInfoTable[cellIndex].freePaths)
						> freePath_Length) {
						
					/* Only use necessary amount of cell */
					distanceTracked += (freePath_Length - *freePathsUsed)/
						phoTrkCellInfoTable[cellIndex].attenuation;
					
					/* Indicate we used all free paths */
					*freePathsUsed = freePath_Length;
		
					/* Calculate new position and object indexes */
					{
						PhoTrkProject(&startingPosition, &direction, distanceTracked, newPosition);
							
						trackingPhotonPtr->xIndex = phoTrkCellInfoTable[cellIndex].xIndex;
						trackingPhotonPtr->yIndex = phoTrkCellInfoTable[cellIndex].yIndex;
						trackingPhotonPtr->sliceIndex = phoTrkCellInfoTable[cellIndex].sliceIndex;
						
					}

					/* We are done */
					break;
						
				}
				else {
					distanceTracked += phoTrkCellInfoTable[cellIndex].distance;
					*freePathsUsed += phoTrkCellInfoTable[cellIndex].freePaths;
				}
			}
			/* If we got here without running out of cells, we interacted so we are done */
			if (cellIndex < phoTrkCellsInUse){
			
				break;
			}
			else {
			
				/*	We went through all cells, so we
					may, or may not be on the target cylinder depending
					on where the cell list was created. We need to test for
					distance to object cylinder to find out.
				*/
			
				/* Calculate distance to object's cylindrical surface */
				CylPosCalcDistanceToObjectSurface(&startingPosition, &direction, &distToObjectSurface);

				/* If we moved to surface, handle as detection, note in this case we may have moved 
					to the end of the object cylinder if the last cell we traveled through was
					in the z direction so this test must include the end of the object cylinder */
				if ((PhgMathRealNumAreEqual(distanceTracked, distToObjectSurface, -7, 
				        0, 0, 0)) || (newPosition->z_position <= CYLPOSGetObjZMin()) ||
						(newPosition->z_position >= CYLPOSGetObjZMax())){
				

					/* Calculate new position (this may or may not be redundant based on 
						previous path above)
					*/
					{
						PhoTrkProject(&startingPosition, &direction, distanceTracked, newPosition);
							
						trackingPhotonPtr->xIndex = phoTrkCellInfoTable[phoTrkCellsInUse-1].xIndex;
						trackingPhotonPtr->yIndex = phoTrkCellInfoTable[phoTrkCellsInUse-1].yIndex;
						trackingPhotonPtr->sliceIndex = phoTrkCellInfoTable[phoTrkCellsInUse-1].sliceIndex;
						
					}
					
					/* See if we are beyond end of limit cylinder */
					if ((newPosition->z_position <= limitZMin) ||
							(newPosition->z_position >= limitZMax)) {
						
						action = PhoTrkDiscard;
						break;
					}
					
					/* See if we are outside of our acceptance angle */
					if (fabs(direction.cosine_z) > PHGGetSineOfAccAngle()) {
						
						/* Discard it */
						action = PhoTrkDiscard;
						break;
					}
						
					/* Project to target cylinder if it is not equal to object cylinder */
					if (CylPosGetTargetRadius() <= CylPosGetObjectRadius()) {

						/* Set action to detect */
						action = PhoTrkDetect;
												
						/* Get out of routine */
						break;
					}
					else if (!CylPosProjectToTargetCylinder(newPosition, &direction,
							 &distToTargetSurface)) {
						action = PhoTrkDiscard;
						break;
					}
					else {
					
						/* Set action to detect */
						action = PhoTrkDetect;
						
						/* Update distance tracked */
						distanceTracked += distToTargetSurface;
						
						/* Get out of routine */
						break;
					}
				}
				else {
				
					/*	Our critical zone is smaller than the object cylinder,
						so we have run out of cells before interacting or reaching
						the surface.
						We need to update our starting position. We then use 
						the originaltracking method.
					*/
					PhoTrkProject(&startingPosition, &direction, distanceTracked, newPosition);
						
					trackingPhotonPtr->xIndex = phoTrkCellInfoTable[phoTrkCellsInUse-1].xIndex;
					trackingPhotonPtr->yIndex = phoTrkCellInfoTable[phoTrkCellsInUse-1].yIndex;
					trackingPhotonPtr->sliceIndex = phoTrkCellInfoTable[phoTrkCellsInUse-1].sliceIndex;
					
					/* Now it is possible that we have reached the end of the cylinder */

					/* If we moved to surface, handle as detection, note in this case we may have moved 
						to the end of the object cylinder if the last cell we traveled through was
						in the z direction so this test must include the end of the object cylinder */
					if ((newPosition->z_position <= CYLPOSGetObjZMin()) ||
							(newPosition->z_position >= CYLPOSGetObjZMax())){
							
						action = PhoTrkDiscard;
						break;
					}

					/* We made it to here so we'll update our starting position
						and continue tracking the hard way
					*/
					startingPosition = *newPosition;
					
					/* Reset free paths counter */
					freePath_Length -= *freePathsUsed;
					 *freePathsUsed = 0.0;
				}
			}
		}
		
		/* Calculate distance to object's cylindrical surface */
		CylPosCalcDistanceToObjectSurface(&startingPosition, &direction, &distToObjectSurface);
		

		/* It is possible to start out right on the surface so check for that */
		if (distToObjectSurface == 0.0) {

			/* See if we are beyond end of limit cylinder */
			if ((startingPosition.z_position <= limitZMin) ||
					(startingPosition.z_position >= limitZMax)) {

				action = PhoTrkDiscard;
				break;
			}
			
			/* Check our acceptance angle */
			if (fabs(direction.cosine_z) > PHGGetSineOfAccAngle()) {
				
				action = PhoTrkDiscard;
				break;
			}
			
			/* If target = object then this is a detection */
			if (CylPosGetTargetRadius() <= CylPosGetObjectRadius()) {

				/* Set action to detect */
				action = PhoTrkDetect;
										
				/* Get out of routine */
				break;
			}
			/* See if photon projects to target cylinder within limits */
			else if (!CylPosProjectToTargetCylinder(&startingPosition, &direction, &distToTargetSurface)) {

				action = PhoTrkDiscard;
				break;
			}
			else {
			
				/* Compute distance to target surface */
				CylPosCalcDistanceToTargetSurface(&startingPosition,
					&direction, &distanceTracked);
					
				/* Photon hit target within acceptance, it is done */
				action = PhoTrkDetect;
				break;
			}

		}
		
		/* Get initial distance to cell wall */
		SubObjGetInnerCellDistance(&startingPosition, &direction, trackingPhotonPtr->sliceIndex,
			trackingPhotonPtr->xIndex, trackingPhotonPtr->yIndex,
			&initialXDistance, &initialYDistance, &initialZDistance);
		
		/* Avoid div by zero from direction cosines */
		{
			if ((direction.cosine_x <= 0.0000001) && 
				(direction.cosine_x >= -0.0000001)) {				
				direction.cosine_x = 
				(0.0000001 * ((direction.cosine_x < 0) ? -1 : 1));
			}
			if ((direction.cosine_y <= 0.0000001) && 
				(direction.cosine_y >= -0.0000001)) {				
				direction.cosine_y = 
				(0.0000001 * ((direction.cosine_y < 0) ? -1 : 1));
			}
			if ((direction.cosine_z <= 0.0000001) && 
				(direction.cosine_z >= -0.0000001)) {				
				direction.cosine_z = 
				(0.0000001 * ((direction.cosine_z < 0) ? -1 : 1));
			}
		}
		
		/* Calculate initial distances (normalized to positive value) */
		{
			distToNextX = initialXDistance/direction.cosine_x;
			distToNextY = initialYDistance/direction.cosine_y;
			distToNextZ = initialZDistance/direction.cosine_z;
		}
		
		/* Calculate general distances */
		{
			generalDistToX = SUBOBJGetSliceAttVoxelWidth(trackingPhotonPtr->sliceIndex)/
				fabs(direction.cosine_x);
			generalDistToY = SUBOBJGetSliceAttVoxelHeight(trackingPhotonPtr->sliceIndex)/
				fabs(direction.cosine_y);
			generalDistToZ = SUBOBJGetSliceVoxelDepth(trackingPhotonPtr->sliceIndex)/
				fabs(direction.cosine_z);
		}
		
		/* Zero loop variables */
		distanceTracked = 0;
		nextDist = 0;
		nextFreePathsToUse = 0;
		
		/* Start process Loop. NOTE Only exits with break statement */		
		do {
				
			/* Get shortest distance to next cell boundary */
			nextDist = PHOTRKGetMinDeltaT(distToNextX, distToNextY, distToNextZ);
			
			/* See if this will take us out of the object.
				If so, adjust distance to go onto object cylinder
			*/
			if ((nextDist > distToObjectSurface) || (PhgMathRealNumAreEqual(nextDist, distToObjectSurface, -7, 0,0,0) == true)) {

				/* Set distance to put us onto cylinder */
				nextDist = distToObjectSurface;		
			}
			
			/* Get attenutation of current cell, fails if we are out of the object */
			if (!SubObjGetCellAttenuation(trackingPhotonPtr->sliceIndex,
					trackingPhotonPtr->xIndex, trackingPhotonPtr->yIndex, photonEnergy, &attenuation)) {
				PhgAbort("Attempt to get attenuation for cell not in subobject (PhoTrkCalcNewPosition)",
					true);
			}
	
			/* Calculate amount of free paths to be used */
			nextFreePathsToUse = fabs(nextDist - distanceTracked) * 
				attenuation;
			
			/* See if we will interact on this move */
			if ((freePath_Length != -1) && 
					((*freePathsUsed + nextFreePathsToUse) >= freePath_Length)) {

				/* Adjust distance tracked to end of freePath_Length */
				nextDist = distanceTracked + 
					((freePath_Length - *freePathsUsed)/attenuation);
				distanceTracked = nextDist;
				*freePathsUsed = freePath_Length;
	
				/* Calculate new position */
				PhoTrkProject(&startingPosition, &direction, nextDist, newPosition);

				
				/* We are done */
				break;
			}
			else {
			

				/* If we moved to surface, project to target */
				if (PhgMathRealNumAreEqual(nextDist, distToObjectSurface, -7, 
				        0, 0, 0)){

					/* Move the photon to its new location */
					PhoTrkProject(&startingPosition, &direction, nextDist, newPosition);

					/* Save exact value */
					distanceTracked = distToObjectSurface;
					
		
					/* See if we are beyond end of limit cylinder */
					if ((newPosition->z_position <= limitZMin) ||
							(newPosition->z_position >= limitZMax)) {
						
						action = PhoTrkDiscard;
						break;
					}
					
					/* Check our acceptance angle */
					if (fabs(direction.cosine_z) > PHGGetSineOfAccAngle()) {
						
						
						action = PhoTrkDiscard;
						break;
					}
					
					/* If target = object then this is a detection */	
					if (CylPosGetTargetRadius() <= CylPosGetObjectRadius()) {
						
						/* Set action to detect */
						action = PhoTrkDetect;
						
						/* Update free paths */
						*freePathsUsed += nextFreePathsToUse;
						
						/* Update distance tracked */
						distanceTracked += distToTargetSurface;
						
						/* Get out of routine */
						break;
					}
						/* Project to target cylinder */
					else if (!CylPosProjectToTargetCylinder(newPosition, &direction,
							 &distToTargetSurface)) {
						action = PhoTrkDiscard;
						break;
					}
					else {
					
						/* Photon reached target cylinder within limits */
						action = PhoTrkDetect;
						*freePathsUsed += nextFreePathsToUse;
						
						/* Update distance tracked */
						distanceTracked += distToTargetSurface;
						
						break;
					}
				}
			}
			
			/* Update free paths used */
			*freePathsUsed += nextFreePathsToUse;
			
			/* Update distance traveled */
			distanceTracked = nextDist;
			
			/* See which value we used */
			if (nextDist == distToNextX) {
				distToNextX += generalDistToX;
				trackingPhotonPtr->xIndex += ((direction.cosine_x >= 0) ? 1 : -1);
				
				#ifdef PHG_DEBUG
				if ((trackingPhotonPtr->xIndex < 0) || ((LbUsFourByte)trackingPhotonPtr->xIndex >= SubObjObject[trackingPhotonPtr->sliceIndex].attNumXBins)) {
					LbInPrintf("xIndex out of range %d, [%d, %d]\n", trackingPhotonPtr->xIndex, 0, SubObjObject[trackingPhotonPtr->sliceIndex].attNumXBins-1);
				}
				#endif				
				
			}
			else if (nextDist == distToNextY) {
				distToNextY += generalDistToY;
				
				/* Notice that Y goes top to bottom so it is reversed from X */
				trackingPhotonPtr->yIndex += ((direction.cosine_y >= 0) ? -1 : 1);
				
				#ifdef PHG_DEBUG
				if ((trackingPhotonPtr->yIndex < 0) || ((LbUsFourByte)trackingPhotonPtr->yIndex >= SubObjObject[trackingPhotonPtr->sliceIndex].attNumYBins)) {
					LbInPrintf("yIndex out of range %d, [%d, %d]\n", trackingPhotonPtr->yIndex, 0, SubObjObject[trackingPhotonPtr->sliceIndex].attNumYBins-1);
				}
				#endif				
			}
			else {
				
				/* Compute the new position */
				PhoTrkProject(&startingPosition, &direction, distanceTracked, newPosition);

				/* We are going to a new slice so update the values */
				newSliceIndex = trackingPhotonPtr->sliceIndex + 
					((direction.cosine_z >= 0) ? 1 : -1);

				/* See if we went out the end of the object */
				if ((newSliceIndex < 0) || ((LbUsFourByte)newSliceIndex >= SubObjNumSlices)) {
					
					/* See if we are beyond end of limit cylinder */
					if ((newPosition->z_position <= limitZMin) ||
							(newPosition->z_position >= limitZMax)) {

						action = PhoTrkDiscard;
						break;
					}
					
					/* Check our acceptance angle */
					if (fabs(direction.cosine_z) > PHGGetSineOfAccAngle()) {
						
						
						action = PhoTrkDiscard;
						break;
					}
						
					/* Project to target cylinder if it is not equal to object cylinder */
					if (CylPosGetTargetRadius() <= CylPosGetObjectRadius()) {

						/* Set action to detect */
						action = PhoTrkDetect;
												
						/* Get out of routine */
						break;
					}
						/* Project to target cylinder */
					else if (!CylPosProjectToTargetCylinder(newPosition, &direction, &distToTargetSurface)) {
		
						action = PhoTrkDiscard;
						break;
					}
					else {
						
						/* Update distance tracked */
						distanceTracked += distToTargetSurface;

						action = PhoTrkDetect;
						break;
					}
		
				}
	
				/* If we didn't break, we entered a new slice of the object; recalculate slice variables */						
				PhoTrkEnterNewSlice(*newPosition, direction,
					trackingPhotonPtr->sliceIndex, distanceTracked,
					&distToNextX, &distToNextY, &distToNextZ,
			   		&generalDistToX, &generalDistToY, &generalDistToZ,
					&(trackingPhotonPtr->sliceIndex), &(trackingPhotonPtr->xIndex),
					&(trackingPhotonPtr->yIndex));
			}
			
		} while (true);
	} while (false);

	/* Save the distance traveled */
	*distTraveled = distanceTracked;
	
	/* Always clear this counter */
	phoTrkCellsInUse = 0;

	return (action);
}

/*********************************************************************************
*
*			Name:		PhoTrkProject
*
*			Summary:		Project a photon from a given position in a given
*							direction, a given distance to a new position.
*
*
*			Arguments:
*				PHG_Position	*startPos	- The starting position.
*				PHG_Direction	*direction	- The direction of travel.
*				double			distance	- The distance to travel.
*				PHG_Position	*newPos		- The new position.
*			Function return: None.
*
*********************************************************************************/
void PhoTrkProject(PHG_Position *startPos, PHG_Direction *direction, 
			double distance, PHG_Position *newPos)
{
	/* Compute the new position */
	newPos->x_position = startPos->x_position +
		(distance * direction->cosine_x);

	newPos->y_position = startPos->y_position +
		(distance * direction->cosine_y);

	newPos->z_position = startPos->z_position +
		(distance * direction->cosine_z);
	
}

/*********************************************************************************
*
*			Name:		phoTrkCalcScatterAngleSPECT
*
*			Summary:		Calculate a scatter angle that will assure the photon's
*							detection.
*
*
*			Arguments:
*				PHG_TrackingPhoton	*photonPtr	- The tracking photon.
*
*			Function return: True if valid angle computed.
*
*********************************************************************************/
Boolean phoTrkCalcScatterAngleSPECT(PHG_TrackingPhoton	*photonPtr)		
{
	Boolean angleFound = false;
	
	double	radialPos, muOut, deltaRphi, deltaMuOutMin, deltaMuOutMax, targetProb, deltaPhiMin, deltaPhiMax;
	double	xCosOut, yCosOut, zCosOut, absZSinOut, phiOut, cosScatAngle, eOut, weightOut, knValue, eIn;
	LbUsFourByte eIndex, muInIndex, zIndex, rIndex, deltaPhiIndex, deltaMuOutIndex, index;

	do {
	
		/*  If the photon will not enter within the acceptance angle limits break */
		if (phoTrkPositionIsAcceptable(photonPtr, &deltaMuOutMin, &deltaMuOutMax) == false) {
			break;
		}
	
		/* Compute indexes */
		eIndex = (LbUsFourByte) floor((
				(photonPtr->energy - PhgRunTimeParams.PhgMinimumEnergy) *
				CBFD_NUM_E_IN_BINS) / (PhgRunTimeParams.PhgNuclide.photonEnergy_KEV-PhgRunTimeParams.PhgMinimumEnergy));
				
		/* Check boundary */
		if (eIndex == CBFD_NUM_E_IN_BINS)
			eIndex = CBFD_NUM_E_IN_BINS -1;

		muInIndex = (LbUsFourByte) floor((
				(photonPtr->angle.cosine_z - -1.0) *
				CBFD_NUM_MU_IN_BINS) / (2.0));
				
		/* Check boundary */
		if (muInIndex == CBFD_NUM_MU_IN_BINS)
			muInIndex = CBFD_NUM_MU_IN_BINS -1;

		zIndex = (LbUsFourByte) floor((
				(photonPtr->location.z_position - CYLPOSGetCritZMin()) *
				CBFD_NUM_AXIAL_BINS) / (CYLPOSGetCritZMax()-CYLPOSGetCritZMin()));
				
		/* Check boundary */
		if (zIndex == CBFD_NUM_AXIAL_BINS)
			zIndex = CBFD_NUM_AXIAL_BINS -1;
			
		
		radialPos = PHGMATH_RadialPos(photonPtr->location.x_position, photonPtr->location.y_position);
		
		rIndex =  (LbUsFourByte) floor((
				(radialPos) *
				CBFD_NUM_RADIAL_BINS) / (CylPosGetObjectRadius()));
				
		/* Check boundary */
		if (rIndex == CBFD_NUM_RADIAL_BINS)
			rIndex = CBFD_NUM_RADIAL_BINS -1;
		
		
		/* Get target probability */
		targetProb = PhgMathGetRandomNumber() * (*phoTrkCBFDTotalProbAccept)[eIndex][muInIndex][zIndex][rIndex];
		
		/* Find largest deltaPhi, deltaMuOut such that targetProb < CumProbDeltaPhiMu() */
		index = phoTrkLookupProb((double *)(*phoTrkCBFDCumProbDeltaPhiMu)[eIndex][muInIndex][zIndex][rIndex],
					CBFD_NUM_DELTA_MU_OUT_BINS * CBFD_NUM_DELTA_PHI_BINS, targetProb);

		/* Convert index to deltaPhi and deltaMuOut */
		deltaPhiIndex = index / CBFD_NUM_DELTA_MU_OUT_BINS;
		deltaMuOutIndex = index % CBFD_NUM_DELTA_MU_OUT_BINS;

		deltaPhiMin = -PHGMATH_PI + (deltaPhiIndex *(PHGMATH_2PI/CBFD_NUM_DELTA_PHI_BINS));		
		deltaPhiMax = -PHGMATH_PI + ((deltaPhiIndex+1) *(PHGMATH_2PI/CBFD_NUM_DELTA_PHI_BINS));		
		deltaMuOutMin = (*phoTrkCBFDMinDeltaMu)[rIndex][zIndex] + deltaMuOutIndex * (((*phoTrkCBFDMaxDeltaMu)[rIndex][zIndex]
				- (*phoTrkCBFDMinDeltaMu)[rIndex][zIndex])/CBFD_NUM_DELTA_MU_OUT_BINS);
		deltaMuOutMax = (*phoTrkCBFDMinDeltaMu)[rIndex][zIndex] + (deltaMuOutIndex+1) * (((*phoTrkCBFDMaxDeltaMu)[rIndex][zIndex]
				- (*phoTrkCBFDMinDeltaMu)[rIndex][zIndex])/CBFD_NUM_DELTA_MU_OUT_BINS);
		
		phiOut = atan2(photonPtr->angle.cosine_y, photonPtr->angle.cosine_x)+
			(PhgMathGetRandomNumber()*(deltaPhiMax-deltaPhiMin) + deltaPhiMin);
		
		if (UNCCOLIsConeBeam() == false) {
			muOut = PhgMathGetRandomNumber() * (deltaMuOutMax - deltaMuOutMin) + deltaMuOutMin;
		}
		else {
			
			deltaRphi = phiOut - atan2(photonPtr->location.y_position, photonPtr->location.x_position);
			
			muOut = PHGMATH_Sine(phoTrkCalcCBFDOmega(deltaRphi, phoTrkRadiusOfFocalCircle,
								 radialPos, photonPtr->location.z_position)) +		
						PhgMathGetRandomNumber() * (deltaMuOutMax - deltaMuOutMin) + deltaMuOutMin; 

			
		}
		
		zCosOut = muOut;
		
		absZSinOut = PHGMATH_SquareRoot(1.0 - PHGMATH_Square(zCosOut));
		
		xCosOut = PHGMATH_Cosine(phiOut) * absZSinOut;
		yCosOut = PHGMATH_Sine(phiOut) * absZSinOut;
		cosScatAngle = (xCosOut * photonPtr->angle.cosine_x) + (yCosOut * photonPtr->angle.cosine_y) +
			(zCosOut * photonPtr->angle.cosine_z);
		
		eIn = photonPtr->energy/511.0;
		eOut = eIn/(1.0+eIn*(1.0-cosScatAngle));   /* !!!RH */
		if (eOut * 511.0 < PhgRunTimeParams.PhgMinimumEnergy) {
			break;
		}
		else {
		
			knValue = 0.5 * (PHGMATH_Square(eOut/eIn) *
				(eOut/eIn + eIn/eOut - 1 + PHGMATH_Square(cosScatAngle)));

			weightOut = photonPtr->photon_scatter_weight *
				(knValue * (*phoTrkCBFDTotalProbAccept)[eIndex][muInIndex][zIndex][rIndex])/
				((*phoTrkCBFDProbDeltaPhiMu)[eIndex][muInIndex][zIndex][rIndex][deltaPhiIndex][deltaMuOutIndex] *
					phoTrkTotalKN[(LbUsFourByte) photonPtr->energy]);
			
			photonPtr->angle.cosine_x = xCosOut;
			photonPtr->angle.cosine_y = yCosOut;
			photonPtr->angle.cosine_z = zCosOut;
			photonPtr->photon_scatter_weight = weightOut;		
			photonPtr->energy = eOut * 511.0;
		}
		
		angleFound = true;
	} while (false);
	
	return (angleFound);
}
/*********************************************************************************
*
*			Name:		phoTrkCalcScatterAnglePET
*
*			Summary:		Calculate a scatter angle that will assure the photon's
*							detection.
*							This is the pre-conebeam collimator modifications
*
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*
*			Function return: True if valid angle computed.
*
*********************************************************************************/
Boolean phoTrkCalcScatterAnglePET(PHG_TrackingPhoton	*trackingPhotonPtr)
{
	Boolean			angleFound;			/* Did we find an acceptable angle? */
	double			deltaAzimuth;		/* Change in azimuth */
	double			incomingAzimuth;	/* The incoming azimuth angle */
	double			incomingEnergy;		/* The incoming energy */
	double			inclination;		/* The inclination */
	double			minAccSine;			/* Minimum acceptable angle (to be detected) */
	double			maxAccSine;			/* Maximum accetpable angle (to be detected) */
	double			normalizedEnergy;	/* Energy value normalized by 511 */
	double			outgoingAzimuth;	/* The outgoing azimuth angle */
	double			targetDensityRange;	/* Range of possible target densitys */
	double			targetDensity;		/* Our target density */
	double			zCosOutRange;		/* Range of valid z cosine out ranges */
	double			cosScatterAngle;	/* Cosine of the scatter angle */
	double			knValue;			/* Klien Nishina probability */
	double			forcedProb;			/* Probability with which the photon was forced into the detection angle */
	double			fdFactor;			/* Product of phi_out & w_out bin size */
	double			totKNCrossSection;	/* Klien Nishina cross section for given energy */
	LbUsFourByte	inEnIndex;			/* Incoming energy FD Table index */
	LbUsFourByte	inZCosIndex;		/* Incoming z cosine FD Table index */
	LbUsFourByte	minAccAngSineIndex;	/* Minimum acceptance angle sine FD Table index */
	LbUsFourByte	maxAccAngSineIndex;	/* Maximum acceptance angle sine FD Table index */
	LbUsFourByte	deltaAzimuthIndex;	/* Change in azimuth FD Table index */
	LbUsFourByte	zCosOutIndex;		/* Outgoing z cosine FD Table index */
	PHG_Direction	newDir;				/* Our new direction */
	
	angleFound = false;
	
	do { /* Process Loop */
		
		/* See if the photon will enter within the acceptance angle limits */
		if ((phoTrkCalcAcceptanceRange(trackingPhotonPtr,
				&minAccSine, &maxAccSine)) == false) {
			
			break;
		}
				
		#ifdef PHG_DEBUG
			if (trackingPhotonPtr->energy < phoTrkFDTInfo.min_iei) {
			
				sprintf(phoTrkErrString, "Incoming energy below FD Table minimum energy (phoTrkCalcScatterAnglePET).\n"
					"FD Table minimum = %2.4f\t incoming energy = %2.4f\n", phoTrkFDTInfo.min_iei,
					trackingPhotonPtr->energy);
					
				PhgAbort(phoTrkErrString, true);
			}
		#endif
		
		/* Compute incoming energy in units of electron mass */
		incomingEnergy = trackingPhotonPtr->energy/511.0;

		/* Compute the incoming energy index (into the FD table) */
		inEnIndex = (LbUsFourByte) ((trackingPhotonPtr->energy -
			phoTrkFDTInfo.min_iei)/phoTrkFDTInfo.delta_iei);
		
		/* Check boundary (Necessary if energy index == max value in table) */
		if (inEnIndex >= phoTrkFDTInfo.num_iei)
			inEnIndex = phoTrkFDTInfo.num_iei - 1;
			
		/* Compute incoming z cosine index (into FD table) */
		inZCosIndex = (LbUsFourByte) ((trackingPhotonPtr->angle.cosine_z+1)/
		   (phoTrkFDTInfo.delta_iwi));
	
		/* Check boundary (Necessary if index == max value in table) */
		if (inZCosIndex >= phoTrkFDTInfo.num_iwi)
			inZCosIndex = phoTrkFDTInfo.num_iwi - 1;
			
		/* Compute minimum acceptance angle range index. */
		minAccAngSineIndex = (minAccSine+phoTrkFDTInfo.range_normalizer)/
			(phoTrkFDTInfo.delta_iwo);
		
		/* Check boundary */
		if (minAccAngSineIndex >= phoTrkFDTInfo.num_iwo)
			minAccAngSineIndex = phoTrkFDTInfo.num_iwo - 1;
			
		/* Compute maximum acceptance angle range index. */
		maxAccAngSineIndex = (maxAccSine+phoTrkFDTInfo.range_normalizer)/
			(phoTrkFDTInfo.delta_iwo);
		
		/* Check boundary */
		if (maxAccAngSineIndex >= phoTrkFDTInfo.num_iwo)
			maxAccAngSineIndex = phoTrkFDTInfo.num_iwo - 1;
					
		/* Compute target density range */
		targetDensityRange = 
			phoTrkFDTable[inEnIndex][inZCosIndex].iwoTable[maxAccAngSineIndex]
		   -
			phoTrkFDTable[inEnIndex][inZCosIndex].ipoCumTable[minAccAngSineIndex * phoTrkFDTInfo.num_ipo];
			
		/* Pick a target density */
		targetDensity = 
			phoTrkFDTable[inEnIndex][inZCosIndex].ipoCumTable[minAccAngSineIndex * phoTrkFDTInfo.num_ipo] +
			(PhgMathGetRandomNumber() * targetDensityRange);
			
		/* Lookup target density to get azimuth delta index */
		deltaAzimuthIndex = phoTrkLookupProb(
			phoTrkFDTable[inEnIndex][inZCosIndex].ipoCumTable +
			(minAccAngSineIndex * phoTrkFDTInfo.num_ipo),
			((maxAccAngSineIndex - minAccAngSineIndex)+1) * phoTrkFDTInfo.num_ipo,
			targetDensity);
		
		/* Compute zCosOutIndex knowing deltaAzimuthIndex */
		zCosOutIndex = minAccAngSineIndex + (deltaAzimuthIndex/phoTrkFDTInfo.num_ipo);

		/* Compute zCosOut range */
		zCosOutRange = phoTrkFDTInfo.delta_iwo;
		
		/* Pick an outgoing z cosine */
		newDir.cosine_z = (phoTrkFDTInfo.min_iwo + 
			(phoTrkFDTInfo.delta_iwo * zCosOutIndex)) +
			(zCosOutRange * PhgMathGetRandomNumber());
		
		
		/* Convert deltaAzimuthIndex from relative to global coordinates */
		deltaAzimuthIndex = (minAccAngSineIndex * phoTrkFDTInfo.num_ipo) + 
			deltaAzimuthIndex;
		
		/* Pick azimuth out from range of azimuth change bin */
		deltaAzimuth = 	
			(phoTrkFDTInfo.min_ipo + 
			(phoTrkFDTInfo.delta_ipo * (deltaAzimuthIndex % phoTrkFDTInfo.num_ipo)))
			+
			(PhgMathGetRandomNumber() * phoTrkFDTInfo.delta_ipo);
		
		/* Compute incoming azimuth and inclination */
		phoTrkCalcDirToAzimuth(&(trackingPhotonPtr->angle), &inclination,
			&incomingAzimuth);
			
		/* Compute azimuth out from azimuth in and azimuth change */
		outgoingAzimuth = incomingAzimuth + deltaAzimuth;
		
		/* Convert cosines and azimuths to direction cosine */
		phoTrkCalcAzimuthToDir(newDir.cosine_z, outgoingAzimuth, &newDir);
		
		/* Compute the cosine of the scatter angle for KN */
		cosScatterAngle = (trackingPhotonPtr->angle.cosine_x *
	        newDir.cosine_x) + (trackingPhotonPtr->angle.cosine_y *
			newDir.cosine_y) + (trackingPhotonPtr->angle.cosine_z *
            newDir.cosine_z);

		/* Save the photon's new direction */
		trackingPhotonPtr->angle = newDir;
		
		/* Compute the outgoing energy */
		{
			/* Normalize energy */
			normalizedEnergy = trackingPhotonPtr->energy/511.0;
			
			trackingPhotonPtr->energy =
				(normalizedEnergy /
					(1 + 
						(normalizedEnergy * 
							(1 - cosScatterAngle)
						)
					)
				)
				* 511.0;
		}
		
		/* See if the new energy is below our threshold */
		if (trackingPhotonPtr->energy < PhgRunTimeParams.PhgMinimumEnergy)
			break;
			
		/* Compute Klein-Nishina for given photon */
		normalizedEnergy = trackingPhotonPtr->energy/511.0;

		knValue = 0.5 * 
		    (PHGMATH_Square(normalizedEnergy/incomingEnergy)) *
			(normalizedEnergy/incomingEnergy +
			 	 incomingEnergy/normalizedEnergy - 1 +
			 	PHGMATH_Square(cosScatterAngle));

		/* Compute azimuth/inclination bin size */
		fdFactor = phoTrkFDTInfo.delta_ipo * phoTrkFDTInfo.delta_iwo;
		
		/* Compute the probability used to pick the chosen azimuth/inclination
			and divide by fdFactor 
		*/
		forcedProb = 
			(phoTrkFDTable[inEnIndex][inZCosIndex].ipoTable[deltaAzimuthIndex]/
			(targetDensityRange * fdFactor));

		
		/* Lookup Integral of KN formula at incoming energy */
		totKNCrossSection = phoTrkTotalKN[(LbUsFourByte) (incomingEnergy * 511.0)];
		
		/* Adjust photon weight according to probability of chosen scatter */
		trackingPhotonPtr->photon_scatter_weight = 
			trackingPhotonPtr->photon_scatter_weight *
			(knValue/(totKNCrossSection * forcedProb));
		
		angleFound = true;
	} while (false);
	
	return (angleFound);
}

/*********************************************************************************
*
*			Name:		phoTrkDoWeightWindow
*
*			Summary:	Analyze photon's weight, adjust if necessary,
*						and determine number of splits/roulettes to do.
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*				LbTwoByte			*numSplitsPtr		- Number of split attempts.
*				LbTwoByte			*numRoulettesPtr	- Number of roulette attempts.
*				Boolean				*discardItPtr		- Discard flag.
*
*			Function return: None.
*
*********************************************************************************/
void phoTrkDoWeightWindow(PHG_TrackingPhoton *trackingPhotonPtr,
		LbTwoByte *numSplitsPtr, LbTwoByte *numRoulettesPtr, Boolean *discardItPtr)
{
	double	estDetectedWeight;	/* The estimated detected weight */
	double	cellProductivity;	/* Productivity of given cell */
	double	randomNumber;		/* Random number used for russian roulette */
	
	/* Assume no attempts */
	*numSplitsPtr = -1;
	*numRoulettesPtr = -1;
	*discardItPtr = false;
	
	/* Compute the estimated detected weight */
	cellProductivity = PRODTBLGetScatCellProductivity(trackingPhotonPtr->sliceIndex,
		trackingPhotonPtr->angleIndex);
		
	estDetectedWeight = trackingPhotonPtr->photon_scatter_weight *
		cellProductivity;
	
	/* See if we are below threshold */
	if (estDetectedWeight < (PhgRunTimeParams.PhgMinWWRatio * trackingPhotonPtr->scatter_target_weight)) {
	
		/* This is Russian Roulette */
		*numRoulettesPtr = 0;
		
		/* Get random value */
		randomNumber = PhgMathGetRandomNumber();
		
		/* See if we should raise the weight */
		if (randomNumber < (estDetectedWeight/trackingPhotonPtr->scatter_target_weight)) {
			
			/* Adjust photon's weight up */
			trackingPhotonPtr->photon_scatter_weight = 
				trackingPhotonPtr->photon_scatter_weight *
				trackingPhotonPtr->scatter_target_weight/estDetectedWeight;
		}
		else {
			/* Russian Roulette says discard it */
			*discardItPtr = true;
			*numRoulettesPtr = 1;
		}
	}
	else if (estDetectedWeight > (PhgRunTimeParams.PhgMaxWWRatio * 
				trackingPhotonPtr->scatter_target_weight)) { /* See if we are above threshold */
	
		/* This is splitting */
		*numSplitsPtr = (LbUsFourByte) estDetectedWeight/
			trackingPhotonPtr->scatter_target_weight;
	}
}


/*********************************************************************************
*
*			Name:		phoTrkLookupProb
*
*			Summary:		Lookup a given probability value in the FD table.
*
*
*			Arguments:
*				double			*phiArray	- The data to search for the probability.
*				LbUsFourByte	numElements	- The number of elements in the array.
*				double			phi			- The value we are searching for.
*
*			Function return: Index within array containing phi.
*
*********************************************************************************/
LbUsFourByte phoTrkLookupProb(double *phiArray, LbUsFourByte numElements, double phi)
{
	LbUsFourByte	phiIndex;		/* The index we will find */
	LbUsFourByte	minIndex;		/* Lower bound of search segment */
	LbUsFourByte	maxIndex;		/* Upper bound of search segment */
	LbUsFourByte	testIndex;		/* The index we will test */

	/* Do simple sequential search if number of elements is small */
	if (numElements < 20) {
	
		/* Find first element strictly greater than */
		for (phiIndex = 0; phiIndex < numElements; phiIndex++) {
			 if (phiArray[phiIndex] > phi)
				break;
		}
		
		/* Walk back through indexes that are equal to this one to
			get to the last bin less than or equal to our value
		*/
		
		while((phiIndex > 0)
				&& (phiArray[phiIndex-1] == phiArray[phiIndex])) {
			phiIndex--;
		}
		
		/* Now select bin that is one less */
		if (phiIndex > 0)
			phiIndex--;
	}
	
	else {	/* Do binary search */
	
		/* Initialize loop variables */
		minIndex = 0;
		maxIndex = numElements;
		
		/* Find first element strictly greater than */
		while ((maxIndex - minIndex) != 1){
		
			/* Calculate index to test */
			testIndex = minIndex + (maxIndex - minIndex)/2;
			
			if (phi > phiArray[testIndex])
				minIndex = testIndex;
			else
				maxIndex = testIndex;
		}
	
		phiIndex = minIndex;
	}

	return (phiIndex);
}

/*********************************************************************************
*
*			Name:		PhoTrkInitialize
*
*			Summary:	Initialize the photon tracking module.
*
*			Arguments:
*
*			Function return: TRUE unless an error occurs.
*
*********************************************************************************/
Boolean PhoTrkInitialize()
{
	double			a;					/* See Haynor */
	double			a1;					/* See Haynor */
	double			a2;					/* See Haynor */
	double			a3;					/* See Haynor */
	LbUsFourByte	ieiIndex;			/* Index for iei tables */
	
	
	do { /* Process Loop */
	
		/* If we're not doing cone beam then just return */
		if (PHOTRKIsConeBeamForcedDetection() == true) {
			ErStGeneric("Forced Detection is not supported with SPECT Collimation right now, (PhoTrkInitialize)\n");
			break;
		}
	
		/* Clear cells in use variable */
		phoTrkCellsInUse = 0;
		
		/* Clear table pointers so termination will know to free their memory or not */
		phoTrkCBFDMaxDeltaMu = 0;
		phoTrkCBFDMinDeltaMu = 0;
		phoTrkCBFDProbDeltaPhiMu = 0;
		phoTrkCBFDCumProbDeltaPhiMu = 0;
		phoTrkCBFDTotalProbAccept = 0;
		
		/* Read in FD Table if requested */
		if (PHG_IsForcedDetection() == true) {
		
			/* Pre initialize variables in case of failure */
			phoTrkFDTable = 0;
		
			
			/* Call phoTrkInitializeFDTbl to create the table on the fly */
			if (phoTrkInitializeFDTbl(&phoTrkFDTInfo) == false){
				break;
			}
			
			/* If doing cone beam forced detection then initialize those tables */
			if (PHOTRKIsConeBeamForcedDetection()) {
				if (phoTrkInitializeCBFDTbl() == false)
					break;
			}
				

			/* Create table = Integral of KN formula by incoming energy
				(Zero entry can't be computed, loop runs 
				from 1 to 999
			*/
			phoTrkTotalKN[0] = 0;
			for (ieiIndex = 1; ieiIndex < 1000; ieiIndex++) {
				/* Initialize total KN table */
				a = (double) ieiIndex/511.0;
				a1 = a+1;
				a2 = 2.0 * a + 1;
				a3 = 3.0 * a + 1;
				phoTrkTotalKN[ieiIndex] = PHGMATH_2PI *
					((a1/pow(a,3.0)) * (2.*a*a1/a2 - log(a2)) + log(a2)/(2.*a) -
					a3/pow(a2,2.0));
			}
		}
		phoTrkIsInitialized = true;
	} while (false);
	
	return (phoTrkIsInitialized);
}

/*********************************************************************************
*
*			Name:		phoTrkInitializeFDTbl
*
*			Summary:	Initialize the forced detection data.
*
*			Arguments:
*				PhoTrkFDTInfoTy	*fdInfoPtr	- Our FD information.
*
*			Function return: TRUE unless an error occurs.
*
*********************************************************************************/
Boolean phoTrkInitializeFDTbl(PhoTrkFDTInfoTy	*fdInfoPtr)
{
#define			NUM_ZCOS_IN 	20
#define			NUM_E_IN 		10
#define			NUM_AZIMUTH_OUT 20
#define			NUM_ZCOS_OUT 	10

	Boolean			okay = false;		/* Process Flag */
	double			a;					/* See Haynor */
	double			a1;					/* See Haynor */
	double			a2;					/* See Haynor */
	double			a3;					/* See Haynor */
	LbUsFourByte	ieiIndex;			/* Index for iei tables */
	LbUsFourByte	iwiIndex;			/* Index for iwi tables */
	LbUsFourByte	curIpo;				/* Overall current ipo */

	Boolean			zero_cells = true;
	double			cosScatAngle;
	double			e_min;
	double			e_max;
	double			e_min511;
	double			energyIn;
	double			zCosIn;
	double			xCosIn;
	double			yCosIn;
	double			energyOut;
	double			zCosOut;
	double			xCosOut;
	double			yCosOut;
	double			azimuthOut;
	double			absZSinOut;
	double			knDens;
	double			deltaZCosIn;
	double			zCosInTbl[20];
	double			deltaEnergyIn;
	double			deltaAzimuthOut;
	double			deltaZCosOut;
	double			maxZCosOut;
	double			energyInTbl[NUM_E_IN];
	double			azimuthOutTbl[NUM_AZIMUTH_OUT];
	double			zCosOutTbl[NUM_ZCOS_OUT];
	LbUsFourByte	safety_factor = 5;
	LbUsFourByte	numEnergyIn;
	LbUsFourByte	numZCosOut;
	LbUsFourByte	numZCosIn;
	LbUsFourByte	numAzimuthOut;
	LbUsFourByte	tableIndex;
	LbUsFourByte	zCosInIndex;
	LbUsFourByte	energyInIndex;
	LbUsFourByte	azimuthOutIndex;
	LbUsFourByte	zCosOutIndex;

	typedef			double	probTbl1Ty[NUM_AZIMUTH_OUT][NUM_ZCOS_OUT][NUM_ZCOS_IN][NUM_E_IN];
	typedef			double	probTbl2Ty[NUM_ZCOS_OUT][NUM_ZCOS_IN][NUM_E_IN];

	probTbl1Ty		*knDensTbl=0;			/* Klein-Nishina densities */
	probTbl1Ty		*knDensCumTbl=0;		/* For each incoming energy and z cosine cumulates the
											   knDensTbl over the outgoing z cosine and azimuth change
											*/
	probTbl2Ty		*knDensCumTotalTbl=0;	/* For each incoming energy and z cosine cumulates the
											   knDensTbl, giving the total for all outgoing z cosines
											   less than or equal to the current z cosine bin
											*/
	
	do { /* Process Loop */
		
		/* Read in FD Table if requested */
		if (PHG_IsForcedDetection() == true) {

			/* Pre initialize variables in case of failure */
			phoTrkFDTable = 0;
						
			/* Initialize header parameters to default values */
			{
				/* Set num_ipo to max, 20, and compute remaining values */
				fdInfoPtr->num_ipo = 20;
				fdInfoPtr->delta_ipo = (PHGMATH_2PI)/fdInfoPtr->num_ipo;
				fdInfoPtr->min_ipo = (1 - 0.5) * fdInfoPtr->delta_ipo;

				/* Convert mid-point of bin to minimum value of bin */
				fdInfoPtr->min_ipo = fdInfoPtr->min_ipo - (fdInfoPtr->delta_ipo/2);
				
				/* Set num_iwo to max, 10, anc compute remaining values */
				fdInfoPtr->num_iwo = 10;
				fdInfoPtr->delta_iwo = 2.0 * PHGGetSineOfAccAngle()/fdInfoPtr->num_iwo;
				fdInfoPtr->delta_iwo = 2.0 * PHGGetSineOfAccAngle()/fdInfoPtr->num_iwo;
				fdInfoPtr->min_iwo = -PHGGetSineOfAccAngle() + (1 - 0.5) * fdInfoPtr->delta_iwo;
				fdInfoPtr->min_iwo = -PHGGetSineOfAccAngle() + (1 - 0.5) * fdInfoPtr->delta_iwo;
				
				/* Convert mid-point of bin to minimum value of bin */
				fdInfoPtr->min_iwo = fdInfoPtr->min_iwo - (fdInfoPtr->delta_iwo/2);

				/* Calculate maximum value */
				fdInfoPtr->max_iwo = fdInfoPtr->min_iwo + 
					(fdInfoPtr->delta_iwo * fdInfoPtr->num_iwo);
							
				/* Calculate range normalization */
				fdInfoPtr->range_normalizer = 
					(fdInfoPtr->max_iwo - fdInfoPtr->min_iwo)/2;
				
				/* Set num_iwi to max, 20, and compute remaining values */	
				fdInfoPtr->num_iwi = 20;
				fdInfoPtr->delta_iwi = 2.0/fdInfoPtr->num_iwi;
				fdInfoPtr->min_iwi = -1.0+(1 - 0.5) * fdInfoPtr->delta_iwi;

				/* Set num_iei to max, 10, and compute remaining values */
				fdInfoPtr->num_iei = 10;
				fdInfoPtr->min_iei = 0;
				fdInfoPtr->delta_iei = 
					(PhgRunTimeParams.PhgNuclide.photonEnergy_KEV - fdInfoPtr->min_iei)/
					fdInfoPtr->num_iei;
				fdInfoPtr->min_iei = fdInfoPtr->delta_iei;
				
				/* Verify at this point that the table's minimum energy is not above
					the user specified minimum energy threshold.
				*/
				if (phoTrkFDTInfo.min_iei > PhgRunTimeParams.PhgMinimumEnergy) {
					sprintf(phoTrkErrString, "\nUser supplied minimum energy (%3.3f) is less than"
						" Forced Detection minimum (%3.3f).\nPlease edit your parameter file"
						" and increase the minimum energy.\n",
						PhgRunTimeParams.PhgMinimumEnergy, phoTrkFDTInfo.min_iei);
					ErStGeneric(phoTrkErrString);
					break;
				}
			}
			
			/* Initialize local variables. This is just to conform to the method
				used in build_fdt. These variables should be replaced with a straight
				use of the fdInfoPtr fields AFTER those fields have been renamed to
				conform to the names of these locals 
			*/
			{
				
				numZCosIn = fdInfoPtr->num_iwi;
				deltaZCosIn = fdInfoPtr->delta_iwi;

				/* Initialize zCosInTbl table */
				for (tableIndex = 1; tableIndex <= numZCosIn; tableIndex++) {
					/* Notice that zCosInTbl = middle of i'th interval */
					zCosInTbl[tableIndex-1] = -1.0+(tableIndex - 0.5) * deltaZCosIn;
				}
				
				numEnergyIn = fdInfoPtr->num_iei;
				e_min = 0;
				e_max = PhgRunTimeParams.PhgNuclide.photonEnergy_KEV;
				
				/* Calculate deltaEnergyIn */
				deltaEnergyIn = (e_max - e_min)/numEnergyIn;
				
				/* Initialize energyInTbl array */
				for (tableIndex = 1; tableIndex <= numEnergyIn; tableIndex++) {
					energyInTbl[tableIndex-1] = e_min + (tableIndex * deltaEnergyIn);
				}
				/* Calculate e_min511 */
				e_min511 = (e_min - safety_factor)/511.0;
				
				numAzimuthOut = fdInfoPtr->num_ipo;
				
				/* Calculate deltaAzimuthOut */
				deltaAzimuthOut = fdInfoPtr->delta_ipo;
		
				/* Initialize azimuthOutTbl table */
				for (tableIndex = 1; tableIndex <= numAzimuthOut; tableIndex++) {
					/* Exiting phi's from 0 to 2 pi */
					azimuthOutTbl[tableIndex-1] = (tableIndex - 0.5) * deltaAzimuthOut;
				}
		
				numZCosOut = fdInfoPtr->num_iwo;
				maxZCosOut = PHGGetSineOfAccAngle();
				deltaZCosOut = fdInfoPtr->delta_iwo;
				
				/* Initialize zCosOutTbl table */
				for (tableIndex = 1; tableIndex <= numZCosOut; tableIndex++) {
					/* zCosOutTbl both pos and neg */
					zCosOutTbl[tableIndex-1] = -maxZCosOut + (tableIndex - 0.5) * deltaZCosOut;
				}
				
				
			}

			/* Allocate iei table */
			if ((phoTrkFDTable = (PhoTrkFDTIeiTy) LbMmAlloc(sizeof(PhoTrkFDTIeiTy) *
					fdInfoPtr->num_iei)) == 0) {
					
				break;
			}
			
			/* Initialize the table to zero */
			for (ieiIndex = 0; ieiIndex < fdInfoPtr->num_iei; ieiIndex++) {
				phoTrkFDTable[ieiIndex] = 0;
			}
			
			/* This loop allocates all of the sub-tables */
			for (ieiIndex = 0; ieiIndex < fdInfoPtr->num_iei; ieiIndex++) {
				
				/* Allocate an iwi table */
				if ((phoTrkFDTable[ieiIndex] = (PhoTrkFDTIwiTy *) LbMmAlloc(sizeof(PhoTrkFDTIwiTy) *
						fdInfoPtr->num_iwi)) == 0) {
					
					goto FAIL;
				}
				
				/* Pre initialize to zero */
				for (iwiIndex = 0; iwiIndex < fdInfoPtr->num_iwi; iwiIndex++) {
					phoTrkFDTable[ieiIndex][iwiIndex].iwoTable = 0;
					phoTrkFDTable[ieiIndex][iwiIndex].ipoTable = 0;
					phoTrkFDTable[ieiIndex][iwiIndex].ipoCumTable = 0;
				}

				/* For each iwi element */
				for (iwiIndex = 0; iwiIndex < fdInfoPtr->num_iwi; iwiIndex++) {
					
					/* Allocate an iwo table */
					if ((phoTrkFDTable[ieiIndex][iwiIndex].iwoTable = (double *)
							LbMmAlloc(sizeof(double) * fdInfoPtr->num_iwo)) == 0) {
							
						goto FAIL;
					}
								
					/* Allocate an ipo table */
					if ((phoTrkFDTable[ieiIndex][iwiIndex].ipoTable = (double *)
							LbMmAlloc(sizeof(double) * fdInfoPtr->num_iwo *
							fdInfoPtr->num_ipo)) == 0) {
							
						goto FAIL;
					}
						
					/* Allocate an ipo_cum table */
					if ((phoTrkFDTable[ieiIndex][iwiIndex].ipoCumTable = (double *)
							LbMmAlloc(sizeof(double) * fdInfoPtr->num_iwo *
							fdInfoPtr->num_ipo)) == 0) {
							
						goto FAIL;
					}
						
				}
					
			}
			
			/* Create table = Integral of KN formula by incoming energy
				(Zero entry can't be computed, loop runs 
				from 1 to 999
			*/
			phoTrkTotalKN[0] = 0;
			for (ieiIndex = 1; ieiIndex < 1000; ieiIndex++) {
				/* Initialize total KN table */
				a = (double) ieiIndex/511.0;
				a1 = a+1;
				a2 = 2.0 * a + 1;
				a3 = 3.0 * a + 1;
				phoTrkTotalKN[ieiIndex] = PHGMATH_2PI *
					((a1/pow(a,3.0)) * (2.*a*a1/a2 - log(a2)) + log(a2)/(2.*a) -
					a3/pow(a2,2.0));
			}

			/*** This next section is the code from build_fdt which creates the
				FD tables. The statements which used to write values to the
				FD file have been replaced with assignments into our dynamic
				FD tables.
			***/

			/* Allocate memory for knDensTbl table */
			if ((knDensTbl = (probTbl1Ty *) malloc(sizeof(probTbl1Ty))) == 0){
				ErStGeneric("Unable to allocate memory for probability table.");
				break;
			}
	
			/* Allocate memory for knDensCum table */
			if ((knDensCumTbl = (probTbl1Ty *) malloc(sizeof(probTbl1Ty))) == 0){
				ErStGeneric("Unable to allocate memory for probability table.");
				break;
			}
	
			/* Allocate memory for knDensCumTotal table */
			if ((knDensCumTotalTbl = (probTbl2Ty *) malloc(sizeof(probTbl2Ty))) == 0){
				ErStGeneric("Unable to allocate memory for probability table.");
				break;
			}
			
			/* Loop through all energies to create the tables */
			for (energyInIndex = 1; energyInIndex <= numEnergyIn; energyInIndex++){
				
				/* Normalize photon energy */
				energyIn = energyInTbl[energyInIndex-1]/511.0;
				
				for (zCosInIndex = 1; zCosInIndex <= numZCosIn; zCosInIndex++) {
				
					zCosIn = zCosInTbl[zCosInIndex-1];
					
					/* The tables rely on the change in azimuth, deltaAzimuth, rather than
					   on a fixed azimuthal value.  Thus we pick an arbitrary azimuth (i.e.
					   x and y cosine) for our incoming azimuth, with the condition that
					   xCosIn*xCosIn + yCosIn*yCosIn +  zCosIn*zCosIn = 1.  We chose
					   yCosIn = 0.0.
					*/
					{
						xCosIn = sqrt(1.0 - (zCosIn * zCosIn));
						yCosIn = 0.0;
					}
									
					knDens = 0.0;
					
					for (zCosOutIndex = 1; zCosOutIndex <= numZCosOut; zCosOutIndex++) {
						zCosOut = zCosOutTbl[zCosOutIndex-1];
						
						absZSinOut = sqrt(1.0 - (zCosOut * zCosOut));
						
						for ( azimuthOutIndex = 1; azimuthOutIndex <= numAzimuthOut; azimuthOutIndex++) {
							
							azimuthOut = azimuthOutTbl[azimuthOutIndex-1];
							
							xCosOut = cos(azimuthOut) * absZSinOut;
							
							yCosOut = sin(azimuthOut) * absZSinOut;
							
							/* Now have in and out angles; calculate cosine of scattering angle */
							cosScatAngle = (xCosOut * xCosIn) + (yCosOut * yCosIn) + (zCosOut * zCosIn);
							
							/* Conservation of momentum yields out energy */
							energyOut = energyIn/(1.0 + energyIn*(1.0 - cosScatAngle));
							
							if (azimuthOutIndex == 1) {
								if (zCosOutIndex == 1) {
									(*knDensCumTbl)[0][0][zCosInIndex-1][energyInIndex-1] = 0.0;
								}
								else {
									(*knDensCumTbl)[0][zCosOutIndex-1][zCosInIndex-1][energyInIndex-1] = 
										knDens + (*knDensCumTbl)[numAzimuthOut-1][zCosOutIndex-2][zCosInIndex-1][energyInIndex-1];
								}
							}
							else {
								(*knDensCumTbl)[azimuthOutIndex-1][zCosOutIndex-1][zCosInIndex-1][energyInIndex-1] =
									knDens + (*knDensCumTbl)[azimuthOutIndex-2][zCosOutIndex-1][zCosInIndex-1][energyInIndex-1];
							}
							
							/* Now calculate value of knDens for the next iteration */
							if ((energyOut < e_min511) && (zero_cells == true)) {
								/* Blank cells leading to too-low energy */
								knDens = 0.0;
							}
							else {
								knDens = 0.5 * ((energyOut/energyIn) * (energyOut/energyIn)) * 
									(energyOut/energyIn + energyIn/energyOut - 1.0 + 
									(cosScatAngle*cosScatAngle));
							}
							
							(*knDensTbl)[azimuthOutIndex-1][zCosOutIndex-1][zCosInIndex-1][energyInIndex-1] = knDens;
						}
						
						/* Each bin in knDensCumTotalTbl gives the total relative probability that
						   the outgoing z cosine will be <= the corresponding z cosine bin
						*/
						(*knDensCumTotalTbl)[zCosOutIndex-1][zCosInIndex-1][energyInIndex-1] = 
							(*knDensCumTbl)[numAzimuthOut-1][zCosOutIndex-1][zCosInIndex-1][energyInIndex-1] + knDens;
	
					} /* ENDFOR zCosOutIndex */
				
				
					/* Reset our curIpo counter */
					curIpo = 0;
					
					/* Write out values of knDensCumTotalTbl and values of knDensTbl for this value of zCosInTbl, energyInTbl */
					for (zCosOutIndex = 0; zCosOutIndex < numZCosOut; zCosOutIndex++) {
						
						/* Get our knDensCumTotalTbl value */
						phoTrkFDTable[energyInIndex-1][zCosInIndex-1].iwoTable[zCosOutIndex] = 
							(*knDensCumTotalTbl)[zCosOutIndex][zCosInIndex-1][energyInIndex-1];
						
						/* Loop through the outgoing azimuth values */
						for (azimuthOutIndex = 0; azimuthOutIndex < numAzimuthOut; azimuthOutIndex++) {
						
							/* Get our knDensTbl value */
							phoTrkFDTable[energyInIndex-1][zCosInIndex-1].ipoTable[curIpo] =
								(*knDensTbl)[azimuthOutIndex][zCosOutIndex][zCosInIndex-1][energyInIndex-1];
							
							/* Get our knDensCumTbl value */
							phoTrkFDTable[energyInIndex-1][zCosInIndex-1].ipoCumTable[curIpo] = 
								(*knDensCumTbl)[azimuthOutIndex][zCosOutIndex][zCosInIndex-1][energyInIndex-1];
							
							/* Increment ipo counter */
							curIpo++;
						}
					}
				} /* ENDFOR zCosInIndex */
						
			} /* ENDFOR energyInIndex */
		
		if (knDensTbl != 0)
			free(knDensTbl);
			
		if (knDensCumTbl != 0)
			free(knDensCumTbl);
			
		if (knDensCumTotalTbl != 0)
			free(knDensCumTotalTbl);

		}
		okay = true;
		FAIL:;
	} while (false);
	
	/* Free memory if error occured */
	if (!okay) {
		if (phoTrkFDTable != 0) {
			/* Free the memory used by the FDT tables */
			for (ieiIndex = 0; ieiIndex < phoTrkFDTInfo.num_iei; ieiIndex++) {
				
				if (phoTrkFDTable[ieiIndex] != 0) {
					/* For each iwi element */
					for (iwiIndex = 0; iwiIndex < phoTrkFDTInfo.num_iwi; iwiIndex++) {
						
						/* Free an iwo table */
						if (phoTrkFDTable[ieiIndex][iwiIndex].iwoTable != 0)
							LbMmFree((void **) &(phoTrkFDTable[ieiIndex][iwiIndex].iwoTable));
									
						/* Free an ipo table */
						if (phoTrkFDTable[ieiIndex][iwiIndex].ipoTable != 0)
							LbMmFree((void **) &(phoTrkFDTable[ieiIndex][iwiIndex].ipoTable));
							
						/* Free an ipo_cum table */
						if (phoTrkFDTable[ieiIndex][iwiIndex].ipoCumTable != 0)
							LbMmFree((void **) &(phoTrkFDTable[ieiIndex][iwiIndex].ipoCumTable));
							
					}
					LbMmFree((void **) &(phoTrkFDTable[ieiIndex]));
				}
			}
			
			free(phoTrkFDTable);
		}
			
		if (knDensTbl != 0)
			free(knDensTbl);
			
		if (knDensCumTbl != 0)
			free(knDensCumTbl);
			
		if (knDensCumTotalTbl != 0)
			free(knDensCumTotalTbl);
	}
	return (okay);
}

/*********************************************************************************
*
*			Name:		phoTrkInitializeCBFDTbl
*
*			Summary:	Initialize the cone beam forced detection tables.
*
*			Arguments:
*
*			Function return: TRUE unless an error occurs.
*
*********************************************************************************/
Boolean phoTrkInitializeCBFDTbl()		
{
	Boolean			okay = false;		/* Process Flag */

	double	omega;
	double	sineOmega;

	double	radialPos;
	double	radiusOfTarget;
	double	targetMinZ;
	double	targetMaxZ;
	double	thetaAccept;
	double	thetaMaxDev;				/* Maximum allowable deviation from the central cone-beam angle */
	double	deviationSafetyFactor;		/* Make sure we don't throw away any good angles */

	LbUsFourByte	rIndex, zIndex;
	double			zMin, zMax, rMin, rMax, tanMax, tanMin;
	double			radialBinSize, axialBinSize;
	double			omegaMin;
	Boolean			targetIntersection;
	double	muOut;
	double	e_min511, eOut;
	double	safety_factor = 0.7;
	double	eInCent, muInCent, zCent, rCent, deltaMuOutCent, deltaPhiCent, deltaPhiMin, deltaPhiMax;
	double	energyBinSize, deltaMuBinSize, deltaPhiBinSize, muInBinSize;
	double	deltaMuOutMin, deltaMuOutMax;
	double absZSinOut, xrCosOut, cosScatAngle;
	LbUsFourByte eIndex, muInIndex, deltaPhiIndex, deltaMuOutIndex;
				
		
		do { /* Process Loop */
			
			/* Clear for error checking purposes */
			phoTrkCBFDMaxDeltaMu = 0;
			phoTrkCBFDMinDeltaMu = 0;
			phoTrkCBFDProbDeltaPhiMu = 0;
			phoTrkCBFDTotalProbAccept = 0;
			phoTrkCBFDCumProbDeltaPhiMu = 0;
			
			
			/* Allocate table */
			if ((phoTrkCBFDMaxDeltaMu = (PhoTrkCBFDDeltaMuTy*)LbMmAlloc(sizeof(PhoTrkCBFDDeltaMuTy))) == 0) {
				break;
			}
			
			/* Allocate table */
			if ((phoTrkCBFDMinDeltaMu = (PhoTrkCBFDDeltaMuTy*)LbMmAlloc(sizeof(PhoTrkCBFDDeltaMuTy))) == 0) {
				break;
			}
			
			/* Allocate table */
			if ((phoTrkCBFDProbDeltaPhiMu = (PhoTrkCBFDProbTy*)LbMmAlloc(sizeof(PhoTrkCBFDProbTy))) == 0) {
				break;
			}
			
			/* Allocate table */
			if ((phoTrkCBFDTotalProbAccept = (PhoTrkCBFDTotProbTy*)LbMmAlloc(sizeof(PhoTrkCBFDTotProbTy))) == 0) {
				break;
			}
			
			/* Allocate table */
			if ((phoTrkCBFDCumProbDeltaPhiMu = (PhoTrkCBFDProbTy*)LbMmAlloc(sizeof(PhoTrkCBFDProbTy))) == 0) {
				break;
			}


	/* Initialize Delta Mu Tables */
	{
		
		/* Get some limits */
		CylPosGetTargetZRange(&targetMinZ, &targetMaxZ);
		radiusOfTarget = CylPosGetTargetRadius();
		
		phoTrkRadiusOfFocalCircle = UNCCOLGetFocalLength() - ColRunTimeParams[ColCurParams].UNCSPECTCol.RadiusOfRotation;	
		phoTrkFocalLength = UNCCOLGetFocalLength();	
		phoTrkCollimatorZMin = UNCCOLGetCollZMin();	
		phoTrkCollimatorZMax = UNCCOLGetCollZMax();	 
		
		
		thetaAccept = PHGMATH_RadiansFromDegrees(PHGGetAccAngle());
		
		
		if (UNCCOLIsConeBeam() == false){			
			thetaMaxDev = thetaAccept;				
		} else {									
			deviationSafetyFactor = 1.05;	/* To make sure we don't throw away good angles */			
			thetaMaxDev = deviationSafetyFactor *
					atan( 2.0*ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleRadius / ColRunTimeParams[ColCurParams].UNCSPECTCol.Thickness );		
		}											
		
		/* Compute the bin sizes */
		radialBinSize = CylPosGetObjectRadius()/CBFD_NUM_RADIAL_BINS;
		axialBinSize = (CYLPOSGetCritZMax()-CYLPOSGetCritZMin())/CBFD_NUM_AXIAL_BINS;
		
		/* Clear loop incrementers */
		/* the following two lines are unnecessary--these values are set in every loop:
		rMin = 0.0;
		zMin = CYLPOSGetCritZMin();	*/

		for (rIndex = 0; rIndex < CBFD_NUM_RADIAL_BINS; rIndex++){
		
			/* Get min/max values of current bins */
			rMin =  (rIndex*radialBinSize);
			rMax = rMin + radialBinSize;
			
			radialPos = (rMax-rMin)/2.0;
			
			for (zIndex = 0; zIndex < CBFD_NUM_AXIAL_BINS; zIndex++){
			
				/* Get min/max values of current bins */
				zMin =  CYLPOSGetCritZMin() +(zIndex* axialBinSize);
				zMax = zMin + axialBinSize;
				
				if (UNCCOLIsConeBeam() == false){
					
					omegaMin = 0.0;
					
					tanMax = tan(thetaAccept);
					tanMin = -tanMax;
					if ((zMin > 0.0) && ((zMin + ((radiusOfTarget-rMax)*tanMin) > targetMaxZ))){	
						targetIntersection = false;
					}
					else if ((zMax < 0.0) && ((zMax + ((radiusOfTarget-rMax)*tanMax) < targetMinZ))) {	
						targetIntersection = false;
					}
					else {
						targetIntersection = true;
					}
				}
				else if (zMin > 0.0) {
					
					omegaMin = atan(zMin/fabs(phoTrkRadiusOfFocalCircle+rMax));
					
					if ((zMin + ((radiusOfTarget-rMax)*tan(omegaMin - thetaMaxDev))) > targetMaxZ) {	
						targetIntersection = false;
					}
					else {
						targetIntersection = true;
					}
				}
				else if (zMax < 0.0) {
				
					omegaMin = atan(zMax/fabs(phoTrkRadiusOfFocalCircle+rMax));
					
					
					if ((zMax + ((radiusOfTarget-rMax)*tan(omegaMin + thetaMaxDev))) < targetMinZ) {	 
						targetIntersection = false;
					}
					else {
						targetIntersection = true;
					}
				}
				else {
					omegaMin = 0.0;
					targetIntersection = true;
				}
				
				if (targetIntersection == false) {
					(*phoTrkCBFDMaxDeltaMu)[rIndex][zIndex] = 0.0;
					(*phoTrkCBFDMinDeltaMu)[rIndex][zIndex] = 0.0;
				}
				else if (omegaMin >= (thetaMaxDev/2.0)) {
					(*phoTrkCBFDMaxDeltaMu)[rIndex][zIndex] = PHGMATH_Sine(omegaMin+thetaMaxDev) - PHGMATH_Sine(omegaMin);
					(*phoTrkCBFDMinDeltaMu)[rIndex][zIndex] = PHGMATH_Sine(omegaMin-thetaMaxDev) - PHGMATH_Sine(omegaMin);
				}
				else if (omegaMin <= -(thetaMaxDev/2.0)) {
					(*phoTrkCBFDMaxDeltaMu)[rIndex][zIndex] = PHGMATH_Sine(omegaMin+thetaMaxDev) - PHGMATH_Sine(omegaMin);
					(*phoTrkCBFDMinDeltaMu)[rIndex][zIndex] = PHGMATH_Sine(omegaMin-thetaMaxDev) - PHGMATH_Sine(omegaMin);
				}
				else  {
					(*phoTrkCBFDMaxDeltaMu)[rIndex][zIndex] = PHGMATH_Sine(thetaMaxDev);
					(*phoTrkCBFDMinDeltaMu)[rIndex][zIndex] = -PHGMATH_Sine(thetaMaxDev);
				}
					
			}
		}
	}
	/* Initialize remaining tables */
	{

		/* Calculate e_min511 */
		e_min511 = safety_factor * (PhgRunTimeParams.PhgMinimumEnergy/511.0);	 
		
		/* Get some limits */
		CylPosGetTargetZRange(&targetMinZ, &targetMaxZ);
		radiusOfTarget = CylPosGetTargetRadius();
		
		/* Compute the bin sizes */
		radialBinSize = CylPosGetObjectRadius()/CBFD_NUM_RADIAL_BINS;
		axialBinSize = (CYLPOSGetCritZMax()-CYLPOSGetCritZMin())/CBFD_NUM_AXIAL_BINS;
		energyBinSize = (PhgRunTimeParams.PhgNuclide.photonEnergy_KEV -
			PhgRunTimeParams.PhgMinimumEnergy)/CBFD_NUM_E_IN_BINS;
		deltaPhiBinSize = ((PHGMATH_2PI))/CBFD_NUM_DELTA_PHI_BINS;
		muInBinSize =  (2.0)/CBFD_NUM_MU_IN_BINS;
		
		eInCent = (PhgRunTimeParams.PhgMinimumEnergy + (energyBinSize/2))/511.0;
		energyBinSize /= 511.0;
		
		
		for (eIndex = 0; eIndex < CBFD_NUM_E_IN_BINS; eIndex++) {
			muInCent = -1.0 + (muInBinSize/2.0);
		for (muInIndex = 0; muInIndex < CBFD_NUM_MU_IN_BINS; muInIndex++){
			zCent = CYLPOSGetCritZMin() + (axialBinSize/2.0);
		for (zIndex = 0; zIndex < CBFD_NUM_AXIAL_BINS; zIndex++) {
			rCent = radialBinSize/2.0;
		for (rIndex = 0; rIndex < CBFD_NUM_RADIAL_BINS; rIndex++){
			
			deltaPhiMin = -PHGMATH_PI;
			(*phoTrkCBFDTotalProbAccept)[eIndex][muInIndex][zIndex][rIndex] = 0.0;
			deltaPhiMax = deltaPhiMin + deltaPhiBinSize;
			deltaPhiCent = deltaPhiMin + ((deltaPhiMax-deltaPhiMin)/2.0);

			deltaMuBinSize = ((*phoTrkCBFDMaxDeltaMu)[rIndex][zIndex]-(*phoTrkCBFDMinDeltaMu)[rIndex][zIndex]) /
								CBFD_NUM_DELTA_MU_OUT_BINS;    

			omega = phoTrkCalcCBFDOmega(deltaPhiCent, phoTrkRadiusOfFocalCircle, rCent, zCent); 
			sineOmega = PHGMATH_Sine(omega);
			for (deltaPhiIndex = 0; deltaPhiIndex < CBFD_NUM_DELTA_PHI_BINS; deltaPhiIndex++){
				
				deltaMuOutMin = (*phoTrkCBFDMinDeltaMu)[rIndex][zIndex];    
				deltaMuOutMax = deltaMuOutMin + deltaMuBinSize;   
				deltaMuOutCent = deltaMuOutMin + ((deltaMuOutMax - deltaMuOutMin)/2.0);   

				if (UNCCOLIsConeBeam() == false) {
					muOut = deltaMuOutCent;
				}
				else {
						muOut = sineOmega + deltaMuOutCent;    
				}
				
				for (deltaMuOutIndex = 0; deltaMuOutIndex < CBFD_NUM_DELTA_MU_OUT_BINS; deltaMuOutIndex++) {
				
					
					absZSinOut = PHGMATH_SquareRoot(1.0 - PHGMATH_Square(muOut));
					xrCosOut = PHGMATH_Cosine(deltaPhiCent) * absZSinOut;
					cosScatAngle = (xrCosOut * PHGMATH_SquareRoot(1 - PHGMATH_Square(muInCent))) +
						(muOut * muInCent);
					
					eOut = eInCent/(1.0+eInCent *(1.0 - cosScatAngle));
					
					if ((eOut < e_min511) || ((*phoTrkCBFDMaxDeltaMu)[rIndex][zIndex]
												<= (*phoTrkCBFDMinDeltaMu)[rIndex][zIndex])) {   
					
						(*phoTrkCBFDProbDeltaPhiMu)[eIndex][muInIndex][zIndex][rIndex][deltaPhiIndex][deltaMuOutIndex] 
							= 0.0;
					}
					else {
						(*phoTrkCBFDProbDeltaPhiMu)[eIndex][muInIndex][zIndex][rIndex][deltaPhiIndex][deltaMuOutIndex] 
							= 0.5 * ((eOut/eInCent) * (eOut/eInCent)) *
									(eOut/eInCent + eInCent/eOut - 1.0 +
										(cosScatAngle*cosScatAngle));
					}
					
					(*phoTrkCBFDTotalProbAccept)[eIndex][muInIndex][zIndex][rIndex] +=
						(*phoTrkCBFDProbDeltaPhiMu)[eIndex][muInIndex][zIndex][rIndex][deltaPhiIndex][deltaMuOutIndex] 
						*
						(deltaPhiMax - deltaPhiMin) *
						deltaMuBinSize;   
						
					(*phoTrkCBFDCumProbDeltaPhiMu)[eIndex][muInIndex][zIndex][rIndex][deltaPhiIndex][deltaMuOutIndex] =
						(*phoTrkCBFDTotalProbAccept)[eIndex][muInIndex][zIndex][rIndex];    

					deltaMuOutMin += deltaMuBinSize;   
					deltaMuOutMax = deltaMuOutMin + deltaMuBinSize;    
					deltaMuOutCent += deltaMuBinSize;   
						
				}
				
				/* Increment bin centers */
				deltaPhiMin += deltaPhiBinSize;   
				deltaPhiMax = deltaPhiMin + deltaPhiBinSize;
				deltaPhiCent += deltaPhiBinSize;
				
			}
				
			/* Increment bin centers */
			rCent = rCent + radialBinSize;
		
		}
			/* Increment bin centers */
			zCent = zCent + axialBinSize;
		}
			/* Increment bin centers */
			muInCent = muInCent + muInBinSize;
		}
			/* Increment bin centers */
			eInCent = eInCent +  energyBinSize;
		}

	}
			okay = true;
			FAIL:;
		} while (false);
		
		/* Free memory if error occured */
		if (!okay) {
			
			/* Free memory that was allocated */
			if (phoTrkCBFDMaxDeltaMu != 0) {
				LbMmFree((void **)&phoTrkCBFDMaxDeltaMu);
			}
			
			/* Free memory that was allocated */
			if (phoTrkCBFDMinDeltaMu != 0) {
				LbMmFree((void **)&phoTrkCBFDMinDeltaMu);
			}
			
			/* Free memory that was allocated */
			if (phoTrkCBFDProbDeltaPhiMu != 0) {
				LbMmFree((void **)&phoTrkCBFDProbDeltaPhiMu);
			}
			
			/* Free memory that was allocated */
			if (phoTrkCBFDTotalProbAccept != 0) {
				LbMmFree((void **)&phoTrkCBFDTotalProbAccept);
			}
			
			/* Free memory that was allocated */
			if (phoTrkCBFDCumProbDeltaPhiMu != 0) {
				LbMmFree((void **)&phoTrkCBFDCumProbDeltaPhiMu);
			}
		}
				
		return (okay);
}

/*********************************************************************************
*
*			Name:		phoTrkCalcCBFDOmega
*
*			Summary:	Calculate 'omega' for cone beam forced detection.
*
*			Arguments:
*				double	deltaRphi			- Delta R Phi
*				double	phoTrkRadiusOfFocalCircle	
*				double	radialPos			- radial position
*				double	zPos				- axial position
*				double	phoTrkRadiusOfFocalCircle - self explanatary
*			Function return: double - omega.
*
*********************************************************************************/
double phoTrkCalcCBFDOmega(double deltaRphi, double phoTrkRadiusOfFocalCircle, double radialPos,
		double zPos)		
{
	double 			omega;
	LbUsFourByte	numRoots;
	double			minRoot, maxRoot;
	double			a,b,c;
	double			Xf, Yf;
	
	#ifdef PHG_DEBUG
			double solution1, solution2;
	#endif
	
	/* Solve for Xf */
	if (PHGMATH_Cosine(deltaRphi) < -0.1) {		
		
		a = 1/(PHGMATH_Square(PHGMATH_Cosine(deltaRphi)));
		b = (-2.0*radialPos)*((1/PHGMATH_Square(PHGMATH_Cosine(deltaRphi))-1.0));
		c = (PHGMATH_Square(radialPos)*((1/PHGMATH_Square(PHGMATH_Cosine(deltaRphi))-1.0)))-PHGMATH_Square(phoTrkRadiusOfFocalCircle);
		
		/* Note, unless in debug mode, we know that there will be one real root solution here */
		numRoots = PhgMathSolveQuadratic(a, b, c, &minRoot, &maxRoot);

		#ifdef PHG_DEBUG
		{
			
			if (numRoots != 2) {
				PhgAbort("Did not get 2 roots from solution to quadratic (phoTrkCalcCBFDOmega)", false);
			}
			
			/* Verify the roots solve the quadratic equation */
			solution1 = ((a*PHGMATH_Square(minRoot)) + (b*minRoot) + c);
			
			if (PhgMathRealNumAreEqual(solution1, 0.0, 
		        -4, 0, 0, 0) == false) {
				PhgAbort("Root from quadratic do not solve properly (minRoot)", false);
			}
			solution2 = ((a*PHGMATH_Square(maxRoot)) + (b*maxRoot) + c);
			
			if (PhgMathRealNumAreEqual(solution2, 0.0, 
		        -4, 0, 0, 0) == false) {
				PhgAbort("Root from quadratic do not solve properly (maxRoot)", false);
			}
		}
		#endif
		
		Xf = maxRoot;
	
		#ifdef PHG_DEBUG
		if ((maxRoot < radialPos) || (minRoot > radialPos)) {
			PhgAbort("Roots from quadratic do not staddle radialPos in cos(deltaRphi) < 0.1 (phoTrkCalcCBFDOmega)", false);
		}
		#endif
		
		/* must use equation 6a here => cos not sin:  omega = atan(zPos*fabs(PHGMATH_Sine(deltaRphi)/(Xf-radialPos))); */
		omega = atan(zPos*fabs(PHGMATH_Cosine(deltaRphi)/(Xf-radialPos)));		
		
	}
	else if (PHGMATH_Cosine(deltaRphi) > 0.1) {
		
		a = 1/(PHGMATH_Square(PHGMATH_Cosine(deltaRphi)));
		b = (-2.0*radialPos)*((1/PHGMATH_Square(PHGMATH_Cosine(deltaRphi))-1.0));
		c = (PHGMATH_Square(radialPos)*((1/PHGMATH_Square(PHGMATH_Cosine(deltaRphi))-1.0)))-PHGMATH_Square(phoTrkRadiusOfFocalCircle);
		
		/* Note, unless in debug mode, we know that there will be one real root solution here */
		numRoots = PhgMathSolveQuadratic(a, b, c, &minRoot, &maxRoot);


		#ifdef PHG_DEBUG
		{
			
			if (numRoots != 2) {
				PhgAbort("Did not get 2 roots from solution to quadratic (phoTrkCalcCBFDOmega)", false);
			}
			
			/* Verify the roots solve the quadratic equation */
			solution1 = ((a*PHGMATH_Square(minRoot)) + (b*minRoot) + c);
			
			if (PhgMathRealNumAreEqual(solution1, 0.0, 
		        -4, 0, 0, 0) == false) {
				PhgAbort("Root from quadratic do not solve properly (minRoot)", false);
			}
			solution2 = ((a*PHGMATH_Square(maxRoot)) + (b*maxRoot) + c);
			if (PhgMathRealNumAreEqual(solution2, 0.0, 
		        -4, 0, 0, 0) == false) {
				PhgAbort("Root from quadratic do not solve properly (maxRoot)", false);
			}
		}
		#endif
		
		Xf = minRoot;
		
		
		#ifdef PHG_DEBUG
		if ((maxRoot < radialPos) || (minRoot > radialPos)) {
			PhgAbort("Roots from quadratic do not staddle radialPos in cos(deltaRphi) > 0.1 (phoTrkCalcCBFDOmega)", false);
		}
		#endif
				
		omega = atan(zPos*fabs(PHGMATH_Cosine(deltaRphi)/(Xf-radialPos)));		
	}
	else  if (PHGMATH_Sine(deltaRphi) < 0.0) {
		
		a = 1/(PHGMATH_Square(PHGMATH_Sine(deltaRphi)));
		b = (-2.0*radialPos)*((1/PHGMATH_Square(PHGMATH_Sine(deltaRphi))-1.0));		
		c = PHGMATH_Square(radialPos) - PHGMATH_Square(phoTrkRadiusOfFocalCircle);
		
		/* Note, unless in debug mode, we know that there will be one real root solution here */
		numRoots = PhgMathSolveQuadratic(a, b, c, &minRoot, &maxRoot);

		#ifdef PHG_DEBUG
		{
			if (numRoots != 2) {
				PhgAbort("Did not get 2 roots from solution to quadratic (phoTrkCalcCBFDOmega)", false);
			}
			
			/* Verify the roots solve the quadratic equation */
			solution1 = ((a*PHGMATH_Square(minRoot)) + (b*minRoot) + c);
			
			
			if (PhgMathRealNumAreEqual(solution1, 0.0, 
		        -4, 0, 0, 0) == false) {
				PhgAbort("Root from quadratic do not solve properly (minRoot)", false);
			}
			solution2 = ((a*PHGMATH_Square(maxRoot)) + (b*maxRoot) + c);
			if (PhgMathRealNumAreEqual(solution2, 0.0, 
		        -4, 0, 0, 0) == false) {
				PhgAbort("Root from quadratic do not solve properly (maxRoot)", false);
			}
		}
		#endif

		Yf = maxRoot;
		
		#ifdef PHG_DEBUG
		if ((maxRoot < 0.0) || (minRoot > 0.0)) {
			PhgAbort("Roots from quadratic do not staddle 0.0 in sin(deltaRphi) > 0.1 (phoTrkCalcCBFDOmega)", false);
		}
		#endif
		
		omega = atan(zPos*fabs(PHGMATH_Sine(deltaRphi)/Yf));
	}
	else  {
		
		a = 1/(PHGMATH_Square(PHGMATH_Sine(deltaRphi)));
		b = (-2.0*radialPos)*((1/PHGMATH_Square(PHGMATH_Sine(deltaRphi))-1.0));		
		c = PHGMATH_Square(radialPos) - PHGMATH_Square(phoTrkRadiusOfFocalCircle);
		
		/* Note, unless in debug mode, we know that there will be one real root solution here */
		numRoots = PhgMathSolveQuadratic(a, b, c, &minRoot, &maxRoot);

		#ifdef PHG_DEBUG
		{
			if (numRoots != 2) {
				PhgAbort("Did not get 2 roots from solution to quadratic (phoTrkCalcCBFDOmega)", false);
			}
			
			
			/* Verify the roots solve the quadratic equation */
			solution1 = ((a*PHGMATH_Square(minRoot)) + (b*minRoot) + c);
			
			if (PhgMathRealNumAreEqual(solution1, 0.0, 
		        -4, 0, 0, 0) == false) {
				PhgAbort("Root from quadratic do not solve properly (minRoot)", false);
			}
			solution2 = ((a*PHGMATH_Square(maxRoot)) + (b*maxRoot) + c);
			if (PhgMathRealNumAreEqual(solution2, 0.0, 
		        -4, 0, 0, 0) == false) {
				PhgAbort("Root from quadratic do not solve properly (maxRoot)", false);
			}
		#endif

		Yf = minRoot;
		
		#ifdef PHG_DEBUG
		if ((maxRoot < 0.0) || (minRoot > 0.0)) {
			PhgAbort("Roots from quadratic do not staddle 0.0 in sin(deltaRphi) < 0.1 (phoTrkCalcCBFDOmega)", false);
		}
		}
		#endif
		
		omega = atan(zPos*fabs(PHGMATH_Sine(deltaRphi)/Yf));
	}

	
	return(omega);
}

/*********************************************************************************
*
*			Name:		PhoTrkTerminate
*
*			Summary:	Terminate the photon tracking module.
*
*			Arguments:
*
*			Function return: None.
*
*********************************************************************************/
void PhoTrkTerminate()
{
	LbUsFourByte	ieiIndex;			/* Index for iei tables */
	LbUsFourByte	iwiIndex;			/* Index for iwi tables */

	do { /* Process Loop */
	
		/* Only deal with this if we've been initialized */
		if (phoTrkIsInitialized == true) {
			if (PHG_IsForcedDetection() == true) {
			
				/* Free the memory used by the FDT tables */
				for (ieiIndex = 0; ieiIndex < phoTrkFDTInfo.num_iei; ieiIndex++) {
					
					/* For each iwi element */
					for (iwiIndex = 0; iwiIndex < phoTrkFDTInfo.num_iwi; iwiIndex++) {
						
						/* Free an iwo table */
						if (phoTrkFDTable[ieiIndex][iwiIndex].iwoTable != 0)
							LbMmFree((void **) &(phoTrkFDTable[ieiIndex][iwiIndex].iwoTable));
									
						/* Free an ipo table */
						if (phoTrkFDTable[ieiIndex][iwiIndex].ipoTable != 0)
							LbMmFree((void **) &(phoTrkFDTable[ieiIndex][iwiIndex].ipoTable));
							
						/* Free an ipo_cum table */
						if (phoTrkFDTable[ieiIndex][iwiIndex].ipoCumTable != 0)
							LbMmFree((void **) &(phoTrkFDTable[ieiIndex][iwiIndex].ipoCumTable));
							
					}
					/* Free an iwi table */
					if (phoTrkFDTable[ieiIndex] != 0)
						LbMmFree((void **) &(phoTrkFDTable[ieiIndex]));
				}
				
				/* Free the iei table */
				if (phoTrkFDTable != 0) {
					LbMmFree((void **) &phoTrkFDTable);
				}
				
				/* Free the cone beam tables */
					if (phoTrkCBFDMaxDeltaMu != 0) {
						LbMmFree((void **) &phoTrkCBFDMaxDeltaMu);
					}
					if (phoTrkCBFDMinDeltaMu != 0) {
						LbMmFree((void **) &phoTrkCBFDMinDeltaMu);
					}
					if (phoTrkCBFDProbDeltaPhiMu != 0) {
						LbMmFree((void **) &phoTrkCBFDProbDeltaPhiMu);
					}
					if (phoTrkCBFDCumProbDeltaPhiMu != 0) {
						LbMmFree((void **) &phoTrkCBFDCumProbDeltaPhiMu);
					}
					if (phoTrkCBFDTotalProbAccept != 0) {
						LbMmFree((void **) &phoTrkCBFDTotalProbAccept);
					}

			}
			
			phoTrkIsInitialized = false;
		}
	} while (false);
}
	
/*********************************************************************************
*
*			Name:		PhoTrkUpdatePhotonPosition
*
*			Summary:	Update the tracking_photon to new_position.
*
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*				PHG_Position		*newPositionPtr		- The photon's new position.
*
*			Function return: None.
*
*********************************************************************************/
void PhoTrkUpdatePhotonPosition(PHG_TrackingPhoton *trackingPhotonPtr,
		PHG_Position  *newPositionPtr)
{
	double	distanceTraveled;	/* The distance traveled */
	
	/* Calculate the change in distance traveled to the new position */
	distanceTraveled = PHGMATH_SquareRoot(
		PHGMATH_Square(newPositionPtr->x_position - trackingPhotonPtr->location.x_position) +
					      PHGMATH_Square(newPositionPtr->y_position - trackingPhotonPtr->location.y_position) +
		PHGMATH_Square(newPositionPtr->z_position - trackingPhotonPtr->location.z_position));
		
	/* Update total distance traveled */
	trackingPhotonPtr->travel_distance += distanceTraveled;
	
	/* Set new position */
	trackingPhotonPtr->location = *newPositionPtr;
	
}
