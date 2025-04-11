/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1998-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		DetGeometric.c
*			Revision Number:	2.2
*			Date last revised:	20 June 2013
*			Programmer:			Steven Gillispie, Steven Vannoy
*			Date Originated:	21 November 2005
*
*			Module Overview:	Simulates general geometric detectors.
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
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s);		Steven Gillispie
*
*			Revision date:		20 June 2013
*
*			Revision description:	Moved detGeomDecidePhotonAction contents from 
*									here to PhoTrk.c as PhoTrkDecidePhotonAction.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		26 January 2012
*
*			Revision description:	Changed form of user functions to pointers
*
*********************************************************************************/

#define DETECTOR_GEOMETRIC


#include <stdio.h>
#include <memory.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhoTrk.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "PhgMath.h"
#include "CylPos.h"
#include "ProdTbl.h"
#include "SubObj.h"
#include "EmisList.h"
#include "PhoHStat.h"
#include "PhoHFile.h"
#include "ColUsr.h"
#include "UNCCollimator.h"
#include "Collimator.h"
#include "ColSlat.h"
#include "Detector.h"
#include "DetGeometric.h"
#include "DetPlanar.h"
#include "DetCylinder.h"
#include "DetBlock.h"
#include "DetUsr.h"
#include "phg.h"
#include "PhgBin.h"


void			detGeomInitDetectors(DetEn_DetectorTypeTy detectorType);
void			detGeomInitPhotons(DetEn_DetectorTypeTy detectorType, 
					PHG_Decay *decayPtr, PHG_TrackingPhoton *photonPtr);
void			detGeomEndDetection(DetEn_DetectorTypeTy detectorType, 
					PHG_TrackingPhoton *photonPtr);
Boolean			detGeomProjectToDetector(DetEn_DetectorTypeTy detectorType, 
					PHG_TrackingPhoton *photonPtr, 
					LbFourByte *ringNum, LbFourByte *detNum);
double			detGeomGtTruncatedFreePaths(DetEn_DetectorTypeTy detectorType, 
					PHG_TrackingPhoton *photonPtr, Boolean firstTime,
					LbUsFourByte ringNum, LbUsFourByte detNum, double *weight);
void			detGeomFindNextInteraction(DetEn_DetectorTypeTy detectorType, 
					PHG_TrackingPhoton *photonPtr, 
					Boolean firstTime, LbFourByte interactionNum, 
					LbFourByte *ringNum, LbFourByte *detNum, 
					detEn_ActionTy *actionType, double *fpToGo, 
					LbUsFourByte *detMaterial, Boolean *isActive);
#ifdef PHG_DEBUG
void			detGeomDumpInteractions(LbFourByte numInteractions,
					detInteractionInfoTy *interactions);
#endif
detEn_ActionTy	detGeomDecidePhotonAction(LbUsFourByte detMaterial, double photonEnergy, 
					Boolean modelingAbsorption, Boolean modelingCohScatter);
void			detGeomFindDetPosition(DetEn_DetectorTypeTy detectorType, 
					PHG_TrackingPhoton *photonPtr, LbUsFourByte interactionIndex);



/*********************************************************************************
*
*		Name:			detGeomInitDetectors
*
*		Summary:		Initialize the geometric detectors based on detector type.
*
*		Arguments:
*			DetEn_DetectorTypeTy	detectorType		- The detector type.
*
*		Function return:	None.
*
*********************************************************************************/

void detGeomInitDetectors(DetEn_DetectorTypeTy detectorType)

{
	switch (detectorType) {
		case DetEn_Planar:
		case DetEn_DualHeaded:
		case DetEn_Block:
			/* These detectors don't need to do anything */
			break;
		
		case DetEn_Cylindrical:
			DetCylInitCylinders();
			break;
		
		default:
			/* Do nothing */
			break;
	}
}


/*********************************************************************************
*
*		Name:			detGeomInitPhotons
*
*		Summary:		Initialize the photons based on detector type.
*
*		Arguments:
*			DetEn_DetectorTypeTy	detectorType	- The detector type.
*			PHG_Decay				*decayPtr		- The decayPtr that started the process.
*			PHG_TrackingPhoton 		*photonPtr		- The detected photon.
*
*		Function return:	None.
*
*********************************************************************************/

void detGeomInitPhotons(DetEn_DetectorTypeTy detectorType, 
							PHG_Decay *decayPtr, PHG_TrackingPhoton *photonPtr)

{
	switch (detectorType) {
		case DetEn_Planar:
		case DetEn_DualHeaded:
			DetPlnrInitPhotons(decayPtr, photonPtr);
			break;
		
		case DetEn_Cylindrical:
		case DetEn_Block:
			/* These detectors don't need to do anything */
			break;
		
		default:
			/* Do nothing */
			break;
	}
}


/*********************************************************************************
*
*		Name:			detGeomEndDetection
*
*		Summary:		Set the final detected position.
*
*		Arguments:
*			DetEn_DetectorTypeTy	detectorType	- The detector type.
*			PHG_TrackingPhoton 		*photonPtr		- The detected photon.
*
*		Function return:	None.
*
*********************************************************************************/

void detGeomEndDetection(DetEn_DetectorTypeTy detectorType, 
							PHG_TrackingPhoton *photonPtr)

{
	switch (detectorType) {
		case DetEn_Planar:
		case DetEn_DualHeaded:
			/* Convert centroid location to tomo coordinates */
			photonPtr->location.x_position = 
				(DetRunTimeParams[DetCurParams].PlanarDetector.InnerRadius*PHGMATH_Cosine(photonPtr->detectorAngle))+
				(photonPtr->detLocation.x_position *
				PHGMATH_Cosine(photonPtr->detectorAngle)) -
				(photonPtr->detLocation.y_position * PHGMATH_Sine(photonPtr->detectorAngle));

			photonPtr->location.y_position =
				(DetRunTimeParams[DetCurParams].PlanarDetector.InnerRadius*PHGMATH_Sine(photonPtr->detectorAngle))+
				(photonPtr->detLocation.x_position *
				PHGMATH_Sine(photonPtr->detectorAngle)) +
				(photonPtr->detLocation.y_position * PHGMATH_Cosine(photonPtr->detectorAngle));
			
			photonPtr->location.z_position = photonPtr->detLocation.z_position;
			break;
		
		
		case DetEn_Cylindrical:
		case DetEn_Block:
			/* Copy the detected position */
			photonPtr->location.x_position = photonPtr->detLocation.x_position;
			photonPtr->location.y_position = photonPtr->detLocation.y_position;
			photonPtr->location.z_position = photonPtr->detLocation.z_position;
			break;
		
		
		default:
			/* Do nothing */
			break;
	}
}


/*********************************************************************************
*
*			Name:			DetGeometric
*
*			Summary:		Perform detection for geometric detectors.
*
*			Arguments:
*				DetEn_DetectorTypeTy detectorType	- The type of geometric detector.
*				PHG_Decay			*decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton	*bluePhotons	- The blue photons detected.
*				LbUsFourByte 		numBluePhotons	- The number of blue photons.
*				PHG_TrackingPhoton	*pinkPhotons	- The pink photons detected.
*				LbUsFourByte		numPinkPhotons	- The number of pink photons.
*				DetectedPhotonsTy	*detPhotonsPtr	- The accepted detected photons.
*
*			Function return:	None
*
*********************************************************************************/

void DetGeometric(DetEn_DetectorTypeTy detectorType, 
		PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		DetectedPhotonsTy *detPhotonsPtr)

{
	PHG_TrackingPhoton 		*photons;		/* Used to access both blues and pinks */
	LbUsFourByte			numPhotons;		/* Used for both blues and pinks */
	LbUsFourByte			loopCount;		/* Current time through the loop once for pinks and once for blues */
	LbUsFourByte			pIndex;			/* Current photon */
	
	
	/* Clear the counters */
	detPhotonsPtr->NumDetectedBluePhotons = 0;
	detPhotonsPtr->NumDetectedPinkPhotons = 0;
	detData[DetCurParams].detDetectedBluePhotonIndex = 0;
	detData[DetCurParams].detDetectedPinkPhotonIndex = 0;
	#ifdef PHG_DEBUG
	detData[DetCurParams].DetDiscardBlue = true;
	detData[DetCurParams].DetDiscardPink = true;
	#endif
	
	/* Initialize photon pointer and loop control to be blue photons */
	photons = bluePhotons;
	numPhotons = numBluePhotons;
	loopCount = 1;
	
	do { /* Loop for blue and pink photons */
		
		/* Loop through all same-colored photons */
		for (pIndex = 0; pIndex < numPhotons; pIndex++) {
			
			/* Initialize photons as needed */
			detGeomInitPhotons(detectorType, decayPtr, &(photons[pIndex]));
			
			/* Let user modify and/or reject photons */
			if (PHG_IsSPECT()) {
				if (DetUsrModSPECTPhotonsFPtr && 
						(*DetUsrModSPECTPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
							&(photons[pIndex])) == false) {
					
					/* They rejected it so go to next photon */
					continue;
				}
			}
			else {
				if (DetUsrModPETPhotonsFPtr && 
						(*DetUsrModPETPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
							&(photons[pIndex])) == false) {
					
					/* They rejected it so go to next photon */
					continue;
				}
			}
			
			/* Track the photon */
			if (DetGeomTrack(detectorType, decayPtr, &(photons[pIndex]))) {
				
				/* We made it to here so save the detected photons */
				if (photons == bluePhotons) {
					/* Save a blue photon */
					#ifdef PHG_DEBUG
						detData[DetCurParams].DetDiscardBlue = false;
					#endif
					detPhotonsPtr->DetectedTrkngBluePhotons[detPhotonsPtr->NumDetectedBluePhotons]
						= photons[pIndex];
					
					/* Increment the counter */
					detPhotonsPtr->NumDetectedBluePhotons++;
				}
				else {
					/* Save a pink photon */
					#ifdef PHG_DEBUG
						detData[DetCurParams].DetDiscardPink = false;
					#endif
					detPhotonsPtr->DetectedTrkngPinkPhotons[detPhotonsPtr->NumDetectedPinkPhotons]
						= photons[pIndex];
					
					/* Increment the counter */
					detPhotonsPtr->NumDetectedPinkPhotons++;
				}
				
				/* Update our detected photon block if doing history file */
				if  (DET_IsDoHistory()) {
					DetUpdateDetectedPhotonBlock(&photons[pIndex]);
				}
			}
		}
		
		/* if no blue photons were found and this is a coincidence-only
			simulation, break out of the loop */
		if ( PHG_IsPETCoincidencesOnly() && (detPhotonsPtr->NumDetectedBluePhotons == 0) ) {
			
			break;
			
		}
		
		/* Now set up for pink photons */
		photons = pinkPhotons;
		numPhotons = numPinkPhotons;
		loopCount++;
		
	} while (loopCount <= 2);
	
	/* Both blue and (if appropriate) pink photons all now detected */
	
}


/*********************************************************************************
*
*			Name:			DetGeomTrack
*
*			Summary:		Perform tracking for geometric detectors.
*
*			Arguments:
*				DetEn_DetectorTypeTy detectorType	- The type of geometric detector.
*				PHG_Decay			*decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton 	*photonPtr		- The photon to track and detect.
*
*			Function return: True if photon should be accepted.
*
*********************************************************************************/

Boolean DetGeomTrack(DetEn_DetectorTypeTy detectorType, 
						PHG_Decay *decayPtr, PHG_TrackingPhoton *photonPtr)

{
	Boolean				acceptPhoton = false;	/* Flag for acceptance of photon */
	Boolean				modelingAbsorption;		/* Whether absorption is being modeled or not */
	Boolean				modelingCohScatter;		/* Whether coherent scatter is modeled or not */
	LbFourByte			ringNum;				/* Detector ring at interaction point */
	LbFourByte			detNum;					/* Detector number at interaction point */
	#ifdef PHG_DEBUG
	double				depositedEnergy;		/* Total energy deposited by photon */
	#endif
	double				depositedActiveEnergy;	/* Energy deposited in active layers */
	LbFourByte			curInteraction;			/* Current interaction number */
	double				fpToGo;					/* Free paths to travel */
	Boolean				firstTime;				/* Flag for first round of tracking */
	detEn_ActionTy		action;					/* Type of interaction of photon */
	LbUsFourByte		detMaterial;			/* Material of detector at interaction point */
	Boolean				isActive;				/* Whether the interaction layer is active */
	double				initialPhotonEnergy;	/* Photon energy before interaction */
	double				energyDepositedFromInteraction;
												/* Energy deposited by interaction */
	
	
	/* Eliminate compiler warning about unused parameter */
	{
		PHG_Decay			*localDecayPtr;			/* Local copy of decayPtr */
		
		localDecayPtr = decayPtr;
	}
	
	/* Collect the modeling parameters for simpler reference */
	modelingAbsorption = true;
	modelingCohScatter = PHG_IsModelCoherentInTomo();
	
	/* Initialize the geometric detector */
	detGeomInitDetectors(detectorType);
	
	/* The target cylinder of the PHG and the inner-most surface of the
		detector may be different.  Hence, we will first project to
		the inner-most surface.
	*/
	if (! detGeomProjectToDetector(detectorType, photonPtr, &ringNum, &detNum)) {
		/* Photon exited without being detected */
		goto REJECT;
	}
	
	/* Increment count of photons that reach the detector */
	detData[DetCurParams].detTotReachingCrystal++;
	
	
	/*	Now track the photon through the detector until it is absorbed
		or it escapes.  At each interaction point tally the energy deposited
		and record the position of the event (if the interaction occurs in an
		active layer (see below))
	*/
	
	/* Clear tracking variables */
	#ifdef PHG_DEBUG
	depositedEnergy = 0.0;
	#endif
	depositedActiveEnergy = 0.0;
	curInteraction = -1;
	
	/* Get free paths to travel */  
	if (DetRunTimeParams[DetCurParams].DoForcedInteraction) {
		fpToGo = detGeomGtTruncatedFreePaths(detectorType, photonPtr, true, 
					ringNum, detNum, &(photonPtr->photon_current_weight));
	} else {
		PhgMathGetTotalFreePaths(&fpToGo);
	}
	
	/* Enter the tracking loop */
	firstTime = true;
	do {
		
		/* Find and move to the next interaction point */
		detGeomFindNextInteraction(detectorType, photonPtr, firstTime, curInteraction, 
							&ringNum, &detNum, &action, &fpToGo, &detMaterial, &isActive);
		
		
		/* If the next interaction really was an interaction (photon didn't exit), 
			determine whether it should be an absorption or scatter, and act accordingly
		*/
		if (action == detEnAc_Interact) {
			/* Determine the photon's type of interaction */
			action = detGeomDecidePhotonAction(detMaterial, photonPtr->energy, 
							modelingAbsorption, modelingCohScatter);
			
			/* Increment counts of interactions */
			curInteraction++;
			photonPtr->num_det_interactions = curInteraction + 1;	/* start from 1, not 0 */
			#ifdef PHG_DEBUG
			detData[DetCurParams].detCylCountInteractions++;	/* ### Make universal */
			#endif
			
			if (curInteraction == MAX_DET_INTERACTIONS) {
				/* Absorb the photon immediately: 
					set the action to absorption */
				action = detEnAc_Absorb;
			}
			
			/* Record the new interaction position */
			photonPtr->det_interactions[curInteraction].pos = photonPtr->location;
			photonPtr->det_interactions[curInteraction].posIndices.ringNum = ringNum;
			photonPtr->det_interactions[curInteraction].posIndices.blockNum = detNum;
			photonPtr->det_interactions[curInteraction].posIndices.layerNum = 0;
			photonPtr->det_interactions[curInteraction].posIndices.elementNum = 0;
			
			
			/* Save the current energy */
			initialPhotonEnergy = photonPtr->energy;
			
			/* Proceed according to the resultant action type */
			switch (action) {
				case detEnAc_Absorb:
				{
					/* Absorb the photon */
					photonPtr->energy = 0.0;
					
					break;
				}
				
				
				case detEnAc_CohScatter:
				{
					/* Coherent scatter */
					/*	Direction may be changed; energy is not */
					
					EmisListDoCoherent(photonPtr, detMaterial);
					#ifdef PHG_DEBUG
					detData[DetCurParams].detCylCountCohInteractions++;	/* ### Make universal */
					#endif
					
					break;
				}
				
				
				case detEnAc_ComptonScatter:
				{
					/* Compton scatter */
					/*	Both direction and energy are likely to change */
					
					EmisListDoComptonInteraction(photonPtr);
					
					break;
				}
				
				
				default:
					/* Ignore other values (they should not occur) */
					break;
			}
			
			if (action != detEnAc_Absorb) {
				/* If the energy is now too small, an absorption is forced */
				if (photonPtr->energy < PHG_MIN_PHOTON_ENERGY) {
					
					/* Change action to absorption */
					action = detEnAc_Absorb;
					
					/* Absorb the photon */
					photonPtr->energy = 0.0;
					
					detData[DetCurParams].detTotForcedAbsorptions++;
				}
			}
			
			/* Compute the deposited energy */
			energyDepositedFromInteraction = initialPhotonEnergy - photonPtr->energy;
			
			
			/* Record the interaction energy if it occured in an active detector */
			if (isActive) {
				/* Mark interaction as in an active layer */
				photonPtr->det_interactions[curInteraction].isActive = true;
				
				/* Record the deposited energy */
				photonPtr->det_interactions[curInteraction].energy_deposited =
					energyDepositedFromInteraction;
				
				/* Sum the deposited energy */
				#ifdef PHG_DEBUG
				depositedEnergy += energyDepositedFromInteraction;
				#endif
				depositedActiveEnergy += energyDepositedFromInteraction;
			}
			else {
				photonPtr->det_interactions[curInteraction].isActive = false;
				photonPtr->det_interactions[curInteraction].energy_deposited = 0.0;
			}
			
			if (action != detEnAc_Absorb) {
				/*  Choose free paths to the next interaction */
				PhgMathGetTotalFreePaths(&fpToGo);
			}
		}
		
		/* Clear "first time" flag */
		firstTime = false;
		
	} while ((action != detEnAc_Discard) && (action != detEnAc_Absorb));
	
	
	/* If photon was absorbed, update statistics */
	if (action == detEnAc_Absorb) {
		double			weight;		/* Local copy for efficiency */
		
		/* Increment counters */
		detData[DetCurParams].detTotPhotonsAbsorbed++;
		if ((curInteraction+1) == MAX_DET_INTERACTIONS) {
			detData[DetCurParams].detNumReachedMaxInteractions++;
		}
		
		weight = photonPtr->photon_current_weight;
		
		/* Increment weight counter */
		detData[DetCurParams].detTotWtAbsorbed += (weight*photonPtr->decay_weight);
		
		/* Update absorbed array, as long as interaction count isn't too high */
		if (curInteraction < MAX_DET_INTERACTIONS) {
			detData[DetCurParams].detWeightAbsorbedBins[curInteraction+1] += (weight*photonPtr->decay_weight);
		}
		
		/* If it is the first interaction, update statistics for 
			absorption on first interaction
		*/
		if (curInteraction == 0) {
			detData[DetCurParams].detTotFirstTimeAbsorptions++;
			detData[DetCurParams].detTotWtFirstTimeAbsorbed += (weight*photonPtr->decay_weight);
		}
	}
	
	
	/* If energy was deposited then compute centroid,
		possibly blur energy, and add the photon to
		the detected photon's list
	 */
	/* Compute centroid if energy deposited */
	if (depositedActiveEnergy > 0.0){	
		
		/* Increment counter */
		detData[DetCurParams].detTotPhotonsDepositingEnergy++;
		
		/* Determine the detected position and save it */
		detGeomFindDetPosition(detectorType, photonPtr, curInteraction);
		detGeomEndDetection(detectorType, photonPtr);
		
		/* Clear out the direction vector, as it is no longer meaningful */
		photonPtr->angle.cosine_x = 0.0;
		photonPtr->angle.cosine_y = 0.0;
		photonPtr->angle.cosine_z = 0.0;
		
		/* Blur energy if it was requested */
		if (DET_DoEnergyBlur() == true) {
			photonPtr->energy = DetGaussEnergyBlur(photonPtr->energy);				
		}
		
		/* Blur time if it was requested */
		if ( DET_DoTofBlur() == true ) {
			photonPtr->travel_distance = DetGaussTimeBlur(photonPtr->travel_distance);				
		}
		
		/* Mark photon for acceptance */
		acceptPhoton = true;
	}
	
	
	/* Increment pass-through counter if appropriate */
	if (curInteraction < 0) {
		/* Increment pass through counter */
		detData[DetCurParams].detTotPhotonsPassingThrough++;
	}
	
	
REJECT:;
	return(acceptPhoton);
}


/*********************************************************************************
*
*			Name:			detGeomProjectToDetector
*
*			Summary:		Project to the inner-most surface of the detector.
*
*			Arguments:
*				DetEn_DetectorTypeTy	detectorType	- The detector type.
*				PHG_TrackingPhoton 		*photonPtr		- The photon to track and detect.
*				LbUsFourByte			*ringNum		- The ring projected to.
*				LbUsFourByte			*detNum			- The detector projected to.
*
*			Function return:	True if photon not rejected (still valid).
*
*********************************************************************************/

Boolean detGeomProjectToDetector(DetEn_DetectorTypeTy detectorType, 
									PHG_TrackingPhoton *photonPtr, 
									LbFourByte *ringNum, LbFourByte *detNum)

{
	Boolean		result;				/* True if photon still valid */
	
	
	result = true;	/* ### */
	switch (detectorType) {
		case DetEn_Planar:
		case DetEn_DualHeaded:
			result = DetPlnrProjectToDetector(photonPtr, ringNum, detNum);
			break;
		
		case DetEn_Cylindrical:
			result = DetCylProjectToDetector(photonPtr, ringNum, detNum);
			break;
		
		case DetEn_Block:
			result = DetBlocProjectToDetector(photonPtr, ringNum, detNum);
			break;
		
		default:
			/* Do nothing */
			break;
	}
	
	return(result);
}


/*********************************************************************************
*
*		Name:			detGeomGtTruncatedFreePaths
*
*		Summary:		Return the free paths to travel for a forced first interaction.
*
*		Arguments:		
*			DetEn_DetectorTypeTy	detectorType	- The detector type.
*			PHG_TrackingPhoton		*photonPtr		- The photon.
*			Boolean					firstTime		- Optimizer for tracking.
*			LbUsFourByte			ringNum			- The current ring.
*			LbUsFourByte			detNum			- The current detector.
*			double					*weight			- The photon weight.
*
*		Function return: Free-paths to travel.
*
*********************************************************************************/

double detGeomGtTruncatedFreePaths(DetEn_DetectorTypeTy detectorType, 
						PHG_TrackingPhoton *photonPtr, Boolean firstTime,
						LbUsFourByte ringNum, LbUsFourByte detNum, double *weight)

{
	double			fpToGo;			/* Returned value of free paths */
	
	
	fpToGo = 0.0;
	
	switch (detectorType) {
		case DetEn_Planar:
		case DetEn_DualHeaded:
			fpToGo = DetPlnrGtTruncatedFreePaths(photonPtr, weight);
			break;
		
		case DetEn_Cylindrical:
			fpToGo = DetCylGtTruncatedFreePaths(photonPtr, firstTime, 
												ringNum, detNum, weight);
			break;
		
		case DetEn_Block:
			fpToGo = DetBlocGtTruncatedFreePaths(photonPtr, ringNum, detNum, weight);
			break;
		
		default:
			/* Do nothing */
			break;
	}
	
	return (fpToGo);
}


/*********************************************************************************
*
*		Name:			detGeomFindNextInteraction
*
*		Summary:		Find and move the photon to the next interaction point.
*
*		Arguments:
*			DetEn_DetectorTypeTy	detectorType	- The detector type.
*			PHG_TrackingPhoton 		*photonPtr		- The photon to track and detect.
*			Boolean					firstTime		- Optimizer for tracking.
*			LbFourByte				interactionNum	- Current interaction count.
*			LbFourByte				*ringNum		- The ring projected to.
*			LbFourByte				*detNum			- The detector projected to.
*			detEn_ActionTy			*actionType		- The type of interaction that occurred.
*			double					*fpToGo			- Remaining free paths to travel.
*			LbUsFourByte			*detMaterial	- The interaction material of the detector.
*			Boolean					*isActive		- Active status of the interaction material.
*
*		Function return:	None.
*
*********************************************************************************/

void detGeomFindNextInteraction(DetEn_DetectorTypeTy detectorType, 
									PHG_TrackingPhoton *photonPtr, 
									Boolean firstTime, LbFourByte interactionNum, 
									LbFourByte *ringNum, LbFourByte *detNum, 
									detEn_ActionTy *actionType, double *fpToGo, 
									LbUsFourByte *detMaterial, Boolean *isActive)

{
	switch (detectorType) {
		case DetEn_Planar:
		case DetEn_DualHeaded:
			DetPlnrFindNextInteraction(photonPtr, interactionNum, detNum, 
										actionType, fpToGo, detMaterial, isActive);
			break;
		
		case DetEn_Cylindrical:
			DetCylFindNextInteraction(photonPtr, firstTime, ringNum, detNum, 
										actionType, fpToGo, detMaterial, isActive);
			break;
		
		case DetEn_Block:
			DetBlocFindNextInteraction(photonPtr, ringNum, detNum, 
										actionType, fpToGo, detMaterial, isActive);
			break;
		
		default:
			/* Do nothing */
			break;
	}
}


#ifdef PHG_DEBUG
/*********************************************************************************
*
*		Name:			detGeomDumpInteractions
*
*		Summary:		Prints out interaction list for debugging purposes.
*
*		Arguments:
*			LbFourByte					numInteractions	- The number of interactions.
*			detGeomInteractionInfoTy	interactions	- The interactions.
*
*		Function return: None.
*
*********************************************************************************/
void detGeomDumpInteractions(LbFourByte numInteractions, detInteractionInfoTy *interactions)
{
	LbFourByte curInteraction = 0;
	
	/* Loop through interactions printing info for each one */
	while (curInteraction < numInteractions) {
		LbInPrintf("\nFor interaction %d", curInteraction);
		LbInPrintf("\n\tpos.x = %3.2f", interactions[curInteraction].pos.x_position);
		LbInPrintf("\n\tpos.y = %3.2f", interactions[curInteraction].pos.y_position);
		LbInPrintf("\n\tpos.z = %3.2f", interactions[curInteraction].pos.z_position);
		LbInPrintf("\n\tenergy = %3.2f", interactions[curInteraction].energy_deposited);
		
		curInteraction++;
	}
}
#endif


/*********************************************************************************
*
*		Name:			detGeomDecidePhotonAction
*
*		Summary:		Decide, probabilistically, what action the photon will undergo.
*							
*		Arguments:
*			LbUsFourByte	detMaterial			- The interaction material of the detector.
*			double			photonEnergy		- The energy of the photon.
*			Boolean			modelingAbsorption	- If absorption is simulated.
*			Boolean			modelingCohScatter	- If coherent scatter is simulated.
*
*		Function return: The detEn_ActionTy type of photon action.
*
*********************************************************************************/

detEn_ActionTy detGeomDecidePhotonAction(LbUsFourByte detMaterial, double photonEnergy, 
							Boolean modelingAbsorption, Boolean modelingCohScatter)

{
	phtrkEn_ActionTy	phtrkAction;			/* Obtained result */
	detEn_ActionTy		detAction;				/* Returned result */
	
	
	/* Use PhoTrkDecidePhotonAction */
	phtrkAction = PhoTrkDecidePhotonAction(detMaterial, photonEnergy, 
										modelingAbsorption, modelingCohScatter);
	
	/* Switch to detector types for returned result */
	switch (phtrkAction) {
		case phtrkEnAc_Interact:		detAction = detEnAc_Interact; break;
		case phtrkEnAc_Absorb:			detAction = detEnAc_Absorb; break;
		case phtrkEnAc_ComptonScatter:	detAction = detEnAc_ComptonScatter; break;
		case phtrkEnAc_CohScatter:		detAction = detEnAc_CohScatter; break;
		default:						detAction = detEnAc_Null; break;
	}
	
	return(detAction);
}


/*********************************************************************************
*
*		Name:			detGeomFindDetPosition
*
*		Summary:		Find the photon's detected position.
*
*		Arguments:
*			DetEn_DetectorTypeTy	detectorType		- The detector type.
*			PHG_TrackingPhoton		*photonPtr			- The photon.
*			LbUsFourByte 			interactionIndex	- The number of interactions.
*
*		Function return: None.
*
*********************************************************************************/

void detGeomFindDetPosition(DetEn_DetectorTypeTy detectorType, 
							PHG_TrackingPhoton *photonPtr, LbUsFourByte interactionIndex)

{
	Boolean				interactsSet[MAX_DET_INTERACTIONS];	/* Which interactions to use */
	LbUsFourByte		i;									/* Index through interactions */
	
	
	/*	Update the photon to reflect its status after passing through
		the detector.  Note that it would not be valid to continue to track
		this photon through any remaining space.
	*/
	if (detectorType == DetEn_Block) {
		/* Special handling required for block detectors */
		DetBlocFindDetPosition(photonPtr, interactionIndex);
	}
	else {
		/* Standard computation */
		
		/* Indicate the requested set of interactions */
		for (i=0; i<=interactionIndex; i++) {
			/* Request all */
			interactsSet[i] = true;
		}
		
		/* Compute the centroid and deposited energy for the requested set */
		DetGeomCompCentroid(photonPtr, 
								&(photonPtr->detLocation), &(photonPtr->energy), 
								interactionIndex, interactsSet);
		
		/* Save the y position to *temporarily* make Planar computing easier */
		photonPtr->transaxialPosition = photonPtr->detLocation.y_position;
	}
}


/*********************************************************************************
*
*		Name:			DetGeomCompCentroid
*
*		Summary:		Compute the energy-weighted centroid position for the 
*							supplied set of interactions.
*
*		Arguments:
*			PHG_TrackingPhoton	*photonPtr			- The photon.
*			PHG_Position		*centroidPos		- Computed centroid position.
*			double				*depositedEnergy	- The deposited energy.
*			LbUsFourByte 		interactionIndex	- The number of interactions.
*			Boolean				*interactsSet		- Which interactions to use.
*
*		Function return: None.
*
*********************************************************************************/

void DetGeomCompCentroid(PHG_TrackingPhoton *photonPtr, 
							PHG_Position *centroidPos, double *depositedEnergy, 
							LbUsFourByte interactionIndex, Boolean *interactsSet)

{
	PHG_Position			centPos;					/* The computed centroid position */
	double					depEnergy;					/* The amount of energy deposited */
	LbFourByte				curInteraction;				/* Index through interactions */
	detInteractionInfoTy	*curInteractPtr;			/* Pointer to current interaction */
	double					photonEnergyDeposited;		/* Single photon energy deposited */
	
	
	/* Clear the calculated variables */
	centPos.x_position = 0.0;
	centPos.y_position = 0.0;
	centPos.z_position = 0.0;
	depEnergy = 0.0;
	
	
	/* Loop through interactions */
	curInteraction = interactionIndex;
	while (curInteraction >= 0) {
		if (interactsSet[curInteraction]) {
			/* Only look at interactions in the requested set */
			
			curInteractPtr = &(photonPtr->det_interactions[curInteraction]);
			if (curInteractPtr->isActive) {
				/* Only contribute to centroid if interaction was in an active layer */
				
				photonEnergyDeposited = curInteractPtr->energy_deposited;
				depEnergy += photonEnergyDeposited;
				
				centPos.x_position += 
					(curInteractPtr->pos.x_position * photonEnergyDeposited);
				centPos.y_position += 
					(curInteractPtr->pos.y_position * photonEnergyDeposited);
				centPos.z_position += 
					(curInteractPtr->pos.z_position * photonEnergyDeposited);
			}
		}
		
		curInteraction--;
	}
	
	if (depEnergy > 0.0) {
		/* Complete computation of centPos */
		centPos.x_position /= depEnergy;
		centPos.y_position /= depEnergy;
		centPos.z_position /= depEnergy;
	}
	
	
	/* Return the computed values */
	*centroidPos = centPos;
	*depositedEnergy = depEnergy;
}
