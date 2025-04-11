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
*			Module Name:		DetCylinder.c
*			Revision Number:	2.3
*			Date last revised:	4 June 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	January 20, 1998
*
*			Module Overview:	Simulates cylindrical detector functionality.
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
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		26 January 2012
*
*			Revision description:	Changed form of user functions to pointers
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
*						- support for simulate_PET_coincidences_only option
*						- support for DetGeometric
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		17 July 2004
*
*			Revision description:
*						- Added time-of-flight blurring.
*						- Corrected some precision problems that caused
*						history file processing to fail.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		21 Nov 2005
*
*			Revision description:
*						- Fixed bug in energy for coherent scatter that
*						lead to an inaccurate photon energy.
*
*********************************************************************************/

#define DETECTOR_CYLINDER


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
#include "DetCylinder.h"
#include "Detector.h"
#include "DetUsr.h"
#include "DetPlanar.h"
#include "phg.h"
#include "PhgBin.h"

typedef enum { detCylEnAc_Null,
				detCylEnAc_Detect,
				detCylEnAc_Discard,
				detCylEnAc_Escape,
				detCylEnAc_Interact,
				detCylEnAc_Absorb,
				detCylEnAc_AxialCross,
				detCylEnAc_LayerCross} detCylEn_ActionTy;

#ifdef PHG_DEBUG
static char					detCylErrStr[1024];				/* Storage for creating error strings */
#endif

void				detCylComputeFreePathsToExit(PHG_Position pos,
						PHG_Direction dir, double energy,
						double *fpPtr);
Boolean				detTrackCylinder(PHG_TrackingPhoton *photonPtr);
void 				detDHDoRandomPos(PHG_Decay *decayPtr,
						PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
						PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
						DetectedPhotonsTy *detPhotonsPtr);

detCylEn_ActionTy	detCylProject(PHG_TrackingPhoton *photonPtr,
						DetCylnRingTy *curRingInfo, Boolean firstTime, 
						PHG_Position *newPosPtr, LbFourByte *curRingIndex,
						LbFourByte *curLayerPtr, double *distPtr);

void				detCylCompCentroid(PHG_TrackingPhoton *photonPtr, LbUsFourByte interactionIndex,
						double depositedEnergy);
void 				detDualHeaded(PHG_Decay *decayPtr,
						PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
						PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
						DetectedPhotonsTy *detPhotonsPtr);
void				detDumpCylInteractions(LbFourByte numInteractions,
						detInteractionInfoTy *interactions);

void				detCylInitCylinders(LbUsFourByte curRing, LbUsFourByte curLayer);
double				detCylFindInnerDist(PHG_Position *posPtr,
						PHG_Direction *dirPtr);
LbUsFourByte		detCylFindRing(double zPos);
double				detCylGtFreePaths(PHG_TrackingPhoton *photonPtr, Boolean firstTime,
						LbUsFourByte curRingIndex, LbUsFourByte curLayerIndex);
void				detCylCompFreePathsToExit(PHG_TrackingPhoton *photonPtr,
						Boolean firstTime, 
						LbFourByte curRingIndex,
						LbFourByte curLayerIndex, double *fpToGo);
void 				detCylFindLayer(PHG_Position *pos, LbFourByte ringIndex, LbFourByte *theLayer);


static double				detCylInnerRadius;					/* The inner most radius of the detector */
static double				detCylOuterRadius;					/* The outer most radius */
static double				detCylLayerDepth;					/* The layer depth of the current ring */
static CylPosCylinderTy		detCylInBoundCyl;					/* The inner bounding cylinder for the detector */
static CylPosCylinderTy		detCylOutBoundCyl;					/* The outer bounding cylinder for the detector */


/*********************************************************************************
*
*			Name:			DetCylGtTruncatedFreePaths
*
*			Summary:		Returns the free paths to travel for a forced first interaction.
*
*			Arguments:
*							PHG_TrackingPhotonPtr	photonPtr		- The photon.
*							Boolean					firstTime		- Optimizer for tracking.
*							LbUsFourByte			curRingIndex	- The current ring.
*							LbUsFourByte			curLayerIndex	- The current layer.
*							double					*weight			- The photon weight.
*			Function return: Free-paths to travel.
*
*********************************************************************************/
double	DetCylGtTruncatedFreePaths(PHG_TrackingPhoton *photonPtr, Boolean firstTime,
			LbUsFourByte curRingIndex, LbUsFourByte curLayerIndex, double *weight)
{

	double 			newWeight;
	double 			freePathsToExit;
	double			fpToGo;
	double			randFromExp;
	LbUsFourByte	wholeNum;				/* Whole number temp variable */
	
	/* Compute the free paths to exit */
	detCylCompFreePathsToExit(photonPtr,
				firstTime, 
				curRingIndex,
				curLayerIndex, &freePathsToExit);
						
	/* Adjust the photon's weight by the probability that an interaction would occur */
	newWeight = (*weight) * (1 - exp(-freePathsToExit));
		
	/* Increment global counter for reporting */
	detData[DetCurParams].detWeightAdjusted += ((*weight) - newWeight);
	
	/* Update weight */
	*weight = newWeight;

	
	/* Pick a free path to go based on a truncated exponential distribution */
	{
		/* Start with a random from the exponential distribution */
		PhgMathGetTotalFreePaths(&randFromExp);
		
		/* Truncate to desired range */
		{

			/* Now compute the free paths to go */
			/* Truncate to desired range */
			{
				/* Compute whole number part NOTE I AM IGNORING THE
				   POSSIBILITY OF INTEGER OVERFLOW
				*/
				wholeNum = (LbUsFourByte) (randFromExp/
					freePathsToExit);

				fpToGo = ((randFromExp/
			   		freePathsToExit) -
					wholeNum)
					* freePathsToExit;
					
				if (fpToGo > freePathsToExit) {
					#ifdef HEAVY_PHG_DEBUG
						/* send alert, but don't terminate simulation */
						ErAlert("Invalid calculation of fpToGo for forced interaction (detCylinder)--PLEASE REPORT", false);
					#endif
					/* this should be an infrequent numerical condition, set to boundary value */
					fpToGo = freePathsToExit;
				}
			}

		}	
	}				
	
	return (fpToGo);
}


/*********************************************************************************
*
*			Name:			detCylCompCentroid
*
*			Summary:		Compute the centroid.
*
*			Arguments:
*				PHG_TrackingPhoton	*photonPtr			- The photon.
*				LbUsFourByte 		interactionIndex	- The number of interactions.
*				double				depositedEnergy		- The amount of energy deposited.
*
*			Function return: None.
*
*********************************************************************************/
void detCylCompCentroid(PHG_TrackingPhoton *photonPtr, LbUsFourByte interactionIndex,
		double depositedEnergy)
{
	PHG_Position		centroidPos;						/* The computed centroid position */
	LbFourByte			curInteraction = interactionIndex;	/* LCV */

	
	/* Clear the centroid position */
	centroidPos.x_position = 0.0;
	centroidPos.y_position = 0.0;
	centroidPos.z_position = 0.0;


	/* Loop through interactions */
	while (curInteraction >= 0) {

		/* Only contribute to centroid if interaction was in an active layer */
		if (photonPtr->det_interactions[curInteraction].isActive) {
			centroidPos.x_position += (photonPtr->det_interactions[curInteraction].pos.x_position *
				photonPtr->det_interactions[curInteraction].energy_deposited);
			
			centroidPos.y_position += (photonPtr->det_interactions[curInteraction].pos.y_position *
				photonPtr->det_interactions[curInteraction].energy_deposited);
			
			centroidPos.z_position += (photonPtr->det_interactions[curInteraction].pos.z_position *
				photonPtr->det_interactions[curInteraction].energy_deposited);
		}
								
		curInteraction--;
	}

	/* Complete computation of centroidPos */
	centroidPos.x_position /= depositedEnergy;
	centroidPos.y_position /= depositedEnergy;
	centroidPos.z_position /= depositedEnergy;


	
	/*	We'll update the photon to reflect its status after passing through
		the detector. Note that it would not be valid to continue to track
		this photon through any remaining space.
	*/
	photonPtr->detLocation = centroidPos;
	photonPtr->energy = depositedEnergy;
}

/*********************************************************************************
*
*			Name:			DetCylinder
*
*			Summary:		Perform detection for cylindrical detectors
*
*			Arguments:
*				PHG_Decay			decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton *bluePhotons		- The blue photons detected.
*				LbUsFourByte 		numBluePhotons	- The number of blue photons.
*				PHG_TrackingPhoton *pinkPhotons		- The  pink photons detected.
*				LbUsFourByte		numPinkPhotons	- The number of pink photons.
*				DetectedPhotonsTy *detPhotonsPtr	- The accepted detected photons
*
*			Function return: None.
*
*********************************************************************************/
void DetCylinder(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		DetectedPhotonsTy *detPhotonsPtr)
{
	Boolean						firstTime;		/* Flag for first round of tracking */
	PHG_TrackingPhoton 			*photons;		/* Used to access both blue and pinks */
	LbUsFourByte				loopCount;		/* Current time through the loop once for pinks and once for blues */
	LbUsFourByte				pIndex;			/* Current photon */
	LbUsFourByte				numPhotons;		/* Used for both blue and pinks */
	LbFourByte					curInteraction;	/* Current interaction */
	double						energyDepositedFromInteraction;		/* Energy deposited in due to a scatter */
	double						depositedEnergy;					/* Energy deposited at current location */
	
	/* Clear the counters */
	detPhotonsPtr->NumDetectedBluePhotons = 0;
	detPhotonsPtr->NumDetectedPinkPhotons = 0;
	detData[DetCurParams].detDetectedBluePhotonIndex = 0;
	detData[DetCurParams].detDetectedPinkPhotonIndex = 0;
	#ifdef PHG_DEBUG
	detData[DetCurParams].DetDiscardBlue = true;
	detData[DetCurParams].DetDiscardPink = true;
	#endif
	
	/* Initialize photon pointer and loop control */
	photons = bluePhotons;
	numPhotons = numBluePhotons;
	loopCount = 1;

	do { /* Loop for blue and pink photons */
	curInteraction = -1;
	firstTime = true;
	energyDepositedFromInteraction=0.0;
	depositedEnergy = 0.0;
	
	/* Loop through all photons */
	for (pIndex = 0; pIndex < numPhotons; pIndex++) {

		/* Update our cylinders */
		detCylInitCylinders(0, 0);

		/* Let user modify and/or reject photons */
		if (DetUsrModPETPhotonsFPtr && 
				(*DetUsrModPETPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
					&(photons[pIndex])) == false) {
			
			/* They rejected it so go to next photon */
			continue;
		}
		
		/* Track the photon */
		if (DetCylTrack(decayPtr, &(photons[pIndex]))){
			/* We made it to here so save the detected photons */
			if (photons == bluePhotons) {
#ifdef PHG_DEBUG
detData[DetCurParams].DetDiscardBlue = false;
#endif
				detPhotonsPtr->DetectedTrkngBluePhotons[detPhotonsPtr->NumDetectedBluePhotons]
					= photons[pIndex];

				/* Increment the counters */
				detPhotonsPtr->NumDetectedBluePhotons++;
			}
			else {
#ifdef PHG_DEBUG
detData[DetCurParams].DetDiscardPink = false;
#endif
				detPhotonsPtr->DetectedTrkngPinkPhotons[detPhotonsPtr->NumDetectedPinkPhotons]
					= photons[pIndex];

				/* Increment the counters */
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
	
	/* Now setup for pink photons */
	photons = pinkPhotons;
	numPhotons = numPinkPhotons;
	loopCount++;
	
	} while (loopCount <= 2);


}

/*********************************************************************************
*
*			Name:			DetCylTrack
*
*			Summary:		Perform tracking for cylindrical detectors
*
*			Arguments:
*				PHG_Decay			*decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton 	*photonPtr			- The blue photons detected.
*
*			Function return: True if photon should be accepted.
*
*********************************************************************************/
Boolean DetCylTrack(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *photonPtr)
{
	PHG_Decay*					dummyPtr;		/* Removes compiler warning */
	Boolean						acceptPhoton = false;	/* Flag for acceptence of photon */
	Boolean						firstTime;		/* Flag for first round of tracking */
	PHG_Position				newPos;			/* Projected position */
	LbFourByte					curLayer=0;		/* Current layer */
	LbFourByte					curRing;		/* Current ring */
	LbFourByte					curInteraction;	/* Current interaction */
	LbFourByte					newLayer;		/* Potential next layer */
	double						attenuation;	/* Attenuation of current material */
	double						fpToGo;			/* Free paths to travel */
	double						rSquared;		/* Temp for position calculation */
	double						distance;		/* Distance variable, may be used multiple ways */
	double						energyDepositedFromInteraction;		/* Energy deposited in due to a scatter */
	double						depositedEnergy;					/* Energy deposited at current location */
	double						depositedActiveEnergy;
	double						weight;
	double						comptonToScatterProbability;	/* Ratio of prob of compton to prob of scatter */
	double						scatterProbability;				/* Probability of a scatter */
	double						interactionProbability;			/* Interaction-type probability */
	DetCylnRingTy				*ringPtr;						/* Current ring info */
	detCylEn_ActionTy			action;							/* What to do in response to photon travel */
	
	
	/* Avoid unused parameter compiler warning */
	dummyPtr = decayPtr;
	
	/* Clear the counters */
	curInteraction = -1;
	firstTime = true;
	energyDepositedFromInteraction=0.0;
	depositedEnergy = 0.0;
	
	/* Update our cylinders */
	detCylInitCylinders(0, 0);
			
	/* Because we use the weight so much we select it and use a local variable */
	weight = photonPtr->photon_current_weight;

	/* The target cylinder of the PHG and the inner-most radius of the
		detector may be different. Hence, we will first check to
		see if we are on the inner-most surface. If not, we will
		project to there
	*/
	{
		/* Compute distance of photon from origin */
		rSquared = PHGMATH_Square(photonPtr->location.x_position) +
			PHGMATH_Square(photonPtr->location.y_position);
			
		
		/* If not on the surface of the bounding cylinder, project onto it */
		if (rSquared < (PHGMATH_Square(detCylInBoundCyl.radius) -0.000000001))  {  
			
			/* Project photon to cylinder */
			if (CylPosProjectToCylinder(&(photonPtr->location),
					&(photonPtr->angle), &detCylInBoundCyl,
					&newPos, &distance) == false) {
			
				/*	Go to next photon */
				goto REJECT;
			}
			
			/* See if we fall outside of axial limits */
			if ((newPos.z_position >= detCylInBoundCyl.zMax) ||
					(newPos.z_position <= detCylInBoundCyl.zMin)) {
				
				/*	Go to next photon */
				goto REJECT;
			}
			
			/* We made it here so update our position on the cylinder */
			photonPtr->location = newPos;
			photonPtr->travel_distance += distance;
		}
		else {
			
			/* See if we fall outside of axial limits */
			if ((photonPtr->location.z_position >= detCylInBoundCyl.zMax) ||
					(photonPtr->location.z_position <= detCylInBoundCyl.zMin)) {
				
				/*	Go to next photon */
				goto REJECT;
			}
		}		
	}
	/* Determine which ring the photon is in */
	for (curRing = 0; 
			curRing < (LbFourByte)(DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings); 
			curRing++){
		if (photonPtr->location.z_position <=
				(DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].MaxZ)){
			
			break;
		}
	}
	#ifdef PHG_DEBUG
		/* Verify we are within this ring */
		if ((photonPtr->location.z_position <
				DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].MinZ) ||
				(curRing == (LbFourByte)(DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings))){
				
			PhgAbort("Invalid computation of ring (detCylnPET)", false);
		}
	#endif


	/*	Now track the photon through the block layers until it is absorbed
		or it escapes.  At each interaction point tally the energy deposited
		and record the position of the event (if the interaction occurs in an
		active layer (see below))
	*/
	/* Clear loop counters */
	depositedEnergy = 0.0;
	depositedActiveEnergy = 0.0;
	curInteraction = -1;
	curLayer = 0;
	
	/* Update our cylinders */
	detCylInitCylinders(curRing, curLayer);
	
	/* Specify short-cut to ring information */
	ringPtr = &DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing];
	
	/* Get free paths to travel */  
	if (DetRunTimeParams[DetCurParams].DoForcedInteraction) {
		fpToGo = DetCylGtTruncatedFreePaths(photonPtr, firstTime, curRing, curLayer, &weight);
	} else {
		PhgMathGetTotalFreePaths(&fpToGo);
	}

	/* Begin tracking loop */
	do {
	
		/* Get the attenuation of the crystal at the current energy */
		SubObjGetAttenuationInTomo(
			ringPtr->LayerInfo[curLayer].LayerMaterial,
			photonPtr->energy, &attenuation);
		
		/* Compute the distance to travel, considering
		   no boundary intersections
		*/
		 distance = fpToGo/attenuation;

		 /* Project the photon the selected distance */
		newLayer = curLayer;

		action = detCylProject(photonPtr, ringPtr, firstTime, &newPos, &curRing,
			&newLayer, &distance);
		
		/* See if photon crossed a layer */
		if (action == detCylEnAc_LayerCross) {

			/* Check layer crossing for boundaries. If outer crossing
				is outer most layer, it escapes. If inner crossing is
				inner most layer it escapes. Otherwise the loop
				continues.
			*/
			if (newLayer == (LbFourByte)(DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].NumLayers)) {
				action = detCylEnAc_Escape;
				break;
			}
			else if (newLayer < 0) {
				action = detCylEnAc_Escape;
				break;
			}
			else {
				/* Set new layer */
				curLayer = newLayer;
				
				/*  Update fpToGo to reflect the free paths traveled */
				fpToGo = fpToGo - (distance * attenuation);		
				
				/* Update our cylinders */
				detCylInitCylinders(curRing, curLayer);

			}
		}
		
		/* See if photon crossed an axial segment */
		if (action == detCylEnAc_AxialCross) {
		
			/* Check ring crossing for boundaries. If crossing
				is boundary ring, then discard it.
			*/
			if (curRing == (LbFourByte)(DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings)) {
				action = detCylEnAc_Escape;
				break;
			}
			
			/* If crossing out through first ring, discard it */
			if (curRing < 0) {
				action = detCylEnAc_Escape;
				break;
			}
			
			
			/* Process photon within this new axial ring */
			/* Determine which layer we are in within this ring */
			detCylFindLayer(&newPos, curRing, &curLayer);
				
			/* Specify short-cut to ring information */
			ringPtr = &DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing];
			
			/* Update fpToGo to reflect the free paths traveled */
			fpToGo = fpToGo - (distance * attenuation);	
			
			/* Update our cylinders */
			detCylInitCylinders(curRing, curLayer);
			
		}
		
		/* If we made it to here, we want to update the 
			position and distance traveled.
		*/
		photonPtr->location = newPos;
		photonPtr->travel_distance += distance;
		
		
		/* If photon should interact, pick a scatter angle and
			adjust the angle
		*/
		if (action == detCylEnAc_Interact) {
			/* Increment count of interactions */
			curInteraction += 1;
			photonPtr->num_det_interactions = curInteraction + 1;
			
			#ifdef PHG_DEBUG
			detData[DetCurParams].detCylCountInteractions++;
			#endif

			/* Save the current position */
			photonPtr->det_interactions[curInteraction].pos = newPos;
			
			/* Get probability of compton to probability scatter ratio */
			comptonToScatterProbability =
				SubObjGetProbComptToScatInTomo2(
				ringPtr->LayerInfo[curLayer].LayerMaterial,
				photonPtr->energy);
			
				
			/* Get probability of scatter */
			scatterProbability = SubObjGetProbScatterInTomo2(
				ringPtr->LayerInfo[curLayer].LayerMaterial,
				photonPtr->energy);
			
			/* Sample for probability of photo-electric absorbtion */
			interactionProbability = PhgMathGetRandomNumber();
			
			/* See if photon is absorbed */
			if ((interactionProbability > scatterProbability) ||
					((curInteraction+1) == MAX_DET_INTERACTIONS)) {
				
				/* Change action status from interaction to absorption */
				action = detCylEnAc_Absorb;
					
				/* Break out */
				break;
			}

			/* We made it to here, so model a scatter */
			{
				
				/* Store the current energy (it gets changed in EmisListDoComptonInteraction) */
				energyDepositedFromInteraction = photonPtr->energy;
				

				/* Determine which type of interaction to model, if it should be coherent
					but coherent modelling is turned off, then no interaction is modeled
				*/
				if (interactionProbability > (scatterProbability * comptonToScatterProbability	)) {
					if (PHG_IsModelCoherentInTomo()) {
						EmisListDoCoherent(photonPtr,
							ringPtr->LayerInfo[curLayer].LayerMaterial);
					}
					else {

						/*  Select a free paths for the next interaction */
						PhgMathGetTotalFreePaths(&fpToGo);
					
						/* Decrement interactions because in this case we decide not to do one */
						curInteraction--;

						#ifdef PHG_DEBUG
						detData[DetCurParams].detCylCountInteractions--;
						detData[DetCurParams].detCylCountCohInteractions++;							
						#endif

						continue; 
					}
				}
				else {

					/* Compute scatter angle and new energy */
					EmisListDoComptonInteraction(photonPtr);
					
					/* If we fall below our supported energy due to the scatter, we'll force an absorption */
					if (photonPtr->energy < PHG_MIN_PHOTON_ENERGY) {
					
						/* Increment number of photons absorbed */
						detData[DetCurParams].detTotPhotonsAbsorbed++;
						detData[DetCurParams].detTotForcedAbsorptions++;
						
						/* Change action status from interaction to absorption */
						action = detCylEnAc_Absorb;
						
						/* Break out */
						break;
						
					}
				}
				
				/* 	If we made it here it was a normal scatter.
					Compute the deposited energy
				*/
				energyDepositedFromInteraction -= photonPtr->energy;
				
				/*  Select a free paths for the next interaction */
				PhgMathGetTotalFreePaths(&fpToGo);	
			
			} /* End of Comment block - We made it here so model a scatter */
			
			/* Store the interaction energy */
			if (ringPtr->LayerInfo[curLayer].IsActive){
				
				/* Mark interaction as from active layer */
				photonPtr->det_interactions[curInteraction].isActive = true;

				/* Save the current energy */
				photonPtr->det_interactions[curInteraction].energy_deposited =
					energyDepositedFromInteraction;
				
				/* Sum the deposited energy */
				depositedEnergy += energyDepositedFromInteraction;

				/* Sum active layer energy */
				depositedActiveEnergy += energyDepositedFromInteraction;

			}
			else {
				photonPtr->det_interactions[curInteraction].isActive = false;
				photonPtr->det_interactions[curInteraction].energy_deposited = 0.0;
			}
			
			/* Update the position */
			photonPtr->location = newPos;					
		}
		
		/* Clear "first time" flag */
		firstTime = false;
			
	} while ((action != detCylEnAc_Detect) && (action != detCylEnAc_Discard) &&
		(action != detCylEnAc_Absorb) && (action != detCylEnAc_Escape));
	
	/* If photon was absorbed, update statistics */
	if (action == detCylEnAc_Absorb) {
				
		/* If we reached maximum interactions, increment the counter */
		if ((curInteraction+1) == MAX_DET_INTERACTIONS)
			detData[DetCurParams].detNumReachedMaxInteractions++;
			
		/* Save the current energy */
		if (ringPtr->LayerInfo[curLayer].IsActive){
			
			/* Mark interaction as from active layer */
			photonPtr->det_interactions[curInteraction].isActive = true;

			/* Save the current energy */
			photonPtr->det_interactions[curInteraction].energy_deposited =
				photonPtr->energy;
			
			/* Sum the deposited energy */
			depositedEnergy += photonPtr->det_interactions[curInteraction].energy_deposited;

			/* Sum active layer energy */
			depositedActiveEnergy += photonPtr->det_interactions[curInteraction].energy_deposited;
		}
		else {
			photonPtr->det_interactions[curInteraction].isActive = false;
			photonPtr->det_interactions[curInteraction].energy_deposited = 0.0;
		}
		/* Increment counter */
		detData[DetCurParams].detTotPhotonsAbsorbed++;
		
		/* Increment weight counter */
		detData[DetCurParams].detTotWtAbsorbed += (weight*photonPtr->decay_weight);
		
		/* Update absorbed array, as long as interaction count isn't too hight */
		if (curInteraction < MAX_DET_INTERACTIONS)
			detData[DetCurParams].detWeightAbsorbedBins[curInteraction+1] += (weight*photonPtr->decay_weight);
		
		/* If it is the first interaction, update statistics for 
			absorbtion on first interaction
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
	if (depositedEnergy != 0.0){	
		
		/* Increment counter */
		detData[DetCurParams].detTotPhotonsDepositingEnergy++;
		
		photonPtr->photon_current_weight = weight;

		detCylCompCentroid(photonPtr, curInteraction, depositedActiveEnergy);
		
		photonPtr->location.x_position = photonPtr->detLocation.x_position;
		photonPtr->location.y_position = photonPtr->detLocation.y_position;
		photonPtr->location.z_position = photonPtr->detLocation.z_position;
		
		
		/* Clear out the direction vector, it's no longer meaningful */
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

#ifdef PHG_DEBUG
/*********************************************************************************
*
*			Name:			detDumpCylInteractions
*
*			Summary:		Prints out interaction list for debugging purposes.
*
*			Arguments:
*				LbFourByte					numInteractions	- The number of interactions.
*				detCylInteractionInfoTy	interactions	- The interactions.
*
*			Function return: None.
*
*********************************************************************************/
void detDumpCylInteractions(LbFourByte numInteractions,
		detInteractionInfoTy *interactions)
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
*			Name:			detCylInitCylinders
*
*			Summary:		Initialize our cylinders based on current
*							ring and layer.
*
*			Arguments:
*				LbUsFourByte	curRing			- The current ring
*				LbUsFourByte	curLayer		- The current layer
*			Function return: None.
*
*********************************************************************************/
void detCylInitCylinders(LbUsFourByte curRing, LbUsFourByte curLayer)
{
	/*
	!!SV Robert, you are right about this routine the values should be set according to the incoming
	!!SV parameters and/or the overall min and max ring number.  It would be best to make this two routines
	!!SV since two types of initializations are being done and there is some redundancy here. I'll make a
	!!SV note of it and do it later.  There is no separate initialization routine for this module at this time
	!!SV so that would need to be added.
	*/
	
	/* Start by setting inner and outer radius to that of requested segment */ 
	detCylInnerRadius = DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].LayerInfo[curLayer].InnerRadius;
	detCylOuterRadius = DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].LayerInfo[curLayer].OuterRadius;
	detCylLayerDepth = detCylOuterRadius - detCylInnerRadius;
	
	/* Now initialize the cylinders */
	detCylInBoundCyl.radius = detCylInnerRadius;
	detCylInBoundCyl.zMin = DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[0].MinZ;
	detCylInBoundCyl.zMax = DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings-1].MaxZ;
	detCylInBoundCyl.centerX = 0.0;
	detCylInBoundCyl.centerY = 0.0;
	
	detCylOutBoundCyl.radius = detCylOuterRadius;
	detCylOutBoundCyl.zMin = DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[0].MinZ;
	detCylOutBoundCyl.zMax = DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings-1].MaxZ;
	detCylOutBoundCyl.centerX = 0.0;
	detCylOutBoundCyl.centerY = 0.0;
}

/*********************************************************************************
*
*			Name:			detCylProject
*
*			Summary:		Project a photon through the detector.
*							Truncate the projection based on the free paths to travel,
*							the inner/outer cylinder boundaries, the current
*							segment boundary.
*
*			Arguments:
*				PHG_TrackingPhoton				*photonPtr		- The photon.
*				DetCylnRingTy					*curRingInfo	- The current ring.
*				Boolean							firstTime		- Is this the first projection?
*				PHG_Position					*newPosPtr		- The new photon position.
*				LbFourByte						*curRingIndex	- The current ring index.
*				LbFourByte						*curLayerIndex	- The current layer index.
*				double							*distPtr		- The distance traveled.
*
*			Function return: The action to take based on the results of tracking.
*
*********************************************************************************/
detCylEn_ActionTy detCylProject(PHG_TrackingPhoton *photonPtr,
					DetCylnRingTy *curRingInfo, Boolean firstTime, 
					PHG_Position *newPosPtr, LbFourByte *curRingIndex,
					LbFourByte *curLayerPtr, double *distPtr)
{
	double				distToInner = 0.0;			/* Distance to inner cylinder */
	double				distToOuter = 0.0;			/* Distance to inner cylinder */
	double				distToWall = 0.0;			/* Distance to donut wall */
	double				distToWall1 = -1.0;			/* Distance to donut wall */
	double				distToWall2 = -1.0;			/* Distance to donut wall */
	detCylEn_ActionTy	action = detCylEnAc_Null;	/* What action will result? */
	
	
	do { /* Process Loop */
	
		
		/*	If this is the first projection, we know the photon is
			traveling outwards so we don't worry about intersections
			with the inner cylinder.
			
			However, if not, if the photon is going to hit the inner cylinder, it will
			do it before hitting the outer cylinder. So we will start by
			getting the distance to the inner cylinder. If it is non-
			negative, and not "equal" to zero,
			then the photon is traveling towards the inner
			cylinder.
		*/
		if ((firstTime == false) && ((distToInner = detCylFindInnerDist((&photonPtr->location),
				&(photonPtr->angle))) > 0)){
		
			/* Find intersection with donut walls */
			if (photonPtr->angle.cosine_z == 0.0) {
				/* If cosine_z is zero then the wall intersection will not happen */
				distToWall1 = distToWall2 = distToInner+1;
			}
			else {
				distToWall1 = ((curRingInfo->MaxZ) - photonPtr->location.z_position)/photonPtr->angle.cosine_z;
				distToWall2 = (curRingInfo->MinZ - photonPtr->location.z_position)/photonPtr->angle.cosine_z;

				/* Save distToWall as maximum positive (non-zero) of 1 & 2 */
				/* NOTE: if both are negative it is weeded out below */
				if ((distToWall1 > 0.0) && ((distToWall2 < 0.0) || (PhgMathRealNumAreEqual(distToWall2, 0.0, -5, 0, 0, 0) == true))) {
					distToWall = distToWall1;
				} else if ((distToWall2 > 0.0) && ((distToWall1 < 0.0) || (PhgMathRealNumAreEqual(distToWall1, 0.0, -5, 0, 0, 0) == true))) {
					distToWall = distToWall2;
				} else {
					#ifdef PHG_DEBUG
							PhgAbort("Faulty distances to z ring boundaries [1] (detCylProject).", true);
					#endif
					distToWall =  ((distToWall1 > distToWall2) 
							? distToWall1 : distToWall2);
				}
			}
							
			/* Pick the smallest positive solution (we know distToInner is > 0) */
			if ((distToWall > 0) && (distToWall < distToInner)) {
				
				/* Compare distance to wall with given distance to travel,
					determines if photon will interact
				 */
				if (distToWall > *distPtr) {
				
					/* Photon will interact at distance given */
					action = detCylEnAc_Interact;
					break;
				
				}
				else {
				
					/* Photon will interact at cell boundary */
					*distPtr = distToWall; 
					
					/* Inc/dec segment index depending on which wall was intersected */
					if (distToWall == distToWall1)
						*curRingIndex += 1;
					else
						*curRingIndex -= 1;
						
					/* Set our action */
					action = detCylEnAc_AxialCross;
					break;
				}
			}
			else {
				
				/* Compare distance to inner cylinder with given distance to travel,
					determines if photon will interact
				 */
				if (distToInner > *distPtr) {
				
					/* Photon will interact at given distance */
					action = detCylEnAc_Interact;
					break;
				
				}
				else {
					
					/* The photon will reach the inner cylinder surface
						without interacting. 
					*/
					*distPtr = distToInner;
					action = detCylEnAc_LayerCross;
					*curLayerPtr -= 1;
					break;
				}
			}
			
		}
		else {
		
			/* The photon is going to intersect with the outer cylinder */
			if (CylPosProjectToCylinder(&(photonPtr->location),
					&(photonPtr->angle), &detCylOutBoundCyl, newPosPtr, &distToOuter)
					== false) {

				distToOuter = LBDOUBLE_MAX;
			}
			
		
			/* Find intersection with donut wall */
			if (photonPtr->angle.cosine_z == 0.0) {
				/* If cosine_z is zero then the wall intersection will not happen */
				distToWall1 = distToWall2 = distToOuter+1;
			}
			else {
				distToWall1 = ((curRingInfo->MaxZ) - photonPtr->location.z_position)/photonPtr->angle.cosine_z;
				distToWall2 = (curRingInfo->MinZ - photonPtr->location.z_position)/photonPtr->angle.cosine_z;

				/* Save distToWall as minimum positive of 1 & 2 */
				/* NOTE: if both are negative it is weeded out below */
				if ((distToWall1 > 0.0) && ((distToWall2 < 0.0) || (PhgMathRealNumAreEqual(distToWall2, 0.0, -5, 0, 0, 0) == true))) {
					distToWall = distToWall1;
				} else if ((distToWall2 > 0.0) && ((distToWall1 < 0.0) || (PhgMathRealNumAreEqual(distToWall1, 0.0, -5, 0, 0, 0) == true))) {
					distToWall = distToWall2;
				} else {
					#ifdef PHG_DEBUG
							PhgAbort("Faulty distances to z ring boundaries [2] (detCylProject).", true);
					#endif
					distToWall =  ((distToWall1 > distToWall2) 
							? distToWall1 : distToWall2);
				}
			}

			/* If distToWall is less than distToOuter, that is the maximum distance */
			if ((distToWall > 0) && (distToWall < distToOuter)) {
				
				/* Compare distance to wall with given distance to travel,
					determines if photon will interact
				 */
				if (distToWall > *distPtr) {
				
					/* Photon will interact at given distance */
					action = detCylEnAc_Interact;
					break;
				
				}
				else {
				
					/* Photon will interact at wall */
					*distPtr = distToWall;
					
					/* Inc/dec segment index depending on which wall was intersected */
					if (distToWall == distToWall1)
						*curRingIndex += 1;
					else
						*curRingIndex -= 1;
						
					/* Set our action */
					action = detCylEnAc_AxialCross;
					break;
				}
			}
			else {
				
				/* Compare distance to inner cylinder with given distance to travel,
					determines if photon will interact
				 */
				if (distToOuter > *distPtr) {
				
					/* Photon will interact at given distance */
					action = detCylEnAc_Interact;
					break;
				
				}
				else {
					
					/* The photon will reach the outer cylinder surface
						without interacting. 
					*/
					*distPtr = distToOuter;
					action = detCylEnAc_LayerCross;
					*curLayerPtr += 1;
					
					break;
				}
			}
		}
				
	} while (false);
	
	#ifdef PHG_DEBUG
		/* Make sure we set some action */
		if (action == detCylEnAc_Null)
			PhgAbort("Failed to set an action in tracking through cylindrical detector (detCylProject).", true);
	#endif
	
	/* Project to new location */
	newPosPtr->x_position = photonPtr->location.x_position +
		(photonPtr->angle.cosine_x * *distPtr);
			
	newPosPtr->y_position = 	photonPtr->location.y_position +
		(photonPtr->angle.cosine_y * *distPtr);
			
	newPosPtr->z_position = 	photonPtr->location.z_position +
		(photonPtr->angle.cosine_z * *distPtr);
			
	return(action);
}

/*********************************************************************************
*
*			Name:			detCylCompFreePathsToExit
*
*			Summary:		Compute the free-paths necessary for the photon
*							to exit the detector cylinder.
*
*			Arguments:
*				PHG_TrackingPhoton				*photonPtr		- The photon.
*				Boolean							firstTime		- Is this the first projection?
*				LbFourByte						curRingIndex	- The current ring index.
*				LbFourByte						curLayerIndex	- The current layer index.
*				double							*fpToGo			- The free paths to exit.
*
*			Function return: None.
*
*********************************************************************************/
void detCylCompFreePathsToExit(PHG_TrackingPhoton *photonPtr,
					Boolean firstTime, 
					LbFourByte curRingIndex,
					LbFourByte curLayerIndex, double *fpToGo)
{
	double				distToInner = 0.0;			/* Distance to inner cylinder */
	double				distToOuter = 0.0;			/* Distance to inner cylinder */
	double				distToWall = 0.0;			/* Distance to donut wall */
	double				distToWall1 = -1.0;			/* Distance to donut wall */
	double				distToWall2 = -1.0;			/* Distance to donut wall */
	double				distance = 0.0;				/* The incremental distance */
	double				attenuation = 0.0;			/* The attenuation of the current material */
	detCylEn_ActionTy	action = detCylEnAc_Null;	/* What action will result? */
	DetCylnRingTy		*ringPtr;					/* The current ring */
	LbUsFourByte		startingRingIndex = curRingIndex;		/* The starting value used to reset module globals */
	LbUsFourByte		startingLayerIndex = curLayerIndex;		/* The starting value used to reset module globals */
	PHG_Position 		position = photonPtr->location;			/* Temp for computing new position */	
	PHG_Position 		tempPos;								/* Temp for computing new position */	
	 
	/* Clear loop variables */
	*fpToGo = 0.0;

	do { /* Tracking Loop */
		
		/* Set our current ring */
		ringPtr = &DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRingIndex];

		/* Get attenuation for current location */
		SubObjGetAttenuationInTomo(
			ringPtr->LayerInfo[curLayerIndex].LayerMaterial,
				photonPtr->energy, &attenuation);

		
		/*	If this is the first projection, we know the photon is
			traveling outwards so we don't worry about intersections
			with the inner cylinder.
			
			However, if not, if the photon is going to hit the inner cylinder, it will
			do it before hitting the outer cylinder. So we will start by
			getting the distance to the inner cylinder. If it is non-
			negative, and not "equal" to zero,
			then the photon is traveling towards the inner
			cylinder.
		*/
		if ((firstTime == false) && ((distToInner = detCylFindInnerDist((&position),
				&(photonPtr->angle))) > 0)){
		
			/* Find intersection with donut walls */
			if (photonPtr->angle.cosine_z == 0.0) {
				/* If cosine_z is zero then the wall intersection will not happen */
				distToWall1 = distToWall2 = distToOuter+1;
			}
			else {
				distToWall1 = ((ringPtr->MaxZ) - position.z_position)/photonPtr->angle.cosine_z;
				distToWall2 = (ringPtr->MinZ - position.z_position)/photonPtr->angle.cosine_z;
	
				/* Save distToWall as maximum positive (non-zero) of 1 & 2 */
				/* NOTE: if both are negative it is weeded out below */
				if ((distToWall1 > 0.0) && ((distToWall2 < 0.0) || (PhgMathRealNumAreEqual(distToWall2, 0.0, -5, 0, 0, 0) == true))) {
					distToWall = distToWall1;
				} else if ((distToWall2 > 0.0) && ((distToWall1 < 0.0) || (PhgMathRealNumAreEqual(distToWall1, 0.0, -5, 0, 0, 0) == true))) {
					distToWall = distToWall2;
				} else {
					#ifdef PHG_DEBUG
							PhgAbort("Faulty distances to z ring boundaries [3] (detCylCompFreePathsToExit).", true);
					#endif
					distToWall =  ((distToWall1 > distToWall2) 
							? distToWall1 : distToWall2);
				}
			}
							
			/* Pick the smallest positive solution (we know distToInner is > 0) */
			if ((distToWall > 0) && (distToWall < distToInner)) {
			
				/* Photon will interact at cell boundary */
				distance = distToWall;
				
				/* Inc/dec segment index depending on which wall was intersected */
				if (distToWall == distToWall1)
					curRingIndex += 1;
				else
					curRingIndex -= 1;
					
				/* Set our action */
				action = detCylEnAc_AxialCross;
			}
			else {

				/* The photon will reach the inner cylinder surface
					without interacting. 
				*/
				distance = distToInner;
				action = detCylEnAc_LayerCross;
				curLayerIndex -= 1;
			}
			
		}
		else {
		
			/* The photon is going to intersect with the outer cylinder */
			if (CylPosProjectToCylinder(&position,
					&(photonPtr->angle), &detCylOutBoundCyl, &tempPos, &distToOuter)
					== false) {

				distToOuter = LBDOUBLE_MAX;
			}
			
		
			/* Find intersection with donut wall */
			if (photonPtr->angle.cosine_z == 0.0) {
				/* If cosine_z is zero then the wall intersection will not happen */
				distToWall1 = distToWall2 = distToOuter+1;
			}
			else {
				distToWall1 = ((ringPtr->MaxZ) - position.z_position)/photonPtr->angle.cosine_z;
				distToWall2 = (ringPtr->MinZ - position.z_position)/photonPtr->angle.cosine_z;
				/* Save distToWall as minimum positive of 1 & 2 */
				/* NOTE: if both are negative it is weeded out below */
				if ((distToWall1 > 0.0) && ((distToWall2 < 0.0) || (PhgMathRealNumAreEqual(distToWall2, 0.0, -5, 0, 0, 0) == true))) {
					distToWall = distToWall1;
				} else if ((distToWall2 > 0.0) && ((distToWall1 < 0.0) || (PhgMathRealNumAreEqual(distToWall1, 0.0, -5, 0, 0, 0) == true))) {
					distToWall = distToWall2;
				} else {
					#ifdef PHG_DEBUG
							PhgAbort("Faulty distances to z ring boundaries [4](detCylCompFreePathsToExit).", true);
					#endif
					distToWall =  ((distToWall1 > distToWall2) 
							? distToWall1 : distToWall2);
				}
			}

			/* If distToWall is less than distToOuter, that is the minimum distance */ 
			if ((distToWall > 0) && (distToWall < distToOuter)) {
				
				/* Photon will interact at wall */
				distance = distToWall;
				
				/* Inc/dec segment index depending on which wall was intersected */
				if (distToWall == distToWall1)
					curRingIndex += 1;
				else
					curRingIndex -= 1;
					
				/* Set our action */
				action = detCylEnAc_AxialCross;
			}
			else {
					
				/* The photon will reach the outer cylinder surface
					without intersecting a boundary. 
				*/
				distance = distToOuter;
				action = detCylEnAc_LayerCross;
				curLayerIndex += 1;
			}
		}
				
		/* Update free paths to go based on current material */
		*fpToGo += attenuation * distance;

		/* See if photon crossed a layer */
		if (action == detCylEnAc_LayerCross) {

			/* Check layer crossing for boundaries. If outer crossing
				is outer layer, it is detected. If inner crossing is
				minimum layer it is discarded. Otherwise the loop
				continues.
			*/
			if (curLayerIndex == (LbFourByte)(ringPtr->NumLayers)) {
				action = detCylEnAc_Escape;
				break;
			}
			else if (curLayerIndex < 0) {
				action = detCylEnAc_Escape;
				break;
			}
			else {
				/* Update our cylinders */
				detCylInitCylinders(curRingIndex, curLayerIndex);
				
			}
		}
		
		/* See if photon crossed a ring */
		if (action == detCylEnAc_AxialCross) {
			/* Check segment crossing for boundaries. If crossing
				is boundary segment, then discard it.
			*/
			if (curRingIndex == (LbFourByte)(DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings)) {
				action = detCylEnAc_Escape;
				break;
			}
			
			/* Check for minimum boundary */
			if (curRingIndex < 0) {
				action = detCylEnAc_Escape;
				break;
			}
		}
		
		/* If still in the loop then update position and travel through next section */
		position.x_position = position.x_position +
			(photonPtr->angle.cosine_x * distance);
				
		position.y_position = 	position.y_position +
			(photonPtr->angle.cosine_y * distance);
				
		position.z_position = 	position.z_position +
			(photonPtr->angle.cosine_z * distance);
		
		/* If we crossed an axial boundary, layer and cylinder updating must follow
			updating the position so it will be done now if appropriate.
		*/
		if (action == detCylEnAc_AxialCross) {
					
			/* Determine which layer we are in within this ring */
			detCylFindLayer(&position, curRingIndex, &curLayerIndex);
			
			/* Update our cylinders */
			detCylInitCylinders(curRingIndex, curLayerIndex);
		}
	} while (true);
	
	/* Reset module globals */
	detCylInitCylinders(startingRingIndex, startingLayerIndex);
}

/*********************************************************************************
*
*			Name:			detCylFindInnerDist
*
*			Summary:		Determine the distance to the inner cylinder.
*
*			Arguments:
*				PHG_Position	*posPtr	- The position of the photon.
*				PHG_Direction	*dirPtr	- The direction of the photon.
*
*			Function return: Distance to inner cylinder, negative if
*								it doesn't intersect.
*
*********************************************************************************/
double detCylFindInnerDist(PHG_Position *posPtr,
				PHG_Direction *dirPtr)
{
	double		xCord = 0.0;				/* X coordinate normalized to center of cylinder */	
	double		yCord = 0.0;				/* Y coordinate normalized to center of cylinder */
	double		distToInner = -1.0;			/* Distance to inner cylinder */
	double		aa, bb, cc, minRoot, maxRoot;
	LbFourByte	numRoots;
	

	 do { /* Process Loop */
	 
		/* Compute 'a' for quadratic */
		aa = 1 - PHGMATH_Square(dirPtr->cosine_z);

		/* Comput b and c for quadratic */		
		xCord = posPtr->x_position - detCylInBoundCyl.centerX;
		yCord = posPtr->y_position - detCylInBoundCyl.centerY;

		bb = 2 * ((dirPtr->cosine_x * xCord) + (dirPtr->cosine_y * yCord));
		cc = PHGMATH_Square(xCord) + PHGMATH_Square(yCord) - PHGMATH_Square(detCylInBoundCyl.radius);

		/* Solve quadratic */
		numRoots = PhgMathSolveQuadratic(aa, bb, cc,
						&minRoot, &maxRoot);

		/* Check for no intersection with inner cylinder */				
		if (numRoots < 2) {
			break;
		}
	
		/* Check for negative-only intersection */				
		if ((maxRoot < 0) || (PhgMathRealNumAreEqual(maxRoot, 0.0, -5, 0, 0, 0) == true)) {
			break;
		}
	
		/* Distance to inner intersection is the minimum root from the quadratic */
		distToInner = minRoot;
		
		
		#ifdef PHG_DEBUG
			if (distToInner < 0.0)
				PhgAbort("Faulty distances to inner cylinder:  dist1 < 0 < dist2 (detCylFindInnerDist).", true);
		#endif
				
	} while (false);			

	return(distToInner);
}

/*********************************************************************************
*
*			Name:			detCylFindRing
*
*			Summary:		Figure out which axial ring the photon is in.
*							
*			Arguments:
*				double			zPos		- The z position of the photon.
*				LbUsFourByte	curLayer	- The current layer.
*
*			Function return: None.
*
*********************************************************************************/
LbUsFourByte detCylFindRing(double zPos)
{
	LbUsFourByte	ringIndex;		/* The photon's segment */
	
	/* Loop  segments to find the one that contains the photon */
	for (ringIndex = 0; ringIndex < DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings; ringIndex++) {
	
		if ((DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[ringIndex].MinZ <= zPos)
				&&
				(DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[ringIndex].MaxZ > zPos)) {
			break;
		}
	}
	
	#ifdef PHG_DEBUG
		/* Verify we found a segment */
		if (ringIndex == DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings) {
			sprintf(detCylErrStr, "Photon starting out on detector ring out of z-bounds (detCylFindRing)"
				"\nzPos = %3.5lf, zMin = %3.5lf, zMax = %3.5lf\n",
				zPos, DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[ringIndex].MinZ,
				DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[ringIndex].MaxZ);
				
			PhgAbort(detCylErrStr, true);
		}
		
	#endif
	
	
	return(ringIndex);
}

/*********************************************************************************
*
*			Name:			detCylFindLayer
*
*			Summary:		Figure out which layer the photon is in.
*							
*			Arguments:
*				PHG_Position	*pos		- The position of the photon.
*				LbUsFourByte	ringIndex	- The ring.
*				LbUsFourByte	*theLayer	- The found layer.
*			Function return: None.
*
*********************************************************************************/
void detCylFindLayer(PHG_Position *pos, LbFourByte rngIdx, LbFourByte *theLayer)
{
	LbUsFourByte	layIdx;		/* The photon's layer */
	
	/* Loop  segments to find the one that contains the photon */
	for (layIdx = 0; layIdx < DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[rngIdx].NumLayers; layIdx++) {
	
		/* Since we are starting from the inside, as soon as we are within a layer's outer diameter, 
			we are in that layer */
		if ((PHGMATH_SquareRoot(PHGMATH_Square(pos->x_position)+PHGMATH_Square(pos->y_position))) <
				DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[rngIdx].LayerInfo[layIdx].OuterRadius) {
			break;
		}
	
	}
	
	#ifdef PHG_DEBUG
		/* Verify we found a layer */
		if (layIdx == DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[rngIdx].NumLayers) {
			PhgAbort("Photon moved into a new ring and the current layer is invalid (detCylFindLayer)", true);
		}
	#endif
	
	*theLayer = layIdx;
}


/*********************************************************************************
*
*			Name:			DetCylInitCylinders
*
*			Summary:		Initialize our cylinders based on current
*							ring and layer.
*
*			Arguments:
*				LbFourByte		curRing			- The current ring
*				LbFourByte		curLayer		- The current layer
*
*			Function return: 	None.
*
*********************************************************************************/

void DetCylInitCylinders(void)/* ### */
{
	/* ### */
	detCylInitCylinders(0, 0);
}


/*********************************************************************************
*
*			Name:			DetCylProjectToDetector
*
*			Summary:		Project the photon to the detector edge.
*
*			Arguments:
*				PHG_TrackingPhoton 	*photonPtr	- The photon to track and detect.
*				LbFourByte			*ringNum	- The ring projected to.
*				LbFourByte			*detNum		- The detector projected to.
*
*			Function return: 	True if photon not rejected (still valid).
*
*********************************************************************************/

Boolean DetCylProjectToDetector(PHG_TrackingPhoton *photonPtr, 
								LbFourByte *ringNum, LbFourByte *detNum)

{
	Boolean			valid;					/* True if photon still valid */
	double			rSquared;				/* Temp for position calculation */
	PHG_Position	newPos;					/* Projected position */
	double			distance;				/* Projected distance */
	LbFourByte		curRing;				/* Ring index */
	
	
	/* The target cylinder of the PHG and the inner-most radius of the
		detector may be different. Hence, we will first check to
		see if we are on the inner-most surface. If not, we will
		project to there
	*/
	
	valid = true;
	
	/* Compute distance of photon from origin */
	rSquared = PHGMATH_Square(photonPtr->location.x_position) +
		PHGMATH_Square(photonPtr->location.y_position);
	
	/* If not on the surface of the bounding cylinder, project onto it */
	if (rSquared < (PHGMATH_Square(detCylInBoundCyl.radius) - 0.000000001)) {  
		
		/* Project photon to cylinder */
		if (CylPosProjectToCylinder(&(photonPtr->location),
				&(photonPtr->angle), &detCylInBoundCyl,
				&newPos, &distance) == false) {
		
			/* Couldn't project the photon */
			valid = false;
		}
		
		if (valid) {
			/* See if we fall outside of axial limits */
			if ((newPos.z_position >= detCylInBoundCyl.zMax) ||
					(newPos.z_position <= detCylInBoundCyl.zMin)) {
				
				/* Outside axial limits */
				valid = false;
			}
		}
		
		if (valid) {
			/* Update the photon position on the cylinder */
			photonPtr->location = newPos;
			photonPtr->travel_distance += distance;
		}
	}
	else {
		/* Radially beyond inside detector edge already */
		
		/* See if we fall outside of axial limits */
		if ((photonPtr->location.z_position >= detCylInBoundCyl.zMax) ||
				(photonPtr->location.z_position <= detCylInBoundCyl.zMin)) {
			
			/* Outside axial limits */
			valid = false;
		}
	}
	
	if (valid) {
		/* Determine which ring the photon is in */
		*ringNum = DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings;	/* invalid number */
		*detNum = 0;	/* always */
		for (curRing = 0; curRing < (LbFourByte)(DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings); curRing++) {
			if (photonPtr->location.z_position <=
					(DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].MaxZ)){
				*ringNum = curRing;
				break;
			}
		}
		#ifdef PHG_DEBUG
			/* Validate the ring number */
			if ((*ringNum < 0) || 
						(*ringNum >= (LbFourByte)(DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings))) {
				
				PhgAbort("Invalid ring number (DetCylProjectToDetector)", false);
			}
			
			/* Verify we are within this ring */
			if ((photonPtr->location.z_position <
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[*ringNum].MinZ) ||
					(*ringNum == (LbFourByte)(DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings))){
				
				PhgAbort("Invalid computation of ring (DetCylProjectToDetector)", false);
			}
		#endif
	}
	
	return(valid);
}


/*********************************************************************************
*
*			Name:			DetCylFindNextInteraction
*
*			Summary:		Find and move the photon to the next interaction point.
*
*			Arguments:
*				PHG_TrackingPhoton 		*photonPtr		- The photon to track and detect.
*				Boolean					firstTime		- Optimizer for tracking.
*				LbFourByte				*curRingNum		- The current ring number.
*				LbFourByte				*curDetNum		- The current detector number.
*				detEn_ActionTy			*actionType		- The type of interaction that occurred.
*				double					*fpToGo			- Remaining free paths to travel.
*				LbUsFourByte			*detMaterial	- The interaction material of the detector.
*				Boolean					*isActive		- Active status of the interaction material.
*
*			Function return:	None.
*
*********************************************************************************/

void DetCylFindNextInteraction(PHG_TrackingPhoton *photonPtr, Boolean firstTime, 
								LbFourByte *curRingNum, LbFourByte *curDetNum, 
								detEn_ActionTy *actionType, double *fpToGo, 
								LbUsFourByte *detMaterial, Boolean *isActive)

{
	LbFourByte					curRing;		/* Current ring */
	LbFourByte					curLayer;		/* Current layer */
	DetCylnRingTy				*ringPtr;		/* Current ring info */
	double						attenuation;	/* Attenuation of current material */
	double						distance;		/* Distance variable, may be used multiple ways */
	LbFourByte					newLayer;		/* Potential next layer */
	detCylEn_ActionTy			action;			/* What to do in response to photon travel */
	detEn_ActionTy				detAction;		/* Photon action as detector type */
	PHG_Position				newPos;			/* Projected position */
	
	
	/* Set the ring the photon is in */
	curRing = *curRingNum;
	curLayer = *curDetNum;
	
	
	/*	Now track the photon through the cylinder layers until an interaction occurs.
		At the interaction point record the position and distance of the event.
	*/
	
	/* Update our cylinders */
	detCylInitCylinders(curRing, curLayer);
	
	/* Specify short-cut to ring information */
	ringPtr = &DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing];
	
	/* Get the attenuation of the crystal at the current energy */
	SubObjGetAttenuationInTomo(
		ringPtr->LayerInfo[curLayer].LayerMaterial,
		photonPtr->energy, &attenuation);
	
	/* Compute the distance to travel, considering
	   no boundary intersections */
	distance = *fpToGo / attenuation;
	
	/* Project the photon the selected distance */
	newLayer = curLayer;
	action = detCylProject(photonPtr, ringPtr, firstTime, &newPos, &curRing,
							&newLayer, &distance);
	
	/* Respond to different actions */
	if (action == detCylEnAc_LayerCross) {
		/* Photon crossed a layer */
		
		/* Check layer crossing for boundaries. If outer crossing
			is outer most layer, it escapes. If inner crossing is
			inner most layer it escapes. Otherwise the loop
			continues.
		*/
		if (newLayer == (LbFourByte)(DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].NumLayers)) {
			detAction = detEnAc_Discard;
		}
		else if (newLayer < 0) {
			detAction = detEnAc_Discard;
		}
		else {
			detAction = detEnAc_LayerCross;
			
			/* Set new layer */
			curLayer = newLayer;
			
			/* Update fpToGo to reflect the free paths traveled */
			*fpToGo = *fpToGo - (distance * attenuation);		
			
			/* Update our cylinders */
			detCylInitCylinders(curRing, curLayer);
		}
	}
	else if (action == detCylEnAc_AxialCross) {
		/* Photon crossed an axial segment */
		
		/* Check ring crossing for boundaries */
		if (curRing == (LbFourByte)(DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings)) {
			/* If crossing is boundary ring, then discard it */
			detAction = detEnAc_Discard;
		}
		else if (curRing < 0) {
			/* If crossing out through first ring, discard it */
			detAction = detEnAc_Discard;
		}
		else {
			detAction = detEnAc_AxialCross;
			
			/* Determine which layer we are now in within this ring */
			detCylFindLayer(&newPos, curRing, &curLayer);
			
			/* Update fpToGo to reflect the free paths traveled */
			*fpToGo = *fpToGo - (distance * attenuation);	
			
			/* Update our cylinders */
			detCylInitCylinders(curRing, curLayer);
		}
	}
	else {
		/* Only one other action is possible at this point */
		detAction = detEnAc_Interact;
	}
	
	/* Update the detector position */
	*curRingNum = curRing;
	*curDetNum = curLayer;
	
	/* Update the position and distance traveled */
	photonPtr->location = newPos;
	photonPtr->travel_distance += distance;
	
	/* Record the action, detector material, and layer active status */
	*actionType = detAction;
	if (detAction != detEnAc_Discard) {
		ringPtr = &DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing];
		*detMaterial = ringPtr->LayerInfo[curLayer].LayerMaterial;
		*isActive = ringPtr->LayerInfo[curLayer].IsActive;
	}
}
