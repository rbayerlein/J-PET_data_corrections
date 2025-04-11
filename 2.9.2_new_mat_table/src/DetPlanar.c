/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1996-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		DetPlanar.c
*			Revision Number:	2.2
*			Date last revised:	2 January 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	5 December 1996
*
*			Module Overview:	Simulates planar detector functionality.
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
**********************************************************************************
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
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		July 19, 2004
*
*			Revision description:	Added time-of-flight blurring for PET.
*
*********************************************************************************/

#define DETECTOR_PLANAR


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


#ifdef PHG_DEBUG
static char			detPnrErrStr[1024];				/* Storage for creating error strings */
#endif
static PHG_Decay	*detPlnrLastDecayPtr = NULL;	/* Current decay when initializing photons */

void	detPlnrComputeFreePathsToExit(PHG_Position pos,
			PHG_Direction dir, double energy,
			double *fpPtr);
Boolean	detTrackPlanar(PHG_TrackingPhoton *photonPtr);
void 			detDHDoRandomPos(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					DetectedPhotonsTy *detPhotonsPtr);
void			detDHDoOptimalPos(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					DetectedPhotonsTy *detPhotonsPtr);			

void			detPlnrProjectToNewPos(PHG_Position pos, 
					PHG_Direction dir, double distance, double frontOfLayer,
					double backOfLayer, PHG_Position *newPos, double *distTraveled,
					Boolean *exitXFront, Boolean *exitXBack, Boolean *exitY, Boolean *exitZ);

void			detPlnrCompCentroid(PHG_TrackingPhoton *photonPtr, LbUsFourByte interactionIndex,
					double depositedEnergy);
LbUsFourByte	detGtDetectableAngle(PHG_Direction bDir, PHG_Position bPos,
					PHG_Direction pDir, PHG_Position pPos, double *detPosition);
void 			detDualHeaded(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					DetectedPhotonsTy *detPhotonsPtr);
void			detDumpPlnrInteractions(LbFourByte numInteractions,
					detInteractionInfoTy *interactions);


/*********************************************************************************
*
*			Name:			DetPlanarSPECT
*
*			Summary:		Perform planar detection.
*
*			Arguments:
*				PHG_Decay			decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton *photons			- The blue photons detected.
*				LbUsFourByte 		numPhotons	- The number of blue photons.
*				DetectedPhotonsTy *detPhotonsPtr	- The accepted detected photons
*
*			Function return: None.
*
*********************************************************************************/
void DetPlanarSPECT(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
		DetectedPhotonsTy *detPhotonsPtr)
{
	LbUsFourByte				pIndex;					/* Index for current photon */
	double						detPosition;
		
	/* Clear the counters */
	detPhotonsPtr->NumDetectedBluePhotons = 0;
	detPhotonsPtr->NumDetectedPinkPhotons = 0;
	detData[DetCurParams].detDetectedBluePhotonIndex = 0;
	detData[DetCurParams].detDetectedPinkPhotonIndex = 0;


	do { /* Process Loop */
	
		/* Compute the detector position, either incremental or continuous (If collimation is done on the fly, then 
			the detector position is already established 
		*/
		if (PHG_IsCollimateOnTheFly() == false){
			if (DetGtNumViews() > 0) {
				detPosition = 	DetRunTimeParams[DetCurParams].PlanarDetector.MinAngle +
					((detData[DetCurParams].detPlnrDelta/2) + (floor(PhgMathGetRandomNumber() * DetGtNumViews()) * detData[DetCurParams].detPlnrDelta));
			}
			else {
				detPosition = 	PhgMathGetRandomNumber() * PHGMATH_2PI;
			}
		}
		
		/* Loop through all photons */
		for (pIndex = 0; pIndex < numPhotons; pIndex++) {

			/* If we've requested a fixed direction and there is no collimator on then
				we need to specify the detector angle
			*/
			if (PHG_IsCollimateOnTheFly() == false) {
				/* Pick a random detector position */
				photons[pIndex].detectorAngle = detPosition;
			}	

			/* Let user modify and/or reject photons */
			if (DetUsrModSPECTPhotonsFPtr && 
					(*DetUsrModSPECTPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
						&(photons[pIndex])) == false) {
				
				/* They rejected it so go to next photon */
				continue;
			}
			
			/* Track through the planar detector */
			if (detTrackPlanar(&photons[pIndex])) {
			

				/* Blur the energy  if requested */
				if (DET_DoEnergyBlur() == true) {
							
					photons[pIndex].energy = DetGaussEnergyBlur(photons[pIndex].energy);				

				}

				/* Save the photon */
				detPhotonsPtr->DetectedTrkngBluePhotons[detPhotonsPtr->NumDetectedBluePhotons]
					= photons[pIndex];

				/* Increment the counters */
				detPhotonsPtr->NumDetectedBluePhotons++;

				/* Update our detected photon block if doing history file */
				if  (DET_IsDoHistory()) {
					DetUpdateDetectedPhotonBlock(&photons[pIndex]);
				}
			}
		}

	} while (false);
}

/*********************************************************************************
*
*			Name:			DetDualHeaded
*
*			Summary:		Perform detection for dual headed detectors
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
void DetDualHeaded(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		DetectedPhotonsTy *detPhotonsPtr)
{

	/* Clear the counters */
	detPhotonsPtr->NumDetectedBluePhotons = 0;
	detPhotonsPtr->NumDetectedPinkPhotons = 0;
	detData[DetCurParams].detDetectedBluePhotonIndex = 0;
	detData[DetCurParams].detDetectedPinkPhotonIndex = 0;
	
	detDHDoRandomPos(decayPtr, bluePhotons, numBluePhotons, pinkPhotons, numPinkPhotons, detPhotonsPtr);
}

/*********************************************************************************
*
*			Name:			DetGetRandomDetPosition
*
*			Summary:		Compute a random detector angle
*
*			Arguments:
*				double				*detPos		- The position of the detector.
*
*			Function return: None.
*
*********************************************************************************/
void DetGetRandomDetPosition(double *detPos)
{
	/* Pick a random detector position */
	if (DetGtNumViews() > 0) {
		*detPos = 	DetRunTimeParams[DetCurParams].PlanarDetector.MinAngle +
			((detData[DetCurParams].detPlnrDelta/2) + (floor(PhgMathGetRandomNumber() * DetGtNumViews()) * detData[DetCurParams].detPlnrDelta));
	}
	else {
		*detPos = PHGMATH_PI * PhgMathGetRandomNumber();
	}
}

/*********************************************************************************
*
*			Name:			DetGetDetAngle
*
*			Summary:		Compute a random detector angle
*
*			Arguments:
*				PHG_Position		*location		- The location of the photon.
*				PHG_Direction		*angle			- The angle of travel.
*				double				dPos			- Detector position.
*				double				*detAngle		- The angle of the detector.
*
*			Function return: None.
*
*********************************************************************************/
void DetGetDetAngle(PHG_Position *location, PHG_Direction *angle,
		double dPos, double *detAngle)
{


	double				bBeta;
	double				difference;
	double				distance;
	PHG_Position		bigPos;
	
	/* Project blue photon to "big" cylinder */
	(void) CylPosProjectToCylinder(location,
		angle, &detData[DetCurParams].detPlnrBigCylinder, &bigPos, &distance);
			
	/* Compute polar coordinate angle  */
	bBeta = atan2(bigPos.y_position,bigPos.x_position);
		
	/* Determine which head the photon will hit */
	if (dPos-bBeta < -PHGMATH_PI) {
		difference = (dPos-bBeta) + (PHGMATH_2PI);
	}
	else if (dPos-bBeta > PHGMATH_PI){
		difference = (dPos-bBeta) - (PHGMATH_2PI);
	}
	else {
		difference = (dPos-bBeta);
	}
	
	if (fabs(difference) > (PHGMATH_PI_DIV2)) {
		if (dPos < PHGMATH_PI) {
			*detAngle = dPos + PHGMATH_PI;
		}
		else {
			*detAngle = dPos - PHGMATH_PI;
		}
	}
	else {
		*detAngle = dPos;
	}
}

/*********************************************************************************
*
*			Name:			detDHDoRandomPos
*
*			Summary:		Perform detection for dual headed detectors using
*							random detector head position
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
void detDHDoRandomPos(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		DetectedPhotonsTy *detPhotonsPtr)
{

	double				dPos;
	LbUsFourByte		b,p;
	
		
	/* Clear the counters */
	detPhotonsPtr->NumDetectedBluePhotons = 0;
	detPhotonsPtr->NumDetectedPinkPhotons = 0;
	detData[DetCurParams].detDetectedBluePhotonIndex = 0;
	detData[DetCurParams].detDetectedPinkPhotonIndex = 0;

	do { /* Process Loop */
	
	/* Compute the detector position, either incremental or continuous
		if collimation is done on the fly, then the detector position
		is already established
	 */
	if (PHG_IsCollimateOnTheFly() == false) {
		if (DetGtNumViews() > 0) {
			dPos = 	DetRunTimeParams[DetCurParams].PlanarDetector.MinAngle +
				((detData[DetCurParams].detPlnrDelta/2) + (floor(PhgMathGetRandomNumber() * DetGtNumViews()) * detData[DetCurParams].detPlnrDelta));
		}
		else {
			dPos = 	PhgMathGetRandomNumber() * PHGMATH_2PI;
		}
	}
	
	/* Loop through all blue photons */
	for (b = 0; b < numBluePhotons; b++) {
		
		/* If collimation has not been done then we need the detector position */
		if (PHG_IsCollimateOnTheFly() == false) {
		
			DetGetDetAngle(&bluePhotons[b].location, &bluePhotons[b].angle,
				dPos, &bluePhotons[b].detectorAngle);
		}
				
		/* Let user modify and/or reject photons */
		if (DetUsrModPETPhotonsFPtr && 
				(*DetUsrModPETPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
					&(bluePhotons[b])) == false) {
			
			/* They rejected it so go to next photon */
			continue;
		}
		
		/* Track the photon through the planar detector */
		if (detTrackPlanar(&bluePhotons[b]) == false)
			continue;	

		/* Blur the energy  if requested */
		if (DET_DoEnergyBlur() == true) {
					
			bluePhotons[b].energy = DetGaussEnergyBlur(bluePhotons[b].energy);				

		}
				
		/* Blur the energy  if requested */
		if ( DET_DoTofBlur() == true ) {
					
			bluePhotons[b].travel_distance = DetGaussTimeBlur(bluePhotons[b].travel_distance);				

		}
				
		/* We made it to here so save the detected photons */
		{
			
			detPhotonsPtr->DetectedTrkngBluePhotons[detPhotonsPtr->NumDetectedBluePhotons]
				= bluePhotons[b];

			/* Increment the counters */
			detPhotonsPtr->NumDetectedBluePhotons++;

			/* Update our detected photon block if doing history file */
			if  (DET_IsDoHistory()) {
				DetUpdateDetectedPhotonBlock(&bluePhotons[b]);
			}
		}
	}
	
	/* if no blue photons were found and this is a coincidence-only
	simulation, break out of the loop */
	if ( PHG_IsPETCoincidencesOnly() && (detPhotonsPtr->NumDetectedBluePhotons == 0) ) {
		
		break;
		
	}
	
	/* Loop through all pink photons */
	for (p = 0; p < numPinkPhotons; p++) {

		
		/* If collimation has not been done then we need the detector position */
		if (PHG_IsCollimateOnTheFly() == false) {
		
			DetGetDetAngle(&pinkPhotons[p].location, &pinkPhotons[p].angle,
				dPos, &pinkPhotons[p].detectorAngle);
		}

		/* Let user modify and/or reject photons */
		if (DetUsrModPETPhotonsFPtr && 
				(*DetUsrModPETPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
					&(pinkPhotons[p])) == false) {
			
			/* They rejected it so go to next pink */
			continue;
		}

		/* Track the photon through the planar detector */
		if (detTrackPlanar(&pinkPhotons[p]) == false)
			continue;	

		/* Blur the energy  if requested */
		if (DET_DoEnergyBlur() == true) {
					
			pinkPhotons[p].energy = DetGaussEnergyBlur(pinkPhotons[p].energy);				

		}
					
		/* Blur time if it was requested */
		if ( DET_DoTofBlur() == true ) {
		
			pinkPhotons[p].travel_distance = DetGaussTimeBlur(pinkPhotons[p].travel_distance);	
						
		}

		/* We made it to here so save the detected photons */
		{

			detPhotonsPtr->DetectedTrkngPinkPhotons[detPhotonsPtr->NumDetectedPinkPhotons]
				= pinkPhotons[p];

			/* Increment the counters */
			detPhotonsPtr->NumDetectedPinkPhotons++;

			/* Update our detected photon block if doing history file */
			if  (DET_IsDoHistory()) {
				DetUpdateDetectedPhotonBlock(&pinkPhotons[p]);
			}
		}
	}
	} while (false);
}

/*********************************************************************************
*
*			Name:			detPlnrCompCentroid
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
void detPlnrCompCentroid(PHG_TrackingPhoton *photonPtr, LbUsFourByte interactionIndex,
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
	photonPtr->transaxialPosition = centroidPos.y_position;
	photonPtr->energy = depositedEnergy;
}
/*********************************************************************************
*
*			Name:			detPlnrProjectToNewPos
*
*			Summary:		Project photon proposed distance, truncating to any
*							layer boundaries and flagging exit.
*
*			Arguments:
*				PHG_Position		pos				- The photon's position.
*				PHG_Direction		dir				- The photon's direction.
*				double				distance		- Proposed distance to travel.
*				double				frontOfLayer	- Proposed distance to travel.
*				double				backOfLayer,	- Proposed distance to travel.
*				PHG_Position		*newPos			- The photon's new position.
*				double				*distTraveled	- The distance traveled.
*				Boolean				*exitXFront		- Exits in the X direction.
*				Boolean				*exitXBack		- Exits in the X direction.
*				Boolean				*exitY			- Exits in the Y direction.
*				Boolean				*exitZ			- Exits in the Z direction.
*
*			Function return: None.
*
*********************************************************************************/
void detPlnrProjectToNewPos(PHG_Position pos, 
		PHG_Direction dir, double distance, double frontOfLayer,
		double backOfLayer, PHG_Position *newPos, double *distTraveled,
		Boolean *exitXFront, Boolean *exitXBack, Boolean *exitY, Boolean *exitZ)
{
	double			distX;				/* Distance to travel with respect to X axis */
	double			distY;				/* Distance to travel with respect to Y axis */
	double			distZ;				/* Distance to travel with respect to Z axis */
	
	/* Clear variables */
	*exitXFront = *exitXBack = *exitY = *exitZ = false;
	*distTraveled = 0.0;
	
	/*	Check which direction we are going with respect to Y and compute
		distance via that projection
	*/
	if (dir.cosine_y > 0.0){
		distY = (detData[DetCurParams].detPlnTransLimit - pos.y_position)/dir.cosine_y;
	}
	else if (dir.cosine_y < 0.0){
		distY = (-detData[DetCurParams].detPlnTransLimit - pos.y_position)/dir.cosine_y;
	}
	else {
		distY = MAXFLOAT;
	}
		
	/*	Check which direction we are going with respect to Z and compute
		distance via that projection
	*/
	if (dir.cosine_z > 0.0){
		distZ = (detData[DetCurParams].detInBoundCyl.zMax - pos.z_position)/dir.cosine_z;
	}
	else if (dir.cosine_z < 0.0) {
		distZ = (detData[DetCurParams].detInBoundCyl.zMin - pos.z_position)/dir.cosine_z;
	}
	else {
		distZ = MAXFLOAT;
	}
	
	if (dir.cosine_x > 0.0) {
		distX = (backOfLayer
			- pos.x_position)/dir.cosine_x;
	}
	else if (dir.cosine_x < 0.0) {
		distX = (frontOfLayer - pos.x_position)/dir.cosine_x;
	}
	else {
		distX = MAXFLOAT;
	}
	/* Now choose the minimum value */
	if ((distX < distY) && (distX < distZ) && (distX < distance)){
		
		*distTraveled = distX;
		if (dir.cosine_x > 0) {
			*exitXBack = true;
		}
		else {
			*exitXFront = true;
		}
	}
	else if ((distY < distZ) && (distY < distance)) {
		*distTraveled = distY;
		*exitY = true;
	}
	else if (distZ < distance){
		*distTraveled = distZ;
		*exitZ = true;
	}
	else {
		*distTraveled = distance;
	}	
		
	/* Project to position, do this even if the photon will exit so that it's
		final position is at the edge of the detector
	*/
	newPos->x_position = pos.x_position + (*distTraveled * dir.cosine_x);
	newPos->y_position = pos.y_position + (*distTraveled * dir.cosine_y);
	newPos->z_position = pos.z_position + (*distTraveled * dir.cosine_z);
}


/*********************************************************************************
*
*			Name:			detPlnrComputeFreePathsToExit
*
*			Summary:		Compute the free paths "necessary" for the
*							photon to exit the planar detector.
*
*			Arguments:
*				PHG_Position		pos			- The photon's position.
*				PHG_Direction		dir			- The photon's direction.
*				double				energy		- The photon's energy
*				double				*fpToGo		- The computed free paths
*
*			Function return: None.
*
*********************************************************************************/
void detPlnrComputeFreePathsToExit(PHG_Position pos,
		PHG_Direction dir, double energy,
		double *fpToGo)
{
	Boolean			exit;				/* Loop control variable */
	double			distance;			/* Distance to travel, used more than once */
	double			distanceTracked;	/* Distance already accounted for */
	double			attenuation;		/* Attenuation of crystal at current energy */
	double			distX;				/* Distance to travel with respect to X axis */
	double			distY;				/* Distance to travel with respect to Y axis */
	double			distZ;				/* Distance to travel with respect to Z axis */
	LbUsFourByte	curLayer = 0;		/* Index for current layer */

	/* Clear counters and accumulators */
	*fpToGo = 0.0;
	distanceTracked = 0.0;
	distance = 0.0;
	exit = false;
					
	/*	Check which direction we are going with respect to Y and compute
		distance via that projection
	*/
	if (dir.cosine_y > 0.0){
		distY = (detData[DetCurParams].detPlnTransLimit - pos.y_position)/dir.cosine_y;
	}
	else if (dir.cosine_y < 0.0){
		distY = (-detData[DetCurParams].detPlnTransLimit - pos.y_position)/dir.cosine_y;
	}
	else {
		distY = MAXFLOAT;
	}
		
	/*	Check which direction we are going with respect to Z and compute
		distance via that projection
	*/
	if (dir.cosine_z > 0.0){
		distZ = (detData[DetCurParams].detInBoundCyl.zMax - pos.z_position)/dir.cosine_z;
	}
	else if (dir.cosine_z < 0.0) {
		distZ = (detData[DetCurParams].detInBoundCyl.zMin - pos.z_position)/dir.cosine_z;
	}
	else {
		distZ = MAXFLOAT;
	}
		
	do { /* Loop through each layer */
			
			/* Get the attenuation of the crystal at the current energy */
			SubObjGetAttenuationInTomo(
				DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerMaterial,
				energy, &attenuation);
	
			/* Compute the distance to the back of the block
				Note this is called as the photon enters the detector for
				forced-interaction, the photon is sitting at x = 0
			*/
			distX = (DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerDepth
				/dir.cosine_x) + distanceTracked;
							
			/* Now choose the minimum value */
			if ((distX < distY) && (distX < distZ)){
				distance = distX;
			}
			else if (distY < distZ) {
				distance = distY;
				exit = true;
			}
			else {
				distance = distZ;
				exit = true;
			}
			/* Compute the free paths to go that distance */
			*fpToGo += ((distance - distanceTracked) * attenuation);
			
			
			#ifdef PHG_DEBUG
			
				/* Verify we have computed a positive free paths to exit */
				if (*fpToGo <= 0.0) {
					sprintf(detPnrErrStr, "\nInvalid computation of free paths to exit\n"
						"\t*fpPtr = %3.3f\tdistance = %3.3f\t distanceTracked = %3.3f\t"
						"attenuation = %3.3f\n"
						"distX = %3.3f\tdistY = %3.3f\tdistZ=%3.3f\n cosX = %3.3f\tcosY = %3.3f\t cosZ = %3.3f\n",
						*fpToGo, distance, distanceTracked, attenuation,
						distX, distY, distZ,
						dir.cosine_x, dir.cosine_y, dir.cosine_z);
					PhgAbort(detPnrErrStr, true);
				}
			#endif

			/* Updated distance tracked */
			distanceTracked = distance;
			
			/* Go to next layer */
			curLayer++;
			
	} while ((curLayer < DetRunTimeParams[DetCurParams].PlanarDetector.NumLayers) && (exit == false));
}
/*********************************************************************************
*
*			Name:			detTrackPlanar
*
*			Summary:		Determine if photon is detected by planar detector.
*
*			Arguments:
*				PHG_TrackingPhoton *photons			- The blue photons detected.
*				LbUsFourByte 		numPhotons	- The number of blue photons.
*				DetectedPhotonsTy *detPhotonsPtr	- The accepted detected photons
*
*			Function return: True if detection occurs.
*
*********************************************************************************/
Boolean detTrackPlanar(PHG_TrackingPhoton *photonPtr)
{
	Boolean 					photonDetected = false;	/* Successful detection or not? */
	Boolean						layerIsActive;			/* Flags current layer as active or not */
	Boolean						exitXFront;				/* Flags exit through front of layer */
	Boolean						exitXBack;				/* Flags exit through back of layer */
	Boolean						exitY;					/* Flags exit through y limit */
	Boolean						exitZ;					/* Flags exit through z limit */
	LbFourByte					curLayer = 0;			/* Index for current layer */
	LbUsFourByte				wholeNum;				/* Whole number temp variable */
	LbFourByte					curInteraction = -1;	/* Index for current interaction */
	double						t;						/* Projection distance */
	PHG_Position				pos;					/* The photon position */
	PHG_Position				newPos;					/* The new position after projection */
	PHG_Direction				dir;					/* The photon direction of travel */
	double						distance;				/* Distance to travel, used more than once */
	double						distTraveled;			/* Truncated distance for boundaries */
	double						energyDepositedFromInteraction=0.0;			/* Energy deposited in due to a scatter */
	double						attenuation;			/* Attenuation of crystal at current energy */
	double						fpToGo;					/* Free paths to travel in crystal */
	double						comptonToScatterProbability;		/* Ratio of prob of compton to prob of scatter */
	double						scatterProbability;		/* Probability of a scatter */
	double						interactionProbability;	/* Interaction type probability */
	double						depositedEnergy = 0.0;	/* Energy deposited at current location */
	double						depositedActiveEnergy = 0.0;	/* Energy deposited in active layers at current location */
	double						backOfLayer = 0.0;		/* "Back" of the block in local coordinates */
	double						frontOfLayer = 0.0;		/* "Front" of the block in local coordinates */
	double						freePathsToExit;		/* Free Paths to exit the block (for FD) */
	double						randFromExp;			/* Random number from exponential districution. */
	double						newWeight;				/* Temp for forced interaction adjustment */
	double						weight;					/* Incoming weight */
	
	do { /* Process Loop */
	
		/* Because we use the weight so much we use a local variable */
		weight = photonPtr->photon_current_weight;


		/* Get the current position/direction to have a working copy */
		pos = photonPtr->location;
		dir = photonPtr->angle;
		
		if ( (PHG_IsCollimateOnTheFly() == false) ||
			 (ColRunTimeParams[ColCurParams].ColType == ColEn_unc_spect) ) {
			/* Rotate the position/direction into detector coordinates */
			{
				pos.x_position = (photonPtr->location.x_position *
					PHGMATH_Cosine(photonPtr->detectorAngle)) +
					(photonPtr->location.y_position * PHGMATH_Sine(photonPtr->detectorAngle)) -
					DetRunTimeParams[DetCurParams].PlanarDetector.InnerRadius;

				pos.y_position = (-photonPtr->location.x_position *
					PHGMATH_Sine(photonPtr->detectorAngle)) +
					(photonPtr->location.y_position * PHGMATH_Cosine(photonPtr->detectorAngle));
					
				dir.cosine_x = (photonPtr->angle.cosine_x *
					PHGMATH_Cosine(photonPtr->detectorAngle)) +
					(photonPtr->angle.cosine_y * PHGMATH_Sine(photonPtr->detectorAngle));
		
				dir.cosine_y = (-photonPtr->angle.cosine_x *
					PHGMATH_Sine(photonPtr->detectorAngle)) +
					(photonPtr->angle.cosine_y * PHGMATH_Cosine(photonPtr->detectorAngle));
			}
		} else {
			/*	The slat collimator leaves the photon in the collimator coordinate system.
			 *	The only difference between this and the detector coordinate system is
			 *	in the zero point for the x-coordinate. */
			/* Change translated coordinates from x = 0.0 at collimator face
				to x = 0.0 at detector face
			*/
			pos.x_position = photonPtr->location.x_position -
			(DetRunTimeParams[DetCurParams].PlanarDetector.InnerRadius -
			ColGtInsideRadius());
		}
		
		/* If photon does not intersect the detector, finished */
		if (dir.cosine_x < 0.0)
			break;
			
		/*	Now compute how far it is from the back of the collimator to the face
			of the detector
		*/
		/* t = pos.x_position/dir.cosine_x; */
		
		t = (-pos.x_position/dir.cosine_x);
		
		/* Project to the face of the detector */
		pos.x_position = pos.x_position + (t * dir.cosine_x);
		pos.y_position = pos.y_position + (t * dir.cosine_y);
		pos.z_position = pos.z_position + (t * dir.cosine_z);

		/* The x position should be zero, due to small round off issues it is often close
			to zero but not exactly zero. The big problem is that it is often negative and
			we can't have that. So, a debug check will verify it is close to zero, while 
			always setting it to exactly zero afterwards.
		*/
		#ifdef PHG_DEBUG
		if (PhgMathRealNumAreEqual(pos.x_position, 0.0, -3, 0, 0, 0) == false) {
			PhgAbort("Invalid positioning of photon onto detector face (detTrackPlanar)", false);
		}
		#endif
	
		
		/* Set x to positive zero */
		pos.x_position = 0.0;
		
		/* Verify we haven't gone out of bounds in axial direction */
		if ((pos.z_position > detData[DetCurParams].detInBoundCyl.zMax) ||
				(pos.z_position < detData[DetCurParams].detInBoundCyl.zMin)) {
	
			/*	Finished */
			break;
		}
		
		/* Verify we haven't gone out of bounds in the transaxial direction */
		if ((pos.y_position > detData[DetCurParams].detPlnTransLimit) ||
				(pos.y_position < -detData[DetCurParams].detPlnTransLimit)) {
			
			/*	Finished */
			break;
		}
		
		/* Increment count of photons that reach the crystal */
		detData[DetCurParams].detTotReachingCrystal++;

		/* If we are doing FI then we compute the free paths to exit the block,
			pick an interaction point, and force the interaction to occur there. Otherwise
			we just sample for a free path value.
		*/
		if (DetRunTimeParams[DetCurParams].DoForcedInteraction) {
		
			/* Compute the free paths to exit */
			detPlnrComputeFreePathsToExit(pos, dir, photonPtr->energy, &freePathsToExit);
			
			/* Adjust the photon's weight by the probability that an interaction would occur */
			newWeight = weight * (1 - exp(-freePathsToExit));
				
			detData[DetCurParams].detWeightAdjusted += (weight - newWeight);
			
			weight = newWeight;
		
			/* Update weight */
			photonPtr->photon_current_weight = weight;
			
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
							
						#ifdef PHG_DEBUG
							if (fpToGo > freePathsToExit)
								PhgAbort("Invalid calculation of fpToGo for forced interaction (detTrackPlanar)", false);
						#endif
					}

				}	
			}

		}
		else {
			/*	Sample free paths to travel */
			PhgMathGetTotalFreePaths(&fpToGo);
		}
		
	
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
		backOfLayer = DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerDepth;
		frontOfLayer = 0.0;
		
		/*	Assign local variable to active status of layer. */
		layerIsActive = DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].IsActive;

		
		/* Begin tracking loop */
		do {
			/* Get the attenuation of the crystal at the current energy */
			SubObjGetAttenuationInTomo(
				DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerMaterial,
				photonPtr->energy, &attenuation);
			
			/* Compute the distance to travel, considering
			   no boundary intersections
			*/
			 distance = fpToGo/attenuation;

			 /* Project the photon the selected distance */
			detPlnrProjectToNewPos(pos, dir, distance, frontOfLayer, backOfLayer,
				&newPos, &distTraveled, &exitXFront, &exitXBack, &exitY, &exitZ);
			
			/* Test to see if we are out of the block */
			{
									
				/* See if the photon is going into the next layer of block */
				if (exitXBack == true) {
					
					/* Increment to next layer */
					curLayer++;
				
					/* If there are no more layers we are done */
					if (curLayer == (LbFourByte)(DetRunTimeParams[DetCurParams].PlanarDetector.NumLayers)) {
							
						/* Increment histogram of escaped weight */
						if (curInteraction < MAX_DET_INTERACTIONS)
							detData[DetCurParams].detWeightEscapedBins[curInteraction+1] += (weight*photonPtr->decay_weight);


						#ifdef PHG_DEBUG
							if ((curInteraction == -1) && (DetRunTimeParams[DetCurParams].DoForcedInteraction)) {
								sprintf(detPnrErrStr, "Escaping photon in x direction with FI on (detTrackPlanar)"
									"photon number is %lld", photonPtr->number);
									
								PhgAbort(detPnrErrStr, false);
							}
						#endif
						
						/* Break out */
						break;
					}
					
					/* We made it here so we are going into a new layer.
						First we adjust the free paths to account for those
						already traveled and then just continue to top of loop
					*/
										
					/* Update starting position to this location */
					pos = newPos;
					
					/* Adjust free paths to go */
					fpToGo -= (distTraveled * attenuation);
					
					#ifdef PHG_DEBUG
						if (fpToGo <= 0)
							PhgAbort("Invalid adjustment to fpToGo (detTrackPlanar)", false);
					#endif
					
					/* Compute new front and back of layer */
					frontOfLayer = backOfLayer;
					backOfLayer += DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerDepth;
					layerIsActive = DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].IsActive;
					
					/* Now "continue" to compute next position */
					continue;
					
				}

				/* See if the photon has "turned around" and headed out the face of the layer */
				if (exitXFront == true) {
					
					/* Decrement to previous layer */
					curLayer--;
				
					/* If there are no more layers we are done */
					if (curLayer == -1) {
							
						/* Increment histogram of escaped weight */
						if (curInteraction < MAX_DET_INTERACTIONS)
							detData[DetCurParams].detWeightEscapedBins[curInteraction+1] += (weight*photonPtr->decay_weight);


						#ifdef PHG_DEBUG
							if ((curInteraction == -1) && (DetRunTimeParams[DetCurParams].DoForcedInteraction)) {
								sprintf(detPnrErrStr, "Escaping photon in x direction with FI on (detTrackPlanar)"
									"photon number is %lld", photonPtr->number);
									
								PhgAbort(detPnrErrStr, false);
							}
						#endif
						
						/* Break out */
						break;
					}
					
					/* We made it here so we are going into a new layer.
						First we adjust the free paths to account for those
						already traveled and then start the loop over
					*/
					
					/* Update starting position to this location */
					pos = newPos;
					
					/* Adjust free paths to go */
					fpToGo -= (distTraveled * attenuation);
					
					#ifdef PHG_DEBUG
						if (fpToGo <= 0)
							PhgAbort("Invalid adjustment to fpToGo (detTrackPlanar)", false);
					#endif
					
					/* Compute new back of layer */
					backOfLayer = frontOfLayer;
					frontOfLayer -= DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerDepth;
					layerIsActive = DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].IsActive;
					
					/* Now "continue" to compute next position */
					continue;
				}
				
				
				/* Check y dimension */
				if (exitY == true) {
					
					/* Increment histogram of escaped weight */
					if (curInteraction < MAX_DET_INTERACTIONS)
						detData[DetCurParams].detWeightEscapedBins[curInteraction+1] += (weight*photonPtr->decay_weight);

					#ifdef PHG_DEBUG_YE_IMAGES
						if ((curInteraction == -1) && (DetRunTimeParams[DetCurParams].DoForcedInteraction)) {
							sprintf(detPnrErrStr, "Escaping photon in y direction with FI on (detTrackPlanar)"
								"photon number is %lld", photonPtr->number);
								
							PhgAbort(detPnrErrStr, false);
						}
						
						if (curInteraction == -1) {
							detPlnrExitYWeight += (weight*photonPtr->decay_weight);
							detPlnrExitYWeightSq += PHGMATH_Square(weight*photonPtr->decay_weight);							
							detPlnrExitYCount++;
						}
						
					#endif

					/* Break out */
					break;
				}
			
				/* Verify photon hasn't gone out of bounds in axial direction */
				if (exitZ == true) {

					/* Increment histogram of escaped weight */
					if (curInteraction < MAX_DET_INTERACTIONS)
						detData[DetCurParams].detWeightEscapedBins[curInteraction+1] += (weight*photonPtr->decay_weight);
										

					#ifdef PHG_DEBUG_YE_IMAGES
						if ((curInteraction == -1) && (DetRunTimeParams[DetCurParams].DoForcedInteraction)) {
							sprintf(detPnrErrStr, "Escaping photon in z direction with FI on (detTrackPlanar)"
								"photon number is %lld", photonPtr->number);
								
							PhgAbort(detPnrErrStr, false);
						}
						
						if (curInteraction == -1) {
							detPlnrExitZWeight += (weight*photonPtr->decay_weight);
							detPlnrExitZWeightSq += PHGMATH_Square(weight*photonPtr->decay_weight);							
							detPlnrExitZCount++;
						}
					#endif
	
					/* Break out */
					break;
				}
			}
			
			#ifdef PHG_DEBUG
				/* Verify positive x value */
				if (newPos.x_position < 0.0) {
					sprintf(detPnrErrStr, "Invalid x pos calculated, x = %3.2f, dir.cosine_x = %3.2f\n",
						newPos.x_position, dir.cosine_x);
					ErAlert(detPnrErrStr, false);
				}
			#endif
			
			/* We are here so do an interaction */
			{
				/* Increment count of active interactions */
				curInteraction++;
				
				#ifdef PHG_DEBUG
					if (curInteraction == MAX_DET_INTERACTIONS) {
						sprintf(detPnrErrStr, "curInteraction incremented to MAX_DET_INTERACTIONS at line %d in (detTrackPlanar)\n", __LINE__);
						PhgAbort(detPnrErrStr, false);
					}
				#endif
							
				
				/* Save the current position */
				photonPtr->det_interactions[curInteraction].pos = newPos;
				
				/* Save the active status */
				photonPtr->det_interactions[curInteraction].isActive = layerIsActive;
				
				/* Get probability of  scatter */
				scatterProbability = SubObjGetProbScatterInTomo2(
					DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerMaterial,
					photonPtr->energy);
				
				/* Get probability of compton scatter */
				comptonToScatterProbability =
					SubObjGetProbComptToScatInTomo2(
					DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerMaterial,
					photonPtr->energy);
				
				/* Sample for probability of photo-electric absorbtion */
				interactionProbability = PhgMathGetRandomNumber();
				
				#ifdef PHG_DEBUG
					if (PHGDEBUG_DetAbsorbOnly())
						interactionProbability = 1.0;
				#endif
				
				/* See if photon is absorbed */
				if (interactionProbability > scatterProbability) {
					
					/* Save the current energy */
					photonPtr->det_interactions[curInteraction].energy_deposited =
						photonPtr->energy;
					
					photonPtr->energy = 0.0;
					
					/* Sum the deposited energy */
					depositedEnergy += photonPtr->det_interactions[curInteraction].energy_deposited;
					
					/* Sum active layer energy */
					if (layerIsActive)
						depositedActiveEnergy += photonPtr->det_interactions[curInteraction].energy_deposited;
						
					/* Increment counter */
					detData[DetCurParams].detTotPhotonsAbsorbed++;
					
					/* Increment weight counter */
					detData[DetCurParams].detTotWtAbsorbed += (weight*photonPtr->decay_weight);
					if (curInteraction < MAX_DET_INTERACTIONS)
						detData[DetCurParams].detWeightAbsorbedBins[curInteraction+1] += (weight*photonPtr->decay_weight);
					
					/* If it is the first interaction, update statistics for 
						absorbtion on first interaction
					*/
					if (curInteraction == 0) {
						detData[DetCurParams].detTotFirstTimeAbsorptions++;
						detData[DetCurParams].detTotWtFirstTimeAbsorbed += (weight*photonPtr->decay_weight);
					}
					
					/* Break out of tracking loop */
					break;
				}
				
				/* We made it to here, so compute scatter angle and new energy */
				{
					/* First store the current energy (it gets changed in EmisListDoComptonInteraction) */
					energyDepositedFromInteraction = photonPtr->energy;

					/* Compute the interaction */
					{
						/* First update the photon's direction to the rotated coordinates */
						photonPtr->angle = dir;
						
						if (interactionProbability > (scatterProbability * comptonToScatterProbability)) {
						
							if (PHG_IsModelCoherentInTomo()) {
								EmisListDoCoherent(photonPtr,
									DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerMaterial);
							}
							else {
								/* Get the next free paths to travel */
								PhgMathGetTotalFreePaths(&fpToGo);
							
								/* Decrement interactions because in this case we decide not to do one */
								curInteraction--;

								/* Update the position */
								pos = newPos;

								continue; 
							}
							
						}
						else  {
							/* Compute scatter angle and new energy */
							EmisListDoComptonInteraction(photonPtr);
						}
						
						/* Now update the local, working copy of the direction */
						dir = photonPtr->angle;
					}
					
					
					/* If the photon energy is below supported minimum energy, force an absorption */
					if (photonPtr->energy < PHG_MIN_PHOTON_ENERGY) {
					
						/* Save the energy deposited */
						depositedEnergy += energyDepositedFromInteraction;

						/* Sum active layer energy */
						if (layerIsActive)
							depositedActiveEnergy += energyDepositedFromInteraction;
							
						if (curInteraction < MAX_DET_INTERACTIONS)
							detData[DetCurParams].detWeightAbsorbedBins[curInteraction+1] += (weight*photonPtr->decay_weight);
						
						/* Save the energy  deposited */
						photonPtr->det_interactions[curInteraction].energy_deposited =
							energyDepositedFromInteraction;
							
						/* Increment counters */
						detData[DetCurParams].detTotPhotonsAbsorbed++;
						detData[DetCurParams].detTotForcedAbsorptions++;
						detData[DetCurParams].detTotWtForcedAbsorbed += (weight*photonPtr->decay_weight);
						
						/* Break out */
						break;
							
					}
						
					/* 	If we made it here it was a normal scatter.
						Compute the deposited energy (irreleavant if inactive layer)
					*/
					energyDepositedFromInteraction -= photonPtr->energy;
				}
							
				/* Save the current energy */
				photonPtr->det_interactions[curInteraction].energy_deposited =
					energyDepositedFromInteraction;
				
				/* Sum the deposited energy */
				depositedEnergy += energyDepositedFromInteraction;
				
				/* Sum active layer energy */
				if (layerIsActive)
					depositedActiveEnergy += energyDepositedFromInteraction;

				/* Update the position */
				pos = newPos;
				
			} /* End of "We are here so do an interaction" */
			
			/* Get the next free paths to travel */
			PhgMathGetTotalFreePaths(&fpToGo);

		} while (true);
		
		/* Store number of interactions */
		photonPtr->num_det_interactions = curInteraction+1;
		
		/* Compute centroid if energy deposited */
		if (depositedActiveEnergy > 0.0) {	
			
			/* Increment counter */
			detData[DetCurParams].detTotPhotonsDepositingEnergy++;
			
			photonPtr->photon_current_weight = weight;

			detPlnrCompCentroid(photonPtr, curInteraction, depositedActiveEnergy);
			
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
			
			
			/* Clear out the direction vector, it's no longer meaningful */
			photonPtr->angle.cosine_x = 0.0;
			photonPtr->angle.cosine_y = 0.0;
			photonPtr->angle.cosine_z = 0.0;
			
			/* Mark this as a successful detection */
			photonDetected = true;
		}
		
		if (photonPtr->num_det_interactions == 0) {
			
			/* Increment counter for number of photons that don't interact */
			detData[DetCurParams].detTotPhotonsPassingThrough++;
		}

	} while (false);
	
	return (photonDetected);
}

/*********************************************************************************
*
*			Name:			detGtDetectableAngle
*
*			Summary:		Selects a detector angle for the given photon direction.
*
*			Arguments:
*					PHG_Direction	dir				- Direction of travel.
*					PHG_Position	pos				- Initial position.
*					float			*detectorPos	- Selected detector position.
*			Function return: True if photon hits the detector.
*
*********************************************************************************/
LbUsFourByte detGtDetectableAngle(PHG_Direction bDir, PHG_Position bPos, 
				PHG_Direction pDir, PHG_Position pPos, double *detectorPos)
{
	double				detAngle;
	double				bBeta, pBeta;
	double				distance;
	double				detPositions[1440]; /* 360*4 = 1/4 degree sampling this should be dynamically allocated via numViews but I'm lazy */			
	double				r1Min, r1Max, r2Min, r2Max;
	double				posSize;
	double				dPos;
	LbFourByte			numViews = DetGtNumViews();
	LbFourByte			i;
	LbUsFourByte		j, k;
	PHG_Position		bigPos;

	do { /* Process Loop */

/** NOTE: This routine has been modified until the complete algorithm can be implemented.
	There is a RETURN statement at the end of this block which is always executed.
***/
{
	
	/* Project blue photon to "big" cylinder */
	(void) CylPosProjectToCylinder(&bPos, &bDir, &detData[DetCurParams].detPlnrBigCylinder, &bigPos, &distance);
		
	/* Compute polar coordinate angle (shift into range [0, 2pi]) */
	bBeta = atan2(bigPos.y_position,bigPos.x_position);
	if (bBeta < 0)
		bBeta += (PHGMATH_2PI);
	
	/* Pick a random detector position */
	dPos = PhgMathGetRandomNumber() * PHGMATH_PI;
	
	/* Determine which head the photon will hit */
	if (fabs(dPos-bBeta) > (PHGMATH_PI_DIV2)) {
		detectorPos[0] = dPos + PHGMATH_PI;
	}
	else {
		detectorPos[0] = dPos;
	}
	
	/* Project pink photon to big cylinder */
	(void) CylPosProjectToCylinder(&pPos, &pDir, &detData[DetCurParams].detPlnrBigCylinder, &bigPos, &distance);
		
	/* Compute polar coordinate angle (shift into range [0, 2pi]) */
	pBeta = atan2(bigPos.y_position,bigPos.x_position);
	if (pBeta < 0)
		pBeta += (PHGMATH_2PI);
	
	/* Determine which head the photon will hit */
	if (fabs(dPos-pBeta) > (PHGMATH_PI_DIV2)) {
		detectorPos[1] = dPos + PHGMATH_PI;
	}
	else {
		detectorPos[1] = dPos;
	}

	return(2);
}	
	
	/* Project pink photon to big cylinder */
	(void) CylPosProjectToCylinder(&pPos, &pDir, &detData[DetCurParams].detPlnrBigCylinder, &bigPos, &distance);
		
	/* Compute polar coordinate angle (shift into range [0, 2pi]) */
	pBeta = atan2(bigPos.y_position,bigPos.x_position);
	if (pBeta < 0)
		pBeta += (PHGMATH_2PI);

	
	/* Initialize */
	i = j = k = 0;
	
	
	
	/* Compute ranges */
	{
		/* Initialize */
		r1Min = r1Max = r2Min = r2Max = 0;

		/* Set range */
		if (bBeta >= pBeta) {
			if ((pBeta + (2 * detData[DetCurParams].detPlnrAngularCoverage)) > bBeta) {
				r1Max = pBeta + detData[DetCurParams].detPlnrAngularCoverage;
				r1Min = bBeta - detData[DetCurParams].detPlnrAngularCoverage;
			}
		}
		else if ((bBeta + (2 * detData[DetCurParams].detPlnrAngularCoverage)) > pBeta) {
			r1Max = bBeta + detData[DetCurParams].detPlnrAngularCoverage;
			r1Min = pBeta - detData[DetCurParams].detPlnrAngularCoverage;
		}
	
		/* Shift into 0->180 */
		if (r1Max >= 360) {
			r2Min = r1Min - 180;
			r2Max = 180;
			r1Min = 0;
			r1Max = r1Max - 360;
		}
		else if (r1Min >= 180) {
			r1Min = r1Min - 180;
			r1Max = r1Max - 180;
		}
		else if (r1Max >= 180) {
			r2Min = r1Min;
			r2Max = 180;
			r1Min = 0;
			r1Max = r1Max - 180;
		}
		else if (r1Min < 0) {
			r2Min = r1Min + 180;
			r2Max = 180;
			r1Min = 0;
		}
	}
	
	/* Setup loop variables */
	posSize = detData[DetCurParams].detViewSize/numViews;
	detAngle = posSize/2;
	
	/* Loop over viable detector angles */
	for (i = 0; i < numViews; i++) {
		
		/* If detector position is in range store it */
		if ((detAngle > r1Min) && (detAngle < r1Max)) {
			detPositions[j] = detAngle;
			j++;
		}
		else if ((detAngle > r2Min) && (detAngle < r2Max)) {
			detPositions[j] = detAngle;
			j++;
		}
		
		detAngle += posSize;
	}

	/* Pick a detector position position */
	if (j > 0) {
	
		k = j * PhgMathGetRandomNumber();
		
		*detectorPos = detPositions[k];
		
		/* See if we are looking at a blue or pink position */
		if ((pBeta >= (*detectorPos-detData[DetCurParams].detPlnrAngularCoverage)) && (pBeta <= (*detectorPos+detData[DetCurParams].detPlnrAngularCoverage)))
			*detectorPos += 180;
			
		/* Convert to radians */
		*detectorPos = PHGMATH_RadiansFromDegrees(*detectorPos);
	}
	
	} while (false);
	return(k);	
}

/*********************************************************************************
*
*			Name:			detDHDoOptimalPos
*
*			Summary:		Perform detection for dual headed detectors using
*							optimal detector head position
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
void detDHDoOptimalPos(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		DetectedPhotonsTy *detPhotonsPtr)
{
	double sumWeight;
	double	pairProbs[PHG_MAX_DETECTED_PHOTONS*PHG_MAX_DETECTED_PHOTONS];
	double	luckyProb;
	double	detectorPos[2];
	LbUsFourByte	b, p;
	LbUsFourByte	numDetPositions;
	
	do { /* Process Loop */
	/* First, sum the weight products */
	sumWeight = 0;
	for (b = 0; b < numBluePhotons; b++) {
		for (p = 0; p < numPinkPhotons; p++) {
			sumWeight += bluePhotons[b].photon_current_weight *
				pinkPhotons[p].photon_current_weight;
			pairProbs[(b*numPinkPhotons)+p] = sumWeight;
		}
	}
	
	/* Build the cumulative probability list */
	for (b = 0; b < numBluePhotons; b++) {
		for (p = 0; p < numPinkPhotons; p++) {
			pairProbs[(b*numPinkPhotons)+p] /= sumWeight;
		}
	}

	/* Pick the lucky probability */
	luckyProb = PhgMathGetRandomNumber();
	
	/* Pick the lucky photon pair */
	for (b = 0; b < numBluePhotons; b++) {
		for (p = 0; p < numPinkPhotons; p++) {
		
			/* If our lucky probability is less than this bin's probability, we are done */
			if (luckyProb <= pairProbs[(b*numPinkPhotons)+p])
				goto LUCKY_FOUND;
		}
	}
	LUCKY_FOUND:;
	#ifdef PHG_DEBUG
		if ((b == numBluePhotons) || (p == numPinkPhotons)) {
			PhgAbort("Error in finding lucky probability (detDualHeaded)", false);
		}
	#endif

	/*	Select position for the detector */
	numDetPositions = detGtDetectableAngle(bluePhotons[b].angle, bluePhotons[b].location,
		pinkPhotons[p].angle, pinkPhotons[p].location, detectorPos);
	
	/* if detection not possible, we're done */
	if (numDetPositions == 0)
		break;
	
#ifdef FIGURED_IT_OUT		
	/* Adjust the weight 
	decayPtr->eventWeight *= numDetPositions/DetGtNumViews(); */
#endif

	/* Loop through all blue photons */
	for (b = 0; b < numBluePhotons; b++) {
		
		/* Set the decay position */			
		bluePhotons[b].detectorAngle = detectorPos[0];
		
		/* Let user modify and/or reject photons */
		if (DetUsrModPETPhotonsFPtr && 
				(*DetUsrModPETPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
					&(bluePhotons[b])) == false) {
			
			/* They rejected it so go to next photon */
			continue;
		}
		
		/* Track the photon through the planar detector */
		if (detTrackPlanar(&bluePhotons[b]) == false)
			continue;	
					
		/* We made it to here so save the detected photons */
		{
			detPhotonsPtr->DetectedTrkngBluePhotons[detPhotonsPtr->NumDetectedBluePhotons]
				= bluePhotons[b];

			/* Increment the counters */
			detPhotonsPtr->NumDetectedBluePhotons++;

			/* Update our detected photon block if doing history file */
			if  (DET_IsDoHistory()) {
				DetUpdateDetectedPhotonBlock(&bluePhotons[b]);
			}
		}
	}
	
	/* Loop through all pink photons */
	for (p = 0; p < numPinkPhotons; p++) {

		/* Set the decay position */			
		pinkPhotons[p].detectorAngle = detectorPos[1];

		/* Let user modify and/or reject photons */
		if (DetUsrModPETPhotonsFPtr && 
				(*DetUsrModPETPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
					&(pinkPhotons[p])) == false) {
			
			/* They rejected it so go to next pink */
			continue;
		}

		/* Track the photon through the planar detector */
		if (detTrackPlanar(&pinkPhotons[p]) == false)
			continue;	
			
		/* We made it to here so save the detected photons */
		{
			detPhotonsPtr->DetectedTrkngPinkPhotons[detPhotonsPtr->NumDetectedPinkPhotons]
				= pinkPhotons[p];

			/* Increment the counters */
			detPhotonsPtr->NumDetectedPinkPhotons++;

			/* Update our detected photon block if doing history file */
			if  (DET_IsDoHistory()) {
				DetUpdateDetectedPhotonBlock(&pinkPhotons[p]);
			}
		}
	}
} while (false);
}

/*********************************************************************************
*
*			Name:			detDumpPlnrInteractions
*
*			Summary:		Prints out interaction list for debugging purposes.
*
*			Arguments:
*				LbFourByte					numInteractions	- The number of interactions.
*				detInteractionInfoTy	interactions	- The interactions.
*
*			Function return: None.
*
*********************************************************************************/
void detDumpPlnrInteractions(LbFourByte numInteractions,
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


/*********************************************************************************
*
*		Name:			DetPlnrInitPhotons
*
*		Summary:		Initialize the photons.
*
*		Arguments:
*			PHG_Decay				*decayPtr		- The decay that started the process.
*			PHG_TrackingPhoton 		*photonPtr		- The detected photon.
*
*		Function return:	None.
*
*********************************************************************************/

void DetPlnrInitPhotons(PHG_Decay *decayPtr, PHG_TrackingPhoton *photonPtr)

{
	static double		detPosition;	/* Randomly calculated detector position */
	
	
	/* If this is a new decay, then...  (do this only once per decay) */
	if (decayPtr != detPlnrLastDecayPtr) {
		detPlnrLastDecayPtr = decayPtr;
		
		/* Compute the detector position, either incremental or continuous.
			if collimation is done on the fly, then the detector position
			is already established.
		 */
		if (PHG_IsCollimateOnTheFly() == false) {
			if (DetGtNumViews() > 0) {
				detPosition = DetRunTimeParams[DetCurParams].PlanarDetector.MinAngle +
					((detData[DetCurParams].detPlnrDelta/2) + 
					(floor(PhgMathGetRandomNumber() * DetGtNumViews()) * 
							detData[DetCurParams].detPlnrDelta));
			}
			else {
				detPosition = 	PhgMathGetRandomNumber() * PHGMATH_2PI;
			}
		}
	}
	
	/* If collimation has not been done then we need the detector position */
	if (PHG_IsCollimateOnTheFly() == false) {
		if (DetRunTimeParams[DetCurParams].DetectorType == DetEn_DualHeaded) {
			/* Dual-headed coincidence detectors */
			DetGetDetAngle(&(photonPtr->location), &(photonPtr->angle),
				detPosition, &(photonPtr->detectorAngle));
		}
		else {
			/* SPECT */
			photonPtr->detectorAngle = detPosition;
		}
	}
}


/*********************************************************************************
*
*		Name:			DetPlnrGtTruncatedFreePaths
*
*		Summary:		Returns the free paths to travel for a forced first interaction.
*
*		Arguments:
*			PHG_TrackingPhoton		*photonPtr		- The photon.
*			double					*weight			- The photon weight.
*
*		Function return: Free-paths to travel.
*
*********************************************************************************/

double DetPlnrGtTruncatedFreePaths(PHG_TrackingPhoton *photonPtr, double *weight)

{
	double				freePathsToExit;	/* Free paths to exit the block */
	double				localWeight;		/* Incoming weight */
	double				newWeight;			/* Temp for forced interaction adjustment */
	double				randFromExp;		/* Random number from exponential distribution */
	LbUsFourByte		wholeNum;			/* Whole number temp variable */
	double				fpToGo;				/* Free paths to travel in crystal */
	
	
	/* Compute the free paths to exit */
	detPlnrComputeFreePathsToExit(photonPtr->location, photonPtr->angle, photonPtr->energy, 
									&freePathsToExit);
	
	/* Adjust the photon's weight by the probability that an interaction would occur */
	localWeight = *weight;
	newWeight = localWeight * (1 - exp(-freePathsToExit));
	detData[DetCurParams].detWeightAdjusted += (localWeight - newWeight);
	*weight = newWeight;
	
	/* Pick a free path to go based on a truncated exponential distribution */
	{
		/* Start with a random from the exponential distribution */
		PhgMathGetTotalFreePaths(&randFromExp);
		
		/* Truncate to desired range */
		{
			/* Remove whole number part NOTE I AM IGNORING THE
			   POSSIBILITY OF INTEGER OVERFLOW
			*/
			wholeNum = (LbUsFourByte) (randFromExp/freePathsToExit);
			
			fpToGo = ((randFromExp/freePathsToExit) - wholeNum) * freePathsToExit;
			
			#ifdef PHG_DEBUG
				if (fpToGo > freePathsToExit)
					PhgAbort("Invalid calculation of fpToGo for forced interaction (DetPlnrGtTruncatedFreePaths)", false);
			#endif
		}	
	}
	
	return(fpToGo);
}


/*********************************************************************************
*
*		Name:			DetPlnrProjectToDetector
*
*		Summary:		Project the photon to the detector edge.
*						Also, convert its coordinates to detector ones.
*
*		Arguments:
*			PHG_TrackingPhoton 	*photonPtr	- The photon to track and detect.
*			LbFourByte			*ringNum	- The ring projected to.
*			LbFourByte			*detNum		- The detector projected to.
*
*		Function return: 	True if photon not rejected (still valid).
*
*********************************************************************************/

Boolean DetPlnrProjectToDetector(PHG_TrackingPhoton *photonPtr, 
									LbFourByte *ringNum, LbFourByte *detNum)

{
	Boolean				valid;				/* True if photon still valid */
	PHG_Position		pos;				/* The photon position */
	PHG_Direction		dir;				/* The photon direction of travel */
	double				t;					/* Projection distance */
	
	
	valid = false;
	
	do {
		/* Get the current position/direction to have a working copy */
		pos = photonPtr->location;
		dir = photonPtr->angle;
		
		/* Transform the photon coordinates to detector coordinates */
		if ( (PHG_IsCollimateOnTheFly() == false) || 
				(ColRunTimeParams[ColCurParams].ColType == ColEn_unc_spect) ) {
			/* Rotate the position/direction into detector coordinates */
			{
				pos.x_position = (photonPtr->location.x_position *
					PHGMATH_Cosine(photonPtr->detectorAngle)) +
					(photonPtr->location.y_position * PHGMATH_Sine(photonPtr->detectorAngle)) -
					DetRunTimeParams[DetCurParams].PlanarDetector.InnerRadius;
				
				pos.y_position = (-photonPtr->location.x_position *
					PHGMATH_Sine(photonPtr->detectorAngle)) +
					(photonPtr->location.y_position * PHGMATH_Cosine(photonPtr->detectorAngle));
				
				dir.cosine_x = (photonPtr->angle.cosine_x *
					PHGMATH_Cosine(photonPtr->detectorAngle)) +
					(photonPtr->angle.cosine_y * PHGMATH_Sine(photonPtr->detectorAngle));
				
				dir.cosine_y = (-photonPtr->angle.cosine_x *
					PHGMATH_Sine(photonPtr->detectorAngle)) +
					(photonPtr->angle.cosine_y * PHGMATH_Cosine(photonPtr->detectorAngle));
			}
		}
		else {
			/*	The slat collimator leaves the photon in the collimator coordinate system.
			 *	The only difference between this and the detector coordinate system is
			 *	in the zero point for the x-coordinate.  */
			/* Change translated coordinates from x = 0.0 at collimator face
				to x = 0.0 at detector face.
			*/
			pos.x_position = photonPtr->location.x_position -
								(DetRunTimeParams[DetCurParams].PlanarDetector.InnerRadius -
									ColGtInsideRadius());
		}
		
		/* If photon does not intersect the detector, finished */
		if (dir.cosine_x < 0.0) {
			break;
		}
		
		
		/*	Now compute how far it is from the back of the collimator to the face
			of the detector */
		t = (-pos.x_position/dir.cosine_x);
		
		/* Project to the face of the detector */
		pos.x_position = pos.x_position + (t * dir.cosine_x);
		pos.y_position = pos.y_position + (t * dir.cosine_y);
		pos.z_position = pos.z_position + (t * dir.cosine_z);
		
		/* The x position should be zero, but due to small round off issues it is often close
			to zero but not exactly zero.  The big problem is that it is often negative and
			we can't have that.  So, a debug check will verify it is close to zero, while 
			always setting it to exactly zero afterwards.
		*/
		#ifdef PHG_DEBUG
			if (PhgMathRealNumAreEqual(pos.x_position, 0.0, -3, 0, 0, 0) == false) {
				PhgAbort("Invalid positioning of photon onto detector face (DetPlnrProjectToDetector)", false);
			}
		#endif
		
		
		/* Set x to positive zero */
		pos.x_position = 0.0;
		
		/* Verify we haven't gone out of bounds in the axial direction */
		if ((pos.z_position > detData[DetCurParams].detInBoundCyl.zMax) ||
				(pos.z_position < detData[DetCurParams].detInBoundCyl.zMin)) {
			
			/*	Finished */
			break;
		}
		
		/* Verify we haven't gone out of bounds in the transaxial direction */
		if ((pos.y_position > detData[DetCurParams].detPlnTransLimit) ||
				(pos.y_position < -detData[DetCurParams].detPlnTransLimit)) {
			
			/*	Finished */
			break;
		}
		
		/* Photon will reach detector */
		valid = true;
		
		/* Change photon to detector coordinates */
		photonPtr->location = pos;
		photonPtr->angle = dir;
		
	} while (false);
	
	
	/* Always both 0 for planar detectors */
	*ringNum = 0;
	*detNum = 0;
	
	return(valid);
}


/*********************************************************************************
*
*			Name:			DetPlnrFindNextInteraction
*
*			Summary:		Find and move the photon to the next interaction point.
*
*			Arguments:
*				PHG_TrackingPhoton 		*photonPtr		- The photon to track and detect.
*				LbFourByte				interactionNum	- Current interaction count.
*				LbFourByte				*curDetNum		- The current detector number.
*				detEn_ActionTy			*actionType		- The type of interaction that occurred.
*				double					*fpToGo			- Remaining free paths to travel.
*				LbUsFourByte			*detMaterial	- The interaction material of the detector.
*				Boolean					*isActive		- Active status of the interaction material.
*
*			Function return:	None.
*
*********************************************************************************/

void DetPlnrFindNextInteraction(PHG_TrackingPhoton *photonPtr, 
								LbFourByte interactionNum, LbFourByte *curDetNum, 
								detEn_ActionTy *actionType, double *fpToGo, 
								LbUsFourByte *detMaterial, Boolean *isActive)

{
	detEn_ActionTy				detAction;		/* Photon action as detector type */
	LbFourByte					curLayer;		/* Current layer */
	double						frontOfLayer;	/* "Front" of the block in local coordinates */
	double						backOfLayer;	/* "Back" of the block in local coordinates */
	LbFourByte					i;				/* Index through layers */
	Boolean						layerIsActive;	/* Flags current layer as active or not */
	double						weight;			/* Incoming weight */
	double						attenuation;	/* Attenuation of current material */
	double						distance;		/* Distance variable */
	PHG_Position				newPos;			/* Projected position */
	double						distTraveled;	/* Truncated distance for boundaries */
	Boolean						exitXFront;		/* Flags exit through front of layer */
	Boolean						exitXBack;		/* Flags exit through back of layer */
	Boolean						exitY;			/* Flags exit through y limit */
	Boolean						exitZ;			/* Flags exit through z limit */
	
	
	/* Initialize values */
	
	detAction = detEnAc_Null;
	
	curLayer = *curDetNum;
	frontOfLayer = 0.0;
	backOfLayer = DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[0].LayerDepth;
	for (i=1; i<=curLayer; i++) {
		frontOfLayer = backOfLayer;
		backOfLayer += DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[i].LayerDepth;
	}
	layerIsActive = DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].IsActive;
	
	weight = photonPtr->photon_current_weight;
	
	
	/* Get the attenuation of the crystal at the current energy */
	SubObjGetAttenuationInTomo(
		DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerMaterial,
		photonPtr->energy, &attenuation);
	
	/* Compute the distance to travel, considering
	   no boundary intersections
	*/
	distance = *fpToGo/attenuation;
	
	/* Project the photon the selected distance */
	detPlnrProjectToNewPos(photonPtr->location, photonPtr->angle, distance, frontOfLayer, backOfLayer,
		&newPos, &distTraveled, &exitXFront, &exitXBack, &exitY, &exitZ);
	
	/* Test to see if we are out of the block */
	do {
		
		/* See if the photon is going into the next layer of block */
		if (exitXBack == true) {
			
			/* Increment to next layer */
			curLayer++;
			
			/* If there are no more layers we are done */
			if (curLayer == (LbFourByte)(DetRunTimeParams[DetCurParams].PlanarDetector.NumLayers)) {
				
				detAction = detEnAc_Discard;
				
				/* Increment histogram of escaped weight */
				if (interactionNum < MAX_DET_INTERACTIONS)
					detData[DetCurParams].detWeightEscapedBins[interactionNum+1] += 
						(weight*photonPtr->decay_weight);
				
				#ifdef PHG_DEBUG
					if ((interactionNum == -1) && (DetRunTimeParams[DetCurParams].DoForcedInteraction)) {
						sprintf(detPnrErrStr, "Escaping photon in x direction with FI on (DetPlnrFindNextInteraction)"
							"photon number is %lld", photonPtr->number);
						
						PhgAbort(detPnrErrStr, false);
					}
				#endif
				
				/* Break out */
				break;
			}
			
			/* We made it here so we are going into a new layer.
				First we adjust the free paths to account for those
				already traveled and then just continue
			*/
			
			detAction = detEnAc_LayerCross;
			
			/* Adjust free paths to go */
			*fpToGo -= (distTraveled * attenuation);
			
			#ifdef PHG_DEBUG
				if (*fpToGo <= 0)
					PhgAbort("Invalid adjustment to fpToGo (DetPlnrFindNextInteraction)", false);
			#endif
			
			/* Compute new front and back of layer */
			frontOfLayer = backOfLayer;
			backOfLayer += DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerDepth;
			layerIsActive = DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].IsActive;
			
			/* Now "continue" to compute next position */
			continue;
		}
		
		
		/* See if the photon has "turned around" and headed out the face of the layer */
		if (exitXFront == true) {
			
			/* Decrement to previous layer */
			curLayer--;
			
			/* If there are no more layers we are done */
			if (curLayer == -1) {
				
				detAction = detEnAc_Discard;
				
				/* Increment histogram of escaped weight */
				if (interactionNum < MAX_DET_INTERACTIONS)
					detData[DetCurParams].detWeightEscapedBins[interactionNum+1] += (weight*photonPtr->decay_weight);
				
				#ifdef PHG_DEBUG
					if ((interactionNum == -1) && (DetRunTimeParams[DetCurParams].DoForcedInteraction)) {
						sprintf(detPnrErrStr, "Escaping photon in x direction with FI on (DetPlnrFindNextInteraction)"
							"photon number is %lld", photonPtr->number);
						
						PhgAbort(detPnrErrStr, false);
					}
				#endif
				
				/* Break out */
				break;
			}
			
			/* We made it here so we are going into a new layer.
				First we adjust the free paths to account for those
				already traveled and then continue
			*/
			
			detAction = detEnAc_LayerCross;
			
			/* Adjust free paths to go */
			*fpToGo -= (distTraveled * attenuation);
			
			#ifdef PHG_DEBUG
				if (*fpToGo <= 0)
					PhgAbort("Invalid adjustment to fpToGo (DetPlnrFindNextInteraction)", false);
			#endif
			
			/* Compute new back of layer */
			backOfLayer = frontOfLayer;
			frontOfLayer -= DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerDepth;
			layerIsActive = DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].IsActive;
			
			/* Now "continue" to compute next position */
			continue;
		}
		
		
		/* Check y dimension */
		if (exitY == true) {
			
			detAction = detEnAc_Discard;
			
			/* Increment histogram of escaped weight */
			if (interactionNum < MAX_DET_INTERACTIONS)
				detData[DetCurParams].detWeightEscapedBins[interactionNum+1] += (weight*photonPtr->decay_weight);
			
			#ifdef PHG_DEBUG_YE_IMAGES
				if ((interactionNum == -1) && (DetRunTimeParams[DetCurParams].DoForcedInteraction)) {
					sprintf(detPnrErrStr, "Escaping photon in y direction with FI on (DetPlnrFindNextInteraction)"
						"photon number is %lld", photonPtr->number);
					
					PhgAbort(detPnrErrStr, false);
				}
				
				if (interactionNum == -1) {
					detPlnrExitYWeight += (weight*photonPtr->decay_weight);
					detPlnrExitYWeightSq += PHGMATH_Square(weight*photonPtr->decay_weight);							
					detPlnrExitYCount++;
				}
			#endif
			
			/* Break out */
			break;
		}
		
		
		/* Verify photon hasn't gone out of bounds in axial direction */
		if (exitZ == true) {
			
			detAction = detEnAc_Discard;
			
			/* Increment histogram of escaped weight */
			if (interactionNum < MAX_DET_INTERACTIONS)
				detData[DetCurParams].detWeightEscapedBins[interactionNum+1] += (weight*photonPtr->decay_weight);
			
			#ifdef PHG_DEBUG_YE_IMAGES
				if ((interactionNum == -1) && (DetRunTimeParams[DetCurParams].DoForcedInteraction)) {
					sprintf(detPnrErrStr, "Escaping photon in z direction with FI on (DetPlnrFindNextInteraction)"
						"photon number is %lld", photonPtr->number);
					
					PhgAbort(detPnrErrStr, false);
				}
				
				if (interactionNum == -1) {
					detPlnrExitZWeight += (weight*photonPtr->decay_weight);
					detPlnrExitZWeightSq += PHGMATH_Square(weight*photonPtr->decay_weight);							
					detPlnrExitZCount++;
				}
			#endif
			
			/* Break out */
			break;
		}
		
		
		/* Otherwise, an interaction occurs here */
		detAction = detEnAc_Interact;
		
		#ifdef PHG_DEBUG
			/* Verify positive x value */
			if (newPos.x_position < 0.0) {
				sprintf(detPnrErrStr, "Invalid x pos calculated, x = %3.2f, dir.cosine_x = %3.2f\n",
					newPos.x_position, photonPtr->angle.cosine_x);
				ErAlert(detPnrErrStr, false);
			}
		#endif
		
	} while (false);
	
	
	/* Update the detector position */
	*curDetNum = curLayer;
	
	/* Update the position and distance traveled */
	photonPtr->location = newPos;
	photonPtr->travel_distance += distTraveled;
	
	/* Record the action, detector material, and layer active status */
	*actionType = detAction;
	if (detAction != detEnAc_Discard) {
		*detMaterial = DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerMaterial;
		*isActive = layerIsActive;
	}
}
