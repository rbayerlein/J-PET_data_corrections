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
*			Module Name:		ColSlat.c
*			Revision Number:	2.3
*			Date last revised:	2 January 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	June 30 1999
*
*			Module Overview:	Model slat collimators. Requires planar or dual
*								headed detector module to be on.
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
*			Revision date:		23 January 2012
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
*						- Support for simulate_PET_coincidences_only option
*
*********************************************************************************/

#define COLLIMATOR_SLAT


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
#include "PhoHFile.h"
#include "Detector.h"
#include "DetPlanar.h"
#include "ProdTbl.h"
#include "SubObj.h"
#include "EmisList.h"
#include "PhoHStat.h"
#include "PhoHFile.h"
#include "ColUsr.h"
#include "UNCCollimator.h"
#include "Collimator.h"
#include "ColSlat.h"
#include "phg.h"
#include "PhgBin.h"

#ifdef PHG_DEBUG
static char					colSltErrStr[1024];				/* Storage for creating error strings */
#endif


void	colSltComputeFreePathsToExit(PHG_Position pos,
			PHG_Direction dir, double energy,
			double *fpPtr);
Boolean	colTrackSlat(PHG_TrackingPhoton *photonPtr);
void 			colDHDoRandomPos(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					CollimatedPhotonsTy *colPhotonsPtr);
void			colDHDoOptimalPos(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					CollimatedPhotonsTy *colPhotonsPtr);			

void			colSltProjectToNewPos(PHG_Position pos, 
					PHG_Direction dir, double distance, double frontOfLayer,
					double backOfLayer, PHG_Position *newPos, double *distTraveled,
					Boolean *exitXFront, Boolean *exitXBack, Boolean *exitY, Boolean *exitZ, double axLow, double axHigh);
void	 		colSlatFindAxialSeg(double zPos, LbFourByte curLayer, LbFourByte *theSeg, Boolean *segIsSlat);

/*********************************************************************************
*
*			Name:			ColSlatSPECT
*
*			Summary:		Perform slat collimation.
*
*			Arguments:
*				PHG_Decay			decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton *photons			- The blue photons to be collimated.
*				LbUsFourByte 		numPhotons		- The number of blue photons.
*				CollimatedPhotonsTy *colPhotonsPtr	- The accepted, collimated photons
*
*			Function return: None.
*
*********************************************************************************/
void ColSlatSPECT(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
		CollimatedPhotonsTy *colPhotonsPtr)
{
	LbUsFourByte	pIndex;					/* Index for current photon */
	double 			detPos;
	
	/* Clear the counters */
	colPhotonsPtr->NumCollimatedBluePhotons = 0;
	colPhotonsPtr->NumCollimatedPinkPhotons = 0;

	/* Bolt if there is nothing to do */
	if (numPhotons == 0) {
		return;
	}
	
	do { /* Process Loop */
	
		/* Get a random detector position for the collimator */
		DetGetRandomDetPosition(&detPos);
		
		/* Loop through all photons */
		for (pIndex = 0; pIndex < numPhotons; pIndex++) {

			/* Get the actual detector/photon incident angle */
			DetGetDetAngle(&photons[pIndex].location, &photons[pIndex].angle,
				detPos, &photons[pIndex].detectorAngle);

			/* Let user modify and/or reject photons */
			if (ColUsrModSPECTPhotonsFPtr && 
					(*ColUsrModSPECTPhotonsFPtr)(&ColRunTimeParams[ColCurParams], decayPtr,
						&(photons[pIndex])) == false) {
				
				/* They rejected it so go to next photon */
				continue;
			}
			
			/* Track through the slat collimator */
			if (colTrackSlat(&photons[pIndex])) {
			
				/* Save the photon */
				colPhotonsPtr->CollimatedTrkngBluePhotons[colPhotonsPtr->NumCollimatedBluePhotons]
					= photons[pIndex];

				/* Increment the counters */
				colPhotonsPtr->NumCollimatedBluePhotons++;
				
				if (photons[pIndex].scatters_in_col == 0) {
					if (photons[pIndex].num_of_scatters > 0) {
						colData[ColCurParams].colTotScatWtPassedThroughCollimator += 
							photons[pIndex].photon_current_weight * decayPtr->startWeight;
					}
					else {
						colData[ColCurParams].colTotPrimWtPassedThroughCollimator += 
							photons[pIndex].photon_current_weight * decayPtr->startWeight;
					}	
				}

				/* Update our collimated photon block if doing history file */
				if  (COL_IsDoHistory()) {
					ColUpdateCollimatedPhotonBlock(&photons[pIndex]);
				}
			}
		}

	} while (false);
}

/*********************************************************************************
*
*			Name:			colSltProjectToNewPos
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
*				double				axLow
*				double				axHigh
*
*			Function return: None.
*
*********************************************************************************/
void colSltProjectToNewPos(PHG_Position pos, 
		PHG_Direction dir, double distance, double frontOfLayer,
		double backOfLayer, PHG_Position *newPos, double *distTraveled,
		Boolean *exitXFront, Boolean *exitXBack, Boolean *exitY, Boolean *exitZ,
		double axLow, double axHigh)
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
		distZ = (axHigh - pos.z_position)/dir.cosine_z;
	}
	else if (dir.cosine_z < 0.0) {
		distZ = (axLow - pos.z_position)/dir.cosine_z;
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
		final position is at the edge of the collimator
	*/
	newPos->x_position = pos.x_position + (*distTraveled * dir.cosine_x);
	newPos->y_position = pos.y_position + (*distTraveled * dir.cosine_y);
	newPos->z_position = pos.z_position + (*distTraveled * dir.cosine_z);
}

/*********************************************************************************
*
*			Name:			colTrackSlat
*
*			Summary:		Track through the slat collimator
*			Arguments:
*				PHG_TrackingPhoton *photons			- The blue photons to collimate.
*				LbUsFourByte 		numPhotons		- The number of blue photons.
*				CollimatedPhotonsTy *colPhotonsPtr	- The accepted, collimated photons
*
*			Function return: True if collimation occurs.
*
*********************************************************************************/
Boolean colTrackSlat(PHG_TrackingPhoton *photonPtr)
{
	Boolean 					photonCollimated = false;	/* Successful collimation or not? */
	Boolean						segIsSlat;					/* Flags position in slat or in gap */
	Boolean						exitXFront;				/* Flags exit through front of layer */
	Boolean						exitXBack;				/* Flags exit through back of layer */
	Boolean						exitY;					/* Flags exit through y limit */
	Boolean						exitZ;					/* Flags exit through z limit */
	LbFourByte					curInteraction = -1;	/* Index for current interaction */
	double						t;						/* Projection distance */
	PHG_Position				pos;					/* The photon position */
	PHG_Position				newPos;					/* The new position after projection */
	PHG_Direction				dir;					/* The photon direction of travel */
	double						distance;				/* Distance to travel, used more than once */
	double						distTraveled;			/* Truncated distance for boundaries */
	double						totDistTraveledInCol;	/* A debugging variable */
	double						totDistTraveledInSlat;	/* A debugging variable */
	double						totDistTraveledInGap;	/* A debugging variable */
	double						attenuation;			/* Attenuation of crystal at current energy */
	double						fpToGo;					/* Free paths to travel in crystal */
	double						comptonToScatterProbability;		/* Ratio of prob of compton to prob of scatter */
	double						scatterProbability;		/* Probability of a scatter */
	double						interactionProbability;	/* Interaction type probability */
	double						backOfLayer = 0.0;		/* "Back" of the block in local coordinates */
	double						frontOfLayer = 0.0;		/* "Front" of the block in local coordinates */
	double						weight;					/* Incoming weight */
	
	LbFourByte		curLayer, curSeg;
	
do { /* Process Loop */

	/* Because we use the weight so much we use a local variable */
	weight = photonPtr->photon_current_weight;

	/* Get the current position/direction to have a working copy */
	pos = photonPtr->location;
	dir = photonPtr->angle;
	totDistTraveledInCol = 0.0;
	totDistTraveledInSlat = 0.0;
	totDistTraveledInGap = 0.0;
	curLayer = 0;
	
	/* Rotate the position/direction into collimator coordinates */
	{
		pos.x_position = (photonPtr->location.x_position *
			PHGMATH_Cosine(photonPtr->detectorAngle)) +
			(photonPtr->location.y_position * PHGMATH_Sine(photonPtr->detectorAngle)) -
			ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].InnerRadius;

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

	/* If photon does not intersect the collimator, finished */
	if (dir.cosine_x < 0.0)
		break;
		
	/*	Now compute how far it is from the back of the object to the face
		of the collimator
	*/

	t = (-pos.x_position/dir.cosine_x);

	/* Project to the face of the collimator */
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
		PhgAbort("Invalid positioning of photon onto collimator face (colTrackSlat)", false);
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

	/* Increment count of photons that reach the collimator */
	colData[ColCurParams].colTotReachingCollimator++;
	
	/*	Sample free paths to travel */
	PhgMathGetTotalFreePaths(&fpToGo);


	/*	Now track the photon through the collimator until it is absorbed
		or it escapes.
	*/

	/* Clear loop counters */
	curInteraction = -1;
	backOfLayer = ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].Depth;
	frontOfLayer = 0.0;
	colSlatFindAxialSeg(pos.z_position, curLayer, &curSeg, &segIsSlat);

	/* NOTE COORDS ARE NOW BACK AT ZERO = CENTER AXIALLY */
	backOfLayer = ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].Depth;
	frontOfLayer = 0.0;
	
	/* Begin tracking through slat loop */
	do {

		/* Get the attenuation of the slat (or air) at the current energy */
		SubObjGetAttenuationInTomo(colData[ColCurParams].colSlatSegs[curLayer][curSeg].Material,
			photonPtr->energy, &attenuation);

		/* Compute the distance to travel, considering
		   no boundary intersections
		*/
		distance = fpToGo/attenuation;

		/* Project the photon the selected distance */
		colSltProjectToNewPos(pos, dir, distance, frontOfLayer, backOfLayer,
			       &newPos, &distTraveled, &exitXFront, &exitXBack, &exitY, &exitZ,
			       colData[ColCurParams].colSlatSegs[curLayer][curSeg].Start, colData[ColCurParams].colSlatSegs[curLayer][curSeg].End);
		
		/* Increment distance counter */
		totDistTraveledInCol += distTraveled;
		if (segIsSlat)
			totDistTraveledInSlat += distTraveled;
		else
			totDistTraveledInGap += distTraveled;

		
		/* Test to see if we are out of the collimator */
		{
			/* See if the photon is going out the back (a successful collimation) */
			if (exitXBack == true) {

				/* Increment the layer and test for exit */
				curLayer += 1;
				if (curLayer == (LbFourByte)(ColRunTimeParams[ColCurParams].SlatCol.NumLayers)) {
					photonCollimated = true;
					
					/* Update starting position to this location */
					photonPtr->location = newPos;
					photonPtr->angle = dir;
					
					/* Now break to go beyond collimator */
					break;
				}
				
				/* We made it here so we are entering a new layer */
				colSlatFindAxialSeg(newPos.z_position, curLayer, &curSeg, &segIsSlat);
				frontOfLayer = backOfLayer;
				backOfLayer = frontOfLayer + ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].Depth;
				
				/* Update the position */
				pos = newPos;
				
				/* Adjust free paths to go */
				fpToGo -= (distTraveled * attenuation);
					
				/* Go back to top and track some more */
				if (fpToGo > 0)
					continue;	  
			}

			/* See if the photon has "turned around" and headed out the face of the layer */
			if (exitXFront == true) {

				/* Test for exit */
				if (curLayer == 0) {

					/* Now break to go beyond collimator */
					goto PHOTON_REJECTED;
				}

				/* Decrement layer */
				curLayer -= 1;
				
				/* We made it here so we are entering a new layer */
				colSlatFindAxialSeg(newPos.z_position, curLayer, &curSeg, &segIsSlat);
				backOfLayer = frontOfLayer;
				frontOfLayer = backOfLayer - ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].Depth;
				
				/* Adjust free paths to go */
				fpToGo -= (distTraveled * attenuation);
				
				/* Update the position */
				pos = newPos;
					
				/* Go back to top and track some more */
				if (fpToGo > 0)
					continue;	  
			}


			/* Check y dimension */
			if (exitY == true) {

				/* Break out */
				goto PHOTON_REJECTED;
			}

			/* Check to see if photon has changed 'layers' in axial direction */
			if (exitZ == true) {

				/* Check to see if up against axial limits */
				if((newPos.z_position <= ColRunTimeParams[ColCurParams].SlatCol.MinZ)|| (newPos.z_position >= ColRunTimeParams[ColCurParams].SlatCol.MaxZ) ||
					(PhgMathRealNumAreEqual(newPos.z_position, ColRunTimeParams[ColCurParams].SlatCol.MinZ, -5, 0, 0, 0) == true) ||
					(PhgMathRealNumAreEqual(newPos.z_position, ColRunTimeParams[ColCurParams].SlatCol.MaxZ, -5, 0, 0, 0) == true)){

					/* Break out */
					goto PHOTON_REJECTED;
				}
				else {
					/* HAVE SWITCHED FROM GAP TO SLAT OR SLAT TO GAP */
					/* NEED TO FIND NEW COORDS FOR BEGIN AND END OF SLAT/GAP */
					/* Update starting position to this location */
					pos = newPos;

					/* SDW add a nudge factor to avoid roundoff error */
					if(dir.cosine_z < 0) {
						curSeg -= 1;
						#ifdef PHG_DEBUG
							if (curSeg < 0) {
								LbInPrintf("Hmmm...");
							}
						#endif
						pos.z_position -= 0.000001;
					}
					else {
						pos.z_position += 0.000001;
						curSeg += 1;
					}
					
					/* Adjust free paths to go */
					fpToGo -= (distTraveled * attenuation);
						
					if (fpToGo > 0)
						continue;	  
				}
			}
		} /* If we run past here we are not out of the collimator */
		
		/* We are here so do an interaction */

		/* Increment count of active interactions */
		curInteraction++;
	

		/* Get probability of scatter */
		comptonToScatterProbability = SubObjGetProbComptToScatInTomo2(colData[ColCurParams].colSlatSegs[curLayer][curSeg].Material, photonPtr->energy);
		scatterProbability = SubObjGetProbScatterInTomo2(colData[ColCurParams].colSlatSegs[curLayer][curSeg].Material, photonPtr->energy);
		interactionProbability = PhgMathGetRandomNumber();
		
		/* See if photon is absorbed */
		if (interactionProbability > scatterProbability) {

			photonPtr->energy = 0.0;

			/* Break out of tracking loop */
			goto PHOTON_REJECTED;
		}

		/* We made it to here, so compute scatter angle and new energy */
		{
			/* Compute the interaction */
			{
				/* First update the photon's direction to the rotated coordinates */
				photonPtr->angle = dir;

				if (interactionProbability > (scatterProbability * comptonToScatterProbability	)) {
					EmisListDoCoherent(photonPtr, colData[ColCurParams].colSlatSegs[curLayer][curSeg].Material);
				}
				else {		
					/* Model Compton scatter (Using Klein-Nishina) */
					EmisListDoComptonInteraction(photonPtr);
				}
				/* Now update the local, working copy of the direction */
				dir = photonPtr->angle;
			}

			/* Update scatter count */
			photonPtr->scatters_in_col += 1;

			/* If the photon energy is below supported minimum energy, force an absorption */
			if (photonPtr->energy < PhgRunTimeParams.PhgMinimumEnergy) {

				/* Break out */
				goto PHOTON_REJECTED;
			}

			/* Update the position */
			pos = newPos;

		} /* End of "We are here so do an interaction" */

		/* Get the next free paths to travel */
		/* LbInPrintf("Update fp to travel...\n");*/
		PhgMathGetTotalFreePaths(&fpToGo);

	} while (true);	/* End of loop for tracking through slats */
} while (false);	/* End of process loop */
	
	PHOTON_REJECTED:;
	return (photonCollimated);
}

/*********************************************************************************
*
*			Name:			ColSlatDualHeaded
*
*			Summary:		Perform Slat Collimation for dual headed detectors
*
*			Arguments:
*				PHG_Decay			decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton *bluePhotons		- The blue photons to collimate.
*				LbUsFourByte 		numBluePhotons	- The number of blue photons.
*				PHG_TrackingPhoton *pinkPhotons		- The  pink photons to collimate.
*				LbUsFourByte		numPinkPhotons	- The number of pink photons.
*				CollimatedPhotonsTy *colPhotonsPtr	- The accepted, collimated photons
*
*			Function return: None.
*
*********************************************************************************/
void ColSlatDualHeaded(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		CollimatedPhotonsTy *colPhotonsPtr)
{

	/* Clear the counters */
	colPhotonsPtr->NumCollimatedBluePhotons = 0;
	colPhotonsPtr->NumCollimatedPinkPhotons = 0;	


	colDHDoRandomPos(decayPtr, bluePhotons, numBluePhotons, pinkPhotons, numPinkPhotons, colPhotonsPtr);


}

/*********************************************************************************
*
*			Name:			ColSlatSetParamsFromDetector
*
*			Summary:		Initializes internal parameters that come from the
*							detector.
*
*			Arguments:
*				double				minAngle		- The starting angle of rotation.
*				double				minZ			- The minimum Z value
*				double				maxZ			- The maximum Z value
*
*			Function return: None.
*
*********************************************************************************/
void ColSlatSetParamsFromDetector(double minAngle, double minZ, double maxZ)
{
	LbUsFourByte	i,j;					/* Loop variables */
	LbUsFourByte	thisSlat, thisSegment;	/* Temps for computing segments */
	double			lastSegmentEnd;			/* Temp for computing segments */
	
	/* Set parameters which come from minimum angle */
	ColRunTimeParams[ColCurParams].SlatCol.MinAngle = minAngle;
	ColRunTimeParams[ColCurParams].SlatCol.MinZ = minZ;
	ColRunTimeParams[ColCurParams].SlatCol.MaxZ = maxZ;

	

	/* Now allocate the segment table */
	if ((colData[ColCurParams].colSlatSegs = (Col_Slat_Slat_Ty **) LbMmAlloc(sizeof(Col_Slat_Slat_Ty *) *
			ColRunTimeParams[ColCurParams].SlatCol.NumLayers)) == 0) {
		
		PhgAbort("Unable to allocate memory for slat table (ColSlatSetParamsFromDetector)", false);
	}
	

	/* Now allocate the segment count table */
	if ((colData[ColCurParams].colSlatNumSegs = (LbFourByte *) LbMmAlloc(sizeof(LbFourByte) *
			ColRunTimeParams[ColCurParams].SlatCol.NumLayers)) == 0) {
		
		PhgAbort("Unable to allocate memory for seg count table (ColSlatSetParamsFromDetector)", false);
	}

	/* Once this is done we can compute our segment table and do some validation */
	/* First count segments and then basically repeat the process, initializing the segments */
	j = 0;
	for (i = 0; i < ColRunTimeParams[ColCurParams].SlatCol.NumLayers; i++) {
		thisSlat = 0;
		thisSegment = 0;
		lastSegmentEnd = ColRunTimeParams[ColCurParams].SlatCol.MinZ;

		do {
			if (ColRunTimeParams[ColCurParams].SlatCol.Layers[i].Slats[thisSlat].Start == lastSegmentEnd) {
				lastSegmentEnd = ColRunTimeParams[ColCurParams].SlatCol.Layers[i].Slats[thisSlat].End;
				thisSlat += 1;
			} 
			else {
				lastSegmentEnd = ColRunTimeParams[ColCurParams].SlatCol.Layers[i].Slats[thisSlat].Start;
			}
			thisSegment += 1;
		} while  (thisSlat < ColRunTimeParams[ColCurParams].SlatCol.Layers[i].NumSlats);
		if (lastSegmentEnd != ColRunTimeParams[ColCurParams].SlatCol.MaxZ) {
			thisSegment += 1;
		}

		colData[ColCurParams].colSlatNumSegs[i] = thisSegment;

		if ((colData[ColCurParams].colSlatSegs[i] = (Col_Slat_Slat_Ty *) LbMmAlloc(sizeof(Col_Slat_Slat_Ty) *
				colData[ColCurParams].colSlatNumSegs[i])) == 0) {
			
			PhgAbort("Unable to allocate memory for slat table (ColSlatSetParamsFromDetector)", false);
		}
	}

	/* Now initialize segments */
	for (i = 0; i < ColRunTimeParams[ColCurParams].SlatCol.NumLayers; i++) {
		thisSegment = 0;
		thisSlat = 0;
		lastSegmentEnd = ColRunTimeParams[ColCurParams].SlatCol.MinZ;

		do {

			if (ColRunTimeParams[ColCurParams].SlatCol.Layers[i].Slats[thisSlat].Start == lastSegmentEnd) {
				colData[ColCurParams].colSlatSegs[i][thisSegment].Start = ColRunTimeParams[ColCurParams].SlatCol.Layers[i].Slats[thisSlat].Start;
				colData[ColCurParams].colSlatSegs[i][thisSegment].End = ColRunTimeParams[ColCurParams].SlatCol.Layers[i].Slats[thisSlat].End;
				colData[ColCurParams].colSlatSegs[i][thisSegment].Material = ColRunTimeParams[ColCurParams].SlatCol.Layers[i].Slats[thisSlat].Material;
				lastSegmentEnd = colData[ColCurParams].colSlatSegs[i][thisSegment].End;
				thisSlat += 1;
			} 
			else {
				colData[ColCurParams].colSlatSegs[i][thisSegment].Start = lastSegmentEnd;
				colData[ColCurParams].colSlatSegs[i][thisSegment].End = ColRunTimeParams[ColCurParams].SlatCol.Layers[i].Slats[thisSlat].Start;
				colData[ColCurParams].colSlatSegs[i][thisSegment].Material = 0;
				lastSegmentEnd = colData[ColCurParams].colSlatSegs[i][thisSegment].End;
			}
			thisSegment += 1;
		} while  (thisSlat < ColRunTimeParams[ColCurParams].SlatCol.Layers[i].NumSlats);
		
		if (lastSegmentEnd != ColRunTimeParams[ColCurParams].SlatCol.MaxZ) {
			colData[ColCurParams].colSlatSegs[i][thisSegment].Start = lastSegmentEnd;
			colData[ColCurParams].colSlatSegs[i][thisSegment].End = ColRunTimeParams[ColCurParams].SlatCol.MaxZ;
			colData[ColCurParams].colSlatSegs[i][thisSegment].Material = 0;
			thisSegment += 1;
		}
	}
	
	/* Now Verify the segments */
	for (i = 0; i < ColRunTimeParams[ColCurParams].SlatCol.NumLayers; i++) {
		for (j = 0; j < (LbUsFourByte)(colData[ColCurParams].colSlatNumSegs[i]); j++) {
			if (colData[ColCurParams].colSlatSegs[i][j].Start >= colData[ColCurParams].colSlatSegs[i][j].End) {
				PhgAbort("Got collimator slat start >= slat end (ColSlatSetParamsFromDetector)", false);
			}
		}
	}
}

/*********************************************************************************
*
*			Name:			colDHDoRandomPos
*
*			Summary:		Perform collimation for dual headed detectors using
*							random detector head position
*
*			Arguments:
*				PHG_Decay			decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton *bluePhotons		- The blue photons to collimate.
*				LbUsFourByte 		numBluePhotons	- The number of blue photons.
*				PHG_TrackingPhoton *pinkPhotons		- The  pink photons to collimate.
*				LbUsFourByte		numPinkPhotons	- The number of pink photons.
*				CollimatedPhotonsTy *colPhotonsPtr	- The accepted, collimated photons
*
*			Function return: None.
*
*********************************************************************************/
void colDHDoRandomPos(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		CollimatedPhotonsTy *colPhotonsPtr)
{

	LbUsFourByte		b,p;
	double				detPos;
	
	do { /* Process Loop */
	
	/* Get a random detector position */
	DetGetRandomDetPosition(&detPos);
	
	/* Loop through the incoming blue and pink photons to determine
	 * the total incoming coincidence weight */
	for (b = 0; b < numBluePhotons; b++) {			
		for (p = 0; p < numPinkPhotons; p++) {			
			colData[ColCurParams].colTotInCoincWt += bluePhotons[b].photon_current_weight
				* pinkPhotons[p].photon_current_weight * decayPtr->startWeight;
		}
	}

	/* Loop through all blue photons */
	for (b = 0; b < numBluePhotons; b++) {
		
		/* increment the incoming blue weight */
		colData[ColCurParams].colTotBluePhotonWt +=  bluePhotons[b].photon_current_weight
			* decayPtr->startWeight;
		
		/* Get the actual detector/photon incident angle */
		DetGetDetAngle(&bluePhotons[b].location, &bluePhotons[b].angle,
			detPos, &bluePhotons[b].detectorAngle);

		/* Let user modify and/or reject photons */
		if (ColUsrModPETPhotonsFPtr && 
				(*ColUsrModPETPhotonsFPtr)(&ColRunTimeParams[ColCurParams], decayPtr,
					&(bluePhotons[b])) == false) {
			
			colData[ColCurParams].colTotRejBluePhotonWt += bluePhotons[b].photon_current_weight
				* decayPtr->startWeight;

			/* They rejected it so go to next photon */
			continue;
		}
		
		/* Track the photon through the slat collimator */
		if (colTrackSlat(&bluePhotons[b]) == false) {

			/* photon rejected, increment rejected weight */
			colData[ColCurParams].colTotRejBluePhotonWt += bluePhotons[b].photon_current_weight
				* decayPtr->startWeight;
			/* go to next photon */
			continue;	
		}
						
		/* We made it to here so save the collimated photons */
		{
			
			colPhotonsPtr->CollimatedTrkngBluePhotons[colPhotonsPtr->NumCollimatedBluePhotons]
				= bluePhotons[b];

			/* Increment the counters */
			colPhotonsPtr->NumCollimatedBluePhotons++;
			
			colData[ColCurParams].colTotAccBluePhotonWt += bluePhotons[b].photon_current_weight
				* decayPtr->startWeight;
				
			if (bluePhotons[b].scatters_in_col == 0) {
				if (bluePhotons[b].num_of_scatters > 0) {
					colData[ColCurParams].colTotScatWtPassedThroughCollimator += 
						bluePhotons[b].photon_current_weight * decayPtr->startWeight;
				}
				else {
					colData[ColCurParams].colTotPrimWtPassedThroughCollimator += 
						bluePhotons[b].photon_current_weight * decayPtr->startWeight;
				}
			}
			
			/* Update our collimated photon block if doing history file */
			if  (COL_IsDoHistory()) {
				ColUpdateCollimatedPhotonBlock(&bluePhotons[b]);
			}
		}
	}
	
	/* if no blue photons were found and this is a coincidence-only
	simulation, break out of the loop */
	if ( PHG_IsPETCoincidencesOnly() && (colPhotonsPtr->NumCollimatedBluePhotons == 0) ) {
		
		break;
		
	}
	
	/* Loop through all pink photons */
	for (p = 0; p < numPinkPhotons; p++) {
		
		/* increment the incoming pink weight */
		colData[ColCurParams].colTotPinkPhotonWt +=  pinkPhotons[p].photon_current_weight
			* decayPtr->startWeight;

		/* Get the actual detector/photon incident angle */
		DetGetDetAngle(&pinkPhotons[p].location, &pinkPhotons[p].angle,
			detPos, &pinkPhotons[p].detectorAngle);

		/* Let user modify and/or reject photons */
		if (ColUsrModPETPhotonsFPtr && 
				(*ColUsrModPETPhotonsFPtr)(&ColRunTimeParams[ColCurParams], decayPtr,
					&(pinkPhotons[p])) == false) {
			
			/* increment rejected pink weight counter */
			colData[ColCurParams].colTotRejPinkPhotonWt += pinkPhotons[p].photon_current_weight
				* decayPtr->startWeight;
			
			/* They rejected it so go to next pink */
			continue;
		}

		/* Track the photon through the slate collimator */
		if (colTrackSlat(&pinkPhotons[p]) == false) {
			
			/* increment rejected pink weight counter */
			colData[ColCurParams].colTotRejPinkPhotonWt += pinkPhotons[p].photon_current_weight
				* decayPtr->startWeight;

			/* rejected so go to next pink */
			continue;	
		}			
	
		/* We made it to here so save the collimated photons */
		{

			colPhotonsPtr->CollimatedTrkngPinkPhotons[colPhotonsPtr->NumCollimatedPinkPhotons]
				= pinkPhotons[p];

			/* Increment the counters */
			colPhotonsPtr->NumCollimatedPinkPhotons++;
		
			/* pink weight clearing collimator */
			colData[ColCurParams].colTotAccPinkPhotonWt += pinkPhotons[p].photon_current_weight
				* decayPtr->startWeight;

			if (pinkPhotons[p].scatters_in_col == 0) {
				/* increment pink weight clearing collimator without interacting */
				if (pinkPhotons[p].num_of_scatters > 0) {
					colData[ColCurParams].colTotScatWtPassedThroughCollimator += 
						pinkPhotons[p].photon_current_weight * decayPtr->startWeight;
				}
				else {
				colData[ColCurParams].colTotPrimWtPassedThroughCollimator += 
					pinkPhotons[p].photon_current_weight * decayPtr->startWeight;
				}
			}

			/* Update our collimated photon block if doing history file */
			if  (COL_IsDoHistory()) {
				ColUpdateCollimatedPhotonBlock(&pinkPhotons[p]);
			}
		}
	}
	
	/* Loop through the outgoing blue and pink photons to determine
	 * the total outgoing collimated coincidence weight */
	for (b = 0; b < colPhotonsPtr->NumCollimatedBluePhotons; b++) {
		for (p = 0; p < colPhotonsPtr->NumCollimatedPinkPhotons; p++) {
			colData[ColCurParams].colTotAccCoincWt += 
				colPhotonsPtr->CollimatedTrkngBluePhotons[b].photon_current_weight * 
				colPhotonsPtr->CollimatedTrkngPinkPhotons[p].photon_current_weight *
				decayPtr->startWeight;
		}
	}


	} while (false);
}

/*********************************************************************************
*
*			Name:			colMCPETFindAxialSeg
*
*			Summary:		Figure out which axial segment the photon is in.
*							
*			Arguments:
*				double			zPos		- The z position of the photon.
*				LbUsFourByte	curLayer	- The current layer.
*				LbFourByte		*theSeg		- Which seg the photon is in.
*				Boolean			*segIsSlat	- If the segment is a slat or air.
*
*			Function return: None.
*
*********************************************************************************/
void colSlatFindAxialSeg(double zPos, LbFourByte curLayer, LbFourByte *theSeg, Boolean *segIsSlat)
{
	LbFourByte		segIndex;		/* The photon's segment */
	
	/* Loop  segments to find the one that contains the photon */
	for (segIndex = 0; segIndex < colData[ColCurParams].colSlatNumSegs[curLayer]; segIndex++) {
	
		if ((colData[ColCurParams].colSlatSegs[curLayer][segIndex].Start <= zPos)
				&&
				(colData[ColCurParams].colSlatSegs[curLayer][segIndex].End >= zPos)) {
				
				*theSeg = segIndex;
				*segIsSlat = (colData[ColCurParams].colSlatSegs[curLayer][segIndex].Material != 0);
			break;
		}
	}
	
	#ifdef PHG_DEBUG
		/* Verify we found a segment */
		if (segIndex == colData[ColCurParams].colSlatNumSegs[curLayer]) {
			sprintf(colSltErrStr, "Photon starting out on collimator segment out of z-bounds (colSlatFindAxialSeg)"
				"\nzPos = %3.5lf, zMin = %3.5lf, zMax = %3.5lf\n",
				zPos, colData[ColCurParams].colSlatSegs[curLayer][0].Start,
				colData[ColCurParams].colSlatSegs[curLayer][segIndex-1].End);
				
			PhgAbort(colSltErrStr, true);
		}
		
	#endif
}



