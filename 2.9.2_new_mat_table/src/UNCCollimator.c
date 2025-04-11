/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1995-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		UNCCollimator.c
*			Revision Number:	1.5
*			Date last revised:	23 July 2013
*			Programmer:			Dave Lewis, Eric Fry, Steven Vannoy
*			Date Originated:	1 January 1995
*
*			Module Overview:	This module provides geometric modelling of
*								a SPECT collimator.  The algorithms come from
*								the work of Dr. Tsui and Dr. Gullberg (see the
*								manual for references).  Much of the code in
*								this module was translated to C and adapted to
*								the PHG by Dave Lewis, using an original
*								FORTRAN implementation by Dr. Frey at the 
*								University of North Carolina.  We would like
*								acknowledge this contribution to the SimSET
*								software package.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:
*				UNCColPrintReport
*				UNCColPrintParams
*				UNCColInitialize
*				UNCCollimate
*				UNCColTerminate
*
*			Global variables defined:		none
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s): Steven Vannoy, Mitch Kaplan
*
*			Revision date: 1/3/96
*
*			Revision description:	The original implementation binned 
*			photons at every detector position at which there was a
*			probability of being detected.  This probability was a function
*			of acceptance angle and the number of views specified.  This 
*			introduced a couple of problems, including the fact that the
*			total weight detected would vary with the number of views, and
*			that  azimuthal angle dimension would be over sampled relative
*			to any other dimension.
*
*			We have revised the algorithm to only bin a photon at a single
*			detector position.  The probability of detection at all possible
*			positions is first computed and then a specific angle for detection
*			is selected from the possibilities.			
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
#define UNC_COLLIMATOR

#include <stdio.h>

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
#include "DetTypes.h"
#include "DetParams.h"
#include "PhoHFile.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "Collimator.h"
#include "UNCCollimator.h"
#include "phg.h"
#include "PhgBin.h"

/* Local types */
typedef struct {
	double				probability;				/* Probability of detection */
	double				transaxialPosition;
	LbUsFourByte		azimuthalAngleIndex;
	double				axialPosition;
	double				detectorAngle;		/* detector angle */
	LbUsFourByte		aaIndex;				/* Index for angle bin (positive integer) */
} colCollimatedParamsTy;
	
/* Local Prototypes */
void grfsetup(void);
double geomrsp(PHG_Position photonPosition, PHG_Direction photonDirection, double detectorAngle, double *y_int, double *z_int);
void xform(PHG_Position *photonPosition, PHG_Direction *photonDirection, double detectorAngle);

/* Global variables */
/*
##
##
	Declare your globals within the PhgUsrBin struct. This will prevent any
	problems with naming and scope.
##
##
*/

/* Maximum number of detector views */
#define	MAX_NUM_VIEWS	500

static char				colErrStr[1024];				/* 	For creating error messages */

static CylPosCylinderTy	colInBoundCyl[PHG_MAX_PARAM_FILES];		/*	The inner bounding cylinder for the collimator */



/*********************************************************************************
*
*			Name:			UNCColPrintParams
*
*			Summary:		Prints the collimation parameters.
*
*			Arguments:
*
*			Function return: None.
*
*********************************************************************************/
void UNCColPrintParams()
{

	switch(ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleGeometry) {
		case PARALLEL:
			LbInPrintf("\n\tHole geometry = parallel");
			break;
			
		case FAN:
			LbInPrintf("\n\tHole geometry = fan");
			LbInPrintf("\n\tFocal length = %3.2f", ColRunTimeParams[ColCurParams].UNCSPECTCol.FocalLength);
			break;
			
		case CONE:
			LbInPrintf("\n\tHole geometry = cone");
			LbInPrintf("\n\tFocal length = %3.2f", ColRunTimeParams[ColCurParams].UNCSPECTCol.FocalLength);
			break;
		
		case HoleTy_NULL:
			LbInPrintf("\n\tHole geometry = NULL");
			LbInPrintf("\n\tInvalid vaule -- parameter probably undefined");
			PhgAbort("UNCColPrintParams finds hole geometry set to NULL", true);
			break;
	}
	
	LbInPrintf("\n\tRadius of rotation = \t%3.2f", ColRunTimeParams[ColCurParams].UNCSPECTCol.RadiusOfRotation);
	LbInPrintf("\n\tCollimator thickness = \t%3.2f", ColRunTimeParams[ColCurParams].UNCSPECTCol.Thickness);
	LbInPrintf("\n\tHole radius = \t%3.3e", ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleRadius);
	LbInPrintf("\n\tSeptal thickness = \t%3.3e", ColRunTimeParams[ColCurParams].UNCSPECTCol.SeptalThickness);
	LbInPrintf("\n\tMin Z = \t%3.2f", ColRunTimeParams[ColCurParams].UNCSPECTCol.MinZ);
	LbInPrintf("\n\tMax Z = \t%3.2f", ColRunTimeParams[ColCurParams].UNCSPECTCol.MaxZ);
	LbInPrintf("\n\tStart angle = \t%3.2f", ColRunTimeParams[ColCurParams].UNCSPECTCol.StartAngle);
	LbInPrintf("\n\tStop angle = \t%3.2f", ColRunTimeParams[ColCurParams].UNCSPECTCol.StopAngle);
	LbInPrintf("\n\tNumber of views = \t%d", ColRunTimeParams[ColCurParams].UNCSPECTCol.NumViews);
	LbInPrintf("\n\tAcceptance Angle = \t%3.3e radians", colData[ColCurParams].colAccAngle);
}

/*********************************************************************************
*
*			Name:			UNCColPrintReport
*
*			Summary:		Prints a report of the final statistics.
*
*			Arguments:
*
*			Function return: None.
*
*********************************************************************************/
void UNCColPrintReport()
{
	LbInPrintf("\n\nSum of accepted primary weight = %3.2e\nSum of accepted scatter weight = %3.2e\n",
			colData[ColCurParams].colAccPrimWeightSum,colData[ColCurParams].colAccScatWeightSum);
			
	if ((colData[ColCurParams].colAccPrimWeightSum != 0.0) && (colData[ColCurParams].colAccScatWeightSum != 0.0))
		LbInPrintf("\nScatter-to-primary ratio = %3.2e\n", 
			colData[ColCurParams].colAccScatWeightSum/colData[ColCurParams].colAccPrimWeightSum);

}
/*********************************************************************************
*
*			Name:			UNCColInitialize
*
*			Summary:		Initialize the UNC Collimator module.
*
*			Arguments:
*
*			Function return: TRUE unless an error occurs.
*
*********************************************************************************/
Boolean UNCColInitialize()
{
	Boolean	okay = false;
	char	errStr[1024];
	
	do { /* Process Loop */
		
		/* Verify that our number of views is less than MAX_NUM_VIEWS */
		if (ColRunTimeParams[ColCurParams].UNCSPECTCol.NumViews >= MAX_NUM_VIEWS) {
			sprintf(errStr, "Your number of views, %ld, is greater than the maximum allowed, %d.\n"
				"Either reduce the number of views in the collimator parameters\n"
				"or increase MAX_NUM_VIEWS in the file UNCCollimator.c and rebuild\n",
				(unsigned long)ColRunTimeParams[ColCurParams].UNCSPECTCol.NumViews, MAX_NUM_VIEWS);
			ErStGeneric(errStr);
			break;
		}
		
		/* Clear our counters */
		colData[ColCurParams].colAccPrimWeightSum = 0.0;
		colData[ColCurParams].colAccScatWeightSum = 0.0;
	
		/* Setup our geometry */
		grfsetup();
	
		/* Compute the acceptance angle depending on the hole geometry */
		if (ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleGeometry == PARALLEL) {
		
			colData[ColCurParams].colAccAngle =
				atan((2*ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleRadius)/
				ColRunTimeParams[ColCurParams].UNCSPECTCol.Thickness);
	
		} else {
	
		    colData[ColCurParams].colAccAngle = 
		    	atan((2*ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleRadius)/
		    	ColRunTimeParams[ColCurParams].UNCSPECTCol.Thickness)
				+ atan((ColGtOutsideRadius())/
				ColRunTimeParams[ColCurParams].UNCSPECTCol.FocalLength);
		}
	
		/* Convert angles to radians */
		ColRunTimeParams[ColCurParams].UNCSPECTCol.StartAngle *= PHGMATH_PI/180.0;
		ColRunTimeParams[ColCurParams].UNCSPECTCol.StopAngle *= PHGMATH_PI/180.0;
		
		/* Compute range of detected angles */
		colData[ColCurParams].colRangeOfDetAngles = ColRunTimeParams[ColCurParams].UNCSPECTCol.StopAngle
			- ColRunTimeParams[ColCurParams].UNCSPECTCol.StartAngle;
	
		/* Compute the inner collimator cylinder */
		colInBoundCyl[ColCurParams].radius = ColRunTimeParams[ColCurParams].UNCSPECTCol.RadiusOfRotation;
		colInBoundCyl[ColCurParams].zMin = ColRunTimeParams[ColCurParams].UNCSPECTCol.MinZ;
		colInBoundCyl[ColCurParams].zMax = ColRunTimeParams[ColCurParams].UNCSPECTCol.MaxZ;
		colInBoundCyl[ColCurParams].centerX = 0.0;
		colInBoundCyl[ColCurParams].centerY = 0.0;

		okay = true;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			UNCCollimate
*
*			Summary:		"Track" the photons through the UNC style collimator.
*
*
*			Arguments:
*				PHG_Decay			decayPtr			- The decayPtr that started the process.
*				PHG_TrackingPhoton *photons			- The photons detected.
*				LbUsFourByte 		numPhotons		- The number of blue photons.
*				CollimatedPhotonsTy *colPhotonsPtr	- The accepted collimated photons
*			Function return: None.
*
*********************************************************************************/
void UNCCollimate(PHG_Decay *decayPtr, PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
		CollimatedPhotonsTy *colPhotonsPtr)
{
	LbUsFourByte			i;						/* counter (used for energy array and for atten coeff of crystal */
	double					prob_det_loc;			/* probability that the photon will be detected at location calculated in geomrsp */
	double					z_intercept;			/* vertical coordinate of photon at detector */
	double					angle_of_photon;		/* azimuthal angle of photon */

	double					distance;				/* horizontal coordinate of photon at detector; what they call "transaxial distance" */
	double					angle_of_detector;		/* detector angle */
	/*double					detectionWeight = 0;*/	/* */
	double					angleIndexInitial;		/* Initial index for angle bin looping */
	double					angleIndexFinal;		/* Final index for angle bin looping */
	double					Index;					/* "Index" for angle bin */
	double					angleIndex;				/* "Index" for angle bin */
	double					rSquared;				/* Distance of photon from origin */
	LbUsFourByte			aaIndex;				/* Index for angle bin (positive integer) */
	LbUsFourByte			numColPhotons;			/* Temp var for convenience */
	
	colCollimatedParamsTy	cumProb[MAX_NUM_VIEWS];						/* Cumulated probability of detection */
	LbUsFourByte			cumulativeIndex;							/* Index for cumulative array */
	double					detPosProb;									/* Probability for which detector position to select */
	PHG_Position			newPosition;								/* Projected location */
	Col_UNC_SPECT_Ty		*colParams = &ColRunTimeParams[ColCurParams].UNCSPECTCol;	/* Just for brevity */
	LbUsFourByte			numViews = colParams->NumViews;				/* Just for brevity */
	double					aRandom;

	do {
	
		/* Use our short variable for number of collimated photons */
		numColPhotons = colPhotonsPtr->NumCollimatedBluePhotons;
		
		/* Loop through all photons */
		for (i = 0; i < numPhotons; i++) {
		
			/* The target cylinder of the PHG and the inner-most radius of the
				collimator may be different. Hence, we will first check to
				see if we are on the inner-most surface. If not, we will
				project to there
			*/
			{
				rSquared = PHGMATH_Square(photons[i].location.x_position) +
					PHGMATH_Square(photons[i].location.y_position);
					
			
				if (rSquared < PHGMATH_Square(colParams->RadiusOfRotation)) {
					
					/* Project photon to cylinder */
					if (CylPosProjectToCylinder(&(photons[i].location),
								&(photons[i].angle), &colInBoundCyl[ColCurParams],
						&newPosition, &distance) == false) {
					
						/*	Go to next photon */
						continue;
					}
					
					/* See if we fall outside of axial limits */
					if ((newPosition.z_position > colParams->MaxZ) ||
							(newPosition.z_position < colParams->MinZ)) {
						
						/*	Go to next photon */
						continue;
					}
					
					/* We made it here so update our position on the cylinder */
					photons[i].location = newPosition;
					photons[i].travel_distance += distance;
				}
			}
		
			/* Compute azimuthal angle */
			angle_of_photon = atan2(photons[i].angle.cosine_y, photons[i].angle.cosine_x);
			if (angle_of_photon < 0.0) 	/* Put photon angle in range (0,2PI) */
				angle_of_photon += PHGMATH_2PI;
	
			/* From azimuthal angle of photon and acceptance angle of collimator, calculate the detector positions */
			/* where the photon has a chance of being detected. Then loop over these positions for binning. */
			/* Note that most of the indexes corresponding to detector angles are double, not integer. */
			
			/* First compute initial and final angle indexes for angle looping */
			if (colParams->NumViews != 1) {
			
				/* More than one view (normal), compute first angle index */
				angleIndexInitial =
					((angle_of_photon - colData[ColCurParams].colAccAngle - colParams->StartAngle) *
					colParams->NumViews)/colData[ColCurParams].colRangeOfDetAngles;
	
				/* Compute final angle index */
				angleIndexFinal =
					((angle_of_photon + colData[ColCurParams].colAccAngle - colParams->StartAngle) *
					colParams->NumViews)/colData[ColCurParams].colRangeOfDetAngles;
			
			} else {
				/* Single view */
				angleIndexInitial = 0.0;
				angleIndexFinal = 0.0;
			}
				
			
			/* Clear index for cumulative probabilities */
			cumulativeIndex = 0;
			cumProb[0].probability = 0.0;

			/* Loop over viable detector angles, note that these are doubles */
			for (Index = angleIndexInitial; Index <= angleIndexFinal; Index++) {
		
				/* Put into temp for transformation */
				angleIndex = Index;
				
				/* See if we are within range minimum, "continue" to next angle if not */
				if (angleIndex <= -1.0) {
					angleIndex += (PHGMATH_2PI * numViews) / colData[ColCurParams].colRangeOfDetAngles;
					
					/* If photon out of range continue to next */
					if (angleIndex <= -1.0 || angleIndex > numViews-1.)
						continue;
				}
				
				
				/* See if we are within range maximum, "continue" to next angle if not */
				if (angleIndex > (numViews - 1.0)) {
					angleIndex -= (PHGMATH_2PI * numViews) / colData[ColCurParams].colRangeOfDetAngles;
	
					/* If photon out of range continue to next */
					if (angleIndex <= -1.0 || angleIndex > numViews-1.0)
						continue;
				}
		
				/* Convert index to integer note, this gets saved for the binning process */
				aaIndex = (LbFourByte) ceil(angleIndex);
				
				/* Compute the angle of the collimator (detector for now) */
				angle_of_detector = colParams->StartAngle + ((colData[ColCurParams].colRangeOfDetAngles * aaIndex) / numViews);


				/* Compute where photon hits detector, and compute value of psf */  			
				prob_det_loc = geomrsp(photons[i].location, photons[i].angle, angle_of_detector, &distance, &z_intercept);
				
				/* If probability is -1 it missed the detector and we ignore it */
				if (prob_det_loc == -1.0)
					continue;
		
				/* See if outside Z range, continue to next photon if so */
				if ((z_intercept < colParams->MinZ) || (z_intercept > colParams->MaxZ))
					continue;
				
				/* Increment the cumulative probability */
				if (cumulativeIndex == 0)
					cumProb[cumulativeIndex].probability = prob_det_loc;
				else
					cumProb[cumulativeIndex].probability = cumProb[cumulativeIndex-1].probability + prob_det_loc;
				
										
				/* Save the collimator dependent parameters */
				{
					/* Save the transaxial position (called distance for poor reasons on my part) */
					cumProb[cumulativeIndex].transaxialPosition = distance;

					/* Save the axial position */
					cumProb[cumulativeIndex].axialPosition = z_intercept;
					
					/* Save the axial angleIndex */
					cumProb[cumulativeIndex].azimuthalAngleIndex = aaIndex;
						
					/* Save the angle of the detector */
					cumProb[cumulativeIndex].detectorAngle = angle_of_detector;
				}
					
				/* Increment cumulative index */
				cumulativeIndex++;
							
			}	/* End of for-loop for angles */
			
			
			/* If probability is zero skip the photon */
			if ((cumulativeIndex == 0) || (cumProb[cumulativeIndex-1].probability == 0.0))
				continue;
							
			/* Pick a random number for selecting the lucky detector position */
			aRandom = PhgMathGetRandomNumber();
			detPosProb = aRandom * cumProb[cumulativeIndex-1].probability;
			
			/* Loop through detector positions to find the lucky one */
			for (aaIndex = 0; aaIndex < cumulativeIndex; aaIndex++){
			
				/* If our lucky probability is less than this bin's probability, we are done */
				if (detPosProb <= cumProb[aaIndex].probability)
					break;
			}
			
			/* If we didn't find a probability , something is wrong */
			#ifdef PHG_DEBUG
				if (aaIndex == cumulativeIndex) {
					sprintf(colErrStr,"\nFailed to find a cumulative probability (UNCCollimate)\n"
						"detPosProb = %3.2f, cumulativeIndex = %ld, cumProb[cumulativeIndex-1] = %3.2f\n",
						detPosProb, (unsigned long)cumulativeIndex, 
						cumProb[cumulativeIndex-1].probability);
					PhgAbort(colErrStr, true);
				}
			#endif

			/*  Update the collimated photons table */
			{
				/* Add this photon to the table */	
				colPhotonsPtr->CollimatedTrkngBluePhotons[numColPhotons]
					= photons[i];

				/* Adjust the weight and update statistics */
				if (colPhotonsPtr->CollimatedTrkngBluePhotons[numColPhotons].num_of_scatters > 0) {
					
					/* Adjust the weight as a function of cumulative probability and number views */
					colPhotonsPtr->CollimatedTrkngBluePhotons[numColPhotons].photon_current_weight *= 
						(cumProb[cumulativeIndex-1].probability/numViews);
						
					/* Sum the detected weight, accounting for the decay weight */
					colData[ColCurParams].colAccScatWeightSum += 
						(colPhotonsPtr->CollimatedTrkngBluePhotons[numColPhotons].photon_current_weight *
						decayPtr->startWeight);
				}
				else {

					/* Adjust the weight as a function of cumulative probability and number views */
					colPhotonsPtr->CollimatedTrkngBluePhotons[numColPhotons].photon_current_weight *=
						(cumProb[cumulativeIndex-1].probability/numViews);

					/* Sum the detected weight, accounting for the decay weight */
					colData[ColCurParams].colAccPrimWeightSum += 
						(colPhotonsPtr->CollimatedTrkngBluePhotons[numColPhotons].photon_current_weight *
						decayPtr->startWeight);
				}
				
				/* Project the photon to the back of the collimator */
				{
					newPosition.x_position = photons[i].location.x_position +
						(photons[i].angle.cosine_x * colParams->Thickness); 
	
					newPosition.y_position = photons[i].location.y_position +
						(photons[i].angle.cosine_y * colParams->Thickness); 
	
					newPosition.z_position = photons[i].location.z_position +
						(photons[i].angle.cosine_z * colParams->Thickness); 
				}

				/* Set the current position to the projected location */
				colPhotonsPtr->CollimatedTrkngBluePhotons[numColPhotons].location = newPosition;
				
				/* Save the transaxial position (called distance for poor reasons on my part) */
				colPhotonsPtr->CollimatedTrkngBluePhotons[numColPhotons].transaxialPosition =
					cumProb[aaIndex].transaxialPosition;

				colPhotonsPtr->CollimatedTrkngBluePhotons[numColPhotons].axialPosition = 
					cumProb[aaIndex].axialPosition;
				
				/* Save the axial angleIndex, the debug check */
				colPhotonsPtr->CollimatedTrkngBluePhotons[numColPhotons].azimuthalAngleIndex =
					cumProb[aaIndex].azimuthalAngleIndex;
					
				/* Save the angle of the detector */
				colPhotonsPtr->CollimatedTrkngBluePhotons[numColPhotons].detectorAngle =
					cumProb[aaIndex].detectorAngle;

				#ifdef PHG_DEBUG
					if (numColPhotons == PHG_MAX_DETECTED_PHOTONS)
						PhgAbort("\nToo many collimated photons (UNCCollimate).\n", true);
				#endif
				
				/* Update our detected photon block if doing history file */
				if  (COL_IsDoHistory()) {
					ColUpdateCollimatedPhotonBlock(
						&colPhotonsPtr->CollimatedTrkngBluePhotons[numColPhotons]);
				}

				/* Increment the counters */
				numColPhotons++;
			}

		} /* End of for-loop for each photon */
		
		/* Update count of collimated photons */
		colPhotonsPtr->NumCollimatedBluePhotons = numColPhotons;
		
	} while (false);
}


/*********************************************************************************
*
*			Name:			UNCColTerminate
*
*			Summary:		Performing any final data processing.
*
*			Arguments:
*
*			Function return: None.
*
*********************************************************************************/
void UNCColTerminate()
{
}

/************

	Steven Vannoy
	1/3/95
	
	
	The rest of this file contains code that originated at The University of North
	Carolina. It is primarily a product of Dr. Eric Frey, and Dave Lewis. We have
	left it in its original form until it can be transformed into the standard SimSET
	format. Note, local global variables have been changed to SimSET format so the
	code below is not exactly as we received.
	
	Changes that have been made to the function of the code have been the elimination
	of detector effects. These issues are dealt with in the detector module. In some
	cases the detector specific code may have been eliminated, and in some cases it
	may have been neutralized by setting specific parameters to neutral values.
*************/

/*************************************************************************
*
*   The following was taken from Detector.c. I am including it here
*   so that we won't have to change the make file(s) in order to compile
*   another file.
*
*************************************************************************/

/***************************************************************************************
* Dependencies: global variables are set up by grfsetup
* Author:  Dave Lewis
* Date: 2/94
* Latest Revision:
*    modified from geomrsp.f by Eric Frey (in ~frey/src/simind)
*
* Procedure:  
*	Calculates the value of the point spread function and the intercept with the collimator backplane.
*	It is assumed that the ray connecting the focal point of the cone beam collimator and the position (0,0) 
*	on the detector passes through the origin of the object volume.
*
* Outputs: 
*      probability that the photon would be detected.
*
* Notes:
*	The projected area seen by the source is r^2(theta-sin(theta)). Divide this by the area of a unit cell.
*	To find sin(theta) we use sin(theta) = 2 sin(theta/2)*cos(theta/2)
*		cos_half_theta is cos(theta/2) in notation of larry's note
*		sin_half_theta = sin(theta/2) = sqrt( 1-cos(theta/2)^2 )
*	In the calculation of rt, the factor 1. / (colDistOriginToColBack - x0) is 1/(z+L+B) from the tsui/metz papers
*	colCellUnitArea = area of collimator unit for circular holes in a hexagonal closed pack arrangement
*
***************************************************************************************/

double geomrsp(PHG_Position photonPosition, PHG_Direction photonDirection, double detectorAngle, double *y_int, double *z_int)
{
	double	x0, y0, z0;				/* Initial position of photon on target cylinder, in transformed coordinates */
	double	cos_x, cos_y, cos_z;	/* Direction cosines describing direction of photon, in transformed coordinates */
	double	cos_half_theta;
	double	sin_half_theta;
	double	x_dist_to_coll;
	double	rt, rty, rtz;			/* Distance between projected position of holes projected on detection plane */
	double	weight;

	/* Transform coordinate systems--this routine (geomrsp) assumes that the collimator is perpendicular to the x-axis */
	xform(&photonPosition, &photonDirection, detectorAngle);

	/* Copy position */
	x0=photonPosition.x_position;
	y0=photonPosition.y_position;
	z0=photonPosition.z_position;
	
	/* Copy angle */
	cos_x=photonDirection.cosine_x;
	cos_y=photonDirection.cosine_y;
	cos_z=photonDirection.cosine_z;

	/* Clear the weight value */
	weight = 0.0;
	
	do {
		
		/* Check to see if our direction is out of bounds */
		if (cos_x < 1.0e-5) {
		
			/* Set weight to flag value and exit */
			weight = -1.;
			break;
		}

		/* Compute distance and z intercept */
		*y_int = cos_y/cos_x*(colData[ColCurParams].colDistOriginToColBack-x0)+y0;
		*z_int = cos_z/cos_x*(colData[ColCurParams].colDistOriginToColBack-x0)+z0;

		x_dist_to_coll = ColRunTimeParams[ColCurParams].UNCSPECTCol.RadiusOfRotation - x0;

		rty=(colData[ColCurParams].k1y - colData[ColCurParams].k2y * x_dist_to_coll)* (*y_int) - colData[ColCurParams].k3y*y0;
		rtz=(colData[ColCurParams].k1z - colData[ColCurParams].k2z * x_dist_to_coll)* (*z_int) - colData[ColCurParams].k3z*z0;
		
		errno = 0;		/* ##rh These changes are to help trap a bug report that I have been unable to reproduce */

		/*	The k* parameters are global variables and were calculated by grfsetup */
		#ifdef WINNT
			rt = _hypot(rty,rtz) / (colData[ColCurParams].colDistOriginToColBack - x0);
		#else
			rt = hypot(rty,rtz) / (colData[ColCurParams].colDistOriginToColBack - x0);
		#endif
		
		if (errno != 0) {		/* ##rh */
			sprintf(colErrStr,"\nerrno non-zero, %d (UNCCollimate geomrsp 1)\n", errno);
			PhgAbort(colErrStr, true);
		}

		cos_half_theta = rt / (2.0*ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleRadius);
		if ( fabs(cos_half_theta) > 1.0)	/* Exit with weight of 0 */		/* ##rh */
			break;

		sin_half_theta = PHGMATH_SquareRoot(1.0 - PHGMATH_Square(cos_half_theta));
		
		if (errno != 0) {		/* ##rh */
			sprintf(colErrStr,"\nerrno non-zero, %d (UNCCollimate geomrsp 2)\n", errno);
			PhgAbort(colErrStr, true);
		}

		weight = PHGMATH_Square(ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleRadius) *
			(2 * PHGMATH_ArcCosine(cos_half_theta) - 2 * cos_half_theta * sin_half_theta)/
			colData[ColCurParams].colCellUnitArea;		
			
		if (errno != 0) {		/* ##rh */
			sprintf(colErrStr,"\nerrno non-zero, %d (UNCCollimate geomrsp 3)\n", errno);
			PhgAbort(colErrStr, true);
		}

	} while (false);

	return(weight);
}

/*****************************************************************************************************************
*     grfsetup: this subroutine calculates various constants used by the routine geomrsp from input parameters

	MODIFIED 1/7/95 SDV
		I removed the uses of mu_crystal, and the dependencies on the existence of the detector. The constants
		are now computed with the back of the collimator as the target, not the face of the detector.
		
*****************************************************************************************************************/

void grfsetup()
{
	
	Col_UNC_SPECT_Ty	*colParams = &ColRunTimeParams[ColCurParams].UNCSPECTCol;	/* Just for brevity */


	/* Compute module globals */
	colData[ColCurParams].colDistOriginToColBack = colParams->RadiusOfRotation+colParams->Thickness;
	colData[ColCurParams].colCellUnitArea = 2*PHGMATH_SquareRoot(3.0) *
		PHGMATH_Square(colParams->HoleRadius + colParams->SeptalThickness);

	switch (colParams->HoleGeometry) {
	
	case PARALLEL: 		/*     parallel */
		colData[ColCurParams].k1y = colData[ColCurParams].k1z = colParams->Thickness;
		colData[ColCurParams].k2y = colData[ColCurParams].k2z = 0;
		colData[ColCurParams].k3y = colData[ColCurParams].k3z = colParams->Thickness;
		break;
		
	case FAN: 		/*     fan */
		colData[ColCurParams].k1z = colParams->Thickness;
		colData[ColCurParams].k2z = 0;
		colData[ColCurParams].k3z = colParams->Thickness;
		
		colData[ColCurParams].k1y = colParams->FocalLength * ColRunTimeParams[ColCurParams].UNCSPECTCol.Thickness / 
			(colParams->Thickness+colParams->FocalLength);
			
		colData[ColCurParams].k2y = colParams->Thickness/(colParams->FocalLength+colParams->Thickness);
		
		colData[ColCurParams].k3y = colParams->Thickness * (colParams->Thickness+colParams->FocalLength)/
			(colParams->FocalLength+colParams->Thickness);
		break;
		
	case CONE: 		/*     cone */
		colData[ColCurParams].k1y = colData[ColCurParams].k1z = colParams->FocalLength
			*colParams->Thickness/(colParams->Thickness + colParams->FocalLength);
			
		colData[ColCurParams].k2y = colData[ColCurParams].k2z = colParams->Thickness/(colParams->FocalLength+colParams->Thickness);
		
		colData[ColCurParams].k3y = colData[ColCurParams].k3z = colParams->Thickness*(colParams->Thickness+colParams->FocalLength)/
			(colParams->FocalLength+colParams->Thickness);
		break;
		
	default:
		ErAbort("Exiting: illegal collimator type in grfsetup\n");
	}
}

/*******************************************************************************************************
*   xform:
*	This function rotates the coordinate system so the x-axis is perpendicular to the detector. 
*	Then it transforms the coordinates of the energy vector. This makes other calculations easier. 
*	In this geometry, z is the axis of rotation, x is perpendicular to the collimator, and y is parallel to 
*	the collimator.
*
*******************************************************************************************************/
void xform(PHG_Position *photonPosition, PHG_Direction *photonDirection, double detectorAngle)
{
	double xtemp, ytemp;
	
	errno = 0;		/* ##rh */
	
	xtemp =   photonPosition->x_position * cos(detectorAngle) + photonPosition->y_position * sin(detectorAngle);
	ytemp = - photonPosition->x_position * sin(detectorAngle) + photonPosition->y_position * cos(detectorAngle);
	photonPosition->x_position = xtemp;
	photonPosition->y_position = ytemp;

	xtemp =   photonDirection->cosine_x * cos(detectorAngle) + photonDirection->cosine_y * sin(detectorAngle);
	ytemp = - photonDirection->cosine_x * sin(detectorAngle) + photonDirection->cosine_y * cos(detectorAngle);
	photonDirection->cosine_x = xtemp;
	photonDirection->cosine_y = ytemp;
	
}

#undef	MAX_NUM_VIEWS	
