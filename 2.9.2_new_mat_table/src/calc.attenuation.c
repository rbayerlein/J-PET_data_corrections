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
*			Module Name:		calc.attenuation.c
*			Revision Number:	1.5
*			Date last revised:	6 June 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	2 May, 1995
*
*			Module Overview:	Computes the "attenuation" of an attenuation image.
*
*			References:			None.
*
**********************************************************************************
*
*			Global functions defined:
*
*			Global macros defined:
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

#define PHG_CALC_ATTEN_MAIN

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbFile.h"
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
#include "CylPos.h"
#include "PhgMath.h"
#include "PhoHFile.h"
#include "PhgHdr.h"
#include "ProdTbl.h"
#include "PhoTrk.h"
#include "SubObj.h"
#include "EmisList.h"
#include "Collimator.h"
#include "Detector.h"
#include "phg.h"
#include "PhgBin.h"

#include "calc.attenuation.h"

#ifdef MPW
	#pragma segment PHG_CALC_ATTEN_MAIN
#endif

/* LOCAL CONSTANTS */
#define	CLCATN_NUM_DIMENIONS	4				/* Number of dimensions supported */

/* LOCAL TYPES */
typedef	struct {	
	PHG_BinEnDimensionsTy	whichDimension;		/* Which dimension is represented */
	LbFourByte				numBins;			/* Number of bins for dimension */
	LbFourByte				dimIndex;			/* Current index for this dimension */
	LbFourByte				*virtualIndexPtr;	/* Reference to the non-subsampled index */
	double					*dimValuePtr;		/* Reference to dimension value, to be reset to zero after each iteration */
} clcAtnDimensionInfoTy;

/* LOCAL GLOBALS */
static Boolean	canceled;								/* Global cancelation flag */
static Boolean	clcAtnIsInitialized = false;			/* Initialization flag */

/* The following local globals are all used due to optimization issues */
#ifdef DO_YOU_REALLY_WANT_THIS
	static	char		clcAtnErrStr[1024];			/* General storage for error strings */
#endif
static	double		clcAtnCurAngle;				/* Angle for current projection */
static	double		clcAtnCurDistance;			/* Distance for current projection */
static	double		clcAtnCurZ1;				/* Z1 for current projection */
static	double		clcAtnCurZ2;				/* Z2 for current projection */
static	double		clcAtnAngleIncr;			/* Amount to increment angle each time */
static	double		clcAtnDistIncr;				/* Amount to increment distance each time */
static	double		clcAtnZIncr;				/* Amount to increment Z values each time */
static	double		clcAtnTargRadius;			/* The target cylinder radius */
static	double		clcAtnObjRadius;			/* The object cylinder radius */
static	double		clcAtnObjZMin;				/* Minimum Z of object cylinder */
static	double		clcAtnObjZMax;				/* Maximum Z of object cylinder */
static	double		clcAtnObjCenterX;			/* Center X of object cylinder */
static	double		clcAtnObjCenterY;			/* Center Y of object cylinder */
static	double		clcAtnAngleRange;			/* Range of angles, varies due to PET/SPECT */
static	LbFourByte	clcAtnCurProjectionNumber;	/* Current projection number */
static	LbFourByte	clcAtnAngleIndex = -1;		/* The current angle index */
static	LbFourByte	clcAtnDistanceIndex = 0;	/* The current distance index */
static	LbFourByte	clcAtnZ1Index = -1;			/* The current Z1 index */
static	LbFourByte	clcAtnZ2Index = -1;			/* The curren Z2 index */
static	void		*clcAtnDAArray = 0;			/* The distance angle data */
static	LbFourByte	clcAtnNumSubBins = 1;		/* Number of sub samples to do for each dimension */
static	LbFourByte	clcAtnNumSamples = 1;		/* Number of samples for each bin */

/* PROTOTYPES */
Boolean			CalcAttenuation(int argc, char *argv[]);
Boolean			clcAtnProjectToCylinder(PHG_Position *posPtr, PHG_Direction *dirPtr,
					double radius, double zMin, double zMax, double centerX,
					double centerY, double *distPtr);
void			clcAtnUpdateDimension(PHG_BinEnDimensionsTy	whichDimension, LbFourByte binNum);
void			clcAtnCalcProjection(Boolean doAverage);
void			clcAtnResetMinimums(PHG_BinEnDimensionsTy whichDimension);

/* FUNCTIONS */

/**********************
*	clcAtnResetMinimums
*
*	Purpose:	Reset dimension minimums after loop execution.
*
*	Result:	Always returns zero indicating no error.
***********************/
void clcAtnResetMinimums(PHG_BinEnDimensionsTy whichDimension)
{
		
	/* Decide which dimension should be updated here */
	switch (whichDimension){
		case PhgBinEn_TD:
			
			/* Compute current distance */
			clcAtnCurDistance = PhgBinParams[0].minTD - (clcAtnDistIncr * .5);
			break;
	
		case PhgBinEn_AA:
		
			/* Compute current angle */
			clcAtnCurAngle = PhgBinParams[0].minAA - (clcAtnAngleIncr * .5);
			break;
		
		case PhgBinEn_Z1:
			clcAtnCurZ1 = PhgBinParams[0].minZ - (clcAtnZIncr * .5);
			break;
			
		case PhgBinEn_Z2:
			clcAtnCurZ2 = PhgBinParams[0].minZ - (clcAtnZIncr * .5);
			break;
		
		/* Ignore other cases */
		case PhgBinEn_Null:
		case PhgBinEn_Energy1:
		case PhgBinEn_Energy2:
		case PhgBinEn_Scatter1:
		case PhgBinEn_Scatter2:
		case PhgBinEn_TOF:
		case PhgBinEn_THETA:
		case PhgBinEn_PHI:
		case PhgBinEn_XR:
		case PhgBinEn_YR:
		case PhgBinEn_Crystal1:
		case PhgBinEn_Crystal2:
			break;
	}
}


/**********************
*	CalcAttenuation
*
*	Purpose:	Execute the program.
*
*	Result:	Always returns zero indicating no error.
***********************/
Boolean CalcAttenuation(int argc, char *argv[])
{
	Boolean			okay = false;			/* Process Loop */
	char			imageFileName[1024];	/* Name of output file */
	FILE			*imageFile = 0;			/* The output image file */



	do { /* Process Loop */
		
		/* Initialize this library */
		if (ClcAtnInitialize(argc, argv) == false)
			break;
		
		/* Get some input values */
		/* Ask for the file name */
		LbInAsk("Enter name of output file", 0, false,
				&canceled, 0, 0, 0, 0,
				imageFileName);
	
		/* Bolt if we canceled */
		if (canceled) {
			ErStCancel("User canceled.");
			goto CANCEL;
		}

		/* Open the output file */
		if ((imageFile = LbFlFileOpen(imageFileName, "w+b")) == 0) {
			ErStFileError("Unable to open image file.");
			break;
		}
		
		/* Compute the attenuation */
		if (ClcAtnCalcAttenuation(&clcAtnDAArray) == false) {
			break;
		}
		
		/* Write out the image file */
		if (fwrite(clcAtnDAArray, PhgBinParams[0].weightImageSize, 1, imageFile) 	!= 1){
			ErStFileError("Unable to write image file.");
			break;
		}
		
		okay = true;
		CANCEL:;
		FAIL:;
	} while (false);
	

		
	/* Handle error situation if one exists */
	if (!okay && canceled) {
		ErHandle("User canceled CalcAttenuation.", false);
		okay = true;
	}

	/* Free local memory */
	if (clcAtnDAArray != 0)
		LbMmFree((void **)&clcAtnDAArray);
	
	/* Close local file */
	if (imageFile != 0)
		fclose(imageFile);
			
	/* Terminate  modules */
	PhoTrkTerminate();
	SubObjTerminate();
	PhgMathTerminate();
	EmisListTerminate();
	

	/* Quit the program */
	return (okay);
}

/**********************
*	ClcAtnCalcAttenuation
*
*	Purpose:	Calculate the attenuation of the attenuation image.
*
*	Arguments:
*		void	**atnAry	- The array data, allocated in this routine
*
*	Result:	true unless an error occurs.
***********************/
Boolean ClcAtnCalcAttenuation(void	**atnAry)
{
	Boolean					okay = false;						/* Process flag */
	Boolean					canceled = false;					/* Cancelation flag */
	LbFourByte				dimIndex;							/* Index for traversing dimensions */
	LbFourByte				index;								/* Generic loop index */
	clcAtnDimensionInfoTy	dimensions[CLCATN_NUM_DIMENIONS];
	
	do { /* Process Loop */
		
		/* See how many sub samples they want */
		clcAtnNumSubBins = LbInAskFourByte("Enter the number of sub-samples for each dimension",
			1, false, false, 1, &canceled, 10, 0, 0, 0);
		
		if (canceled == true)
			goto CANCEL;
	
		/* Verify we are initialized */
		if (clcAtnIsInitialized == false) {
			ErStGeneric("You must call ClcAtnInitialize before any other routine (ClcAtnCalcAttenuation)");
			break;
		}
		
		/* Allocate d/a array */
		if ((*atnAry = LbMmAlloc(PhgBinParams[0].weightImageSize))
				== 0){
			break;
		}
		
		/* Set our local global for convenience */
		clcAtnDAArray = *atnAry;
		
		
		/* Compute geometric values */
		{
			/* Set our radius to the target cylinder */
			clcAtnTargRadius = CylPosGetTargetRadius();
			
			/* Compute angle range, differs due to PET or SPECT and is not in the parameter struct */
			clcAtnAngleRange = PhgBinParams[0].maxAA - PhgBinParams[0].minAA;
			
			
			/* Get object cylinder information */
			SubObjGetObjCylinder(&clcAtnObjRadius, &clcAtnObjZMin,
				&clcAtnObjZMax, &clcAtnObjCenterX, &clcAtnObjCenterY);
		}
		
		/* Setup our dimension variables so that we can sequentially index
			through them while still supporting variable ordering.
		*/
		index = 0;
		for (dimIndex = PHGBIN_NUM_DIMENSIONS-1; dimIndex >= 0; dimIndex--){
		
			switch (PhgBinParams[0].PhgBinDimensions[dimIndex]){
			
				case PhgBinEn_TD:
					dimensions[index].whichDimension = PhgBinEn_TD;
					dimensions[index].virtualIndexPtr = &clcAtnDistanceIndex;
					dimensions[index].dimValuePtr = &clcAtnCurDistance;
					
					/* Special case if there is only 1 bin for this dimension */
					if (PhgBinParams[0].numTDBins == 1) {
						dimensions[index].numBins = 1;
						clcAtnDistIncr = PhgBinParams[0].tdRange/2;
						clcAtnCurDistance = PhgBinParams[0].minTD;
					}
					else {
					
						/* Set the dimensions equal number of bins * sub-binnng */
						dimensions[index].numBins = PhgBinParams[0].numTDBins * clcAtnNumSubBins;
						clcAtnNumSamples *= clcAtnNumSubBins;
						
						/* Compute increment for each projection */
						clcAtnDistIncr = PhgBinParams[0].tdRange/dimensions[index].numBins;
						
						/* Initialize the current distance to the center of the first bin -1,
							then in the general loop it will be incremented one bin size. This
							will start it out at the center of the first bin and have it increment
							to the center of each following bin
						 */
						clcAtnCurDistance = PhgBinParams[0].minTD - (clcAtnDistIncr * .5);
					}

					index++;
					break;
			
				case PhgBinEn_AA:
					dimensions[index].whichDimension = PhgBinEn_AA;
					dimensions[index].virtualIndexPtr = &clcAtnAngleIndex;
					dimensions[index].dimValuePtr = &clcAtnCurAngle;

					/* Set the dimensions equal number of bins */
					dimensions[index].numBins = PhgBinParams[0].numAABins;
				
					/* Compute increment for each projection */
					clcAtnAngleIncr = clcAtnAngleRange/dimensions[index].numBins;
					
					/* Initialize the current distance to the center of the first bin -1,
						then in the general loop it will be incremented one bin size. This
						will start it out at the center of the first bin and have it increment
						to the center of each following bin
					 */
					clcAtnCurAngle = PhgBinParams[0].minAA - (clcAtnAngleIncr * .5);

					index++;
					break;
				
				
				case PhgBinEn_Z1:
					dimensions[index].whichDimension = PhgBinEn_Z1;
					dimensions[index].virtualIndexPtr = &clcAtnZ1Index;
					dimensions[index].dimValuePtr = &clcAtnCurZ1;

					/* Special case if there is only 1 bin for this dimension */
					if (PhgBinParams[0].numZBins == 1) {
						dimensions[index].numBins = 1;
						clcAtnZIncr = PhgBinParams[0].zRange/2;
						clcAtnCurZ1 = PhgBinParams[0].minZ;
					}
					else {
					
						/* Set the dimensions equal number of bins * sub-binnng */
						dimensions[index].numBins = PhgBinParams[0].numZBins * clcAtnNumSubBins;
						clcAtnNumSamples *= clcAtnNumSubBins;
					
						/* Compute increment for each projection */
						clcAtnZIncr = PhgBinParams[0].zRange/dimensions[index].numBins;
						
						/* Initialize the current distance to the center of the first bin -1,
							then in the general loop it will be incremented one bin size. This
							will start it out at the center of the first bin and have it increment
							to the center of each following bin
						 */
						clcAtnCurZ1 = PhgBinParams[0].minZ - (clcAtnZIncr * .5);
					}

					index++;
					break;

				case PhgBinEn_Z2:
					dimensions[index].whichDimension = PhgBinEn_Z2;
					dimensions[index].virtualIndexPtr = &clcAtnZ2Index;
					dimensions[index].dimValuePtr = &clcAtnCurZ2;

					/* Special case if there is only 1 bin for this dimension */
					if (PhgBinParams[0].numZBins == 1) {
						dimensions[index].numBins = 1;
						clcAtnZIncr = PhgBinParams[0].zRange/2;
						clcAtnCurZ2 = PhgBinParams[0].minZ;
					}
					else {
					
						/* Set the dimensions equal number of bins * sub-binnng */
						dimensions[index].numBins = PhgBinParams[0].numZBins * clcAtnNumSubBins;
						clcAtnNumSamples *= clcAtnNumSubBins;
					
						/* Compute increment for each projection */
						clcAtnZIncr = PhgBinParams[0].zRange/dimensions[index].numBins;
						
						/* Initialize the current distance to the center of the first bin -1,
							then in the general loop it will be incremented one bin size. This
							will start it out at the center of the first bin and have it increment
							to the center of each following bin
						 */
						clcAtnCurZ2 = PhgBinParams[0].minZ - (clcAtnZIncr * .5);
					}

					index++;
					break;
					
				case PhgBinEn_THETA:
					break;
				case PhgBinEn_PHI:
					break;
				case PhgBinEn_XR:
					break;
				case PhgBinEn_YR:
					break;
				case PhgBinEn_Null:
					break;
				case PhgBinEn_Energy1:			
					break;
				case PhgBinEn_Energy2:
					break;
				case PhgBinEn_Scatter1:
					break;
				case PhgBinEn_Scatter2:
					break;
					
				default:
					ErStGeneric("Invalid binning dimension in array (ClcAtnCalcAttenuation)");
					goto FAIL;
			}
			
		}
		
		/* Initialize loop variables */
		clcAtnCurProjectionNumber = 1;
		
		/* Loop through slowest varying dimension */
		for (dimensions[0].dimIndex = 0; dimensions[0].dimIndex < dimensions[0].numBins; dimensions[0].dimIndex++){
			
			clcAtnUpdateDimension(dimensions[0].whichDimension, dimensions[0].dimIndex);
			
			for (dimensions[1].dimIndex = 0; dimensions[1].dimIndex < dimensions[1].numBins; dimensions[1].dimIndex++){
			
				clcAtnUpdateDimension(dimensions[1].whichDimension, dimensions[1].dimIndex);
			
				
				for (dimensions[2].dimIndex = 0; dimensions[2].dimIndex < dimensions[2].numBins; dimensions[2].dimIndex++){
				
					clcAtnUpdateDimension(dimensions[2].whichDimension, dimensions[2].dimIndex);
				
					
					for (dimensions[3].dimIndex = 0; dimensions[3].dimIndex < dimensions[3].numBins; dimensions[3].dimIndex++){
					
						clcAtnUpdateDimension(dimensions[3].whichDimension, dimensions[3].dimIndex);

						clcAtnCalcProjection(false);
						
						clcAtnCurProjectionNumber++;
					}
					*(dimensions[3].virtualIndexPtr) = -1;
					clcAtnResetMinimums(dimensions[3].whichDimension);
				}
				*(dimensions[2].virtualIndexPtr) = -1;
					clcAtnResetMinimums(dimensions[2].whichDimension);
			}
			*(dimensions[1].virtualIndexPtr) = -1;
			clcAtnResetMinimums(dimensions[1].whichDimension);
		}
		
		okay = true;
		CANCEL:;
		FAIL:;
	} while (false);
	
	return (okay);
}

/**********************
*	clcAtnCalcProjection
*
*	Purpose:	Calculate the attenuation of the current projection.
*
*	Arguments:
*		Boolean	doAverage - Should we average the bin on this call?
*	Result:	None.
***********************/
void clcAtnCalcProjection(Boolean doAverage)
{
	double					attenuation;			/* The computed attenuation */
	double					k;						/* Temp for computing angles */
	double					x0, y0, s, d1, d2, x1, x2, y1, y2;
	double					tempDist;				/* Temporary distance */
	double					cosAlpha, sineAlpha;
	LbFourByte				imageIndex;				/* Computed index for image */
	PHG_TrackingPhoton		trackingPhoton;			/* The tracking photon */
	PHG_Position			tempPos;				/* Temp position for testing */
	
	if (doAverage) {};		/* Avoid unused parameter compiler warning */
	
	do { /* Process Loop */
		
		/* Compute cosine/sine of alpha because it is used a lot */
		/* NOTE THAT THE NAMES HERE ARE WRONG, BUT THE ORIENTATION IS CORRECT */
		cosAlpha = PHGMATH_Sine(clcAtnCurAngle);
		sineAlpha = -PHGMATH_Cosine(clcAtnCurAngle);

		/* Compute 2d alpha/beta cosines */
		trackingPhoton.angle.cosine_x = -cosAlpha;
		trackingPhoton.angle.cosine_y = -sineAlpha;

		/* Compute some intermediates */
		x0 = (clcAtnTargRadius * cosAlpha)
			- (clcAtnCurDistance * sineAlpha);

		y0 = (clcAtnTargRadius * sineAlpha)
			+ (clcAtnCurDistance * cosAlpha);
	
		/* Find 2d intersection with the target radius */
		if (CylPosFind2dIntersection(x0, y0, trackingPhoton.angle.cosine_x,
				trackingPhoton.angle.cosine_y, clcAtnTargRadius,
				&d1, &d2) == false) {
				
			break;
		}
		
		/* Do "Debug" type check */
	#ifdef DO_YOU_REALLY_WANT_THIS
		if (d2 > d1) {
			sprintf(clcAtnErrStr, "From CylPosFind2dIntersection, d2 > d1\n"
				"d1 = %2.3f and d2 = %3.2f (clcAtnCalcProjection)");
			ErAlert(clcAtnErrStr, false);
		}
	#endif
	
		/* Compute intermediates for z cosine */
		x2 = x0 - d1*cosAlpha;
		y2 = y0 - d1*sineAlpha;
		x1 = x0 - d2*cosAlpha;
		y1 = y0 - d2*sineAlpha;
		
		/* Compute the direction vector, perpendicular to d/a projection */
		trackingPhoton.angle.cosine_z = (clcAtnCurZ1 - clcAtnCurZ2)
			/PHGMATH_SquareRoot(PHGMATH_Square(x1-x2) + PHGMATH_Square(y1-y2)
			+PHGMATH_Square(clcAtnCurZ1 - clcAtnCurZ2));

		/* Compute k */
		k = PHGMATH_SquareRoot(1-PHGMATH_Square(trackingPhoton.angle.cosine_z));
		trackingPhoton.angle.cosine_x = -(k*cosAlpha);
		trackingPhoton.angle.cosine_y = -(k*sineAlpha);
		
		/* Compute our position and direction for traveling onto the target cylinder */
		trackingPhoton.location.x_position = x2;
		trackingPhoton.location.y_position = y2;
		trackingPhoton.location.z_position = clcAtnCurZ2;
	
		/* See if we project onto the object cylinder */
		if (clcAtnProjectToCylinder(&trackingPhoton.location,
				&trackingPhoton.angle, clcAtnObjRadius, clcAtnObjZMin,
				clcAtnObjZMax, clcAtnObjCenterX,clcAtnObjCenterY, &s)
				!= true) {
				
			/* Go to next angle since projection failed */
			break;
		}

		/* See if the photon reaches the opposite surface of the target cylinder within
			the z boundaries
		*/
		tempPos = trackingPhoton.location;
		if (CylPosProjectToTargetCylinder(&tempPos, &trackingPhoton.angle, &tempDist) == false)
			break;
			
		/* Get the indexes */
		SubObjGtPositionIndexes(&trackingPhoton.location,
			&trackingPhoton.sliceIndex,
			&trackingPhoton.xIndex,
			&trackingPhoton.yIndex);
			
		/* Initialize other fields */
		trackingPhoton.flags = 0;
		trackingPhoton.angleIndex = 0;
		trackingPhoton.origSliceIndex = trackingPhoton.sliceIndex;
		trackingPhoton.num_of_scatters = 0;
		trackingPhoton.scatters_in_col = 0;
		trackingPhoton.photon_scatter_weight = 0.0;
		trackingPhoton.photon_primary_weight = 0.0;
		trackingPhoton.photon_current_weight = 0.0;
		trackingPhoton.scatter_target_weight = 0.0;
		trackingPhoton.decay_weight = 0.0;
		trackingPhoton.energy = PhgRunTimeParams.PhgNuclide.photonEnergy_KEV;
		trackingPhoton.number = clcAtnCurProjectionNumber;
		
		/* Compute the attenuation */
		PhoTrkCalcAttenuation(&trackingPhoton, &attenuation);
						
		/* Now compute the actual image index NOTE THAT AA and TD's are NOT reversed */
		imageIndex = ((clcAtnDistanceIndex) * PhgBinParams[0].tdCIsize) + 
			((clcAtnAngleIndex) * PhgBinParams[0].aaCIsize) +
			(clcAtnZ1Index * PhgBinParams[0].z1CIsize) +
			(clcAtnZ2Index * PhgBinParams[0].z2CIsize);
				
		/* Update the image */
		switch(PhgBinParams[0].weight_image_type){
				
			case PHG_BIN_WEIGHT_TYPE_R4:
					
				/* Update the image */
				((float *)clcAtnDAArray)[imageIndex] += exp(attenuation)/(float)clcAtnNumSamples;
				
				break;

			case PHG_BIN_WEIGHT_TYPE_R8:
			
				/* Update the image */
				((double *)clcAtnDAArray)[imageIndex] += exp(attenuation)/(double)clcAtnNumSamples;
				
				break;
		}
		
		
	} while (false);
}

/**********************
*	clcAtnUpdateDimension
*
*	Purpose:	.Update a dimension value.
*
*	Arguments:
*		PHG_BinEnDimensionsTy	whichDimension	- Which dimension to update
*		LbUsFourByte			binIndex		- Which index are we on
*
*	Result:	None.
***********************/
void clcAtnUpdateDimension(PHG_BinEnDimensionsTy whichDimension,
		LbFourByte binIndex)
{
		
	/* Decide which dimension should be updated here */
	switch (whichDimension){
		case PhgBinEn_TD:
			
			/* Compute current distance */
			clcAtnCurDistance += clcAtnDistIncr;
			
			/* See if we need to increment the Distance Index */
			if ((binIndex != 0) && ((binIndex) % clcAtnNumSubBins) == 0)
				clcAtnDistanceIndex++;
				
			break;
	
		case PhgBinEn_AA:
		
			/* Compute current angle */
			clcAtnCurAngle += clcAtnAngleIncr;
			
			/* See if we need to increment the Distance Index */
		/*	if ((binIndex != 0) && ((binIndex) % clcAtnNumSubBins) == 0) */
				clcAtnAngleIndex++;
				
			break;
		
		case PhgBinEn_Z1:
			clcAtnCurZ1 += clcAtnZIncr;
			
			/* See if we need to increment the Distance Index */
	/*		if ((binIndex != 0) && ((binIndex) % clcAtnNumSubBins) == 0) */
				clcAtnZ1Index++;

			break;
			
		case PhgBinEn_Z2:
			clcAtnCurZ2 += clcAtnZIncr;
			
			/* See if we need to increment the Distance Index */
	/*		if ((binIndex != 0) && ((binIndex) % clcAtnNumSubBins) == 0) */
				clcAtnZ2Index++;

			break;
		
		/* Ignore other cases */
		case PhgBinEn_Null:
		case PhgBinEn_Energy1:
		case PhgBinEn_Energy2:
		case PhgBinEn_Scatter1:
		case PhgBinEn_Scatter2:
		case PhgBinEn_TOF:
		case PhgBinEn_THETA:
		case PhgBinEn_PHI:
		case PhgBinEn_XR:
		case PhgBinEn_YR:
		case PhgBinEn_Crystal1:
		case PhgBinEn_Crystal2:
			
			break;
	}
}

/**********************
*	clcAtnProjectToCylinder
*
*	Purpose:	Project photon to cylinder boundary
*
*	Arguments:
*		PHG_Position	*posPtr 	- The starting position
*		PHG_Direction	*dirPtr		- The direction vector
*		double			radius		- The object cylinder radius
*		double			zMin		- The z minimum
*		double			zMax		- The z maximum
*		double			centerX		- The center of the object cylinder
*		double			centerY		- The center of the object cylinder
*		double			*distPtr	- The distance to travel
*
*	Result:	True if the photon intersects the cylinder.
***********************/
Boolean clcAtnProjectToCylinder(PHG_Position *posPtr, PHG_Direction *dirPtr,
			double radius, double zMin, double zMax, double centerX, double centerY,
			double *distPtr)
{
	Boolean			hasSolution = false;	/* Assume no solution */
	double			temp1, temp3, sqrBSqr_4ac, sol1, sol2;
	double			sigma1, sigma2;
	double			oneMinusCosSquGamma;
	double			xCord, yCord;			/* Offsets for cylinder center */
	PHG_Position	newPos;					/* The new position */
	
	
	do { /* Solution loop */
	
		if ((dirPtr->cosine_x > -0.0000001) && (dirPtr->cosine_x < .0000001)) {
			if (dirPtr->cosine_x < 0.0)
				dirPtr->cosine_x = -0.000001;
			else 
				dirPtr->cosine_x = 0.000001;
		}

		if ((dirPtr->cosine_y > -0.0000001) && (dirPtr->cosine_y < .0000001)) {			
			if (dirPtr->cosine_y < 0.0)
				dirPtr->cosine_y = -0.000001;
			else 
				dirPtr->cosine_y = 0.000001;
		}
		
		if ((dirPtr->cosine_z > -0.0000001) && (dirPtr->cosine_z < .0000001)) {
			if (dirPtr->cosine_z < 0.0)
				dirPtr->cosine_z = -0.000001;
			else 
				dirPtr->cosine_z = 0.000001;
		}

		/* Fudge the object radius a tiny amount to make sure we are
			in the cylinder */
		radius -= .001;
		
		/* Normalize position for center of object cylinder */
		xCord = posPtr->x_position - centerX;
		yCord = posPtr->y_position - centerY;

		/* Calculate first temp */
		temp1 = -((xCord * dirPtr->cosine_x) +
			(yCord * dirPtr->cosine_y));
			
		/* Calculate second temp */
		oneMinusCosSquGamma = 1 - PHGMATH_Square(dirPtr->cosine_z);

		
		/* Calculate value under the square root in quadratic formula */
		temp3 = (oneMinusCosSquGamma * PHGMATH_Square(radius)) -
			PHGMATH_Square((dirPtr->cosine_x * yCord) -
			(dirPtr->cosine_y * xCord));
		
		/* If it is negative, there is no solution so bolt */
		if (temp3 < 0.0)
			break;
			
		
		/* Compute value under the radical */
		sqrBSqr_4ac = PHGMATH_SquareRoot(temp3);
		
		/* Compute two possible solutions */
		sol1 = (temp1 + sqrBSqr_4ac)/oneMinusCosSquGamma;
		sol2 = (temp1 - sqrBSqr_4ac)/oneMinusCosSquGamma;
		
		/* Project both solutions in z direction */
		sigma1 = clcAtnCurZ2 + (sol1 * dirPtr->cosine_z);
		sigma2 = clcAtnCurZ2 + (sol2 * dirPtr->cosine_z);
		
		/* See if we don't intersect at all */
		if (((sigma1 >= zMax) && (sigma2 >= zMax)) ||
				((sigma2 <= zMin) && (sigma2 <= zMin)))
			break;
		else if (sigma1 >= zMax) {
			*distPtr = (zMax - clcAtnCurZ2)/dirPtr->cosine_z;	
		}
		else if (sigma1 <= zMin) {
			*distPtr = (zMin - clcAtnCurZ2)/dirPtr->cosine_z;
		}
		else if (sol1 < sol2)
			*distPtr = sol1;
		else
			*distPtr = sol2;
				
		/* Project to cyliner */
		newPos.x_position = posPtr->x_position + (*distPtr * dirPtr->cosine_x);
		newPos.y_position = posPtr->y_position + (*distPtr * dirPtr->cosine_y);
		newPos.z_position = posPtr->z_position + (*distPtr * dirPtr->cosine_z);
		
		/* Check for boundaries (this shouldn't happen due to code above) */
		if ((newPos.z_position < zMin) || (newPos.z_position > zMax))
			break;
			
		/* We made it, so update the position */
		*posPtr = newPos;
			
		hasSolution = true;
	} while (false);
	
	return (hasSolution);
}

/**********************
*	ClcAtnInitialize
*
*	Purpose:	Initialize this module.
*
*	Result:	Always returns zero indicating no error.
***********************/
Boolean ClcAtnInitialize(int argc, char *argv[])
{
	char		fileName[1024];				/* Storage for file name */
	LbFourByte	randSeed;					/* Seed for random generator */
	
	do { /* Process Loop */
		
		/* If there are no command line arguments, prompt user for param
			file. Only one file is allowed.
		   Else, they may specify multiple param files which will all be
		   processed according to the binning parameters of the first file
		 */
		if (argc == 1) {
			/* Ask for the file name */
			LbInAsk("Enter name of param file", 0, false,
					&canceled, 0, 0, 0, 0,
					fileName);
		
			/* Bolt if we canceled */
			if (canceled) {
				ErStCancel("User canceled.");
				goto CANCEL;
			}
			strcpy(PhgRunTimeParams.PhgParamFilePath, fileName);
		}
		else {
		
			/* Get first param file and save number of param files to process */
			strcpy(PhgRunTimeParams.PhgParamFilePath,argv[1]);
		}
		
		/* Get our parameters */
		if (!PhgGetRunTimeParams())
			break;
				
		/* Initialize the math library */
		randSeed = PhgRunTimeParams.PhgRandomSeed;
		if (!PhgMathInit(&randSeed))
			break;
		
		/* Initialize the emission list manager */
		if (!EmisListInitialize())
			break;
				
		/* Initialize the productivity module */
		if (!ProdTblInitialize()) {
			break;
		}
		
		/* Initialize the sub-object manager */
		if (!SubObjInitialize()) {
			break;
		}
		
		/* Setup the Subobject modules */
		if (!SubObjCreate()) {
			break;
		}
			/* Initialize the Cylinder Positions */
		{
			/* Set object cylinder */
			if (!CylPosInitObjectCylinder()) {
				break;
			}
			
			/* Set the criticial zone */
			if (!CylPosInitCriticalZone(PhgRunTimeParams.Phg_AcceptanceAngle)) {
				break;
			}
			
			/* Set the limit cylinder */
			CylPosInitLimitCylinder();
		}
		
	
		/* Initialize the photon tracking module */
		if (!PhoTrkInitialize()) {
			break;
		}
	
		/* Now we should be binning on the fly here, or there is an error */
		if (PHG_IsBinOnTheFly() == false) {
			ErStGeneric("You must supply a binning parameter file for this utility.");
			break;
		}
		
		/* Get the binning parameters */
		if (PhgBinInitParams(PhgRunTimeParams.PhgBinParamsFilePath[0], &PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0]) == false) {
			break;
		}

		clcAtnIsInitialized = true;
		CANCEL:;
	} while (false);
		
	return (clcAtnIsInitialized);
}

#undef PHG_CALC_ATTEN_MAIN
