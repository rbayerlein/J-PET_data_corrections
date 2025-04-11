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
*			Module Name:		atten.correct.c
*			Revision Number:	1.5
*			Date last revised:	4 June 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	17 May, 1995
*
*			Module Overview:	Performs an attenuation correction by multiplying.
*								a PHG image by it's computed attenuation map.
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

/* LOCAL TYPES */
typedef	struct {	
	PHG_BinEnDimensionsTy	whichDimension;	/* Which dimension is represented */
	LbFourByte				numBins;		/* Number of bins for dimension */
	LbFourByte				dimIndex;		/* Current index for this dimension */
} atnCorDimensionInfoTy;

/* LOCAL GLOBALS */
static Boolean	canceled;								/* Global cancelation flag */

/* The following local globals are all used due to optimization issues */
/*static	double			atnCorCurAngle;*/				/* Angle for current projection */
/*static	double			atnCorCurDistance;*/			/* Distance for current projection */
/*static	double			atnCorCurZ1;*/				/* Z1 for current projection */
/*static	double			atnCorCurZ2;*/				/* Z2 for current projection */
/*static	double			atnCorTargRadius;*/			/* The target cylinder radius */
/*static	double			atnCorObjRadius;*/			/* The object cylinder radius */
/*static	double			atnCorObjZMin;*/				/* Minimum Z of object cylinder */
/*static	double			atnCorObjZMax;*/				/* Maximum Z of object cylinder */
/*static	double			atnCorObjCenterX;*/			/* Center X of object cylinder */
/*static	double			atnCorObjCenterY;*/			/* Center Y of object cylinder */
/*static	double			atnCorAngleRange;*/			/* Range of angles, varies due to PET/SPECT */
/*static	LbFourByte		atnCorCurProjectionNumber;*/	/* Current projection number */
static	LbFourByte		*atnCorAngleIndex;			/* The current angle index */
static	LbFourByte		*atnCorDistanceIndex;		/* The current distance index */
static	LbFourByte		*atnCorZ1Index;				/* The current Z1 index */
static	LbFourByte		*atnCorZ2Index;				/* The current Z2 index */
static	LbFourByte		*atnCorS1Index;				/* The current S1 index */
static	LbFourByte		*atnCorS2Index;				/* The current S2 index */
static	LbFourByte		*atnCorE1Index;				/* The current E1 index */
static	LbFourByte		*atnCorE2Index;				/* The current E2 index */

static	void			*atnCorAttenuationAry = 0;	/* The attenuation data */
static	void			*atnCorImgWtAry = 0;		/* The image weight data */
static	void			*atnCorImgWtSqAry = 0;		/* The image weight squared data */
static	FILE			*atnCorImgWtFile = 0;		/* The weight image file */
static	FILE			*atnCorImgWtSqFile = 0;		/* The weight sqared image file */
	
static	PhoHFileHdrTy	atnWeightImgHdr;			/* Header for weight image file */
static	PhoHFileHdrTy	atnWeightSquImgHdr;			/* Header for weight squared image file */
static	LbHdrHkTy		atnWeightImgHdrHk;			/* Hook for weight image file header */
static	LbHdrHkTy		atnWeightSquImgHdrHk;		/* Hook for count image file header */

/* PROTOTYPES */
Boolean			AttnCorrect(int argc, char *argv[]);
Boolean			atnCorInitialize(int argc, char *argv[]);
Boolean			atnCorGtImgData(void);

/* FUNCTIONS */

/**********************
*	AttnCorrect
*
*	Purpose:	Performs attenuation correction.
*
*	Result:	Always returns zero indicating no error.
***********************/
Boolean AttnCorrect(int argc, char *argv[])
{
	Boolean					okay = false;				/* Process Loop */
	LbFourByte				dimIndex;					/* Index for traversing dimensions */
	LbFourByte				attDataIndex;				/* Index into attenuation data */
	LbFourByte				imgDataIndex;				/* Index into image data */
	atnCorDimensionInfoTy	imgDimensions[PHGBIN_NUM_DIMENSIONS];



	do { /* Process Loop */
		
		/* Initialize this library */
		if (ClcAtnInitialize(argc, argv) == false)
			break;
		
		/* Read in the image data */
		if (atnCorGtImgData() == false){
			break;
		}
		
		/* Compute the attenuation map */
		if (ClcAtnCalcAttenuation(&atnCorAttenuationAry) != true){
			break;
		}
		
		/* Setup dimension index for traversing all dimensions */
		for (dimIndex = PHGBIN_NUM_DIMENSIONS-1; dimIndex >= 0; dimIndex--){
		
			switch (PhgBinParams[0].PhgBinDimensions[dimIndex]){
			
				case PhgBinEn_TD:
					atnCorDistanceIndex = &imgDimensions[dimIndex].dimIndex;
					imgDimensions[dimIndex].numBins = PhgBinParams[0].numTDBins;
					imgDimensions[dimIndex].whichDimension = PhgBinEn_TD;
			
					break;
			
				case PhgBinEn_AA:
					atnCorAngleIndex = &imgDimensions[dimIndex].dimIndex;
					imgDimensions[dimIndex].numBins = PhgBinParams[0].numAABins;
					imgDimensions[dimIndex].whichDimension = PhgBinEn_AA;

					break;
				
				
				case PhgBinEn_Z1:
					atnCorZ1Index = &imgDimensions[dimIndex].dimIndex;
					imgDimensions[dimIndex].numBins = PhgBinParams[0].numZBins;
					imgDimensions[dimIndex].whichDimension = PhgBinEn_Z1;

					break;

				case PhgBinEn_Z2:
					atnCorZ2Index = &imgDimensions[dimIndex].dimIndex;
					imgDimensions[dimIndex].numBins = PhgBinParams[0].numZBins;
					imgDimensions[dimIndex].whichDimension = PhgBinEn_Z2;

					break;
					
	
				case PhgBinEn_Energy1:			

					atnCorE1Index = &imgDimensions[dimIndex].dimIndex;
					imgDimensions[dimIndex].numBins = PhgBinParams[0].numE1Bins;
					imgDimensions[dimIndex].whichDimension = PhgBinEn_Energy1;
					break;
					
				case PhgBinEn_Energy2:

					atnCorE2Index = &imgDimensions[dimIndex].dimIndex;
					imgDimensions[dimIndex].numBins = PhgBinParams[0].numE1Bins;
					imgDimensions[dimIndex].whichDimension = PhgBinEn_Energy2;
					break;
					
				case PhgBinEn_Scatter1:

					atnCorS1Index = &imgDimensions[dimIndex].dimIndex;
					imgDimensions[dimIndex].numBins = PhgBinParams[0].numS1Bins;
					imgDimensions[dimIndex].whichDimension = PhgBinEn_Scatter1;
					break;
					
				case PhgBinEn_Scatter2:

					atnCorS2Index = &imgDimensions[dimIndex].dimIndex;
					imgDimensions[dimIndex].numBins = PhgBinParams[0].numS2Bins;
					imgDimensions[dimIndex].whichDimension = PhgBinEn_Scatter2;
					break;
					
				default:
					ErStGeneric("Invalid binning dimension in image array (CalcAttenuation)");
					goto FAIL;
			}
			
		}


		/* Loop image data and attenuation data starting with slowest varying dimension */
		for (imgDimensions[0].dimIndex = 0; imgDimensions[0].dimIndex < imgDimensions[0].numBins; imgDimensions[0].dimIndex++){
			
		for (imgDimensions[1].dimIndex = 0; imgDimensions[1].dimIndex < imgDimensions[1].numBins; imgDimensions[1].dimIndex++){
			
		for (imgDimensions[2].dimIndex = 0; imgDimensions[2].dimIndex < imgDimensions[2].numBins; imgDimensions[2].dimIndex++){
				
		for (imgDimensions[3].dimIndex = 0; imgDimensions[3].dimIndex < imgDimensions[3].numBins; imgDimensions[3].dimIndex++){
				
		for (imgDimensions[4].dimIndex = 0; imgDimensions[4].dimIndex < imgDimensions[4].numBins; imgDimensions[4].dimIndex++){
			
		for (imgDimensions[5].dimIndex = 0; imgDimensions[5].dimIndex < imgDimensions[5].numBins; imgDimensions[5].dimIndex++){
			
		for (imgDimensions[6].dimIndex = 0; imgDimensions[6].dimIndex < imgDimensions[6].numBins; imgDimensions[6].dimIndex++){
			
		for (imgDimensions[7].dimIndex = 0; imgDimensions[7].dimIndex < imgDimensions[7].numBins; imgDimensions[7].dimIndex++){


			/* Compute the index into the attenuation data */
			attDataIndex = (*atnCorDistanceIndex * PhgBinParams[0].tdCIsize) + 
							(*atnCorAngleIndex * PhgBinParams[0].aaCIsize) +
							(*atnCorZ1Index * PhgBinParams[0].z1CIsize) +
							(*atnCorZ2Index * PhgBinParams[0].z2CIsize);
		
			/* Compute the index into the image data */
			imgDataIndex = (*atnCorDistanceIndex * PhgBinParams[0].tdCIsize) + 
							(*atnCorAngleIndex * PhgBinParams[0].aaCIsize) +
							(*atnCorE1Index * PhgBinParams[0].energy1CIsize) +
							(*atnCorE2Index * PhgBinParams[0].energy2CIsize) +
							(*atnCorS1Index * PhgBinParams[0].scatter1CIsize) +
							(*atnCorS2Index * PhgBinParams[0].scatter2CIsize) +
							(*atnCorZ1Index * PhgBinParams[0].z1CIsize) +
							(*atnCorZ2Index * PhgBinParams[0].z2CIsize);
							
			/* Multiply the image by the attenuation data */
				switch(PhgBinParams[0].weight_image_type){
						
					case PHG_BIN_WEIGHT_TYPE_R4:
					
						/* We don't check for floating point overflow because
							the hardware will get it
						*/
						((float *)atnCorImgWtAry)[imgDataIndex] *= 
							((float *)atnCorAttenuationAry)[attDataIndex];

						((float *)atnCorImgWtSqAry)[imgDataIndex] *= 
							PHGMATH_Square(((float *)atnCorAttenuationAry)[attDataIndex]);
						break;

					case PHG_BIN_WEIGHT_TYPE_R8:
					
						/* We don't check for floating point overflow because
							the hardware will get it
						*/
						((double *)atnCorImgWtAry)[imgDataIndex] *= 
							((double *)atnCorAttenuationAry)[attDataIndex];

						((double *)atnCorImgWtSqAry)[imgDataIndex] *= 
							PHGMATH_Square(((double *)atnCorAttenuationAry)[attDataIndex]);
						break;
				}
			
		} } } } } } } }
		
		/* Write out the corrected weights file */
		if (fwrite(atnCorImgWtAry, PhgBinParams[0].weightImageSize, 1, atnCorImgWtFile) 	!= 1){
			ErStFileError("Unable to write image file.");
			break;
		}
		
		/* Set the header flag indicating that this file has been attenuation corrected */
		if (PhgHdrStAttenuationCorrected(&atnWeightImgHdrHk) == false)
			break;
		
		/* Write out the corrected weights squared file */
		if ((fwrite(atnCorImgWtSqAry, PhgBinParams[0].weightSquImageSize, 1, atnCorImgWtSqFile)) != 1) {
			ErStFileError("Unable to write header to corrected weights squared file.");
			break;
		}
		
		/* Set the header flag indicating that this file has been attenuation corrected */
		if (PhgHdrStAttenuationCorrected(&atnWeightSquImgHdrHk) == false)
			break;
		
		okay = true;
		CANCEL:;
		FAIL:;
	} while (false);
	

		
	/* Handle error situation if one exists */
	if (!okay && canceled) {
		ErHandle("User canceled CalcAttenuation.", false);
		okay = true;
	}

	/* Free memory */
	if (atnCorAttenuationAry != 0)
		LbMmFree(&atnCorAttenuationAry);
	

	/* Free  memory */
	if (atnCorImgWtAry != 0)
		LbMmFree(&atnCorImgWtAry);

	if (atnCorImgWtSqAry != 0)
		LbMmFree(&atnCorImgWtSqAry);

	/* Close local file */
	if (atnCorImgWtFile != 0)
		fclose(atnCorImgWtFile);

	/* Close local file */
	if (atnCorImgWtSqFile != 0)
		fclose(atnCorImgWtSqFile);
			
	/* Terminate  modules */
	PhoTrkTerminate();
	SubObjTerminate();
	PhgMathTerminate();
	EmisListTerminate();
	

	/* Quit the program */
	return (okay);
}


/**********************
*	atnCorInitialize
*
*	Purpose:	Initialize this module.
*
*	Result:	True unless an error occurs.
***********************/
Boolean atnCorInitialize(int argc, char *argv[])
{
	Boolean		okay = false;				/* Process Loop */
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

		okay = true;
		CANCEL:;
	} while (false);
		
	return (okay);
}

/**********************
*	atnCorGtImgData
*
*	Purpose:	Read in the image data.
*
*	Result:	True unless an error occurs.
***********************/
Boolean atnCorGtImgData()
{
	Boolean			okay = false;			/* Process Flag */
	char			imageFileName[1024];	/* Name of output file */
	FILE			*dataFile;
	
	do { /* Process Loop */
	
		/* Ask for the file name for the adjusted weight image */
		LbInAskNew(LBINFg_Unique, "Enter name of corrected weights file", 1, false,
				&canceled, PhgBinParams[0].weightImgFilePath, 0, 0, 0,
				imageFileName);
	
		/* Bolt if we canceled */
		if (canceled) {
			ErStCancel("User canceled.");
			goto CANCEL;
		}

		/* Open the input file */
		if ((dataFile = LbFlFileOpen(PhgBinParams[0].weightImgFilePath, "r+b")) == 0) {
			ErStFileError("Unable to open weights file.");
			break;
		}

		/* Read in the header */
		if (PhgHdrGtParams(dataFile, &atnWeightImgHdr, &atnWeightImgHdrHk) == false){
			break;
		}

		/* Allocate the weight image data array */
		if ((atnCorImgWtAry = LbMmAlloc(PhgBinParams[0].weightImageSize)) == 0) {
			break;
		}
		
		/* Read  the weight image data */
		if ((fread(atnCorImgWtAry, PhgBinParams[0].weightImageSize, 1, dataFile)) != 1) {
			ErStFileError("Unable to read weight data.");
			break;
		}

		/* Close the input file */
		fclose(dataFile);
		
		/* Open the output file (this comes after closing the input file, because it
			may be the same.
		 */
		if ((atnCorImgWtFile = LbFlFileOpen(imageFileName, "w+b")) == 0) {
			ErStFileError("Unable to create corrected weights file.");
			break;
		}

		/* Set the header hook to the new file */
		if (PhgHdrStFile(&atnWeightImgHdrHk, atnCorImgWtFile) == false) {
			break;
		}
				
		/* Ask for the file name for the adjusted weight squared image */
		LbInAskNew(LBINFg_Unique, "Enter name of corrected weights squared file", 1, false,
				&canceled, PhgBinParams[0].weightSquImgFilePath, 0, 0, 0,
				imageFileName);
	
		/* Bolt if we canceled */
		if (canceled) {
			ErStCancel("User canceled.");
			goto CANCEL;
		}

		/* Open the input file */
		if ((dataFile = LbFlFileOpen(PhgBinParams[0].weightSquImgFilePath, "r+b")) == 0) {
			ErStFileError("Unable to open weights squared file.");
			break;
		}
		
		/* Allocate the weight image data array */
		if ((atnCorImgWtSqAry = LbMmAlloc(PhgBinParams[0].weightSquImageSize)) == 0) {
			break;
		}

		/* Read in the header */
		if (PhgHdrGtParams(dataFile, &atnWeightSquImgHdr, &atnWeightSquImgHdrHk) == false){
			break;
		}
		
		/* Read  the weight image data */
		if ((fread(atnCorImgWtSqAry, PhgBinParams[0].weightSquImageSize, 1, dataFile)) != 1) {
			ErStFileError("Unable to read weight squared data.");
			break;
		}

		/* Close the input file */
		fclose(dataFile);
	
		/* Open the output file */
		if ((atnCorImgWtSqFile = LbFlFileOpen(imageFileName, "w+b")) == 0) {
			ErStFileError("Unable to create corrected weights file.");
			break;
		}

		/* Set the header hook to the new file */
		if (PhgHdrStFile(&atnWeightSquImgHdrHk, atnCorImgWtSqFile) == false) {
			break;
		}
		
	okay = true;
	CANCEL:;
	} while (false);
	
	
	/* Do some error handling */
	if (okay == false) {
		/* Free memory if allocated */
		if (atnCorImgWtAry != 0) {
			LbMmFree(&atnCorImgWtAry);
		}
		if (atnCorImgWtSqAry != 0) {
			LbMmFree(&atnCorImgWtSqAry);
		}

		/* Close files if necessary */
		if (dataFile != 0) {
			fclose(dataFile);
			dataFile = 0;
		}
		if (atnCorImgWtFile != 0) {
			fclose(atnCorImgWtFile);
			atnCorImgWtFile = 0;
		}
		if (atnCorImgWtSqFile != 0) {
			fclose(atnCorImgWtSqFile);
			atnCorImgWtSqFile = 0;
		}
	}
	
	return (okay);
}
#undef PHG_CALC_ATTEN_MAIN
