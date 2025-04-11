/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1997-2013 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			line3d.c
*     Revision Number:		1.6
*     Date last revised:	23 July 2013
*     Programmer:			Steven Vannoy
*     Date Originated:		Thursday, March 27, 1997
*
*     Module Overview:	Computes statistics and extracts an "axial" line in a
*						3D data set.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*     Global variables defined:
*
*	  Global macros defined:
*
*********************************************************************************/

#define LINE_3D


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

/* Local Constants */

/* Local Globals */
static	char			errString[1024];		/* For building error messages */
static	Boolean			canceled = false;		/* Cancelation flag */
static	LbUsFourByte	line3dRunTimeOptions=0;	/* The runtime options specified */

#define Is_Float()			LbFgIsSet(line3dRunTimeOptions, LBFlag0)	/* Did user specify float data? */
#define Is_Double()			LbFgIsSet(line3dRunTimeOptions, LBFlag1)	/* Did user specify double data? */
#define Is_LbFourByte()		LbFgIsSet(line3dRunTimeOptions, LBFlag2)	/* Did user specify LbFourByte data? */
#define Is_PhgParams()		LbFgIsSet(line3dRunTimeOptions, LBFlag3)	/* Did user supply a PHG parameter file? */

/* Prototypes */
Boolean	line3d(int argc, char *argv[]);
Boolean line3DInitializePHG(char *fileName);
Boolean line3DGtImgData(void);
void	line3dTerminateForPHG(void);	
Boolean	line3dGetIndexesForPHG(LbUsFourByte *numSlices, LbUsFourByte *sliceSize, LbUsFourByte *numColumns);	
Boolean	line3dInitForPHG(char *fileName, LbUsFourByte *numSlices, LbUsFourByte *sliceSize, LbUsFourByte *numColumns);
Boolean	line3dInitForDataFile(char *fileName, void **data, LbUsFourByte *numSlices, LbUsFourByte *sliceSize, LbUsFourByte *numColumns);
void	line3dPrintForPHG(LbUsFourByte imageIndex);
void	line3dPrintForDataFile(LbUsFourByte imageIndex, void *data);


/**********************
*	line3d
*
*	Purpose:	Execute the program.
*
*	Result:	None.
***********************/

Boolean	line3d(int argc, char *argv[])
{
	Boolean			okay = false;			/* Process flag */
	
	LbUsFourByte	sliceSize;				/* Size of a slice for weights */
	LbUsFourByte	index;					/* Loop Control */
/*	double			wtSum=0;*/				/* Cumulative counter */
/*	double			wtSqSum=0;*/				/* Cumulative counter */
/*	LbUsFourByte	ctSum=0;*/				/* Cumulative counter */
/*	double			sumSq =0.0;*/				/* Sum Squared */
/*	double			minWt=FLT_MAX;*/			/* Minimum value */
/*	double			maxWt=FLT_MIN;*/			/* Maximum  value */
/*	LbUsFourByte	minCt=LBUSFOURBYTE_MAX;*/	/* Minimum value */
/*	LbUsFourByte	maxCt=0;*/				/* Maximum  value */
	LbUsFourByte	rowNumber;					/* x bin */
	LbUsFourByte	colNumber;					/* y bin */
	LbUsFourByte	imageIndex;				/* The index into the value of interest */
	LbUsFourByte	numColumns;				/* The number of columns */
	LbUsFourByte	numSlices;				/* The number of slices */
	void			*data;					/* The data */
	
	/* The following variables are for getting run time options from
		the command line 
	*/
	#define	NUM_FLAGS	4
	
	char				*knownOptions[] = {"fdip"};
	
	char				optArgs[NUM_FLAGS][LBEnMxArgLen];
	LbUsFourByte		optArgFlags = LBFlag0 + LBFlag1 + LBFlag2 + LBFlag3;
	LbUsFourByte		argIndex;
	
	do	 {
			
		/* Get our runtime options */
		if (!LbEnGetOptions(argc, argv, knownOptions,
				&line3dRunTimeOptions, optArgs, optArgFlags, &argIndex)) {
	
			break;
		}
		
		/* Verify command line contains an argument */
		if (argc < 2) {
			/* Tell them an image file is required */
			ErAbort("\nThis program requires an input-file name (either data or phg parameter file).\n");
		}
		
		/* Verify the command line indicates the type of data to collapse */
		if (line3dRunTimeOptions == 0) {
			ErAbort("\nThis program requires a switch to indicate the type of data being collapsed.\n"
				"Usage: collapse -{pdfi} input.name output.name\n"
				"Where: -p indicates the input file is a phg parameter file\n"
				"       -d = double precision reals -f single precision reals -i = 4 byte integer data\n");
		}

		if (Is_PhgParams()) {
			if (line3dInitForPHG(argv[argIndex], &numSlices, &sliceSize, &numColumns) == false)
				break;
		}
		else {
			if (line3dInitForDataFile(argv[argIndex], &data, &numSlices, &sliceSize, &numColumns) == false)
				break;
		}
		
		/* Get the x index */
		rowNumber = LbInAskFourByte("Enter row number (starting at zero)",
			1, false, false, 1, &canceled, 0, 0, 0, 0) - 1;

		colNumber = LbInAskFourByte("Enter column number (zero is at the bottom)",
			1, false, false, 1, &canceled, 0, 0, 0, 0) - 1;
			
		
		/* Verify user selected values are within range of given data */
		if ((rowNumber > (sliceSize/numColumns)) || (colNumber > numColumns)) {
			sprintf(errString, "Your row and/or column selection is out of range\n"
				"you said row %ld and column %ld but there are only %ld rows and %ld columns",
				(unsigned long)rowNumber, (unsigned long)colNumber, 
				(unsigned long)(sliceSize/numColumns), (unsigned long)numColumns);
			ErStGeneric(errString);
			goto FAIL;
		}
		
		/* Compute index into desired value */
		imageIndex = (rowNumber * numColumns) + colNumber;

		LbInPrintf("\nLine values\n");
			
		/* Loop through and get the data */
		for (index = 0; index < numSlices; index++){
			if (Is_PhgParams()) {
				line3dPrintForPHG(imageIndex);
			}
			else {
				line3dPrintForDataFile(imageIndex, data);
			}
			
			LbInPrintf("\n");
			imageIndex += sliceSize;
		}	

		
	CANCELED:;
	okay = true;
	FAIL:;
	} while (false);

	/* Close the input files */
	if (Is_PhgParams()){
		if (PhgBinFields[0].CountFile != 0) {
			fclose(PhgBinFields[0].CountFile);
			PhgBinFields[0].CountFile = 0;
		}
		if (PhgBinFields[0].WeightFile != 0) {
			fclose(PhgBinFields[0].WeightFile);
			PhgBinFields[0].WeightFile = 0;
		}
		if (PhgBinFields[0].WeightSquFile != 0) {
			fclose(PhgBinFields[0].WeightSquFile);
			PhgBinFields[0].WeightSquFile = 0;
		}
	}
	
	/* Free Memory */
	if (Is_PhgParams()){
		if (PhgBinData[0].countImage != 0)
			LbMmFree((void **) &PhgBinData[0].countImage);
			
		if (PhgBinData[0].weightImage != 0)
			LbMmFree((void **) &PhgBinData[0].weightImage);
			
		if (PhgBinData[0].weightSquImage != 0)
			LbMmFree((void **) &PhgBinData[0].weightSquImage);
	}
	else {
		LbMmFree(&data);
	}
	
	if (Is_PhgParams()){
		line3dTerminateForPHG();
	}
	
	return(okay);
}

/**********************
*	line3DInitializePHG
*
*	Purpose:	Initialize this module.
*
*	Result:	True unless an error occurs.
***********************/
Boolean line3DInitializePHG(char *fileName)
{
	Boolean		okay = false;				/* Process Loop */
	
	do { /* Process Loop */
		
		
		/* Get first param file and save number of param files to process */
		strcpy(PhgRunTimeParams.PhgParamFilePath, fileName);
		

		/* Get our parameters */
		if (!PhgGetRunTimeParams())
			break;
				
		
		/* Initialize the math library */
		if (!PhgMathInit((LbFourByte *)&PhgRunTimeParams.PhgRandomSeed))
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
*	line3DGtImgData
*
*	Purpose:	Read in the image data.
*
*	Result:	True unless an error occurs.
***********************/
Boolean line3DGtImgData()
{
	Boolean			okay = false;			/* Process Flag */
	FILE			*dataFile;
	
	do { /* Process Loop */

		/* NOTE: We override addToExisting so the existing files get opened */
		PhgBinParams[0].addToExistingImg = true;
		
		
		/* See if they are binning counts */
		if (PhgBinParams[0].doCounts == true) {

			/* Verify the image exists already */
			if ((dataFile = LbFlFileOpen(PhgBinParams[0].countImgFilePath, "rb")) == 0) {
					ErStGeneric("The count image doesn't exist, yet the binning parameters has 'do counts = true', this test can't be run");
					goto FAIL;
			}
			else {
				fclose(dataFile);
			}
		
			/* Allocate the image buffer */
			if ((PhgBinData[0].countImage = LbMmAlloc(PhgBinParams[0].countImageSize))
					== 0) {
				break;
			}			
			
			/* Open the count image file */
			if (PhgBinOpenImage(&PhgBinParams[0], &PhgBinFields[0], PhoHFileEn_BIN_CT,
					PhgBinParams[0].countImgFilePath,
					&PhgBinFields[0].CountFile, &PhgBinFields[0].CountImgHdr,
					&PhgBinFields[0].CountImgHdrHk,PhgBinParams[0].countImageSize,
					(void *)PhgBinData[0].countImage) == false) {
				goto FAIL;
			}
		}

		/* See if they are binning weights */
		if (PhgBinParams[0].doWeights == true) {

			/* Verify the image exists already */
			if ((dataFile = LbFlFileOpen(PhgBinParams[0].weightImgFilePath, "rb")) == 0) {
					ErStGeneric("The weight image doesn't exist, yet the binning parameters has 'do weights = true', this test can't be run");
					goto FAIL;
			}
			else {
				fclose(dataFile);
			}
		
			/* Allocate the image buffer */
			if ((PhgBinData[0].weightImage = LbMmAlloc(PhgBinParams[0].weightImageSize))
					== 0) {
				break;
			}
						
			/* Open the weight image file */
			if (PhgBinOpenImage(&PhgBinParams[0], &PhgBinFields[0], PhoHFileEn_BIN_WT,
					PhgBinParams[0].weightImgFilePath,
					&PhgBinFields[0].WeightFile, &PhgBinFields[0].WeightImgHdr,
					&PhgBinFields[0].WeightImgHdrHk,
							(((PhgBinParams[0].sumAccordingToType == false) && (PhgBinParams[0].weight_image_type == PHG_BIN_WEIGHT_TYPE_R4))
						? PhgBinParams[0].weightImageSize/2 : PhgBinParams[0].weightImageSize),
					(void *)PhgBinData[0].weightImage) == false) {
				break;
			}
		}

		/* See if they are binning weights squared */
		if (PhgBinParams[0].doWeightsSquared == true) {

			/* Verify the image exists already */
			if ((dataFile = LbFlFileOpen(PhgBinParams[0].weightSquImgFilePath, "rb")) == 0) {
					ErStGeneric("The weight square image doesn't exist, yet the binning parameters has 'do weights squared = true', this test can't be run");
					goto FAIL;
			}
			else {
				fclose(dataFile);
			}
		
			/* Allocate the image buffer */
			if ((PhgBinData[0].weightSquImage = LbMmAlloc(PhgBinParams[0].weightSquImageSize))
					== 0) {
				break;
			}
			
			/* Open the weight squared image file */
			if (PhgBinOpenImage(&PhgBinParams[0], &PhgBinFields[0], PhoHFileEn_BIN_WTSQ,
					PhgBinParams[0].weightSquImgFilePath,
					&PhgBinFields[0].WeightSquFile, &PhgBinFields[0].WeightSquImgHdr,
					&PhgBinFields[0].WeightSquImgHdrHk,
							(((PhgBinParams[0].sumAccordingToType == false) && (PhgBinParams[0].weight_image_type == PHG_BIN_WEIGHT_TYPE_R4))
						? PhgBinParams[0].weightSquImageSize/2 : PhgBinParams[0].weightSquImageSize),
					(void *)PhgBinData[0].weightSquImage) == false) {
				break;
			}

		}
		
	okay = true;
	CANCEL:;
	FAIL:;
	} while (false);
	
	
	/* Do some error handling */
	if (okay == false) {
		/* Free memory if allocated */
		if (PhgBinData[0].weightSquImage != 0) {
			LbMmFree(&PhgBinData[0].weightSquImage);
		}
		if (PhgBinData[0].weightImage != 0) {
			LbMmFree(&PhgBinData[0].weightImage);
		}
		if (PhgBinData[0].countImage != 0) {
			LbMmFree(&PhgBinData[0].countImage);
		}
	}
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			line3dGetIndexesForPHG
*
*			Summary:	Terminate the various modules.
*			Arguments:
*				LbUsFourByte	*numSlices	- The number of slices
*				LbUsFourByte	*sliceSize	- The size of the slices
*				LbUsFourByte	*numColumns	- The number of columns
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean line3dGetIndexesForPHG(LbUsFourByte *numSlices, LbUsFourByte *sliceSize, LbUsFourByte *numColumns)	
{
	Boolean		okay = false;
	LbFourByte	index;	/* LCV */
	LbFourByte				curDim = 0;
	
	
		/* First, verify there are only three dimensions, and identify where they are */
		for (index = 0; index < PHGBIN_NUM_DIMENSIONS; index++) {
	
			/* Figure out size of slices */		
			switch (PhgBinParams[0].PhgBinDimensions[index]){
			
				case PhgBinEn_TD:

					if (PhgBinParams[0].numTDBins > 1) {

						if (curDim == 0) {
							*numColumns = PhgBinParams[0].numTDBins;
						}
						else if (curDim == 2) {
							*numSlices = PhgBinParams[0].numTDBins;
							*sliceSize = PhgBinParams[0].tdCIsize;
						}
						else if (curDim > 2) { 
							ErStGeneric("This data has more than three dimensions");
						}
						curDim++;
					}
						
					break;
			
				case PhgBinEn_AA:

					if (PhgBinParams[0].numAABins > 1) {

						if (curDim == 0) {
							*numColumns = PhgBinParams[0].numAABins;
						}
						else if (curDim == 2) {
							*numSlices = PhgBinParams[0].numAABins;
							*sliceSize = PhgBinParams[0].aaCIsize;
						}
						else if (curDim > 2) { 
							ErStGeneric("This data has more than three dimensions");
						}
						curDim++;
					}
			
					break;
				
				
				case PhgBinEn_Z1:

					if (PhgBinParams[0].numZBins > 1) {

						if (curDim == 0) {
							*numColumns = PhgBinParams[0].numZBins;
						}
						else if (curDim == 2) {
							*numSlices = PhgBinParams[0].numZBins;
							*sliceSize = PhgBinParams[0].z1CIsize;
						}
						else if (curDim > 2) { 
							ErStGeneric("This data has more than three dimensions");
						}
						curDim++;
					}
					break;

				case PhgBinEn_Z2:

					if ((PhgBinParams[0].numZBins > 1) && 
							( !(PHG_IsSPECT() || PhgBinParams[0].isBinPETasSPECT) ) ) {

						if (curDim == 0) {
							*numColumns = PhgBinParams[0].numZBins;
						}
						else if (curDim == 2) {
							*numSlices = PhgBinParams[0].numZBins;
							*sliceSize = PhgBinParams[0].z2CIsize;
						}
						else if (curDim > 2) { 
							ErStGeneric("This data has more than three dimensions");
						}
						curDim++;
					}
			
					break;
					

				case PhgBinEn_Energy1:			

					if (PhgBinParams[0].numE1Bins > 1) {

						if (curDim == 0) {
							*numColumns = PhgBinParams[0].numE1Bins;
						}
						else if (curDim == 2) {
							*numSlices = PhgBinParams[0].numE1Bins;
							*sliceSize = PhgBinParams[0].energy1CIsize;
						}
						else if (curDim > 2) { 
							ErStGeneric("This data has more than three dimensions");
						}
						curDim++;
					}			
					break;
					
				case PhgBinEn_Energy2:

					if (PhgBinParams[0].numE2Bins > 1) {

						if (curDim == 0) {
							*numColumns = PhgBinParams[0].numE2Bins;
						}
						else if (curDim == 2) {
							*numSlices = PhgBinParams[0].numE2Bins;
							*sliceSize = PhgBinParams[0].energy2CIsize;
						}
						else if (curDim > 2) { 
							ErStGeneric("This data has more than three dimensions");
						}
						curDim++;
					}
					break;
					
				case PhgBinEn_Scatter1:
	
					if (PhgBinParams[0].numS1Bins > 1) {

						if (curDim == 0) {
							*numColumns = PhgBinParams[0].numS1Bins;
						}
						else if (curDim == 2) {
							*numSlices = PhgBinParams[0].numS1Bins;
							*sliceSize = PhgBinParams[0].scatter1CIsize;
						}
						else if (curDim > 2) { 
							ErStGeneric("This data has more than three dimensions");
						}
						curDim++;
					}
					break;
					
				case PhgBinEn_Scatter2:

					if (PhgBinParams[0].numS2Bins > 1) {

						if (curDim == 0) {
							*numColumns = PhgBinParams[0].numS2Bins;
						}
						else if (curDim == 2) {
							*numSlices = PhgBinParams[0].numS2Bins;
							*sliceSize = PhgBinParams[0].scatter2CIsize;
						}
						else if (curDim > 2) { 
							ErStGeneric("This data has more than three dimensions");
						}
						curDim++;
					}
					break;

				case PhgBinEn_THETA:

					if (PhgBinParams[0].numThetaBins > 1) {

						if (curDim == 0) {
							*numColumns = PhgBinParams[0].numThetaBins;
						}
						else if (curDim == 2) {
							*numSlices = PhgBinParams[0].numThetaBins;
							*sliceSize = PhgBinParams[0].thetaCIsize;
						}
						else if (curDim > 2) { 
							ErStGeneric("This data has more than three dimensions");
						}
						curDim++;
					}
					break;

				case PhgBinEn_PHI:

					if (PhgBinParams[0].numPHIBins > 1) {

						if (curDim == 0) {
							*numColumns = PhgBinParams[0].numPHIBins;
						}
						else if (curDim == 2) {
							*numSlices = PhgBinParams[0].numPHIBins;
							*sliceSize = PhgBinParams[0].phiCIsize;
						}
						else if (curDim > 2) { 
							ErStGeneric("This data has more than three dimensions");
						}
						curDim++;
					}
					break;
					

				case PhgBinEn_XR:

					if (PhgBinParams[0].numXRBins > 1) {

						if (curDim == 0) {
							*numColumns = PhgBinParams[0].numXRBins;
						}
						else if (curDim == 2) {
							*numSlices = PhgBinParams[0].numXRBins;
							*sliceSize = PhgBinParams[0].xrCIsize;
						}
						else if (curDim > 2) { 
							ErStGeneric("This data has more than three dimensions");
						}
						curDim++;
					}
					break;
					

				case PhgBinEn_YR:

					if (PhgBinParams[0].numYRBins > 1) {

						if (curDim == 0) {
							*numColumns = PhgBinParams[0].numYRBins;
						}
						else if (curDim == 2) {
							*numSlices = PhgBinParams[0].numYRBins;
							*sliceSize = PhgBinParams[0].yrCIsize;
						}
						else if (curDim > 2) { 
							ErStGeneric("This data has more than three dimensions");
						}
						curDim++;
					}
					break;
				
				
				/* Ignore others */
				case PhgBinEn_Null:
				case PhgBinEn_TOF:
				case PhgBinEn_Crystal1:
				case PhgBinEn_Crystal2:
					
					break;
				
			}
		}		
	
	/* Verify there were three dimensions */
	if (curDim < 3) {
		ErStGeneric("This data does not have three dimensions!");
		goto FAIL;
	}
	
	/* Completion */
	okay = true;
	FAIL:;
	return (okay);
	}

/*********************************************************************************
*
*			Name:			line3dTerminateForPHG
*
*			Summary:	Terminate the various modules.
*			Arguments:
*
*			Function return: none.
*
*********************************************************************************/
void line3dTerminateForPHG()	
{
				
	/* Terminate the collimation module if initialized */
	if (PHG_IsCollimateOnTheFly()) {
		ColTerminate();
	}
	/* Terminate the detection module if initialized */
	if (PHG_IsDetectOnTheFly()) {
		DetTerminate();
	}
	
	/* Terminate the photon tracking module */
	PhoTrkTerminate();
	
	/* Terminate the emission list manager */
	EmisListTerminate();
	
	/* Terminate the Productivity table manager */
	ProdTblTerminate(0);
	
	/* Terminate the sub object module */
	SubObjTerminate();
	
}

/*********************************************************************************
*
*			Name:			line3dInitForDataFile
*
*			Summary:	Terminate the various modules.
*			Arguments:
*
*			Function return: true - unless an error occurs.
*
*********************************************************************************/
Boolean line3dInitForDataFile(char *fileName, void **data, LbUsFourByte *numSlices,
			LbUsFourByte *sliceSize, LbUsFourByte *numColumns)	
{
	Boolean 		okay = false;
	FILE			*imageFile = 0;
	LbUsFourByte	hdrSize = 0;
	LbUsFourByte	numRows = 0;
	LbUsFourByte	fileSize = 0;
	LbUsFourByte	dataSize = 0;
	LbUsFourByte	numRead = 0;
	
	do { /* Process Loop */
		
		/* Open the image file */
		if ((imageFile = LbFlFileOpen(fileName, "r+b")) == 0) {
			sprintf(errString, "Unable to open image file\n'%s' (line3dInitForDataFile).", fileName);
			ErStFileError(errString);
			break;
		}
		
		/* See if they want To skip a header */
		hdrSize = LbInAskFourByte("Enter size of header to skip",
			1, false, false, 1, &canceled, 32768, 0, 0, 0);
		
		/* Get the number of rows */
		numRows = LbInAskFourByte("Enter the number of rows per slice",
			1, false, false, 1, &canceled, 128, 0, 0, 0);
		
		/* Get the number of slices */
		*numColumns = LbInAskFourByte("Enter the number of columns per slice",
			1, false, false, 1, &canceled, numRows, 0, 0, 0);
		
		/* Compute slice size */
		*sliceSize = numRows * *numColumns;
		
		/* Turn buffering off, it shouldn't be necessary but bugs have been found in the
			DG i/o library that make it so.
		*/
		setbuf(imageFile, 0);
	
		/* Seek to the end of the file */
		if (fseek(imageFile, 0, SEEK_END) != 0) {
			ErStFileError("Unable to seek to end of file (line3dInitForDataFile)");
			break;
		}

		/* Compute the file size */
		fileSize = ftell(imageFile);
		dataSize = fileSize - hdrSize;
	
		/* Seek to the beginning of the data */
		if (fseek(imageFile, hdrSize, SEEK_SET) != 0) {
			ErStFileError("Unable to seek to beginning of file (line3dInitForDataFile)");
			break;
		}
		
		/* Compute num slices based on bytes size, will need to be divided by element
			size below
		*/
		*numSlices = dataSize/(numRows * *numColumns);
		
		
		/* Allocate memory buffer for first file, using Numerical Recipes routine */
		if (Is_Float()){
		
			*numSlices /= sizeof(float);
			*data = LbMmAlloc(dataSize);
			if (*data == 0) {
				goto FAIL;
			}
			
			/* Read in the data */
			if ((numRead = fread(*data, dataSize, 1, imageFile)) != 1) {
				sprintf(errString, "Unable to read data for image file '%s'", fileName);
				ErStFileError(errString);
				break;
			}
		}
		else if (Is_Double()){
		
			*numSlices /= sizeof(double);
			
			*data = LbMmAlloc(dataSize);
			if (*data == 0) {
				goto FAIL;
			}
			
			/* Read in the data */
			if ((numRead = fread(*data, ((numRows * *numColumns) *sizeof(double)), 1, imageFile)) != 1) {
				sprintf(errString, "Unable to read data for image file '%s'", fileName);
				ErStFileError(errString);
				break;
			}
		}
		else {
		
			*numSlices /= sizeof(LbFourByte);
			
			*data = LbMmAlloc(dataSize);
			if (*data == 0) {
				goto FAIL;
			}
			
			/* Read in the data */
			if ((numRead = fread(*data, (dataSize), 1, imageFile)) != 1) {
				sprintf(errString, "Unable to read data for image file '%s'", fileName);
				ErStFileError(errString);
				break;
			}
		}
		
		
	okay = true;
	FAIL:;
	}while (false);
	
	if (imageFile != 0)
		fclose(imageFile);
		
	return (okay);			
	
}

/*********************************************************************************
*
*			Name:			line3dInitForPHG
*
*			Summary:	Terminate the various modules.
*			Arguments:
*
*			Function return: true - unless an error occurs.
*
*********************************************************************************/
Boolean line3dInitForPHG(char *fileName, LbUsFourByte *numSlices,
			LbUsFourByte *sliceSize, LbUsFourByte *numColumns)	
{
	Boolean okay = false;
	
	do { /* Process Loop */
		
		/* Initialize the PHG modules */
		if (line3DInitializePHG(fileName) == false) {
			goto FAIL;
		}
				
		/* Get the data */
		if (line3DGtImgData() == false) {
			goto FAIL;
		}
	
		/* Get index values from params */
		if (line3dGetIndexesForPHG(numSlices, sliceSize, numColumns) == false) {
			goto FAIL;
		}

	okay = true;
	FAIL:;
	}while (false);
	
	return (okay);			
	
}

/*********************************************************************************
*
*			Name:			line3dPrintForPHG
*
*			Summary:	Print values for PHG data.
*			Arguments:
*
*			Function return: none.
*
*********************************************************************************/
void line3dPrintForPHG(LbUsFourByte imageIndex)	
{
	
	/* Update the weight statistics if they created one */
	if (PhgBinParams[0].doWeights == true) {

		if ((PhgBinParams[0].sumAccordingToType == true) && (PhgBinParams[0].weight_image_type == PHG_BIN_WEIGHT_TYPE_R4)) {

			/* Update the statistics */
			LbInPrintf("%3.5e\t", ((float *)PhgBinData[0].weightImage)[imageIndex]);

		}
		else {
			
			/* Update the statistics */
			LbInPrintf("%3.5e\t", ((double *)PhgBinData[0].weightImage)[imageIndex]);
		}
	}

	/* Update the weight statistics if they created one */
	if (PhgBinParams[0].doWeightsSquared == true) {

		if ((PhgBinParams[0].sumAccordingToType == true) && (PhgBinParams[0].weight_image_type == PHG_BIN_WEIGHT_TYPE_R4)) {

			/* Update the statistic */
			LbInPrintf("%3.5e\t", ((float *)PhgBinData[0].weightSquImage)[imageIndex]);

		}
		else {
			
			/* Update the statistic */
			LbInPrintf("%3.5e\t", ((double *)PhgBinData[0].weightSquImage)[imageIndex]);
		}
	}
	
	/* Update the counts statistics if they created one */
	if (PhgBinParams[0].doCounts == true) {
	
		/* Update the counts statistics */
		switch(PhgBinParams[0].count_image_type){
			case PHG_BIN_COUNT_TYPE_I1:
				
				LbInPrintf("%d\t", ((LbUsOneByte *)PhgBinData[0].countImage)[imageIndex]);
				break;
				
			case PHG_BIN_COUNT_TYPE_I2:
					
				/* Update the statistic */
				LbInPrintf("%d\t", ((LbUsTwoByte *)PhgBinData[0].countImage)[imageIndex]);
				break;
				
			case PHG_BIN_COUNT_TYPE_I4:

				/* Update the statistic */
				LbInPrintf("%d\t", ((LbUsFourByte *)PhgBinData[0].countImage)[imageIndex]);
				break;
		}
		
	}
	
}

/*********************************************************************************
*
*			Name:			line3dPrintForDataFile
*
*			Summary:	Print values from a data.
*			Arguments:
*
*			Function return: none.
*
*********************************************************************************/
void line3dPrintForDataFile(LbUsFourByte imageIndex, void *data)	
{
	
	/* Allocate memory buffer for first file, using Numerical Recipes routine */
	if (Is_Float()){

		/* Update the statistics */
		LbInPrintf("%3.5e\t", ((float *)data)[imageIndex]);
	}
	else if (Is_Double()){

		/* Update the statistics */
		LbInPrintf("%3.5e\t", ((double *)data)[imageIndex]);

	}
	else {

		/* Update the statistics */
		LbInPrintf("%d\t", ((LbUsFourByte *)data)[imageIndex]);

	}
}

#undef LINE_3D
