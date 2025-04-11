/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1993-2013 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			makeindexfile.c
*     Revision Number:		1.8
*     Date last revised:	23 July 2013
*     Programmer:			Steven Vannoy
*     Date Originated:		March 16, 1993
*
*     Module Overview:	This program creates a file of indexes for the phg. It
*						asks for a file name, the number columns, and the
*						number of rows. It then request an index for each
*						voxel. The output file is in binary format compatible
*						with the phg.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*     Global variables defined:   
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
*			Revision date:		13 December 2011
*
*			Revision description:	Added 'box' option to object creation options
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		16 August 2011
*
*			Revision description:	Improved bounds checking and user experience
*
*********************************************************************************/

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbError.h"
#include "LbFile.h"
#include "LbInterface.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbHeader.h"

#include "Photon.h"
#include "ProdTbl.h"
#include "PhgMath.h"
#include "PhgParams.h"
#include "SubObj.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhoHFile.h"
#include "phg.h"
#include "PhgBin.h"

/* LOCAL TYPES */
typedef struct {
	double	x;		/* X-axis coordinate */
	double	y;		/* Y-axis coordinate */
	double	z;		/* Z-axis coordinate */
} point;

typedef struct {
	point	topLeft;		/* Top-Left point */
	point	bottomRight;	/* Bottom-Right point */
	double	length;			/* Length */
	double	width;			/* Width */
	double	height;			/* Height*/
} boxTy;

/* LOCAL CONSTANTS */
#define	CYLINDER			1
#define SPHERE				2
#define VOXEL				3
#define BOX					4
#define INDEX_SIZE			2

#define SUBOBJ_NUM_TISSUE_TRANSLATIONS	256		/* Entries in the translation tables */
#define SUBOBJ_NUM_ENERGY_BINS			1000								/* Number of energy bins */

/* The SIOW environment doesn't set argc properly. Hence,
	we have a special constant for checking if there 
	are any arguments to the program or not.
*/

#ifdef MPW
	#define NO_ARGS	0
#else
	#define NO_ARGS 1
#endif
/* LOCAL GLOBALS */
#define					NUM_TISSUE_NAMES_LETTERS	80
Boolean					canceled = false;								/* Did user cancel */
char					mifTissueNames[SUBOBJ_NUM_TISSUE_TRANSLATIONS][NUM_TISSUE_NAMES_LETTERS];
LbUsFourByte			*attenuationTransTbl = 0;	/* Attenuation index translation table */

/* LOCAL PROTOTYPES */
Boolean					makeindexfile(int argc, char *argv[]);
Boolean					allocIndexArrays(void);
Boolean					freeSubObject(Boolean didActivity, Boolean didAttenuation);
void					createObjects(Boolean doActivity);
Boolean					writeIndexFile(FILE *indexFile, Boolean doActivity);
Boolean					initializeIndexArrays(Boolean doActivity);
void					getPoint(point *pointPtr, Boolean *canceledPtr);
void					getBox(boxTy *boxPtr, Boolean *canceledPtr);
void					getCylinder(point *pointPtr, double *lengthPtr,
							double *x_radiusPtr, double *y_radiusPtr,
							Boolean *canceledPtr);
void					getSphere(point *pointPtr, double *radiusPtr,
							Boolean *canceledPtr);
Boolean					pointInSphere(point *pointPtr, point *sphereCntrPtr,
							double radius);
Boolean					pointInCylinder(point *pointPtr, point *cylinderCntrPtr, double xAxis_Radius,
							double yAxis_Radius, double length);
Boolean					pointInBox(point *pointPtr, boxTy *boxPtr);
Boolean					objectFilter(char *inputString, Boolean *isValid,
							char *errorPrompt, LbUserData userData);
void printMaterials(void);
							

/* FUNCTIONS */

/* FUNCTIONS */
/**********************
*	makeindexfile
*
*	Purpose:	Make the index files.
*
*	Result:	True unless an error occurs.
***********************/

Boolean makeindexfile(int argc, char *argv[])
{
	Boolean					okay = false;				/* Process flag */
	Boolean					didActivity;				/* Process flag */
	Boolean					didAttenuation;				/* Process flag */
	char					fileName[256];				/* Name of param file */
	char					prompt[512];				/* Prompt for user */
	char					errString[512];				/* Prompt for user */
	FILE					*indexFile;					/* The file we will create */
	
	do	/* Abort Loop */
	{
		
		
		/* Get input file path */
		if (argc == NO_ARGS) {
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
		else 
			strcpy(PhgRunTimeParams.PhgParamFilePath, argv[1]);
			
						
		/* Get our run-time parameters */
		if (!PhgGetRunTimeParams())
			break;
	
		/* Print out materials to help user */
		printMaterials();
		
		/* Clear processing flag */
		didActivity = false;
		didAttenuation = false;
		
		/* Repeat initialization for activity and attenuation files */
		{
			/* Build prompt string */
			sprintf(prompt, "\n\nBuild indexes for %s", PhgRunTimeParams.PhgSubObjActivityIndexFilePath);
			
			/* See if user wants to do activity file */
			if (LbInAskYesNo(prompt, LBINYes) == LBINYes) {
			
				/* Set processing flag */
				didActivity = true;
				
				/* Open the output file */
				if ((indexFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjActivityIndexFilePath, "wb")) == 0) {
					sprintf(errString, "Unable to open index file named '%s'.",
						PhgRunTimeParams.PhgSubObjActivityIndexFilePath);

					ErStFileError(errString);
					break;
				}
				
				/* Allocate the index arrays */
				if (!allocIndexArrays())
					break;
				
				/* Print out the list of materials */
				{
					LbUsFourByte	tissueIndex;
					
					LbInPrintf("\n\tYour material indexes and names are:\n");
					for (tissueIndex = 0; 
								tissueIndex < SUBOBJ_NUM_TISSUE_TRANSLATIONS; tissueIndex++) {
						if ((attenuationTransTbl[tissueIndex] != 0) || (tissueIndex == 0)) {
							if (strncmp(mifTissueNames[attenuationTransTbl[tissueIndex]], 
										"temp", 4)) {
								LbInPrintf("\t\tIndex = %d, name = %s\n", tissueIndex,
									mifTissueNames[attenuationTransTbl[tissueIndex]]);
							}
						}
					}
				}
				
				/* Initialize the indexes */
				if (!initializeIndexArrays(true))
					break;
				
				/* Create the "objects" */
				createObjects(true);
					
				/* Write the output file */
				if (!writeIndexFile(indexFile, true))
					break;
					
				/* Close the file */
				fclose(indexFile);
			}
			
			/* Build prompt string */
			sprintf(prompt, "\n\nBuild indexes for %s", PhgRunTimeParams.PhgSubObjAttenIndexFilePath);

			/* See if user wants to do attenuation file */
			if (LbInAskYesNo(prompt, LBINYes) == LBINYes) {
				
				/* Allocate the index arrays if necessary */
				if (!didActivity) {
					if (!allocIndexArrays())
						break;
				}
				
				/* Set flag */
				didAttenuation = true;
				
				
				/* Open the output file */
				if ((indexFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjAttenIndexFilePath, "wb")) == 0) {
					sprintf(errString, "Unable to open index file named '%s'.",
						PhgRunTimeParams.PhgSubObjAttenIndexFilePath);

					ErStFileError(errString);
					break;
				}
				
				/* Print out the list of materials again */
				{
					LbUsFourByte	tissueIndex;
					
					LbInPrintf("\n\tYour material indexes and names are:\n");
					for (tissueIndex = 0; 
								tissueIndex < SUBOBJ_NUM_TISSUE_TRANSLATIONS; tissueIndex++) {
						if ((attenuationTransTbl[tissueIndex] != 0) || (tissueIndex == 0)) {
							if (strncmp(mifTissueNames[attenuationTransTbl[tissueIndex]], 
										"temp", 4)) {
								LbInPrintf("\t\tIndex = %d, name = %s\n", tissueIndex,
									mifTissueNames[attenuationTransTbl[tissueIndex]]);
							}
						}
					}
				}
				
				/* Initialize the indexes */
				if (!initializeIndexArrays(false))
					break;
				
				/* Create the "objects" */
				createObjects(false);
					
				/* Write the output file */
				if (!writeIndexFile(indexFile, false))
					break;
					
				/* Close the file */
				fclose(indexFile);
			}
		}
		
		/* Free the subobject if necessary */
		if (didActivity || didAttenuation) {
			if (!freeSubObject(didActivity, didAttenuation))
				break;
		}
		
		okay = true;
		CANCEL:;
	} while (false);

	if (attenuationTransTbl != 0)
		LbMmFree((void **)&attenuationTransTbl);
	
	return (okay);
}


/*********************************************************************************
*
*			Name:		allocIndexArrays
*
*			Summary:	Allocate memory for subobject index arrays.
*
*			Arguments:
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean	allocIndexArrays()
{
	Boolean					okay = false;					/* Process flag */
	char					errString[1024];				/* Error string */
	LbUsFourByte			sliceIndex;						/* Current slice */
	
	do { /* Process Loop */

		/* Verify that the slices are contiguous and of the same x/y dimension */
		for (sliceIndex = 0; sliceIndex < (SubObjNumSlices-1); sliceIndex++) {
		
			/* Verify xMax of current slice is the same as next slice */
			if (SubObjObject[sliceIndex].xMax != SubObjObject[sliceIndex+1].xMax){
				
				/* Build error string */
				sprintf(errString, "Slices must have same x/y dimensions!\nSlice #%ld, has xMax = %3.2f and slice %ld has xMax = %3.2f",
					(unsigned long)sliceIndex, SubObjObject[sliceIndex].xMax, 
					(unsigned long)sliceIndex+1, SubObjObject[sliceIndex+1].xMax);
				ErStGeneric(errString);
				goto FAIL;
			}
			/* Verify xMin of current slice is the same as next slice */
			if (SubObjObject[sliceIndex].xMin != SubObjObject[sliceIndex+1].xMin){
				/* Build error string */
				sprintf(errString, "Slices must have same x/y dimensions!\nSlice #%ld, has xMin = %3.2f and slice %ld has xMin = %3.2f",
					(unsigned long)sliceIndex, SubObjObject[sliceIndex].xMin, 
					(unsigned long)sliceIndex+1, SubObjObject[sliceIndex+1].xMin);
				ErStGeneric(errString);
				goto FAIL;
			}
			/* Verify yMax of current slice is the same as next slice */
			if (SubObjObject[sliceIndex].yMax != SubObjObject[sliceIndex+1].yMax){
				/* Build error string */
				sprintf(errString, "Slices must have same x/y dimensions!\nSlice #%ld, has yMax = %3.2f and slice %ld has yMax = %3.2f",
					(unsigned long)sliceIndex, SubObjObject[sliceIndex].yMax, 
					(unsigned long)sliceIndex+1, SubObjObject[sliceIndex+1].yMax);
				ErStGeneric(errString);
				goto FAIL;
			}
			/* Verify yMin of current slice is the same as next slice */
			if (SubObjObject[sliceIndex].yMin != SubObjObject[sliceIndex+1].yMin){
				/* Build error string */
				sprintf(errString, "Slices must have same x/y dimensions!\nSlice #%ld, has yMin = %3.2f and slice %ld has yMin = %3.2f",
					(unsigned long)sliceIndex, SubObjObject[sliceIndex].yMin, 
					(unsigned long)sliceIndex+1, SubObjObject[sliceIndex+1].yMin);
				ErStGeneric(errString);
				goto FAIL;
			}
			
			/* Verify slices are contiguous */
			if (SubObjObject[sliceIndex].zMax != SubObjObject[sliceIndex+1].zMin){
				/* Build error string */
				sprintf(errString, "Slices must be contiguous in axial direction!\nSlice #%ld, has zMax = %3.2f and slice %ld has zMin = %3.2f",
					(unsigned long)sliceIndex, SubObjObject[sliceIndex].zMax, 
					(unsigned long)sliceIndex+1, SubObjObject[sliceIndex+1].zMin);
				ErStGeneric(errString);
				goto FAIL;
			}
				
		}

		/* Verify xMax of current slice is the same as yMax */
		if (SubObjObject[0].xMax != SubObjObject[0].yMax){
			sprintf(errString,"\nCurrently, the object cylinder must be round.\n\t"
				"In slice %d, you have specified x-max as %3.2f and y-max as %3.2f\n\t"
				"These values must be equal.\n", 0, SubObjObject[0].xMax,
				SubObjObject[0].yMax);
			ErStGeneric(errString);
			goto FAIL;
		}

		/* Verify xMin of current slice is the same as yMin */
		if (SubObjObject[0].xMin != SubObjObject[0].yMin){
			sprintf(errString,"\nCurrently, the object cylinder must be round.\n\t"
				"In slice %d, you have specified x-min as %3.2f and y-min as %3.2f\n\t"
				"These values must be equal.\n", 0, SubObjObject[0].xMin,
				SubObjObject[0].yMin);
			ErStGeneric(errString);
			goto FAIL;
		}

		/* Verify object cylinder is round */
		if ((SubObjObject[0].xMin != -SubObjObject[0].xMax)){
			sprintf(errString,"\nCurrently, the object cylinder must be round.\n\t"
				"In slice %d, you have specified x-min as %3.2f and x-max as %3.2f\n\t"
				"These values must be symmetric.\n", 0, SubObjObject[0].xMin,
				SubObjObject[0].xMax);
			ErStGeneric(errString);
			goto FAIL;
		}

		/* Verify object cylinder is centered */
		if ((SubObjObject[0].xMin + SubObjObject[0].yMax) != 0){
			sprintf(errString,"\nCurrently, the object cylinder must be centered.\n\t"
				"In slice %d, you have specified x-min as %3.2f and x-max as %3.2f\n\t"
				"These values must be symmetric.\n", 0, SubObjObject[0].xMin,
				SubObjObject[0].xMax);
			ErStGeneric(errString);
			goto FAIL;
		}
		
		/* Loop through the parameters at this level */
		for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {

			/* Compute slice width, height, and depth */
			SubObjObject[sliceIndex].sliceWidth = fabs(
				SubObjObject[sliceIndex].xMax - 
				SubObjObject[sliceIndex].xMin);
				
			SubObjObject[sliceIndex].sliceHeight = fabs(
				SubObjObject[sliceIndex].yMax - 
				SubObjObject[sliceIndex].yMin);
			
			SubObjObject[sliceIndex].sliceDepth = fabs(
				SubObjObject[sliceIndex].zMax - 
				SubObjObject[sliceIndex].zMin);

			SubObjObject[sliceIndex].actVoxelWidth = SubObjObject[sliceIndex].sliceWidth/
				SubObjObject[sliceIndex].actNumXBins;

			SubObjObject[sliceIndex].attVoxelWidth = SubObjObject[sliceIndex].sliceWidth/
				SubObjObject[sliceIndex].attNumXBins;

			SubObjObject[sliceIndex].actVoxelHeight = SubObjObject[sliceIndex].sliceHeight/
				SubObjObject[sliceIndex].actNumYBins;

			SubObjObject[sliceIndex].attVoxelHeight = SubObjObject[sliceIndex].sliceHeight/
				SubObjObject[sliceIndex].attNumYBins;
		
			/* Allocate array for indexes */
			if ((SubObjObject[sliceIndex].activityArray = (SubObjActVoxelTy *) LbMmAlloc(
					(sizeof(SubObjActVoxelTy) * 
					(SubObjObject[sliceIndex].actNumXBins * SubObjObject[sliceIndex].actNumYBins)))) == 0) {
				
				goto FAIL;
			}
		
			/* Allocate array for indexes */
			if ((SubObjObject[sliceIndex].attenuationArray = (LbUsFourByte *) LbMmAlloc(
					(sizeof(LbUsFourByte) * 
					(SubObjObject[sliceIndex].attNumXBins * SubObjObject[sliceIndex].attNumYBins)))) == 0) {
				
				goto FAIL;
			}
		}
		okay = true;
		FAIL:;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:		freeSubObject
*
*			Summary:	Free memory for subobject index arrays.
*
*			Arguments:
*					Boolean	didActivity 	- Activity object was allocated.
*					Boolean	didAttenuation	- Attenuation object was allocated.
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean	freeSubObject(Boolean didActivity, Boolean didAttenuation)
{
	Boolean					okay = false;					/* Process flag */
	LbUsFourByte			sliceIndex;						/* Current slice */
	
	do { /* Process Loop */
	
		/* Loop through the parameters at this level */
		for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
			
			if (didActivity)
				LbMmFree((void **) &(SubObjObject[sliceIndex].activityArray));
		
			if (didAttenuation)
				LbMmFree((void **) &(SubObjObject[sliceIndex].attenuationArray));
		}
		
		LbMmFree((void **) &SubObjObject);
		okay = true;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:		initializeIndexArrays
*
*			Summary:	Initialize subobject index arrays.
*
*			Arguments:
*					Boolean	doActivity	- We're doing activity voxels
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean	initializeIndexArrays(Boolean doActivity)
{
	Boolean					okay = false;					/* Process flag */
	char					prompt[256];					/* String for prompting the user */
	char					fileName[256];					/* String for storing the file name */
	char					*indexBuffer;					/* Index buffer */
	FILE					*indexFile;						/* Our index file */
	LbUsFourByte			indexOffset;					/* Offset into file for data */
	LbUsFourByte			interSliceOffset;				/* Offset between slices for data */
	LbUsFourByte			indexSize;						/* Size of index */
	LbUsFourByte			numVoxels;						/* Number of voxels for this slice */
	LbUsFourByte			voxelIndex;						/* Current voxel */
	LbUsFourByte			fillValue;						/* Value to fill indexes with */
	LbUsFourByte			sliceIndex;						/* Current slice */
	
	do { /* Process Loop */

		/* See if they want to fill the entire object */
		if (LbInAskYesNo("\nWould you like to fill the entire object with a single value?", LBINYes) == LBINYes) {
	
			/* Get fill value for voxels */
			fillValue = LbInAskFourByteInRange("Enter fill value for object",
				false, &canceled, 0, 255);
			/* Loop through the parameters at this level */
			for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
				
				/* Compute number of voxels */
				if (doActivity)
					numVoxels = SubObjObject[sliceIndex].actNumXBins *
						SubObjObject[sliceIndex].actNumYBins;
				else
					numVoxels = SubObjObject[sliceIndex].attNumXBins *
						SubObjObject[sliceIndex].attNumYBins;
		
				/* Fill in the voxels with the given fill value */
				for (voxelIndex = 0; voxelIndex < numVoxels; voxelIndex++) {
					if (doActivity)
						SubObjObject[sliceIndex].activityArray[voxelIndex] = fillValue;
					else
						SubObjObject[sliceIndex].attenuationArray[voxelIndex] = fillValue;
					
				}
			}
			
			/*if (!doActivity)*/
			LbInPrintf("Object filled with '%s'\n", mifTissueNames[attenuationTransTbl[fillValue]]);
			
		}
		
		/* See if they want to initialize the entire object from a file */
		if (LbInAskYesNo("\nWould you like to fill the entire object from a single data file?", LBINNo) == LBINYes) {
		
				/* Ask for the file name */
				LbInAsk("Enter name of data file", 0, false,
						&canceled, 0, 0, 0, 0,
						fileName);
			
				/* Bolt if we canceled */
				if (canceled) {
					ErStCancel("User canceled.");
					goto CANCEL;
				}
		
				/* Open the data file */
				if ((indexFile = LbFlFileOpen(fileName, "rb")) == 0) {
					ErStFileError("Unable to open index file.\n");
					goto FAIL;
				}
				
				/* Get the offset to the data */
				indexOffset = LbInAskFourByte("Enter offset to start of data (in bytes)", 0,
					false, false, 0, &canceled, 0, 0, 0, 0);

				/* Bolt if we canceled */
				if (canceled) {
					ErStCancel("User canceled.");
					goto CANCEL;
				}
				
				/* Get the offset to the data */
				if (SubObjNumSlices > 1) {
					interSliceOffset = LbInAskFourByte("Enter offset between slices of data (in bytes)", 0,
						false, false, 0, &canceled, 0, 0, 0, 0);
				}
				else {
					interSliceOffset = 0;
				}
				
				/* Bolt if we canceled */
				if (canceled) {
					ErStCancel("User canceled.");
					goto CANCEL;
				}
				
				/* Get the size of the indexes */
				indexSize = LbInAskTwoByte("Enter size of indexes (as in 1 byte, 2 bytes, 4 bytes = 1,2,or 4)", 0,
					false, false, 0, &canceled, 0, 0, 0, 0);

				/* Bolt if we canceled */
				if (canceled) {
					ErStCancel("User canceled.");
					goto CANCEL;
				}
				
				/* Loop through the parameters at this level */
				for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
					
					/* Compute number of voxels */
					if (doActivity)
						numVoxels = SubObjObject[sliceIndex].actNumXBins *
							SubObjObject[sliceIndex].actNumYBins;
					else
						numVoxels = SubObjObject[sliceIndex].attNumXBins *
							SubObjObject[sliceIndex].attNumYBins;
					
					/* Allocate a buffer for the indexes */
					if ((indexBuffer = (char *) LbMmAlloc(numVoxels * indexSize)) == 0)
						goto FAIL;
					
					/* Position to data */
					if (indexOffset > 0) {
						if (fseek(indexFile, indexOffset, SEEK_SET) != 0) {
							ErStFileError("Unable to position to index data.\n");
							goto FAIL;
						}
					}
					
					/* If there will be more slices then set indexOffset to interSliceOffset
						for subsequent seeks
					 */
					if (SubObjNumSlices > 1) {
					
						indexOffset = interSliceOffset;
					}
					
					/* Read in the indexes */
					if (fread(indexBuffer, indexSize, numVoxels, indexFile)
							!= numVoxels) {
							
						ErStFileError("Unable to read index data.");
						goto FAIL;
					}
			
					/* Fill in the voxels with the given fill value */
					for (voxelIndex = 0; voxelIndex < numVoxels; voxelIndex++) {
					
						switch(indexSize) {
						
							case 1:
								fillValue = (LbUsFourByte) *(indexBuffer+voxelIndex);
								break;
								
							case 2:
								fillValue = *((LbUsTwoByte *)indexBuffer+voxelIndex);
								break;
								
							case 4:
								fillValue = *((LbUsFourByte *)indexBuffer+voxelIndex);
								break;
								
							default:
								ErAbort("\nYou gave me an invalid index size.");
								break;
						}
						
						/* Fill the array */
						if (doActivity)
							SubObjObject[sliceIndex].activityArray[voxelIndex] = fillValue;
						else
							SubObjObject[sliceIndex].attenuationArray[voxelIndex] = fillValue;
			
					}
					
					/* Free the index buffer */
					LbMmFree((void **) &indexBuffer);
					
				}
				
				/* Close the file */
				fclose(indexFile);
		}
		
		/* See if they want to initialize specific slices */
		if (LbInAskYesNo("\nWould you like to initialize specific slice values?", LBINNo) == LBINNo) {
			goto FINISHED;
		}
		
		/* Loop through the parameters at this level */
		for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
			
			/* Compute number of voxels */
			if (doActivity)
				numVoxels = SubObjObject[sliceIndex].actNumXBins *
					SubObjObject[sliceIndex].actNumYBins;
			else
				numVoxels = SubObjObject[sliceIndex].attNumXBins *
					SubObjObject[sliceIndex].attNumYBins;
			
			/* See if they have a prepared file for current slice */
			sprintf(prompt, "Do you have a prepared file for slice %ld? ", (unsigned long)sliceIndex);
			if (LbInAskYesNo(prompt, LBINNo) == LBINYes) {
			
					/* Ask for the file name */
					LbInAsk("Enter data file for slice ", 0, false,
							&canceled, 0, 0, 0, 0,
							fileName);
				
					/* Bolt if we canceled */
					if (canceled) {
						ErStCancel("User canceled.");
						goto CANCEL;
					}
			
					/* Open the data file */
					if ((indexFile = LbFlFileOpen(fileName, "rb")) == 0) {
						ErStFileError("Unable to open index file.\n");
						goto FAIL;
					}
					
					/* Get the offset to the data */
					indexOffset = LbInAskFourByte("Enter offset to data", 0,
						false, false, 0, &canceled, 0, 0, 0, 0);
	
	
					/* Bolt if we canceled */
					if (canceled) {
						ErStCancel("User canceled.");
						goto CANCEL;
					}
					
					/* Get the size of the indexes */
					indexSize = LbInAskTwoByte("Enter size of indexes (as in 1 byte, 2 bytes, 4 bytes = 1,2,or 4)", 0,
						false, false, 0, &canceled, 0, 0, 0, 0);
	
	
					/* Bolt if we canceled */
					if (canceled) {
						ErStCancel("User canceled.");
						goto CANCEL;
					}
					
					/* Allocate a buffer for the indexes */
					if ((indexBuffer = (char *) LbMmAlloc(numVoxels * indexSize)) == 0)
						goto FAIL;
					
					/* Position to data */
					if (fseek(indexFile, indexOffset, SEEK_SET) != 0) {
						ErStFileError("Unable to position to index data.\n");
						goto FAIL;
					}
					
					/* Read in the indexes */
					if (fread(indexBuffer, indexSize, numVoxels, indexFile)
							!= numVoxels) {
							
						ErStFileError("Unable to read index data.");
						goto FAIL;
					}
			
					/* Fill in the voxels with the given fill value */
					for (voxelIndex = 0; voxelIndex < numVoxels; voxelIndex++) {
					
						switch(indexSize) {
						
							case 1:
								fillValue = (LbUsFourByte) *(indexBuffer+voxelIndex);
								break;
								
							case 2:
								fillValue = *((LbUsTwoByte *)indexBuffer+voxelIndex);
								break;
								
							case 4:
								fillValue = *((LbUsFourByte *)indexBuffer+voxelIndex);
								break;
								
							default:
								ErAbort("\nYou gave me an invalid index size.");
								break;
						}
						
						/* Fill the array */
						if (doActivity)
							SubObjObject[sliceIndex].activityArray[voxelIndex] = fillValue;
						else
							SubObjObject[sliceIndex].attenuationArray[voxelIndex] = fillValue;
			
					}
					
					/* Free the index buffer */
					LbMmFree((void **) &indexBuffer);
					
					/* Close the file */
					fclose(indexFile);
			}
			else {
				/* Get fill value for voxels */
				fillValue = LbInAskFourByteInRange("\nEnter fill value for slice",
					false, &canceled, 0, 255);
				LbInPrintf("\n");
		
				/* Bolt if we canceled */
				if (canceled) {
					ErStCancel("User canceled.");
					goto CANCEL;
				}
		
				/* Fill in the voxels with the given fill value */
				for (voxelIndex = 0; voxelIndex < numVoxels; voxelIndex++) {
					if (doActivity)
						SubObjObject[sliceIndex].activityArray[voxelIndex] = fillValue;
					else
						SubObjObject[sliceIndex].attenuationArray[voxelIndex] = fillValue;
					
				}
				/*if (!doActivity)*/
					LbInPrintf("Slice filled with '%s'\n", mifTissueNames[attenuationTransTbl[fillValue]]);
			}
		}
		
		FINISHED:;
		okay = true;
		FAIL:;
		CANCEL:;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:		createObjects
*
*			Summary:	Allow user to create "objects".
*
*			Arguments:
*				Boolean	doActivity	- Are we doing the activity or the attenuation object
*
*			Function return: None.
*
*********************************************************************************/
void	createObjects(Boolean doActivity)
{
	double					objectRadius;					/* Radius of objects */
	double					xAxis_Radius;					/* Radius in x axis */
	double					yAxis_Radius;					/* Radius in y axis */
	double					objectLength;					/* Length of object */
	double					voxelWidth;						/* Width of voxel */
	double					voxelHeight;					/* Height of voxel */
	LbTwoByte				createAnObject;					/* Flag for creating more objects */
	LbUsFourByte			objectType;						/* Type of object they want */
	LbUsFourByte			xIndex;							/* Current x index */
	LbUsFourByte			yIndex;							/* Current y index */
	LbUsFourByte			fillValue;						/* Value to fill indexes with */
	LbUsFourByte			sliceIndex;						/* Current slice */
	LbUsFourByte			fillCount;						/* Count of voxels affected by object fill */
	LbUsFourByte			numXBins;						/* Number of X bins */
	LbUsFourByte			numYBins;						/* Number of Y bins */
	LbUsFourByte			*array;							/* The activity or attenuation */
	point					centerPoint;					/* Center point of object */
	point					voxelCenter;					/* Center of voxel */
	char					userStr[128];					/* String result for user choice */
	LbTwoByte				result;							/* Result from prompt */
	boxTy					box;							/* The object box */
	
	do { /* Process Loop */
	
			
		/* See if they want to create an object */
		createAnObject = LbInAskYesNo("\nWould you like to create an 'object'? ",
			LBINYes);
		
		/* Clear default object type */
		objectType = CYLINDER;
		
		/* Loop while the want to create an object */
		while (createAnObject == LBINYes) {
			
			/* See what type of object they want to create */
			result = LbInAskFiltered(0, "Select a type ", objectType, &canceled,
				"cylinder", "sphere", "voxel", "box",
				userStr, objectFilter, 0, (LbUserData) &objectType);
			
			if (result != LBINUserChoice)
				objectType = result;
			
			
			/* Get remaining info based on object type */
			switch (objectType) {
				case CYLINDER:
					/* They chose a cylinder so get its parameters */
					getCylinder(&centerPoint, &objectLength, &xAxis_Radius,
						&yAxis_Radius, &canceled);
						
					/* See if the user canceled */
					if (canceled == true)
						goto CANCEL;


					break;
					
				case SPHERE:
					/* They chose a sphere so get its parameters */
					getSphere(&centerPoint, &objectRadius, &canceled);
						
					/* See if the user canceled */
					if (canceled == true)
						goto CANCEL;

					
					break;
					
				case VOXEL:
					/* They chose a voxel so get its parameters */
					LbInPrintf("\n");
					getPoint(&centerPoint, &canceled);
						
					/* See if the user canceled */
					if (canceled == true)
						goto CANCEL;

					break;
					
				case BOX:
					/* They chose a voxel so get its parameters */
					getBox(&box, &canceled);
						
					/* See if the user canceled */
					if (canceled == true)
						goto CANCEL;

					break;
					
				default:
					/* This shouldn't happen */
					ErAbort("(createObjects): Got invalid response from LbInAsk for object type.");
					break;
			}
			
			
			/* Get fill value for voxels */
			fillValue = LbInAskFourByteInRange("\nEnter fill value for object",
				false, &canceled, 0, 255);
				

			/* Bolt if we canceled */
			if (canceled) {
				ErStCancel("User canceled.");
				goto CANCEL;
			}

			/* If object is  a point, compute the x/y/z indexes and fill it */
			if (objectType == VOXEL) {
				/* Get z axis index */
				if ((centerPoint.z < SubObjObject[0].zMin) || 
						(centerPoint.z > SubObjObject[SubObjNumSlices-1].zMax)) {
					LbInPrintf( "\n\nVoxel z coordinate outside z-axis boundaries.\n");
					goto CANCEL;
				}
				else {
					for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
						
						/* See if the position is in this slice */
						if (centerPoint.z < SubObjObject[sliceIndex].zMax) {
							if (centerPoint.z < SubObjObject[sliceIndex].zMin) {
								LbInPrintf( 
									"\n\nVoxel coordinate outside z boundary (z-axis).\n");
								goto CANCEL;
							}
							break;
						}
					}
				}
				
				if (doActivity) {
					array = SubObjObject[sliceIndex].activityArray;
					numXBins = SubObjObject[sliceIndex].actNumXBins;
					numYBins = SubObjObject[sliceIndex].actNumYBins;
				}
				else {
					array = SubObjObject[sliceIndex].attenuationArray;
					numXBins = SubObjObject[sliceIndex].attNumXBins;
					numYBins = SubObjObject[sliceIndex].attNumYBins;
				}
					
				/* Calculate the x index */
				if ((centerPoint.x < SubObjObject[sliceIndex].xMin) || 
						(centerPoint.x > SubObjObject[sliceIndex].xMax)) {
					LbInPrintf( "\n\nVoxel x coordinate outside x-axis boundaries.\n");
					goto CANCEL;
				}
				else {
					if (centerPoint.x == SubObjObject[sliceIndex].xMax){
							xIndex = numXBins - 1;
					}
					else {
						xIndex = (LbUsFourByte) 
							(((centerPoint.x - SubObjObject[sliceIndex].xMin)
							/SubObjObject[sliceIndex].sliceWidth) *
							numXBins);
					}
					
					/* Verify we are within bounds */
					if (xIndex >= numXBins) {
						LbInPrintf( "\n\nVoxel coordinate outside object (x-axis).\n");
						goto CANCEL;
					}
				}
				
				/* Calculate the y index */
				if ((centerPoint.y < SubObjObject[sliceIndex].yMin) || 
						(centerPoint.y > SubObjObject[sliceIndex].yMax)) {
					LbInPrintf( "\n\nVoxel y coordinate outside y-axis boundaries.\n");
					goto CANCEL;
				}
				else {
					if (centerPoint.y == SubObjObject[sliceIndex].yMin){
						yIndex = numYBins - 1;
					}
					else {
						yIndex = (LbUsFourByte) 
							((((-1 * centerPoint.y) + SubObjObject[sliceIndex].yMax)
							/SubObjObject[sliceIndex].sliceHeight) *
							numYBins);
					}
					
					/* Verify we are within bounds */
					if (yIndex >= numYBins) {
						LbInPrintf( "\n\nVoxel coordinate outside object (y-axis).");
						goto CANCEL;
					}
				}
				
				
				/* Set the voxel value */
				array[(yIndex *
					numXBins) + xIndex] = fillValue;
				
				/*if (!doActivity)*/
					LbInPrintf( "Voxel filled at (x, y, z) = (%ld, %ld, %ld) with %s.\n",
						(unsigned long)xIndex, (unsigned long)yIndex, (unsigned long)sliceIndex, 
						mifTissueNames[attenuationTransTbl[fillValue]]);
			}
			else {
				/* Clear fill count */
				fillCount = 0;
				
				/* Fill in the voxels with the given fill value */
				for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
				
					/* Compute voxel center in z axis */
					voxelCenter.z = SubObjObject[sliceIndex].zMin + 
						(SubObjObject[sliceIndex].sliceDepth/2);
					
					if (doActivity) {
						array = SubObjObject[sliceIndex].activityArray;
						voxelWidth = SubObjObject[sliceIndex].actVoxelWidth;
						voxelHeight = SubObjObject[sliceIndex].actVoxelHeight;
						numXBins = SubObjObject[sliceIndex].actNumXBins;
						numYBins = SubObjObject[sliceIndex].actNumYBins;
					}
					else {
						array = SubObjObject[sliceIndex].attenuationArray;
						voxelWidth = SubObjObject[sliceIndex].attVoxelWidth;
						voxelHeight = SubObjObject[sliceIndex].attVoxelHeight;
						numXBins = SubObjObject[sliceIndex].attNumXBins;
						numYBins = SubObjObject[sliceIndex].attNumYBins;
					}
					
					/* Traverse y indexes */
					for (yIndex = 0; yIndex < numYBins; yIndex++) {
					
						/* Compute center of voxel in y axis */
						voxelCenter.y = SubObjObject[sliceIndex].yMax -
							(voxelHeight * yIndex) -
							(voxelHeight/2);
							
						for (xIndex = 0; xIndex < numXBins; xIndex++) {
					
							/* Compute center of voxel in x axis */
							voxelCenter.x = SubObjObject[sliceIndex].xMin +
								(voxelWidth * xIndex) +
								(voxelWidth/2);
							
							/* See if voxel is in object */
							switch (objectType) {
							
								case CYLINDER:
									/* See if current voxel is in the cylinder */
									if (pointInCylinder(&voxelCenter, &centerPoint,
											xAxis_Radius,
											yAxis_Radius, objectLength)) {
										
										/* Set the voxel value */
										array[(yIndex *
											numXBins) + xIndex] = fillValue;											
											
										/* Increment our fill count */
										fillCount++;
									}
									break;
									
								case SPHERE:
									/* See if current voxel is in the shere */
									if (pointInSphere(&voxelCenter, &centerPoint,
										objectRadius)) {
											
										/* Set the voxel value */
										array[(yIndex *
											numXBins) + xIndex] = fillValue;											
											
										/* Increment out fill count */
										fillCount++;
									}
									break;
									
								case BOX:
									/* See if current voxel is in the box */
									if (pointInBox(&voxelCenter, &box)) {
											
										/* Set the voxel value */
										array[(yIndex *
											numXBins) + xIndex] = fillValue;											
											
										/* Increment out fill count */
										fillCount++;
									}
									break;
									
								default:
									break;
							}
						}
					}
				}
				
				/* Tell the user how many voxels were filled */
				/*if (!doActivity)*/
					LbInPrintf("Filled %ld voxels with %s.\n", 
					(unsigned long)fillCount, 
					mifTissueNames[attenuationTransTbl[fillValue]]);
			}
			CANCEL:;

			/* See if they want to create an object */
			createAnObject = LbInAskYesNo("\nWould you like to create another 'object'? ",
				LBINNo);


		}
	} while (false);
}

/*********************************************************************************
*
*			Name:		writeIndexFile
*
*			Summary:	Write the index data to an index file.
*
*			Arguments:
*				FILE	*outputFile - The index file to write.
*				Boolean	doActivity	- Are we doing activity or attenuation indexes
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean	writeIndexFile(FILE *outputFile, Boolean doActivity)
{
	Boolean					okay = false;					/* Process flag */
	LbUsFourByte			sliceIndex;						/* Current slice */
	LbUsFourByte			numVoxels;						/* Number of voxels in this slice */
	LbUsFourByte			*array;
	
	do { /* Process Loop */
	
			

		/* Loop through the slices */
		for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {

			/* Point to the write data */
			if (doActivity) {
				array = SubObjObject[sliceIndex].activityArray;
				
				numVoxels = SubObjObject[sliceIndex].actNumXBins *
					SubObjObject[sliceIndex].actNumYBins;
			}
			else {
				array = SubObjObject[sliceIndex].attenuationArray;
				
				numVoxels = SubObjObject[sliceIndex].attNumXBins *
					SubObjObject[sliceIndex].attNumYBins;
			}
				
			/* Write out the array */	
			if (fwrite(array,
					sizeof(LbUsFourByte) * numVoxels, 1, outputFile) != 1) {
				
				ErStFileError("\nFailed to write activity voxel indexes.\n");
				goto FAIL;
			}
		}
		
		okay = true;
		FAIL:;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:		getPoint
*
*			Summary:	Ask the user for three-dimensional coordinates of a
*						"point".
*
*			Arguments:
*				point		*pointPtr	- Point coordinates.
*				Boolean		*canceledPtr	- Did user cancel.
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
void	getPoint(point *pointPtr, Boolean *canceledPtr)
{
	/* Get x-axis coordinate */
	pointPtr->x = LbInAskDouble("Enter x-axis coordinate ", 0, true,
		2, false, 0, canceledPtr, 0, 0, 0, 0);

	/* Get y-axis coordinate */
	pointPtr->y = LbInAskDouble("Enter y-axis coordinate ", 0, true,
		2, false, 0, canceledPtr, 0, 0, 0, 0);
		
	/* Get z-axis coordinate */
	pointPtr->z = LbInAskDouble("Enter z-axis coordinate ", 0, true,
		2, false, 0, canceledPtr, 0, 0, 0, 0);
}

/*********************************************************************************
*
*			Name:		getCylinder
*
*			Summary:	Ask the user for three-dimensional coordinates of a
*						"cylinder".
*
*			Arguments:
*				point			*pointPtr		- Center point coordinates.
*				double			*lengthPtr		- Length of the cylinder.
*				double			*x_radiusPtr	- Radius along x-axis.
*				double			*y_radiusPtr	- Radius along y-axis.
*				Boolean			*canceledPtr	- Did user cancel.
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
void	getCylinder(point *pointPtr, double *lengthPtr,
			double *x_radiusPtr, double *y_radiusPtr,
			Boolean *canceledPtr)
{
	do { /* Process Loop */
		/* Get center point */
		{
			/* Prompt for point coordinates */
			LbInPrintf( "\nPlease select a center point for the cylinder\n");
			
			/* Get the point */
			getPoint(pointPtr, canceledPtr);
			if (*canceledPtr == true)
				break;
		}				
		
		/* Get the length */
		*lengthPtr = LbInAskDouble("Enter length along z-axis ", 0, true,
			2, false, 0, canceledPtr, 0, 0, 0, 0);
		
		/* Get radius in x axis */
		*x_radiusPtr = LbInAskDouble("Enter x-axis radius ", 0, true,
			2, false, 0, canceledPtr, 0, 0, 0, 0);
		
		/* Get radius in y axis */
		*y_radiusPtr = LbInAskDouble("Enter y-axis radius ", 0, true,
			2, false, 0, canceledPtr, 0, 0, 0, 0);
	} while (false);
}

/*********************************************************************************
*
*			Name:		getBox
*
*			Summary:	Ask the user for three-dimensional coordinates of a
*						"boxTy".
*
*			Arguments:
*				boxTy			*boxPtr			- Box point coordinates.
*				Boolean			*canceledPtr	- Did user cancel.
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
void	getBox(boxTy *boxPtr, Boolean *canceledPtr)

{
	point		minPoint;		/* Corner of box with smallest coordinates */
	point		maxPoint;		/* Corner of box with largest coordinates */
	
	do { /* Process Loop */
		/* Prompt for minimum point coordinates */
		LbInPrintf( "\nEnter minimum corner...\n");
		
		/* Get the point */
		getPoint(&minPoint, canceledPtr);
		if (*canceledPtr == true)
			break;
		
		/* Prompt for maximum point coordinates */
		LbInPrintf( "\nEnter maximum corner...\n");
		
		/* Get the point */
		getPoint(&maxPoint, canceledPtr);
		if (*canceledPtr == true)
			break;
		
		/* Fill in the boxTy values */
		boxPtr->topLeft.x = minPoint.x;
		boxPtr->topLeft.y = maxPoint.y;
		boxPtr->topLeft.z = minPoint.z;
		boxPtr->bottomRight.x = maxPoint.x;
		boxPtr->bottomRight.y = minPoint.y;
		boxPtr->bottomRight.z = minPoint.z;
		boxPtr->length = maxPoint.z - minPoint.z;
		boxPtr->width = maxPoint.x - minPoint.x;
		boxPtr->height = maxPoint.y - minPoint.y;
		
	} while (false);
}

/*********************************************************************************
*
*			Name:		getSphere
*
*			Summary:	Ask the user for three-dimensional coordinates of a
*						"sphere".
*
*			Arguments:
*				point			*pointPtr		- Center point coordinates.
*				double			*radiusPtr		- Radius.
*				Boolean			*canceledPtr	- Did user cancel.
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
void	getSphere(point *pointPtr, double *radiusPtr,
			Boolean *canceledPtr)
{
	do { /* Process Loop */
		/* Get center point */
		{
			/* Prompt for point coordinates */
			LbInPrintf( "\nPlease select a center point for the sphere\n");
			
			/* Get the point */
			getPoint(pointPtr, canceledPtr);
			if (*canceledPtr == true)
				break;
		}				
		/* Get radius  */
		*radiusPtr = LbInAskDouble("Enter radius ", 0, true,
			2, false, 0, canceledPtr, 0, 0, 0, 0);
	} while (false);
}

/*********************************************************************************
*
*			Name:		pointInCylinder
*
*			Summary:	See if the given point falls within the given cylinder.
*
*			Arguments:
*				point			*pointPtr			- Point in question.
*				point			*cylinderCntrPtr	- Center point of cylinder.
*				double			xAxis_Radius		- Radius in x axis.
*				double			yAxis_Radius		- Radius in y axis.
*				double			length				- Length of the cylinder.
*
*			Function return: True if point lies within sphere.
*
*********************************************************************************/
Boolean	pointInCylinder(point *pointPtr, point *cylinderCntrPtr, double xAxis_Radius,
			double yAxis_Radius, double length)
{
	Boolean	pointIsWithin = false;	/* Assume its not within object */
	double	distFromCenter;			/* Compute distance from center of sphere */
	
	do { /* Process Loop */
	
		/* If we are outside the length boundary bolt */
		if ((pointPtr->z < (cylinderCntrPtr->z - (length/2))) ||
				(pointPtr->z > (cylinderCntrPtr->z + (length/2)))) {
			
			break;
		}
		
		/* Compute distance from center of cylinder (x/y axis) */
		distFromCenter = PHGMATH_SquareRoot(
			(PHGMATH_Square(cylinderCntrPtr->x - pointPtr->x)/
			PHGMATH_Square(xAxis_Radius)) +
			(PHGMATH_Square(cylinderCntrPtr->y - pointPtr->y)/
			PHGMATH_Square(yAxis_Radius)));
			
		/* If we are less than one unit away, we are within */
		if (distFromCenter <= 1.0)
			pointIsWithin = true;
			
	} while (false);
	
	return (pointIsWithin);
}


/*********************************************************************************
*
*			Name:		pointInSphere
*
*			Summary:	See if the given point falls within the given sphere.
*
*			Arguments:
*				point			*pointPtr		- Point in question.
*				point			*objCenterPtr	- Center point of sphere
*				double			radius			- Radius.
*
*			Function return: True if point lies within sphere.
*
*********************************************************************************/
Boolean	pointInSphere(point *pointPtr, point *sphereCntrPtr, double radius)
{
	Boolean	pointIsWithin = false;	/* Assume its not within object */
	double	distFromCenter;			/* Compute distance from center of sphere */
	
	do { /* Process Loop */
	
		/* Compute distance from the center of the sphere */
		distFromCenter = PHGMATH_SquareRoot(
			PHGMATH_Square(sphereCntrPtr->x - pointPtr->x) +
			PHGMATH_Square(sphereCntrPtr->y - pointPtr->y) +
			PHGMATH_Square(sphereCntrPtr->z - pointPtr->z));

		
		/* We are in the sphere if we are within the radius */
		if (distFromCenter <= radius)
			pointIsWithin = true;
			
	} while (false);
	
	return (pointIsWithin);
}

/*********************************************************************************
*
*			Name:		pointInBox
*
*			Summary:	See if the given point falls within the given boxTy.
*
*			Arguments:
*				point			*pointPtr		- Point in question.
*				boxTy			*boxPtr			- The boxTy.
*
*			Function return: True if point lies within boxTy.
*
*********************************************************************************/
Boolean	pointInBox(point *pointPtr, boxTy *boxPtr)
{
	Boolean	pointIsWithin = false;	/* Assume its not within object */
	
	do { /* Process Loop */
	
		if (pointPtr->x < boxPtr->topLeft.x)
			break;
			
		if (pointPtr->x > boxPtr->bottomRight.x)
			break;
			
		if (pointPtr->y < boxPtr->bottomRight.y)
			break;
			
		if (pointPtr->y > boxPtr->topLeft.y)
			break;
			
		if (pointPtr->z < boxPtr->topLeft.z)
			break;
			
		if (pointPtr->z > (boxPtr->topLeft.z + boxPtr->length))
			break;
			
		pointIsWithin = true;
			
	} while (false);
	
	return (pointIsWithin);
}


/*****************************************
*		objectFilter
*
*	Arguments:
*				char			*inputString		- The input string.
*				Boolean			*isValid			- Result of filter.
*				char			*errorPrompt		- String to chastise user.
*				LbUserData		userData	- Our user data.
*
*	Purpose:	Checks to see if they want to create a BOX which is an
*				object that is not displayed in the prompt.
*
*	Returns:	True unless an error occurs.
*
******************************************/
Boolean	objectFilter(char *inputString, Boolean *isValid, char *errorPrompt,
				LbUserData userData)
{
	Boolean		okay = false;
	
	if (errorPrompt) {
		/* Do nothing, 
			but stops compiler from complaining about unused variable */
	}
	
	do	{ /* Process Loop */
	
		if ((*isValid = LbInMatchStr(inputString, "boxTy")) == true) {
			*((LbUsFourByte *)userData) = BOX;
		}
		okay = true;
	} while (false);

	return (okay);
}

/*****************************************
*		printMaterials
*
*	Arguments:
*
*	Purpose:	Prints out material names and their indexes to help
*				the user remember what material goes with what number.
*
*	Returns:	True unless an error occurs.
*
******************************************/
void printMaterials()
{
	FILE			*attenuationFile=0, *attenuationTransFile=0;
	char			errString[1024];
	char			inputBuffer[1024];
	LbUsFourByte	numTissues;
	LbUsFourByte	tissueIndex;
	unsigned long	tissueIndexUSL;
	LbUsFourByte	loopIndex;
	LbUsFourByte	attenIndex;
	LbUsFourByte	tissueIndexTranslation;
	unsigned long	tissueIndexTranslationUSL;
	LbUsFourByte	letterIndex;
	
	/* Open the file */
	if ((attenuationFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjAttenTableFilePath, "r")) == 0) {
		sprintf(errString, "Unable to open attenuation description file '%s'.",
				PhgRunTimeParams.PhgSubObjAttenTableFilePath);
		ErStFileError(errString);
		goto LEAVE_FUNC;
	}

	/* Get the number of tissue types */
	if (LbFlFGetS(inputBuffer, sizeof(inputBuffer), attenuationFile) == 0) {
		ErStGeneric("Error from LbFlFGetS while reading number of tissue indexes.");
		goto LEAVE_FUNC;
	}

	/* Set table size to supported count */
	numTissues = atol(inputBuffer);

	/* Get input a line at a time and convert it */
	for (tissueIndex = 0; tissueIndex < numTissues; tissueIndex++) {

		/* Get tissue name */
		if (LbFlFGetS(mifTissueNames[tissueIndex], sizeof(inputBuffer), attenuationFile) == 0) {
			ErStGeneric("Error from LbFlFGetS while reading tissue seperator in tissue attenuation file.");
			goto LEAVE_FUNC;
		}
		mifTissueNames[tissueIndex][strlen(mifTissueNames[tissueIndex])-1] = '\0';
		
		/* Read through and knock off positron range data if it exists */
		letterIndex = NUM_TISSUE_NAMES_LETTERS+1;
		do {
			if (letterIndex > NUM_TISSUE_NAMES_LETTERS) {
				/* First time through */
				letterIndex = 0;
			}
			else {
				letterIndex += 1;
			}
		} while ((letterIndex < strlen(mifTissueNames[tissueIndex])) && (mifTissueNames[tissueIndex][letterIndex] == ' '));

		do {
			letterIndex += 1;
		} while ((letterIndex < strlen(mifTissueNames[tissueIndex])) && (mifTissueNames[tissueIndex][letterIndex] != ' '));
		mifTissueNames[tissueIndex][letterIndex] = '\0';
		
		/* Now read past the numbers */
		for (attenIndex = 0; attenIndex < SUBOBJ_NUM_ENERGY_BINS; attenIndex++) {
			
			/* Get energy specific line */
			if (LbFlFGetS(inputBuffer, sizeof(inputBuffer), attenuationFile) == 0) {
				sprintf(errString, "Error from LbFlFGetS while reading tissue values in tissue attenuation file.\n Tissue #%ld, Tissue name '%s', energy (line#) = %ld\n string is '%s'",
					(unsigned long)tissueIndex, mifTissueNames[tissueIndex], 
					(unsigned long)attenIndex, inputBuffer);
				ErStGeneric(errString);
				goto LEAVE_FUNC;
			}
		}
	}

	/* Open the attenuation translation file */
	if ((attenuationTransFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjAttIndexTransFilePath, "r")) == 0) {
		sprintf(errString, "Unable to open attenuation translation file '%s'.",
			PhgRunTimeParams.PhgSubObjAttIndexTransFilePath);
		ErStFileError(errString);
		goto LEAVE_FUNC;
	}

	/* Allocate space for the table */
	if ((attenuationTransTbl = (LbUsFourByte *) LbMmAlloc(
			SUBOBJ_NUM_TISSUE_TRANSLATIONS * sizeof(LbUsFourByte))) == 0){
	
		goto LEAVE_FUNC;
	}
	
	/* Clear the table */
	for (tissueIndex = 0; tissueIndex < SUBOBJ_NUM_TISSUE_TRANSLATIONS; tissueIndex++)
		attenuationTransTbl[tissueIndex] = 0;
		
	/* Read in the file header */
	do {
		
		/* Get header lines */
		if (LbFlFGetS(inputBuffer, sizeof(inputBuffer), attenuationTransFile) == 0) {
			ErStFileError("Error while reading header from attenuation index translation file.");
			goto LEAVE_FUNC;
		}
	} while (inputBuffer[0] == '#');
	
	/* Read in lines and convert them (first line already read from above, its the first non-comment line */
	loopIndex = 0;
	do {
		/* Convert line to index values */
		if (sscanf(inputBuffer, " %ld %ld", 
				&tissueIndexUSL, &tissueIndexTranslationUSL) != 2){
				
				ErStGeneric("Error from sscanf on attenuation index translation value.\n");
				goto LEAVE_FUNC;
		}
		tissueIndex = (LbUsFourByte)tissueIndexUSL;
		tissueIndexTranslation = (LbUsFourByte)tissueIndexTranslationUSL;
		
		/* Verify that value is not too large */
		if (tissueIndex > SUBOBJ_NUM_TISSUE_TRANSLATIONS) {
			sprintf(errString, "Attenuation index out of range, value is %ld.\n"
				"This occurred on translation pair #%ld of the file named %s.\n"
				"The out of range value is the number in the first column.\n",
				(unsigned long)tissueIndex, (unsigned long)loopIndex+1, 
				PhgRunTimeParams.PhgSubObjAttIndexTransFilePath);
				
			ErStGeneric(errString);
			goto LEAVE_FUNC;
		}
		
		/* Verify that neither value is too large */
		if (tissueIndexTranslation > numTissues) {
			sprintf(errString, "Attenuation translation index out of range, value is %ld.\n"
				"This occurred on translation pair #%ld of the file named %s.\n"
				"The out of range value is the number in the second column.\n",
				(unsigned long)tissueIndexTranslation, (unsigned long)loopIndex+1, 
				PhgRunTimeParams.PhgSubObjAttIndexTransFilePath);
			ErStGeneric("Translation tissue index out of range in attenuation index translation table.\n");
			goto LEAVE_FUNC;
		}
		
		/* Store translation in tranlsation table, this affectively performs the translation */
		attenuationTransTbl[tissueIndex] = tissueIndexTranslation;
		
		/* Get next pair of translation values */
		if (tissueIndex < SUBOBJ_NUM_TISSUE_TRANSLATIONS-1) {
			if (LbFlFGetS(inputBuffer, sizeof(inputBuffer), attenuationTransFile) == 0) {
				ErStFileError("Error while reading attenuation index translation value.");
				goto LEAVE_FUNC;
			}
		}
		
		loopIndex++;
	} while (loopIndex <= SUBOBJ_NUM_TISSUE_TRANSLATIONS);
	
	
	/* Now print out the values */
	LbInPrintf("\nFor this file, '%s'\n"
	"\tUsing the attenuation index translation file \n\t\t%s\n",
	PhgRunTimeParams.PhgSubObjAttenTableFilePath,
	PhgRunTimeParams.PhgSubObjAttIndexTransFilePath);
	
	LEAVE_FUNC:;
	
	if (ErIsInError() == true) {
		ErHandle("Had trouble in printing material index names", false);
	}
	
	if (attenuationTransFile != 0)
		fclose(attenuationTransFile);
	
	if (attenuationFile != 0)
		fclose(attenuationFile);
		
					
}

#undef PHG_BIN
