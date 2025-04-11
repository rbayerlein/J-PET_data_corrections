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
*     Module Name: 			convertCoh.c
*     Revision Number:		1.4
*     Date last revised:	23 July 2013
*     Programmer:			Steven Vannoy
*     Date Originated:		Thursday April 4, 1997
*
*     Module Overview:	Converts text based coherent scatter angular distribution
*						files to binary format.
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
#define CONVERTCOH


#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"
#include "LbHeader.h"


/* Local Constants */
#define NUM_ANGLES		100
#define NUM_MATERIALS	22
#define	NUM_COH_ENERGIES	172
/* LOCAL TYPES */
typedef char	labelTy[1024];

typedef	struct {

	LbUsFourByte	energy;
	LbUsFourByte	materialIndex;
	double			angleProbabilities[NUM_ANGLES];
	} materialEntryTy;
	
/* Local Globals */
static	char			convertCohErrString[1024];		/* For building error messages */
static	materialEntryTy	convertCohMaterials[NUM_MATERIALS][NUM_COH_ENERGIES];

static labelTy	convertCohMaterialLabels[] = {
				"air",
				"water",
				"blood",
				"bone",
				"brain",
				"heart",
				"lung",
				"muscle",
				"lead",
				"NaI",
				"BGO",
				"iron",
				"graphite",
				"tin",
				"GI_tract",
				"con_tissue",
				"copper",
				"perfect_absorber",
				"LSO",
				"GSO",
				"aluminum",
				"tungsten"
				};

/* Prototypes */
Boolean	convertCoh(int argc, char *argv[]);
LbFourByte	convertCohGtMaterialIndex(char *theMaterialStr);
void		convertCohTest(FILE *outputFile);

/**********************
*	convertCoh
*
*	Purpose:	Execute the program.
*
*	Result:	None.
***********************/

Boolean	convertCoh(int argc, char *argv[])
{
	Boolean			okay = false;			/* Process flag */
	
	
	char				inputBuffer[1024];
	char				kevStr[1024];
	char				materialStr[1024];
	char				theMaterialStr[1024];
	FILE				*inputFile=0;
	FILE				*outputFile=0;
	double				angleProbability;
	LbFourByte			lineCount;
	LbFourByte			energy;
	long				energyL;
	LbFourByte			materialIndex;
	LbFourByte			numItems;
	LbFourByte			totalWritten = 0;
	LbFourByte			numWritten;
	LbFourByte			nameIndex;
	char				inputName[32];
	char				outputName[32];
	
	/* Remove compiler warnings */
	if (argc) {};
	if (argv) {};
	
	do	 {	/* Process Loop */

	for (nameIndex = 0; nameIndex < NUM_MATERIALS; nameIndex++) {
		sprintf(inputName, "%s.ad", convertCohMaterialLabels[nameIndex]);
		sprintf(outputName, "%s.ad.bin", convertCohMaterialLabels[nameIndex]);
		
		/* Open the image file */
		if ((inputFile = LbFlFileOpen(inputName, "r")) == 0) {
			sprintf(convertCohErrString, "Unable to open image file\n'%s' (convertCoh).", inputName);
			ErStFileError(convertCohErrString);
			goto FAIL;
		}
		
		/* Open the image file */
		if ((outputFile = LbFlFileOpen(outputName, "w+b")) == 0) {
			sprintf(convertCohErrString, "Unable to create/open image file\n'%s' (convertCoh).", outputName);
			ErStFileError(convertCohErrString);
			goto FAIL;
		}		

		lineCount = 0;
		
		while (LbFlFGetS(inputBuffer, 1024, inputFile) != NULL) {
			
			/* If it's the first line, then convert text to header info */
			if ((lineCount == 0) || ((lineCount) % (NUM_ANGLES+1) == 0)) {
				
				/* Convert values */
				if ((numItems = sscanf(inputBuffer, "%ld %s %s = %s",
						&energyL, kevStr, materialStr, theMaterialStr)) != 4){
						
						sprintf(convertCohErrString,"Error from sscanf on header info, string = '%s'\n"
							" expected something like 1 keV, material = air. (convertCoh)",
							inputBuffer);
						ErStGeneric(convertCohErrString);
						goto FAIL;
				}
				energy = (LbFourByte)energyL;
				
				/* Get the material index */
				if ((materialIndex = convertCohGtMaterialIndex(theMaterialStr)) == -1) {
					goto FAIL;
				}
				
				/* Write energy out */
				if ((numWritten = fwrite(&energy, 1, sizeof(energy), outputFile)) != sizeof(energy)) {
					ErStFileError("Error writing energy value (convertCoh)");
					goto FAIL;
				}
				
				totalWritten += numWritten;
				
				/* Write material out */
				if ((numWritten = fwrite(&materialIndex, 1, sizeof(materialIndex), outputFile)) != sizeof(materialIndex)) {
					ErStFileError("Error writing material index value (convertCoh)");
					goto FAIL;
				}
				totalWritten += numWritten;
			}
			else {

				/* Convert values */
				if ((sscanf(inputBuffer, "%lf",
						&angleProbability)) != 1){
						
						sprintf(convertCohErrString,"Error from sscanf on angle probability, string = '%s'\n"
							" expected something like .00874. (convertCoh)",
							inputBuffer);
						ErStGeneric(convertCohErrString);
						goto FAIL;
				}
				
				/* Write material out */
				if ((numWritten = fwrite(&angleProbability, 1, sizeof(angleProbability), outputFile)) != sizeof(angleProbability)) {
					ErStFileError("Error writing angle probability value (convertCoh)");
					goto FAIL;
				}
				totalWritten += numWritten;
			}
			lineCount++;
		}
		
		fclose(inputFile);
		fclose(outputFile);
	}
	okay = true;
	FAIL:;
	} while (0);
	
	#ifdef DO_TESTING
		if (okay == true) {
			fclose(outputFile);

			if ((outputFile = LbFlFileOpen(outputName, "rb")) == 0) {
				sprintf(convertCohErrString, "Unable to create/open image file\n'%s' (convertCoh).", outputName);
				ErStFileError(convertCohErrString);
				return(false);
			}		
			
			convertCohTest(outputFile);
		}
	#endif
	
	if (inputFile != 0)
		fclose(inputFile);
	if (outputFile != 0)
		fclose(outputFile);
	
	return(okay);
}

/**********************
*	convertCohGtMaterialIndex
*
*	Purpose:	Convert text string to material index.
*	Arguments:
*		char	* theMaterialStr - text label of material string
*	Result:	The index of the material.
***********************/

LbFourByte	convertCohGtMaterialIndex(char *theMaterialStr)
{
	LbFourByte	index;
	LbFourByte	materialIndex = -1;
	
	/* Search table for given label */
	for (index = 0; index < NUM_MATERIALS; index++){
	
		/* See if it matches */
		if (strcmp(theMaterialStr, convertCohMaterialLabels[index]) == 0) {
			materialIndex = index;
			break;
		}
	}

	/* If not found, register an error */
	if (materialIndex == -1){
		sprintf(convertCohErrString, "Unable to find '%s' in material table (convertCohGtMaterialIndex)",
			theMaterialStr);
			
		ErStGeneric(convertCohErrString);
	}
	return(materialIndex);
}


/**********************
*	convertCohTest
*
*	Purpose:	Verify binary file is valid.
*	Arguments:
*		FILE	* outoutFile - the file just created
*	Result:	None.
***********************/

void	convertCohTest(FILE *outputFile)
{
	LbUsFourByte	expectedSize = sizeof(convertCohMaterials);
	LbUsFourByte	fileSize = 0;
	
	/* First verify that file is of the correct size */
	if (fseek(outputFile, 0, SEEK_END) != 0){
		ErStFileError("Unable to seek to end binary file (convertCohTest).");
		goto FAIL;
	}
	
	fileSize = ftell(outputFile);
	
	if (fileSize != expectedSize) {
		LbInPrintf("Expected size = %ld, actual size = %ld\n",
			(unsigned long)expectedSize, (unsigned long)fileSize);
	
		goto FAIL;
	}
	
	/* Now read in the data */
	if (fseek(outputFile, 0, SEEK_SET) != 0){
		ErStFileError("Unable to seek to beginning binary file (convertCohTest).");
		goto FAIL;
	}

	if (fread(convertCohMaterials, expectedSize, 1, outputFile) != 1) {
		ErStFileError("Unable to read data: (convertCohTest).");
		goto FAIL;
	}
	
	/* This is basically here so you have a place to break with the debugger */
	if (convertCohMaterials[0][0].energy != 1) {
		LbInPrintf("Exected 1 keV as first entery, got %d\n",
			convertCohMaterials[0][0].energy);
	}
	
	FAIL:;
}
#undef CONVERTCOH
