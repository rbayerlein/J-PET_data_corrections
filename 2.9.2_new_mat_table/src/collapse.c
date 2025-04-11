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
*			Module Name:		collapse.c
*			Revision Number:	1.5
*			Date last revised:	23 July 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	18 December 1995
*
*			Module Overview:	Collapses two dimensional data into one. Prompts
*								the user to determine if they want to collapse
*								the fastest varying or slowest varying dimension.
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
#include "LbInterface.h"
#include "LbParamFile.h"
#include "LbConvert.h"
#include "LbHeader.h"



/* LOCAL CONSTANTS */

/* LOCAL TYPES */

/* LOCAL GLOBALS */
static	Boolean			canceled;				/* Global cancelation flag */
static	char			errStr[1024];			/* Storrage for error messages */
static	LbUsFourByte	corrHdrSize;			/* Size of header on data files */
static	float			*corrDataFlt = 0;		/* Data from image file 1 */
static	float			*corrResultsFlt = 0;	/* Results of correlation test */
static	double			*corrDataDbl = 0;		/* Data from image file 1 */
static	double			*corrResultsDbl = 0;	/* Results of correlation test */
static	LbFourByte		*corrDataLbf = 0;		/* Data from image file 1 */
static	LbFourByte		*corrResultsLbf = 0;	/* Results of correlation test */
static	LbUsFourByte	corrRunTimeOptions=0;	/* The runtime options specified */

#define Is_Float()			LbFgIsSet(corrRunTimeOptions, LBFlag0)	/* Did user specify float data? */
#define Is_Double()			LbFgIsSet(corrRunTimeOptions, LBFlag1)	/* Did user specify double data? */
#define Is_LbFourByte()		LbFgIsSet(corrRunTimeOptions, LBFlag2)	/* Did user specify LbFourByte data? */

/* PROTOTYPES */
Boolean	Collapse(int argc, char *argv[]);
static Boolean readImageData(char *fileName, LbUsFourByte *numBins);
static void collapseSlowestFlt(float *inputData, float *outputData, LbUsFourByte numBinsIn, LbUsFourByte numToSum);
static void collapseFastestFlt(float *inputData, float *outputData, LbUsFourByte numBinsIn, LbUsFourByte numToSum);
static void collapseSlowestDbl(double *inputData, double *outputData, LbUsFourByte numBinsIn, LbUsFourByte numToSum);
static void collapseFastestDbl(double *inputData, double *outputData, LbUsFourByte numBinsIn, LbUsFourByte numToSum);
static void collapseSlowestLbf(LbFourByte *inputData, LbFourByte *outputData, LbUsFourByte numBinsIn, LbUsFourByte numToSum);
static void collapseFastestLbf(LbFourByte *inputData, LbFourByte *outputData, LbUsFourByte numBinsIn, LbUsFourByte numToSum);
/* FUNCTIONS */

/**********************
*	Collapse
*
*	Purpose: Collapse 2 dimensional data into one.
*
*	Result:	True unless an error occurs.
***********************/
Boolean Collapse(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */
	char				errStr[1024];			/* Error string buffer */
	LbUsFourByte		numBins;				/* Number of bins in the data */
	LbUsFourByte		numToSum;				/* Number of bins to sum */
	FILE				*outputFile;			/* The output file */

	/* The following variables are for getting run time options from
		the command line 
	*/
	#define	NUM_FLAGS	3
	
	char				*knownOptions[] = {"dfi"};
	
	char				optArgs[NUM_FLAGS][LBEnMxArgLen];
	LbUsFourByte		optArgFlags = LBFlag0 + LBFlag1 + LBFlag2;
	LbUsFourByte		argIndex;
	
	do { /* Process Loop */
			
		/* Get our runtime options */
		if (!LbEnGetOptions(argc, argv, knownOptions,
				&corrRunTimeOptions, optArgs, optArgFlags, &argIndex)) {
	
			break;
		}
		
		/* Verify command line contains two arguments */
		if (argc < 3) {
			/* Tell them an image file is required */
			ErAbort("\nThis program requires an image-file name and an output file name as input.\n");
		}
		
		/* Verify the command line indicates the type of data to collapse */
		if (corrRunTimeOptions == 0) {
			ErAbort("\nThis program requires a switch to indicate the type of data being collapsed.\n"
				"Usage: collapse -{dfi} input.name output.name\n"
				"Where: -d = double precision reals -f single precision reals -i = 4 byte integer data\n");
		}
		
		/* See if they want To skip a header */
		corrHdrSize = LbInAskFourByte("Enter size of header to skip",
			1, false, false, 1, &canceled, 32768, 0, 0, 0);
	
		/* Get the number of elements to sum */
		numToSum = LbInAskFourByte("Enter number of elements to sum",
			1, false, false, 1, &canceled, 32, 0, 0, 0);

		/* Read image data */
		if (readImageData(argv[argIndex], &numBins) == false) {
			break;
		}
		
		/* Verify they are collapsing an integral amount of the data */
		if ((numBins % numToSum) != 0){
			sprintf(errStr, "You have %ld bins in the data which is not an integral multiple of the number of bins to sum (%ld)",
				(unsigned long)numBins, (unsigned long)numToSum);
			ErStGeneric(errStr);
			break;
		}
		
		/* Allocate the results buffer */
		if (Is_Float()){
			if ((corrResultsFlt = (float *) LbMmAlloc(sizeof(float) * (numBins/numToSum))) == 0){
				ErHandle("Unable to allocate memory for output data", false);
				okay = true;
				break;
			}
		}
		else if (Is_Double()){
			if ((corrResultsDbl = (double *) LbMmAlloc(sizeof(double) * (numBins/numToSum))) == 0){
				ErHandle("Unable to allocate memory for output data", false);
				okay = true;
				break;
			}
		}
		else {
			if ((corrResultsLbf = (LbFourByte *) LbMmAlloc(sizeof(LbFourByte) * (numBins/numToSum))) == 0){
				ErHandle("Unable to allocate memory for output data", false);
				okay = true;
				break;
			}
		}
			
		/* See if they want the fastest or slowest varying dimension collapsed */
		if (LbInAskYesNo("Do you want to collapse the fastest varying dimension", LBINYes)
					== LBINYes) {
					
			if (Is_Float())
				collapseFastestFlt(corrDataFlt, corrResultsFlt, numBins, numToSum);
			else if (Is_Double())
				collapseFastestDbl(corrDataDbl, corrResultsDbl, numBins, numToSum);
			else 
				collapseFastestLbf(corrDataLbf, corrResultsLbf, numBins, numToSum);
		}
		else {
			if (Is_Float())
				collapseSlowestFlt(corrDataFlt, corrResultsFlt, numBins, numToSum);
			else if (Is_Double())
				collapseSlowestDbl(corrDataDbl, corrResultsDbl, numBins, numToSum);
			else 
				collapseSlowestLbf(corrDataLbf, corrResultsLbf, numBins, numToSum);
		}
	
		/* Open the file */
		if ((outputFile = LbFlFileOpen(argv[argIndex+1], "wb")) == 0) {
			sprintf(errStr, "Unable to open file '%s'\n", argv[2]);
			ErStFileError(errStr);
			break;
		}
	
		/* Write the results */
		if (Is_Float()) {
			if (fwrite(corrResultsFlt, ((numBins/numToSum)*sizeof(float)), 1, outputFile) != 1) {
				sprintf(errStr, "Unable to write to file '%s'\n", argv[2]);
				ErStFileError(errStr);
				break;
			}
		}
		else if (Is_Double()){
			if (fwrite(corrResultsDbl, ((numBins/numToSum)*sizeof(double)), 1, outputFile) != 1) {
				sprintf(errStr, "Unable to write to file '%s'\n", argv[2]);
				ErStFileError(errStr);
				break;
			}
		}
		else {
			if (fwrite(corrResultsLbf, ((numBins/numToSum)*sizeof(LbFourByte)), 1, outputFile) != 1) {
				sprintf(errStr, "Unable to write to file '%s'\n", argv[2]);
				ErStFileError(errStr);
				break;
			}
		}
				
		okay = true;
		FAIL:;
		CANCEL:;
	} while (false);
	
	if (corrResultsFlt != 0)
		LbMmFree((void **)&corrResultsFlt);
	
	if (corrResultsDbl != 0)
		LbMmFree((void **)&corrResultsDbl);
	
	if (corrResultsLbf != 0)
		LbMmFree((void **)&corrResultsLbf);
		
	if (corrDataFlt != 0)
		LbMmFree((void **)&corrDataFlt);
		
	if (corrDataDbl != 0)
		LbMmFree((void **)&corrDataDbl);
		
	if (corrDataLbf != 0)
		LbMmFree((void **)&corrDataLbf);
		
		
	/* If error set due to cancellation, handle it, otherwise pass it on */
	if (!okay) {
		if (canceled) {
			ErHandle("User canceled collapse.", false);
			okay = true;
		}
	}
	
	/* Return the status */
	return (okay);
}

/**********************
*	collapseFastestFlt
*
*	Purpose: Collapses the fastest varying dimension of two dimensional data.
*
*	Arguments:
*		float			*inputData 	- The data to collapse
*		float			*outputData	- The data that's been collapsed
*		LbUsFourByte	numBinsIn	- The number of data bins in the input
*		LbUsFourByte	numToSum	- The number of bins to sum
*
*	Result:	None.
***********************/
static void collapseFastestFlt(float *inputData, float *outputData, LbUsFourByte numBinsIn, LbUsFourByte numToSum)
{
		
	double			sum;				/* Running sum */
	LbUsFourByte	i, j;				/* For loop index */
	
	/* Loop through the data summing the values */
	j = 0;
	sum = 0.0;
	for (i = 0; i < numBinsIn; i++){
		
		/* Sum the value */
		sum += inputData[i];
		
		/* If we are a numToSum boundary, store the sum */
		if (((i+1) % numToSum) == 0) {
			outputData[j] = sum;
			sum = 0;
			j++;
		}
	}
}

/**********************
*	collapseFastestDbl
*
*	Purpose: Collapses the fastest varying dimension of two dimensional data.
*
*	Arguments:
*		double			*inputData 	- The data to collapse
*		double			*outputData	- The data that's been collapsed
*		LbUsFourByte	numBinsIn	- The number of data bins in the input
*		LbUsFourByte	numToSum	- The number of bins to sum
*
*	Result:	None.
***********************/
static void collapseFastestDbl(double *inputData, double *outputData, LbUsFourByte numBinsIn, LbUsFourByte numToSum)
{
		
	double			sum;				/* Running sum */
	LbUsFourByte	i, j;				/* For loop index */
	
	/* Loop through the data summing the values */
	j = 0;
	sum = 0.0;
	for (i = 0; i < numBinsIn; i++){
		
		/* Sum the value */
		sum += inputData[i];
		
		/* If we are a numToSum boundary, store the sum */
		if (((i+1) % numToSum) == 0) {
			outputData[j] = sum;
			sum = 0;
			j++;
		}
	}
}

/**********************
*	collapseFastestLbf
*
*	Purpose: Collapses the fastest varying dimension of two dimensional data.
*
*	Arguments:
*		LbFourByte			*inputData 	- The data to collapse
*		LbFourByte			*outputData	- The data that's been collapsed
*		LbUsFourByte	numBinsIn	- The number of data bins in the input
*		LbUsFourByte	numToSum	- The number of bins to sum
*
*	Result:	None.
***********************/
static void collapseFastestLbf(LbFourByte *inputData, LbFourByte *outputData, LbUsFourByte numBinsIn, LbUsFourByte numToSum)
{
		
	LbFourByte		sum;					/* Running sum */
	LbUsFourByte	i, j;					/* For loop index */
	
	/* Loop through the data summing the values */
	j = 0;
	sum = 0;
	for (i = 0; i < numBinsIn; i++){
		
		/* Sum the value */
		sum += inputData[i];
		
		/* If we are a numToSum boundary, store the sum */
		if (((i+1) % numToSum) == 0) {
			outputData[j] = sum;
			sum = 0;
			j++;
		}
	}
}

/**********************
*	collapseSlowestFlt
*
*	Purpose: Collapses the slowest varying dimension of two dimensional data.
*
*	Arguments:
*		float			*inputData 	- The data to collapse
*		float			*outputData	- The data that's been collapsed
*		LbUsFourByte	numBinsIn	- The number of data bins in the input
*		LbUsFourByte	numToSum	- The number of bins to sum
*
*	Result:	None.
***********************/
static void collapseSlowestFlt(float *inputData, float *outputData, LbUsFourByte numBinsIn, LbUsFourByte numToSum)
{

	LbUsFourByte	i, j, k;	/* For loop index */
	double			*sum;
	
	/* Allocate sum buffer */
	if ((sum = (double *) LbMmAlloc((numBinsIn/numToSum)*sizeof(double))) == 0) {
		ErAbort("Unable to allocate memory for sum buffer (collapseSlowestFlt)");
	}
	
	/* Loop through the data summing the values */
	j = 0;
	for (i = 0; i < numBinsIn; i++){
		
		/* Sum the value */
		sum[j] += inputData[i];
		
		/* If we are a numToSum boundary, store the sum */
		if (((j+1) % (numBinsIn/numToSum)) == 0) {
			j = 0;
		}
		else
			j++;
	}

	/* Copy over the double precision sum */
	for (k = 0; k < numBinsIn/numToSum; k++)
		outputData[k] = sum[k];
	
	
	LbMmFree((void **)&sum);
}

/**********************
*	collapseSlowestDbl
*
*	Purpose: Collapses the slowest varying dimension of two dimensional data.
*
*	Arguments:
*		double			*inputData 	- The data to collapse
*		double			*outputData	- The data that's been collapsed
*		LbUsFourByte	numBinsIn	- The number of data bins in the input
*		LbUsFourByte	numToSum	- The number of bins to sum
*
*	Result:	None.
***********************/
static void collapseSlowestDbl(double *inputData, double *outputData, LbUsFourByte numBinsIn, LbUsFourByte numToSum)
{
		
	LbUsFourByte	i, j;				/* For loop index */
	
	/* Loop through the data summing the values */
	j = 0;
	for (i = 0; i < numBinsIn; i++){
		
		/* Sum the value */
		outputData[j] += inputData[i];
		
		/* If we are a numToSum boundary, store the sum */
		if (((j+1) % numBinsIn/numToSum) == 0) {
			j = 0;
		}
		else
			j++;
	}
}

/**********************
*	collapseSlowestLbf
*
*	Purpose: Collapses the slowest varying dimension of two dimensional data.
*
*	Arguments:
*		LbFourByte			*inputData 	- The data to collapse
*		LbFourByte			*outputData	- The data that's been collapsed
*		LbUsFourByte	numBinsIn	- The number of data bins in the input
*		LbUsFourByte	numToSum	- The number of bins to sum
*
*	Result:	None.
***********************/
static void collapseSlowestLbf(LbFourByte *inputData, LbFourByte *outputData, LbUsFourByte numBinsIn, LbUsFourByte numToSum)
{
		
	LbUsFourByte	i, j;				/* For loop index */
	
	/* Loop through the data summing the values */
	j = 0;
	for (i = 0; i < numBinsIn; i++){
		
		/* Sum the value */
		outputData[j] += inputData[i];
		
		/* If we are a numToSum boundary, store the sum */
		if (((j+1) % numBinsIn/numToSum) == 0) {
			j = 0;
		}
		else
			j++;
	}
}

/**********************
*	readImageData
*
*	Purpose: Reads in the image data, verifying that it is four byte real data.
*
*	Arguments:
*		char			*fileName 	- name of the image file
*		LbUsFourByte	*numBins	- The number of data bins
*
*	Result:	True unless an error occurs.
***********************/
static Boolean readImageData(char *fileName, LbUsFourByte *numBins)
{
	Boolean 			okay = false;
	FILE				*imageFile = 0;		/* The image file */
	LbFourByte			numRead;			/* Number of bytes read from file */
	
	do /* Process Loop */
	{
		
		/* Open the image file */
		if ((imageFile = LbFlFileOpen(fileName, "r+b")) == 0) {
			sprintf(errStr, "Unable to open image file\n'%s'.", fileName);
			ErStFileError(errStr);
			break;
		}
		
		/* Turn buffering off, it shouldn't be necessary but bugs have been found in the
			DG i/o library that make it so.
		*/
		setbuf(imageFile, 0);
	
		/* Seek to the end of the file */
		if (fseek(imageFile, 0, SEEK_END) != 0) {
			ErStFileError("Unable to seek to end of file (readImageData)");
			break;
		}
				
		/* Save the number of bins */
		if (Is_Float()){
			*numBins = (ftell(imageFile) - corrHdrSize)/sizeof(float);
		}
		else if (Is_Double()){
			*numBins = (ftell(imageFile) - corrHdrSize)/sizeof(double);
		}
		else {
			*numBins = (ftell(imageFile) - corrHdrSize)/sizeof(LbFourByte);
		}
		
		/* Check for error from ftell */
		if (*numBins < 1) {
			ErStFileError("Error determining size of file (readImageData).");
			break;
		}
			
		/* Seek to beginning of data */
		if (fseek(imageFile, corrHdrSize, SEEK_SET) != 0) {
			ErStFileError("Unable to seek to the beginning of data (readImageData).");
			break;
		}
		
		/* Allocate memory buffer for first file, using Numerical Recipes routine */
		if (Is_Float()){
			corrDataFlt = (float *)LbMmAlloc(*numBins * sizeof(float));
			if (corrDataFlt == 0) {
				goto FAIL;
			}
			
			/* Read in the data */
			if ((numRead = fread(corrDataFlt, (*numBins*sizeof(float)), 1, imageFile)) != 1) {
				sprintf(errStr, "Unable to read data for image file '%s'", fileName);
				ErStFileError(errStr);
				break;
			}
		}
		else if (Is_Double()){
			corrDataDbl = (double *)LbMmAlloc(*numBins * sizeof(double));
			if (corrDataDbl == 0) {
				goto FAIL;
			}
			
			/* Read in the data */
			if ((numRead = fread(corrDataDbl, (*numBins*sizeof(double)), 1, imageFile)) != 1) {
				sprintf(errStr, "Unable to read data for image file '%s'", fileName);
				ErStFileError(errStr);
				break;
			}
		}
		else {
			corrDataLbf = (LbFourByte *)LbMmAlloc(*numBins * sizeof(LbFourByte));
			if (corrDataLbf == 0) {
				goto FAIL;
			}
			
			/* Read in the data */
			if ((numRead = fread(corrDataLbf, (*numBins*sizeof(LbFourByte)), 1, imageFile)) != 1) {
				sprintf(errStr, "Unable to read data for image file '%s'", fileName);
				ErStFileError(errStr);
				break;
			}
		}
		
		okay = true;
		FAIL:;
	} while (false);

	/* Close image file */
	if (imageFile != 0) {
		
		fclose(imageFile);
	}
	/* Free memory if there was an error */
	if (!okay) {
		if (corrDataFlt != 0)
			LbMmFree((void **) &corrDataFlt);
		if (corrDataDbl != 0)
			LbMmFree((void **) &corrDataDbl);
		if (corrDataLbf != 0)
			LbMmFree((void **) &corrDataLbf);
	}
	
	return (okay);
}
