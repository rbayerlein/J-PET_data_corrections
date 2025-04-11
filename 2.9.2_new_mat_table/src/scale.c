/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1996-2012 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		scale.c
*			Revision Number:	1.3
*			Date last revised:	14 December 2012
*			Programmer:			Steven Vannoy
*			Date Originated:	5 December, 1996
*
*			Module Overview:	Scales file 2 to file 1 by multiplying each
*								value in file 2 by sum(file1)/sum(file2).
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

#include <stdio.h>
#include <string.h>
#include <stddef.h>

#include "SystemDependent.h"

#include "LbTypes.h"
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
static	Boolean			canceled;		/* Global cancelation flag */
static	char			errStr[1024];	/* Storrage for error messages */
static	LbUsFourByte	scaleHdrSize1;	/* Size of header on data file 1 */
static	LbUsFourByte	scaleHdrSize2;	/* Size of header on data file 2 */

/* PROTOTYPES */
Boolean	Scale(int argc, char *argv[]);
static Boolean readImageData(char *fileName1, char *fileName2,
		float **imageData1, float **imageData2, LbFourByte *numBins);

/* FUNCTIONS */

/**********************
*	Scale
*
*	Purpose: Scale file 2 to file 1.
*
*	Result:	True unless an error occurs.
***********************/
Boolean Scale(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */


	char				errStr[1024];			/* Error string buffer */
	char				prompt[1024];			/* String for building prompt */
	float				*data1 = 0;				/* Data from image file 1 */
	float				*data2 = 0;				/* Data from image file 2 */
	double				sum1 = 0;				/* Sum of file 1 */
	double				sum2 = 0;				/* Sum of file 2 */
	double				scaleFactor = 0;		/* Scale factor */
	LbFourByte			i;						/* For loop index */
	LbFourByte			numBins;				/* Number of bins in the data */
	FILE				*imageFile = 0;		/* The image file */

	do { /* Process Loop */
		
		/* Verify command line contains two arguments */
		if (argc < 3) {
			/* Tell them 2 image files are required */
			ErAbort("\nThis program requires two image-file names as input.\n");
		}
		
		/* See if they want To skip a header */
		sprintf(prompt, "Enter size of header for '%s'", argv[1]);
		scaleHdrSize1 = LbInAskFourByte(prompt,
			1, false, false, 1, &canceled, 32768, 0, 0, 0);

		sprintf(prompt, "Enter size of header for '%s'", argv[2]);
		scaleHdrSize2 = LbInAskFourByte(prompt,
			1, false, false, 1, &canceled, scaleHdrSize1, 0, 0, 0);
	
		/* Read image data 1 */
		if (readImageData(argv[1], argv[2], &data1, &data2, &numBins) == false) {
			break;
		}
	
		/* Compute the sums */
		for (i = 0; i < numBins; i++){
			sum1 += data1[i];
			sum2 += data2[i];
		}
		
		scaleFactor = sum1/sum2;
		
		/* Print values */
		LbInPrintf("Sum for file 1 is %3.3e,  Sum for file 2 is %3.3e, scale factor is %3.3e\n",
			sum1, sum2, scaleFactor);
			
		/* Give user a chance to abort based on values */
		if (LbInAskYesNo("Do you want still want to scale the data", LBINYes)
					== LBINNo) {

			ErStCancel("They didn't like the scale parameters");
			canceled = true;
			goto CANCEL;
		}
		
		/* Scale the data */
		for (i = 0; i < numBins; i++){
			data2[i] *= scaleFactor;
		}
		
		/* Open the image file */
		if ((imageFile = LbFlFileOpen(argv[2], "r+b")) == 0) {
			sprintf(errStr, "Unable to open image file\n'%s'.", argv[2]);
			ErStFileError(errStr);
			break;
		}
		
		setbuf(imageFile, 0);
		
		/* Seek to beginning of data */
		if (fseek(imageFile, scaleHdrSize2, SEEK_SET) != 0) {
			ErStFileError("Unable to seek to the beginning of data for file 2 (readImageData).");
			break;
		}
		
		/* Write in the data */
		if ((fwrite(data2, (numBins*sizeof(float)), 1, imageFile)) != 1) {
			sprintf(errStr, "Unable to write data for image file '%s'", argv[2]);
			ErStFileError(errStr);
			break;
		}

		fclose(imageFile);
		
		okay = true;
		FAIL:;
		CANCEL:;
	} while (false);
	
		
	if (data1 != 0)
		LbMmFree((void **)&data1);
		
	if (data2 != 0)
		LbMmFree((void **)&data1);

	/* If error set due to cancellation, handle it, otherwise pass it on */
	if (!okay) {
		if (canceled) {
			ErHandle("User canceled scale.", false);
			okay = true;
		}
	}
	
	/* Return the status */
	return (okay);
}

/**********************
*	readImageData
*
*	Purpose: Reads in the image data, verifying that it is four byte real data.
*
*	Arguments:
*		char			*fileName1 	- name of the image file
*		char			*fileName2 	- name of the image file
*		float			**imageData2	- The data
*		float			**imageData2	- The data
*		LbFourByte		*numBins	- The number of data bins
*
*	Result:	True unless an error occurs.
***********************/
static Boolean readImageData(char *fileName1, char *fileName2,
					float **imageData1, float **imageData2, LbFourByte *numBins)
{
	Boolean 			okay = false;
	FILE				*imageFile = 0;		/* The image file */
	LbFourByte			numRead;			/* Number of bytes read from file */
	LbFourByte			numBins2;			/* Number of file 2 bins */
	
	do /* Process Loop */
	{
		/* Clear memory variable right away */
		*imageData1 = 0;
		*imageData2 = 0;
		
		/* Open the image file */
		if ((imageFile = LbFlFileOpen(fileName1, "r+b")) == 0) {
			sprintf(errStr, "Unable to open image file\n'%s'.", fileName1);
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
		*numBins = (ftell(imageFile) - scaleHdrSize1)/sizeof(float);
		
		/* Check for error from ftell */
		if (*numBins < 1) {
			sprintf(errStr, "Error determining size of '%s', possibly incorrect header size (readImageData).",
				fileName1);
			ErStFileError(errStr);
			break;
		}
			
		/* Seek to beginning of data */
		if (fseek(imageFile, scaleHdrSize1, SEEK_SET) != 0) {
			ErStFileError("Unable to seek to the beginning of data (readImageData).");
			break;
		}
		
		/* Allocate memory buffer for first file */
		*imageData1 = (float *)LbMmAlloc(*numBins * sizeof(float));
		if (*imageData1 == 0) {
			ErStGeneric("Failed to allocate memory for the image data (readImageData)");
			goto FAIL;
		}
		
		/* Read in the data */
		if ((numRead = fread(*imageData1, (*numBins*sizeof(float)), 1, imageFile)) != 1) {
			sprintf(errStr, "Unable to read data for image file '%s'", fileName1);
			ErStFileError(errStr);
			break;
		}
		
		fclose(imageFile);
		
		
		/* Open the image file */
		if ((imageFile = LbFlFileOpen(fileName2, "r+b")) == 0) {
			sprintf(errStr, "Unable to open image file\n'%s'.", fileName2);
			ErStFileError(errStr);
			break;
		}
		
		setbuf(imageFile, 0);
	
		/* Seek to the end of the file */
		if (fseek(imageFile, 0, SEEK_END) != 0) {
			ErStFileError("Unable to seek to end of file 2 (readImageData)");
			break;
		}
				
		/* Save the number of bins */
		numBins2 = (ftell(imageFile) - scaleHdrSize2)/sizeof(float);
		
		if (*numBins != numBins2) {
			ErStGeneric("These files are not the same size!");
			break;
		}

		/* Seek to beginning of data */
		if (fseek(imageFile, scaleHdrSize2, SEEK_SET) != 0) {
			ErStFileError("Unable to seek to the beginning of data for file 2 (readImageData).");
			break;
		}
		
		/* Allocate memory buffer for second file */
		*imageData2 = (float *)LbMmAlloc(*numBins * sizeof(float));
		if (*imageData1 == 0) {
			ErStGeneric("Failed to allocate memory for the image data (readImageData)");
			goto FAIL;
		}
		
		/* Read in the data */
		if ((numRead = fread(*imageData2, (*numBins*sizeof(float)), 1, imageFile)) != 1) {
			sprintf(errStr, "Unable to read data for image file '%s'", fileName2);
			ErStFileError(errStr);
			break;
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
		if (*imageData1 != 0)
			LbMmFree((void **)imageData1);
		if (*imageData2 != 0)
			LbMmFree((void **)imageData2);
	}
	
	return (okay);
}

