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
*			Module Name:		extract.c
*			Revision Number:	1.3
*			Date last revised:	23 July 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	18 December, 1995
*
*			Module Overview:	Extracts a given portion of a file and stores it
*								in the output file.
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
#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbMemory.h"
#include "LbInterface.h"
#include "LbParamFile.h"
#include "LbConvert.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhgMath.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhoHFile.h"
#include "PhgHdr.h"
#include "phg.h"
#include "PhgBin.h"


/* LOCAL CONSTANTS */

/* LOCAL TYPES */

/* LOCAL GLOBALS */
static	Boolean			canceled;		/* Global cancelation flag */
static	char			errStr[1024];	/* Storrage for error messages */
static	LbUsFourByte	corrHdrSize;	/* Size of header on data files */

/* PROTOTYPES */
Boolean	Extract(int argc, char *argv[]);
static Boolean readImageData(char *fileName, void **imageData, LbUsFourByte numToExtract, LbUsFourByte firstToKeep, LbUsFourByte elemSize);
/* FUNCTIONS */

/**********************
*	Extract
*
*	Purpose: Extracts variable size elements from a file.
*
*	Result:	True unless an error occurs.
***********************/
Boolean Extract(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */
	char				errStr[1024];			/* Error string buffer */
	void				*data1 = 0;				/* Data from image file 1 */
	LbUsFourByte		firstToKeep;			/* First element to keep */
	LbUsFourByte		numToExtract;				/* Number of elements to keep */
	LbUsFourByte		elemSize;				/* Size of elements */
	FILE				*outputFile;			/* The output file */
	
	do { /* Process Loop */
		
		
		/* Verify command line contains two arguments */
		if (argc < 3) {
			/* Tell them an image file is required */
			ErAbort("\nThis program requires an image-file name and an output file name as input.\n");
		}
		
		
		/* Open the output file, input file get's opened by extraction routine */
		if ((outputFile = LbFlFileOpen(argv[2], "wb")) == 0) {
			sprintf(errStr, "Unable to open file '%s'\n", argv[2]);
			ErStFileError(errStr);
			break;
		}

		/* See if they want To skip a header */
		corrHdrSize = LbInAskFourByte("Enter size of header to skip",
			1, false, false, 1, &canceled, 32768, 0, 0, 0);
	
		/* Enter the size of each element */
		elemSize = LbInAskFourByte("Enter the size of each element to extract (e.g., floats = 4)",
			1, false, false, 1, &canceled, 4, 0, 0, 0);
	
		/* Get the first element to extract */
		firstToKeep = LbInAskFourByte("Enter element number of first element to extract (counting from 0)",
			1, false, false, 1, &canceled, 32, 0, 0, 0);
	
		/* Get the last element to extract */
		numToExtract = LbInAskFourByte("Enter number of elements to extract",
			1, false, false, 1, &canceled, 32, 0, 0, 0);

		/* Read image data 1 */
		if (readImageData(argv[1], &data1, numToExtract, firstToKeep, elemSize) == false) {
			break;
		}
	
		if (fwrite(data1, (numToExtract*elemSize), 1, outputFile) != 1) {
			sprintf(errStr, "Unable to write to file '%s'\n", argv[2]);
			ErStFileError(errStr);
			break;
		}
				
		okay = true;
		FAIL:;
		CANCEL:;
	} while (false);
		
	if (data1 != 0)
		LbMmFree((void **)&data1);
		
		
	/* If error set due to cancellation, handle it, otherwise pass it on */
	if (!okay) {
		if (canceled) {
			ErHandle("User canceled Extract.", false);
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
*		char			*fileName 		- name of the image file
*		float			**imageData		- The data
*		LbUsFourByte	numToExtract	- The number of data elements to extract
*		LbUsFourByte	firstToKeep		- The number of first element to extract
*		LbUsFourByte	elemSize		- The size of the elements
*
*	Result:	True unless an error occurs.
***********************/
static Boolean readImageData(char *fileName, void **imageData, LbUsFourByte numToExtract, LbUsFourByte firstToKeep, LbUsFourByte elemSize)
{
	Boolean 			okay = false;
	FILE				*imageFile = 0;		/* The image file */
	LbFourByte			numRead;			/* Number of bytes read from file */
	
	do /* Process Loop */
	{
		/* Clear memory variable right away */
		*imageData = 0;
		
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
						
		/* Verify the file is as large as the number of elements they asked for */
		if ((numToExtract) > ((ftell(imageFile) - corrHdrSize - (firstToKeep *elemSize))/elemSize)) {
			sprintf(errStr, "File does not have %ld elements starting at %ld, only %ld!\n", (unsigned long)numToExtract,
			(long)ftell(imageFile)/elemSize, (long)((ftell(imageFile) - corrHdrSize- (firstToKeep *elemSize))/elemSize));
			ErStGeneric(errStr);
			break;
		}
			
		/* Seek to beginning of data */
		if (fseek(imageFile, ((firstToKeep * elemSize)+ corrHdrSize), SEEK_SET) != 0) {
			ErStFileError("Unable to seek to the beginning of data (readImageData).");
			break;
		}
		
		/* Allocate memory buffer for data */
		*imageData = (float *)LbMmAlloc(numToExtract * elemSize);
		if (*imageData == 0) {
			break;
		}
		
		/* Read in the data */
		if ((numRead = fread(*imageData, (numToExtract*elemSize), 1, imageFile)) != 1) {
			sprintf(errStr, "Unable to read data for image file '%s'", fileName);
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
		if (*imageData != 0)
			LbMmFree((void **)imageData);
	}
	
	return (okay);
}
