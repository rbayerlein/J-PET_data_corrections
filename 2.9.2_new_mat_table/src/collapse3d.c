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
*			Module Name:		collapse3d.c
*			Revision Number:	1.3
*			Date last revised:	4 June 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	18 December, 1995
*
*			Module Overview:	Collapses three dimensional data into two. The
*								assumption is that we have multiple "slices" of
*								two dimensional data and we want a single two
*								dimensional array containing the sum of the
*								slices.
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
static	LbUsFourByte	corrHdrSize;			/* Size of header on data files */
static	float			*corrDataFlt = 0;		/* Data from image file 1 */
static	float			*corrResultsFlt = 0;	/* Results of correlation test */
static	double			*corrDataDbl = 0;		/* Data from image file 1 */
static	double			*corrResultsDbl = 0;	/* Results of correlation test */
static	LbFourByte		*corrDataLbf = 0;		/* Data from image file 1 */
static	LbFourByte		*corrResultsLbf = 0;	/* Results of correlation test */
static	LbUsFourByte	corrRunTimeOptions=LBFlag0;	/* The runtime options specified, default to float */

#define Is_Float()			LbFgIsSet(corrRunTimeOptions, LBFlag0)	/* Did user specify float data? */
#define Is_Double()			LbFgIsSet(corrRunTimeOptions, LBFlag1)	/* Did user specify double data? */
#define Is_LbFourByte()		LbFgIsSet(corrRunTimeOptions, LBFlag2)	/* Did user specify LbFourByte data? */

/* PROTOTYPES */
Boolean	Collapse3d(int argc, char *argv[]);

/* FUNCTIONS */

/**********************
*	Collapse3d
*
*	Purpose: Collapse three dimensional data into two.
*
*	Result:	True unless an error occurs.
***********************/
Boolean Collapse3d(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */
	char				errStr[1024];			/* Error string buffer */
	LbUsFourByte		numRows;				/* Number of bins to sum */
	LbUsFourByte		numColumns;				/* Number of bins to sum */
	LbUsFourByte		numSlices;				/* Number of bins to sum */
	LbUsFourByte		curSlice;				/* Number of bins to sum */
	LbUsFourByte		curRow;					/* Number of bins to sum */
	LbUsFourByte		curCol;					/* Number of bins to sum */
	LbUsFourByte		numRead;				/* Number of bytes read from the file */
	double				curSum;					/* Running sum */
	FILE				*outputFile;			/* The output file */
	FILE				*imageFile;

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
		
		/* Open the image file */
		if ((imageFile = LbFlFileOpen(argv[argIndex], "r+b")) == 0) {
			sprintf(errStr, "Unable to open image file\n'%s'.", argv[argIndex]);
			ErStFileError(errStr);
			break;
		}
		
		/* See if they want To skip a header */
		corrHdrSize = LbInAskFourByte("Enter size of header to skip",
			1, false, false, 1, &canceled, 32768, 0, 0, 0);
	
		/* Get the number of rows */
		numRows = LbInAskFourByte("Enter number of rows",
			1, false, false, 1, &canceled, 32, 0, 0, 0);
	
		/* Get the number of columns */
		numColumns = LbInAskFourByte("Enter number of columns",
			1, false, false, 1, &canceled, 32, 0, 0, 0);
	
		/* Get the number of slices */
		numSlices = LbInAskFourByte("Enter number of slices",
			1, false, false, 1, &canceled, 32, 0, 0, 0);
		
		/* Turn buffering off, it shouldn't be necessary but bugs have been found in the
			DG i/o library that make it so.
		*/
		setbuf(imageFile, 0);
					
		/* Seek past the header */
		if (corrHdrSize != 0){
			/* Seek to the end of the file */
			if (fseek(imageFile, corrHdrSize, SEEK_SET) != 0) {
				ErStFileError("Unable to seek to beginning of data (Collapse3d)");
				break;
			}
		}	
		/* Allocate the data and results buffers */
		if (Is_Float()){
			if ((corrResultsFlt = (float *) LbMmAlloc(sizeof(float) * (numColumns*numRows))) == 0){
				ErHandle("Unable to allocate memory for output data", false);
				okay = true;
				break;
			}
			
			if ((corrDataFlt = (float *) LbMmAlloc(sizeof(float) * (numColumns*numRows))) == 0){
				ErHandle("Unable to allocate memory for input data", false);
				okay = true;
				break;
			}
		}
		else if (Is_Double()){
			if ((corrResultsDbl = (double *) LbMmAlloc(sizeof(double) * (numColumns*numRows))) == 0){
				ErHandle("Unable to allocate memory for output data", false);
				okay = true;
				break;
			}
			if ((corrDataDbl = (double *) LbMmAlloc(sizeof(double) * (numColumns*numRows))) == 0){
				ErHandle("Unable to allocate memory for output data", false);
				okay = true;
				break;
			}
		}
		else {
			if ((corrResultsLbf = (LbFourByte *) LbMmAlloc(sizeof(LbFourByte) * (numColumns*numRows))) == 0){
				ErHandle("Unable to allocate memory for output data", false);
				okay = true;
				break;
			}
			if ((corrDataLbf = (LbFourByte *) LbMmAlloc(sizeof(LbFourByte) * (numColumns*numRows))) == 0){
				ErHandle("Unable to allocate memory for output data", false);
				okay = true;
				break;
			}
		}
		
		/* Read in the first slice, straight to results buffer */
		if (Is_Float()) {
			if ((numRead = fread(corrResultsFlt, ((numColumns*numRows)*sizeof(float)), 1, imageFile)) != 1) {
				sprintf(errStr, "Unable to read data for image file '%s'", argv[argIndex]);
				ErStFileError(errStr);
				break;
			}
		}
		else if (Is_Double()) {
			if ((numRead = fread(corrResultsDbl, ((numColumns*numRows)*sizeof(double)), 1, imageFile)) != 1) {
				sprintf(errStr, "Unable to read data for image file '%s'",  argv[argIndex]);
				ErStFileError(errStr);
				break;
			}
		}
		else  {
			if ((numRead = fread(corrResultsLbf, ((numColumns*numRows)*sizeof(LbUsFourByte)), 1, imageFile)) != 1) {
				sprintf(errStr, "Unable to read data for image file '%s'",  argv[argIndex]);
				ErStFileError(errStr);
				break;
			}
		}
		
		/* Now for each remaining slice read in the data and sum them together */
		for (curSlice = 1; curSlice < numSlices; curSlice++){
		
			/* Read in the next slice */
			if (Is_Float()) {
				if ((numRead = fread(corrDataFlt, ((numColumns*numRows)*sizeof(float)), 1, imageFile)) != 1) {
					sprintf(errStr, "Unable to read data for image file '%s'",  argv[argIndex]);
					ErStFileError(errStr);
					goto FAIL;
				}
			}
			else if (Is_Double()) {
				if ((numRead = fread(corrDataDbl, ((numColumns*numRows)*sizeof(double)), 1, imageFile)) != 1) {
					sprintf(errStr, "Unable to read data for image file '%s'",  argv[argIndex]);
					ErStFileError(errStr);
					goto FAIL;
				}
			}
			else  {
				if ((numRead = fread(corrDataLbf, ((numColumns*numRows)*sizeof(LbUsFourByte)), 1, imageFile)) != 1) {
					sprintf(errStr, "Unable to read data for image file '%s'",  argv[argIndex]);
					ErStFileError(errStr);
					goto FAIL;
				}
			}
			
			/* Loop through and sum them together */
			for (curRow = 0; curRow < numRows; curRow++){
				
				curSum  = 0.0;
				
				/* Sum the row */
				for (curCol = 0; curCol < numColumns; curCol++) {

					/* Read in the next slice */
					if (Is_Float()) {
						corrResultsFlt[(curRow*numColumns)+curCol] += corrDataFlt[(curRow*numColumns)+curCol];
					}
					else if (Is_Double()) {
						corrResultsDbl[(curRow*numColumns)+curCol] += corrDataDbl[(curRow*numColumns)+curCol];
					}
					else  {
						corrResultsLbf[(curRow*numColumns)+curCol] += corrDataLbf[(curRow*numColumns)+curCol];
					}
				
				}

			}
		}
		
		/* Open the file */
		if ((outputFile = LbFlFileOpen(argv[argIndex+1], "wb")) == 0) {
			sprintf(errStr, "Unable to open file '%s'\n", argv[2]);
			ErStFileError(errStr);
			break;
		}
	
		/* Write the results */
		if (Is_Float()) {
			if (fwrite(corrResultsFlt, ((numColumns*numRows)*sizeof(float)), 1, outputFile) != 1) {
				sprintf(errStr, "Unable to write to file '%s'\n", argv[2]);
				ErStFileError(errStr);
				break;
			}
		}
		else if (Is_Double()){
			if (fwrite(corrResultsDbl, ((numColumns*numRows)*sizeof(double)), 1, outputFile) != 1) {
				sprintf(errStr, "Unable to write to file '%s'\n", argv[2]);
				ErStFileError(errStr);
				break;
			}
		}
		else {
			if (fwrite(corrResultsLbf, ((numColumns*numRows)*sizeof(LbFourByte)), 1, outputFile) != 1) {
				sprintf(errStr, "Unable to write to file '%s'\n", argv[2]);
				ErStFileError(errStr);
				break;
			}
		}
				
		okay = true;
		FAIL:;
		CANCEL:;
	} while (false);
	
	if (corrResultsFlt != 0) {
		LbMmFree((void **)&corrResultsFlt);
	}
	
	if (corrResultsDbl != 0) {
		LbMmFree((void **)&corrResultsDbl);
	}
	
	if (corrResultsLbf != 0) {
		LbMmFree((void **)&corrResultsLbf);
	}
		
	if (corrDataFlt != 0) {
		LbMmFree((void **)&corrDataFlt);
	}
		
	if (corrDataDbl != 0) {
		LbMmFree((void **)&corrDataDbl);
	}
		
	if (corrDataLbf != 0) {
		LbMmFree((void **)&corrDataLbf);
	}
	
	if (imageFile != 0)
		fclose(imageFile);
		
	if (outputFile != 0)
		fclose(outputFile);
		
	/* If error set due to cancellation, handle it, otherwise pass it on */
	if (!okay) {
		if (canceled) {
			ErHandle("User canceled collapse3d.", false);
			okay = true;
		}
	}
	
	/* Return the status */
	return (okay);
}
