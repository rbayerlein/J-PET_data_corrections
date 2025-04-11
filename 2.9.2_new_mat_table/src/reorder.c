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
*     Module Name: 			reorder.c
*     Revision Number:		1.3
*     Date last revised:	4 June 2013
*     Programmer:			Steven Vannoy
*     Date Originated:		Monday, March 31, 1997
*
*     Module Overview:	Reorder 3D data, if data is a[i][j][k], it will end up
*						b[k][j][i].
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
#define REORDER


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

/* Local Globals */
static	char			reorderErrString[1024];		/* For building error messages */
static	Boolean			reorderCanceled = false;	/* Cancelation flag */

/* The variable for run-time options */
static	LbUsFourByte	reorderRunTimeOptions=0;	/* The runtime options specified */

/* Flag macros for accessing run-time options */
#define Is_Integer()			LbFgIsSet(reorderRunTimeOptions, LBFlag0)	/* Did user specify an integer image? */
#define Is_Float()			LbFgIsSet(reorderRunTimeOptions, LBFlag1)	/* Did user specify a float image? */
#define Is_Double()	LbFgIsSet(reorderRunTimeOptions, LBFlag2)			/* Did user specify a double image? */

/* Prototypes */
Boolean	reorder(int argc, char *argv[]);


/**********************
*	reorder
*
*	Purpose:	Execute the program.
*
*	Result:	None.
***********************/

Boolean	reorder(int argc, char *argv[])
{
	Boolean			okay = false;			/* Process flag */
	
	/* The following variables are for getting run time options from
		the command line.  Flag 'd' is not used, but for illustration purpose,  
	*/
	#define	NUM_FLAGS	3
	char				*knownOptions[] = {"ifd"};
	char				optArgs[NUM_FLAGS][LBEnMxArgLen];
	LbUsFourByte		optArgFlags = LBFlag0 + LBFlag1 + LBFlag2;
	LbUsFourByte		argIndex;
	
	char				*newOrderStr[12] = {"slices", "rows", "columns"};
	
	char				*inputName = argv[argc-2];
	char				*outputName = argv[argc-1];
	char				resultStr[12];
	
	FILE				*inputFile=0;
	FILE				*outputFile=0;
	LbUsFourByte		numRead;
	LbUsFourByte		numWritten;
	LbUsFourByte		fileSize;
	LbUsFourByte		dataSize;
	LbUsFourByte		hdrSize;
	LbUsFourByte		numRows;
	LbUsFourByte		numColumns;
	LbUsFourByte		numSlices;
	LbUsFourByte		newRows;
	LbUsFourByte		newColumns;
	LbUsFourByte		newSlices;
	LbUsFourByte		newRowsOrder;
	LbUsFourByte		newColumnsOrder;
	LbUsFourByte		newSlicesOrder;
	LbUsFourByte		*newIndexOrder[3]= {0,0,0};
	LbUsFourByte		i,j,k, ii, jj, kk;
	LbUsFourByte		*intData=0;
	float				*floatData=0;
	double				*doubleData=0;
	LbUsFourByte		*outputInts=0;
	float				*outputFloats=0;
	double				*outputDoubles=0;
	
	do	 {	/* Process Loop */
		
		/* Get our runtime options */
		if (!LbEnGetOptions(argc, argv, knownOptions,
				&reorderRunTimeOptions, optArgs, optArgFlags, &argIndex)) {
	
			goto FAIL;
		}
		
		/* Verify the command line indicates the type of data to collapse */
		if (reorderRunTimeOptions == 0) {
			ErAbort("\nThis program requires a switch to indicate the type of image being processed\n"
				"Usage: readimg -{ifd} input.name output.name\n"
				"Where: -i Indicates an integer image\n"
				"       -f Indicates a float image\n"
				"		-d Indicates a double image\n");
		}


		/* Verify command line contains two arguments */
		if ((argc - argIndex) != 2) {
			ErAbort("\nThis program requires an input and an output file\n"
				"Usage: readimg -{ifd} input.name output.name\n"
				"Where: -i Indicates an integer image\n"
				"       -f Indicates a float image\n"
				"		-d Indicates a double image\n");
		}
		
		/* Open the image file */
		if ((inputFile = LbFlFileOpen(inputName, "r+b")) == 0) {
			sprintf(reorderErrString, "Unable to open image file\n'%s' (reorder).", inputName);
			ErStFileError(reorderErrString);
			break;
		}
		
		/* Open the image file */
		if ((outputFile = LbFlFileOpen(outputName, "w+b")) == 0) {
			sprintf(reorderErrString, "Unable to create/open image file\n'%s' (reorder).", outputName);
			ErStFileError(reorderErrString);
			break;
		}
		
		/* See if they want To skip a header */
		hdrSize = LbInAskFourByte("Enter size of header to skip",
			1, false, false, 1, &reorderCanceled, 32768, 0, 0, 0);

		/* Get the number of slices */
		numSlices = LbInAskFourByte("Enter the number of slices",
			1, false, false, 0, &reorderCanceled, 0, 0, 0, 0);
		
		/* Get the number of rows */
		numRows = LbInAskFourByte("Enter the number of rows per slice",
			1, false, false, 0, &reorderCanceled, 0, 0, 0, 0);
		
		/* Get the number of columns */
		numColumns = LbInAskFourByte("Enter the number of columns per slice",
			1, false, false, 0, &reorderCanceled, 0, 0, 0, 0);
		
		
		/* Get the new order of slices */
		newSlicesOrder = LbInAsk("Enter new position of current slices", 1, true,
			&reorderCanceled, "slices", "rows", "columns", 0,
			resultStr);

		if (reorderCanceled == true)
			goto CANCELED;
		

		do { 
			/* Get the new order of rows */
			newRowsOrder = LbInAsk("Enter new position of current rows", 1, true,
				&reorderCanceled, "slices", "rows", "columns", 0,
				resultStr);

			if (reorderCanceled == true)
				goto CANCELED;
		
			if (newRowsOrder == newSlicesOrder) {
				LbInPrintf("You can't have current slices and current rows mapped to %s\n", resultStr);
			}
		} while (newRowsOrder == newSlicesOrder);
		
		/* Get the new order of columns */
		do {
			/* Get the new order of columns */
			newColumnsOrder = LbInAsk("Enter new position of current columns", 1, true,
				&reorderCanceled, "slices", "rows", "columns", 0,
				resultStr);

			if (reorderCanceled == true)
				goto CANCELED;
		
			if (newColumnsOrder == newSlicesOrder) {
				LbInPrintf("You can't have current slices and current columns mapped to %s\n", resultStr);
			}
		
			if (newColumnsOrder == newRowsOrder) {
				LbInPrintf("You can't have current rows and current columns mapped to %s\n", resultStr);
			}
		} while ((newColumnsOrder == newRowsOrder) || (newColumnsOrder == newSlicesOrder));

		/* Tell user what the new ordering is */
		LbInPrintf("New ordering is current slices = %s, current rows = %s, and current columns = %s\n",
			newOrderStr[newSlicesOrder-1], newOrderStr[newRowsOrder-1], newOrderStr[newColumnsOrder-1]);

		/* Transform ordering into correct size */
		if (newSlicesOrder == 1) {
			newSlices = numSlices;
			newIndexOrder[0] = &ii;
		}
		else if (newSlicesOrder == 2) {
			newRows = numSlices;
			newIndexOrder[0] = &jj;
		}
		else {
			newColumns = numSlices;
			newIndexOrder[0] = &kk;
		}
		
		if (newRowsOrder == 1) {
			newSlices = numRows;
			newIndexOrder[1] = &ii;
		}
		else if (newRowsOrder == 2) {
			newRows = numRows;
			newIndexOrder[1] = &jj;
		}
		else {
			newColumns = numRows;
			newIndexOrder[1] = &kk;
		}
		if (newColumnsOrder == 1) {
			newSlices = numColumns;
			newIndexOrder[2] = &ii;
		}
		else if (newColumnsOrder == 2) {
			newRows = numColumns;
			newIndexOrder[2] = &jj;
		}
		else {
			newColumns = numColumns;
			newIndexOrder[2] = &kk;
		}

		/* Turn buffering off, it shouldn't be necessary but bugs have been found in the
			DG i/o library that make it so.
		*/
		setbuf(inputFile, 0);
		setbuf(outputFile, 0);
	
		/* Seek to the end of the file */
		if (fseek(inputFile, 0, SEEK_END) != 0) {
			ErStFileError("Unable to seek to end of file (reorder)");
			break;
		}

		/* Compute the file size */
		fileSize = ftell(inputFile);
		dataSize = fileSize - hdrSize;
	
		/* Seek to the beginning of the data */
		if (fseek(inputFile, hdrSize, SEEK_SET) != 0) {
			ErStFileError("Unable to seek to beginning of file (reorder)");
			break;
		}
		
		/* Allocate memory buffer for first file, using Numerical Recipes routine */
		if (Is_Integer()){
		
			intData = (LbUsFourByte *) LbMmAlloc(dataSize);
			if (intData == 0) {
				goto FAIL;
			}
			outputInts = (LbUsFourByte *) LbMmAlloc(dataSize);
			if (outputInts == 0) {
				goto FAIL;
			}
			
			/* Read in the data */
			if ((numRead = fread(intData, dataSize, 1, inputFile)) != 1) {
				sprintf(reorderErrString, "Unable to read data for image file '%s'", inputName);
				ErStFileError(reorderErrString);
				break;
			}
		}
		else if (Is_Double()){
		
			
			doubleData = (double *) LbMmAlloc(dataSize);
			if (doubleData == 0) {
				goto FAIL;
			}
			
			outputDoubles = (double *) LbMmAlloc(dataSize);
			if (outputDoubles == 0) {
				goto FAIL;
			}
			
			/* Read in the data */
			if ((numRead = fread(doubleData, ((numRows * numColumns) *sizeof(double)), 1, inputFile)) != 1) {
				sprintf(reorderErrString, "Unable to read data for image file '%s'", inputName);
				ErStFileError(reorderErrString);
				break;
			}
		}
		else {
		
			floatData = (float *) LbMmAlloc(dataSize);
			if (floatData == 0) {
				goto FAIL;
			}
			
			outputFloats = (float *) LbMmAlloc(dataSize);
			if (outputFloats == 0) {
				goto FAIL;
			}
			
			/* Read in the data */
			if ((numRead = fread(floatData, (dataSize), 1, inputFile)) != 1) {
				sprintf(reorderErrString, "Unable to read data for image file '%s'", inputName);
				ErStFileError(reorderErrString);
				break;
			}
		}
		

		/* Now reorder the data */
		ii = jj = kk = 0;
		for (i = 0; i < numSlices; i++) {
			for (j = 0; j < numRows; j++) {
				for (k = 0; k < numColumns; k++) {
				
					if (Is_Integer())
						outputInts[(ii*newColumns*newRows)+(jj*newColumns)+kk] = intData[(i*numRows*numColumns)+(j*numColumns)+k];
					if (Is_Float())
						outputFloats[(ii*newColumns*newRows)+(jj*newColumns)+kk] = floatData[(i*numRows*numColumns)+(j*numColumns)+k];
					if (Is_Double())
						outputDoubles[(ii*newColumns*newRows)+(jj*newColumns)+kk] = doubleData[(i*numRows*numColumns)+(j*numColumns)+k];

					#ifdef HEAVY_DEBUGGING
					LbInPrintf("[%d][%d][%d]  ", (ii*newColumns*newRows),(jj*newColumns),kk);
					#endif
					

					(*(newIndexOrder[2]))++;
				}
				*(newIndexOrder[2]) = 0;
				(*(newIndexOrder[1]))++;
			}
			*(newIndexOrder[1]) = 0;
			(*(newIndexOrder[0]))++;
		}
		
		/* Write the data */
		if (Is_Integer()){
		
			
			if ((numWritten = fwrite(outputInts, dataSize, 1, outputFile)) != 1) {
				sprintf(reorderErrString, "Unable to write data for image file '%s'", outputName);
				ErStFileError(reorderErrString);
				break;
			}
		}
		else if (Is_Double()){
			
			if ((numWritten = fwrite(outputDoubles, dataSize, 1, outputFile)) != 1) {
				sprintf(reorderErrString, "Unable to write data for image file '%s'", outputName);
				ErStFileError(reorderErrString);
				break;
			}
		}
		else {
			
			if ((numWritten = fwrite(outputFloats, dataSize, 1, outputFile)) != 1) {
				sprintf(reorderErrString, "Unable to write data for image file '%s'", outputName);
				ErStFileError(reorderErrString);
				break;
			}
		}
		
	CANCELED:;
	okay = true;
	FAIL:;
	} while (false);
	
	if (intData != 0)
		LbMmFree((void **)&intData);
	if (outputInts != 0)
		LbMmFree((void **)&outputInts);
	if (floatData != 0)
		LbMmFree((void **)&floatData);
	if (outputFloats != 0)
		LbMmFree((void **)&outputFloats);
	if (doubleData != 0)
		LbMmFree((void **)&doubleData);
	if (outputDoubles != 0)
		LbMmFree((void **)&outputDoubles);
	if (inputFile != 0)
		fclose(inputFile);
	if (outputFile != 0)
		fclose(outputFile);
	
	return(okay);
}

#undef REORDER
