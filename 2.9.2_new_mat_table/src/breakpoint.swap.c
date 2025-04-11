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
*			Module Name:		breakpoint.swap.c
*			Revision Number:	1.1
*			Date last revised:	6 September 2012
*			Programmer:			Steven Vannoy
*			Date Originated:	15 December 1994
*
*			Module Overview:	Swaps the bytes in a GE produced "breakpoint"
*								file.  The file has three parts, the first is
*								a 256 byte header which we don't use so we don't
*								 swap it.
*								The second section is 336 floats (scale factors
*								for each projection plane).  Finally 2 byte integer
*								data for the rest of the file.
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

/* This dictates the size of the units we read from disk.
It will be applied to a "generic" buffer which then gets
interpreted based on the type of the binned data. Hence,
it must remain a power of 2 >= 8.
*/
#define BUFF_SIZE (1024*1024)

/* LOCAL TYPES */

/* LOCAL GLOBALS */
static	Boolean	canceled;	/* Global cancelation flag */

/* PROTOTYPES */
Boolean			BreakpointSwap(int argc, char *argv[]);

/* FUNCTIONS */

/**********************
*	BreakpointSwap
*
*	Purpose:	Combine the history files.
*
*	Result:	True unless an error occurs.
***********************/
Boolean BreakpointSwap(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */
	char				errStr[1024];			/* Error string buffer */
	FILE				*finalFile;				/* The output file */
	LbUsFourByte		index;					/* Current element in buffer */
	LbUsFourByte		fileIndex;				/* Current element in buffer */
	LbUsFourByte		numToProcess;			/* Number of files to process */
	LbUsFourByte		elemSize;				/* Size of "elements" per data buffer read */
	LbUsFourByte		numRead;				/* Number of bytes read from final file */
	LbUsFourByte		byteIndex;				/* LCV */
	char				*oldData;				/* Buffer for unswapped data */
	char				*newData;				/* Buffer for swapped data */
	float				scaleData[336];			/* Scaled factor data at start of file (after header) */
	
	
	do { /* Process Loop */
			
		/* Verify command line contains two arguments */
		if (argc < 2) {
			/* Ask for the file name */
			ErAbort("\nThis program requires a file name for input.\n");
		}
		
		/* Loop through all files on the command line */
		numToProcess = argc - 1;
		for (fileIndex = 1; fileIndex <= numToProcess; fileIndex++) {
		
			/* Open the first file, the one that will have the final content */
			if ((finalFile = LbFlFileOpen(argv[fileIndex], "r+b")) == 0) {
				sprintf(errStr, "Unable to open output file\n'%s'.", argv[fileIndex]);
				ErStFileError(errStr);
				goto FAIL;
			}			
		
			/* Allocate memory buffers */
			if ((oldData = (char *) LbMmAlloc(BUFF_SIZE)) == 0) {
				goto FAIL;
			}
			
			if ((newData = (char *) LbMmAlloc(BUFF_SIZE)) == 0) {
				goto FAIL;
			}
			
			/* Seek past the 256 byte header */
			if (fseek(finalFile, 256, SEEK_SET) != 0) {
				sprintf(errStr, "Unable to seek past 256 byte header.\n");
				ErStFileError(errStr);
				goto FAIL;
			}
			
			/* Read in the 336 scale factors */
			if ((numRead = fread((void *)scaleData, sizeof(float), 336, finalFile)) != 336){
				sprintf(errStr, "Unable to read 336 scale factors.\n");
				ErStFileError(errStr);
				goto FAIL;
			}
			
			/* Swap the scale factors */
			for (index = 0; index < 336; index++)
				LbCvSwap((char *)scaleData+(index*sizeof(float)),
					(char *)newData+(index*sizeof(float)), sizeof(float));
			
			/* Seek back to beginnig of 336 scale factors */
			if (fseek(finalFile, -(336*sizeof(float)), SEEK_CUR) != 0) {
				sprintf(errStr, "Unable to seek back to beginning of scale factors.\n");
				ErStFileError(errStr);
				goto FAIL;
			}
			
			/* Write the data to the output file */
			if (fwrite(newData, sizeof(float), 336, finalFile) != 336) {
				ErStFileError("\nUnable to write swapped scale factors to file.");
				goto FAIL;
			}

			/* Compute elements per buffer */
			elemSize = sizeof(short int);
			
			
			/* Loop by reading a buffer from the final file and processing until there
				isn't a full buffer left (the partial buffer is handled next)
			 */
			while ((numRead = fread(oldData, 1, BUFF_SIZE, finalFile)) == BUFF_SIZE) {

				/* Loop through, swapping element bytes */
				for (index = 0; index < BUFF_SIZE; index += elemSize) {
				
					/* Swap the data */
					for (byteIndex = 0; byteIndex < elemSize; byteIndex++) {
						newData[index+((elemSize-1)-byteIndex)] =
							oldData[index+byteIndex];
					}
				}
				
				/* Seek back one buffer in final file */
				if (fseek(finalFile, -BUFF_SIZE, SEEK_CUR) != 0) {
					sprintf(errStr, "Unable to seek back in final file.\n");
					ErStFileError(errStr);
					goto FAIL;
				}
				
				/* Write the data to the output file */
				if (fwrite(newData, BUFF_SIZE, 1, finalFile) != 1) {
					ErStFileError("\nUnable to write updated data to final file.");
					goto FAIL;
				}
			}
			
			/* Process partial buffer if there is  one */
			if (numRead != 0) {
				
				/* Loop through, adding elements together */
				for (index = 0; index < numRead; index += elemSize) {
					
					/* Swap the data */
					for (byteIndex = 0; byteIndex < elemSize; byteIndex++) {
						newData[index+((elemSize-1)-byteIndex)] =
							oldData[index+byteIndex];
						}
				}
				
				/* Seek back one buffer in final file */
				if (fseek(finalFile, -numRead, SEEK_CUR) != 0) {
					sprintf(errStr, "Unable to seek back in final file.\n");
					ErStFileError(errStr);
					goto FAIL;
				}
				
				/* Write the data to the output file */
				if (fwrite(newData, numRead, 1, finalFile) != 1) {
					ErStFileError("\nUnable to write updated data to final file.");
					goto FAIL;
				}
			}
			/* If we are here and haven't read to EOF there is an error */
			else if (feof(finalFile) == 0) {
				sprintf(errStr, "Unable to read from file '%s'.\n", argv[fileIndex]);
				ErStFileError(errStr);
				goto FAIL;
			}
				
			/* Close the final file */
			fclose(finalFile);
		}
		
		okay = true;
		FAIL:;
	} while (false);
	
	/* Free memory and close files */
	if (oldData != 0)
		LbMmFree((void **) &oldData);
		
	if (newData != 0)
		LbMmFree((void **) &newData);
		
	if (finalFile != 0)
		fclose(finalFile);
		
	/* If error set due to cancellation, handle it, otherwise pass it on */
	if (!okay) {
		if (canceled) {
			ErHandle("User canceled breakpoint.swap.", false);
			okay = true;
		}
	}
	
	/* Return the status */
	return (okay);
}
