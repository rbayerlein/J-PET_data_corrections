/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1994-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		reverse.bytes.c
*			Revision Number:	1.6
*			Date last revised:	23 July 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	May 12, 1994
*
*			Module Overview:	"Reverses" the byte ordering of a file. The user
*								specifies the unit size to be reversed, for
*								example a unit size of 4 would cause groups
*								of four bytes to be reversed, while a unit
*								size of eight would cause groups of eight
*								to be reversed.
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
#include "LbMemory.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbError.h"
#include "LbInterface.h"
#include "LbDebug.h"

/* LOCAL CONSTANTS */
/* This dictates the size of the units we read from disk.
It will be applied to a "generic" buffer which then gets
interpreted based on the type of the binned data. Hence,
it must remain a power of 2 >= 8.
*/
#define BUFF_SIZE (1024*1024)

/* LOCAL TYPES */

/* LOCAL GLOBALS */
static	char		inputBuffer[BUFF_SIZE];			/* Input buffer */
static	char		outputBuffer[BUFF_SIZE];		/* Output buffer */

/* PROTOTYPES */
Boolean			ReverseBytes(int argc, char *argv[]);

/* FUNCTIONS */


/**********************
*	ReverseBytes
*
*	Purpose:	Execute the program.
*
*	Result:	Always returns zero indicating no error.
***********************/
Boolean ReverseBytes(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */
	Boolean				canceled = false;		/* Cancelation flag */
	char				errStr[1024];			/* Error string */
	FILE				*inputFile = 0;			/* Our data file */
	FILE				*outputFile = 0;		/* Converted file */
	LbUsFourByte		hdrSize;				/* Bytes Read */
	LbUsFourByte		bytesRead;				/* Bytes Read */
	LbUsFourByte		numBytes;				/* Number of bytes to swap */
	LbUsFourByte		byteIndex;				/* Loop control */
	/*LbUsFourByte		numBytesSwapped = 0;*/	/* For statistics */
	LbUsFourByte		i;						/* LCV */
	
	do { /* Process Loop */
		
		/* Check for proper number of arguments */
		if (argc != 4) {
			ErStFileError("This program requires 3 arguments, \n an input file, an output file, and the number of bytes to swap.");
			break;
		}
		
		
		/* Open data file file */
		if ((inputFile = LbFlFileOpen(argv[1], "rb")) == 0) {
			sprintf(errStr, "Unable to open input file named '%s'", argv[1]);
			ErStFileError(errStr);
			break;
		}
		

		/* Open the output file */
		if ((outputFile = LbFlFileOpen(argv[2], "wb")) == 0) {
			sprintf(errStr, "Unable to create/open output file named '%s'", argv[2]);
			ErStFileError(errStr);
			break;
		}		
		
		/* Save the number of bytes to swap */
		numBytes = atoi(argv[3]);
		if ((numBytes != 2) && (numBytes != 4) && (numBytes != 8)) {
			sprintf(errStr, "The last argument, the number of bytes to swap, must be a multiple of 2 between 2 and 8, you gave %ld\n", (unsigned long)numBytes);
			ErStGeneric(errStr);
			break;
		}
		
		/* See if they want To skip a header */
		hdrSize = LbInAskFourByte("Enter size of header to skip",
			1, false, false, 1, &canceled, 32768, 0, 0, 0);
		
		if (canceled)
			break;

		/* Seek past header if requested */
		if (hdrSize != 0) {

			/* Seek past the header */
			if (fseek(inputFile, hdrSize, SEEK_SET) != 0) {
				ErStFileError("\nUnable to seek past header.");
				break;
			}
		}

		/* Loop through the files */
		do {
			okay = false;
						
			/* Read the input file */
			if ((bytesRead = fread(inputBuffer, 1, BUFF_SIZE, inputFile))
					 == BUFF_SIZE) {
				
				/* Loop through the data swapping the bytes */
				for (i = 0; i < BUFF_SIZE; i += numBytes) {
					
					/* Swap the data */
					for (byteIndex = 0; byteIndex < numBytes; byteIndex++) {
						outputBuffer[i+((numBytes-1)-byteIndex)] =
							inputBuffer[i+byteIndex];
					}
				}
					
				/* Write the output file */
				if (fwrite(outputBuffer, 1, BUFF_SIZE, outputFile) != BUFF_SIZE) {
					ErStFileError("\nUnable to write output file.");
					break;
				}
				
			}
			else if (feof(inputFile) == 0){
			
				/* Error is not end of file so let them know */
				ErStFileError("\nUnable to read image file.");
				break;
			}
			else if (bytesRead != 0) {
				if ((bytesRead % numBytes) == 0) {
					/* Loop through the data swapping the bytes */
					for (i = 0; i < BUFF_SIZE; i += numBytes) {
						
						/* Swap the data */
						for (byteIndex = 0; byteIndex < numBytes; byteIndex++) {
							outputBuffer[i+((numBytes-1)-byteIndex)] =
								inputBuffer[i+byteIndex];
						}
					}
					/* Write the output file */
					if (fwrite(outputBuffer, 1, bytesRead, outputFile) != bytesRead) {
						ErStFileError("\nUnable to write output file.");
						break;
					}
				}
				else {
				
					ErStGeneric("\nFile is not evenly divisible by 'numBytes', processing done upto last buffer\n");
					break;
				}
			}
			else {
				/* We reached the end of the file and are done */
				okay = true;
				break;
			}
		} while (true);
	} while (false);
	
	/* Close the files */
	if (inputFile != 0)
		fclose(inputFile);
		
	if (outputFile != 0)
		fclose(outputFile);
		
	/* Quit the program */
	return (okay);
}
