/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1995-2011 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		strip.header.c
*			Revision Number:	1.1
*			Date last revised:	12 December 2011
*			Programmer:			Steven Vannoy
*			Date Originated:	February 9, 1995
*
*			Module Overview:	Strips an arbitrarily sized header off a file.
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
#define BUFF_SIZE 2048

/* LOCAL TYPES */

/* LOCAL GLOBALS */

/* PROTOTYPES */
Boolean			stripheader(int argc, char *argv[]);

/* FUNCTIONS */


/**********************
*	stripheader
*
*	Purpose:	Execute the program.
*
*	Result:	Always returns zero indicating no error.
***********************/
Boolean stripheader(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */
	Boolean				canceled = false;		/* Cancelation flag */
	char				errStr[1024];			/* Error string */
	FILE				*inputFile=0;				/* Our data file */
	FILE				*outputFile=0;			/* Converted file */
	char				data[BUFF_SIZE];		/* Buffer */
	LbUsFourByte		hdrSize;				/* Size of header */
	LbUsFourByte		numRead;				/* Bytes read */

	do { /* Process Loop */
		
		/* Check for proper number of arguments */
		if (argc != 3) {
			ErStFileError("This program requires 2 arguments, \n an input file, and an output file.\n");
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
		
		/* See if they want To skip a header */
		hdrSize = LbInAskFourByte("Enter size of header to skip",
			1, false, false, 1, &canceled, 31768, 0, 0, 0);
		
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

		/* Loop by reading a buffer from the final file and processing until there
			isn't a full buffer left (the partial buffer is handled next)
		 */
		while ((numRead = fread(data, 1, BUFF_SIZE, inputFile)) == BUFF_SIZE) {
						
			
			/* Write the data to the output file */
			if (fwrite(data, BUFF_SIZE, 1, outputFile) != 1) {
				ErStFileError("\nUnable to write to output file.");
				goto FAIL;
			}
		}
		
		/* Process partial buffer if there is  one */
		if (numRead != 0) {
			
			/* Write the data to the output file */
			if (fwrite(data, numRead, 1, outputFile) != 1) {
				ErStFileError("\nUnable to write to output file.");
				goto FAIL;
			}
		}
		/* If we are here and haven't read to EOF there is an error */
		else if (feof(inputFile) == 0) {
			sprintf(errStr, "Unable to read from file '%s'.\n", argv[1]);
			ErStFileError(errStr);
			goto FAIL;
		}

		okay = true;
		FAIL:;
	} while (false);
	
		
	/* Close the files */
	if (inputFile != 0)
		fclose(inputFile);
		
	if (outputFile != 0)
		fclose(outputFile);
		
	/* Quit the program */
	return (okay);
}
