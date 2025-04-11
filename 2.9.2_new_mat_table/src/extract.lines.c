/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1996-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		extract.lines.c
*			Revision Number:	1.4
*			Date last revised:	23 July 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	23 July, 1996
*
*			Module Overview:	Extracts a user specified number of lines from 
*								a text file.
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



/* LOCAL CONSTANTS */
#define MAXLINE	2048

/* LOCAL TYPES */

/* LOCAL GLOBALS */
static	Boolean			canceled;		/* Global cancelation flag */
static	char			errStr[1024];	/* Storrage for error messages */

/* PROTOTYPES */
Boolean	ExtractLines(int argc, char *argv[]);

/* FUNCTIONS */

/**********************
*	ExtractLine
*
*	Purpose: Extracts variable number of lines from file.
*
*	Result:	True unless an error occurs.
***********************/
Boolean ExtractLines(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */
	LbUsFourByte		localArgc;				/* Unsigned version of argc */
	char				theLine[MAXLINE];		/* Text of current line */
	LbUsFourByte		firstToExtract;			/* First line to keep */
	LbUsFourByte		numToExtract;			/* Number of lines to keep */
	LbUsFourByte		i;						/* The current file we are processing */
	LbUsFourByte		j;						/* The current line we are processing */
	LbUsFourByte		len;					/* The length of the current line we are processing */
	FILE				*inFile = 0;			/* The file we are processing */
	
	
	do { /* Process Loop */
		
		
		/* Verify command line contains two arguments */
		if (argc < 2) {
			/* Tell them a file is required */
			ErAbort("\nThis program requires at least one input file.\n");
		}

		/* Get first line to keep */
		firstToExtract = LbInAskFourByte("Enter line number to start extraction (starting at 1)",
			1, false, false, 1, &canceled, 1, 0, 0, 0);
	
		/* Get number of lines to keep */
		numToExtract = LbInAskFourByte("Enter the number of lines to extract",
			1, false, false, 1, &canceled, 1, 0, 0, 0);
		
		/* Print blank line in case redirecting input/output */
		LbInPrintf("\n");
		
		/* Loop over all input files */
		if (argc < 0) {
			localArgc = 0;
		}
		else {
			localArgc = argc;
		}
		for (i = 1; i < localArgc; i++){

			/* Open the file */
			if ((inFile = LbFlFileOpen(argv[i], "r")) == 0) {
				sprintf(errStr, "Unable to open file '%s'\n", argv[i]);
				ErStFileError(errStr);
				goto FAIL;
			}
			
			/* Loop through until we reach the first line to extract */
			for (j = 1; j < firstToExtract; j++) {
				
				/* Get the current line */
				if (LbFlFGetS(theLine, MAXLINE, inFile) == NULL) {
					sprintf(errStr, "File '%s' has only %ld lines, you requested I start on %ld",
						argv[i], (unsigned long)(j-1), (unsigned long)firstToExtract);
					ErAlert(errStr, false);
					goto NEXTFILE;
				}
				
				/* Verify that the line was terminated */
				if ((len = strlen(theLine)) == (MAXLINE-1)) {
					sprintf(errStr, "Line number %ld in file '%s' is greater than the maximum allowed line length (%d)",
						(unsigned long)j, argv[i], MAXLINE);
					ErAlert(errStr, false);
					goto NEXTFILE;
				}
			}
			
			/* Loop through until we've extracted the lines we want */
			for (j = 0; j < numToExtract; j++) {
				
				/* Get the current line */
				if (LbFlFGetS(theLine, MAXLINE, inFile) == NULL) {
					sprintf(errStr, "File '%s' has only %ld lines, you requested I start on %ld",
						argv[i], (unsigned long)(j-1), (unsigned long)firstToExtract);
					ErAlert(errStr, false);
					goto NEXTFILE;
				}
				
				/* Verify that the line was terminated */
				if ((len = strlen(theLine)) == (MAXLINE-1)) {
					sprintf(errStr, "Line number %ld in file '%s' is greater than the maximum allowed line length (%d)",
						(unsigned long)j, argv[i], MAXLINE);
					ErAlert(errStr, false);
					goto NEXTFILE;
				}
				
				/* "Extract" the line */
				LbInPrintf("%s", theLine);
			}
	
			/* This label is jumped to if the current file needs to be aborted */
			NEXTFILE:;
			
			/* Close the current file */
			if (inFile != 0) {
				fclose(inFile);		
				inFile = 0;
			}
		}
				
		okay = true;
		FAIL:;
		CANCEL:;
	} while (false);
		
	if (inFile != 0)
		fclose(inFile);		
		
	/* If error set due to cancellation, handle it, otherwise pass it on */
	if (!okay) {
		if (canceled) {
			ErHandle("User canceled ExtractLines.", false);
			okay = true;
		}
	}
	
	/* Return the status */
	return (okay);
}

