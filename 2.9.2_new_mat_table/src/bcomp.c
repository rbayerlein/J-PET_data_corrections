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
*			Module Name:		bcomp.c
*			Revision Number:	1.2
*			Date last revised:	23 July 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	4 April 1996
*
*      bcomp -- compare two binary files.  Stops on the first non-matching byte
*
*		SYNOPSIS
*   		   bcomp file1 file2
*
*		DESCRIPTION
*    	  
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

#include 	<stdio.h>
#include 	<string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbMemory.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbError.h"
#include "LbInterface.h"

/* Local Constants */
#define BUFF_SIZE		4096


/* Local Macros */

Boolean	bcomp(int argc, char *argv[]);

	
/* FUNCTIONS */
/**********************
*	bcomp
*
*	Purpose:	Compare two files on a binary basis.
*
*	Result:	Always returns zero indicating no error.
***********************/

Boolean bcomp(int argc, char *argv[]) {

	Boolean			okay = false;		/* Process Flag */
	char			buff1[BUFF_SIZE];	/* Our first input buffer */
	char			buff2[BUFF_SIZE];	/* Our second input buffer */
	char			errStr[1024];		/* For error strings */
	FILE			*file1 = 0;		/* The source file */
	FILE			*file2 = 0;	/* The output file */
	LbUsFourByte	index;				/* Current char index */
	LbUsFourByte	numRead1;			/* Current bytes read from file1 */
	LbUsFourByte	numRead2;			/* Current bytes read from file2 */
	LbUsFourByte	totBytes1 = 0;		/* Total bytes read from file 1 */
	LbUsFourByte	totBytes2 = 0;		/* Total bytes read from file 2 */
	
	do { /* Process Loop */
		
		/* Check for arguments */
		if (argc != 3) {
			ErStGeneric("This program requires two input files to compare as command line arguments");
			goto FAIL;
		}
		
		/* Open input file1 with read privilidge */
		if ((file1 = LbFlFileOpen(argv[1], "rb")) == 0){
			sprintf(errStr,"\nUnable to open 1st input file '%s'\n",argv[1]);
			ErStFileError(errStr);
			goto FAIL;
		}
		
		/* Open input file1 with read privilidge */
		if ((file2 = LbFlFileOpen(argv[2], "rb")) == 0) {
			sprintf(errStr, "\nUnable to open 2nd file '%s'\n",argv[2]);
			ErStFileError(errStr);
			goto FAIL;
		}
		
		/* Process current file */
		do {
		
			/* Read a block of data */
			numRead1 = fread(buff1, 1, BUFF_SIZE, file1);
			numRead2 = fread(buff2, 1, BUFF_SIZE, file2);
			totBytes1 += numRead1;
			totBytes2 += numRead2;
			
			/* Process the block if it was read properly */
			if ((numRead1 != 0) && (numRead1 == numRead2)) {
	
				for (index = 0; index < numRead1; index++) {
					if (buff1[index] != buff2[index]) {
						sprintf(errStr, "File '%s' differs from '%s' at byte %ld\n",
							argv[1], argv[2], 
							(unsigned long)(totBytes1 - (numRead1 - index)));
						ErAlert(errStr, false);
						goto DONE;
					}
				}
			}
			else if (((numRead1 == 0) && (numRead2 != 0)) ||
					((numRead2 == 0) && (numRead1 != 0))) {
				sprintf(errStr, "Files are of different length, they match upto that point\n");
				ErAlert(errStr, false);
				goto DONE;
			}
			
			/* Verify that we reached end of file, and not an error */
			if ((numRead1 == 0) && (feof(file1) == 0)){
			
				/* We were unable to read, but are not at end of file */
				sprintf(errStr, "\nUnable to read from the input file '%s'\n", argv[1]);
				ErStFileError(errStr);
				goto FAIL;
			}
			
			/* Verify that we reached end of file, and not an error */
			if ((numRead2 == 0) && (feof(file2) == 0)){
			
				/* We were unable to read, but are not at end of file */
				sprintf(errStr, "\nUnable to read from the input file '%s'\n", argv[2]);
				ErStFileError(errStr);
				goto FAIL;
			}
			
		} while (feof(file1) == 0);
			
		DONE:;
		okay = true;
		FAIL:;
	} while (false);
		
	/* Close input files */
	{
		if (file1 != 0) {
			fclose(file1);
			file1 = 0;
		}
		
		if (file2 != 0) {
			fclose(file2);
			file2 = 0;
		}
	}
	
	return (okay);
}
