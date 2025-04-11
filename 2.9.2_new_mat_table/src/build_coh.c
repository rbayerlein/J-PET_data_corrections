/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1997-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		build_coh.c
*			Revision Number:	1.1
*			Date last revised:	23 July 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	11 July 1997
*
*			Module Overview:	Combines multiple attenuation and probability
*								files into a single file.
*
*			References:			
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
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbError.h"
#include "LbInterface.h"
#include "LbMemory.h"


/* LOCAL CONSTANTS */
#define NUM_ENERGIES 1000

/* LOCAL GLOBALS */
/* PROTOTYPES */
Boolean	BuildCoh(int argc, char *argv[]);

/* FUNCTIONS */
/**********************
*	BuildCoh
*
*	Purpose:	Execute the program.
*
*	Result:	True unless an error occurs.
***********************/
Boolean BuildCoh(int argc, char *argv[])
{
	Boolean			okay = false;		/* Process Loop */
	char			errStr[1024];
	char			inputBuffer[1024];
	FILE			*masterFile;
	FILE			*inputFile;
	FILE			*outputFile;
	LbUsFourByte	numMaterials;

	do { /* Process Loop */
		
		/* Verify we have enough arguments */
		if (argc != 3) {
			LbInPrintf("This program requires a master-list input file and the name of the combined-data output file as arguments\n");
			ErAbort("Please try again");
		}
		
		/* Open the master-list input */
		if ((masterFile = LbFlFileOpen(argv[1], "r")) == 0) {
			sprintf(errStr, "Unable to open master-list input file named: '%s'. ",
				argv[1]);
			ErAbort(errStr);
		}
		
		/* Create the ouput file */
		if ((outputFile = LbFlFileOpen(argv[2], "w")) == 0) {
			sprintf(errStr, "Unable to open output file named: '%s'. ",
				argv[2]);
			ErAbort(errStr);
		}

		/* Clear counter */
		numMaterials = 0;
		
		/* Loop through and count the names */
		while (LbFlFGetS(inputBuffer, sizeof(inputBuffer), masterFile) != NULL) {
			numMaterials++;
		}
		
		/* Go back to beginning of file */
		if (fseek(masterFile, 0, SEEK_SET) != 0){
			ErStFileError("Unable to seek to beginning of master-file list.");
			goto FAIL;
		}

		/* Put the number of materials at the beginning of the file */
		fprintf(outputFile, "%ld\n", (unsigned long)numMaterials);

		/* Concatenate each file in master-file list into output file */
		while (LbFlFGetS(inputBuffer, sizeof(inputBuffer), masterFile) != NULL) {
			
			/* Put the name of the material in the output file */
			fprintf(outputFile, "%s", inputBuffer);
			
			/* Remove the new-line */
			inputBuffer[strlen(inputBuffer)-1] = '\0';
			
			/* Open the current file */
			if ((inputFile = LbFlFileOpen(inputBuffer, "r")) == 0) {
				sprintf(errStr, "Unable to open data input file named: '%s'. ",
					inputBuffer);
				ErAbort(errStr);
			}

			while (LbFlFGetS(inputBuffer, sizeof(inputBuffer), inputFile) != NULL) {
				fprintf(outputFile, "%s", inputBuffer);
			}
			
			fclose(inputFile);
			
		}
		
		
		/* Close the file */
		fclose(masterFile);
		fclose(outputFile);
		
		okay = true;
		FAIL:;
	} while (false);
	
	/* Quit the program */
	return (okay);
}
