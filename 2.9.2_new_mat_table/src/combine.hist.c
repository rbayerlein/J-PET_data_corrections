/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1990-2012 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		combine.hist.c
*			Revision Number:	1.1
*			Date last revised:	22 October 2012
*			Programmer:			Steven Vannoy
*			Date Originated:	21 September, 1993
*
*			Module Overview:	Combines a series of history files into one.
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
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhgMath.h"
#include "PhoHFile.h"
#include "ProdTbl.h"
#include "PhoTrk.h"
#include "SubObj.h"
#include "EmisList.h"
#include "Collimator.h"
#include "Detector.h"
#include "phg.h"
#include "PhgBin.h"

/* LOCAL CONSTANTS */
#define DATA_BUFF_SIZE 1000000

/* LOCAL TYPES */

/* LOCAL GLOBALS */
static	Boolean	canceled;	/* Global cancelation flag */
static	PhoHFileHdrTy		newHeader;						/* Our Header */
static	PhoHFileHdrTy		header;						/* Our Header */

/* PROTOTYPES */
Boolean			CombineHist(int argc, char *argv[]);

/* FUNCTIONS */
/**********************
*	CombineHist
*
*	Purpose:	Combine the history files.
*
*	Result:	True unless an error occurs.
***********************/
Boolean CombineHist(int argc, char *argv[])
{
	Boolean				okay = false;				/* Process Loop */
	Boolean				updateHdr = false;			/* Will we update the header */
	Boolean				deleteInputs = false;		/* Will we delete input files */
	char				errStr[1024];				/* Error string buffer */
	char				*hdrBuffPtr;				/* Buffer for reading in headers */
	char				*dataBuffPtr;				/* Buffer for reading in data */
	FILE				*inputFile;					/* Current input file */
	FILE				*outputFile;				/* The output file */
	LbUsFourByte		hdrSize;					/* Size of the file's header */
	LbUsFourByte		curFileIndex;				/* Current file index */
	LbUsFourByte		numToProcess;				/* Number of files to process */
	LbUsFourByte		numRead;					/* Number of bytes read */
	
	#ifdef GEN_UNIX
	struct	stat		statBuff;					/* Buffer for file status */
	int					statResult;					/* Result of stat call */
	#endif

	do { /* Process Loop */
			
		/* Verify command line is reasonable */
		if (argc < 3) {
			/* Ask for the file name */
			ErAbort("\nThis program requires at least 2 file names for input.\n");
		}

		#ifdef GEN_UNIX
		/* Check for existance of output file */
		statResult = stat(argv[argc-1], &statBuff);
		
		if (statResult == 0) {
			if (LbInAskYesNo("Output file allready exists, overwrite", LBINYes)
					== LBINNo) {
				ErAbort("\nUser cancelled execution.\n");
			}
		}
		else if (errno != ENOENT) {
				ErAbort("\nUnable to access output file.\n");
		}
		#endif
		
		/* Open output file */
		if ((outputFile = LbFlFileOpen(argv[argc-1], "wb")) == 0) {
			sprintf(errStr, "Unable to open output file\n'%s'.", argv[argc-1]);
			ErStFileError(errStr);
			break;
		}
		
		/* See if we will delete input files on the fly */
		if (LbInAskYesNo("Do you want the input files deleted", LBINYes)
					== LBINYes) {
			deleteInputs = true;
		}
		
		/* Set file count vars */
		numToProcess = argc - 2;
		curFileIndex = 1;
		
		/* Allocate memory for data buffer */
		if ((dataBuffPtr = (char *)LbMmAlloc(DATA_BUFF_SIZE)) == 0) {
			break;
		}
		
		/* Open first input file */
		if ((inputFile = LbFlFileOpen(argv[curFileIndex], "rb")) == 0) {
			sprintf(errStr, "Unable to open input file\n'%s'.", argv[curFileIndex]);
			ErStFileError(errStr);
			break;
		}

		/* Read the header in the first input file */
		{
			/* First, read the header size; it is in the first four bytes */
			if (fread(&hdrSize, sizeof(LbUsFourByte), 1, inputFile) != 1) {
				sprintf(errStr, "\nUnable to read header size from input file '%s'.\n",
					argv[curFileIndex]);
				ErStFileError(errStr);
				break;
			}
			
			/* Verify we know how to deal with this */
			if (hdrSize == sizeof(header)) {
				updateHdr = true;
			}
			
			/* Allocate header buffer if unknown */
			if (updateHdr == false) {
				if ((hdrBuffPtr = (char *) LbMmAlloc(hdrSize)) == 0) {
					break;
				}
			}
			else {
				hdrBuffPtr = (char *) &header;
			}
			
			/* Reset to zero and read in the header */
			if (fseek(inputFile, 0, SEEK_SET) != 0) {
				ErStFileError("\nUnable to reset to beginning of history file.");
				break;
			}

			/* Read in the header */
			if (fread(hdrBuffPtr, hdrSize, 1, inputFile) != 1) {
				ErStFileError("\nUnable to read header from history file.");
				break;
			}

			/* Save the header */
			if (updateHdr == true){
				memcpy(&newHeader, &header, sizeof(header));
			}

		}

		/* Write the header to the output file */
		if (fwrite(hdrBuffPtr, hdrSize, 1, outputFile) != 1) {
			ErStFileError("\nUnable to write header to new history file.");
			break;
		}
			
		/* Process the files */
		for (curFileIndex = 1; curFileIndex <= numToProcess; curFileIndex++){

			/* Read until EOF */
			while ((numRead = fread(dataBuffPtr, 1, DATA_BUFF_SIZE, inputFile)) == DATA_BUFF_SIZE) {
				/* Write the data to the output file */
				if (fwrite(dataBuffPtr, DATA_BUFF_SIZE, 1, outputFile) != 1) {
					ErStFileError("\nUnable to write to new history file.");
					goto FAIL;
				}
			}
			
			/* See if we reached EOF */
			if (feof(inputFile)) {
				/* Write remainder of file */
				if (fwrite(dataBuffPtr, numRead, 1, outputFile) != 1) {
					ErStFileError("\nUnable to write to new history file.");
					goto FAIL;
				}
			}
			else {
				sprintf(errStr, "Unable to read from file '%s'.\n", argv[curFileIndex]);
				ErStFileError(errStr);
				goto FAIL;
			}
			
			/* Close the input file */
			fclose(inputFile);
			
			/* Delete input file if requested */
			if (deleteInputs == true) {
				if (remove(argv[curFileIndex]) != 0) {
					sprintf(errStr, "Unable to delete file '%s'.\n", argv[curFileIndex]);
					ErStFileError(errStr);
					goto FAIL;
				}
			}
			
			/* Clear header buffer if updating = false */
			if (updateHdr == false) {
				LbMmFree((void **)&hdrBuffPtr);
			}
			
			/* Open next file */
			if (curFileIndex < numToProcess) {
				
				/* Open input file */
				if ((inputFile = LbFlFileOpen(argv[curFileIndex+1], "rb")) == 0) {
					sprintf(errStr, "Unable to open input file\n'%s'.", argv[curFileIndex+1]);
					ErStFileError(errStr);
					break;
				}
				/* First, read the header size; it is in the first four bytes */
				if (fread(&hdrSize, sizeof(LbUsFourByte), 1, inputFile) != 1) {
					ErStFileError("\nUnable to read header size from input file.");
					break;
				}
				
				/* Allocate header buffer if unknown */
				if (updateHdr == false) {
					if ((hdrBuffPtr = (char *)LbMmAlloc(hdrSize)) == 0) {
						break;
					}
				}
				
				/* Reset to zero and read in the header */
				if (fseek(inputFile, 0, SEEK_SET) != 0) {
					ErStFileError("\nUnable to reset to beginning of history file.");
					break;
				}
	
				/* Read in the header */
				if (fread(hdrBuffPtr, hdrSize, 1, inputFile) != 1) {
					ErStFileError("\nUnable to read header from history file.");
					break;
				}
				
				/* Update header if requested */
				if (updateHdr == true) {
					newHeader.H.NumPhotons += header.H.NumPhotons;
					newHeader.H.NumDecays += header.H.NumDecays;
					newHeader.H.PhgRunTimeParams.Phg_EventsToSimulate +=
						header.H.PhgRunTimeParams.Phg_EventsToSimulate;
				}
					
			}
		} /* End FOR each input file */

		/* Update header if set up to */
		if (updateHdr == true){
				
	   		/* Reset to zero and write the header */
	   		if (fseek(outputFile, 0, SEEK_SET) != 0) {
	   			ErStFileError("\nUnable to reset to beginning of history file.");
	   			break;
	   		}

			/* Write the header to the output file */
			if (fwrite(&newHeader, hdrSize, 1, outputFile) != 1) {
				ErStFileError("\nUnable to write header to new history file.");
				break;
			}

	   	}

		okay = true;
		FAIL:;
	} while (false);
	
	/* Handle error situation if one exists */
	if (!okay) {
		if (canceled)
			ErHandle("User canceled combine.hist.", false);
			okay = true;
	}
	
	/* Quit the program */
	return (okay);
}
