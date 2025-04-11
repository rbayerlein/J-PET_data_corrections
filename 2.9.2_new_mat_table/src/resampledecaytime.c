/*********************************************************************************
*                                                                                
*                       Source code developed by the                             
*           Imaging Research Laboratory - University of Washington               
*               (C) Copyright 2003-2012 Department of Radiology          	        
*                           University of Washington                             
*                              All Rights Reserved                               
*                                                                                
*********************************************************************************/

/*********************************************************************************
*
*		Module Name:		resampledecaytime.c
*		Revision Number:	1.4
*		Date last revised:	10 October 2012
*		Programmer:			Robert Harrison
*		Date Originated:	15 Aug 2005
*
*		Module Overview:	Scan a previously sorted list mode file and add
*							randoms, delete triples.
*
*		References:			RandomsCodeChanges.doc
*
**********************************************************************************
*
*		Global functions defined:		None.
*
*		Global macros defined:			None.
*				
*		Global variables defined:		None.
*
**********************************************************************************
*
*		Revision Section (Also update version number, if relevant)
*
*		Programmer(s):		
*
*		Revision date:		
*
*		Revision description:
*
*********************************************************************************/

#define RESAMPLE_DECAY_TIME_MAIN	/* Note we are substituting ourselves for phg's main */

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"
#include "LbHeader.h"
#include "LbTiming.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhgMath.h"
#include "PhoHFile.h"
#include "PhgHdr.h"
#include "ProdTbl.h"
#include "PhoTrk.h"
#include "SubObj.h"
#include "EmisList.h"
#include "phg.h"
#include "Collimator.h"
#include "Detector.h"
#include "PhgBin.h"


#ifdef MPW
	#pragma segment RESAMPLE_DECAY_TIME_MAIN
#endif

/* LOCAL CONSTANTS */
#define RESAMPT_IsDeleteFile()		LbFgIsSet(PhgOptions, LBFlag0)		/* Will we delete the original file */

#define	RESAMPT_NumFlags	1		/* Number of flags defined */

/* LOCAL TYPES */
typedef enum  {Null, Decay, Photon} EventTy;


/* LOCAL GLOBALS */
static Boolean				resamptCancelled;				/* Global cancellation flag */
static Boolean				resamptCustomFile = false;		/* True if custom history file used */
static char					resamptErrStr[1024];			/* Error string storage */		
static char					resamptInHistName[1024];		/* Name of input history file */
static char					resamptOutHistName[256];		/* Name of output history file */
static char					resamptHistParamsName[256];		/* Name of history parameters file */
static LbUsFourByte			resamptArgIndex;				/* Index through command line arguments */
static PhoHFileHdrTy		resamptHdrParams;				/* Input header */
static PhoHFileHkTy			resamptHistParamsHk;			/* Hook to custom history file information */
static PhoHFileHkTy			resamptInHistHk;				/* Hook to input history file */
static PhoHFileHkTy			resamptOutHistHk;				/* Hook to output history file for resampled decay times */
static LbUsFourByte			numDecaysRead;					/* number decays read from input */
static LbUsFourByte			numDecaysWritten;				/* number decays written to output */


/* PROTOTYPES */
Boolean			resampt(int argc, char *argv[]);
Boolean 		resamptInitialize(int argc, char *argv[]);
void			resamptTerminate(void);
Boolean			resamptResampleTime(char *argv[]);
Boolean			resamptNewDecayTimes(double scanTime);


/* FUNCTIONS */
/**********************
*	Name:		resampt
*
*	Purpose:	Execute the program.
*
*	Result:		True unless an error occurs.
***********************/
Boolean resampt(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */
					/* Error condition results unless this is changed in body */
	time_t 			curTime;		/* Current time for stamping execution date */
	LbFourByte		argIndex;		/* Index through command line arguments */
	char			executionTimeStr[33];
									/* String for execution date conversion */
	
	/*** NOTE:  Custom history files are not supported and have not been tested. ***/
	
	
	do { /* Process Loop */
	
		/* Perform initialization tasks */
		if (resamptInitialize(argc, argv) == false) {
			sprintf(resamptErrStr, "Resample Decay Time initialization failed.");
			ErStGeneric(resamptErrStr);
			goto FAIL;
		}
		
		/* Print out the command line */
		LbInPrintf("\n\n\nCommand line: ");
		for (argIndex = 0; argIndex < argc; argIndex++)
			LbInPrintf("%s ", argv[argIndex]);
		LbInPrintf("\n");
		
		/* Get current time */
		time(&curTime);
		
		/* Convert time to string format */
		strftime(executionTimeStr, 31, "%H:%M:%S %e-%b-%y", localtime(&curTime));
				
		LbInPrintf("\nExecution of resampledecaytime occurred on %s\n", executionTimeStr);
		
		/* Print the input and output list mode files for resampling */
		LbInPrintf("\nName of input history file to resample: %s.\n", resamptInHistName);
		LbInPrintf("Name of output history file for resampled data: %s.\n", resamptOutHistName);

		/* Scramble the decay times */
		okay = resamptResampleTime(argv);
		if (! okay) {
			/* Error should already have been recorded */
			goto FAIL;
		}
		
		/* write out report for module */
		LbInPrintf("\n\tNumber of decays read in: %d", numDecaysRead);
		LbInPrintf("\n\tNumber of decays written out: %d", numDecaysWritten);
		
		/* Delete the original file, if requested */
		if (RESAMPT_IsDeleteFile()) {
			LbInPrintf("Deleting input history file, as requested.\n");
			if (remove(resamptInHistName) != 0) {
				sprintf(resamptErrStr, "Unable to delete history file\n'%s'.",
					resamptInHistName);
				ErStFileError(resamptErrStr);
				/* But keep going */
			}
		}
		
		/* Get current time */
		time(&curTime);
		
		/* Convert time to string format */
		strftime(executionTimeStr, 31, "%H:%M:%S %e-%b-%y", localtime(&curTime));
			
		/* Print out our final message */
		LbInPrintf("\n\n\nExecution of resampledecaytime finished on %s\n", executionTimeStr);
		
		okay = true;
		FAIL:;
		CANCEL:;
		
	} while (false);
	
	/* Terminate PHG modules */
	resamptTerminate();
	
	/* Handle error situation if one exists */
	if (!okay && resamptCancelled) {
		ErHandle("User cancelled resampledecaytime", false);
		okay = true;
	}
	
	/* Quit the program */
	return (okay);
}


/**********************
*	Name:		resamptInitialize
*
*	Purpose:	Perform initialization tasks.
*
*	Result:		True unless an error occurs.
***********************/
Boolean resamptInitialize(int argc, char *argv[])
{
	Boolean				okay = false;				/* Process Loop */
	char				*knownOptions[] = {"r"};
	char				optArgs[RESAMPT_NumFlags][LBEnMxArgLen];
	LbUsFourByte		optArgFlags = (LBFlag0);
	LbFourByte			randSeed;
	
	
	do { /* Process Loop */
		
		/* Get our options */
		if (!LbEnGetOptions(argc, argv, knownOptions,
				&PhgOptions, optArgs, optArgFlags, &resamptArgIndex)) {
	
			break;
		}
				
		/* See if they supplied two command line arguments */
		if ( (argc - resamptArgIndex) == 2) {
		
			/* Get the input history file */
			strcpy(resamptInHistName,argv[resamptArgIndex]);
			
			/* Get the output history file */
			resamptArgIndex++;
			strcpy(resamptOutHistName,argv[resamptArgIndex]);
			
		}
		else {
			
			/* Report error */					
			ErStGeneric("ERROR - Usage: resampledecaytime inputHistoryFile outputHistoryFile.");
			break;
			
		}
		
		/* Initialize the math library */
		randSeed = PhgRunTimeParams.PhgRandomSeed;
		if (!PhgMathInit(&randSeed))
			break;
		/* Save the seed if it had come from the clock */
		if (PhgRunTimeParams.PhgRandomSeed == 0)
			PhgRunTimeParams.PhgRandomSeed = randSeed;
		
		okay = true;
		CANCEL:;
	} while (false);
	
	return (okay);
}



/**********************
*	Name:		resamptTerminate
*
*	Purpose:	Terminate all of the managers.
*
*	Result:		None.
***********************/
void resamptTerminate()
{
	ProdTblProdTblInfoTy	hstPrdTblInfo;		/* Info for initializing productivity table */
	
	
	/* NOTE:  None of these should have been initialized, anyway */
	
	
	/* Terminate the binning module if initialized */
	if (PHG_IsBinOnTheFly()) {
		PhgBinTerminate(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0]);
	}
	/* Terminate the collimation module if initialized */
	if (PHG_IsCollimateOnTheFly()) {
		for (ColCurParams = 0; ColCurParams < ColNumParams; ColCurParams++) {
			ColTerminate();
		}
		ColCurParams = 0;
	}
	/* Terminate the detection module if initialized */
	if (PHG_IsDetectOnTheFly()) {
		DetTerminate();
	}
	
	/* Terminate the photon tracking module */
	PhoTrkTerminate();
	
	/* Terminate the emission list manager */
	EmisListTerminate();
	
	/* Terminate the Productivity table manager */
	ProdTblTerminate(&hstPrdTblInfo);
	
	/* Terminate the sub object module */
	SubObjTerminate();
}


/*********************************************************************************
*
*	Name:			resamptResampleTime
*
*	Summary:		Give every decay in history file a new decay time randomly
*					sampled from the scan time.
*
*	Arguments:
*		char			*argv[]					- Command line arguments.
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean resamptResampleTime(char *argv[])
{
	Boolean				okay = false;		/* Process flag (function return) */
	FILE				*inHistoryFile;		/* The input history file */
	double				scanTime;			/* Scan time from history file header */
	
	
	/* Remove compiler warning */
	if (argv) {};
	
	do { /* Process Loop */
	
		/* open files and verify parameters */
		{
		
			/* Validate the name of the input history file */
			if (strlen(resamptInHistName) == 0) {
				sprintf(resamptErrStr, "No input history file name supplied.");
				ErStGeneric(resamptErrStr);
				goto FAIL;
			}
			
			/* Validate the name of the resampled decay-time output file */
			if (strlen(resamptOutHistName) == 0) {
				sprintf(resamptErrStr, "No resampled decay-time output file name supplied.");
				ErStGeneric(resamptErrStr);
				goto FAIL;
			}
			
			/* Open input history file */
			if ((inHistoryFile = LbFlFileOpen(resamptInHistName, "rb")) == 0) {
				sprintf(resamptErrStr, "Unable to open input history file\n'%s'.",
					resamptInHistName);
				ErStFileError(resamptErrStr);
				goto FAIL;
			}
			
			resamptInHistHk.histFile = inHistoryFile;
			
			/* Init header hook, read in input header and copy to output header and save it */
			if (PhgHdrGtParams(resamptInHistHk.histFile, &resamptHdrParams, &resamptInHistHk.headerHk) == false){
				sprintf(resamptErrStr, "Unable to read input history file header\n'%s'.",
					resamptInHistName);
				ErStFileError(resamptErrStr);
				goto FAIL;
			}
			
			/* set local variable for scan time to make things read more clearly */
			scanTime = (double)(resamptHdrParams.H.PhgRunTimeParams.Phg_LengthOfScan);
			
			/* Open and initialize the custom history file parameters, if any */
			{
				
				if (resamptCustomFile) {
					/* Setup header hook */
					
					resamptHistParamsHk.doCustom = true;
					
					/* Do custom parameter initialization  */
					if (PhoHFileGetRunTimeParams(resamptHistParamsName, &(resamptHistParamsHk.customParams)) == false) {
						sprintf(resamptErrStr,"Unable to get custom parameters for history file named '%s'",
							resamptHistParamsName);
						ErAlert(resamptErrStr, false);
						goto FAIL;
					}
					
					if ( resamptHistParamsHk.customParams.doDecayTime == false ) {
						sprintf(resamptErrStr, "Randoms processing of custom list mode requires doDecayTime = true.");
						ErStFileError(resamptErrStr);
						goto FAIL;
					}
					
					if ( resamptHistParamsHk.customParams.doDecayType == false ) {
						sprintf(resamptErrStr, "Randoms processing of custom list mode requires doDecayTye = true.");
						ErStFileError(resamptErrStr);
						goto FAIL;
					}
					
					if ( resamptHistParamsHk.customParams.doTravelDistance == false ) {
						sprintf(resamptErrStr, "Randoms processing of custom list mode requires doTravelDistance = true.");
						ErStFileError(resamptErrStr);
						goto FAIL;
					}
					
					strcpy(resamptHistParamsHk.customParamsName, resamptHistParamsName);
					
					resamptHistParamsHk.bluesReceived = 0;
					resamptHistParamsHk.bluesAccepted = 0;
					resamptHistParamsHk.pinksReceived = 0;
					resamptHistParamsHk.pinksAccepted = 0;
				}
				
			}
			
			/* Open the output history file */
			if (PhoHFileCreate(resamptOutHistName, "", 
					resamptHdrParams.H.HdrKind, &resamptOutHistHk) == false) {
				sprintf(resamptErrStr,"Unable to create output history file named:\n"
					"'%s'\n"
					" (resamptResampleTime)",
					resamptOutHistName);
				ErStFileError(resamptErrStr);
				goto FAIL;
			}
			
		}
		
		/* Scramble the decay times and create the new list mode data */
		if ( !resamptNewDecayTimes( scanTime ) ) {
			/* Error should already have been reported */
			goto FAIL;
		}
		
			/* write header to output randoms added file and close */
		{
			
			/* update the file header */
			if (PhgHdrUpHeader(NULL, &resamptHdrParams, &(resamptOutHistHk.headerHk)) == false) {
				sprintf(resamptErrStr, "Unable to write header to output history file\n'%s'.",
					resamptOutHistName);
				ErStFileError(resamptErrStr);
				goto FAIL;
			}
			
			/* close the output file */
			if (fclose(resamptOutHistHk.histFile) != 0) {
				resamptOutHistHk.histFile = NULL;
				sprintf(resamptErrStr,"Unable to close randoms added file named:\n"
					"'%s'\n"
					" (resamptResampleTime)",
					resamptOutHistName);
				ErStFileError(resamptErrStr);
				goto FAIL;
			}
			
		}
		
		okay = true;
		FAIL:;
		CANCEL:;
	} while (false);
	
	
	return (okay);
}


/*********************************************************************************
*
*	Name:			resamptNewDecayTimes
*
*	Summary:		Resample the decay times and write the altered events to
*					the output list mode file.
*
*	Arguments:
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean resamptNewDecayTimes( double scanTime )
{
	
	#define				PROGRESSminutes		1		/* Time interval between progress displays */
	#define				PROGRESScheckNum	1000000	/* number decays between progress displays */

	Boolean				okay = false;			/* Function return value */
	PhoHFileEventType	eventType;				/* Type of current event */
	PHG_Decay			decay;					/* decay info */
	LbUsFourByte		numBluePhotons;			/* number of blue photons in this decay */
	LbUsFourByte		numPinkPhotons;			/* number of pink photons in this decay */
	PHG_DetectedPhoton	bluePhotons[PHG_MAX_DETECTED_PHOTONS];
												/* the blue photons for this decay */
	PHG_DetectedPhoton	pinkPhotons[PHG_MAX_DETECTED_PHOTONS];
												/* the pink photons for this decay */
	PHG_Decay			nextDecay;				/* Current decay data for input or output */
	PHG_DetectedPhoton	nextPhoton;				/* Current photon data for input or output */
	LbUsFourByte		progressDecayCounter;	/* decay count since last progress report check */
	double				accumSecs;				/* Seconds since last progress display */
	double				accumCPUSecs;			/* Total used CPU time */
	LbUsFourByte		nextAccumSecs;			/* Time between progress displays */
	LbTmTimingType		clockTiming;			/* Timing data for progress messages */
	double				realSecs;				/* Wall clock time elapsed */
	double				cpuSecs;				/* CPU time elapsed */
	LbUsEightByte		bytesReadIn;			/* number bytes read from input history file */
	LbUsFourByte		photonSize;				/* number bytes read in for photon */
	LbUsFourByte		decaySize;				/* number bytes read in for decay */
	time_t				curTime;				/* Current wall clock time */
	char				timeStr[32];			/* String for printing the time */
	
	do {
		
		/* Initialize progress report messages */
		accumSecs = 0.0;
		accumCPUSecs = 0.0;
		nextAccumSecs = PROGRESSminutes * 60;	/* = every xx minutes */
		LbTmStartTiming(&clockTiming);
		bytesReadIn = PHG_HDR_HEADER_SIZE;	/* header already read in above */
		decaySize = sizeof( PHG_Decay ) + 1; /* increment for reading decay (+1 for decay/photon flag) */
		photonSize = sizeof( PHG_DetectedPhoton ) + 1; /* increment for reading photon */

		/* Report start of file creation */
		{
			/* Get current time */
			time(&curTime);
			
			/* Convert time to string format */
			strftime(timeStr, 31, "(%H:%M %e-%b-%y)", localtime(&curTime));
			
			/* Display first progress message */
			LbInPrintf("\n\n********* Resampling of history file %s *********\n",
				timeStr);
		}
		
		/* initialize counter statistics */
		numDecaysRead = 0;
		numDecaysWritten = 0;
			
		if (!resamptCustomFile) {
			/* standard history file processing */
			
			/* Read in first decay */
			{
				eventType = PhoHFileReadEvent(resamptInHistHk.histFile, &decay, 
								&nextPhoton);
				if (eventType != Decay) {
					ErStGeneric("Expected first event to be decay, and it wasn't (resamptNewDecayTimes).");
					goto FAIL;
				}
				bytesReadIn += decaySize;
				numDecaysRead = 1;
				progressDecayCounter = 1;
			}
			
			/* Loop until the list mode file is empty */
			while (eventType == PhoHFileDecayEvent) {
				
				/* clear the current photon counts */
				numBluePhotons = 0;
				numPinkPhotons = 0;
				
				/* reassign decay time */
				decay.decayTime = PhgMathGetDPRandomNumber() * scanTime;
				
				while ( (eventType = PhoHFileReadEvent(resamptInHistHk.histFile, &nextDecay, &nextPhoton)) == 
							PhoHFilePhotonEvent ) {
					
					/* increment the number of input bytes read in */
					bytesReadIn += photonSize;
					
					/* copy this photon to our blue or pink photons */
					if ( PHG_IsBlue(&nextPhoton) ) {
						
						if (numBluePhotons >= PHG_MAX_DETECTED_PHOTONS) {
							ErStGeneric("Number blue photons exceeds maximum (resamptNewDecayTimes).");
							goto FAIL;
						}
						
						bluePhotons[numBluePhotons] = nextPhoton;
						numBluePhotons += 1;
						
					} else {
						
						if (numPinkPhotons >= PHG_MAX_DETECTED_PHOTONS) {
							ErStGeneric("Number pink photons exceeds maximum (resamptNewDecayTimes).");
							goto FAIL;
						}
						
						pinkPhotons[numPinkPhotons] = nextPhoton;
						numPinkPhotons += 1;
						
					}
										
				}
				
				if (eventType == PhoHFileDecayEvent) {
				
					/* increment the decays and bytes read in; generate periodic progress messages */
					{
						
						bytesReadIn += decaySize;
						progressDecayCounter++;
						
						if ( progressDecayCounter >= PROGRESScheckNum ) {
							
							progressDecayCounter = 1;	/* reset counter */
							
							/* check to see if it is time for a progress message */
							if (LbTmStopTiming(&clockTiming, &realSecs, &cpuSecs)) {
								accumSecs += realSecs;
								accumCPUSecs += cpuSecs;
								
								/* Restart the timer */
								LbTmStartTiming(&clockTiming);
								
								if ( (LbUsFourByte) accumSecs >= nextAccumSecs) {
									
									accumSecs = 0.0;	/* easier than incrementing to total */
									
									/* Get current time */
									time(&curTime);
									
									/* Convert time to string format */
									strftime(timeStr, 31, "(%H:%M %e-%b-%y)", localtime(&curTime));
									
									/* Display a progress message */
									LbInPrintf("%d MB resampled, CPU seconds = %3.0f %s.\n",
										(int)(bytesReadIn/1024/1024), 
										accumCPUSecs, timeStr);
								}
							}
						
						}
						
					}
					
					/* write last decay to output list mode file */
					PhoHFileRewriteDetections( &resamptOutHistHk, &(decay),
						bluePhotons, numBluePhotons,
						pinkPhotons, numPinkPhotons);
					
					/* increment the number of decays written out */
					numDecaysWritten++;
					
					/* copy new decay to current decay buffer */
					decay = nextDecay;
					
					/* increment decays read in */
					numDecaysRead++;
					
									
				}
			
			}
			
			/* it is probable there is still a decay to write out, flush it before leaving */
			if ( (numBluePhotons > 0) || (numPinkPhotons > 0) ) {

				/* write this decay to output list mode file */
				PhoHFileRewriteDetections( &resamptOutHistHk, &(decay),
					bluePhotons, numBluePhotons,
					pinkPhotons, numPinkPhotons);
				
				/* increment the number of decays written out and the number of randoms created */
				numDecaysWritten++;
				
			}
			
		}

		/* print out final progress message */
		if (LbTmStopTiming(&clockTiming, &realSecs, &cpuSecs)) {
			accumSecs += realSecs;
			accumCPUSecs += cpuSecs;
			
			/* Get current time */
			time(&curTime);
			
			/* Convert time to string format */
			strftime(timeStr, 31, "(%H:%M %e-%b-%y)", localtime(&curTime));
			
			/* Display a progress message */
			LbInPrintf("%d MB read in,  CPU seconds = %3.0f %s.\n",
				(int)(bytesReadIn/1024/1024), 
				accumCPUSecs, timeStr);
		}
					
		/* Report end of file creation */
		{
			/* Get current time */
			time(&curTime);
			
			/* Convert time to string format */
			strftime(timeStr, 31, "(%H:%M %e-%b-%y)", localtime(&curTime));
			
			/* Display first progress message */
			LbInPrintf("******* Resampling of history file complete %s *******\n",
				timeStr);
		}
				
	} while (false);
	
	okay = true;
	
	
	FAIL:;
	/* NOTE:  In case of failure, files may be left unclosed; therefore,
		it is best to completely exit the program on failure and let the 
		system clean up these files. */
		
	return (okay);
}



#undef RESAMPLE_DECAY_TIME_MAIN
