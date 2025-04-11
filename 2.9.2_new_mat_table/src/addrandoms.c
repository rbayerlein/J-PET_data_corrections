/*********************************************************************************
*                                                                                
*                       Source code developed by the                             
*           Imaging Research Laboratory - University of Washington               
*               (C) Copyright 2003-2013 Department of Radiology          	        
*                           University of Washington                             
*                              All Rights Reserved                               
*                                                                                
*********************************************************************************/

/*********************************************************************************
*
*		Module Name:		addrandoms.c
*		Revision Number:	2.4
*		Date last revised:	7 June 2013
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
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		25 June 2012
*
*			Revision description:	Added parameters to the two detection functions
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		10 January 2012
*
*			Revision description:	Changed form of user functions to pointers
*
*********************************************************************************/

#define ADD_RANDOMS_MAIN	/* Note we are substituting ourselves for phg's main */

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
#include "addrandoms.h"
#include "addrandUsr.h"


#ifdef MPW
	#pragma segment ADD_RANDOMS_MAIN
#endif

/* LOCAL CONSTANTS */
#define ADRAND_IsTestWindow()				LbFgIsSet(PhgOptions, LBFlag0)		/* Will we only check whether the file is sorted */
#define ADRAND_IsDeleteFile()				LbFgIsSet(PhgOptions, LBFlag1)		/* Will we delete the original file after sorting */

#define	ADRAND_NumFlags	2						/* Number of flags defined */

/* LOCAL TYPES */
typedef enum  {Null, Decay, Photon} EventTy;

/* LOCAL GLOBALS */
static Boolean				adrandCancelled;				/* Global cancellation flag */
static Boolean				adrandDoRandoms = false;		/* True randoms should be done */
static Boolean				adrandCustomFile = false;		/* True if custom history file used */
static DetEn_TriplesMethodTy	adrandTriplesMethod = DetEn_DeleteTriples;	/* Method to use for triples */
static char					adrandErrStr[1024];				/* Error string storage */		
static LbUsFourByte			adrandNumToProc;				/* Number of histories to process */
static char					adrandHistName[1024];			/* Name of history file */
static char					adrandDetParamsFileName[256];	/* Name of detector parameters file (from command line) */
static char					adrandRandomsFileName[256];		/* Name of output randoms history file */
static char					adrandHistParamsName[256];		/* Name of history parameters file */
static double				adrandTimingWindow = 0.0;		/* the coincidence timing window in seconds */
static LbUsFourByte			adrandArgIndex;					/* Index through command line arguments */
static timeWindowDetectionsTy	timeWindowDetections;		/* The decays, photons in current time window */
static PhoHFileHdrTy		adrandHdrParams;				/* Input header */
/*static LbHdrHkTy			adrandHeaderHk;*/					/* Hook to original history file header */
static PhoHFileHkTy			adrandHistParamsHk;				/* Hook to custom history file information */
static PhoHFileHkTy			adrandHistHk;					/* Hook to input  sorted history file */
static PhoHFileHkTy			adrandRandomsFileHk;			/* Hook to randoms list mode file */
static LbUsFourByte			numDecaysRead;					/* Number of decays read in so far */
static LbUsFourByte			numDecaysWritten;				/* Number of created random decays written out */
static LbUsFourByte			numDecaysUnchanged;				/* Number of decays written out unchanged */
static LbUsFourByte			numDecaysRandom;				/* Number of random decays created */
static LbUsFourByte			numDecaysLostTriples;			/* Number of decays lost to the triples rule */
static LbUsFourByte			histDecaysPerWindow[ADRAND_MaxTWDecays];	/* Histogram of decays per time window */
static LbUsFourByte			numLostCorrectWindow;			/* Number randoms lost to the by-photon windowing */


/* PROTOTYPES */
Boolean			adrand(int argc, char *argv[]);
Boolean 		adrandInitialize(int argc, char *argv[]);
Boolean			adrandGetParams(void);
void			adrandTerminate(void);
Boolean			adrandTest(char *argv[]);
Boolean			adrandAddRandoms(char *argv[]);
Boolean			adrandCreateRandomsList(void);
Boolean 		adrandProcessTimeWindow(void);


/* FUNCTIONS */
/**********************
*	Name:		adrand
*
*	Purpose:	Execute the program.
*
*	Result:		True unless an error occurs.
***********************/
Boolean adrand(int argc, char *argv[])
{
	Boolean				okay = false;	/* Process Loop */
					/* Error condition results unless this is changed in body */
	LbUsFourByte	numDecays;			/* counter */
	time_t 			curTime;			/* Current time for stamping execution date */
	LbFourByte		argIndex;			/* Index through command line arguments */
	char			executionTimeStr[33];
										/* String for execution date conversion */
	
	/*** NOTE:  Custom history files are not supported and have not been tested. ***/
	
	
	do { /* Process Loop */
	
		/* Perform initialization tasks */
		if (adrandInitialize(argc, argv) == false) {
			sprintf(adrandErrStr, "Add Randoms initialization failed.");
			ErStGeneric(adrandErrStr);
			goto FAIL;
		}
		
		/* Get current time */
		time(&curTime);
		
		/* Print out the command line */
		LbInPrintf("\n\n\nCommand line: ");
		for (argIndex = 0; argIndex < argc; argIndex++)
			LbInPrintf("%s ", argv[argIndex]);
		LbInPrintf("\n");
		
		/* Get current time */
		time(&curTime);
		
		/* Convert time to string format */
		strftime(executionTimeStr, 31, "%H:%M:%S %e-%b-%y", localtime(&curTime));
		
		LbInPrintf("\nExecution of addrandoms occurred at %s\n", executionTimeStr);
		
		/* Print out the main parameters for randoms processing */
		LbInPrintf("\nName of input history file: %s.\n", adrandHistName);
		if (adrandCustomFile) {
			
			LbInPrintf("Using custom history parameters: %s.\n", adrandHistParamsName);
			
		}
		LbInPrintf("Coincidence window = %3.3f nanoseconds.\n",
					DetRunTimeParams[DetCurParams].CoincidenceTimingWindowNS);
		LbInPrintf("Name of output randoms-added history file: %s.\n", adrandHistName);
		
		
		if (ADRAND_IsTestWindow()) {
			/* Just check that no two 'decays' are left in any window */
			/* okay = adrandTest(argv); */
			if (! okay) {
				/* Error should already have been recorded */
				goto FAIL;
			}
		}
		else {
			/* Apply the timing parameters */
			okay = adrandAddRandoms(argv);
			if (! okay) {
				/* Error should already have been recorded */
				goto FAIL;
			}
		}
		
		/* write out report for module */
		LbInPrintf("\n\tNumber of decays read in: %d", numDecaysRead);
		
		LbInPrintf("\n\n\tNumber of decays written out: %d", numDecaysWritten);
		LbInPrintf("\n\tOf the decays written out, %d were written out unchanged.", numDecaysUnchanged);
		LbInPrintf("\n\tOf the decays written out, %d were randoms created by addrandoms.", numDecaysRandom);
		
		LbInPrintf("\n\n\tNumber of decays lost as triples: %d", numDecaysLostTriples);
		LbInPrintf("\n\n\tNumber of decays lost to correct windowing: %d", numLostCorrectWindow);
		
		/* write out histogram of decays found per time window */
		LbInPrintf("\n\n\tHistogram of the number of decays found per time window:");
		LbInPrintf("\n\n\tTime windows with 1 decay  =\t\t%d", histDecaysPerWindow[0]);
		for ( numDecays = 2; numDecays < ADRAND_MaxTWDecays; numDecays++ ) {
			
			LbInPrintf("\n\tTime windows with %d decays =\t\t%d", numDecays, histDecaysPerWindow[numDecays-1]);
			
		}
		LbInPrintf("\n\tTime windows with %d or more decays =\t%d\n\n", numDecays, histDecaysPerWindow[numDecays-1]);

		/* Get current time */
		time(&curTime);
		
		/* Convert time to string format */
		strftime(executionTimeStr, 31, "%H:%M:%S %e-%b-%y", localtime(&curTime));
			
		/* Print out our final message */
		LbInPrintf("\n\n\nExecution of addrandoms finished at %s\n", executionTimeStr);
		
		okay = true;
		FAIL:;
		CANCEL:;
		
	} while (false);
	
	
	/* Terminate PHG modules */
	adrandTerminate();
	
	/* Handle error situation if one exists */
	if (!okay && adrandCancelled) {
		ErHandle("User cancelled addrandoms", false);
		okay = true;
	}
	
	/* Quit the program */
	return (okay);
}


/**********************
*	Name:		adrandInitialize
*
*	Purpose:	Perform initialization tasks.
*
*	Result:		True unless an error occurs.
***********************/
Boolean adrandInitialize(int argc, char *argv[])
{
	Boolean				okay = false;				/* Process Loop */
	char				*knownOptions[] = {"tr"};
	char				optArgs[ADRAND_NumFlags][LBEnMxArgLen];
	LbUsFourByte		optArgFlags = (LBFlag0);
	char				fileName[1024];		/* Name of param file */
	LbFourByte			randSeed;
	LbUsFourByte		numDecays;			/* counter */
	
	
	do { /* Process Loop */
		
		/* Get our options */
		if (!LbEnGetOptions(argc, argv, knownOptions,
				&PhgOptions, optArgs, optArgFlags, &adrandArgIndex)) {
	
			break;
		}
				
		/* See if they supplied a command line argument */
		if ((adrandArgIndex != 0) && (argv[adrandArgIndex] != 0)) {
		
			/* Get first param file and save number of param files to process */
			strcpy(adrandDetParamsFileName,argv[adrandArgIndex]);
			adrandNumToProc = (argc - adrandArgIndex);
		}
		else {
			/* Ask for the file name */
			LbInAsk("Enter name of detector parameters file", 0, false,
					&adrandCancelled, 0, 0, 0, 0,
					fileName);
		
			/* Bolt if we canceled */
			if (adrandCancelled) {
				ErStCancel("User cancelled.");
				goto CANCEL;
			}
			adrandNumToProc = 1;
			strcpy(adrandDetParamsFileName, fileName);
		}
		
		/* Get our run-time parameters */
		if (!adrandGetParams())
			break;
		
		if (addrandUsrInitializeFPtr) {
			(*addrandUsrInitializeFPtr)(&DetRunTimeParams[DetCurParams]);
		}
		
		/* initialize counters */
		numDecaysRead = 0;
		numDecaysWritten = 0;
		numDecaysUnchanged = 0;
		numDecaysRandom = 0;
		numDecaysLostTriples = 0;
		numLostCorrectWindow = 0;
		for ( numDecays = 0; numDecays < ADRAND_MaxTWDecays; numDecays++ ) {
			
			histDecaysPerWindow[numDecays] = 0;
			
		}
		
		/* Initialize the math library */
		/* NOTE:  Random numbers not known to be used in this program */
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


/*********************************************************************************
*
*	Name:			adrandGetParams
*
*	Summary:		Read in the addrandom parameters.
*
*	Arguments:		None.
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean adrandGetParams()	
{
	LbUsFourByte			loopV =0;
	Boolean					okay = false;					/* Process flag */
	
	
	do { /* Process Loop */
		
		/* Clear out path name variables in case of error */
		PhgRunTimeParams.PhgRandomSeed = 0.0;
 		memset(PhgRunTimeParams.PhgSubObjActivityIndexFilePath, '\0', PATH_LENGTH);
		memset(PhgRunTimeParams.PhgSubObjActivityTableFilePath, '\0', PATH_LENGTH);
 		memset(PhgRunTimeParams.PhgSubObjActIndexTransFilePath, '\0', PATH_LENGTH);
 		memset(PhgRunTimeParams.PhgSubObjActImgFilePath, '\0', PATH_LENGTH);
 		memset(PhgRunTimeParams.PhgSubObjAttenIndexFilePath, '\0', PATH_LENGTH);
 		memset(PhgRunTimeParams.PhgSubObjAttenTableFilePath, '\0', PATH_LENGTH);
	 	memset(PhgRunTimeParams.PhgSubObjAttIndexTransFilePath, '\0', PATH_LENGTH);
  		memset(PhgRunTimeParams.PhgSubObjAttImgFilePath, '\0', PATH_LENGTH);
  		memset(PhgRunTimeParams.PhgSubObjCohScatTableFilePath, '\0', PATH_LENGTH);
		memset(PhgRunTimeParams.PhgProdTblInputTableFilePath, '\0', PATH_LENGTH);
 		memset(PhgRunTimeParams.PhgProdTblOutputTableFilePath, '\0', PATH_LENGTH);
	    memset(PhgRunTimeParams.PhgPhoTrkForcedDetectionFilePath, '\0', PATH_LENGTH);
	    memset(PhgRunTimeParams.PhgPhoTrkForcedDetectionFilePath, '\0', PATH_LENGTH);
		PhgRunTimeParams.PhgIsBinOnTheFly = false;
		PhgNumBinParams = 0;
		PhgNumTomoFiles = 0;
		ColNumParams = 0;
		DetNumParams = 0;
		
		/* Clear out array of possible tomograph path names */
		for (loopV = 0; loopV < PHG_MAX_PARAM_FILES; loopV++) {
		    memset(PhgRunTimeParams.PhgCollimatorParamsFilePath[loopV], '\0', PATH_LENGTH);
		    memset(PhgRunTimeParams.PhgDetectorParamsFilePath[loopV], '\0', PATH_LENGTH);
		    memset(PhgRunTimeParams.PhgBinParamsFilePath[loopV], '\0', PATH_LENGTH);
			memset(PhgRunTimeParams.PhgTomoParamsFilePath[loopV], '\0', PATH_LENGTH);
		}
		
		/* At first, it was assumed that all parameters would be present since all were included
			in the sample files.  As they are added, I am making sure that new parameters get set
			to default values.  Sometimes this is done in module initialization routines, and sometimes
			it's done here.
		*/
		PhgRunTimeParams.PhgIsModelCoherent = false;
		PhgRunTimeParams.PhgIsModelCoherentInObj = false;
		PhgRunTimeParams.PhgIsModelCoherentInTomo = false;
		PhgRunTimeParams.PhgIsModelPolarization = false;
		PhgRunTimeParams.PhgNuclide.isotope = PhgEn_IsotopType_NULL;
		EmisListIsotopeDataFilePath[0] = '\0';
		
		/* assign the run parameters field used to open the detector parameters file */
		DetCurParams = 0;	/* this may need to be changed if addrandoms ever gets called in main phg loop */
		strcpy(PhgRunTimeParams.PhgDetectorParamsFilePath[DetCurParams], adrandDetParamsFileName);
		
		/* Initialize detector parameters */
		DetRunTimeParams[DetCurParams].DetectorType = DetEn_DetType_NULL;
		
		/* Read the parameters */
		if ( !DetGetRunTimeParams() ) break;
		
		/* Check that a detector history file name is given */
		if ( strlen(DetRunTimeParams[DetCurParams].DetHistoryFilePath) == 0 ) {
				
				ErStGeneric("No detector history file supplied in detector parameters file.");
				break;
				
		}
		
		/* Save the randoms-processing related parameters to local variables. */
		{
			/* input history file name */
			strcpy(adrandHistName, DetRunTimeParams[DetCurParams].DetHistoryFilePath);
			/* custom history params */
			strcpy(adrandHistParamsName, DetRunTimeParams[DetCurParams].DetHistoryParamsFilePath);
			/* should randoms be done? */
			adrandDoRandoms = DetRunTimeParams[DetCurParams].DoRandomsProcessing;
			/* coincidence window - convert from nanoseconds to seconds as all times in tracking are seconds */
			adrandTimingWindow = DetRunTimeParams[DetCurParams].CoincidenceTimingWindowNS * 1.0E-9;
			/*	how do we handle triples? */
			adrandTriplesMethod = DetRunTimeParams[DetCurParams].TriplesMethod;
			/* output randoms-added list data */
			strcpy( adrandRandomsFileName, DetRunTimeParams[DetCurParams].DetRandomsHistoryFilePath);
		}
		
		okay = true;
		FAIL:;
	} while (false);
	
	return (okay);
}


/**********************
*	Name:		adrandTerminate
*
*	Purpose:	Terminate all of the managers.
*
*	Result:		None.
***********************/
void adrandTerminate()
{
	ProdTblProdTblInfoTy	hstPrdTblInfo;		/* Info for initializing productivity table */
	
	/* call the user-add-on termination function */
	if (addrandUsrTerminateFPtr) {
		(*addrandUsrTerminateFPtr)(&DetRunTimeParams[DetCurParams]);
	}
	
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


/**********************
*	Name:		adrandTest
*
*	Purpose:	Test a standard history file to see if it is sorted.
*
*	Result:		True unless an error occurs.
***********************/
Boolean adrandTest(char *argv[])
{
	Boolean				okay = false;				/* Process flag (function return) */
	char**				dummyPtr;					/* Removes compiler warning */
	
	/*LbUsFourByte		curFileIndex;*/				/* Current file index */
	/*FILE				*historyFile;*/				/* The history file we are going to process */
	/*LbHdrHkTy			histHeaderHk;*/				/* Hook to history file header */
	/*EventTy				eventType;*/					/* Type of current event */
	
	dummyPtr = argv;		/* Avoid unused parameter compiler warning */
	
	
		
	
	return (okay);
}


/*********************************************************************************
*
*	Name:			adrandAddRandoms
*
*	Summary:		Combine singles to create randoms in a history file.
*
*	Arguments:
*		char			*argv[]					- Command line arguments.
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean adrandAddRandoms(char *argv[])
{
	Boolean				okay = false;		/* Process flag (function return) */
	LbUsFourByte		curFileIndex;		/* Current file index */
	FILE				*historyFile;		/* The detector history file to add randoms to */
	
	
	do { /* Process Loop */
		
		/* Process the files */
		for (curFileIndex = 1; curFileIndex <= adrandNumToProc; curFileIndex++){
			
			/* open files and verify parameters */
			{
			
				/* Validate the name of the sorted file (supplied by params file) */
				if (strlen(adrandHistName) == 0) {
					sprintf(adrandErrStr, "No randoms-added input file name supplied.");
					ErStGeneric(adrandErrStr);
					goto FAIL;
				}
				
				/* Validate the name of the randoms-added file (supplied by params file) */
				if (strlen(adrandRandomsFileName) == 0) {
					sprintf(adrandErrStr, "No randoms-added output file name supplied.");
					ErStGeneric(adrandErrStr);
					goto FAIL;
				}
				
				/* Check that randoms processing is selected */
				if ( !adrandDoRandoms ) {
					sprintf(adrandErrStr, "(adrandAddRandoms) Randoms processing not specified.");
					ErStGeneric(adrandErrStr);
					goto FAIL;
				}
				
				/* Check a positive coincidence timing window is selected */
				if ( adrandTimingWindow <= 0.0 ) {
					sprintf(adrandErrStr, "(adrandAddRandoms) Non-positive coincidence timing window not allowed.");
					ErStGeneric(adrandErrStr);
					goto FAIL;
				}
				
				/* Open input history file */
				if ((historyFile = LbFlFileOpen(adrandHistName, "rb")) == 0) {
					sprintf(adrandErrStr, "Unable to open history file\n'%s'.",
						adrandHistName);
					ErStFileError(adrandErrStr);
					goto FAIL;
				}
				
				adrandHistHk.histFile = historyFile;
				
				/* Init header hook, read in input header and copy to output header and save it */
				/* This will overwrite the DetRunTimeParams previously read in with those actually
				used to create the input data */
				if (PhgHdrGtParams(adrandHistHk.histFile, &adrandHdrParams, &adrandHistHk.headerHk) == false){
					sprintf(adrandErrStr, "Unable to read history file header\n'%s'.",
						adrandHistName);
					ErStFileError(adrandErrStr);
					goto FAIL;
				}
				
				/* Modify the parameters that are specific to addrandoms */
				adrandHdrParams.H.DetRunTimeParams.DoRandomsProcessing = adrandDoRandoms;
				/* convert timing window for header back to nanoseconds */
				adrandHdrParams.H.DetRunTimeParams.CoincidenceTimingWindowNS = adrandTimingWindow * 1.0E9;
				adrandHdrParams.H.DetRunTimeParams.TriplesMethod = adrandTriplesMethod;
				strcpy( adrandHdrParams.H.DetRunTimeParams.DetRandomsHistoryFilePath, adrandRandomsFileName );
				
				/* Check that all the importance sampling features are off */
				{
					
					if ( adrandHdrParams.H.PhgRunTimeParams.PhgIsCalcEventsToSimulate != true ) {
						sprintf(adrandErrStr, "Data for randoms processing must be generated with num_to_simulate = 0.");
						ErStFileError(adrandErrStr);
						goto FAIL;
					}
						
					if ( adrandHdrParams.H.PhgRunTimeParams.PhgIsForcedDetection ) {
						sprintf(adrandErrStr, "Randoms processing is incompatible with Forced Detection.");
						ErStFileError(adrandErrStr);
						goto FAIL;
					}
						
					if ( adrandHdrParams.H.PhgRunTimeParams.PhgIsStratification ) {
						sprintf(adrandErrStr, "Randoms processing is incompatible with Stratification.");
						ErStFileError(adrandErrStr);
						goto FAIL;
					}
						
					if ( adrandHdrParams.H.PhgRunTimeParams.PhgIsForcedNonAbsorbtion ) {
						sprintf(adrandErrStr, "Randoms processing is incompatible with Forced Non-Absorption.");
						ErStFileError(adrandErrStr);
						goto FAIL;
					}
					
					if ( adrandHdrParams.H.DetRunTimeParams.DoForcedInteraction ) {
						sprintf(adrandErrStr, "Randoms processing is incompatible with Forced Interaction in detector.");
						ErStFileError(adrandErrStr);
						goto FAIL;
					}
					
				}

				/* Check that this is a simulate_PET_coincidences_plus_singles scan */
				if ( adrandHdrParams.H.PhgRunTimeParams.PhgIsPETCoincPlusSingles != true ) {
					sprintf(adrandErrStr, "Randoms processing requires simulate_PET_coincidences_plus_singles.");
					ErStFileError(adrandErrStr);
					goto FAIL;
				}
				
				/* Check that this file has been time sorted */
				if ( adrandHdrParams.H.isTimeSorted != true ) {
					sprintf(adrandErrStr, "Randoms processing requires a time-sorted history file.");
					ErStFileError(adrandErrStr);
					goto FAIL;
				}
				
				/* Open and initialize the custom history file parameters, if any */
				{
					
					if (adrandCustomFile) {
						/* Setup header hook */
						
						adrandHistParamsHk.doCustom = true;
						
						/* Do custom parameter initialization  */
						if (PhoHFileGetRunTimeParams(adrandHistParamsName, &(adrandHistParamsHk.customParams)) == false) {
							sprintf(adrandErrStr,"Unable to get custom parameters for history file named '%s'",
								adrandHistParamsName);
							ErAlert(adrandErrStr, false);
							goto FAIL;
						}
						
						if ( adrandHistParamsHk.customParams.doDecayTime == false ) {
							sprintf(adrandErrStr, "Randoms processing of custom list mode requires doDecayTime = true.");
							ErStFileError(adrandErrStr);
							goto FAIL;
						}
						
						if ( adrandHistParamsHk.customParams.doDecayType == false ) {
							sprintf(adrandErrStr, "Randoms processing of custom list mode requires doDecayType = true.");
							ErStFileError(adrandErrStr);
							goto FAIL;
						}
						
						if ( adrandHistParamsHk.customParams.doTravelDistance == false ) {
							sprintf(adrandErrStr, "Randoms processing of custom list mode requires doTravelDistance = true.");
							ErStFileError(adrandErrStr);
							goto FAIL;
						}
						
						strcpy(adrandHistParamsHk.customParamsName, adrandHistParamsName);
						
						adrandHistParamsHk.bluesReceived = 0;
						adrandHistParamsHk.bluesAccepted = 0;
						adrandHistParamsHk.pinksReceived = 0;
						adrandHistParamsHk.pinksAccepted = 0;
						adrandHistParamsHk.histFile = historyFile;
					}
					
				}
				
				/* Open the randoms-added output file */
				if (PhoHFileCreate(adrandRandomsFileName, "", 
						adrandHdrParams.H.HdrKind, &adrandRandomsFileHk) == false) {
					sprintf(adrandErrStr,"Unable to create output randoms file named:\n"
						"'%s'\n"
						" (adrandAddRandoms)",
						adrandRandomsFileName);
					ErStFileError(adrandErrStr);
					goto FAIL;
				}
				
			}
			
			/* Create the list mode data with randoms */
			if ( !adrandCreateRandomsList() ) {
				/* Error should already have been reported */
				goto FAIL;
			}
			
			/* write header to output randoms added file and close */
			{
				
				/* mark the file as having randoms added */
				adrandHdrParams.H.isRandomsAdded = true;
				if (LbHdrStElem(&(adrandRandomsFileHk.headerHk), HDR_HISTORY_FILE_IS_RANDOMS_ADDED_ID,
						sizeof(adrandHdrParams.H.isRandomsAdded),
						(void *)&(adrandHdrParams.H.isRandomsAdded)) == false){
					
					sprintf(adrandErrStr,"Unable to set header 'randoms added indicator' parameter.");
					ErStFileError(adrandErrStr);
					goto FAIL;
				}
				
				/* update the file header */
				if (PhgHdrUpHeader(NULL, &adrandHdrParams, &(adrandRandomsFileHk.headerHk)) == false) {
					sprintf(adrandErrStr, "Unable to write header to randoms added history file\n'%s'.",
						adrandRandomsFileName);
					ErStFileError(adrandErrStr);
					goto FAIL;
				}
				
				/* close the output file */
				if (fclose(adrandRandomsFileHk.histFile) != 0) {
					adrandRandomsFileHk.histFile = NULL;
					sprintf(adrandErrStr,"Unable to close randoms added file named:\n"
						"'%s'\n"
						" (adrandAddRandoms)",
						adrandRandomsFileName);
					ErStFileError(adrandErrStr);
					goto FAIL;
				}
				
			}
			
			
			/* Delete the original file, if requested */
			if (ADRAND_IsDeleteFile()) {
				if (remove(adrandHistName) != 0) {
					sprintf(adrandErrStr, "Unable to delete history file\n'%s'.",
						adrandHistName);
					ErStFileError(adrandErrStr);
					/* But keep going */
				}
			}
			
			
			/* Open next parameters file */
			if (curFileIndex < adrandNumToProc) {
				
				strcpy(PhgRunTimeParams.PhgParamFilePath,argv[curFileIndex+adrandArgIndex]);
				
				/* Get our run-time parameters */
				if (!adrandGetParams()) {
					goto FAIL;
				}
				
				/* Let them know what is going on */
				printf("\n***********************\n");
				printf("About to process parameter file '%s'\n", PhgRunTimeParams.PhgParamFilePath);
				printf("\n***********************\n");

				if (strlen(adrandHistParamsName) != 0) {
					/* Custom history file */
					adrandCustomFile = true;
				}
				else {
					adrandCustomFile = false;
				}
				
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
*	Name:			adrandCreateRandomsList
*
*	Summary:		Read listmode data in one time window at a time and
*					call processing/writing out function.
*
*	Arguments:
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean adrandCreateRandomsList( )
{
	#define				PROGRESSminutes		1		/* Time interval between progress displays */
	#define				PROGRESScheckNum	1000000	/* number decays between progress displays */
	
	Boolean				okay = false;				/* Function return value */
	double				photonDetectionTime;		/* time between beginning of scan and detection */
	PHG_Decay			nextDecay;					/* Current decay data for input or output */
	PHG_DetectedPhoton	nextPhoton;					/* Current photon data for input or output */
	PhoHFileEventType	eventType;					/* Type of current event */
	LbUsFourByte		curDecayNum;				/* current decay counter */
	LbUsFourByte		curPhotonNum;				/* current photon number */
	LbUsFourByte		progressDecayCounter;		/* decay count since last progress report check */
	double				accumSecs;					/* Seconds since last progress display */
	double				accumCPUSecs;				/* Total used CPU time */
	LbUsFourByte		nextAccumSecs;				/* Time between progress displays */
	LbTmTimingType		clockTiming;				/* Timing data for progress messages */
	double				realSecs;					/* Wall clock time elapsed */
	double				cpuSecs;					/* CPU time elapsed */
	LbUsEightByte		bytesReadIn;				/* number bytes read from input history file */
	LbUsFourByte		photonSize;					/* number bytes read in for photon */
	LbUsFourByte		decaySize;					/* number bytes read in for decay */
	time_t				curTime;					/* Current wall clock time */
	char				timeStr[32];				/* String for printing the time */

	
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
		LbInPrintf("\n\n******* Beginning creation of randoms-added data %s *******\n",
			timeStr);
	}
	
	/* initialize time window */
	timeWindowDetections.numDecays = 0;
	timeWindowDetections.lastDetectionTime = 0.0;
	for (curDecayNum = 0; curDecayNum < ADRAND_MaxTWDecays; curDecayNum++) {
		timeWindowDetections.decays[curDecayNum].numBluePhotons = 0;
		timeWindowDetections.decays[curDecayNum].numPinkPhotons = 0;
	}
			

	if (! adrandCustomFile) {
		/* standard history file processing */
		
		/* Fill the first time window decay buffer */
		{
			eventType = PhoHFileReadEvent(adrandHistHk.histFile, &(timeWindowDetections.decays[0].decay), 
							&nextPhoton);
			if (eventType != Decay) {
				ErStGeneric("Expected first event to be decay, and it wasn't.");
				goto FAIL;
			}
			bytesReadIn += decaySize;
			timeWindowDetections.numDecays++;
			numDecaysRead = 1;
			progressDecayCounter = 1;
			curDecayNum = 0;
		}
		
		/* Loop until the list mode file is empty */
		while (eventType == PhoHFileDecayEvent) {
			
			while ( (eventType = PhoHFileReadEvent(adrandHistHk.histFile, &nextDecay, &nextPhoton)) == 
						PhoHFilePhotonEvent ) {
				
				/* increment the number of input bytes read in */
				bytesReadIn += photonSize;
				
				/* check if this photon extends our time window */
				photonDetectionTime = timeWindowDetections.decays[curDecayNum].decay.decayTime +
									nextPhoton.time_since_creation;
									
				if ( photonDetectionTime > timeWindowDetections.lastDetectionTime ) {
					timeWindowDetections.lastDetectionTime = photonDetectionTime;
				}
				
				if ( PHG_IsBlue(&nextPhoton) ) {
					
					curPhotonNum = timeWindowDetections.decays[curDecayNum].numBluePhotons;
					timeWindowDetections.decays[curDecayNum].bluePhotons[curPhotonNum] = nextPhoton;
					timeWindowDetections.decays[curDecayNum].numBluePhotons += 1;
					
				} else {
					
					curPhotonNum = timeWindowDetections.decays[curDecayNum].numPinkPhotons;
					timeWindowDetections.decays[curDecayNum].pinkPhotons[curPhotonNum] = nextPhoton;
					timeWindowDetections.decays[curDecayNum].numPinkPhotons += 1;
					
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
								LbInPrintf("%d MB read in,  CPU seconds = %3.0f %s.\n",
									(int)(bytesReadIn/1024/1024), 
									accumCPUSecs, timeStr);
							}
						}
					
					}
					
				}
				
				if ( nextDecay.decayTime >= (timeWindowDetections.lastDetectionTime + adrandTimingWindow) ) {
					
					
					/* process last time window */
					adrandProcessTimeWindow();
					
					/* start new time window with the last decay */
					timeWindowDetections.decays[0].decay = nextDecay;
					timeWindowDetections.lastDetectionTime = nextDecay.decayTime;
					timeWindowDetections.decays[0].numBluePhotons = 0;
					timeWindowDetections.decays[0].numPinkPhotons = 0;				
					timeWindowDetections.numDecays = 1;
					numDecaysRead++;
					curDecayNum = 0;
				
				} else if ( timeWindowDetections.numDecays < (ADRAND_MaxTWDecays - 1) ) {
					
					/* add decay to current time window */
					curDecayNum++;
					timeWindowDetections.decays[curDecayNum].decay = nextDecay;
					if ( nextDecay.decayTime > timeWindowDetections.lastDetectionTime )
						timeWindowDetections.lastDetectionTime = nextDecay.decayTime;
					timeWindowDetections.decays[curDecayNum].numBluePhotons = 0;
					timeWindowDetections.decays[curDecayNum].numPinkPhotons = 0;				
					timeWindowDetections.numDecays += 1;
					numDecaysRead++;
					
				} else {
					
					/* place this decay as last decay in current time window, 
					 bouncing previous last decay.  This is okay as long as one doesn't
					 use a window with so many decays in it, e.g. currently all decays
					 in the window will be dumped anyway as we reject triples. */
					timeWindowDetections.decays[curDecayNum].decay = nextDecay;
					if ( nextDecay.decayTime > timeWindowDetections.lastDetectionTime )
						timeWindowDetections.lastDetectionTime = nextDecay.decayTime;
					timeWindowDetections.decays[curDecayNum].numBluePhotons = 0;
					timeWindowDetections.decays[curDecayNum].numPinkPhotons = 0;				
					timeWindowDetections.numDecays = ADRAND_MaxTWDecays;
					numDecaysRead++;
					
				}
								
			}
		
		}
		
		/* it is probable there are still decay(s) in the time window, flush them before leaving */
		if ( timeWindowDetections.numDecays > 0 ) {
			adrandProcessTimeWindow();
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
			LbInPrintf("******* Creation of randoms-added data complete %s *******\n",
				timeStr);
		}
	}
	
	okay = true;
	
	FAIL:;
	/* NOTE:  In case of failure, files may be left unclosed; therefore,
		it is best to completely exit the program on failure and let the 
		system clean up these files. */
		
	return (okay);
}


/*********************************************************************************
*
*	Name:			adrandProcessTimeWindow
*
*	Summary:		Write out prompt and random events from current time window,
*					purge the time window.
*
*	Arguments:
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean adrandProcessTimeWindow( )
{
	
	Boolean				okay = false;			/* Function return value */
	LbUsFourByte		totalPhotons = 0;		/* total photons in this time window */
	LbUsFourByte		decayNum;				/* loop counter for decays */
	
	do {
		
		/* increment the histogram of decays per time window */
		histDecaysPerWindow[timeWindowDetections.numDecays - 1] += 1;
		
		/* give user a chance to modify/process time window decays and change counters 
			before standard processing */
		if (addrandUsrModDecays1FPtr) {
			(*addrandUsrModDecays1FPtr)(&timeWindowDetections, 
				&adrandRandomsFileHk, &numDecaysWritten, &numDecaysUnchanged, &numDecaysRandom, 
				&numLostCorrectWindow, &numDecaysLostTriples);
		}
		
		/* process the time window depending on how many decays and photons there are in it */
		if ( timeWindowDetections.numDecays == 1 ) {
			
			if ( (timeWindowDetections.decays[0].numBluePhotons >= 1) &&
					(timeWindowDetections.decays[0].numPinkPhotons >= 1) ) {
				
				if ((!addrandUsrModDecays2FPtr) || 
					((*addrandUsrModDecays2FPtr)(&timeWindowDetections, &adrandRandomsFileHk, 
						&numDecaysWritten, &numDecaysUnchanged, 
						&numDecaysRandom, &numLostCorrectWindow, &numDecaysLostTriples))) {
					
					/* write this decay to output list mode file */
					PhoHFileRewriteDetections( &adrandRandomsFileHk, &(timeWindowDetections.decays[0].decay),
						timeWindowDetections.decays[0].bluePhotons, timeWindowDetections.decays[0].numBluePhotons,
						timeWindowDetections.decays[0].pinkPhotons, timeWindowDetections.decays[0].numPinkPhotons);
					
					/* increment the number of decays written out, and the number unchnged written out */
					numDecaysWritten++;
					numDecaysUnchanged++;
					
				}
				
			}
			
		} else if ( timeWindowDetections.numDecays > 1 ) {
			
			for ( decayNum = 0; decayNum < timeWindowDetections.numDecays; decayNum++ ) {
				
				totalPhotons += timeWindowDetections.decays[decayNum].numBluePhotons;
				totalPhotons += timeWindowDetections.decays[decayNum].numPinkPhotons;
				
			}
			
			/* if this time window has exactly two photons (from different decays) we will create
			 a random from them.  Otherwise we will discard all the photons/decays in the window.
			 This is where one would insert a scheme for handling triples and other multiples.
			 Currently they are discarded */
			decayNum = 0;
			if ( totalPhotons == 2 ) {
				
				totalPhotons = 0;
				while (totalPhotons < 2) {
					
					if ( timeWindowDetections.decays[decayNum].numBluePhotons == 1 ) {
						
						if ( totalPhotons == 0 ) {
							
							timeWindowDetections.decays[0].bluePhotons[0] =
									timeWindowDetections.decays[decayNum].bluePhotons[0];
							timeWindowDetections.decays[0].bluePhotons[0].flags =
									PHGFg_PhotonBlue;
							timeWindowDetections.decays[0].numBluePhotons = 1;
							totalPhotons++;
							
						} else {
							
							timeWindowDetections.decays[0].pinkPhotons[0] =
									timeWindowDetections.decays[decayNum].bluePhotons[0];
							timeWindowDetections.decays[0].pinkPhotons[0].flags = 0;
							timeWindowDetections.decays[0].numPinkPhotons = 1;
							/* we must also adjust the photon time-since-creation
							 to reflect the difference in decay times */
							timeWindowDetections.decays[0].pinkPhotons[0].time_since_creation +=
									(timeWindowDetections.decays[1].decay.decayTime -
									timeWindowDetections.decays[0].decay.decayTime);
							totalPhotons++;
							
						}
						
					} else if ( timeWindowDetections.decays[decayNum].numPinkPhotons == 1 ) {
						
						if ( totalPhotons == 0 ) {
							
							timeWindowDetections.decays[0].bluePhotons[0] =
									timeWindowDetections.decays[decayNum].pinkPhotons[0];
							timeWindowDetections.decays[0].bluePhotons[0].flags =
									PHGFg_PhotonBlue;
							timeWindowDetections.decays[0].numBluePhotons = 1;
							totalPhotons++;
							
						} else {
							
							timeWindowDetections.decays[0].pinkPhotons[0] =
									timeWindowDetections.decays[decayNum].pinkPhotons[0];
							timeWindowDetections.decays[0].pinkPhotons[0].flags = 0;
							timeWindowDetections.decays[0].numPinkPhotons = 1;
							/* we must also adjust the photon time-since-creation
							 to reflect the difference in decay times */
							timeWindowDetections.decays[0].pinkPhotons[0].time_since_creation +=
									(timeWindowDetections.decays[1].decay.decayTime -
									timeWindowDetections.decays[0].decay.decayTime);
							totalPhotons++;
							
						}						
						
					}
					
					decayNum++;
					
					/* if we have been through all the decays we should have created a random event.
					 If not, we abort with error message. */
					if ( (decayNum == timeWindowDetections.numDecays) && (totalPhotons < 2) ) {
						
						ErAbort("Unable to create random event (adrandProcessTimeWindow).");
						
					}
				}
				
				/* set the decay type to random event */
				timeWindowDetections.decays[0].decay.decayType = PhgEn_PETRandom;
				
				/* check to see if the two photons in the random do pass within time window -
				 to this point we have only compared the decay time of the second photon...*/
				if ( fabs( timeWindowDetections.decays[0].pinkPhotons[0].time_since_creation -
							timeWindowDetections.decays[0].bluePhotons[0].time_since_creation )
							<= adrandTimingWindow ) {
					
					if ((!addrandUsrModDecays2FPtr) || 
						((*addrandUsrModDecays2FPtr)(&timeWindowDetections, &adrandRandomsFileHk, 
							&numDecaysWritten, &numDecaysUnchanged, 
							&numDecaysRandom, &numLostCorrectWindow, &numDecaysLostTriples))) {
						
						/* write this decay to output list mode file */
						PhoHFileRewriteDetections( &adrandRandomsFileHk, &(timeWindowDetections.decays[0].decay),
							timeWindowDetections.decays[0].bluePhotons, timeWindowDetections.decays[0].numBluePhotons,
							timeWindowDetections.decays[0].pinkPhotons, timeWindowDetections.decays[0].numPinkPhotons);
						
						/* increment the number of decays written out and the number of randoms created */
						numDecaysWritten++;
						numDecaysRandom++;
						
					}
					
				} else {
					
					/* keep track of how many randoms are lost to this second windowing */
					numLostCorrectWindow++;
					
				}
				
			} else if ( totalPhotons > 2 ) {
				
				/* Triples - keep track how many decays are lost to triples */
				numDecaysLostTriples += timeWindowDetections.numDecays;
				
			}
			
		}
		
	} while (false);
	
	okay = true;
	
	
	FAIL:;
	/* NOTE:  In case of failure, files may be left unclosed; therefore,
		it is best to completely exit the program on failure and let the 
		system clean up these files. */
		
	return (okay);
}


#undef ADD_RANDOMS_MAIN

