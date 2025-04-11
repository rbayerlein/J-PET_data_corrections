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
*		Module Name:		timesort.c
*		Revision Number:	1.92
*		Date last revised:	23 July 2013
*		Programmer:			Steven Gillispie
*		Date Originated:	12 December 2003
*
*		Module Overview:	Sort a phg history file according to decay time.
*
*		References:			Knuth (1973) vol III, on tape sorting, may be useful.
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
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		12 September 2012
*
*			Revision description:	Corrected sort indicators to be signed
*
*********************************************************************************/

#define TIME_SORT_MAIN	/* Note we are substituting ourselves for phg's main */

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
#include "LbSort.h"
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
	#pragma segment TIME_SORT_MAIN
#endif

/* Debugging option to keep all intermediate sort files */
/*#define TMSORT_KEEP_TEMP_FILES*/

/* Debugging option to study photon count distribution */
/*#define TMSORT_PHO_DIST*/

/* LOCAL CONSTANTS */
#define TMSORT_IsUsePHGHistory()			LbFgIsSet(PhgOptions, LBFlag0)		/* Will we use the PHG history file */
#define TMSORT_IsUseColHistory()			LbFgIsSet(PhgOptions, LBFlag1)		/* Will we use the Collimator history file */
#define TMSORT_IsUseDetHistory()			LbFgIsSet(PhgOptions, LBFlag2)		/* Will we use the Detector history file */
#define TMSORT_IsTestSorted()				LbFgIsSet(PhgOptions, LBFlag3)		/* Will we only check whether the file is sorted */
#define TMSORT_IsDeleteFile()				LbFgIsSet(PhgOptions, LBFlag4)		/* Will we delete the original file after sorting */

#define	TMSORT_NumFlags	5							/* Number of flags defined */

#define		TMSORT_MINKEY	-1.0				/* Minimum sort key; also, must never actually occur */
#define		TMSORT_MAXKEY	LBDOUBLE_MAX		/* Maximum sort key; also, must never actually occur */

#define		TMSORT_TEMP_COUNT_OFFSET	0		/* Offset into decay data to temporarily store photon count */

#define		TMSORT_DIST_CT		4				/* Number of read slot distributions to track */


/* LOCAL TYPES */
typedef enum  {Null, Decay, Photon} EventTy;

typedef struct  {			/* Temporary decay data marker */
	LbUsFourByte		numPhotons;			/* Number of photons in the decay data */
	PHG_Decay 			*decayDataPtr;		/* Pointer into the buffer to the decay data */
} DecaySortData;
struct DecaySortBlock  {	/* DecaySortData block marker */
	DecaySortData			*thisDataBlock;	/* Pointer to a dynamically allocated block of DecaySortData */
	struct DecaySortBlock 	*nextSortBlock;	/* Pointer to the next block marker in the list */
};
typedef struct DecaySortBlock		DecaySortBlock;

typedef double				TmSortKeyType;			/* The timesort key type */
typedef TmSortKeyType		*TmSortKeyPtr;			/* Pointer to a TmSortKeyType */

typedef struct  {			/* Merge decay data marker */
	LbUsFourByte		numPhotons;			/* Number of photons in the decay data */
	PHG_Decay 			decayData;			/* The decay data */
} DecayMergeData;
typedef struct  {			/* Merge input file data */
	FILE				*theFile;			/* The input file info */
	DecayMergeData		*bufStart;			/* Start of the merge data buffer */
	LbUsFourByte		bufSize;			/* Size of the merge data buffer */
	LbUsFourByte		bufDecays;			/* Number of decays in the buffer */
	PHG_Decay 			nextDecay;			/* The next decay to read from the file */
	Boolean				nextDecayPending;	/* Whether nextDecay is valid */
} MergeFileData;


/* LOCAL GLOBALS */
static LbUsOneByte			tmsortDecayFlag;				/* Decay identification flag */
static LbUsOneByte			tmsortPhotonFlag;				/* Photon identification flag */
static Boolean				tmsortCancelled;				/* Global cancellation flag */
static Boolean				tmsortCustomFile = false;		/* True if custom history file used */
static char					tmsortErrStr[1024];				/* Error string storage */		
static LbUsFourByte			tmsortNumToProc;				/* Number of histories to process */
static char					tmsortHistName[1024];			/* Name of history file */
static char					tmsortSortNameBase[1024];		/* Base name of sorted history file */
static char					tmsortSortName[1024];			/* Name of sorted history file */
static char					tmsortHistParamsName[1024];		/* Name of history parameters file */
static LbUsFourByte			tmsortArgIndex;					/* Index through command line arguments */
static PhoHFileHdrTy		tmsortHdrParams;				/* Input header */
static LbHdrHkTy			tmsortHistHeaderHk;				/* Hook to original history file header */
static PhoHFileHkTy			tmsortHistParamsHk;				/* Hook to history file information */
static PhoHFileHkTy			tmsortMergeFileHk;				/* Hook to merged file */
static Boolean				tmsortPreserveHeader;			/* Whether the header is altered or not */
static LbUsFourByte			tmsortDecaySize = 0;			/* Byte size of decay event */
static LbUsFourByte			tmsortMergeDecaySize = 0;		/* Byte size of merge decay event */
static LbUsFourByte			tmsortMergeFileSize = 0;		/* Byte size of merge file data */
static LbUsFourByte			tmsortPhotonSize = 0;			/* Byte size of photon event */
static TmSortKeyType		tmsortMinKey;					/* Minimum sort key; never actually occurs */
static TmSortKeyType		tmsortMaxKey;					/* Maximum sort key; never actually occurs */
static PHG_Decay			tmsortNextDecayEvent;			/* Storage for next partially read deacy */
static Boolean				tmsortDecayPending = false;		/* Status of tmsortNextDecayEvent */
static LbUsFourByte			tmsortDataBufferSizeParam = 0;	/* Size of main data buffer (Mbytes) */
static LbUsFourByte			tmsortDataBufferSize = 0;		/* Size of main data buffer (bytes) */
static PHG_Decay			*tmsortSortBufferPtr = NULL;	/* Main buffer holding decay data */
static PHG_Decay			*tmsortReserveBufferPtr = NULL;	/* Main buffer section for excess data */
static PHG_Decay			*tmsortResBufferStart = NULL;	/* Original start of reserve buffer */
static DecaySortBlock		*tmsortTempDecayBlocksList = NULL;
															/* Dynamic list of DecaySortBlock */
static LbUsFourByte			tmsortBlockSize;				/* Number of DecaySortData in a DecaySortBlock */
static PHG_Decay			**tmsortDecaySlots[PHG_MAX_DETECTED_PHOTONS+1];
															/* Lists of decay data slots in 
																tmsortSortBufferPtr, by 
																number of photons */
static LbUsFourByte			tmsortDecaySlotCounts[PHG_MAX_DETECTED_PHOTONS+1];
															/* Counts of the list sizes of 
																tmsortDecaySlots */
static LbUsFourByte			tmsortDecaySlotMaxs[PHG_MAX_DETECTED_PHOTONS+1];
															/* Maximums for the list sizes of 
																tmsortDecaySlotCounts */
static LbUsFourByte			tmsortDecaySlotDists[TMSORT_DIST_CT][PHG_MAX_DETECTED_PHOTONS+1];
															/* Counts of the list sizes for 
																the last few reads */
static LbUsFourByte			tmsortCurDist;					/* Current read slot distribution */
static DecayMergeData		*tmsortMergeBufferPtr = NULL;	/* Pointer to main buffer during merging */
static MergeFileData		*tmsortMergeFiles = NULL;		/* Array of files to be merged */
#ifdef TMSORT_PHO_DIST
static LbUsFourByte			tmsortPhotonCounts[PHG_MAX_DETECTED_PHOTONS+1];
															/* Count of total decays read, by photons */
#endif


/* PROTOTYPES */
Boolean			tmsort(int argc, char *argv[]);
Boolean 		tmsortInitialize(int argc, char *argv[]);
Boolean			tmsortGetParams(void);
void			tmsortTerminate(void);
Boolean			tmsortTest(char *argv[]);
Boolean			tmsortTimeSort(char *argv[]);
Boolean			tmsortSort(LbUsFourByte *numSortedFiles);
Boolean			tmsortMerge(LbUsFourByte sortFilesToMerge, LbUsFourByte *finalFileNum);
Boolean			tmsortMergeFileBatch(LbUsFourByte numSortFiles);
Boolean			tmsortGetFirstBatch(FILE *historyFile, 
					LbUsFourByte bufferSize, LbUsFourByte *numDecays);
Boolean			tmsortGetNextBatch(FILE *historyFile, LbSortListPtr sortListPtr, 
					LbUsFourByte *reserveSpace, LbUsFourByte *numDecays);
void			tmsortClearReserveBuffer(LbSortListPtr sortListPtr, 
					LbUsFourByte *reserveSpace, LbUsFourByte *clearedDecays);
LbUsFourByte	tmsortEstimateReadDecays(void);
void			tmsortFillMergeBuffer(LbUsFourByte fileIndex);
Boolean			tmsortSortDecay(LbSortListPtr sortListPtr, 
					PHG_Decay *decayPtr, LbUsFourByte numPhotons);
Boolean			tmsortStoreDecay(PHG_Decay *decayPtr, LbUsFourByte numPhotons,
					PHG_Decay **decayLocation);
void			tmsortRemoveDecay(PHG_Decay *decayLocation, LbUsFourByte numPhotons);
Boolean			tmsortGetDecaySlot(LbUsFourByte numPhotons, PHG_Decay **decaySlot);
LbFourByte		tmsortKeyComparator(LbSortKeyPtr keyPtr1, LbSortKeyPtr keyPtr2);
void			tmsortGetDecaySortKey(PHG_Decay *decayPtr, TmSortKeyPtr decaySortKeyPtr);
void			tmsortSetDecaySortKey(PHG_Decay *decayPtr, TmSortKeyType decaySortKey);
Boolean			tmsortReadDecay(FILE *historyFile, 
					PHG_Decay *decayPtr, 
					LbUsFourByte *dataSize,
					LbUsFourByte *numPhotons);
EventTy			tmsortReadEvent(FILE *historyFile,
					PHG_Decay *decayPtr,
					PHG_DetectedPhoton *photonPtr);
Boolean			tmsortWriteDecay(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay,
					PHG_DetectedPhoton *decayPhotons, LbUsFourByte numPhotons);
Boolean			tmsortWriteDetections(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay,
					PHG_DetectedPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_DetectedPhoton *pinkPhotons, LbUsFourByte numPinkPhotons);


/* FUNCTIONS */
/**********************
*	Name:		tmsort
*
*	Purpose:	Execute the program.
*
*	Result:		True unless an error occurs.
***********************/
Boolean tmsort(int argc, char *argv[])
{
	Boolean			okay = false;	/* Process Loop */
	time_t 			curTime;		/* Current time for stamping execution date */
	LbFourByte		argIndex;		/* Index through command line arguments */
	char			executionTimeStr[33];
									/* String for execution date conversion */
									/* Error condition results */
	/*
	Boolean			isInError;
	ErMgCdTy		managerCode;
	ErErCdTy		errorCode;
	*/
	
	
	/*** NOTE:  Custom history files have not been tested. ***/
	
	
	do { /* Process Loop */
	
		/* Perform initialization tasks */
		if (tmsortInitialize(argc, argv) == false) {
			goto FAIL;
		}
		
		/* Get current time */
		time(&curTime);
		
		/* Print out the command line */
		LbInPrintf("\n\n\nCommand line: ");
		for (argIndex = 0; argIndex < argc; argIndex++)
			LbInPrintf("%s ", argv[argIndex]);
		LbInPrintf("\n");
		
		/* Convert time to string format */
		strftime(executionTimeStr, 31, "%e-%b-%Y (%H:%M:%S)", localtime(&curTime));
		/*
		#ifdef __MWERKS__
			strftime(executionTimeStr, 31, "%m %d %Y (%I:%M:%S %p)", localtime(&curTime));
		#elif defined WINNT
			strftime(executionTimeStr, 31, "%x (%X)", localtime(&curTime));
		#else 
			strftime(executionTimeStr, 31, "%m %d %Y (%r)", localtime(&curTime));
		#endif
		*/
		
		LbInPrintf("\nExecution of timesort occurred on %s\n", executionTimeStr);
		
		if (TMSORT_IsUsePHGHistory()) {
			LbInPrintf("Using PHG history files.\n");
		}
		else if (TMSORT_IsUseColHistory()) {
			LbInPrintf("Using collimator history files.\n");
		}
		else if (TMSORT_IsUseDetHistory()) {
			LbInPrintf("Using detector history files.\n");
		}
		
		
		if (TMSORT_IsTestSorted()) {
			/* Just check the sorting */
			okay = tmsortTest(argv);
			if (! okay) {
				/* Error should already have been recorded */
				goto FAIL;
			}
		}
		else {
			/* Do the sorting */
			okay = tmsortTimeSort(argv);
			if (! okay) {
				/* Error should already have been recorded */
				goto FAIL;
			}
		}
		
		
		/* Get current time */
		time(&curTime);
		
		/* Convert time to string format */
		strftime(executionTimeStr, 31, "%e-%b-%Y (%H:%M:%S)", localtime(&curTime));
		/*
		#ifdef __MWERKS__
			strftime(executionTimeStr, 31, "%m %d %Y (%I:%M:%S %p)", localtime(&curTime));
		#elif defined WINNT
			strftime(executionTimeStr, 31, "%x (%X)", localtime(&curTime));
		#else 
			strftime(executionTimeStr, 31, "%m %d %Y (%r)", localtime(&curTime));
		#endif
		*/
		
		/* Print out our final message */
		LbInPrintf("\nExecution of timesort finished on %s\n", executionTimeStr);
		
		okay = true;
		FAIL:;
		
		/* Report any errors (locally) */
		/*
		if (ErIsInError()) {
			ErWhatError(&isInError, &managerCode, &errorCode, tmsortErrStr);
			sprintf(tmsortErrStr, "Failed timesort (%d,%d)", 
										managerCode, errorCode);
			ErHandle(tmsortErrStr, false);
			okay = true;
		}
		*/
	} while (false);
	
	/* Terminate PHG modules */
	tmsortTerminate();
	
	/* Handle error situation if one exists */
	if (!okay && tmsortCancelled) {
		ErHandle("User cancelled timesort", false);
		okay = true;
	}
	
	/* Quit the program */
	return (okay);
}


/**********************
*	Name:		tmsortInitialize
*
*	Purpose:	Perform initialization tasks.
*
*	Result:		True unless an error occurs.
***********************/
Boolean tmsortInitialize(int argc, char *argv[])
{
	Boolean				okay = false;				/* Process Loop */
	char				*knownOptions[] = {"pcdtr"};
	char				optArgs[TMSORT_NumFlags][LBEnMxArgLen];
	LbUsFourByte		optArgFlags = (LBFlag0);
	LbUsFourByte		typeCount;			/* Count of file types requested */
	char				fileName[1024];		/* Name of param file */
	LbFourByte			randSeed;			/* Seed for random generator */
	
	
	/* Set our type flags, used later on for marking what type of event is being written */
	PHG_SetIsADecay(tmsortDecayFlag);
	PHG_SetIsAPhoton(tmsortPhotonFlag);
	
	do { /* Process Loop */
		
		/* Get our options */
		if (!LbEnGetOptions(argc, argv, knownOptions,
				&PhgOptions, optArgs, optArgFlags, &tmsortArgIndex)) {
	
			break;
		}
		
		/* Make sure they specified exactly one history file */
		typeCount = 0;
		if (TMSORT_IsUsePHGHistory())
			typeCount++;
		if (TMSORT_IsUseColHistory())
			typeCount++;
		if (TMSORT_IsUseDetHistory())
			typeCount++;
		if (typeCount > 1) {
			ErStGeneric("You can only specify one type of history file, -p, -c, or -d.");
			break;
		}
		if (typeCount == 0) {
			ErStGeneric("You must specify the use of a PHG history file (-p)\n"
				" or a collimator history file (-c)\n"
				" or a detector history file (-d).");
			break;
		}
		
		/* See if they supplied a command line argument */
		if ((tmsortArgIndex != 0) && (argv[tmsortArgIndex] != 0)) {
		
			/* Get first param file and save number of param files to process */
			strcpy(PhgRunTimeParams.PhgParamFilePath,argv[tmsortArgIndex]);
			tmsortNumToProc = (argc - tmsortArgIndex);
		}
		else {
			/* Ask for the file name */
			LbInAsk("Enter name of param file", 0, false,
					&tmsortCancelled, 0, 0, 0, 0,
					fileName);
		
			/* Bolt if we canceled */
			if (tmsortCancelled) {
				ErStCancel("User cancelled.");
				goto CANCEL;
			}
			tmsortNumToProc = 1;
			strcpy(PhgRunTimeParams.PhgParamFilePath, fileName);
		}
		
		/* Get our run-time parameters */
		if (!tmsortGetParams())
			break;
		
		/* Clear the file name parameters */
		tmsortHistParamsName[0] = '\0';
		tmsortHistName[0] = '\0';
		
		/* Let them know what is going on */
		LbInPrintf("\n***********************\n");
		LbInPrintf("About to process parameter file '%s'\n", PhgRunTimeParams.PhgParamFilePath);
		LbInPrintf("\n***********************\n");
		
		
		/* If user requested to sort PHG history file, use the one specified in
			the param file.
		*/
		if (TMSORT_IsUsePHGHistory() &&
				(strlen(PhgRunTimeParams.PhgPhoHFileHistoryFilePath) == 0)) {
			
			ErStGeneric("No history file supplied in run time parameters.");
			break;
		}
		else {
			strcpy(tmsortHistName, PhgRunTimeParams.PhgPhoHFileHistoryFilePath);
			strcpy(tmsortHistParamsName, PhgRunTimeParams.PhgPhoHParamsFilePath);
		}
		if (strlen(tmsortHistParamsName) != 0) {
			/* Custom history file */
			tmsortCustomFile = true;
		}
		else {
			tmsortCustomFile = false;
		}
		
		/* Set decay and photon sizes */
		
		/* Both standard and custom files use PHG_Decay */
		tmsortDecaySize = sizeof(PHG_Decay);
		/* Both standard and custom files use PHG_Decay */
		tmsortMergeDecaySize = sizeof(DecayMergeData);
		/* Both standard and custom files use PHG_Decay */
		tmsortMergeFileSize = sizeof(MergeFileData);
		if (! tmsortCustomFile) {
			/* Standard files use PHG_DetectedPhoton */
			tmsortPhotonSize = sizeof(PHG_DetectedPhoton);
		}
		else {
			/* Custom files use PHG_TrackingPhoton */
			tmsortPhotonSize = sizeof(PHG_TrackingPhoton);
		}
		
		/* Initialize the default minimum and maximum key values */
		tmsortMinKey = TMSORT_MINKEY;
		tmsortMaxKey = TMSORT_MAXKEY;
		
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
*	Name:			tmsortGetParams
*
*	Summary:		Read in the timesort parameters.
*
*	Arguments:		None.
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean tmsortGetParams(void)	
{
	LbUsFourByte			loopV =0;
	LbPfHkTy				phgParamFileHk;					/* The parameter file */
	double					paramBuffer[LBPF_PARAM_LEN];
	Boolean					okay = false;					/* Process flag */
	Boolean					isEOF;							/* End of file flag */
	Boolean					switchOkay;						/* Switch flag */
	LbPfEnPfTy				paramType;
	LbUsFourByte			paramSize;
	char					paramLabel[LBPF_LABEL_LEN];
	PhgEn_RunTimeParamsTy	whichParam;
	
	
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
		
		
		/* Attempt to open parameter file */
		if (!LbPfOpen(PhgRunTimeParams.PhgParamFilePath, 0, &phgParamFileHk)) {
			sprintf(tmsortErrStr, "An error occurred opening the command line input parameters file '%s'.",
				PhgRunTimeParams.PhgParamFilePath);
			ErAlert(tmsortErrStr, false);
			break;
		}
		
		/* Read the parameters */
		while(LbPfGetParam(&phgParamFileHk, (void *)paramBuffer,
				&paramType, &paramSize, paramLabel, &isEOF)) {
			
			/* Find the runtime parameter */
			whichParam = PhgLookupRunTimeParamLabel(paramLabel);
			switchOkay = false;
			if (whichParam == PhgEn_NULL) {
				/* Check local list */
				do {
					/* Sorted output file name */
					if (strcmp(paramLabel, "sorted_history_file") == 0) {
						strcpy(tmsortSortNameBase, (char *)paramBuffer);
						switchOkay = true;
						break;
					}
					
					/* Buffer memory size */
					if (strcmp(paramLabel, "buffer_size") == 0) {
						tmsortDataBufferSizeParam = *((LbUsFourByte *)paramBuffer);
						switchOkay = true;
						break;
					}
				} while (false);
			}
			
			if (!switchOkay) {
				switchOkay = true;
				switch (whichParam) {
					
					case PhoHFileEn_history_file:
							strcpy(PhgRunTimeParams.PhgPhoHFileHistoryFilePath,
								(char *) paramBuffer);
							
							PhgRunTimeParams.PhgIsHistoryFile = 
								(PhgRunTimeParams.PhgPhoHFileHistoryFilePath[0] != '\0');
						break;
					
					case PhoHFileEn_history_params_file:
							strcpy(PhgRunTimeParams.PhgPhoHParamsFilePath,
								(char *) paramBuffer);
							
							PhgRunTimeParams.PhgIsHistoryParamsFile = 
								(PhgRunTimeParams.PhgPhoHParamsFilePath[0] != '\0');
						break;
					
					case PhgEn_NULL:
						sprintf(tmsortErrStr, "(tmsortGetParams) Unknown (hence unused) parameter (%s).",
	                        paramLabel);
						ErAlert(tmsortErrStr, false);
						break;
					
					default:
						sprintf(tmsortErrStr, "(tmsortGetParams) Unknown (hence unused) parameter (%s).",
	                        paramLabel);
						ErAlert(tmsortErrStr, false);
						break;
				}
			}
			if (!switchOkay)
				break;
		}
		
		/* Close the parameter file */
		LbPfClose(&phgParamFileHk);
		
		/* See if we quit due to error */
		if (!isEOF)
			break;
		
		okay = true;
	} while (false);
	
	return (okay);
}


/**********************
*	Name:		tmsortTerminate
*
*	Purpose:	Terminate all of the managers.
*
*	Result:		None.
***********************/
void tmsortTerminate()
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


/**********************
*	Name:		tmsortTest
*
*	Purpose:	Test a standard history file to see if it is sorted.
*
*	Result:		True unless an error occurs.
***********************/
Boolean tmsortTest(char *argv[])
{
	Boolean				okay = false;				/* Process flag (function return) */
	LbUsFourByte		curFileIndex;				/* Current file index */
	FILE				*historyFile;				/* The history file we are going to process */
	LbHdrHkTy			histHeaderHk;				/* Hook to history file header */
	EventTy				eventType;					/* Type of current event */
	TmSortKeyType		minimumKey;					/* Smallest sort key for outputs */
	Boolean				dataSorted;					/* Indicator of result */
	LbUsFourByte		numDecaysChecked;			/* Number of decays checked */
	LbUsFourByte		dataSize;					/* Size of single decay data */
	LbUsFourByte		numPhotons;					/* Number of photons in the decay */
	TmSortKeyType		decaySortKey;				/* Current decay time sort key */
	
	
	do { /* Process Loop */
		
		/* Process the files */
		for (curFileIndex = 1; curFileIndex <= tmsortNumToProc; curFileIndex++){
			
			/* Report the input file name(s) */
			LbInPrintf("\n");
			if (tmsortCustomFile) {
				LbInPrintf("Using custom parameters file:\n"
							"  '%s'.\n", 
							tmsortHistParamsName);
			}
			LbInPrintf("Checking whether file:\n"
						"  '%s'\n"
						"  is sorted...\n", 
						tmsortHistName);
			
			/* Open history file file */
			if ((historyFile = LbFlFileOpen(tmsortHistName, "rb")) == 0) {
				sprintf(tmsortErrStr, "Unable to open history file\n'%s'.",
					tmsortHistName);
				ErStFileError(tmsortErrStr);
				goto FAIL;
			}
			
			/* Read in the original header */
			if (PhgHdrGtParams(historyFile, &tmsortHdrParams, &histHeaderHk) == false){
				sprintf(tmsortErrStr, "Unable to read history file header\n'%s'.",
					tmsortHistName);
				ErStFileError(tmsortErrStr);
				goto FAIL;
			}
			
			/* Verify that requested file is of the right type */
			{
				if (TMSORT_IsUsePHGHistory() && (tmsortHdrParams.H.HdrKind != PhoHFileEn_PHG)) {
					ErStGeneric("File specified as PHG history file is not valid.");
					goto FAIL;
				}
				
				if (TMSORT_IsUseColHistory() && (tmsortHdrParams.H.HdrKind != PhoHFileEn_COL)) {
					ErStGeneric("File specified as Collimator history file is not valid.");
					goto FAIL;
				}
				
				if (TMSORT_IsUseDetHistory() && (tmsortHdrParams.H.HdrKind != PhoHFileEn_DET)) {
					ErStGeneric("File specified as Detector history file is not valid.");
					goto FAIL;
				}
			}
			
			if (tmsortCustomFile) {
				/* Setup header hook */
				
				tmsortHistParamsHk.doCustom = true;
				
				/* Do custom parameter initialization  */
				if (PhoHFileGetRunTimeParams(tmsortHistParamsName, &(tmsortHistParamsHk.customParams)) == false) {
					sprintf(tmsortErrStr,"Unable to get custom parameters for history file named '%s'",
						tmsortHistParamsName);
					ErAlert(tmsortErrStr, false);
					goto FAIL;
				}
				
				strcpy(tmsortHistParamsHk.customParamsName, tmsortHistParamsName);
				
				tmsortHistParamsHk.bluesReceived = 0;
				tmsortHistParamsHk.bluesAccepted = 0;
				tmsortHistParamsHk.pinksReceived = 0;
				tmsortHistParamsHk.pinksAccepted = 0;
				tmsortHistParamsHk.histFile = historyFile;
			}
			
			/* Allocate a small decay data buffer */
			tmsortDataBufferSize = 2*tmsortDecaySize + PHG_MAX_DETECTED_PHOTONS*tmsortPhotonSize;
			if ((tmsortSortBufferPtr = 
					(PHG_Decay *) LbMmAlloc(tmsortDataBufferSize)) == NULL) {
				sprintf(tmsortErrStr, "Unable to allocate memory:  tmsortSortBufferPtr.");
				ErStGeneric(tmsortErrStr);
				goto FAIL;
			}
			
			if (! tmsortCustomFile) {
				/* Fill the first decay event buffer (decay reading assumes decay 
					event always already read in from file */
				{
					tmsortDecayPending = false;
					eventType = tmsortReadEvent(historyFile, &tmsortNextDecayEvent, 
									(PHG_DetectedPhoton *)tmsortSortBufferPtr);
					if (eventType != Decay) {
						ErStGeneric("Expected first event to be decay, and it wasn't.");
						goto FAIL;
					}
					else {
						tmsortDecayPending = true;
					}
				}
			}
			
			/* Initialize the smallest sort key */
			minimumKey = tmsortMinKey;
			
			
			/* Do the checking */
			dataSorted = true;		/* Hope for the best */
			numDecaysChecked = 0;
			while (tmsortReadDecay(historyFile, tmsortSortBufferPtr, &dataSize, &numPhotons)) {
				numDecaysChecked++;
				
				/* Get the decay sort key */
				tmsortGetDecaySortKey(tmsortSortBufferPtr, &decaySortKey);
				
				if (decaySortKey >= minimumKey) {
					/* Still in sorted order; update the key */
					minimumKey = decaySortKey;
				}
				else {
					/* Not sorted */
					dataSorted = false;
					break;
				}
			}
			
			/* Close the history file file */
			fclose(historyFile);
			
			/* Free the buffer */
			if (tmsortSortBufferPtr != NULL) {
				LbMmFree((void **) &tmsortSortBufferPtr);
			}
			
			/* Report the results */
			if (dataSorted) {
				LbInPrintf("This file is already in sorted order.\n");
			}
			else {
				LbInPrintf("This file is not in its desired sorted order.\n");
			}
			
			
			/* Open next file */
			if (curFileIndex < tmsortNumToProc) {
				
				strcpy(PhgRunTimeParams.PhgParamFilePath,argv[curFileIndex+tmsortArgIndex]);
				
				/* Get our run-time parameters */
				if (!tmsortGetParams()) {
					goto FAIL;
				}
				
				/* Let them know what is going on */
				LbInPrintf("\n***********************\n");
				LbInPrintf("About to process parameter file '%s'\n", PhgRunTimeParams.PhgParamFilePath);
				LbInPrintf("\n***********************\n");

				/* Set the next history file to the next collimator history */
				if (TMSORT_IsUseColHistory() &&
						(strlen(ColRunTimeParams[ColCurParams].ColHistoryFilePath) == 0)){
		
					/* Let the user know they must supply the collimator history file */
					ErStGeneric("You specified using a collimatory history file (-c)\n"
						"but didn't supply one in the collimator parameter file");
					break;
				}
				else {
					/* Save the current history file name */
					strcpy(tmsortHistName, ColRunTimeParams[ColCurParams].ColHistoryFilePath);
				}
				
				/* Or, set the next history file to the PHG history file */
				if (TMSORT_IsUsePHGHistory() &&
						(strlen(PhgRunTimeParams.PhgPhoHFileHistoryFilePath) == 0)) {
						
						ErStGeneric("No history file supplied in run time parameters.");
						break;
				}
				else {
					strcpy(tmsortHistName, PhgRunTimeParams.PhgPhoHFileHistoryFilePath);
				}
				if (strlen(tmsortHistParamsName) != 0) {
					/* Custom history file */
					tmsortCustomFile = true;
				}
				else {
					tmsortCustomFile = false;
				}
				
				/* Set decay and photon sizes */
				
				/* Both standard and custom files use PHG_Decay */
				tmsortDecaySize = sizeof(PHG_Decay);
				/* Both standard and custom files use PHG_Decay */
				tmsortMergeDecaySize = sizeof(DecayMergeData);
				/* Both standard and custom files use PHG_Decay */
				tmsortMergeFileSize = sizeof(MergeFileData);
				if (! tmsortCustomFile) {
					/* Standard files use PHG_DetectedPhoton */
					tmsortPhotonSize = sizeof(PHG_DetectedPhoton);
				}
				else {
					/* Custom files use PHG_TrackingPhoton */
					tmsortPhotonSize = sizeof(PHG_TrackingPhoton);
				}
			}
		}		

		okay = true;
		FAIL:;
	} while (false);
	
	
	/* Free the buffer */
	if (tmsortSortBufferPtr != NULL) {
		LbMmFree((void **) &tmsortSortBufferPtr);
	}
	
	return (okay);
}


/*********************************************************************************
*
*	Name:			tmsortTimeSort
*
*	Summary:		Sort a history file.
*
*	Arguments:
*		char			*argv[]					- Command line arguments.
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean tmsortTimeSort(char *argv[])
{
	Boolean				okay = false;				/* Process flag (function return) */
	LbUsFourByte		curFileIndex;				/* Current file index */
	LbUsFourByte		sortFilesToMerge;			/* Number of smaller sort files produced */
	LbUsFourByte		finalFileNumber;			/* Final index of sorted output file */
	
	
	do { /* Process Loop */
		
		/* Process the files */
		for (curFileIndex = 1; curFileIndex <= tmsortNumToProc; curFileIndex++){
			
			/* Validate the name of the sorted files (supplied by params file) */
			if (strlen(tmsortSortNameBase) == 0) {
				sprintf(tmsortErrStr, "No sorted output file name supplied.");
				ErStGeneric(tmsortErrStr);
				goto FAIL;
			}
			
			/* Report the parameters */
			if (tmsortCustomFile) {
				LbInPrintf("Using custom parameters file:\n"
							"  '%s'.\n", 
							tmsortHistParamsName);
			}
			LbInPrintf("Sorting input file:\n"
						"  '%s'\n"
						"  into output file:\n"
						"  '%s'.\n", 
						tmsortHistName, 
						tmsortSortNameBase);
			if (tmsortDataBufferSizeParam == 1) {
				LbInPrintf("Buffer size is %ld megabyte.\n", 
							(unsigned long)tmsortDataBufferSizeParam);
			}
			else {
				LbInPrintf("Buffer size is %ld megabytes.\n", 
							(unsigned long)tmsortDataBufferSizeParam);
			}
			
			/* Check size of main data buffer (supplied by params file) */
			/* Must be large enough for 3% reserve buffer to hold some full decays */
			tmsortDataBufferSize = tmsortDataBufferSizeParam * 1048576;
			if (tmsortDataBufferSize < 200000) {	/* 3000000 needed for custom files */
				sprintf(tmsortErrStr, "Supplied buffer memory size (%ld) is too small.", 
										(unsigned long)tmsortDataBufferSizeParam);
				ErStGeneric(tmsortErrStr);
				goto FAIL;
			}
			
			/* Allocate main decay data buffer */
			if ((tmsortSortBufferPtr = 
					(PHG_Decay *) LbMmAlloc(tmsortDataBufferSize)) == NULL) {
				sprintf(tmsortErrStr, "Unable to allocate memory:  tmsortSortBufferPtr.");
				ErStGeneric(tmsortErrStr);
				goto FAIL;
			}
			
			
			LbInPrintf("\nBeginning creation of intermediate subfiles...\n");
			
			/* Sort the history file into multiple smaller sort files (phase I) */
			if (! tmsortSort(&sortFilesToMerge)) {
				/* Error should already have been reported */
				goto FAIL;
			}
			
			if (sortFilesToMerge == 1) {
				LbInPrintf("Created a single sorted subfile; no merging required.\n");
			}
			else {
				LbInPrintf("Created %ld sorted subfiles, now to be merged...\n", 
							(unsigned long)sortFilesToMerge);
			}
			
			
			/* Delete the original file, if requested */
			if (TMSORT_IsDeleteFile()) {
				LbInPrintf("Deleting input history file, as requested.\n");
				if (remove(tmsortHistName) != 0) {
					sprintf(tmsortErrStr, "Unable to delete history file\n'%s'.",
						tmsortHistName);
					ErStFileError(tmsortErrStr);
					/* But keep going */
				}
			}
			
			
			/* Merge the separate sorted files (phase II) */
			if (sortFilesToMerge == 1) {
				/* Only one sort file produced; nothing to merge */
				finalFileNumber = 1;
			}
			else {
				/* Merge the multiple sorted files */
				if (! tmsortMerge(sortFilesToMerge, &finalFileNumber)) {
					/* Error should already have been reported */
					goto FAIL;
				}
			}
			
			
#ifndef TMSORT_KEEP_TEMP_FILES
			/* Rename the final file to the desired output file name */
			sprintf(tmsortSortName, "%s%ld", 
					tmsortSortNameBase, (unsigned long)finalFileNumber);
			remove(tmsortSortNameBase);
			if (rename(tmsortSortName, tmsortSortNameBase) != 0) {
				sprintf(tmsortErrStr, "Unable to rename last sorted file\n'%s'\n"
					"to final sorted file\n'%s'.", 
					tmsortSortName, tmsortSortNameBase);
				ErStFileError(tmsortErrStr);
				/* But keep going */
			}
#endif
			
			LbInPrintf("Final file and sorting process now completed.\n");
			
			/* Free the main buffer */
			if (tmsortSortBufferPtr != NULL) {
				LbMmFree((void **) &tmsortSortBufferPtr);
			}
			
			
			/* Open next parameters file */
			if (curFileIndex < tmsortNumToProc) {
				
				strcpy(PhgRunTimeParams.PhgParamFilePath,argv[curFileIndex+tmsortArgIndex]);
				
				/* Get our run-time parameters */
				if (!tmsortGetParams()) {
					goto FAIL;
				}
				
				/* Let them know what is going on */
				LbInPrintf("\n***********************\n");
				LbInPrintf("About to process parameter file '%s'\n", PhgRunTimeParams.PhgParamFilePath);
				LbInPrintf("\n***********************\n");

				/* Set the next history file to the next collimator history */
				if (TMSORT_IsUseColHistory() &&
						(strlen(ColRunTimeParams[ColCurParams].ColHistoryFilePath) == 0)){
		
					/* Let the user know they must supply the collimator history file */
					ErStGeneric("You specified using a collimatory history file (-c)\n"
						"but didn't supply one in the collimator parameter file");
					break;
				}
				else {
					/* Save the current history file name */
					strcpy(tmsortHistName, ColRunTimeParams[ColCurParams].ColHistoryFilePath);
				}
				
				/* Or, set the next history file to the PHG history file */
				if (TMSORT_IsUsePHGHistory() &&
						(strlen(PhgRunTimeParams.PhgPhoHFileHistoryFilePath) == 0)) {
						
						ErStGeneric("No history file supplied in run time parameters.");
						break;
				}
				else {
					strcpy(tmsortHistName, PhgRunTimeParams.PhgPhoHFileHistoryFilePath);
				}
				if (strlen(tmsortHistParamsName) != 0) {
					/* Custom history file */
					tmsortCustomFile = true;
				}
				else {
					tmsortCustomFile = false;
				}
				
				/* Set decay and photon sizes */
				
				/* Both standard and custom files use PHG_Decay */
				tmsortDecaySize = sizeof(PHG_Decay);
				/* Both standard and custom files use PHG_Decay */
				tmsortMergeDecaySize = sizeof(DecayMergeData);
				/* Both standard and custom files use PHG_Decay */
				tmsortMergeFileSize = sizeof(MergeFileData);
				if (! tmsortCustomFile) {
					/* Standard files use PHG_DetectedPhoton */
					tmsortPhotonSize = sizeof(PHG_DetectedPhoton);
				}
				else {
					/* Custom files use PHG_TrackingPhoton */
					tmsortPhotonSize = sizeof(PHG_TrackingPhoton);
				}
			}
		}		
		
		okay = true;
		FAIL:;
	} while (false);
	
	
	/* Free memory */
	if (tmsortSortBufferPtr != NULL) {
		LbMmFree((void **) &tmsortSortBufferPtr);
	}
	
	return (okay);
}


/*********************************************************************************
*
*	Name:			tmsortSort
*
*	Summary:		Sort a large history file into multiple, smaller sorted files.
*
*	Arguments:
*		LbUsFourByte	 	*numSortedFiles			- The # of files produced.
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean tmsortSort(LbUsFourByte *numSortedFiles)
{
	#define		kTSProgressMinutes		10	/* Time interval between progress displays */
	#define		kTSProgressFiles		10	/* File interval between progress displays */
	
	Boolean				okay = false;				/* Function return value */
	LbUsFourByte		i;							/* Generic index variable */
	LbUsFourByte		j;							/* Index through decay distributions */
	FILE				*historyFile;				/* The history file to sort */
	/*double				histFileSize;*/				/* Size of the history file */
	LbUsFourByte		accumSecs;					/* Seconds since last progress display */
	double				accumCPUSecs;				/* Total used CPU time */
	LbUsFourByte		nextAccumSecs;				/* Time between progress displays */
	LbUsEightByte		sortedFilesSize;			/* Cumulative size of the sorted output files */
	LbUsFourByte		nextFileDisplay;			/* Next file display point */
	/*double				nextSizeDisplay;*/			/* Next cumulative size for a report message */
	LbTmTimingType		clockTiming;				/* Timing data for progress messages */
	LbUsFourByte		maxReserveSpace;			/* Original size of reserve buffer */
	LbUsFourByte		reserveSpace;				/* Current size of reserve buffer */
	EventTy				eventType;					/* Type of current event */
	LbUsFourByte		numDecaysRead;				/* Number of decays read in so far */
	LbUsFourByte		numDecaysWritten;			/* Number of decays written out so far */
	LbUsFourByte		numInitialDecays;			/* Number of initial buffer decays */
	LbSortItemPtr		listItemPtr;				/* Pointer to sorted list item */
	LbSortListPtr		sortListPtr;				/* Pointer to sort sorted list */
	DecaySortBlock		*curBlockPtr;				/* Pointer to the current block */
	DecaySortData		*sortDataPtr;				/* Pointer to the current data */
	LbUsFourByte		ioDecayAmount;				/* Usual number of buffer decays in i/o ops */
	LbUsFourByte		sortFileNumber;				/* Index to current sorted output file */
	Boolean				allDone;					/* File has been sorted into separate files */
	PhoHFileHkTy		sortFileHk;					/* Hook to sorted file */
	TmSortKeyType		minimumKey;					/* Smallest sort key for outputs */
	LbUsFourByte		ioNumDecays;				/* Actual number of buffer decays in an i/o op */
	LbSortItemPtr		minSortItem;				/* Item containing last minimumKey */
	Boolean				foundDecay;					/* Good decay found in sorted list */
	LbUsFourByte		numPhotons;					/* Number of photons in the decay */
	PHG_Decay			*nextDecayPtr;				/* Current output decay data */
	LbUsFourByte		clearedDecays;				/* Count of decays cleared from reserve buffer */
	LbUsFourByte		estNumDecays;				/* Estimated count of decays to read in */
	Boolean				gotData;					/* Decays were read in from file */
	double				realSecs;					/* Wall clock time elapsed */
	double				cpuSecs;					/* CPU time elapsed */
	time_t				curTime;					/* Current wall clock time */
	char				timeStr[32];				/* String for printing the time */
	
	
	/* Initialize variables */
	tmsortTempDecayBlocksList = NULL;
	for (i=0; i<=PHG_MAX_DETECTED_PHOTONS; i++) {
		tmsortDecaySlots[i] = NULL;
		
		#ifdef TMSORT_PHO_DIST
			tmsortPhotonCounts[i] = 0;
		#endif
	}
	for (j=0; j<TMSORT_DIST_CT; j++) {
		for (i=0; i<=PHG_MAX_DETECTED_PHOTONS; i++) {
			tmsortDecaySlotDists[j][i] = 0;
		}
	}
	sortListPtr = NULL;
	
	/* Open history file file */
	if ((historyFile = LbFlFileOpen(tmsortHistName, "rb")) == 0) {
		sprintf(tmsortErrStr, "Unable to open history file\n'%s'.",
			tmsortHistName);
		ErStFileError(tmsortErrStr);
		goto FAIL;
	}
	
	/* Determine the size of the history file */
	/*
	if (fseek(historyFile, 0, SEEK_END) == 0) {
		histFileSize = (double)(ftello(historyFile));
	}
	else {
		histFileSize = -1.0;
	}
	fseek(historyFile, 0, SEEK_SET);
	*/
	
	/* Set up for progress display messages */
	accumSecs = 0;
	accumCPUSecs = 0.0;
	nextAccumSecs = kTSProgressMinutes*60;	/* = every xx minutes */
	sortedFilesSize = 0;
	nextFileDisplay = kTSProgressFiles;	/* = every xx sorted output files */
	LbTmStartTiming(&clockTiming);
	
	/* Read in the original header and save it */
	if (PhgHdrGtParams(historyFile, &tmsortHdrParams, &tmsortHistHeaderHk) == false){
		sprintf(tmsortErrStr, "Unable to read history file header\n'%s'.",
			tmsortHistName);
		ErStFileError(tmsortErrStr);
		goto FAIL;
	}
	
	/* Verify that requested file is of the right type */
	{
		if (TMSORT_IsUsePHGHistory() && (tmsortHdrParams.H.HdrKind != PhoHFileEn_PHG)) {
			ErStGeneric("File specified as PHG history file is not valid.");
			goto FAIL;
		}
		
		if (TMSORT_IsUseColHistory() && (tmsortHdrParams.H.HdrKind != PhoHFileEn_COL)) {
			ErStGeneric("File specified as Collimator history file is not valid.");
			goto FAIL;
		}
		
		if (TMSORT_IsUseDetHistory() && (tmsortHdrParams.H.HdrKind != PhoHFileEn_DET)) {
			ErStGeneric("File specified as Detector history file is not valid.");
			goto FAIL;
		}
	}
	
	if (tmsortCustomFile) {
		/* Setup header hook */
		
		tmsortHistParamsHk.doCustom = true;
		
		/* Do custom parameter initialization  */
		if (PhoHFileGetRunTimeParams(tmsortHistParamsName, &(tmsortHistParamsHk.customParams)) == false) {
			sprintf(tmsortErrStr,"Unable to get custom parameters for history file named '%s'",
				tmsortHistParamsName);
			ErAlert(tmsortErrStr, false);
			goto FAIL;
		}
		
		strcpy(tmsortHistParamsHk.customParamsName, tmsortHistParamsName);
		
		tmsortHistParamsHk.bluesReceived = 0;
		tmsortHistParamsHk.bluesAccepted = 0;
		tmsortHistParamsHk.pinksReceived = 0;
		tmsortHistParamsHk.pinksAccepted = 0;
		tmsortHistParamsHk.histFile = historyFile;
	}
	
	/* If the history file has been sorted already then we preserve its header */
	/* If not, then the header is modified to show that it has been sorted. */
	/* The file is sorted anyway, no matter what the header indicates. */
	tmsortPreserveHeader = tmsortHdrParams.H.isTimeSorted;
	
	
	/* The tmsortMaxKey may be improved upon by making it closer to the upper 
		limit of possible decay time sort values */
	if (tmsortHdrParams.H.PhgRunTimeParams.Phg_LengthOfScan > 0) {
		/* The scan time exists, and no decay can occur after the end of scan */
		tmsortMaxKey = tmsortHdrParams.H.PhgRunTimeParams.Phg_LengthOfScan + 1.0;
	}
	
	/* Set aside a reserve buffer at the end of the main buffer for reading 
		new decay data and for decay misfits */
	maxReserveSpace = 0.03 * tmsortDataBufferSize;	/* 3% sounds good */
	reserveSpace = maxReserveSpace;
	/*reserveSpace = 32 * */		/* Arbitrary choice of number */
	/*	(tmsortDecaySize + PHG_MAX_DETECTED_PHOTONS*tmsortPhotonSize);*/
	tmsortReserveBufferPtr = (PHG_Decay *)((LbUsOneByte *)tmsortSortBufferPtr + 
									tmsortDataBufferSize-maxReserveSpace);
	tmsortResBufferStart = tmsortReserveBufferPtr;
	
	/* Allocate a start to the temporary, first batch decay list */
	tmsortBlockSize = 1000;		/* Arbitrary choice */
	if ((tmsortTempDecayBlocksList = (DecaySortBlock *) 
			LbMmAlloc(sizeof(DecaySortBlock))) == NULL) {
		sprintf(tmsortErrStr, "Unable to allocate memory:  tmsortTempDecayBlocksList.");
		ErStGeneric(tmsortErrStr);
		goto FAIL;
	}
	tmsortTempDecayBlocksList->nextSortBlock = NULL;
	if ((tmsortTempDecayBlocksList->thisDataBlock = (DecaySortData *) 
			LbMmAlloc(tmsortBlockSize*sizeof(DecaySortData))) == NULL) {
		sprintf(tmsortErrStr, "Unable to allocate memory:  DecaySortData.");
		ErStGeneric(tmsortErrStr);
		goto FAIL;
	}
	
	if (! tmsortCustomFile) {
		/* Fill the first decay event buffer (decay reading assumes decay 
			event always already read in from file */
		tmsortDecayPending = false;
		eventType = tmsortReadEvent(historyFile, &tmsortNextDecayEvent, 
						(PHG_DetectedPhoton *)tmsortSortBufferPtr);
		if (eventType != Decay) {
			ErStGeneric("Expected first event to be decay, and it wasn't.");
			goto FAIL;
		}
		else {
			tmsortDecayPending = true;
		}
	}
	
	/* Initialize decay counts */
	numDecaysRead = 0;
	numDecaysWritten = 0;
	
	/* Fill the sort buffer from the input file */
	if (! tmsortGetFirstBatch(historyFile, 
								tmsortDataBufferSize-reserveSpace, 
								&numInitialDecays)) {
		/* Unable to allocate enough memory; error message already sent */
		goto FAIL;
	}
	numDecaysRead = numInitialDecays;
	
	/* Initialize the read distribution with the one just read */
	tmsortCurDist = 0;
	for (i=0; i<=PHG_MAX_DETECTED_PHOTONS; i++) {
		tmsortDecaySlotDists[tmsortCurDist][i] = tmsortDecaySlotMaxs[i];
	}
	
	/* Make a new sorted list */
	if ((sortListPtr = LbSortNewList(numInitialDecays+1, 
						(LbSortKeyPtr)&tmsortMinKey, sizeof(TmSortKeyType), 
						tmsortKeyComparator)) == NULL) {
		/* Severe problems */
		sprintf(tmsortErrStr, "Got failure creating sorted list.");
		ErStGeneric(tmsortErrStr);
		goto FAIL;
	}
	
	/* Add an initial tmsortMinKey item (required by output algorithm) */
	if (! LbSortInsert(sortListPtr, (LbSortKeyPtr)&tmsortMinKey, 0, NULL) ) {
		/* Severe problems */
		sprintf(tmsortErrStr, "Got failure adding first sorted list item.");
		ErStGeneric(tmsortErrStr);
		goto FAIL;
	}
	
	/* Add the first batch of decays to the sorted list */
	/* NOTE:  Can't be done at time of reading because sorted list isn't 
		created until its size is known (after reading) */
	curBlockPtr = tmsortTempDecayBlocksList;
	while (curBlockPtr != NULL) {
		for (i=0; i<tmsortBlockSize; i++) {
			sortDataPtr = curBlockPtr->thisDataBlock + i;
			if (sortDataPtr->decayDataPtr != NULL) {
				/* Add the data to the sorted list */
				if (! tmsortSortDecay(sortListPtr, 
						sortDataPtr->decayDataPtr, sortDataPtr->numPhotons)) {
					sprintf(tmsortErrStr, "Unable to sort data; all results are invalid.");
					ErStGeneric(tmsortErrStr);
					goto FAIL;
				}
			}
			else {
				/* End of data */
				curBlockPtr = NULL;
				break;
			}
		}
		if (curBlockPtr != NULL) {
			/* Go to next block */
			curBlockPtr = curBlockPtr->nextSortBlock;
		}
	}
	
	/* Delete the temporary, first batch decay list */
	curBlockPtr = tmsortTempDecayBlocksList;
	while (curBlockPtr != NULL) {
		if (curBlockPtr->thisDataBlock != NULL) {
			/* Delete the data block */
			LbMmFree((void **) &curBlockPtr->thisDataBlock);
		}
		tmsortTempDecayBlocksList = curBlockPtr->nextSortBlock;
		LbMmFree((void **) &curBlockPtr);
		curBlockPtr = tmsortTempDecayBlocksList;
	}
	
	/* Allocate the decay slot lists */
	for (i=0; i<=PHG_MAX_DETECTED_PHOTONS; i++) {
		if (tmsortDecaySlotMaxs[i] > 0) {
			/* Allocate the memory */
			if ((tmsortDecaySlots[i] = (PHG_Decay **) 
					LbMmAlloc(tmsortDecaySlotMaxs[i]*sizeof(PHG_Decay *))) == NULL) {
				sprintf(tmsortErrStr, "Unable to allocate memory:  tmsortDecaySlots[%ld].", 
					(unsigned long)i);
				ErStGeneric(tmsortErrStr);
				goto FAIL;
			}
		}
	}
	
	/* Set the proportion of the sort buffer to be read/written */
	/* This size is chosen to have a large portion of decay records 
		in the sorted list, yet still allow bulk i/o to minimize disk 
		seeks.  Note:  the size has not been optimally determined.  */
	ioDecayAmount = numInitialDecays / 16;
	if (ioDecayAmount <= 0) {
		ioDecayAmount = numInitialDecays;
	}
	
	
	/* Sort the input file into as many output files as needed */
	sortFileNumber = 0;
	allDone = false;
	do {
		/* Create the next sortFileNumber numbered output file */
		sortFileNumber++;
		sprintf(tmsortSortName, 
				"%s%ld", tmsortSortNameBase, (unsigned long)sortFileNumber);
		if (PhoHFileCreate(tmsortSortName, "", 
				tmsortHdrParams.H.HdrKind, &sortFileHk) == false) {
			sprintf(tmsortErrStr,"Unable to create sorted file named:\n"
				"'%s'\n"
				" (tmsortSort)",
				tmsortSortName);
			ErStFileError(tmsortErrStr);
			goto FAIL;
		}
		
		/* Initialize the smallest sort key */
		minimumKey = tmsortMinKey;
		
		do {
			/* Write out the next batch of sorted output records */
			ioNumDecays = ioDecayAmount;
			if (LbSortFind(sortListPtr, (LbSortKeyPtr)&minimumKey, &minSortItem)) {
				for (i=0; i<ioNumDecays; i++) {
					foundDecay = LbSortNext(sortListPtr, 
						minSortItem, &listItemPtr, &numPhotons, (void **)&nextDecayPtr);
					if (foundDecay) {
						/* Update the minimum key */
						tmsortGetDecaySortKey(nextDecayPtr, &minimumKey);
						
						/* Write out the decay */
						if (! tmsortWriteDecay(&sortFileHk, nextDecayPtr,
								(PHG_DetectedPhoton *)(nextDecayPtr+1), numPhotons)) {
							sprintf(tmsortErrStr, "Got failure writing sort decay to disk.");
							ErStFileError(tmsortErrStr);
							goto FAIL;
						}
						
						/* Remove the decay from the sorted list */
						if (! LbSortDelete(sortListPtr, listItemPtr)) {
							/* Severe problems */
							sprintf(tmsortErrStr, "Got failure removing sort decay from sorted list.");
							ErStGeneric(tmsortErrStr);
							goto FAIL;
						}
						
						/* Clear the decay data's buffer spot */
						tmsortRemoveDecay(nextDecayPtr, numPhotons);
					}
					else {
						/* No more data >= minimumKey */
						
						ioNumDecays = i;	/* Decays found */
						
						break;
					}
				}
			}
			numDecaysWritten += ioNumDecays;
			
			/* Clear out the reserve buffer as much as possible */
			tmsortClearReserveBuffer(sortListPtr, &reserveSpace, &clearedDecays);
			ioNumDecays -= clearedDecays;
			
			if (reserveSpace == maxReserveSpace) {
				/* Reserve buffer is unused -- check if can read in more */
				estNumDecays = tmsortEstimateReadDecays();
				ioNumDecays = PHGMATH_Max(ioNumDecays, estNumDecays);
			}
			
			/* Refill the sort buffer */
			gotData = tmsortGetNextBatch(historyFile, sortListPtr, 
											&reserveSpace, &ioNumDecays);
			numDecaysRead += ioNumDecays;
		} while (gotData);
		
		if (LbSortNext(sortListPtr, minSortItem, &listItemPtr, 
								&numPhotons, (void **)&nextDecayPtr)) {
			/* Write out the remaining decays in the sorted list */
			do {
				/* Write out the decay */
				if (! tmsortWriteDecay(&sortFileHk, nextDecayPtr,
						(PHG_DetectedPhoton *)(nextDecayPtr+1), numPhotons)) {
					sprintf(tmsortErrStr, "Got failure writing sort decay to disk.");
					ErStFileError(tmsortErrStr);
					goto FAIL;
				}
				numDecaysWritten++;
				
				/* Remove the decay from the sorted list */
				if (! LbSortDelete(sortListPtr, listItemPtr)) {
					/* Severe problems */
					sprintf(tmsortErrStr, "Got failure removing sort decay from sorted list.");
					ErStGeneric(tmsortErrStr);
					goto FAIL;
				}
				
				/* Clear the decay data's buffer spot */
				tmsortRemoveDecay(nextDecayPtr, numPhotons);
				
			} while (LbSortNext(sortListPtr, minSortItem, &listItemPtr, 
												&numPhotons, (void **)&nextDecayPtr));
			
			/* Clear out the reserve buffer as much as possible */
			tmsortClearReserveBuffer(sortListPtr, &reserveSpace, &clearedDecays);
		}
		
		if (LbSortFind(sortListPtr, (LbSortKeyPtr)&tmsortMinKey, &minSortItem)) {
			if (! LbSortNext(sortListPtr, minSortItem, &listItemPtr, 
										&numPhotons, (void **)&nextDecayPtr)) {
				/* No data items remain in list */
				
				/* Input file has been completely sorted (for phase I) */
				allDone = true;
			}
			else {
				/* Need to start a new file */
			}
		}
		else {
			/* Sorted list is incorrectly empty!! (missing null item) */
			/* Shouldn't ever happen */
			
			/* Assume file has been completely sorted (for phase I) */
			allDone = true;
		}
		
		/* Copy the original file header data to the new file (it is desirable 
			that a sorted file have an identical header as its unsorted copy) */
		memcpy(sortFileHk.headerHk.headerData,
				tmsortHistHeaderHk.headerData, tmsortHistHeaderHk.headerSize);
		
		if (! tmsortPreserveHeader) {
			/* Add the timesorted indicator */
			tmsortHdrParams.H.isTimeSorted = true;
			if (LbHdrStElem(&(sortFileHk.headerHk), HDR_HISTORY_FILE_IS_SORTED_ID,
					sizeof(tmsortHdrParams.H.isTimeSorted),
					(void *)&(tmsortHdrParams.H.isTimeSorted)) == false){
				
				sprintf(tmsortErrStr,"Unable to set header 'timesorted indicator' parameter.");
				ErStFileError(tmsortErrStr);
				goto FAIL;
			}
		}
		
		/* Save the header to the sorted output file */
		if (PhgHdrUpHeader(NULL, &tmsortHdrParams, &(sortFileHk.headerHk)) == false) {
			sprintf(tmsortErrStr, "Unable to write header to sorted history file\n'%s'.",
				tmsortSortName);
			ErStFileError(tmsortErrStr);
			goto FAIL;
		}
		
		/* Accumulate the size of the sorted output file */
		if (fseek(sortFileHk.histFile, 0, SEEK_END) == 0) {
			if (ftell(sortFileHk.histFile) > 0) {
				sortedFilesSize += ftell(sortFileHk.histFile);
			}
		}
		
		/* Close the output file */
		/* NOTE:  Can't use PhoHFileClose because it restores the new header */
		PhgHdrFrHeader(&sortFileHk.headerHk);
		if (fclose(sortFileHk.histFile) != 0) {
			sortFileHk.histFile = NULL;
			sprintf(tmsortErrStr,"Unable to close sorted file named:\n"
				"'%s'\n"
				" (tmsortSort)",
				tmsortSortName);
			ErStFileError(tmsortErrStr);
			goto FAIL;
		}
		sortFileHk.histFile = NULL;
		
		/* Check for a progress message display */
		if (LbTmStopTiming(&clockTiming, &realSecs, &cpuSecs)) {
			accumSecs += (LbUsFourByte) realSecs;
			accumCPUSecs += cpuSecs;
			
			/* Restart the timer */
			LbTmStartTiming(&clockTiming);
			
			if ((accumSecs >= nextAccumSecs) && 	/* use longer of two */
						(sortFileNumber >= nextFileDisplay)) {
				
				accumSecs = 0;	/* easier than incrementing to total */
				nextFileDisplay += kTSProgressFiles;	/* every xx sort files */
				
				/* Get current time */
				time(&curTime);
				
				/* Convert time to string format */
				strftime(timeStr, 31, "(%H:%M  %e-%b-%y)", localtime(&curTime));
				
				/* Display a progress message */
				LbInPrintf(" %d subfiles created (%d MB),  CPU seconds = %3.0f  %s.\n",
					sortFileNumber, (int)(sortedFilesSize/1024/1024), 
					accumCPUSecs, timeStr);
			}
		}
	} while (! allDone);
	
	*numSortedFiles = sortFileNumber;
	
	/* Display the 100% progress message */
	if (LbTmStopTiming(&clockTiming, &realSecs, &cpuSecs)) {
		accumCPUSecs += cpuSecs;
		
		/* Get current time */
		time(&curTime);
		
		/* Convert time to string format */
		strftime(timeStr, 31, "(%H:%M  %e-%b-%y)", localtime(&curTime));
		
		/* Display a progress message */
		LbInPrintf(" %d subfiles created (%d MB),  CPU seconds = %3.0f  %s.\n",
			sortFileNumber, (int)(sortedFilesSize/1024/1024), 
			accumCPUSecs, timeStr);
	}
	
	/* Close the history file file */
	fclose(historyFile);
	
	okay = true;
	
	
	FAIL:;
	/* NOTE:  In case of failure, files may be left unclosed; therefore,
		it is best to completely exit the program on failure and let the 
		system clean up these files. */
	
	/* Free data structures */
	curBlockPtr = tmsortTempDecayBlocksList;
	while (curBlockPtr != NULL) {
		if (curBlockPtr->thisDataBlock != NULL) {
			/* Delete the data block */
			LbMmFree((void **) &curBlockPtr->thisDataBlock);
		}
		tmsortTempDecayBlocksList = curBlockPtr->nextSortBlock;
		LbMmFree((void **) &curBlockPtr);
	}
	for (i=0; i<=PHG_MAX_DETECTED_PHOTONS; i++) {
		if (tmsortDecaySlots[i] != NULL) {
			LbMmFree((void **) &tmsortDecaySlots[i]);
		}
	}
	LbSortDisposeList(&sortListPtr);
	
	return (okay);
}


/*********************************************************************************
*
*	Name:			tmsortMerge
*
*	Summary:		Merge multiple, smaller sorted files into a single sorted history file.
*
*	Arguments:
*		LbUsFourByte	 	sortFilesToMerge		- The # of files to merge.
*		LbUsFourByte	 	*finalFileNum			- The index of the final file.
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean tmsortMerge(LbUsFourByte sortFilesToMerge, LbUsFourByte *finalFileNum)
{
	Boolean				okay = false;				/* Function return value */
	LbUsFourByte		maxSimulSortFiles;			/* Maximum number of open sort files */
	LbUsFourByte		mergeFileNumber;			/* File index of a merged output file */
	LbUsFourByte		simulSortFiles;				/* Number of files to merge at one time */
	LbUsFourByte		i;							/* Generic index variable */
	PhoHFileHdrTy		sortHdrParams;				/* Input header */
	LbHdrHkTy			sortHeaderHk;				/* Hook to sort file header */
	LbUsFourByte		mergedFiles;				/* Current count of files merged */
	LbUsFourByte		filesToMerge;				/* Current total count of files to merge */
	LbUsFourByte		curNumSortFiles;			/* Current count of files to merge */
	time_t				curTime;					/* Current wall clock time */
	char				timeStr[32];				/* String for printing the time */
	
	
	/* Set maximum number of open files */
	/* Since allowed main buffer minimum is based on a minimal size for 3% of the 
		buffer, allow 100/3 = 33 maximum files--the reason is that the buffer 
		will be divided into this many sub-buffers. */
	maxSimulSortFiles = 33;
	
	/* Create an array of files for each file to be merged */
	if ((tmsortMergeFiles = (MergeFileData *) 
			LbMmAlloc(maxSimulSortFiles*tmsortMergeFileSize)) == NULL) {
		sprintf(tmsortErrStr, "Unable to allocate memory:  tmsortMergeFiles.");
		ErStGeneric(tmsortErrStr);
		goto FAIL;
	}
	
	/* Open an output file */
	mergeFileNumber = sortFilesToMerge + 1;
	sprintf(tmsortSortName, 
			"%s%ld", tmsortSortNameBase, (unsigned long)mergeFileNumber);
	if (PhoHFileCreate(tmsortSortName, "", 
			tmsortHdrParams.H.HdrKind, &tmsortMergeFileHk) == false) {
		sprintf(tmsortErrStr,"Unable to create sorted file named:\n"
			"'%s'\n"
			" (tmsortMerge)",
			tmsortSortName);
		ErStFileError(tmsortErrStr);
		goto FAIL;
	}
	
	/* Open as many simultaneous input files as possible, up to the maximum */
	simulSortFiles = 0;
	for (i=0; i<sortFilesToMerge; i++) {
		/* Try to open the file */
		sprintf(tmsortSortName, "%s%ld", tmsortSortNameBase, (unsigned long)(i+1));
		if ((tmsortMergeFiles[i].theFile = LbFlFileOpen(tmsortSortName, "rb")) != NULL) {
			simulSortFiles++;
			
			/* Read in the header (but ignore it) */
			if (PhgHdrGtParams(tmsortMergeFiles[i].theFile, &sortHdrParams, &sortHeaderHk) == false){
				sprintf(tmsortErrStr, "Unable to read input sort file header\n'%s'.",
					tmsortSortName);
				ErStFileError(tmsortErrStr);
				goto FAIL;
			}
		}
		else {
			/* No more can be opened */
			break;
		}
		
		if (simulSortFiles == maxSimulSortFiles) {
			/* Have reached the maximum allowable */
			break;
		}
	}
	
	/* Convert the main buffer to a merge buffer */
	tmsortMergeBufferPtr = (DecayMergeData *)tmsortSortBufferPtr;
	
	/* Merge the files */
	mergedFiles = 0;
	filesToMerge = sortFilesToMerge;
	curNumSortFiles = simulSortFiles;
	while (curNumSortFiles > 1) {
		/* Get current time */
		time(&curTime);
		
		/* Convert time to string format */
		strftime(timeStr, 31, "(%H:%M  %e-%b-%y)", localtime(&curTime));
		
		/* Display a progress message */
		LbInPrintf(" Merging subfiles %ld-%ld,  %s.\n", 
						(unsigned long)(mergedFiles+1), 
						(unsigned long)(mergedFiles+curNumSortFiles), timeStr);
		
		/* Merge the current set of files */
		if (! tmsortMergeFileBatch(curNumSortFiles)) {
			/* Error in merging files; should already have been reported */
			goto FAIL;
		}
		
		/* Close up the output file */
		{
			/* Copy the original file header data to the new file (it is desirable 
				that a sorted file have an identical header as its unsorted copy) */
			memcpy(tmsortMergeFileHk.headerHk.headerData,
					tmsortHistHeaderHk.headerData, tmsortHistHeaderHk.headerSize);
			
			if (! tmsortPreserveHeader) {
				/* Add the timesorted indicator */
				tmsortHdrParams.H.isTimeSorted = true;
				if (LbHdrStElem(&(tmsortMergeFileHk.headerHk), HDR_HISTORY_FILE_IS_SORTED_ID,
						sizeof(tmsortHdrParams.H.isTimeSorted),
						(void *)&(tmsortHdrParams.H.isTimeSorted)) == false){
					
					sprintf(tmsortErrStr,"Unable to set header 'timesorted indicator' parameter.");
					ErStFileError(tmsortErrStr);
					goto FAIL;
				}
			}
			
			/* Save the original header to the merged output file */
			if (PhgHdrUpHeader(NULL, &tmsortHdrParams, &(tmsortMergeFileHk.headerHk)) == false) {
				sprintf(tmsortSortName, 
					"%s%ld", tmsortSortNameBase, (unsigned long)mergeFileNumber);
				sprintf(tmsortErrStr, "Unable to write header to merged history file\n'%s'.",
					tmsortSortName);
				ErStFileError(tmsortErrStr);
				goto FAIL;
			}
			
			/* Close the output file */
			/* NOTE:  Can't use PhoHFileClose because it restores the new header */
			PhgHdrFrHeader(&tmsortMergeFileHk.headerHk);
			if (fclose(tmsortMergeFileHk.histFile) != 0) {
				tmsortMergeFileHk.histFile = NULL;
				sprintf(tmsortSortName, 
					"%s%ld", tmsortSortNameBase, (unsigned long)mergeFileNumber);
				sprintf(tmsortErrStr,"Unable to close merged file named:\n"
					"'%s'\n"
					" (tmsortMerge)",
					tmsortSortName);
				ErStFileError(tmsortErrStr);
				goto FAIL;
			}
			tmsortMergeFileHk.histFile = NULL;
		}
		
		/* Close and delete the merged input files */
		for (i=0; i<curNumSortFiles; i++) {
			/* Close the file (ignore errors) */
			fclose(tmsortMergeFiles[i].theFile);
			tmsortMergeFiles[i].theFile = NULL;
			
#ifndef TMSORT_KEEP_TEMP_FILES
			/* Delete the file (report errors, but keep going) */
			sprintf(tmsortSortName, 
				"%s%ld", tmsortSortNameBase, (unsigned long)(mergedFiles+i+1));
			if (remove(tmsortSortName) != 0) {
				sprintf(tmsortErrStr, "Unable to delete sort file\n'%s'.",
					tmsortSortName);
				ErStFileError(tmsortErrStr);
			}
#endif
		}
		
		/* Update the counts of files to merge */
		mergedFiles += curNumSortFiles;
		filesToMerge = filesToMerge - curNumSortFiles + 1;
		curNumSortFiles = PHGMATH_Min(simulSortFiles, filesToMerge);
		
		if (curNumSortFiles > 1) {
			/* Set up for the next merge batch */
			
			/* Open an output file */
			mergeFileNumber++;
			sprintf(tmsortSortName, 
				"%s%ld", tmsortSortNameBase, (unsigned long)mergeFileNumber);
			if (PhoHFileCreate(tmsortSortName, "", 
					tmsortHdrParams.H.HdrKind, &tmsortMergeFileHk) == false) {
				sprintf(tmsortErrStr,"Unable to create sorted file named:\n"
					"'%s'\n"
					" (tmsortMerge)",
					tmsortSortName);
				ErStFileError(tmsortErrStr);
				goto FAIL;
			}
			
			/* Open the input files */
			for (i=0; i<curNumSortFiles; i++) {
				/* Try to open the file */
				sprintf(tmsortSortName, 
					"%s%ld", tmsortSortNameBase, (unsigned long)(mergedFiles+i+1));
				if ((tmsortMergeFiles[i].theFile = LbFlFileOpen(tmsortSortName, "rb")) != NULL) {
					/* Read in the header (but ignore it) */
					if (PhgHdrGtParams(tmsortMergeFiles[i].theFile, &sortHdrParams, &sortHeaderHk) == false){
						sprintf(tmsortErrStr, "Unable to read input sort file header\n'%s'.",
							tmsortSortName);
						ErStFileError(tmsortErrStr);
						goto FAIL;
					}
				}
			}
		}
	}
	
	*finalFileNum = mergeFileNumber;
	
	/* Free memory */
	if (tmsortMergeFiles != NULL) {
		LbMmFree((void **) &tmsortMergeFiles);
	}
	
	okay = true;
	
	
	FAIL:;
	/* NOTE:  In case of failure, files may be left unclosed; therefore,
		it is best to completely exit the program on failure and let the 
		system clean up these files.  This also requires tmsortMergeFiles 
		be left unfreed. */
	
	return (okay);
}


/*********************************************************************************
*
*	Name:			tmsortMergeFileBatch
*
*	Summary:		Merge numSortFiles sorted files into a single sorted file.
*
*	Arguments:
*		LbUsFourByte	 	numSortFiles			- The # of files to merge.
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean tmsortMergeFileBatch(LbUsFourByte numSortFiles)
{
	Boolean				okay = false;				/* Function return value */
	LbSortListPtr		mergeListPtr;				/* Pointer to merge sorted list */
	LbUsFourByte		i;							/* Generic index variable */
	LbUsFourByte		bufferUsed;					/* Amount of main buffer used for current sort file */
	LbUsFourByte		bufOffset;					/* Start of buffer for current sort file */
	DecayMergeData		*mergeBufferPtr;			/* Pointer into current sorted merge buffer */
	EventTy				eventType;					/* Type of current event */
	LbUsFourByte		numDecaysRead;				/* Number of decays read in so far */
	LbUsFourByte		numDecaysWritten;			/* Number of decays written out so far */
	TmSortKeyType		minSortKey;					/* Minimum sort key */
	Boolean				allDone;					/* Files have been merged */
	LbSortItemPtr		minSortItem;				/* Item containing last minimum key */
	LbSortItemPtr		mergeItemPtr;				/* Item containing current minimum decay */
	LbUsFourByte		numPhotons;					/* Number of photons in the decay */
	PHG_Decay			*nextDecayPtr;				/* Current output decay data */
	TmSortKeyType		nextSortKey;				/* Next minimum sort key */
	LbSortItemPtr		nextItemPtr;				/* Item containing next minimum decay */
	LbUsFourByte		fileIndex;					/* Current sort file index number */
	PHG_Decay			*mergeDecayPtr;				/* Ptr to current minimum decay data */
	LbUsFourByte		remainingDecays;			/* Decays still in current sort buffer */
	LbUsFourByte		dataSize;					/* Size of current decay data */
	
	
	/* Make a new sorted list */
	if ((mergeListPtr = LbSortNewList(numSortFiles+1, 
							(LbSortKeyPtr)&tmsortMinKey, sizeof(TmSortKeyType), 
							tmsortKeyComparator)) == NULL) {
		/* Severe problems */
		sprintf(tmsortErrStr, "Got failure creating sorted list.");
		ErStGeneric(tmsortErrStr);
		goto FAIL;
	}
	
	/* Add an initial tmsortMinKey item (to increase access speed) */
	if (! LbSortInsert(mergeListPtr, (LbSortKeyPtr)&tmsortMinKey, 0, NULL) ) {
		/* Severe problems */
		sprintf(tmsortErrStr, "Got failure adding first sorted list item.");
		ErStGeneric(tmsortErrStr);
		goto FAIL;
	}
	
	/* Divide the main buffer into numSortFiles sections */
	bufferUsed = 0;
	tmsortMergeFiles[0].bufStart = tmsortMergeBufferPtr;
	for (i=1; i<numSortFiles; i++) {
		bufOffset = (LbUsFourByte)
			((float)i/numSortFiles * tmsortDataBufferSize);
		tmsortMergeFiles[i].bufStart = (DecayMergeData *)
			((LbUsOneByte *)tmsortMergeBufferPtr + bufOffset);
		tmsortMergeFiles[i-1].bufSize = bufOffset - bufferUsed;
		bufferUsed = bufOffset;
	}
	tmsortMergeFiles[numSortFiles-1].bufSize = 
		tmsortDataBufferSize - bufferUsed;
	
	/* Read in the first decay for each file */
	for (i=0; i<numSortFiles; i++) {
		if (! tmsortCustomFile) {
			/* But only for standard files */
			mergeBufferPtr = tmsortMergeFiles[i].bufStart;
			tmsortDecayPending = false;
			eventType = tmsortReadEvent(tmsortMergeFiles[i].theFile, 
							&(tmsortMergeFiles[i].nextDecay), 
							(PHG_DetectedPhoton *) &(mergeBufferPtr->decayData));
			if (eventType != Decay) {
				ErStGeneric("Expected first event to be decay, and it wasn't.");
				goto FAIL;
			}
		}
		
		/* Mark the file as ready for reading */
		tmsortMergeFiles[i].nextDecayPending = true;
	}
	
	/* Initialize decay counts */
	numDecaysRead = 0;
	numDecaysWritten = 0;
	
	/* Fill the buffers from each file */
	for (i=0; i<numSortFiles; i++) {
		tmsortFillMergeBuffer(i);
		numDecaysRead += tmsortMergeFiles[i].bufDecays;
	}
	
	/* Insert the start of each buffer into the sorted list */
	for (i=0; i<numSortFiles; i++) {
		if (tmsortMergeFiles[i].bufDecays > 0) {
			tmsortGetDecaySortKey(&(tmsortMergeFiles[i].bufStart->decayData), 
									&minSortKey);
			if (! LbSortInsert(mergeListPtr, 
						(LbSortKeyPtr)&minSortKey, 
						i,	/* File index */
						(PHG_Decay *)(tmsortMergeFiles[i].bufStart)
							/* Buffer start */
					)) {
				
				/* Severe problems; but just keep going */
				sprintf(tmsortErrStr, "Got failure inserting decay into sorted list.");
				ErStGeneric(tmsortErrStr);
			}
		}
	}
	
	/* Find the special, first item in the sorted list */
	if (! LbSortFirst(mergeListPtr, &minSortItem, &numPhotons, (void **)&nextDecayPtr)) {
		/* Severe problems; but should never occur */
		sprintf(tmsortErrStr, "Couldn't find null merge file sort item.");
		ErStGeneric(tmsortErrStr);
		goto FAIL;
	}
	
	/* Do the merge */
	allDone = false;
	while (! allDone) {
		/* Find the minimal decay data buffer in the sorted list */
		if (LbSortNext(mergeListPtr, minSortItem, &mergeItemPtr, 
									&numPhotons, (void **)&nextDecayPtr)) {
			if (LbSortNext(mergeListPtr, mergeItemPtr, &nextItemPtr, 
										&numPhotons, (void **)&nextDecayPtr)) {
				/* Note the next smallest buffer data */
				LbSortGetItemSortKey(nextItemPtr, 
						(LbSortKeyPtr)&nextSortKey, sizeof(TmSortKeyType));
			}
			else {
				/* Only one file left to merge */
				nextSortKey = tmsortMaxKey;
			}
		}
		else {
			/* Nothing left in the list--all files merged */
			allDone = true;
		}
		
		if (! allDone) {
			/* Set up at the minimum decay */
			fileIndex = LbSortGetItemIndex(mergeItemPtr);
			mergeBufferPtr = (DecayMergeData *)LbSortGetItemDataPtr(mergeItemPtr);
			mergeDecayPtr = &(mergeBufferPtr->decayData);
			remainingDecays = tmsortMergeFiles[fileIndex].bufDecays;
			
			/* Write out all minimal decays in this buffer */
			do {
				/* Write out the minimum decay */
				if (! tmsortWriteDecay(&tmsortMergeFileHk, mergeDecayPtr,
						(PHG_DetectedPhoton *)(mergeDecayPtr+1), 
						mergeBufferPtr->numPhotons)) {
					sprintf(tmsortErrStr, "Got failure writing merge decay to disk.");
					ErStFileError(tmsortErrStr);
					goto FAIL;
				}
				numDecaysWritten++;
				dataSize = tmsortMergeDecaySize + 
					mergeBufferPtr->numPhotons*tmsortPhotonSize;
				mergeBufferPtr = (DecayMergeData *)
					((LbUsOneByte *)mergeBufferPtr + dataSize);
				mergeDecayPtr = &(mergeBufferPtr->decayData);
				remainingDecays--;
				
				if (remainingDecays > 0) {
					/* Get the next sort key in the buffer */
					tmsortGetDecaySortKey(mergeDecayPtr, &minSortKey);
				}
			} while ((remainingDecays > 0) && (minSortKey <= nextSortKey));
			
			if (remainingDecays > 0) {
				/* Save this buffer's decay count */
				tmsortMergeFiles[fileIndex].bufDecays = remainingDecays;
			}
			else {
				/* Reload the buffer from the input sort file */
				tmsortFillMergeBuffer(fileIndex);
				mergeBufferPtr = tmsortMergeFiles[fileIndex].bufStart;
				remainingDecays = tmsortMergeFiles[fileIndex].bufDecays;
				numDecaysRead += remainingDecays;
				mergeDecayPtr = &(mergeBufferPtr->decayData);
				tmsortGetDecaySortKey(mergeDecayPtr, &minSortKey);
			}
			
			/* Readjust the sorted list */
			if (! LbSortDelete(mergeListPtr, mergeItemPtr)) {
				/* Severe problems; but just keep going */
				sprintf(tmsortErrStr, "Got failure removing merge decay from sorted list.");
				ErStGeneric(tmsortErrStr);
			}
			if (remainingDecays > 0) {
				if (! LbSortInsert(mergeListPtr, 
							(LbSortKeyPtr)&minSortKey, 
							fileIndex, (PHG_Decay *)mergeBufferPtr)) {
					/* Severe problems; but just keep going */
					sprintf(tmsortErrStr, "Got failure adding merge decay to sorted list.");
					ErStGeneric(tmsortErrStr);
				}
			}
		}
	}
	
	okay = true;
	
	
	FAIL:;
	
	/* Free memory */
	LbSortDisposeList(&mergeListPtr);
	
	return (okay);
}


/*********************************************************************************
*
*	Name:			tmsortGetFirstBatch
*
*	Summary:		Fill the data buffer as much as possible with new events 
*						from the given input file.
*					Record the locations and photon counts for later sorting.
*					Record the distribution of decay slots by photon counts.
*
*	Arguments:
*		FILE					*historyFile	- History file to read from.
*		LbUsFourByte			bufferSize		- Space available in buffer.
*		LbUsFourByte			*numDecays		- Number of decays read in.
*
*	Function return: True if successful; false if not.
*
*********************************************************************************/
Boolean tmsortGetFirstBatch(FILE *historyFile, 
							LbUsFourByte bufferSize, LbUsFourByte *numDecays)
{
	Boolean				result;				/* Function return value */
	PHG_Decay			*bufPtr;			/* Local pointer into data buffer */
	LbUsFourByte		availSpace;			/* Size of remaining buffer space */
	LbUsFourByte		minSpaceReqd;		/* Guaranteed empty space */
	DecaySortBlock		*curBlockPtr;		/* Pointer to the current block */
	LbUsFourByte		curBlockCount;		/* Data count of current block */
	LbUsFourByte		i;					/* Index through number of photons */
	LbUsFourByte		dataSize;			/* Size of single decay data */
	LbUsFourByte		numPhotons;			/* Photons associated with a decay */
	DecaySortData		*sortDataPtr;		/* Pointer to the current data */
	
	
	result = true;
	bufPtr = tmsortSortBufferPtr;
	availSpace = bufferSize;
	minSpaceReqd = 2*tmsortDecaySize + PHG_MAX_DETECTED_PHOTONS*tmsortPhotonSize;
	curBlockPtr = tmsortTempDecayBlocksList;
	curBlockCount = 0;
	
	/* Zero out the decay counts (for use as temporary counters) */
	for (i=0; i<=PHG_MAX_DETECTED_PHOTONS; i++) {
		tmsortDecaySlotCounts[i] = 0;
	}
	
	/* Fill the buffer as much as possible */
	while (result && (availSpace >= minSpaceReqd)) {
		if (tmsortReadDecay(historyFile, bufPtr, &dataSize, &numPhotons)) {
			/* A decay was read; reduce the available buffer space */
			availSpace -= dataSize;
			
			/* Record the data for later sorting */
			if (curBlockCount >= tmsortBlockSize) {
				/* Need to allocate a new block */
				if ((curBlockPtr->nextSortBlock = (DecaySortBlock *) 
						LbMmAlloc(sizeof(DecaySortBlock))) == NULL) {
					sprintf(tmsortErrStr, "Unable to allocate memory:  DecaySortBlock.");
					ErStGeneric(tmsortErrStr);
					result = false;
					break;
				}
				curBlockPtr = curBlockPtr->nextSortBlock;
				curBlockPtr->nextSortBlock = NULL;
				if ((curBlockPtr->thisDataBlock = (DecaySortData *) 
						LbMmAlloc(tmsortBlockSize*sizeof(DecaySortData))) == NULL) {
					sprintf(tmsortErrStr, "Unable to allocate memory:  DecaySortData.");
					ErStGeneric(tmsortErrStr);
					result = false;
					break;
				}
				curBlockCount = 0;
			}
			sortDataPtr = curBlockPtr->thisDataBlock + curBlockCount;
			sortDataPtr->numPhotons = numPhotons;
			sortDataPtr->decayDataPtr = bufPtr;
			curBlockCount++;
			
			/* Record the decay data by photon count */
			tmsortDecaySlotCounts[numPhotons]++;
			
			bufPtr = (PHG_Decay *)((LbUsOneByte *)bufPtr + dataSize);
		}
		else {
			/* Ran out of decay data */
			
			/* Record the end of data, if necessary */
			if (curBlockCount < tmsortBlockSize) {
				sortDataPtr = curBlockPtr->thisDataBlock + curBlockCount;
				sortDataPtr->decayDataPtr = NULL;
			}
			
			break;
		}
	}
	
	/* Add up the total decays and record the maximum photon counts */
	*numDecays = 0;
	for (i=0; i<=PHG_MAX_DETECTED_PHOTONS; i++) {
		*numDecays += tmsortDecaySlotCounts[i];
		
		tmsortDecaySlotMaxs[i] = tmsortDecaySlotCounts[i];
		tmsortDecaySlotCounts[i] = 0;	/* Restore to 0 */
	}
	
	return (result);
}


/*********************************************************************************
*
*	Name:			tmsortGetNextBatch
*
*	Summary:		Fill the sort buffer as much as possible with numDecays 
*						from the given input file.
*					Record the decay positions in the record array.
*
*	Arguments:
*		FILE				*historyFile		- History file to read from.
*		LbSortListPtr		sortListPtr			- The sort list.
*		LbUsFourByte	 	*reserveSpace		- Size of reserve buffer space.
*		LbUsFourByte	 	*numDecays			- The # of decays to read, 
*													and actually read.
*
*	Function return: True if data read, false if none remained or could be read.
*
*********************************************************************************/
Boolean tmsortGetNextBatch(FILE *historyFile, LbSortListPtr sortListPtr, 
							LbUsFourByte *reserveSpace, LbUsFourByte *numDecays)
{
	Boolean				dataRead;			/* Function return value */
	LbUsFourByte		availSpace;			/* Size of remaining buffer space */
	LbUsFourByte		minSpaceReqd;		/* Minimum reserve buffer size */
	LbUsFourByte		decaysRead;			/* Number of decays read in */
	PHG_Decay			*localBufferPtr;	/* Local ptr into reserve buffer */
	LbUsFourByte		i;					/* Index through slot counts */
	LbUsFourByte		d;					/* Index through numDecays */
	LbUsFourByte		dataSize;			/* Size of decay */
	LbUsFourByte		numPhotons;			/* Number of photons in decay */
	PHG_Decay			*decayLocation;		/* Where decay is stored */
	
	
	dataRead = false;
	availSpace = *reserveSpace;
	minSpaceReqd = 2*tmsortDecaySize + PHG_MAX_DETECTED_PHOTONS*tmsortPhotonSize;
	
	decaysRead = 0;
	if ((availSpace >= minSpaceReqd) && (*numDecays > 0)) {
		localBufferPtr = tmsortReserveBufferPtr;
		
		/* Clear the next read distribution */
		tmsortCurDist = (tmsortCurDist + 1) % TMSORT_DIST_CT;
		for (i=0; i<=PHG_MAX_DETECTED_PHOTONS; i++) {
			tmsortDecaySlotDists[tmsortCurDist][i] = 0;
		}
		
		for (d=0; d<*numDecays; d++) {
			if (tmsortReadDecay(historyFile, localBufferPtr, &dataSize, &numPhotons)) {
				dataRead = true;
				decaysRead++;
				tmsortDecaySlotDists[tmsortCurDist][numPhotons]++;
				
				/* Try to fit the data into the buffer memory */
				if (! tmsortStoreDecay(localBufferPtr, numPhotons, &decayLocation)) {
					/* Couldn't find a spot, so use reserve buffer space for it */
					localBufferPtr = (PHG_Decay *)((LbUsOneByte *)localBufferPtr + dataSize);
					availSpace -= dataSize;
				}
				
				/* Add the data to the sorted list */
				if (! tmsortSortDecay(sortListPtr, decayLocation, numPhotons)) {
					sprintf(tmsortErrStr, "Unable to sort data; all results are invalid.");
					ErStGeneric(tmsortErrStr);
					break;
				}
				
				/* Check for low memory */
				if (availSpace < minSpaceReqd) {
					break;
				}
			}
			else {
				/* Ran out of decay data */
				break;
			}
		}
		
		/* Update reserve space */
		*reserveSpace = availSpace;
		
		/* Update reserve buffer ptr */
		tmsortReserveBufferPtr = localBufferPtr;
	}
	
	/* Update numDecays */
	*numDecays = decaysRead;
	
	
	return (dataRead);
}


/*********************************************************************************
*
*	Name:			tmsortClearReserveBuffer
*
*	Summary:		Clear out the reserve buffer as much as possible.
*
*	Arguments:
*		LbSortListPtr		sortListPtr		- The sort list.
*		LbUsFourByte	 	*reserveSpace	- Size of reserve buffer.
*		LbUsFourByte	 	*clearedDecays	- Number of decays removed from reserve buffer.
*
*	Function return: None.
*
*********************************************************************************/
void tmsortClearReserveBuffer(LbSortListPtr sortListPtr, 
						LbUsFourByte *reserveSpace, LbUsFourByte *clearedDecays)
{
	LbUsFourByte		clearedCount;		/* Number of decays cleared from reserve buffer */
	PHG_Decay			*usedResBufferPtr;	/* Pointer to end of used reserve buffer */
	PHG_Decay			*curResBufferPtr;	/* Pointer into previous reserve buffer */
	TmSortKeyType		decaySortKey;		/* Decay time sort key */
	LbSortItemPtr		theSortItemPtr;		/* Sorted list item for current decay */
	Boolean				foundDecay;			/* Result of sorted list data search */
	LbUsFourByte		numPhotons;			/* Number of photons in decay */
	PHG_Decay			*decayLocation;		/* Where decay is stored */
	TmSortKeyType		foundSortKey;		/* Sort key of theSortItemPtr */
	LbUsFourByte		dataSize;			/* Size of decay data */
	PHG_Decay			*storedLocation;	/* New location of decay data */
	
	
	/* Scan through the data in the reserve buffer */
	clearedCount = 0;
	usedResBufferPtr = tmsortResBufferStart;
	curResBufferPtr = tmsortResBufferStart;
	while (curResBufferPtr < tmsortReserveBufferPtr) {
		/* Get the decay sort key */
		tmsortGetDecaySortKey(curResBufferPtr, &decaySortKey);
		
		if (decaySortKey != tmsortMinKey) {
			/* Find the data in the sorted list */
			if (! LbSortFind(sortListPtr, (LbSortKeyPtr)&decaySortKey, &theSortItemPtr)) {
				/* Data isn't in sorted list!!  Severe error. */
				foundDecay = false;
				sprintf(tmsortErrStr, "Sort data missing; all results are invalid.");
				ErStGeneric(tmsortErrStr);
				break;
			}
			
			/* Make sure the sorted list item is the correct one */
			foundDecay = true;
			while (((PHG_Decay *)LbSortGetItemDataPtr(theSortItemPtr)) != curResBufferPtr) {
				/* Get the previous one */
				foundDecay = LbSortPrev(sortListPtr, 
					theSortItemPtr, &theSortItemPtr, &numPhotons, (void **)&decayLocation);
				if (! foundDecay) {
					/* Data isn't in sorted list!!  Severe error. */
					sprintf(tmsortErrStr, "Sort data missing; all results are invalid.");
					ErStGeneric(tmsortErrStr);
					break;
				}
				LbSortGetItemSortKey(theSortItemPtr, 
						(LbSortKeyPtr)&foundSortKey, sizeof(TmSortKeyType));
				if (foundSortKey != decaySortKey) {
					/* Data isn't in sorted list!!  Severe error. */
					foundDecay = false;
					sprintf(tmsortErrStr, "Sort data missing; all results are invalid.");
					ErStGeneric(tmsortErrStr);
					break;
				}
			}
			if (! foundDecay) {
				break;
			}
			
			/* Try to store the decay in the main buffer */
			numPhotons = LbSortGetItemIndex(theSortItemPtr);
			dataSize = tmsortDecaySize + numPhotons*tmsortPhotonSize;
			if (tmsortStoreDecay(curResBufferPtr, numPhotons, &storedLocation)) {
				/* Decay data now no longer in reserve buffer */
				LbSortSetItemDataPtr(theSortItemPtr, storedLocation);
				clearedCount++;
				*reserveSpace += dataSize;
			}
			else {
				/* Couldn't find a spot, so it stays in the reserve buffer */
				if (usedResBufferPtr < curResBufferPtr) {
					/* Compress the reserve buffer */
					memcpy(usedResBufferPtr, curResBufferPtr, dataSize);
					
					/* Update the recorded data location */
					LbSortSetItemDataPtr(theSortItemPtr, usedResBufferPtr);
				}
				usedResBufferPtr = (PHG_Decay *)((LbUsOneByte *)usedResBufferPtr + dataSize);
			}
		}
		else {
			/* Decay has been written out already; data is invalid except for photon count */
			dataSize = tmsortDecaySize + 
				(*((LbUsFourByte *)((LbUsOneByte *)curResBufferPtr + TMSORT_TEMP_COUNT_OFFSET))) * 
				tmsortPhotonSize;
			*reserveSpace += dataSize;
		}
		
		curResBufferPtr = (PHG_Decay *)((LbUsOneByte *)curResBufferPtr + dataSize);
	}
	*clearedDecays = clearedCount;
	
	/* Update the start of the reserve buffer */
	tmsortReserveBufferPtr = usedResBufferPtr;
}


/*********************************************************************************
*
*	Name:			tmsortEstimateReadDecays
*
*	Summary:		Make an estimate of the maximum number of decays to read.
*
*	Arguments:		None.
*
*	Function return: The estimated number of decays that can be read.
*
*********************************************************************************/
LbUsFourByte tmsortEstimateReadDecays(void)
{
	LbUsFourByte		curReadTotal;		/* Count of current distribution decays */
	LbUsFourByte		i;					/* Index through photon counts */
	LbUsFourByte		j;					/* Index through decay distributions */
	LbUsFourByte		curReadDist[PHG_MAX_DETECTED_PHOTONS+1];
											/* Running average distribution */
	LbUsFourByte		histoCell;			/* Count in histogram single cell */
	float				bestMult;			/* Best scale multiplier of current distribution */
	LbUsFourByte		estDecays;			/* The estimated decay count */
	
	
	/* Compute the running average distribution currently being read */
	curReadTotal = 0;
	for (i=0; i<=PHG_MAX_DETECTED_PHOTONS; i++) {
		histoCell = tmsortDecaySlotDists[0][i];
		curReadDist[i] = histoCell;
		curReadTotal += histoCell;
	}
	for (j=1; j<TMSORT_DIST_CT; j++) {
		for (i=0; i<=PHG_MAX_DETECTED_PHOTONS; i++) {
			histoCell = tmsortDecaySlotDists[j][i];
			curReadDist[i] += histoCell;
			curReadTotal += histoCell;
		}
	}
	
	/* Find the largest multiplier of the current distribution that will fit 
		in the open slots distribution */
	for (i=0; i<=PHG_MAX_DETECTED_PHOTONS; i++) {
		if ((tmsortDecaySlotMaxs[i] > 0) && (curReadDist[i] > 0)) {
			/* Set best multiplier to start with */
			bestMult = (float)tmsortDecaySlotCounts[i] / 
						curReadDist[i];
			break;
		}
	}
	for (i=0; i<=PHG_MAX_DETECTED_PHOTONS; i++) {
		if ((tmsortDecaySlotMaxs[i] > 0) && (curReadDist[i] > 0)) {
			if (bestMult*curReadDist[i] > tmsortDecaySlotCounts[i]) {
				/* Have to reduce the multiplier */
				bestMult = (float)tmsortDecaySlotCounts[i] / 
							curReadDist[i];
			}
		}
	}
	
	/* Compute the estimated number to read in */
	estDecays = (LbUsFourByte)(bestMult*curReadTotal + 0.5);
	
	/* Return the estimate */
	return (estDecays);
}


/*********************************************************************************
*
*	Name:			tmsortFillMergeBuffer
*
*	Summary:		Fill the data buffer memory for the supplied sort file.
*
*	Arguments:
*		LbUsFourByte	 		fileIndex		- Index of the input sort file.
*
*	Function return: None.
*
*********************************************************************************/
void tmsortFillMergeBuffer(LbUsFourByte fileIndex)
{
	LbUsFourByte		bufDecays;			/* Number of decays read in */
	LbUsFourByte		minSpaceReqd;		/* Empty space required at buffer end */
	DecayMergeData		*bufPtr;			/* Start of fileIndex buffer */
	LbUsFourByte		curBufferSpace;		/* Remaining space in the buffer */
	FILE				*theFile;			/* Local copy of sort file */
	LbUsFourByte		dataSize;			/* Size of the decay data */
	LbUsFourByte		numPhotons;			/* Number of photons in the decay */
	LbUsFourByte		dataSpace;			/* Buffer space used by the decay */
	
	
	if (! tmsortMergeFiles[fileIndex].nextDecayPending) {
		/* Can't read anything in */
		bufDecays = 0;
	}
	else {
		if (! tmsortCustomFile) {
			/* Set up the pending decay */
			memcpy(&tmsortNextDecayEvent, 
						&(tmsortMergeFiles[fileIndex].nextDecay), tmsortDecaySize);
			tmsortDecayPending = true;
		}
		
		/* Fill the buffer */
		minSpaceReqd = tmsortMergeDecaySize + tmsortDecaySize + 
						PHG_MAX_DETECTED_PHOTONS*tmsortPhotonSize;
		bufDecays = 0;
		bufPtr = tmsortMergeFiles[fileIndex].bufStart;
		curBufferSpace = tmsortMergeFiles[fileIndex].bufSize;
		theFile = tmsortMergeFiles[fileIndex].theFile;
		while (curBufferSpace > minSpaceReqd) {
			if (tmsortReadDecay(theFile, &(bufPtr->decayData), &dataSize, &numPhotons)) {
				bufDecays++;
				bufPtr->numPhotons = numPhotons;
				dataSpace = tmsortMergeDecaySize + dataSize - tmsortDecaySize;
				bufPtr = (DecayMergeData *)((LbUsOneByte *)bufPtr + dataSpace);
				curBufferSpace -= dataSpace;
			}
			else {
				/* Nothing left in file */
				tmsortMergeFiles[fileIndex].nextDecayPending = false;
				break;
			}
		}
		
		if (! tmsortCustomFile) {
			/* Save the pending decay (may not actually exist if eof) */
			memcpy(&(tmsortMergeFiles[fileIndex].nextDecay), 
						&tmsortNextDecayEvent, tmsortDecaySize);
		}
	}
	
	/* Record the decays read */
	tmsortMergeFiles[fileIndex].bufDecays = bufDecays;
}


/*********************************************************************************
*
*	Name:			tmsortSortDecay
*
*	Summary:		Record the decay into its sorted position in the sorted list.
*
*	Arguments:
*		LbSortListPtr		sortListPtr		- The sorted list.
*		PHG_Decay			*decayPtr		- The decay data.
*		LbUsFourByte	 	numPhotons		- The number of photons in the decay.
*
*	Function return: True if successful.
*
*********************************************************************************/
Boolean tmsortSortDecay(LbSortListPtr sortListPtr, 
							PHG_Decay *decayPtr, LbUsFourByte numPhotons)
{
	Boolean				result = false;		/* Function return */
	TmSortKeyType		decaySortKey;		/* Decay time sort key */
	
	
	/* Get the decay sort key */
	tmsortGetDecaySortKey(decayPtr, &decaySortKey);
	
	/* Add the decay to the sorted list */
	if (! LbSortInsert(sortListPtr, (LbSortKeyPtr)&decaySortKey, numPhotons, decayPtr)) {
		/* Severe problems */
		sprintf(tmsortErrStr, "Got failure inserting decay into sorted list.");
		ErStGeneric(tmsortErrStr);
		result = false;
	}
	else {
		result = true;
	}
	
	return (result);
}


/*********************************************************************************
*
*	Name:			tmsortStoreDecay
*
*	Summary:		Store the decay record in the buffer memory.
*
*	Arguments:
*		PHG_Decay				*decayPtr		- The decay data to store.
*		LbUsFourByte	 		numPhotons		- Number of photons in the decay.
*		PHG_Decay				**decayLocation	- Returned decay data location.
*
*	Function return: True if decay data stored in buffer;
*						false if no place for it.
*
*********************************************************************************/
Boolean tmsortStoreDecay(PHG_Decay *decayPtr, LbUsFourByte numPhotons,
							PHG_Decay **decayLocation)
{
	Boolean				found;				/* If proper memory location found */
	PHG_Decay			*decaySlot;			/* Proper decay location */
	
	
	/* Find a memory location for the decay data */
	found = tmsortGetDecaySlot(numPhotons, &decaySlot);
	
	if (found) {
		/* Move the decay data to the located slot */
		memcpy(decaySlot, decayPtr, 
				tmsortDecaySize+numPhotons*tmsortPhotonSize);
		*decayLocation = decaySlot;
	}
	else {
		/* No available slot, so decay must remain in current buffer */
		*decayLocation = decayPtr;
	}
	
	/* Return the result */
	return (found);
}


/*********************************************************************************
*
*	Name:			tmsortRemoveDecay
*
*	Summary:		Remove the decay record from the buffer memory.
*
*	Arguments:
*		PHG_Decay				*decayLocation	- The decay data buffer location.
*		LbUsFourByte	 		numPhotons		- Number of photons in the decay.
*
*	Function return: None.
*
*********************************************************************************/
void tmsortRemoveDecay(PHG_Decay *decayLocation, LbUsFourByte numPhotons)
{
	LbUsFourByte		index;			/* Photon index into decay slot array */
	LbUsFourByte		slotCount;		/* Available decay slots */
	PHG_Decay			**slotPtr;		/* Pointer to decay data pointer */
	
	
	if (decayLocation < tmsortResBufferStart) {
		/* Record the main buffer decay location as an empty slot */
		index = numPhotons;
		slotCount = tmsortDecaySlotCounts[index];
		slotPtr = tmsortDecaySlots[index] + slotCount;
		*slotPtr = decayLocation;
		tmsortDecaySlotCounts[index]++;
	}
	else {
		/* Decay was stored in the reserve buffer; mark it as invalid */
		tmsortSetDecaySortKey(decayLocation, tmsortMinKey);
		
		/* Save the number of photons in the invalid data location */
		*((LbUsFourByte *)((LbUsOneByte *)decayLocation + TMSORT_TEMP_COUNT_OFFSET)) 
			= numPhotons;
	}
}


/*********************************************************************************
*
*	Name:			tmsortGetDecaySlot
*
*	Summary:		Find and return a memory slot for the decay.
*					Record the used slot in the open slot list.
*
*	Arguments:
*		LbUsFourByte 		numPhotons		- The number of decay photons.
*		PHG_Decay			**decaySlot		- Located memory slot.
*
*	Function return: True if slot located.
*
*********************************************************************************/
Boolean tmsortGetDecaySlot(LbUsFourByte numPhotons, PHG_Decay **decaySlot)
{
	Boolean				found;			/* Function result */
	LbUsFourByte		index;			/* Photon index into decay slot array */
	LbUsFourByte		slotCount;		/* Available decay slots */
	PHG_Decay			**slotPtr;		/* Pointer to decay data slot */
	
	
	found = false;
	
	index = numPhotons;
	slotCount = tmsortDecaySlotCounts[index];
	if (slotCount > 0) {
		/* A spot for it is available */
		slotPtr = tmsortDecaySlots[index] + slotCount-1;
		*decaySlot = *slotPtr;
		tmsortDecaySlotCounts[index]--;
		found = true;
	}
	else {
		/* No spot available for that photon count */
		found = false;
	}
	
	return (found);
}


/*********************************************************************************
*
*	Name:			tmsortKeyComparator
*
*	Summary:		Compare the two sort keys, key1 and key2.
*
*	Arguments:
*		LbSortKeyPtr		keyPtr1			- Pointer to first sort key.
*		LbSortKeyPtr		keyPtr2			- Pointer to second sort key.
*
*	Function return:  -1 if key1<key2; 0 if key1=key2; 1 if key1>key2.
*
*********************************************************************************/
LbFourByte tmsortKeyComparator(LbSortKeyPtr keyPtr1, LbSortKeyPtr keyPtr2)
{
	TmSortKeyType		sortKey1;		/* First sort key */
	TmSortKeyType		sortKey2;		/* Second sort key */
	LbFourByte			result;			/* Result of function */
	
	
	/* Compare the two sort keys and return the result */
	
	sortKey1 = *((TmSortKeyPtr) keyPtr1);
	sortKey2 = *((TmSortKeyPtr) keyPtr2);
	
	if (sortKey1 < sortKey2) {
		result = -1;
	}
	else if (sortKey1 == sortKey2) {
		result = 0;
	}
	else {
		result = 1;
	}
	
	return (result);
}


/*********************************************************************************
*
*	Name:			tmsortGetDecaySortKey
*
*	Summary:		Return the sort key for the supplied decay.
*
*	Arguments:
*		PHG_Decay			*decayPtr		- Pointer to decay data.
*		TmSortKeyPtr		decaySortKeyPtr	- Pointer to returned decay sort key.
*
*	Function return: None.
*
*********************************************************************************/
void tmsortGetDecaySortKey(PHG_Decay *decayPtr, TmSortKeyPtr decaySortKeyPtr)
{
	/* Find and return the sort key */
	/* NOTE:  If you change this be sure to also change the Set Key function */
	/* 		Also check tmsortMinKey, tmsortMaxKey, and tmsortKeyComparator */
	
	/* ### Correct as appropriate */
	/**decaySortKeyPtr = ((PHG_DetectedPhoton *)(decayPtr+1))->time_since_creation;*/
	*decaySortKeyPtr = decayPtr->decayTime;
}


/*********************************************************************************
*
*	Name:			tmsortSetDecaySortKey
*
*	Summary:		Set the sort key for the supplied decay.
*
*	Arguments:
*		PHG_Decay			*decayPtr		- Pointer to decay data.
*		TmSortKeyType		decaySortKey	- New decay sort key.
*
*	Function return: None.
*
*********************************************************************************/
void tmsortSetDecaySortKey(PHG_Decay *decayPtr, TmSortKeyType decaySortKey)
{
	/* Set the new sort key */
	/* NOTE:  If you change this be sure to also change the Get Key function */
	/* 		Also check tmsortMinKey, tmsortMaxKey, and tmsortKeyComparator */
	
	/* ### Correct as appropriate */
	/*((PHG_DetectedPhoton *)(decayPtr+1))->time_since_creation = decaySortKey;*/
	decayPtr->decayTime = decaySortKey;
}


/**********************
*
*	Name:		tmsortReadDecay
*
*	Summary:	Read the next decay event data from the input file.
*
*	Arguments:
*		FILE				*historyFile	- The history file.
*		PHG_Decay 			*decayPtr		- Storage for decay data.
*			NOTE:  Must be guaranteed to hold complete decay data + decay event.
*		LbUsFourByte 		*dataSize		- Returned byte size of decay data.
*		LbUsFourByte 		*numPhotons		- Returned number of decay photons.
*
*	Result:		True if data read in.
*
***********************/
Boolean tmsortReadDecay(FILE *historyFile, 
						PHG_Decay *decayPtr, 
						LbUsFourByte *dataSize,
						LbUsFourByte *numPhotons)
{
	LbUsOneByte			*bufPtr;		/* Buffer to be filled (as bytes) */
	Boolean				allDone;		/* Indicates all decay photons read */
	EventTy				eventType;		/* Type of each event read in */
	PHG_TrackingPhoton	*photonPtr;		/* Pointer to custom photons */
	LbUsOneByte			numBlue;		/* Number of blue photons for this decay */
	LbUsFourByte		bIndex;			/* Current blue photon for this decay */
	LbUsOneByte			numPink;		/* Number of pink photons for this decay */
	LbUsFourByte		pIndex;			/* Current pink photon for this decay */
	Boolean				photonAccepted;	/* Ignored */
	
	
	if (! tmsortCustomFile) {
		/* Standard history file */
		
		if (tmsortDecayPending) {
			/* Copy over previous decay */
			memcpy(decayPtr, &tmsortNextDecayEvent, tmsortDecaySize);
			bufPtr = (LbUsOneByte *)decayPtr + tmsortDecaySize;
			*numPhotons = 0;
			
			/* Continue reading subsequent photon events */
			allDone = false;
			while (! allDone) {
				eventType = tmsortReadEvent(historyFile, 
											(PHG_Decay *)bufPtr, 
											(PHG_DetectedPhoton *)bufPtr);
				
				if (eventType == Photon) {
					/* A photon associated with the decay was read */
					bufPtr += tmsortPhotonSize;
					(*numPhotons)++;
				}
				else if (eventType == Decay) {
					/* New decay; copy it to the storage buffer */
					memcpy(&tmsortNextDecayEvent, bufPtr, tmsortDecaySize);
					allDone = true;
				}
				else {
					/* Some other sort of result, most likely end of file */
					allDone = true;
					tmsortDecayPending = false;
				}
			}
			
			*dataSize = bufPtr - (LbUsOneByte *)decayPtr;
		}
		else {
			*dataSize = 0;
		}
	}
	else {
		/* Custom history file */
		
		*dataSize = 0;
		*numPhotons = 0;
		photonPtr = (PHG_TrackingPhoton *)((LbUsOneByte *)decayPtr + tmsortDecaySize);
		while (feof(historyFile) == false) {
			
			/* Get the number of blue photons */
			if ((fread(&numBlue, sizeof(LbUsOneByte), 1, historyFile)) != 1) {
				if (feof(historyFile) == false) {
					/* Some kind of error, but just ignore it */
				}
				
				/* Stop reading data */
				break;
			}
			
			/* Loop through the blue photons */
			for (bIndex=0; bIndex<numBlue; bIndex++) {
				
				/* Read the photon information */
				if (PhoHFileReadFields(&tmsortHistParamsHk, decayPtr, 
										photonPtr, &photonAccepted) == false) {
					/* Problem */
					break;
				}
				
				/* Record the photon */
				LbFgSet(photonPtr->flags, PHGFg_PhotonBlue);
				photonPtr++;
				(*numPhotons)++;
			}
			
			/* If we are doing PET then read the pink photons */
			if (PHG_IsPET()) {
				/* Get the number of pink photons */
				if ((fread(&numPink, sizeof(LbUsOneByte), 1, historyFile)) != 1) {
					if (feof(historyFile) == false) {
						/* Some kind of error, but just ignore it */
					}
					
					/* Stop reading data */
					break;
				}
				
				/* Loop through the pink photons */
				for (pIndex=0; pIndex<numPink; pIndex++) {
					
					/* Read the photon information */
					if (PhoHFileReadFields(&tmsortHistParamsHk, decayPtr, 
											photonPtr, &photonAccepted) == false) {
						/* Problem */
						break;
					}
					
					/* Record the photon */
					LbFgClear(photonPtr->flags, PHGFg_PhotonBlue);
					photonPtr++;
					(*numPhotons)++;
				}
			}
		}
		
		*dataSize = (LbUsOneByte *)photonPtr - (LbUsOneByte *)decayPtr;
	}
	
	#ifdef TMSORT_PHO_DIST
		if (*dataSize != 0) {
			tmsortPhotonCounts[*numPhotons]++;
		}
	#endif
	
	return (*dataSize != 0);
}


/**********************
*
*	Name:		tmsortReadEvent
*
*	Purpose:	Read the next event from the input file.
*
*	Arguments:
*		FILE				*historyFile	- The history file.
*		PHG_Decay 			*decayPtr		- Storage for decay.
*		PHG_DetectedPhoton	*photonPtr		- Storage for photon.
*
*	Result:		The type of the event read.
*
***********************/
EventTy tmsortReadEvent(FILE *historyFile, PHG_Decay *decayPtr, PHG_DetectedPhoton * photonPtr)
{
	EventTy			eventType = Null;	/* The event we read */
	LbUsOneByte		flag;				/* Storage for type of event to read */
	LbUsFourByte	decaySize;			/* Size of a decay */
	LbUsFourByte	photonSize;			/* Size of a photon */
	
	
	do { /* Process Loop */

		/* See what type of event we have */
		if ((fread(&flag, sizeof(LbUsOneByte), 1, historyFile)) != 1) {
		
			/* See if we are not at end of file */
			if (feof(historyFile) == 0) {
				ErAbort("Unable to read event type.");
			}
			else {
				/* We are end of file so just break */
				break;
			}
		}
	
		/* See if we have a decay or a photon */
		if (PHG_IsADecay((LbUsFourByte) flag)) {
			
			/* Set our event type */
			eventType = Decay;
			
			/* Read the decay */
			decaySize = tmsortDecaySize;
			
			if ((fread(decayPtr, decaySize, 1, historyFile)) != 1) {
				ErAbort("Unable to read decay.");
			}
			
		}
		else if (PHG_IsAPhoton((LbUsFourByte) flag)){
			/* Set our event type */
			eventType = Photon;
			
			/* Read the photon */
			photonSize = tmsortPhotonSize;
			
			if ((fread(photonPtr, photonSize, 1, historyFile)) != 1) {
				ErAbort("Unable to read photon.");
			}
		}
		else {
			/* Bad flag event type -- this probably means the file is using an older 
				format with different length event types than the current ones, thus 
				causing the program to end up trying to read an event starting in the 
				middle of another event.  */
			ErAbort("This data file appears to be stored in an outdated format.");
		}
		
	} while (false);
	
	return (eventType);
}


/*********************************************************************************
*
*	Name:			tmsortWriteDecay
*
*	Summary:		Write the decay and photons to the file.
*
*	Arguments:
*		PhoHFileHkTy 		*hdrHkTyPtr		- The hook to the file.
*		PHG_Decay			decay			- The decay that started the process.
*		PHG_DetectedPhoton *decayPhotons	- The photons detected.
*		LbUsFourByte 		numPhotons		- The number of photons.
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean tmsortWriteDecay(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay,
			PHG_DetectedPhoton *decayPhotons, LbUsFourByte numPhotons)
{
	PHG_TrackingPhoton	*photonPtr;		/* Pointer to custom photons */
	PHG_TrackingPhoton	*bluePhotons;	/* Pointer to start of blue photons */
	LbUsFourByte		numBluePhotons;	/* Count of blue photons */
	PHG_TrackingPhoton	*pinkPhotons;	/* Pointer to start of pink photons */
	LbUsFourByte		numPinkPhotons;	/* Count of pink photons */
	LbUsFourByte		index;			/* Index through custom photons */
	
	
	if (! tmsortCustomFile) {
		/* Standard history file */
		return (PhoHFileRewriteDetections(hdrHkTyPtr, decay, 
					decayPhotons, numPhotons, NULL, 0));
	}
	else {
		/* Custom history file */
		
		photonPtr = (PHG_TrackingPhoton *)((LbUsOneByte *)decay + tmsortDecaySize);
		bluePhotons = photonPtr;
		numBluePhotons = 0;
		pinkPhotons = NULL;
		numPinkPhotons = 0;
		
		/* Count the blue and pink photons */
		for (index=0; index<numPhotons; index++) {
			if (LbFgIsSet(photonPtr->flags, PHGFg_PhotonBlue)) {
				/* Blue photon */
				numBluePhotons++;
			}
			else {
				/* Pink photon */
				numPinkPhotons++;
				if (pinkPhotons == NULL) {
					pinkPhotons = photonPtr;
				}
			}
			photonPtr++;
		}
		
		return (PhoHFileWriteDetections(hdrHkTyPtr, decay, 
					bluePhotons, numBluePhotons, pinkPhotons, numPinkPhotons));
	}
}


#undef TIME_SORT_MAIN
