/*********************************************************************************
 *
 *                       Source code developed by the
 *           Imaging Research Laboratory - University of Washington
*               (C) Copyright 1997-2013 Department of Radiology          	        
 *                           University of Washington
 *                              All Rights Reserved
 *
 *********************************************************************************/

/*********************************************************************************
 *
 *			Module Name:		readhist.c
*			Revision Number:	1.7
*			Date last revised:	10 June 2013
 *			Programmer:			Steven Vannoy
 *			Date Originated:	8 August 1997
 *
 *			Module Overview:	Reads a custom phg history file.
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
 **********************************************************************************
 *
 *			Revision Section (Also update version number, if relevant)
 *
 *			Programmer(s):      Robert Harrison
 *
 *			Revision date:      1 Jan 2013
 *
 *			Revision description:   Write out the detector crystal number to
 *                           the history file.
 *                                  Changes the photon direction so its norm
 *                          is acceptably (to other SimSET routines) close to
 *                          one (fixes a float -> double problem).
 *
**********************************************************************************
 *
 *			Revision Section (Also update version number, if relevant)
 *
 *			Programmer(s):		Robert Harrison
 *
 *			Revision date:		Mar 2010
 *			Revision description:
 *						- Changed initialization of trackingPhoton.sliceIndex
 *						to 0 to get rid of segmentation error.
 *
 **********************************************************************************
 *
 *			Revision Section (Also update version number, if relevant)
 *
 *			Programmer(s):		Robert Harrison
 *
 *			Revision date:		Jan-Feb 2005
 *			Revision description:
 *						- Added code for new custom fields: detector angle,
 *						decay position, and decay time.
 *						- Changed custom history file read function to call
 *						the processing option for each decay, rather than
 *						at the end of file processing.  Previously the
 *						processing lost all decay-specific info.
 *						- phgbinProcessPhotons now prints out whether blue
 *						or pink photons are currently being processed for PET.
 *
 **********************************************************************************
 *
 *			Revision Section (Also update version number, if relevant)
 *
 *			Programmer(s):		Robert Harrison
 *
 *			Revision date:		19 July 2004
 *			Revision description:
 *						- Fixed bugs in the processing of history files,
 *						including failure to zero out some counters,
 *						failure to account for scatter in the collimator.
 *						- Added logic so that old history files could
 *						still be processed as necessary.
 *
 *********************************************************************************/
#define PHG_BIN_MAIN	/* Note we are substituting ourselves for phg's main */

#include <stdio.h>
#include <string.h>

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
#include "Collimator.h"
#include "Detector.h"
#include "phg.h"
#include "PhgBin.h"


#ifdef MPW
#pragma segment PHG_BIN_MAIN
#endif

/* LOCAL CONSTANTS */
#define PHGRDHST_IsUsePHGHistory()		LbFgIsSet(PhgOptions, LBFlag0)		/* Will we use the PHG history file */
#define PHGRDHST_IsUseColHistory()		LbFgIsSet(PhgOptions, LBFlag1)		/* Will we use the Collimator history file */
#define PHGRDHST_IsUseDetHistory()		LbFgIsSet(PhgOptions, LBFlag2)		/* Will we use the Detector history file */

#define	PHGRDHST_NumFlags	3													/* Number of flags defined */

/* LOCAL TYPES */
typedef enum  {Null, Decay, Photon} EventTy;

/* LOCAL GLOBALS */
static Boolean				phgrdhstCanceled;				/* Global cancelation flag */
static CollimatedPhotonsTy	phgrdhstColPhotons;				/* These are the successfully collimated photons */
static DetectedPhotonsTy	phgrdhstDetPhotons;				/* These are the successfully detected photons */
static char					phgrdhstErrStr[1024];			/* Error string storage */
static LbUsFourByte			phgrdhstNumToProc;				/* Number of histories to process */
static char					phgrdhstHistName[1024];			/* Name of history file */
static char					phgrdhstHistParamsName[1024];	/* Name of history parameters file */
static LbUsFourByte			phgrdhstArgIndex;
static ProdTblProdTblInfoTy	phgrdhstPrdTblInfo;				/* Info for initializing productivity table */
static PhoHFileHdrTy		phgrdhstHdrParams;				/* Input header */

/* PROTOTYPES */
Boolean 		phgrdhstInitialize(int argc, char *argv[]);
void			phgrdhstTerminate(void);
Boolean			phgbin(int argc, char *argv[]);
EventTy			readEvent(FILE *historyFile,
						  PHG_Decay *decayPtr,
						  PHG_DetectedPhoton *photonPtr);
EventTy			oldReadEvent(FILE *historyFile,
                             PHG_Decay *decayPtr,
                             PHG_DetectedPhoton *photonPtr,
                             Boolean isOldPhotons1,
                             Boolean isOldPhotons2);
Boolean			phgrdhstStandard(char *argv[]);
Boolean			phgrdhstCustom(char *argv[]);
void			phgbinProcessPhotons(PhoHFileHkTy *histHk, PHG_Decay *decayPtr,
                                     PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBlues,
                                     PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinks);


/* FUNCTIONS */
/**********************
 *	phgrdhstInitialize
 *
 *	Purpose:	Perform initialization tasks.
 *
 *	Result:		True unless an error occurs.
 ***********************/
Boolean phgrdhstInitialize(int argc, char *argv[])
{
	Boolean				okay = false;				/* Process Loop */
    char			*knownOptions[] = {"pcd"};
	char				optArgs[PHGRDHST_NumFlags][LBEnMxArgLen];
	LbUsFourByte		optArgFlags = (LBFlag0);
	char				fileName[1024];				/* Name of param file */
	LbFourByte			randSeed;			/* Seed for random generator */
	
	do { /* Process Loop */
        
		/* Get our options */
		if (!LbEnGetOptions(argc, argv, knownOptions,
                            &PhgOptions, optArgs, optArgFlags, &phgrdhstArgIndex)) {
            
			break;
		}
		
		/* Make sure the didn't specify more than one history file */
		if (PHGRDHST_IsUsePHGHistory() && (PHGRDHST_IsUseColHistory() || PHGRDHST_IsUseDetHistory())) {
			ErStGeneric("You can only specify one type of history file, -p, -c, or -d.");
			break;
		}
		
		/* Make sure they specified a history file */
		if (!PHGRDHST_IsUsePHGHistory() && !PHGRDHST_IsUseColHistory() && !PHGRDHST_IsUseDetHistory()) {
			ErStGeneric("You must specify the use of a PHG history file (-p)\n"
                        " or a collimator history file (-c)\n"
                        " or a detector history file (-d)\n");
			break;
		}
        
		/* See if they supplied a command line argument */
		if ((phgrdhstArgIndex != 0) && (argv[phgrdhstArgIndex] != 0)) {
            
			/* Get first param file and save number of param files to process */
			strcpy(PhgRunTimeParams.PhgParamFilePath,argv[phgrdhstArgIndex]);
			phgrdhstNumToProc = (argc - phgrdhstArgIndex);
		}
		else {
			/* Ask for the file name */
			LbInAsk("Enter name of param file", 0, false,
					&phgrdhstCanceled, 0, 0, 0, 0,
					fileName);
            
			/* Bolt if we canceled */
			if (phgrdhstCanceled) {
				ErStCancel("User canceled.");
				goto CANCEL;
			}
			phgrdhstNumToProc = 1;
			strcpy(PhgRunTimeParams.PhgParamFilePath, fileName);
		}
		
		/* Get our run-time parameters */
		if (!PhgGetRunTimeParams())
			break;
        
		/* Clear the file name parameters */
		phgrdhstHistParamsName[0] = '\0';
		phgrdhstHistName[0] = '\0';
		
		/* Let them know what is going on */
		LbInPrintf("\n***********************\n");
		LbInPrintf("About to process parameter file '%s'\n", PhgRunTimeParams.PhgParamFilePath);
		LbInPrintf("\n***********************\n");
		
		
		/* If user requested to bin PHG history file, use the one specified in
         the param file.
         */
		if (PHGRDHST_IsUsePHGHistory() &&
            (strlen(PhgRunTimeParams.PhgPhoHFileHistoryFilePath) == 0)) {
			
			ErStGeneric("No history file supplied in run time parameters.");
			break;
		}
		else {
			strcpy(phgrdhstHistName, PhgRunTimeParams.PhgPhoHFileHistoryFilePath);
			strcpy(phgrdhstHistParamsName, PhgRunTimeParams.PhgPhoHParamsFilePath);
		}
		
		/* Initialize the math library */
		randSeed = PhgRunTimeParams.PhgRandomSeed;
		if (!PhgMathInit(&randSeed))
			break;
		/* Save the seed if it had come from the clock */
		if (PhgRunTimeParams.PhgRandomSeed == 0)
			PhgRunTimeParams.PhgRandomSeed = randSeed;
        
		/* Initialize the emission list manager */
		if (!EmisListInitialize())
			break;
		
		/* Initialize the sub-object manager */
		if (!SubObjInitialize()) {
			break;
  		}
        
		/* Initialize the productivity table manager */
		if (!ProdTblInitialize()) {
			break;
		}
		
		/* Create sub objects */
		if (!SubObjCreate()) {
			break;
		}
        
		/* Initialize the Cylinder Positions */
		{
			/* Set object cylinder */
			if (!CylPosInitObjectCylinder()) {
				break;
			}
			
			/* Set the criticial zone */
			if (!CylPosInitCriticalZone(PhgRunTimeParams.Phg_AcceptanceAngle)) {
				break;
			}
			
			/* Set the limit cylinder */
			CylPosInitLimitCylinder();
		}
		
		/* Setup up the productivity information */
		{
			/* Initialize productivity table info */
			if (PHG_IsStratification() == false) {
				phgrdhstPrdTblInfo.inputFileName = "";
				phgrdhstPrdTblInfo.outputFileName = "";
			}
			else {
				phgrdhstPrdTblInfo.inputFileName =
                PhgRunTimeParams.PhgProdTblInputTableFilePath;
				phgrdhstPrdTblInfo.outputFileName =
                PhgRunTimeParams.PhgProdTblOutputTableFilePath;
			}
            
			phgrdhstPrdTblInfo.acceptanceAngle =
            PhgRunTimeParams.Phg_AcceptanceAngle;
			
            
			/* Create productivity table */
			if (!SubObjGetStartingProdValues(&phgrdhstPrdTblInfo)) {
				break;
			}
            
		}
        
		/* Initialize the photon tracking module */
		if (!PhoTrkInitialize()) {
			break;
		}
        
		/* If we are going to use the PHG history file then we initialize
         the collimator module with "doColHistory" == true, the actual
         creation will depend on the collimator parameters.
         */
		/* Initialize the collimation module if necessary */
		if (PHG_IsCollimateOnTheFly()) {
			ColCurParams = 0;
			if (!ColInitialize(PHGRDHST_IsUsePHGHistory()))
				break;
		}
		
		/* Initialize the detection module if necessary */
		if (PHG_IsDetectOnTheFly()) {
			if (!DetInitialize(PHGRDHST_IsUsePHGHistory()))
				break;
		}
		
		/* make sure that slat collimation is only done with dual_headed detectors */
		if (PHG_IsCollimateOnTheFly()) {
			for (ColCurParams = 0; ColCurParams < ColNumParams; ColCurParams++) {
				DetCurParams = ColCurParams;
				if ( ( ColIsSlat() ) && ( !DetIsDualHeaded() ) ) {
					PhgAbort("Slat collimation is only available with dual-headed detectors (phgrdhstInitialize)", false);
				}
			}
			ColCurParams = 0;
		}
        
		/* Initialize the binning module if necessary */
		/* If a collimator history file was created in the simulation,
         we'll use it instead of the PHG history file.
         */
		if (PHGRDHST_IsUseColHistory() &&
            (strlen(ColRunTimeParams[ColCurParams].ColHistoryFilePath) == 0)){
            
			/* Let the user know they must supply the collimator history file */
			ErStGeneric("You specified using a collimator history file (-c)\n"
                        "but didn't supply one in the collimator parameter file");
			break;
		}
		else if (PHGRDHST_IsUseColHistory() == true){
			/* Save the current history file name */
			strcpy(phgrdhstHistName, ColRunTimeParams[ColCurParams].ColHistoryFilePath);
			strcpy(phgrdhstHistParamsName, ColRunTimeParams[ColCurParams].ColHistoryParamsFilePath);
		}
		if (PHGRDHST_IsUseDetHistory() &&
            (strlen(DetRunTimeParams[DetCurParams].DetHistoryFilePath) == 0)){
            
			/* Let the user know they must supply the detector history file */
			ErStGeneric("You specified using a detector history file (-d)\n"
                        "but didn't supply one in the detector parameter file");
			break;
		}
		else if (PHGRDHST_IsUseDetHistory() == true){
			/* Save the current history file name */
			strcpy(phgrdhstHistName, DetRunTimeParams[DetCurParams].DetHistoryFilePath);
			strcpy(phgrdhstHistParamsName, DetRunTimeParams[DetCurParams].DetHistoryParamsFilePath);
		}
        
		okay = true;
    CANCEL:;
	} while (false);
	
	return(okay);
    
}

/**********************
 *	phgrdhstTerminate
 *
 *	Purpose:	Terminate all of the managers.
 *
 *	Result:		None.
 ***********************/
void phgrdhstTerminate()
{
	
	/* Print the collimation report if initialized */
	if (!ErIsInError()){
		if (PHG_IsCollimateOnTheFly()&& !PHGRDHST_IsUseColHistory()) {
			ColPrintReport();
		}
        
		/* Print the detection report if initialized */
		if (PHG_IsDetectOnTheFly() && !PHGRDHST_IsUseDetHistory()) {
			DetPrintReport();
		}
        
		/* Print the binning report if initialized */
		if (PHG_IsBinOnTheFly()) {
			PhgBinPrintReport(&PhgBinParams[0], &PhgBinFields[0]);
		}
	}
	
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
	ProdTblTerminate(&phgrdhstPrdTblInfo);
	
	/* Terminate the sub object module */
	SubObjTerminate();
	
}

/**********************
 *	phgrdhstStandard
 *
 *	Purpose:	Process a standard history file.
 *
 *	Result:		True unless an error occurs.
 ***********************/
Boolean phgrdhstStandard(char *argv[])
{
	Boolean				okay = false;				/* Process flag */
	EventTy				eventType;					/* Type of current event */
	FILE				*historyFile;				/* The history file we are going to process */
	LbHdrHkTy			headerHk;					/* Hook to header */
	LbUsFourByte		numBluePhotons;				/* Number of blue photons for this decay */
	LbUsFourByte		numPinkPhotons;				/* Number of pink photons for this decay */
	LbUsFourByte		curFileIndex;				/* Current file index */
	LbUsFourByte		numPhotonsProcessed;		/* Number of photons processed */
	LbUsFourByte		numDecaysProcessed;			/* Number of decays processed */
	PHG_Decay		   	decay;						/* The decay */
	PHG_Decay		   	nextDecay;					/* The decay */
	PHG_DetectedPhoton	detectedPhoton;				/* The detected photon */
	PHG_TrackingPhoton	trackingPhoton;				/* The tracking photon */
	PHG_TrackingPhoton	*bluePhotons = 0;			/* Blue photons for current decay */
	PHG_TrackingPhoton	*pinkPhotons = 0;			/* Pink photon for current decay*/
	double 				angle_norm;					/* for normalizing photon direction */
	Boolean				isOldPhotons1;				/* is this a very old history file--must be
                                                     read using oldReadEvent */
	Boolean				isOldPhotons2;				/* is this a moderately old history file--must be
                                                     read using oldReadEvent */
	Boolean				isOldDecays;				/* is this an old history file--must be
                                                     read using oldReadEvent */
	Boolean				isPHGList;					/* this is PHG history file */
	Boolean				isColList;					/* this is collimator history file */
	Boolean				isDetList;					/* this is detector history file */
    
	do { /* Process Loop */
        
		/* Complete initialization for processing standard history file */
		{
            
			/* Setup for binning */
			if (!PhgBinInitialize(PhgRunTimeParams.PhgBinParamsFilePath[0], &PhgBinParams[0],
                                  &PhgBinData[0], &PhgBinFields[0]))
				break;
		}
		
		/* Attempt to allocate memory for photon arrays */
		if ((bluePhotons = (PHG_TrackingPhoton *) LbMmAlloc(sizeof(PHG_TrackingPhoton) *
                                                            PHG_MAX_DETECTED_PHOTONS)) == 0) {
			break;
		}
        
		if ((pinkPhotons = (PHG_TrackingPhoton *) LbMmAlloc(sizeof(PHG_TrackingPhoton) *
                                                            PHG_MAX_DETECTED_PHOTONS)) == 0) {
			break;
		}
        
		/* Clear loop variables */
		numPhotonsProcessed = 0;
		numDecaysProcessed = 0;
		
		/* Process the files */
		for (curFileIndex = 1; curFileIndex <= phgrdhstNumToProc; curFileIndex++){
            
			/* Open history file file */
			if ((historyFile = LbFlFileOpen(phgrdhstHistName,
                                            "rb")) == 0) {
				sprintf(phgrdhstErrStr, "Unable to open history file\n'%s'.",
                        phgrdhstHistName);
				ErStFileError(phgrdhstErrStr);
				goto FAIL;
			}
			
			/* Let them know what is going on */
			LbInPrintf("\n***********************\n");
			LbInPrintf("Read in header.\n");
			LbInPrintf("\n***********************\n");
            
			/* Read in the header and verify it is the right type of file */
			if (PhgHdrGtParams(historyFile, &phgrdhstHdrParams, &headerHk) == false){
				break;
			}
			
			/* Let them know what is going on */
			LbInPrintf("\n***********************\n");
			LbInPrintf("Check for input conflicts.\n");
			LbInPrintf("\n***********************\n");
            
			/* Verify old collimator/detector list mode files are not being used for SPECT/DHCI:
             photons had insufficient information for further processing--no detector angle was saved */
			{
				if ( (PHG_IsSPECT() || (DetRunTimeParams[DetCurParams].DetectorType == DetEn_DualHeaded))
                    && (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_COLOLD) ) {
					ErStGeneric("Old SPECT/DHCI collimator history file is inadequate for further processing.");
					goto FAIL;
				}
                
				if ( (PHG_IsSPECT() || (DetRunTimeParams[DetCurParams].DetectorType == DetEn_DualHeaded))
                    && (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_DETOLD) ) {
					ErStGeneric("Old SPECT/DHCI detector history file is inadequate for further processing.");
					goto FAIL;
				}
			}
			
			/* Set flags for type of list mode file being processed */
			if ( 	(phgrdhstHdrParams.H.HdrKind == PhoHFileEn_PHG) ||
                (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_PHG2625) ||
                (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_PHGOLD) ) {
				
				isPHGList = true;
				isColList = false;
				isDetList = false;
				
			} else if ( (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_COL) ||
                       (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_COL2625) ||
                       (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_COLOLD) ) {
				
				isPHGList = false;
				isColList = true;
				isDetList = false;
				
			} else if ( (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_DET) ||
                       (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_DET2625) ||
                       (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_DETOLD) ) {
				
				isPHGList = false;
				isColList = false;
				isDetList = true;
				
			} else {
				
				isPHGList = false;
				isColList = false;
				isDetList = false;
				
			}
			
			/* Verify that requested file is of the right type */
			{
				if ( PHGRDHST_IsUsePHGHistory() && (!isPHGList) ) {
					ErStGeneric("File specified as PHG history file is not valid.");
					goto FAIL;
				}
				
				if ( PHGRDHST_IsUseColHistory() && (!isColList) ) {
					ErStGeneric("File specified as Collimator history file is not valid.");
					goto FAIL;
				}
				
				if ( PHGRDHST_IsUseDetHistory() && (!isDetList) ) {
					ErStGeneric("File specified as Detector history file is not valid.");
					goto FAIL;
				}
			}
			
			/* Determine if this is an old-style history file--if so we will call oldReadEvent below */
			if ( 	(phgrdhstHdrParams.H.HdrKind == PhoHFileEn_PHGOLD) ||
                (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_COLOLD) ||
                (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_DETOLD) ) {
				
				isOldPhotons1 = true;
				isOldPhotons2 = false;
				isOldDecays = true;
				
			} else {
                
				isOldPhotons1 = false;
				if ( (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_PHG2625) ||
					(phgrdhstHdrParams.H.HdrKind == PhoHFileEn_COL2625) ||
					(phgrdhstHdrParams.H.HdrKind == PhoHFileEn_DET2625) ) {
					
					isOldPhotons2 = true;
					isOldDecays = true;
                    
				} else {
                    
					isOldPhotons2 = false;
					isOldDecays = false;
					
				}
				
			}
			
			/* Do the binning */
            
			/* Let them know what is going on */
			LbInPrintf("\n***********************\n");
			LbInPrintf("Begin decay processing.\n");
			LbInPrintf("\n***********************\n");
            
			{
				/* Read the initial event */
				if (isOldDecays) {
					eventType = oldReadEvent(historyFile, &decay, &detectedPhoton, isOldPhotons1, isOldPhotons2);
				} else {
					eventType = readEvent(historyFile, &decay, &detectedPhoton);
				}
                
				if (eventType != Decay) {
					ErStGeneric("Expected first event to be decay, and it wasn't.");
					goto FAIL;
				}
                
				/* Loop through all decays */
				while (eventType == Decay) {
                    
					numDecaysProcessed++;
					
					/* Clear decay variables */
					numBluePhotons = 0;
					numPinkPhotons = 0;
					
					/* Read next event */
					if (isOldDecays) {
						eventType = oldReadEvent(historyFile, &nextDecay, &detectedPhoton, isOldPhotons1, isOldPhotons2);
					} else {
						eventType = readEvent(historyFile, &nextDecay, &detectedPhoton);
					}
					
					while ( eventType == Photon) {
                        
						/* Convert to a tracking photon */
						{
							
							/* initialize fields not stored in standard history file */
							trackingPhoton.sliceIndex = 0;
							trackingPhoton.angleIndex =  -1;
							trackingPhoton.origSliceIndex = -1;
							trackingPhoton.origAngleIndex = -1;
							trackingPhoton.xIndex = -1;
							trackingPhoton.yIndex = -1;
							trackingPhoton.scatters_in_col = 0; /* We don't know this value--user must
                                                                 create a custom history file if they
                                                                 want it, otherwise these scatters are
                                                                 folded into the num_of_scatters field below  */
							trackingPhoton.scatter_target_weight  = 0;
							trackingPhoton.num_det_interactions = 0;
							
							/* copy decay weight from decay */
							trackingPhoton.decay_weight = decay.startWeight;
                            
							/* the rest of the tracking photon is copied from detectedPhoton */
							{
								
								/* Flags are currently stored in first two bytes, with upper bits
                                 representing number of scatters
                                 */
								trackingPhoton.flags = (detectedPhoton.flags & 3);
								trackingPhoton.num_of_scatters =
                                (detectedPhoton.flags >> 2);
                                
								trackingPhoton.location.x_position = detectedPhoton.location.x_position;
								trackingPhoton.location.y_position = detectedPhoton.location.y_position;
								trackingPhoton.location.z_position = detectedPhoton.location.z_position;
								trackingPhoton.angle.cosine_x = detectedPhoton.angle.cosine_x;
								trackingPhoton.angle.cosine_y = detectedPhoton.angle.cosine_y;
								trackingPhoton.angle.cosine_z = detectedPhoton.angle.cosine_z;
                                
								{
									/* normalize photon direction - it is written out as float but needs to be
                                     double precision unit length for som,e functions */
									angle_norm = sqrt(	(double)detectedPhoton.angle.cosine_x * (double)detectedPhoton.angle.cosine_x +
                                                      (double)detectedPhoton.angle.cosine_y * (double)detectedPhoton.angle.cosine_y +
                                                      (double)detectedPhoton.angle.cosine_z * (double)detectedPhoton.angle.cosine_z );
									trackingPhoton.angle.cosine_x = (double)detectedPhoton.angle.cosine_x / angle_norm;
									trackingPhoton.angle.cosine_y = (double)detectedPhoton.angle.cosine_y / angle_norm;
									trackingPhoton.angle.cosine_z = (double)detectedPhoton.angle.cosine_z / angle_norm;
								}
                                
								trackingPhoton.transaxialPosition = detectedPhoton.transaxialPosition;
								trackingPhoton.azimuthalAngleIndex = detectedPhoton.azimuthalAngleIndex;
								trackingPhoton.detectorAngle = detectedPhoton.detectorAngle;
								trackingPhoton.detCrystal = detectedPhoton.detCrystal;
								if (trackingPhoton.num_of_scatters == 0){
									trackingPhoton.photon_scatter_weight = 0;
									trackingPhoton.photon_primary_weight =
                                    detectedPhoton.photon_weight;
									trackingPhoton.photon_current_weight =
                                    detectedPhoton.photon_weight;
								}
								else {
									trackingPhoton.photon_scatter_weight =
                                    detectedPhoton.photon_weight;
									trackingPhoton.photon_current_weight =
                                    detectedPhoton.photon_weight;
									trackingPhoton.photon_primary_weight = 0;
								}
								trackingPhoton.energy = detectedPhoton.energy;
								trackingPhoton.travel_distance  =
                                detectedPhoton.time_since_creation * PHGMATH_SPEED_OF_LIGHT;
								trackingPhoton.numStarts = 0;
								trackingPhoton.number = (LbUsEightByte) -1;
                                
							}
						}
                        
						/* See if it is blue or pink */
						if (LbFgIsSet(trackingPhoton.flags, PHGFg_PhotonBlue)) {
							bluePhotons[numBluePhotons] = trackingPhoton;
							numBluePhotons++;
						}
						else {
							pinkPhotons[numPinkPhotons] = trackingPhoton;
							numPinkPhotons++;
						}
						
						numPhotonsProcessed++;
                        
						/* Read next event */
						if (isOldDecays) {
							eventType = oldReadEvent(historyFile, &nextDecay, &detectedPhoton, isOldPhotons1,isOldPhotons2);
						} else {
							eventType = readEvent(historyFile, &nextDecay, &detectedPhoton);
						}
                        
					}
					
					/* zero out processed event counters--this keeps us from mistakenly
                     reprocessing photons from the previous decay */
					phgrdhstColPhotons.NumCollimatedBluePhotons = 0;
					phgrdhstColPhotons.NumCollimatedPinkPhotons = 0;
					phgrdhstDetPhotons.NumDetectedBluePhotons = 0;
					phgrdhstDetPhotons.NumDetectedPinkPhotons = 0;
                    
					/* Write the info if there were any detections */
					if (PHG_IsPETCoincidencesOnly()){
						
						/* For PET check for both blue and pink photons to process detections */
						if ((numBluePhotons != 0) &&
                            (numPinkPhotons != 0)) {
							
							/* Collimate them if requested, note that if we are reading
                             a collimator history file then we don't want to do this,
                             even though PHG_IsCollimateOnTheFly is set to true.
							 */
							if ( PHG_IsCollimateOnTheFly() && (isPHGList) ) {
                                
								ColPETPhotons(&decay,
                                              bluePhotons, numBluePhotons,
                                              pinkPhotons, numPinkPhotons,
                                              &phgrdhstColPhotons);
							}
							
							/* Detect them if requested */
							if ( PHG_IsDetectOnTheFly() && (isPHGList || isColList) ) {
                                
								/* Send collimated photons if it was done */
								if ( PHG_IsCollimateOnTheFly() && (isPHGList) ) {
									DetPETPhotons(&decay,
                                                  phgrdhstColPhotons.CollimatedTrkngBluePhotons,
                                                  phgrdhstColPhotons.NumCollimatedBluePhotons,
                                                  phgrdhstColPhotons.CollimatedTrkngPinkPhotons,
                                                  phgrdhstColPhotons.NumCollimatedPinkPhotons,
                                                  &phgrdhstDetPhotons);
								}
								else {
									DetPETPhotons(&decay,
                                                  bluePhotons, numBluePhotons,
                                                  pinkPhotons, numPinkPhotons,
                                                  &phgrdhstDetPhotons);
								}
								
								PhgBinPETPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                 &decay,
                                                 phgrdhstDetPhotons.DetectedTrkngBluePhotons,
                                                 phgrdhstDetPhotons.NumDetectedBluePhotons,
                                                 phgrdhstDetPhotons.DetectedTrkngPinkPhotons,
                                                 phgrdhstDetPhotons.NumDetectedPinkPhotons);
							}
							else if ( PHG_IsCollimateOnTheFly() && (isPHGList) ) {
                                
                                PhgBinPETPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                 &decay,
                                                 phgrdhstColPhotons.CollimatedTrkngBluePhotons,
                                                 phgrdhstColPhotons.NumCollimatedBluePhotons,
                                                 phgrdhstColPhotons.CollimatedTrkngPinkPhotons,
                                                 phgrdhstColPhotons.NumCollimatedPinkPhotons);
							}
							else {
								/* Bin up non-collimated photons */
								PhgBinPETPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                 &decay,
                                                 bluePhotons, numBluePhotons,
                                                 pinkPhotons, numPinkPhotons);
							}
						}
					} else if (PHG_IsPETCoincPlusSingles()){
						
						/* For PET check there are blue or pink photons to process detections */
						if ((numBluePhotons != 0) ||
                            (numPinkPhotons != 0)) {
							
							/* Collimate them if requested, note that if we are reading
                             a collimator history file then we don't want to do this,
                             even though PHG_IsCollimateOnTheFly is set to true.
							 */
							if ( PHG_IsCollimateOnTheFly() && (isPHGList) ) {
                                
								ColPETPhotons(&decay,
                                              bluePhotons, numBluePhotons,
                                              pinkPhotons, numPinkPhotons,
                                              &phgrdhstColPhotons);
							}
							
							/* Detect them if requested */
							if ( PHG_IsDetectOnTheFly() && (isPHGList || isColList) ) {
                                
								/* Send collimated photons if it was done */
								if ( PHG_IsCollimateOnTheFly() && (isPHGList) ) {
									DetPETPhotons(&decay,
                                                  phgrdhstColPhotons.CollimatedTrkngBluePhotons,
                                                  phgrdhstColPhotons.NumCollimatedBluePhotons,
                                                  phgrdhstColPhotons.CollimatedTrkngPinkPhotons,
                                                  phgrdhstColPhotons.NumCollimatedPinkPhotons,
                                                  &phgrdhstDetPhotons);
								}
								else {
									DetPETPhotons(&decay,
                                                  bluePhotons, numBluePhotons,
                                                  pinkPhotons, numPinkPhotons,
                                                  &phgrdhstDetPhotons);
								}
								
								if ( PhgBinParams[0].isBinPETasSPECT ) {
									
									PhgBinSPECTPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                       &decay,
                                                       phgrdhstDetPhotons.DetectedTrkngBluePhotons,
                                                       phgrdhstDetPhotons.NumDetectedBluePhotons);
                                    
									PhgBinSPECTPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                       &decay,
                                                       phgrdhstDetPhotons.DetectedTrkngPinkPhotons,
                                                       phgrdhstDetPhotons.NumDetectedPinkPhotons);
                                    
									
								} else {
                                    
									PhgBinPETPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                     &decay,
                                                     phgrdhstDetPhotons.DetectedTrkngBluePhotons,
                                                     phgrdhstDetPhotons.NumDetectedBluePhotons,
                                                     phgrdhstDetPhotons.DetectedTrkngPinkPhotons,
                                                     phgrdhstDetPhotons.NumDetectedPinkPhotons);
									
								}
								
							}
							else if ( PHG_IsCollimateOnTheFly() && (isPHGList) ) {
                                
								if ( PhgBinParams[0].isBinPETasSPECT ) {
									
									PhgBinSPECTPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                       &decay,
                                                       phgrdhstColPhotons.CollimatedTrkngBluePhotons,
                                                       phgrdhstColPhotons.NumCollimatedBluePhotons);
                                    
									PhgBinSPECTPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                       &decay,
                                                       phgrdhstColPhotons.CollimatedTrkngPinkPhotons,
                                                       phgrdhstColPhotons.NumCollimatedPinkPhotons);
									
								} else {
                                    
									PhgBinPETPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                     &decay,
                                                     phgrdhstColPhotons.CollimatedTrkngBluePhotons,
                                                     phgrdhstColPhotons.NumCollimatedBluePhotons,
                                                     phgrdhstColPhotons.CollimatedTrkngPinkPhotons,
                                                     phgrdhstColPhotons.NumCollimatedPinkPhotons);
									
								}
								
							}
							else {
								/* Bin up non-collimated photons */
								if ( PhgBinParams[0].isBinPETasSPECT ) {
									
									PhgBinSPECTPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                       &decay,
                                                       bluePhotons, numBluePhotons);
                                    
									PhgBinSPECTPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                       &decay,
                                                       pinkPhotons, numPinkPhotons);
									
								} else {
									
									PhgBinPETPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                     &decay,
                                                     bluePhotons, numBluePhotons,
                                                     pinkPhotons, numPinkPhotons);
									
								}
								
							}
						}
					}
					else {
                        
						/* For SPECT check for blue photos to process detections */
						if (numBluePhotons != 0) {
                            
							/* Collimate them if necessary */
							if ( PHG_IsCollimateOnTheFly() && (isPHGList) ) {
								
								ColSPECTPhotons(&decay,
                                                bluePhotons, numBluePhotons,
                                                &phgrdhstColPhotons);
							}
                            
							/* Detect them if requested */
							if ( PHG_IsDetectOnTheFly() && ( isPHGList || isColList ) ) {
                                
								if ( PHG_IsCollimateOnTheFly() && (isPHGList) ) {
									
									DetSPECTPhotons(&decay,
                                                    phgrdhstColPhotons.CollimatedTrkngBluePhotons,
                                                    phgrdhstColPhotons.NumCollimatedBluePhotons,
                                                    &phgrdhstDetPhotons);
								}
								else {
									DetSPECTPhotons(&decay,
                                                    bluePhotons, numBluePhotons,
                                                    &phgrdhstDetPhotons);
								}
                                
								PhgBinSPECTPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                   &decay,
                                                   phgrdhstDetPhotons.DetectedTrkngBluePhotons,
                                                   phgrdhstDetPhotons.NumDetectedBluePhotons);
							}
							else if ( PHG_IsCollimateOnTheFly() && (isPHGList) ) {
                                
								PhgBinSPECTPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                   &decay,
                                                   phgrdhstColPhotons.CollimatedTrkngBluePhotons,
                                                   phgrdhstColPhotons.NumCollimatedBluePhotons);
							}
							else {
								/* phgrdhst up non-collimated photons */
								PhgBinSPECTPhotons(&PhgBinParams[0], &PhgBinData[0], &PhgBinFields[0],
                                                   &decay,
                                                   bluePhotons, numBluePhotons);
							}
						}
					}
					
					/* Update current decay */
					decay = nextDecay;
				}
			}
			/* Open next file */
			if (curFileIndex < phgrdhstNumToProc) {
				
				strcpy(PhgRunTimeParams.PhgParamFilePath,argv[curFileIndex+phgrdhstArgIndex]);
				
				/* Get our run-time parameters */
				if (!PhgGetRunTimeParams()) {
					goto FAIL;
				}
				
				/* Let them know what is going on */
				LbInPrintf("\n***********************\n");
				LbInPrintf("About to process parameter file '%s'\n", PhgRunTimeParams.PhgParamFilePath);
				LbInPrintf("\n***********************\n");
                
				/* Get the collimator parameters */
				if (PHG_IsCollimateOnTheFly() == true) {
					if (ColGetRunTimeParams() == false) {
						goto FAIL;
					}
				}
				
				/* Get the detector parameters */
				if (PHG_IsDetectOnTheFly() == true) {
					if (DetGetRunTimeParams() == false) {
						goto FAIL;
					}
				}
				
				/* Set the next history file to the next collimator history */
				if (PHGRDHST_IsUseColHistory() &&
                    (strlen(ColRunTimeParams[ColCurParams].ColHistoryFilePath) == 0)){
                    
					/* Let the user know they must supply the collimator history file */
					ErStGeneric("You specified using a collimatory history file (-c)\n"
                                "but didn't supply one in the collimator parameter file");
					break;
				}
				else {
					/* Save the current history file name */
					strcpy(phgrdhstHistName, ColRunTimeParams[ColCurParams].ColHistoryFilePath);
				}
				
				/* Or, set the next history file to the PHG history file */
				if (PHGRDHST_IsUsePHGHistory() &&
                    (strlen(PhgRunTimeParams.PhgPhoHFileHistoryFilePath) == 0)) {
                    
                    ErStGeneric("No history file supplied in run time parameters.");
                    break;
				}
				else {
					strcpy(phgrdhstHistName, PhgRunTimeParams.PhgPhoHFileHistoryFilePath);
				}
			}
		}
        
		okay = true;
    FAIL:;
    CANCEL:;
	} while (false);
	
	/* Free memory */
	if (bluePhotons != 0)
		LbMmFree((void **) &bluePhotons);
    
	if (pinkPhotons != 0)
		LbMmFree((void **) &pinkPhotons);
    
	return(okay);
}

/**********************
 *	phgrdhstCustom
 *
 *	Purpose:	Process a custom history file.
 *
 *	Result:		True unless an error occurs.
 ***********************/
Boolean phgrdhstCustom(char *argv[])
{
	char**				dummyPtr;			/* Remove compiler warning */
	Boolean 			okay = false; 		/* Process flag */
	Boolean				photonAccepted;		/* Was the photon accepted */
	FILE				*historyFile=0;		/* The file we are processing */
	PHG_Decay			decay;				/* The decay we are reading */
	PHG_TrackingPhoton	*bluePhotons;		/* The blue photons we read */
	PHG_TrackingPhoton	*pinkPhotons;		/* The pink photon we read */
	PHG_TrackingPhoton	bPhoton;			/* The blue photon we are reading */
	PHG_TrackingPhoton	pPhoton;			/* The pink photon we are reading */
	LbUsOneByte			numBlue;			/* Number of blue photons for this decay */
	LbUsOneByte			numPink;			/* Number of pink photons for this decay */
	LbUsOneByte			acceptedBlues;		/* Number of blue photons for this decay */
	LbUsOneByte			acceptedPinks;		/* Number of pink photons for this decay */
	LbUsOneByte			bIndex;				/* Current blue photon for this decay */
	LbUsOneByte			pIndex;				/* Current pink photon for this decay */
	LbHdrHkTy			headerHk;			/* Hook to header */
	PhoHFileHkTy		histHk;				/* Hook to history file information */
	Boolean				isPHGList;			/* this is PHG history file */
	Boolean				isColList;			/* this is collimator history file */
	Boolean				isDetList;			/* this is detector history file */
	
	dummyPtr = argv;		/* Avoids unused parameter compiler warning */
	
	do { /* Process Loop */
        
		/* Open history file file */
		if ((historyFile = LbFlFileOpen(phgrdhstHistName,
                                        "rb")) == 0) {
			sprintf(phgrdhstErrStr, "Unable to open history file\n'%s'.",
                    phgrdhstHistName);
			ErStFileError(phgrdhstErrStr);
			goto FAIL;
		}
		
		/* Read in the header and verify it is the right type of file */
		if (PhgHdrGtParams(historyFile, &phgrdhstHdrParams, &headerHk) == false){
			break;
		}
		
        /* Set flags for type of list mode file being processed */
        if ( 	(phgrdhstHdrParams.H.HdrKind == PhoHFileEn_PHG) ||
            (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_PHG2625) ||
            (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_PHGOLD) ) {
            
            isPHGList = true;
            isColList = false;
            isDetList = false;
            
        } else if ( (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_COL) ||
                   (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_COL2625) ||
                   (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_COLOLD) ) {
            
            isPHGList = false;
            isColList = true;
            isDetList = false;
            
        } else if ( (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_DET) ||
                   (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_DET2625) ||
                   (phgrdhstHdrParams.H.HdrKind == PhoHFileEn_DETOLD) ) {
            
            isPHGList = false;
            isColList = false;
            isDetList = true;
            
        } else {
            
            isPHGList = false;
            isColList = false;
            isDetList = false;
            
        }
        
        /* Verify that requested file is of the right type */
        {
            if ( PHGRDHST_IsUsePHGHistory() && (!isPHGList) ) {
                ErStGeneric("File specified as PHG history file is not valid.");
                goto FAIL;
            }
            
            if ( PHGRDHST_IsUseColHistory() && (!isColList) ) {
                ErStGeneric("File specified as Collimator history file is not valid.");
                goto FAIL;
            }
            
            if ( PHGRDHST_IsUseDetHistory() && (!isDetList) ) {
                ErStGeneric("File specified as Detector history file is not valid.");
                goto FAIL;
            }
        }
        
		
		/* Setup header hook */
		{
			histHk.doCustom = true;
			
			/* Do custom parameter initialization  */
			if (PhoHFileGetRunTimeParams(phgrdhstHistParamsName, &(histHk.customParams)) == false) {
				sprintf(phgrdhstErrStr,"Unable to get custom parameters for history file named '%s'",
                        PhgRunTimeParams.PhgPhoHParamsFilePath);
				ErAlert(phgrdhstErrStr, false);
				goto FAIL;
			}
			
			strcpy(histHk.customParamsName, phgrdhstHistParamsName);
            
			histHk.bluesReceived = 0;
			histHk.bluesAccepted = 0;
			histHk.pinksReceived = 0;
			histHk.pinksAccepted = 0;
			histHk.histFile = historyFile;
		}
		
		/* Attempt to allocate memory for photon arrays */
		if ((bluePhotons = (PHG_TrackingPhoton *) LbMmAlloc(sizeof(PHG_TrackingPhoton) *
                                                            PHG_MAX_DETECTED_PHOTONS)) == 0) {
			break;
		}
        
		if ((pinkPhotons = (PHG_TrackingPhoton *) LbMmAlloc(sizeof(PHG_TrackingPhoton) *
                                                            PHG_MAX_DETECTED_PHOTONS)) == 0) {
			break;
		}
        
		/* Loop through reading event information */
		while (feof(historyFile) == false) {
            
			/* Clear counters */
			acceptedBlues = 0;
			acceptedPinks = 0;
			
			/* Get the number of blue photons */
			if ((fread(&numBlue, sizeof(LbUsOneByte), 1, historyFile)) != 1) {
				if (feof(historyFile) == false)
					goto FAIL;
				else
					break;
                
			}
			
			/* Loop through blue photons */
			for (bIndex = 0; bIndex < numBlue; bIndex++) {
                
				/* Read the photon information */
				if (PhoHFileReadFields(&histHk, &decay, &bPhoton, &photonAccepted) == false) {
					goto FAIL;
				}
				
				/* Go to next photon if not accepted by filtering */
				if (photonAccepted == false)
					continue;
                
				/* Save the photon */
				bluePhotons[acceptedBlues] = bPhoton;
				
				acceptedBlues++;
			}
            
			/* If we are doing PET then read pink */
			if (PHG_IsPET()) {
				/* Get the number of pink photons */
				if ((fread(&numPink, sizeof(LbUsOneByte), 1, historyFile)) != 1) {
					goto FAIL;
				}
				
				/* Loop through pink photons */
				for (pIndex = 0; pIndex < numPink; pIndex++) {
                    
					/* Read the photon information */
					if (PhoHFileReadFields(&histHk, &decay, &pPhoton, &photonAccepted) == false) {
						goto FAIL;
					}
					
					/* Go to next photon if not accepted by filtering */
					if (photonAccepted == false)
						continue;
                    
					
					/* Save the photon */
					pinkPhotons[acceptedPinks] = pPhoton;
					
					acceptedPinks++;
				}
			}
			
			/*  Process the photons */
			phgbinProcessPhotons(&histHk, &decay, bluePhotons, acceptedBlues,
                                 pinkPhotons, acceptedPinks);
			
			
		}
		
		
        okay = true;
	FAIL:;
	} while (false);
	
	if (historyFile != 0) {
		fclose(historyFile);
	}
	return(okay);
}

/*********************************************************************************
 *
 *			Name:			phgbinProcessPhotons
 *
 *			Summary:		Process the photon information read.
 *
 *			Arguments:
 *				PhoHFileHkTy			*histHk					- Hook to history file information.
 *				PHG_Decay				*decayPtr				- The decay that started the process.
 *				PHG_TrackingPhoton 		*bluePhotons			- The photon detected.
 *				LbUsFourByte	 		numBlues				- The # of photon detected.
 *				PHG_TrackingPhoton 		*pinkPhotons			- The photon detected.
 *				LbUsFourByte	 		numPinks				- The # of photon detected.
 *
 *			Function return: None.
 *
 *********************************************************************************/
void phgbinProcessPhotons(PhoHFileHkTy *histHk, PHG_Decay *decayPtr,
                          PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBlues,
                          PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinks)
{
    
	LbUsFourByte bIndex, pIndex;
    
	/* If we are doing PET, specify that we are printing the blue photons */
	if (PHG_IsPET()) {
		LbInPrintf("\nBlue photons\n");
	}
	
	/* Process the blues */
	for (bIndex = 0; bIndex < numBlues; bIndex++) {
		PhoHFilePrintFields(histHk, decayPtr, &bluePhotons[bIndex]);
	}
    
	/* Process the pinks */
	if (PHG_IsPET()) {
		LbInPrintf("\nPink photons\n");
		for (pIndex = 0; pIndex < numPinks; pIndex++) {
			PhoHFilePrintFields(histHk, decayPtr, &pinkPhotons[pIndex]);
		}
	}
	
}

/**********************
 *	phgbin
 *
 *	Purpose:	Execute the program.
 *
 *	Result:		True unless an error occurs.
 ***********************/
Boolean phgbin(int argc, char *argv[])
{
	Boolean				okay = false;				/* Process Loop */
	
    
	do { /* Process Loop */
        
		/* Perform initialization tasks */
		if (phgrdhstInitialize(argc, argv) == false) {
			break;
		}
		
		/* See if this is a standard or custom history file */
		if (strlen(phgrdhstHistParamsName) == 0) {
			if ( phgrdhstStandard(argv) == false ) {
				goto FAIL;
			}
		}
		else {
			if ( phgrdhstCustom(argv) == false ) {
				goto FAIL;
			}
		}
		okay = true;
    FAIL:;
    CANCEL:;
	} while (false);
	
	/* Terminate PHG modules */
	phgrdhstTerminate();
	
	/* Handle error situation if one exists */
	if (!okay && phgrdhstCanceled) {
		ErHandle("User canceled readhist", false);
		okay = true;
	}
    
	/* Quit the program */
	return (okay);
}

/**********************
 *	readEvent
 *
 *	Purpose:	Read the next event from the input file.
 *
 *	Arguments:
 *		FILE				*historyFile	- The history file.
 *		PHG_Decay 			*decayPtr		- Storage for decay.
 *		PHG_DetectedPhoton	*photonPtr		- Storage for photon.
 *
 *	Result:	Enumerated type corresponding to event type.
 ***********************/
EventTy readEvent(
                  FILE *historyFile,
                  PHG_Decay *decayPtr,
                  PHG_DetectedPhoton * photonPtr	)
{
	EventTy				retEventType = Null;			/* The returned event type we read */
	PhoHFileEventType	eventType = PhoHFileNullEvent;	/* The event type we read */
	
	/* Call PhoHFileReadEvent */
	eventType = PhoHFileReadEvent( historyFile, decayPtr, photonPtr );
	
	/* Convert to the local event type */
	switch ( eventType ) {
		case PhoHFileNullEvent:
			retEventType = Null;
			break;
            
		case PhoHFileDecayEvent:
			retEventType = Decay;
			break;
            
		case PhoHFilePhotonEvent:
			retEventType = Photon;
			break;
            
		default:
			retEventType = Null;
			break;
	}
	
	return (retEventType);
}

/**********************
 *	oldReadEvent
 *
 *	Purpose:	Read the next event from an old input file.
 *				Old list mode files did not have all the
 *				info needed for processing SPECT/DHCI collimator
 *				and detector list mode data.  This routine
 *				is used for reading other old list mode files
 *				(e.g. PHG list mode files) which can be
 *				processed.
 *	Arguments:
 *		FILE				*historyFile		- The history file.
 *		PHG_Decay 			*decayPtr		- Storage for decay.
 *		PHG_DetectedPhoton	*photonPtr		- Storage for photon.
 *
 *	Result:	The index of the first array element >= phi
 ***********************/
EventTy oldReadEvent(
                     FILE *historyFile, 
                     PHG_Decay *decayPtr, 
                     PHG_DetectedPhoton *photonPtr,
                     Boolean isOldPhotons1,
                     Boolean isOldPhotons2 )
{
	EventTy				retEventType = Null;			/* The returned event type we read */
	PhoHFileEventType	eventType = PhoHFileNullEvent;	/* The event type we read */
	PHG_OldDecay		oldDecay;
	
	/* Call PhoHFileReadEvent */
	eventType = PhoHFileOldReadEvent( historyFile, &oldDecay, photonPtr, isOldPhotons1, isOldPhotons2 );
    
	/* Convert to the local event type */
	switch ( eventType ) {
		case PhoHFileNullEvent:
			retEventType = Null;
			break;
            
		case PhoHFileDecayEvent:
			(*decayPtr).location = oldDecay.location;
			(*decayPtr).startWeight = oldDecay.startWeight;
			
			/* set decay time to random time within current time frame
             in picoseconds from beginning */
			(*decayPtr).decayTime = PhgMathGetRandomNumber() * 
            (double)(PHGGetLengthOfScan());
			/* only two types of decay were available in old decays */
			if (PHG_IsPET()) {
				(*decayPtr).decayType = PhgEn_Positron;
			} else {
				(*decayPtr).decayType = PhgEn_SinglePhoton;
			}
			
			retEventType = Decay;
			break;
            
		case PhoHFilePhotonEvent:
			retEventType = Photon;
			break;
            
		default:
			retEventType = Null;
			break;
	}
	
	return (retEventType);
}
#undef PHG_BIN_MAIN
