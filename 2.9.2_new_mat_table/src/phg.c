/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992-2013 Department of Radiology	     		*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			phg.c
*     Revision Number:		1.7
*     Date last revised:	23 July 2013
*     Programmer:			Steven Vannoy
*     Date Originated:		9 September 1992
*
*     Module Overview:	This is the main module for the phg.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*     Global variables defined:
*
*	  Global macros defined:
*		
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		
*
*			Revision date:		
*
*			Revision description:
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*			Revision description:
*						- support for randoms and eight-byte number of decays
*
*********************************************************************************/

#define PHG_MAIN

#include <stdio.h>
#include <signal.h>
#include <time.h>
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
#include "PhgMath.h"
#include "CylPos.h"
#include "ProdTbl.h"
#include "SubObj.h"
#include "EmisList.h"
#include "PhoHStat.h"
#include "PhoHFile.h"
#include "ColUsr.h"
#include "UNCCollimator.h"
#include "Collimator.h"
#include "ColSlat.h"
#include "DetCylinder.h"
#include "Detector.h"
#include "DetUsr.h"
#include "DetPlanar.h"
#include "PhoTrk.h"
#include "phg.h"
#include "PhgBin.h"

/* LOCAL CONSTANTS */
#define PHG_PARAM_FILENAME		"phg.run"			/* Name of parameter file */

/* LOCAL TYPES */

/* LOCAL GLOBALS */
ProdTblProdTblInfoTy	prodTableInfo;				/* Info for initializing productivity table */
static char				phgErrString[1024];			/* General error string */

/* PROTOTYPES */
Boolean 				phg_initialize(void);
void					phg_terminate(void);
Boolean					phgValidateParams(void);	

/* FUNCTIONS */

/*********************************************************************************
*
*			Name:		PhgRun
*
*			Summary:	Begin the PHG run.
*			Arguments:
*
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean PhgRun(int argc, char *argv[])	
{
	Boolean	okay = false;			/* Process Flag */
	Boolean	randFromClock = false;	/* Did we take random seed from the clock */
	#ifdef AOS_VS
		char		*knownOptions[] = {"Debug"
										(char *) 0};
	#else
		char			*knownOptions[] = {"d:"};
	#endif
	char				optArgs[PHG_NumFlags][LBEnMxArgLen];
	LbUsFourByte		optArgFlags = (LBFlag0);
	LbUsFourByte		argIndex;
	time_t 				curTime;							/* Current time for stamping execution date */
	
	do {
	
		/* Get our options */
		if (!LbEnGetOptions(argc, argv, knownOptions,
				&PhgOptions, optArgs, optArgFlags, &argIndex)) {
	
			break;
		}

		/* If they gave us debug options, get them, otherwise set to default 'nil' */
		if (PHG_IsDebugOptions()) {
			
			PhgDebugOptions = atoi(optArgs[0]);
		}
		else {
			PhgDebugOptions = 0;
		}

		/* If they gave us a run file, then get it */
		if ((argIndex != 0) && (argv[argIndex] != 0)) {
			strcpy(PhgRunTimeParams.PhgParamFilePath, argv[argIndex]);

			/* Go to next argument */
			argIndex++;
		}
		else
			strcpy(PhgRunTimeParams.PhgParamFilePath, PHG_PARAM_FILENAME);
		
		#ifdef PHG_DEBUG_DO_DUMP_FILE
			/* Open the debug dump file */
			if (argIndex < argc) {
				if ((PhgDebugDumpFile = LbFlFileOpen(argv[argIndex], "w")) == 0){
					LbInPrintf("\nUnable to open dump file %s.\n", argv[argIndex]);
					break;
				}
				argIndex++;
			}
			else if ((PhgDebugDumpFile = LbFlFileOpen("phg_debug_dump", "w")) == 0){
				LbInPrintf("\nUnable to open dump file.\n");
				break;
			}

		#endif
		
		/* Set our scan length to 60 seconds, this is the default */
		PhgRunTimeParams.Phg_LengthOfScan = 60;

		/* Get our parameters */
		if (PhgGetRunTimeParams() == false)
			break;
		
		/* If the random seed is zero now, then we know it will be set by
			the system clock, we want to know this for the user's sake.
		*/
		randFromClock = (PhgRunTimeParams.PhgRandomSeed <= 0);

		/* Do consistency checks on parameters */
		if (phgValidateParams() == false) {
			break;
		}
		
		/* Initialize the phg */
		if (phg_initialize() == false) {
			break;
		}
		
		/* Print out operation parameters */
		PhgPrintParams(argc, argv, randFromClock);
		
		/* Generate the list */
		if (EmisListCreatePhotonList() == false) {
			break;
		}
	
		okay = true;
	} while (false);

	if (okay) {
		/* Get current time */
		time(&curTime);
		
		/* Convert time to string format */
		#ifdef __MWERKS__
			strftime(PhgExecutionDateStr, 31, "%m %d %Y (%I:%M:%S %p)", localtime(&curTime));
		#elif defined WINNT
			strftime(PhgExecutionDateStr, 31, "%x (%X)", localtime(&curTime));	
		#else 
			strftime(PhgExecutionDateStr, 31, "%m %d %Y (%r)", localtime(&curTime));
		#endif

		/* Print out our final message */
		LbInPrintf("\nExecution of PHG finished on %s\n", PhgExecutionDateStr);
	}
	
	
	/* Terminate the process */
	phg_terminate();
		
	if (okay == false) {
		ErHandle("Unable to perform simulation", false);
		okay = true;
	}
				
	#ifdef PHG_DEBUG
		/* Close the debug dump file */
		if (PhgDebugDumpFile != 0){
			fclose(PhgDebugDumpFile);
			PhgDebugDumpFile = 0;
		}
	#endif

	
	return(okay);
}

/*********************************************************************************
*
*			Name:		PhgPrintParams
*
*			Summary:	Print the simulation parameters.
*			Arguments:
*				int		argc			- The number of arguments to the program
*				char  	*argv[]			- The program arguments
*				Boolean	randFromClock	- Did the random seed come from the clock
*			Function return: none.
*
*********************************************************************************/
void PhgPrintParams(int argc, char *argv[], Boolean randFromClock)	
{
	LbFourByte		argIndex;
	LbUsFourByte	curBinParams;	/* LCV for multiple binning parameters */
	time_t			curTime;		/* Current time for stamping execution date */
	LbUsFourByte	curTomo;
	
	/* Get current time */
	time(&curTime);
	
	/* Print out the command line */
	LbInPrintf("\n\n\nCommand line: ");
	for (argIndex = 0; argIndex < argc; argIndex++)
		LbInPrintf("%s ", argv[argIndex]);
	LbInPrintf("\n");
	
	/* Convert time to string format */
	#ifdef __MWERKS__
		strftime(PhgExecutionDateStr, 31, "%m %d %Y (%I:%M:%S %p)", localtime(&curTime));
	#elif defined WINNT
		strftime(PhgExecutionDateStr, 31, "%x (%X)", localtime(&curTime));
	#else 
		strftime(PhgExecutionDateStr, 31, "%m %d %Y (%r)", localtime(&curTime));
	#endif
	
	
	LbInPrintf("\nExecution of phg occurred on %s", PhgExecutionDateStr);
	
	if ( PHG_IsSPECT() ) {
		LbInPrintf("\nSimulating SPECT.");
	} else if ( PHG_IsPETCoincidencesOnly() ) {
		LbInPrintf("\nSimulating PET (coincidences only).");
	} else if ( PHG_IsPETCoincPlusSingles() ) {
		LbInPrintf("\nSimulating PET (coincidences plus singles).");
	}
		
	LbInPrintf("\nStratification is %s.", (PHG_IsStratification() ? "on" : "off"));
	LbInPrintf("%sForced Detection is %s.", 
		(PHOTRKIsConeBeamForcedDetection() ? "\nCone Beam " : "\n"), (PHG_IsForcedDetection() ? "on" : "off"));
	LbInPrintf("\nForced non-absorption is %s.", (PHG_IsNoForcedNonAbsorbtion() ? "off" : "on"));
	LbInPrintf("\nModelling positron range is %s.", PHG_IsRangeAdjust() ? "on" : "off");
	if (PHG_IsRangeAdjust() == true) {
		LbInPrintf("\n\tIsotope being modelled is %s", phgEn_IsotopeStr[PhgRunTimeParams.PhgNuclide.isotope]);
	}
	LbInPrintf("\nModelling non-collinearity is %s.", PHG_IsNonCollinearityAdjust() ? "on" : "off");
	LbInPrintf("\nModelling coherent scatter in object is %s.", (PHG_IsModelCoherentInObj() ? "on" : "off"));
	LbInPrintf("\nModelling coherent scatter in tomograph is %s.", (PHG_IsModelCoherentInTomo() ? "on" : "off"));
	LbInPrintf("\nModelling polarization is %s.", (PHG_IsModelPolarization() ? "on" : "off"));
	LbInPrintf("\nPhoton energy is            %3.1f keV.",
		PhgRunTimeParams.PhgNuclide.photonEnergy_KEV);
	LbInPrintf("\nMinimum energy threshold is %3.1f", PhgRunTimeParams.PhgMinimumEnergy);
	LbInPrintf("\nAcceptance angle is (+/-) %3.1f degrees.",
		PhgRunTimeParams.Phg_AcceptanceAngle);
	LbInPrintf("\nSine of acceptance angle is %3.3f.", PHGGetSineOfAccAngle());
	
	/* user can specify Phg_EventsToSimulate, or leave it to be calculated by SimSET */
	if ( !PhgRunTimeParams.PhgIsCalcEventsToSimulate ) {
		LbInPrintf("\nNumber of decays to simulate (user-supplied) = %lld.", PhgRunTimeParams.Phg_EventsToSimulate);
	} else {
		LbInPrintf("\nDecays to simulate (calculated from activity and scan time) = %lld.", PhgRunTimeParams.Phg_EventsToSimulate);
	}
	
	LbInPrintf("\nScan time = %9.3f seconds.", PHGGetLengthOfScan());
	
	if (PhgNumTomoFiles > 0) {
		for (curTomo = 0; curTomo < PhgNumTomoFiles; curTomo++) {
			LbInPrintf("\nTomograph files being used are: '%s'",
				PhgRunTimeParams.PhgTomoParamsFilePath[curTomo]);
		}
	}
	LbInPrintf("\nCollimator modelling is %s.", PHG_IsCollimateOnTheFly() ? "on" : "off");
	LbInPrintf("\nDetector modelling is %s.", PHG_IsDetectOnTheFly() ? "on" : "off");
	LbInPrintf("\nBinning is %s.", PHG_IsBinOnTheFly() ? "on" : "off");
	LbInPrintf("\nCreation of target-cylinder history file is %s.", PHG_IsNoHist() ? "off" : "on");
	if (PHG_IsNoHist() == false) {
		LbInPrintf("\nCustomized target-cylinder history file is %s.", PHG_IsHistParams() ? "on" : "off");
	}
	#ifdef PHG_DEBUG
	
		LbInPrintf("\nDebug options = %d.", PhgDebugOptions);
		
		if (PHGDEBUG_FixedDirection()) {
			LbInPrintf("\n\tPhotons started in fixed direction (%7.4f, %7.4f, %7.4f).",
						SUBOBJ_FIXED_DIR_COSINE_X, SUBOBJ_FIXED_DIR_COSINE_Y, SUBOBJ_FIXED_DIR_COSINE_Z);
		}
		if (PHGDEBUG_BinFirstDetPosition()) {
			LbInPrintf("\n\tFirst interaction within detector is being binned (will print to screen).");
		}
		if (PHGDEBUG_ReadFromRandSeedFile()) {
			LbInPrintf("\n\tRandom number seed being read from file.");
		}
		if (PHGDEBUG_WriteToRandSeedFile()) {
			LbInPrintf("\n\tRandom number seed being written to file.");
		}
		if (PHGDEBUG_ResetRandSeedFile()) {
			LbInPrintf("\n\tRandom seed file is being reset.");
		}
		if (PHGDEBUG_BinCentroidPosition()) {
			LbInPrintf("\n\tDetector is binning centroid positions (files named detector.*).");
		}
		if (PHGDEBUG_DetAbsorbOnly()) {
			LbInPrintf("\n\tDetector is forcing all interactions to be absorptions.");
		}
		if (PHGDEBUG_ObjCohOnly()) {
			LbInPrintf("\n\tAll interactions in object are being modeled as coherent.");
		}
		
	#else
		LbInPrintf("\nDebugging is OFF.");
	#endif
	
	LbInPrintf("\nVoxel as point source is %s.", (PHG_IsVoxelPointSource() ? "on" : "off"));
	LbInPrintf("\nVoxel as line source is %s.", (PHG_IsVoxelLineSource() ? "on" : "off"));
	LbInPrintf("\nRandom seed is %ld  (%s).",
			(unsigned long)PhgRunTimeParams.PhgRandomSeed, ((randFromClock == true) ? "from system clock" : "user specified"));
	LbInPrintf("\n");
	
	/* Print out collimator parameters if appropriate */
	if (PHG_IsCollimateOnTheFly())
		ColPrintParams();

	/* Print out detector parameters if appropriate */
	if (PHG_IsDetectOnTheFly())
		DetPrintParams();

	/* Print out binning parameters if appropriate */
	if (PHG_IsBinOnTheFly()) {
		for (curBinParams = 0; curBinParams < PhgNumBinParams; curBinParams++) {
			PhgBinPrintParams(PhgRunTimeParams.PhgBinParamsFilePath[curBinParams],
				&PhgBinParams[curBinParams], &PhgBinFields[curBinParams]);
		}
	}
}

/*********************************************************************************
*
*			Name:		PhgAbort
*
*			Summary:	Abort program execution.
*			Arguments:
*				char	*abortStr		- Message describing abort condition.
*				Boolean	dumpTrkPhoton	- Request to dump current tracking photon
*
*
*			Function return: none.
*
*********************************************************************************/
void PhgAbort(char *abortStr, Boolean	dumpTrkPhoton)	
{

	/* Print the abort string */
	LbInPrintf("\n%s", abortStr);

	/* Do debug type stuff if in debug mode */
	#ifdef PHG_DEBUG
		
		/* Dump the current tracking photon if there is one */
		if (dumpTrkPhoton == true) {
			
			/* Dump the restart info */
			LbInPrintf("\n***** Photon Restart Info ******\n");
			LbInPrintf("\tRestart Slice = %ld, Restart Angle = %ld, Restart Voxel = %ld\n",
				(unsigned long)EmisListRestartSlice, 
				(unsigned long)EmisListRestartAngle, 
				(unsigned long)EmisListRestartVoxel);
				
			LbInPrintf("\n***** ******************* ******\n");
				
			/* Dump the decay */
			PhgDumpDecay(&EmisListNewDecay);
			
			/* See if we are tracking a blue photon */
			if (EmisListNumCreated == (LbEightByte)EmisListBluePhoton.number) {
			
				/* Dump the blue photon */
				PhgDumpPhoton(&EmisListBluePhoton);
			}
			else {
			
				/* Dump the pink photon */
				PhgDumpPhoton(&EmisListPinkPhoton);
			}

			/* If forced detection is on, and we have an initialized FD photon, dump it */
			if (PHG_IsForcedDetection()  && (EmisListFDPhoton.number != 0)) {
				LbInPrintf("\n\n*** Forced Detection Photon **\n");
				
				PhgDumpPhoton(&EmisListFDPhoton);
			}
		}
	#else
		if (dumpTrkPhoton) {};		/* Eliminate compiler warning */
	#endif
	
	/* Call the library abort routine */
	ErAbort(abortStr);
}

/*********************************************************************************
*
*			Name:		phgValidateParams
*
*			Summary:	Does consistency and range checks on parameters.
*			Arguments:
*
*
*			Function return: none.
*
*********************************************************************************/
Boolean	phgValidateParams()	
{
	Boolean 		okay = false;
	LbUsFourByte	typeIndex;
	
	do { /* Process Loop */
	
	/* Verify that the number of detector setups is the same as the number of collimators */
	if ((DetNumParams > 0) && (ColNumParams > 0) && (DetNumParams != ColNumParams)) {
		sprintf(phgErrString, "You must have an equal number of collimators and detectors.");
		ErStGeneric(phgErrString);
		goto FAIL;
	}
	
	/* Verify that they didn't specify history file parameters and no history file */
	if ((PHG_IsHist() == false) && (PHG_IsHistParams() == true)) {
		sprintf(phgErrString,"You specified history-file parameters but no history file name.\n"
			"Use 'history_file' in the PHG parameters file '%s' to specify a history file name.",
			PhgRunTimeParams.PhgParamFilePath);
			 
		ErStGeneric(phgErrString);
			goto FAIL;
	}
	
	/* Verify that all collimators are the same type (no check is made if there is only one specified)  */
	for (ColCurParams = 1; ColCurParams < ColNumParams; ColCurParams++) {
		if (ColRunTimeParams[0].ColType != ColRunTimeParams[ColCurParams].ColType) {
			sprintf(phgErrString, "All collimators must be of the same type, please check your collimator parameter files.");
			ErStGeneric(phgErrString);
			goto FAIL;
		}
	}
	
	/* Verify that all detectors are the same type (no check is made if there is only one specified) */
	for (DetCurParams = 1; DetCurParams < DetNumParams; DetCurParams++) {
		if (DetRunTimeParams[0].DetectorType != DetRunTimeParams[DetCurParams].DetectorType) {
			sprintf(phgErrString, "All detectors must be of the same type, please check your detectors parameter files.");
			ErStGeneric(phgErrString);
			goto FAIL;
		}
	}
	
	/* Verify that the minimum energy is not greater than the incident energy */
	if (PhgRunTimeParams.PhgNuclide.photonEnergy_KEV <
			PhgRunTimeParams.PhgMinimumEnergy) {
		sprintf(phgErrString,"You have specified an incident photon energy of %f which is less than the minimum energy of %f.\n, Please edit your parameter file, '%s' to fix this\n",
			PhgRunTimeParams.PhgNuclide.photonEnergy_KEV,
			PhgRunTimeParams.PhgMinimumEnergy,
			PhgRunTimeParams.PhgParamFilePath);
			
		ErStGeneric(phgErrString);
		goto FAIL;
	}	
	
	/* Verify that they don't try to use Forced Detection and Coherent Scatter */
	if (PHG_IsModelCoherentInObj() && PHG_IsForcedDetection()) {
		sprintf(phgErrString, "You can not use forced detection with coherent scatter modelling in the object, sorry.");
		ErStGeneric(phgErrString);
		goto FAIL;
	}
	
	/* Verify that they don't turn polarization on without PET */
	if ((PHG_IsModelPolarization()) && (PHG_IsPET() == false)) {
		sprintf(phgErrString, "You can not model polarization in SPECT, fix this in your parameters file.");
		ErStGeneric(phgErrString);
		goto FAIL;
	}
	
	/* Verify that isotope is specified if use positron range is selected */
	if ((PHG_IsRangeAdjust() == true) && (PhgRunTimeParams.PhgNuclide.isotope == PhgEn_IsotopType_NULL)) {
		LbInPrintf("In order to model positron range you must specify an isotope in the\n"
			" phg run time parameters file with the parameter 'isotope'\n"
			" Isotope is an enumerated type 'ENUM' with the following types supported\n");
		
		for (typeIndex = 1; typeIndex < (NUM_ISOTOPE_TYPES); typeIndex++){
			LbInPrintf("\n\t%s", phgEn_IsotopeStr[typeIndex]);
		}
		LbInPrintf("\n\nPlease edit your run time parameter file '%s' and try again", PhgRunTimeParams.PhgParamFilePath);
		ErStGeneric("You have selected positron range modelling without specifying an isotope\n");
		goto FAIL;
	}
	
	/* Verify that isotope path is specified if use positron range is selected */
	if ((PHG_IsRangeAdjust() == true) && (EmisListIsotopeDataFilePath[0] == '\0')) {

		LbInPrintf("In order to model positron range you must specify a path to the isotope data in the\n"
			" phg run time parameters file with the parameter 'isotope_data_file'\n"
			" The default table should be supplied in the 'phg.data' directory and should be named 'isotope_positron_range_data'\n");
		
		LbInPrintf("\n\nPlease edit your run time parameter file '%s' and try again", PhgRunTimeParams.PhgParamFilePath);
		ErStGeneric("You have selected positron range modelling without specifying a path to the positron range data\n");
		goto FAIL;
	}
	
	/* Verify that if we are doing slat collimation we have a planar detector */
	if ((ColIsSlat() == true) && !((DetRunTimeParams[DetCurParams].DetectorType == DetEn_Planar) ||
			(DetRunTimeParams[DetCurParams].DetectorType == DetEn_DualHeaded))) {
		ErStGeneric("For slat collimation you must have a planar or DHCI detector module");
		goto FAIL;
	}

	okay = true;
	FAIL:;
	} while (false);
	
	return(okay);			
}

/*********************************************************************************
*
*			Name:		PhgDumpDecay
*
*			Summary:	Perform ASCII dump of decay.
*			Arguments:
*				PHG_Decay	*decayPtr.
*
*
*			Function return: none.
*
*********************************************************************************/
void PhgDumpDecay(PHG_Decay	*decayPtr)	
{
	
	LbInPrintf("\n**** DECAY DUMP ****\n");
	LbInPrintf("\tlocation  \n");
	LbInPrintf("\t\tpos.x = %3.2f\n", decayPtr->location.x_position);
	LbInPrintf("\t\tpos.y = %3.2f\n", decayPtr->location.y_position);
	LbInPrintf("\t\tpos.z = %3.2f\n", decayPtr->location.z_position);
	LbInPrintf("\tdecay_weight = %3.3e\n", decayPtr->startWeight);
	LbInPrintf("\n*******************\n");
			
}

/*********************************************************************************
*
*			Name:		PhgDumpPhoton
*
*			Summary:	Perform ASCII dump of photon.
*			Arguments:
*				PHG_TrackingPhoton	*trkPhotonPtr.
*
*
*			Function return: none.
*
*********************************************************************************/
void PhgDumpPhoton(PHG_TrackingPhoton	*trkPhotonPtr)	
{

	LbUsTwoByte	lcv;	/* Simple loop control */
	
	LbInPrintf("\n**** PHOTON DUMP ****\n");
	LbInPrintf("\tflags = %d\n", trkPhotonPtr->flags);
	LbInPrintf("\tlocation  \n");
	LbInPrintf("\t\tpos.x = %3.2f\n", trkPhotonPtr->location.x_position);
	LbInPrintf("\t\tpos.y = %3.2f\n", trkPhotonPtr->location.y_position);
	LbInPrintf("\t\tpos.z = %3.2f\n", trkPhotonPtr->location.z_position);
	LbInPrintf("\tangle  \n");
	LbInPrintf("\t\tcos.x = %3.2f\n", trkPhotonPtr->angle.cosine_x);
	LbInPrintf("\t\tcos.y = %3.2f\n", trkPhotonPtr->angle.cosine_y);
	LbInPrintf("\t\tcos.z = %3.2f\n", trkPhotonPtr->angle.cosine_z);
	LbInPrintf("\tsliceIndex = %ld\n", (unsigned long)trkPhotonPtr->sliceIndex);
	LbInPrintf("\tangleIndex = %ld\n", (unsigned long)trkPhotonPtr->angleIndex);
	LbInPrintf("\torigSliceIndex = %ld\n", (unsigned long)trkPhotonPtr->origSliceIndex);
	LbInPrintf("\torigAngleIndex = %ld\n", (unsigned long)trkPhotonPtr->origAngleIndex);
	LbInPrintf("\txIndex = %ld\n", (unsigned long)trkPhotonPtr->xIndex);
	LbInPrintf("\tyIndex = %ld\n", (unsigned long)trkPhotonPtr->yIndex);
	LbInPrintf("\tnum_Of_Scatters = %ld\n", (unsigned long)trkPhotonPtr->num_of_scatters);
	LbInPrintf("\tphoton_scatter_weight = %3.3e\n", trkPhotonPtr->photon_scatter_weight);
	LbInPrintf("\tphoton_primary_weight = %3.3e\n", trkPhotonPtr->photon_primary_weight);
	LbInPrintf("\tphoton_current_weight = %3.3e\n", trkPhotonPtr->photon_current_weight);
	LbInPrintf("\tscatter_target_weight = %3.3e\n", trkPhotonPtr->scatter_target_weight);
	LbInPrintf("\tdecay_weight = %3.3e\n", trkPhotonPtr->decay_weight);
	LbInPrintf("\tenergy = %3.2f\n", trkPhotonPtr->energy);
	LbInPrintf("\ttravel_distance = %3.3f\n", trkPhotonPtr->travel_distance);
	LbInPrintf("\tnumStarts = %ld\n", (unsigned long)trkPhotonPtr->numStarts);
	LbInPrintf("\ttransaxialPosition = %3.3f\n", trkPhotonPtr->transaxialPosition);
	LbInPrintf("\tazimuthalAngleIndex = %3.3f\n", trkPhotonPtr->azimuthalAngleIndex);
	LbInPrintf("\taxialPosition = %3.3f\n", trkPhotonPtr->axialPosition);
	LbInPrintf("\tdetectorAngle = %3.3f\n", trkPhotonPtr->detectorAngle);
	
	/* Run through starts list */
	if (trkPhotonPtr->numStarts > 1) {
		LbInPrintf("\tstarts_list\n");
		
		for (lcv = 1; lcv <= trkPhotonPtr->numStarts; lcv++) {
			
			LbInPrintf("\t\tstart[%d] angleIndex = %ld, sliceIndex = %ld\n",
				lcv, (unsigned long)(trkPhotonPtr->starts_list[lcv-1].angleIndex),
				(unsigned long)(trkPhotonPtr->starts_list[lcv-1].sliceIndex));
		}
	}
	LbInPrintf("\tnumber = %lld\n", trkPhotonPtr->number);
	LbInPrintf("\n*******************\n");
			
}
/*********************************************************************************
*
*			Name:			phg_initialize
*
*			Summary:	Initialize the various phg modules.
*			Arguments:
*
*			Function return: true unless an error occurs.
*
*********************************************************************************/
Boolean phg_initialize()	
{
	Boolean			okay = false;	/* Process flag */
	LbFourByte		randSeed;		/* Seed for random generator */
	LbUsFourByte	curBinParams;	/* LCV for multiple binning parameters */
	
	do { /* Process Loop */
		
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
				
		/* Setup the Subobject modules */
		{
			/* Create sub objects */
			if (!SubObjCreate()) {
				break;
			}
		}
		
		/* If the user has selected to calculate the number of decays to simulate
		 from the scan time and activity, fill in the events to simulate
		 field in the header with the computed value. */
		if ( PhgRunTimeParams.PhgIsCalcEventsToSimulate ) {
			PhgRunTimeParams.Phg_EventsToSimulate = SUBOBJGetTotalRealDecays();
		}

		/* Initialize the Cylinder Positions */
		{
			/* Set object cylinder */
			if (!CylPosInitObjectCylinder()) {
				break;
			}
			
			/* Set the limit cylinder */
			CylPosInitLimitCylinder();
		}
		
		/* Setup up the productivity information */
		{
			/* Initialize productivity table info */
			if (PHG_IsStratification() == false) {
				prodTableInfo.inputFileName = "";
				prodTableInfo.outputFileName = "";
			}
			else {
				prodTableInfo.inputFileName = 
					PhgRunTimeParams.PhgProdTblInputTableFilePath;
				prodTableInfo.outputFileName = 
					PhgRunTimeParams.PhgProdTblOutputTableFilePath;
			}
				
			prodTableInfo.acceptanceAngle = 
				PhgRunTimeParams.Phg_AcceptanceAngle;
			

			/* Create productivity table */
			if (!SubObjGetStartingProdValues(&prodTableInfo)) {
				break;
			}
	
		}

		
		/* Initialize the collimation module if necessary */
		if (PHG_IsCollimateOnTheFly()) {
			for (ColCurParams = 0; ColCurParams < ColNumParams; ColCurParams++) {
				if (!ColInitialize(true)) 
					goto FAIL;
			}
			ColCurParams = 0;
		}
		
		/* Set the criticial zone.  Note, if doing cone beam SPECT collimation
			then this initialization is differnet, hence it must be performed
			after potential collimator initialization
		 */
		if (!CylPosInitCriticalZone(PhgRunTimeParams.Phg_AcceptanceAngle)) {
			break;
		}

		/* Initialize the photon tracking module */
		if (!PhoTrkInitialize()) {
			break;
		}
		
		/* Initialize the detection module if necessary */
		if (PHG_IsDetectOnTheFly()) {
			for (DetCurParams = 0; DetCurParams < DetNumParams; DetCurParams++) {
			ColCurParams = DetCurParams;
			if (!DetInitialize(true)) 
				goto FAIL;
			}
			DetCurParams = 0;
			ColCurParams = 0;
		}
		
		/* make sure that slat collimation is only done with dual_headed detectors */
		if (PHG_IsCollimateOnTheFly()) {
			for (ColCurParams = 0; ColCurParams < ColNumParams; ColCurParams++) {
				DetCurParams = ColCurParams;
				if ( ( ColIsSlat() ) && ( !DetIsDualHeaded() ) ) {
					PhgAbort("Slat collimation is only available with dual-headed detectors (phg_initialize)", false);
				}
			}
			ColCurParams = 0;
		}

		/* Initialize the binning module if necessary */
		if (PHG_IsBinOnTheFly()) {
			for (curBinParams = 0; curBinParams < PhgNumBinParams; curBinParams++) {
				if (!PhgBinInitialize(PhgRunTimeParams.PhgBinParamsFilePath[curBinParams],
						&PhgBinParams[curBinParams],
						&PhgBinData[curBinParams], &PhgBinFields[curBinParams])) {
					goto FAIL;
				}
			}
		}
		
		/* If collimation is being done be sure that the radius of the collimator is larger than the
		radius of the target cylinder
		*/
		if (PHG_IsCollimateOnTheFly() == true) {
			if (ColGtInsideRadius() < CylPosGetTargetRadius()) {
				LbInPrintf("Your target cylinder radius, %3.3f, is greater than your collimator radius %3.3f",
					CylPosGetTargetRadius(), ColGtInsideRadius());
				
				PhgAbort("Unable to perform simulation", false);
			}
		}
			
		/* If detection is being done be sure that the radius of the detector is larger than the
		radius of the target (or collimator if one is being used) cylinder
		*/
		if (PHG_IsDetectOnTheFly() == true) {
			if (ColGtOutsideRadius() > DetGtInsideRadius()) {
				LbInPrintf("Your target (or collimator) cylinder radius, %3.3f, is greater than your detector radius %3.3f",
					ColGtOutsideRadius(), DetGtInsideRadius());
				
				PhgAbort("Unable to perform simulation", false);
			}
		}
						
		okay = true;
		FAIL:;
	} while (false);

	/* Cleanup if things didn't work */
	if (!okay) {

		/* Terminate previous modules to clear any memory or file allocations */
		PhoTrkTerminate();
		ProdTblTerminate(&prodTableInfo);
		SubObjTerminate();
		PhgMathTerminate();
		EmisListTerminate();
	}	

	return (okay);
}

#ifdef DGUX
/*********************************************************************************
*
*			Name:		PhgSigHandler
*
*			Summary:	Handle floating point exception signals.
*			Arguments:
*
*			Function return: none.
*
*********************************************************************************/
void PhgSigHandler(int sigNum, siginfo_t *sigInfoPtr, ucontext_t *sigContextPtr)	
{
	char sigMessage[1024];	/* Error string */
	
	/* Build error message */
	sprintf(sigMessage, "\nA floating-point exception has occurred! %s",
		"Use your debugger to track the offending location.");
	
	/* Abort the process */
	PhgAbort(sigMessage, true);
	
}
#endif

/*********************************************************************************
*
*			Name:			phg_terminate
*
*			Summary:	Terminate the various phg modules.
*			Arguments:
*
*			Function return: none.
*
*********************************************************************************/
void phg_terminate()	
{
	LbUsFourByte curBinParams;
	
	/* Print the collimation report if initialized */
	if (!ErIsInError()){
		if (PHG_IsCollimateOnTheFly()) {
			ColPrintReport();
		}
	
		/* Print the detection report if initialized */
		if (PHG_IsDetectOnTheFly()) {
			DetPrintReport();
		}
	
		/* Print the binning report if initialized */
		if (PHG_IsBinOnTheFly()) {
			for (curBinParams = 0; curBinParams < PhgNumBinParams; curBinParams++) {
				PhgBinPrintReport(&PhgBinParams[curBinParams], &PhgBinFields[curBinParams]);
			}
		}
	}
	
	/* Terminate the binning module if initialized */
	if (PHG_IsBinOnTheFly()) {
		for (curBinParams = 0; curBinParams < PhgNumBinParams; curBinParams++) {
			PhgBinTerminate(&PhgBinParams[curBinParams], &PhgBinData[curBinParams],
				&PhgBinFields[curBinParams]);
		}
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
	ProdTblTerminate(&prodTableInfo);
	
	/* Terminate the sub object module */
	SubObjTerminate();
	
}
#undef PHG_MAIN
