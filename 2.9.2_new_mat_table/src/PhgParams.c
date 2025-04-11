/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1993-2013 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			PhgParams.c
*     Revision Number:		1.4
*     Date last revised:	24 July 2013
*     Programmer:			Steven Vannoy
*     Date Originated:		29 July 1993
*
*     Module Overview:	The parameter module.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*			PhgLookupRunTimeParamLabel
*			PhgGetRunTimeParams
*			PhgLookupBinParamLabel
*			PhgGetBinParams
*
*     Global variables defined:
*
*
*	  Global macros defined:
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
*						- binning by random-state and crystal number.
*						- binning PET data as SPECT
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:	Added support for time-of-flight.
*
*********************************************************************************/

#define PHG_PARAMS

#include <stdio.h>
#include <string.h>
#include <time.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "PhgMath.h"
#include "CylPos.h"
#include "ProdTbl.h"
#include "SubObj.h"
#include "EmisList.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "PhoHFile.h"
#include "phg.h"
#include "PhgBin.h"

/* LOCAL CONSTANTS */
/* LOCAL TYPES */
typedef char labelTy[LBPF_LABEL_LEN];

/* LOCAL GLOBALS */
static 	char		phgBinErrString[512];				/* Storage for error string s*/

static	LbPfHkTy	phgParamFileHk;				/* Our parameter file */

/* When changing the following list also change PhgEn_RunTimeParamsTy in PhgParams.h.
The two lists must have the same order */
static	labelTy		phgRunTimeParamLabels[] = {	/* The parameter labels */
					"simulate_stratification",
					"simulate_forced_detection",
					"forced_non_absorbtion",	/* see also forced_non_absorption below */
					"acceptance_angle",
					"num_to_simulate",
					"simulate_SPECT",
					"simulate_PET_coincidences_only",
					"simulate_PET_coincidences_plus_singles",
					"simulate_Multi_Emission_Isotope",
					"adjust_for_positron_range",
					"adjust_for_collinearity",
					"model_coherent_scatter",
					"model_coherent_scatter_in_tomo",
					"model_coherent_scatter_in_obj",
					"minimum_energy",
					"photon_energy",
					"weight_window_ratio",
					"point_source_voxels",
					"line_source_voxels",
					"random_seed",
					"isotope",
					"isotope_data_file",
					"length_of_scan",
					"model_polarization",
					"object",
					"num_slices",
					"slice",
					"slice_number",
					"zMin",
					"zMax",
					"xMin",
					"xMax",
					"yMin",
					"yMax",
					"num_X_bins",
					"num_Y_bins",
					"num_act_X_bins",
					"num_act_Y_bins",
					"num_att_X_bins",
					"num_att_Y_bins",
					"target_cylinder",
					"target_zMin",
					"target_zMax",
					"centerX",
					"centerY",
					"radius",
					"activity_indexes",
					"activity_table",
					"activity_index_trans",
					"activity_image",
					"attenuation_indexes",
					"attenuation_table",
					"attenuation_index_trans",
					"attenuation_image",
					"coherent_scatter_table",
					"productivity_input_table",
					"productivity_output_table",
					"forced_detection_table",
					"statistics_file",
					"bin_params_file",
					"collimator_params_file",
					"detector_params_file",
					"tomograph_params_file",
					"history_file",
					"history_params_file",
					"forced_non_absorption",	/* correct spelling!, replacing forced_non_absorbtion */
					""};

/* When changing the following list also change PhgEn_BinParamsTy in PhgParams.h.
The two lists must have the same order */
static	labelTy	phgBinParamLabels[] = {
					"num_z_bins",
					"num_pa_bins",
					"num_td_bins",
					"num_aa_bins",
					"num_e_bins",
					"num_e1_bins",
					"num_e2_bins",
					"scatter_param",
					"scatter_random_param",
					"accept_randoms",
					"num_theta_bins",
					"num_phi_bins",
					"num_xr_bins",
					"num_yr_bins",
					"num_tof_bins",
					"min_z",
					"max_z",
					"min_td",
					"max_td",
					"min_pa",
					"max_pa",
					"min_e",
					"max_e",
					"min_s",
					"max_s",
					"max_theta",
					"min_xr",
					"max_xr",
					"min_yr",
					"max_yr",
					"min_tof",
					"max_tof",
					"bin_by_crystal",
					"sum_according_to_type",
					"weight_image_type",
					"count_image_type",
					"add_to_existing_img",
					"weight_image_path",
					"weight_squared_image_path",
					"count_image_path",
					"do_ssrb",
					"do_msrb",
					"detector_radius",
					"history_file",
					"history_params_file",
					"binPETasSPECT",
					""};

LbUsFourByte			*dumPtr;

/* PROTOTYPES */
Boolean phgGetTomoParams(void);


/* FUNCTIONS */

/*********************************************************************************
*
*			Name:			phgGetTomoParams
*
*			Summary:	Read in the tomograph parameters. This routine may get called
*						with the initial run time parameters file, or it may get called
*						with the "tomograph" file, depending on whether the user supplies
*						a tomograph file or not.
*			Arguments:
*				FILE	*phgParamFileHk	- Either the run time parameters or the tomograph file
*
*			Function return: true unless an error occurs.
*
*********************************************************************************/
Boolean phgGetTomoParams()	
{
	double					paramBuffer[LBPF_PARAM_LEN];
	Boolean					okay = false;					/* Process flag */
	Boolean					isEOF;							/* End of file flag */
	Boolean					switchOkay;						/* Switch flag */
	LbPfEnPfTy				paramType;
	LbUsFourByte			paramSize;
	char					paramLabel[LBPF_LABEL_LEN];
	PhgEn_RunTimeParamsTy	whichParam;
	LbUsFourByte			loopV =0;
	LbPfHkTy				tomoParamFileHk;
	
	
	/* Process all the tomograph files */
	for (loopV = 0; loopV < PhgNumTomoFiles; loopV++ ) { /* Process Loop */
	
		/* Attempt to open parameter file */
		if (LbPfOpen(PhgRunTimeParams.PhgTomoParamsFilePath[loopV], 0, &tomoParamFileHk) == false) {
			sprintf(phgBinErrString, "An error occurred opening the tomograph parameters file '%s'\n",
				PhgRunTimeParams.PhgTomoParamsFilePath[loopV]);
			ErAlert(phgBinErrString, false);
			break;
		}
	
		while(1) {
			if (LbPfGetParam(&tomoParamFileHk, (void *)paramBuffer,
				&paramType, &paramSize, paramLabel, &isEOF) == false) {
				
				break;
			}
				
				/* Find the runtime parameter */
				whichParam = PhgLookupRunTimeParamLabel(paramLabel);
				
				switchOkay = true;
				switch (whichParam) {
					
					case PhgEn_bin_params_file:
					
							/* Allocate memory for binning records */
							strcpy(PhgRunTimeParams.PhgBinParamsFilePath[PhgNumBinParams],
								(char *) paramBuffer);
							
							if (PhgRunTimeParams.PhgBinParamsFilePath[PhgNumBinParams][0] != '\0') {
								PhgRunTimeParams.PhgIsBinOnTheFly = true;
								PhgNumBinParams++;
							}
							
							/* Verify we don't go over the maximum number allowed */
							if (PhgNumBinParams == PHG_MAX_PARAM_FILES) {
								sprintf(phgBinErrString, "(PhgGetTomoParams) You have exceeded the maximum number of binning parameter files (%d).\n",
		                            PHG_MAX_PARAM_FILES);
								ErStGeneric(phgBinErrString);
								goto FAIL;
							};
							
							
						break;
					
					case PhgEn_collimator_params_file:
							strcpy(PhgRunTimeParams.PhgCollimatorParamsFilePath[ColNumParams],
								(char *) paramBuffer);
							
							if (PhgRunTimeParams.PhgCollimatorParamsFilePath[ColNumParams][0] != '\0') {
								PhgRunTimeParams.PhgIsCollimateOnTheFly = true;
								ColNumParams += 1;
							}
						break;
					
					case PhgEn_detector_params_file:
							strcpy(PhgRunTimeParams.PhgDetectorParamsFilePath[DetNumParams],
								(char *) paramBuffer);
							
							if (PhgRunTimeParams.PhgDetectorParamsFilePath[DetNumParams][0] != '\0') {
								PhgRunTimeParams.PhgIsDetectOnTheFly = true;
								DetNumParams += 1;
							}
						break;


					case PhgEn_NULL:
						sprintf(phgBinErrString, "(PhgGetTomoParams) Unknown (hence unused) parameter (%s).\n",
                            paramLabel);
						ErAlert(phgBinErrString, false);
						break;
					
					default:
						sprintf(phgBinErrString, "(PhgGetTomoParams) Unknown (hence unused) parameter (%s).\n",
                            paramLabel);
						ErAlert(phgBinErrString, false);
						break;
				}
				if (!switchOkay)
					break;
			}
		
		okay = true;
	}
	FAIL:;
	
	return(okay);
}

/*********************************************************************************
*
*			Name:			PhgGetRunTimeParams
*
*			Summary:	Read in the runtime parameters.
*			Arguments:
*
*			Function return: true unless an error occurs.
*
*********************************************************************************/
Boolean PhgGetRunTimeParams()	
{
	double					paramBuffer[LBPF_PARAM_LEN];
	Boolean					okay = false;					/* Process flag */
	Boolean					isEOF;							/* End of file flag */
	Boolean					switchOkay;						/* Switch flag */
	LbPfEnPfTy				paramType;
	LbUsFourByte			paramSize;
	LbUsFourByte			typeIndex;
	char					paramLabel[LBPF_LABEL_LEN];
	PhgEn_RunTimeParamsTy	whichParam;
	LbUsFourByte			loopV =0;

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
		
		/* Set the following options to false: the user must set exactly one to true
			the program will halt with error */
		PhgRunTimeParams.PhgIsSPECT = false;
		PhgRunTimeParams.PhgIsPETCoincidencesOnly = false;
		PhgRunTimeParams.PhgIsPETCoincPlusSingles = false;
		PhgRunTimeParams.PhgIsMultiEmission = false;
		
		/* These fields are set true as appropriate when the above are set */
		PhgRunTimeParams.PhgIsPET = false;
		PhgRunTimeParams.PhgIsPETCoincidences = false;
		PhgRunTimeParams.PhgIsPETSingles = false;
		
		
		
		/* Attempt to open parameter file */
		if (!LbPfOpen(PhgRunTimeParams.PhgParamFilePath, 0, &phgParamFileHk)) {
			sprintf(phgBinErrString, "An error occurred opening the command line input parameters file '%s'\n",
				PhgRunTimeParams.PhgParamFilePath);
			ErAlert(phgBinErrString, false);
			break;
		}
		while(LbPfGetParam(&phgParamFileHk, (void *)paramBuffer,
				&paramType, &paramSize, paramLabel, &isEOF)) {
				
				/* Find the runtime parameter */
				whichParam = PhgLookupRunTimeParamLabel(paramLabel);
				
				switchOkay = true;
				switch (whichParam) {


					case PhgEn_simulate_stratification:
						PhgRunTimeParams.PhgIsStratification =
							*((Boolean *) paramBuffer);
						break;
						
					case PhgEn_simulate_forced_detection:
						PhgRunTimeParams.PhgIsForcedDetection =
							*((Boolean *) paramBuffer);
						break;
						
					case PhgEn_forced_non_absorption:		/* fixed spelling in parameter file, though I've left the wrong spelling in PhgRunTimeParams */
					case PhgEn_forced_non_absorbtion:
						PhgRunTimeParams.PhgIsForcedNonAbsorbtion =
							*((Boolean *) paramBuffer);
						break;
						
					case PhgEn_acceptance_angle:
						PhgRunTimeParams.Phg_AcceptanceAngle = 
							*((double *) paramBuffer);
								
						/* Set sine of acceptance angle */
						PhgRunTimeParams.Phg_SineOfAcceptanceAngle = 
							PHGMATH_SinFromDegrees(PhgRunTimeParams.Phg_AcceptanceAngle);

						break;
						
					case PhgEn_num_to_simulate:
						if (paramSize == 8) {
							PhgRunTimeParams.Phg_EventsToSimulate = 
								*((LbUsEightByte *) paramBuffer);
						} else {
							PhgRunTimeParams.Phg_EventsToSimulate = 
								*((LbUsFourByte *) paramBuffer);
						}
						
						/* if the user specified 0 events to simulate, SimSET will calculate
						 the number of decays to simulate */
						PhgRunTimeParams.PhgIsCalcEventsToSimulate = 
								(PhgRunTimeParams.Phg_EventsToSimulate == 0);
								
						break;
						
					case PhgEn_simulate_SPECT:
						PhgRunTimeParams.PhgIsSPECT = *((Boolean *) paramBuffer);
						
						/* check that no other simulation modalities are set */
						if ( 	PhgRunTimeParams.PhgIsSPECT && (PhgRunTimeParams.PhgIsPETCoincidencesOnly ||
								PhgRunTimeParams.PhgIsPETCoincPlusSingles || 
								PhgRunTimeParams.PhgIsMultiEmission)
							) {
							
							sprintf(phgBinErrString, "User supplied multiple modalities.  Only one of"
								" simulate_SPECT, simulate_PET_coincidences_only,\n"
								"simulate_PET_coincidences_plus_singles, and simulate_PET_singles_only"
								" may be set to true.\n");
								
							ErStGeneric(phgBinErrString);
							goto FAIL;
						
						}
						
						break;
						
					case PhgEn_simulate_PETCoincidencesOnly:
						PhgRunTimeParams.PhgIsPETCoincidencesOnly = *((Boolean *) paramBuffer);
						
						/* check that no other simulation modalities are set */
						if ( 	PhgRunTimeParams.PhgIsPETCoincidencesOnly && (PhgRunTimeParams.PhgIsSPECT || 
								PhgRunTimeParams.PhgIsPETCoincPlusSingles ||
								PhgRunTimeParams.PhgIsMultiEmission)
							) {
							
							sprintf(phgBinErrString, "User supplied multiple modalities.  Only one of"
								" simulate_SPECT, simulate_PET_coincidences_only,\n"
								"simulate_PET_coincidences_plus_singles, and simulate_PET_singles_only"
								" may be set to true.\n");
								
							ErStGeneric(phgBinErrString);
							goto FAIL;
						
						}
						
						/* set other PET coincidence flags */
						if ( PhgRunTimeParams.PhgIsPETCoincidencesOnly ) {
							
							PhgRunTimeParams.PhgIsPET = true;
							PhgRunTimeParams.PhgIsPETCoincidences = true;
						}
						
						break;
						
					case PhgEn_simulate_PETCoincPlusSingles:
						PhgRunTimeParams.PhgIsPETCoincPlusSingles = *((Boolean *) paramBuffer);
						
						/* check that no other simulation modalities are set */
						if ( 	PhgRunTimeParams.PhgIsPETCoincPlusSingles && (PhgRunTimeParams.PhgIsSPECT || 
								PhgRunTimeParams.PhgIsPETCoincidencesOnly || 
								PhgRunTimeParams.PhgIsMultiEmission)
							) {
							
							sprintf(phgBinErrString, "User supplied multiple modalities.  Only one of"
								" simulate_SPECT, simulate_PET_coincidences_only,\n"
								"simulate_PET_coincidences_plus_singles, and simulate_PET_singles_only"
								" may be set to true.\n");
								
							ErStGeneric(phgBinErrString);
							goto FAIL;
						
						}
						
						/* set other PET coincidence and singles flags */
						if ( PhgRunTimeParams.PhgIsPETCoincPlusSingles ) {
							
							PhgRunTimeParams.PhgIsPET = true;
							PhgRunTimeParams.PhgIsPETCoincidences = true;
							PhgRunTimeParams.PhgIsPETSingles = true;
						}
						
						break;
						
					case PhgEn_simulate_MultiEmission:
						PhgRunTimeParams.PhgIsMultiEmission = *((Boolean *) paramBuffer);
						
						/* this option is not yet supported */
						if ( PhgRunTimeParams.PhgIsMultiEmission ) {
							
							sprintf(phgBinErrString, "User set simulate_Multi_Emission_Isotope to true.\n"
								"This option is not yet supported.\n");
								
							ErStGeneric(phgBinErrString);
							goto FAIL;
						
						}						
						
						/* check that no other simulation modalities are set */
						if ( 	PhgRunTimeParams.PhgIsMultiEmission && (PhgRunTimeParams.PhgIsSPECT || 
								PhgRunTimeParams.PhgIsPETCoincidencesOnly || PhgRunTimeParams.PhgIsPETCoincPlusSingles)
							) {
							
							sprintf(phgBinErrString, "User supplied multiple modalities.  Only one of"
								" simulate_SPECT, simulate_PET_coincidences_only,\n"
								"simulate_PET_coincidences_plus_singles, and simulate_PET_singles_only"
								" may be set to true.\n");
								
							ErStGeneric(phgBinErrString);
							goto FAIL;
						
						}
						
						break;
						
					case PhgEn_adjust_for_positron_range:
						PhgRunTimeParams.PhgIsAdjForPosRange = 
							*((Boolean *) paramBuffer);
						break;
						
					case PhgEn_adjust_for_collinearity:
						PhgRunTimeParams.PhgIsAdjForCollinearity = *((Boolean *) paramBuffer);
						break;
						
					case PhgEn_model_coherent:
						PhgRunTimeParams.PhgIsModelCoherent = *((Boolean *) paramBuffer);
						PhgRunTimeParams.PhgIsModelCoherentInTomo = *((Boolean *) paramBuffer);
						PhgRunTimeParams.PhgIsModelCoherentInObj = *((Boolean *) paramBuffer);
						break;
						
					case PhgEn_model_coherent_in_tomo:
						PhgRunTimeParams.PhgIsModelCoherentInTomo = *((Boolean *) paramBuffer);
						break;
						
					case PhgEn_model_coherent_in_obj:
						PhgRunTimeParams.PhgIsModelCoherentInObj = *((Boolean *) paramBuffer);
						break;
						
					case PhgEn_minimum_energy:
						PhgRunTimeParams.PhgMinimumEnergy = *((double *) paramBuffer);
						
						/* Do a verification on the energy */
						if (PHG_MIN_PHOTON_ENERGY > PhgRunTimeParams.PhgMinimumEnergy) {
							sprintf(phgBinErrString, "User supplied minimum energy (%3.3f) is"
								" less than minimum supported energy (%3.3f).\n"
								" Please edit the parameter file to raise the minimum energy.\n",
								PhgRunTimeParams.PhgMinimumEnergy, PHG_MIN_PHOTON_ENERGY);
								
							ErStGeneric(phgBinErrString);
							goto FAIL;
						}
						break;
						
					case PhgEn_photon_energy:
						PhgRunTimeParams.PhgNuclide.photonEnergy_KEV = *((double *) paramBuffer);
						break;
					
					case PhgEn_weight_window_ratio:
						PhgRunTimeParams.PhgMaxWWRatio = *((double *) paramBuffer);
						PhgRunTimeParams.PhgMinWWRatio = 1/PhgRunTimeParams.PhgMaxWWRatio;
						break;
					
					case PhgEn_point_source_voxels:
						PhgRunTimeParams.PhgIsVoxelPointSource = 
							*((Boolean *) paramBuffer);
						break;
					
					case PhgEn_line_source_voxels:
						PhgRunTimeParams.PhgIsVoxelLineSource = 
							*((Boolean *) paramBuffer);
						break;
					
					case PhgEn_random_seed:
						PhgRunTimeParams.PhgRandomSeed = 
							*((LbUsFourByte *) paramBuffer);
						break;

					case SubObjEn_object:
						/* Let subobject code read in the geometry */
						switchOkay = SubObjGetObjGeometry(phgParamFileHk,
							*((LbUsFourByte *) paramBuffer));
						break;

					case CylPosEn_target_cylinder:
						switchOkay = CylPosInitTargetCylinder(phgParamFileHk,
							*((LbUsFourByte *) paramBuffer));
						break;

					case SubObjEn_activity_indexes:
							strcpy(PhgRunTimeParams.PhgSubObjActivityIndexFilePath,
								(char *) paramBuffer);
						break;

					case SubObjEn_activity_table:
							strcpy(PhgRunTimeParams.PhgSubObjActivityTableFilePath,
								(char *) paramBuffer);
						break;

					case SubObjEn_coherent_scatter_table:
							strcpy(PhgRunTimeParams.PhgSubObjCohScatTableFilePath,
								(char *) paramBuffer);
						break;

					case SubObjEn_activity_index_trans:
							strcpy(PhgRunTimeParams.PhgSubObjActIndexTransFilePath,
								(char *) paramBuffer);
						break;

					case SubObjEn_activity_image:
							strcpy(PhgRunTimeParams.PhgSubObjActImgFilePath,
								(char *) paramBuffer);
						break;

					case SubObjEn_attenuation_indexes:
							strcpy(PhgRunTimeParams.PhgSubObjAttenIndexFilePath,
								(char *) paramBuffer);
						break;

					case SubObjEn_attenuation_table:
							strcpy(PhgRunTimeParams.PhgSubObjAttenTableFilePath,
								(char *) paramBuffer);
						break;

					case SubObjEn_attenuation_index_trans:
							strcpy(PhgRunTimeParams.PhgSubObjAttIndexTransFilePath,
								(char *) paramBuffer);
						break;

					case SubObjEn_attenuation_image:
							strcpy(PhgRunTimeParams.PhgSubObjAttImgFilePath,
								(char *) paramBuffer);
						break;

					case ProdTblEn_productivity_input_table:
							strcpy(PhgRunTimeParams.PhgProdTblInputTableFilePath,
								(char *) paramBuffer);
								
							PhgRunTimeParams.PhgIsComputedProductivityTbl = 
								(PhgRunTimeParams.PhgProdTblInputTableFilePath[0] != '\0');
						break;

					case ProdTblEn_productivity_output_table:
							strcpy(PhgRunTimeParams.PhgProdTblOutputTableFilePath,
								(char *) paramBuffer);
						break;

					case PhoTrkEn_forced_detection_table:
							strcpy(PhgRunTimeParams.PhgPhoTrkForcedDetectionFilePath,
								(char *) paramBuffer);
						break;

					
					case PhoHStatEn_statistics_file:
							strcpy(PhgRunTimeParams.PhgPhoHStatFilePath,
								(char *) paramBuffer);
					break;
					
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
					
					case PhgEn_length_of_scan:

						if (paramType == LbPfEnReal) {
							PhgRunTimeParams.Phg_LengthOfScan =
								(float)( *((double *) paramBuffer) );
						} else {
							PhgRunTimeParams.Phg_LengthOfScan =
								(float)( *((LbUsFourByte *) paramBuffer) );
						}
						
						break;
					
					case PhgEn_model_polarization:
							PhgRunTimeParams.PhgIsModelPolarization =
								*((Boolean *) paramBuffer);
						break;
					
					case PhgEn_bin_params_file:
					
							/* Verify a tomograph file hasn't already been specified */
							if ((PhgNumTomoFiles > 0) && (strlen((char *)paramBuffer) > 0)) {
								ErStGeneric("You can not specify a binning paramters file and a tomograph file in your run time parameter file.");
								goto FAIL;
							}
					
							/* Allocate memory for binning records */
							strcpy(PhgRunTimeParams.PhgBinParamsFilePath[PhgNumBinParams],
								(char *) paramBuffer);
							
							PhgRunTimeParams.PhgIsBinOnTheFly |= 
								(strlen(PhgRunTimeParams.PhgBinParamsFilePath[PhgNumBinParams]) != '\0');
								
							if (PhgRunTimeParams.PhgIsBinOnTheFly) {
								PhgNumBinParams++;
							}
							
							/* Verify we don't go over the maximum number allowed */
							if (PhgNumBinParams == PHG_MAX_PARAM_FILES) {
								sprintf(phgBinErrString, "(PhgGetRunTimeParams) You have exceeded the maximum number of binning parameter files (%d).\n",
		                            PHG_MAX_PARAM_FILES);
								ErStGeneric(phgBinErrString);
								goto FAIL;
							};
							
							
						break;
					
					case PhgEn_collimator_params_file:
					
							/* Verify a tomograph file hasn't already been specified */
							if ((PhgNumTomoFiles > 0) && (strlen((char *)paramBuffer) > 0)) {
								ErStGeneric("You can not specify a collimator paramters file and a tomograph file in your run time parameter file.");
								goto FAIL;
							}
							
							strcpy(PhgRunTimeParams.PhgCollimatorParamsFilePath[ColNumParams],
								(char *) paramBuffer);
							
							PhgRunTimeParams.PhgIsCollimateOnTheFly = 
								(PhgRunTimeParams.PhgCollimatorParamsFilePath[0][0] != '\0');
								
							if (PhgRunTimeParams.PhgIsCollimateOnTheFly) {
								ColNumParams += 1;
							}
						break;
					
					case PhgEn_detector_params_file:
					
							/* Verify a tomograph file hasn't already been specified */
							if ((PhgNumTomoFiles > 0) && (strlen((char *)paramBuffer) > 0)) {
								ErStGeneric("You can not specify a detector paramters file and a tomograph file in your run time parameter file.");
								goto FAIL;
							}
							
							strcpy(PhgRunTimeParams.PhgDetectorParamsFilePath[DetNumParams],
								(char *) paramBuffer);
							
							PhgRunTimeParams.PhgIsDetectOnTheFly = 
								(PhgRunTimeParams.PhgDetectorParamsFilePath[0][0] != '\0');
								
							if (PhgRunTimeParams.PhgIsDetectOnTheFly) {
								DetNumParams += 1;
							}
						break;
					
					case PhgEn_tomograph_params_file:
							if (strlen((char *) paramBuffer) > 0) {
							
								/* First verify that old tomograph parameters where not specified
									in the file as well as a tomograph path. If they specify a tomograph file
									then all tomograph parameters must be provided there
									
								*/
								if (strlen(PhgRunTimeParams.PhgCollimatorParamsFilePath[0]) != 0) {
									ErStGeneric("You can not specify a collimator paramters file and a tomograph file in your run time parameter file.");
									goto FAIL;
								}
								if (strlen(PhgRunTimeParams.PhgDetectorParamsFilePath[0]) != 0) {
									ErStGeneric("You can not specify a detector paramters file and a tomograph file in your run time parameter file.");
									goto FAIL;
								}
								if (strlen(PhgRunTimeParams.PhgBinParamsFilePath[0]) != 0) {
									ErStGeneric("You can not specify a binning paramters file and a tomograph file in your run time parameter file.");
									goto FAIL;
								}
								
								/* Since there were no conflicting parameters, copy the tomograph parameters file name */
								strcpy(PhgRunTimeParams.PhgTomoParamsFilePath[PhgNumTomoFiles],
									(char *) paramBuffer);
								
								/* Increment global counter for number of parameter files */
								PhgNumTomoFiles += 1;
							}
						break;
					
					
					case PhgEn_isotope_data_file:
							strcpy(EmisListIsotopeDataFilePath,
								(char *) paramBuffer);
							
						break;
					
					case PhgEn_isotope:
						for (typeIndex = 0; typeIndex < (NUM_ISOTOPE_TYPES); typeIndex++){
						
							/* See if it matches */
							if (strcmp((char *)paramBuffer, phgEn_IsotopeStr[typeIndex]) == 0) {
								PhgRunTimeParams.PhgNuclide.isotope = (PhgEn_IsotopeTypeTy) typeIndex;
								break;
							}
						}
						if (typeIndex == NUM_ISOTOPE_TYPES) {
							LbInPrintf("Invalid isotope type in phg parameters '%s', valid types are:\n",
								(char *)paramBuffer);
							for (typeIndex = 1; typeIndex < (NUM_ISOTOPE_TYPES); typeIndex++){
								LbInPrintf("'%s'\n", phgEn_IsotopeStr[typeIndex]);
							}
							ErStGeneric("An invalid isotope type was supplied");
							switchOkay = false;
						}						
						break;


					case PhgEn_NULL:
						sprintf(phgBinErrString, "(PhgGetRunTimeParams) Unknown (hence unused) parameter (%s).\n",
                            paramLabel);
						ErAlert(phgBinErrString, false);
						break;
					
					default:
						sprintf(phgBinErrString, "(PhgGetRunTimeParams) Unknown (hence unused) parameter (%s).\n",
                            paramLabel);
						ErAlert(phgBinErrString, false);
						break;
								}
				if (!switchOkay)
					break;
		}
		/* Close the parameter file */
		LbPfClose(&phgParamFileHk);
		
		/* See if we quit due to error */
		if (!isEOF)
			break;
			
		/* Make sure they picked a modality */
		if ( !( PHG_IsSPECT() || PHG_IsPETCoincidencesOnly() ||
				PHG_IsPETCoincPlusSingles() || PHG_IsMultiEmission() ) ) {
			
			sprintf(phgBinErrString, "User must choose a modality.  One of"
				" simulate_SPECT, simulate_PET_coincidences_only,\n"
				"and simulate_PET_coincidences_plus_singles must be set to true.\n");
				
			ErStGeneric(phgBinErrString);
			goto FAIL;
		
		}
		
		/* If they supplied a tomograph file, then process it */
		if (PhgNumTomoFiles > 0) {
			if (phgGetTomoParams() == false) {
				goto FAIL;
			}
		}
		
		okay = true;
		FAIL:;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			PhgLookupRunTimeParamLabel
*
*			Summary:	Find the given label in the label table and return
*						its enumeration value.
*			Arguments:
*				char		*label		- The label to find.
*			Function return: Enumerated type corresponding to label.
*
*********************************************************************************/
PhgEn_RunTimeParamsTy  PhgLookupRunTimeParamLabel(char *label)	
{
	PhgEn_RunTimeParamsTy	whichParam;			/* Parameter we find */
	LbUsFourByte			labelIndex;			/* Label we are checking */

	/* Set param to null */
	whichParam = PhgEn_NULL;
	
	/* Search table for given label */
	for (labelIndex = 0; labelIndex < (PhgEn_NULL); labelIndex++){
	
		/* See if it matches */
		if (strcmp(label, phgRunTimeParamLabels[labelIndex]) == 0) {
			whichParam = (PhgEn_RunTimeParamsTy) labelIndex;
			break;
		}
	}
	
	return (whichParam);
}

/*********************************************************************************
*
*			Name:			PhgGetBinParams
*
*			Summary:	Read in the binning parameters.
*			Arguments:
*				char			*ParamsName	- The name of the parameter file to open.
*				PHG_BinparamsTy	*binParams 	- The parameter struct to initialize.
*			Function return: true unless an error occurs.
*
*********************************************************************************/
Boolean PhgGetBinParams(char *ParamsName, PHG_BinParamsTy *binParams)	
{
	double					paramBuffer[LBPF_PARAM_LEN];
	Boolean					okay = false;					/* Process flag */
	Boolean					isEOF;							/* End of file flag */
	Boolean					switchOkay;						/* Switch flag */
	LbPfEnPfTy				paramType;
	LbUsFourByte			paramSize;
	LbUsFourByte			dimensionCount = 0;				/* Current dimension */
	char					paramLabel[LBPF_LABEL_LEN];
	PhgEn_BinParamsTy		whichParam;
	LbUsFourByte			curBinParams;	/* LCV for multiple binning parameters */

	do { /* Process Loop */

		/* Clear out path name variables in case of error */
 		memset(binParams->weightImgFilePath, '\0', PATH_LENGTH);
	    memset(binParams->weightSquImgFilePath, '\0', PATH_LENGTH);
 		memset(binParams->countImgFilePath, '\0', PATH_LENGTH);

		/* Clear new param  in case its not set (user has old param file) */
		binParams->addToExistingImg = false;

		/* Set count image type to 4 byte int in case they are using an old
		   param file that doesn't specify it
		*/
		binParams->count_image_type = 2;
		
		/* In PET random coincidences are rejected by default */
		binParams->acceptRandoms  = false;

		/* Clear various parameters */

		binParams->numZBins = 0;
		binParams->numPABins = 0;
		binParams->numTDBins = 0;
		binParams->numAABins = 0;
		binParams->numTOFBins = 0;
		binParams->numE1Bins = 0;
		binParams->numE2Bins = 0;
		binParams->numS1Bins = 0;
		binParams->numS2Bins = 0;
		binParams->numPHIBins = 0;
		binParams->numThetaBins = 0;
		binParams->numXRBins = 0;
		binParams->numYRBins = 0;
		binParams->numCrystalBins = 0;
		binParams->doCrystal = false;
		binParams->doSSRB = false;
		binParams->doMSRB = false;
		binParams->sumAccordingToType = false;
		binParams->minTheta = PHGMATH_DegreesFromRadians(-PHGMATH_PI_DIV2);
		binParams->maxTheta = PHGMATH_DegreesFromRadians(PHGMATH_PI_DIV2);
		binParams->isHistoryFile = false;
		binParams->isHistoryParamsFile = false;
		binParams->isBinPETasSPECT = false;
		
		/* Attempt to open parameter file */
		if (!LbPfOpen(ParamsName, 0, &phgParamFileHk)) {
			sprintf(phgBinErrString, "An error occurred opening the binning parameters file '%s'\n"
				"Check your PHG parameters file, '%s' for a valid path\n", ParamsName,
				PhgRunTimeParams.PhgParamFilePath);
			ErAlert(phgBinErrString,false);
			 
			break;
		}
		while(LbPfGetParam(&phgParamFileHk, (void *)paramBuffer,
				&paramType, &paramSize, paramLabel, &isEOF)) {
				
				/* Find the runtime parameter */
				whichParam = PhgLookupBinParamLabel(paramLabel);
				
				switchOkay = true;
				switch (whichParam) {


					case PhgBinEn_num_z_bins:
						binParams->numZBins = *((LbUsFourByte *) paramBuffer);
						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_Z1;
						dimensionCount++;
						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_Z2;
						dimensionCount++;
						break;
						
					case PhgBinEn_num_pa_bins:
						binParams->numPABins = *((LbUsFourByte *) paramBuffer);
						break;
						
					case PhgBinEn_num_td_bins:
						binParams->numTDBins = *((LbUsFourByte *) paramBuffer);
						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_TD;
						dimensionCount++;
						break;
						
					case PhgBinEn_num_aa_bins:
						binParams->numAABins = *((LbUsFourByte *) paramBuffer);
						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_AA;
						dimensionCount++;
						break;
						
					case PhgBinEn_num_tof_bins:
						binParams->numTOFBins = *((LbUsFourByte *) paramBuffer);
						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_TOF;
						dimensionCount++;
						break;
						
					case PhgBinEn_num_e_bins:
					case PhgBinEn_num_e1_bins:
						binParams->numE1Bins = *((LbUsFourByte *) paramBuffer);
						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_Energy1;
						dimensionCount++;
						binParams->numE2Bins = *((LbUsFourByte *) paramBuffer);
						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_Energy2;
						dimensionCount++;
						break;
						
					case PhgBinEn_num_e2_bins:
						binParams->numE2Bins = *((LbUsFourByte *) paramBuffer);
						
						/* If they gave us bins for this, and its SPECT, let them know */
						if ( (PHG_IsSPECT() || binParams->isBinPETasSPECT) && (binParams->numE2Bins > 1)) {
							ErAlert("\nYou have specified numE2bins > 1 for SPECT binning, this value is ignored.\n", false);
						}
						break;
					
					/* PhgBinEn_scatter_random_param replaced PhgBinEn_scatter_param, though the latter is still understood */
					case PhgBinEn_scatter_random_param:
					case PhgBinEn_scatter_param:
						binParams->scatterRandomParam = *((LbUsFourByte *) paramBuffer);
						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_Scatter1;
						dimensionCount++;
						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_Scatter2;
						dimensionCount++;
						break;
						
					case PhgBinEn_accept_randoms:
						binParams->acceptRandoms  = 
							*((Boolean *) paramBuffer);
						break;

					case PhgBinEn_min_z:
						binParams->minZ = *((double *) paramBuffer);
						break;
						
					case PhgBinEn_max_z:
						binParams->maxZ = *((double *) paramBuffer);
						break;
						
					case PhgBinEn_min_td:
						binParams->minTD = *((double *) paramBuffer);
						break;
						
					case PhgBinEn_max_td:
						binParams->maxTD = *((double *) paramBuffer);
						break;
					
					case PhgBinEn_min_tof:
						binParams->minTOF = *((double *) paramBuffer);
						break;
					
					case PhgBinEn_max_tof:
						binParams->maxTOF = *((double *) paramBuffer);
						break;
					
					case PhgBinEn_min_pa:
						binParams->minPA = *((double *) paramBuffer);
						break;
					
					case PhgBinEn_max_pa:
						binParams->maxPA = *((double *) paramBuffer);
						break;
					
					case PhgBinEn_min_e:
						binParams->minE = *((double *) paramBuffer);
						break;

					case PhgBinEn_max_e:
						binParams->maxE= *((double *) paramBuffer);
						break;
						
					case PhgBinEn_min_s:
						binParams->minS = *((LbUsFourByte *) paramBuffer);
						break;

					case PhgBinEn_max_s:
						binParams->maxS = *((LbUsFourByte *) paramBuffer);
						break;

					case PhgBinEn_weight_image_type:
						binParams->weight_image_type  = 
							*((LbUsFourByte *) paramBuffer);
							
						/* Verify they didn't specify a non-supported type */
						if ((binParams->weight_image_type != PHG_BIN_WEIGHT_TYPE_R4) &&
							(binParams->weight_image_type != PHG_BIN_WEIGHT_TYPE_R8)) {

							sprintf(phgBinErrString, "You have specified an invalid type for your weight image in the binning parameters file (type = %ld)\n",
								(unsigned long)binParams->weight_image_type);
							ErStGeneric(phgBinErrString);
							switchOkay = false;
						}
						
						break;

					case PhgBinEn_count_image_type:
						binParams->count_image_type  = 
							*((LbUsFourByte *) paramBuffer);
						break;

					case PhgBinEn_sum_according_to_type:
						binParams->sumAccordingToType  = 
							*((Boolean *) paramBuffer);
						break;

					case PhgBinEn_add_to_existing_img:
						binParams->addToExistingImg = *((Boolean *) paramBuffer);
						break;

					case PhgBinEn_weight_image_path:
							strcpy(binParams->weightImgFilePath,
								(char *) paramBuffer);
							
							binParams->doWeights = 
								(binParams->weightImgFilePath[0] != '\0');
						break;

					case PhgBinEn_weight_squared_image_path:
							strcpy(binParams->weightSquImgFilePath,
								(char *) paramBuffer);
							
							binParams->doWeightsSquared = 
								(binParams->weightSquImgFilePath[0] != '\0');
						break;

					case PhgBinEn_count_image_path:
							strcpy(binParams->countImgFilePath,
								(char *) paramBuffer);
							
							binParams->doCounts = 
								(binParams->countImgFilePath[0] != '\0');
						break;
						
					case 	PhgBinEn_num_theta_bins:
					
						binParams->numThetaBins  = 
							*((LbUsFourByte *) paramBuffer);

						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_THETA;
						dimensionCount++;
						break;

					case 	PhgBinEn_max_theta:
						binParams->maxTheta  = 
							*((double *) paramBuffer);
						binParams->minTheta = -binParams->maxTheta;
						break;
						
					case 	PhgBinEn_num_phi_bins:
					
						binParams->numPHIBins  = 
							*((LbUsFourByte *) paramBuffer);

						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_PHI;
						dimensionCount++;
						break;

					case 	PhgBinEn_num_xr_bins:

						binParams->numXRBins  = 
							*((LbUsFourByte *) paramBuffer);

						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_XR;
						dimensionCount++;
						break;
						
					case 	PhgBinEn_num_yr_bins:
						binParams->numYRBins  = 
							*((LbUsFourByte *) paramBuffer);

						binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_YR;
						dimensionCount++;
						break;
						
					case 	PhgBinEn_min_xr:
						binParams->minXR  = 
							*((double *) paramBuffer);
						break;
						
					case 	PhgBinEn_max_xr:
						binParams->maxXR  = 
							*((double *) paramBuffer);
						break;
						
					case 	PhgBinEn_min_yr:
						binParams->minYR  = 
							*((double *) paramBuffer);
						break;
						
					case 	PhgBinEn_max_yr:
						binParams->maxYR  = 
							*((double *) paramBuffer);
						break;
						
					case 	PhgBinEn_bin_by_crystal:
						binParams->doCrystal  = 
							*((Boolean *) paramBuffer);
						if (binParams->doCrystal) {
						
							/* check to make sure that the user has specified a block detector tomograph
							 with at least one active crystal / else set number of crystals for binning */
							if (DetRunTimeParams[DetCurParams].BlockTomoDetector.NumActiveElementsInTomo <= 0) {
								if (DetRunTimeParams[DetCurParams].DetectorType == DetEn_Block) {
									/* report error */
									sprintf(phgBinErrString, "You have selected binning by crystal, but your detector has no active crystals. (PhgBinGetParams)\n");
									ErStGeneric(phgBinErrString);
									switchOkay = false;
								} else {
									/* report error */
									sprintf(phgBinErrString, "Binning by crystal is available only with block detectors. (PhgBinGetParams)\n");
									ErStGeneric(phgBinErrString);
									switchOkay = false;
								}
								
							} else {
								binParams->numCrystalBins = DetRunTimeParams[DetCurParams].BlockTomoDetector.NumActiveElementsInTomo;
							}
							
							binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_Crystal1;
							dimensionCount++;
							binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-dimensionCount-1] = PhgBinEn_Crystal2;
							dimensionCount++;
						}
						break;
						
					case 	PhgBinEn_do_ssrb:
						binParams->doSSRB  = 
							*((Boolean *) paramBuffer);
						break;
						
					case 	PhgBinEn_do_msrb:
						binParams->doMSRB  = 
							*((Boolean *) paramBuffer);
						break;
												
					case 	PhgBinEn_detector_radius:
						binParams->detector_radius  = 
							*((double *) paramBuffer);
						break;
					
					case PhgBinEn_history_file:
							strcpy(binParams->history_path,
								(char *) paramBuffer);
							
							binParams->isHistoryFile = 
								(binParams->history_path[0] != '\0');
						break;
					
					case PhgBinEn_history_params_file:
							strcpy(binParams->history_params_path,
								(char *) paramBuffer);
							
							binParams->isHistoryParamsFile = 
								(binParams->history_params_path[0] != '\0');
						break;

					case PhgBinEn_binPETasSPECT:
						binParams->isBinPETasSPECT  = 
							*((Boolean *) paramBuffer);
						
						/* the BinPETasSPECT option is incompatible with the SPECT
						and PETCoincidencesOnly simulations */
						if ( binParams->isBinPETasSPECT && 
							( PHG_IsPETCoincidencesOnly() || PHG_IsSPECT() ) ) {
							
							sprintf(phgBinErrString, 
								"binPETasSPECT = TRUE is not compatible with simulate_PET_coincidences_only or simulate_SPECT = TRUE.\n");
							ErStGeneric(phgBinErrString);
							switchOkay = false;
							
						}
						
						break;

					case PhgBinEn_NULL:
						sprintf(phgBinErrString, "(PhgGetBinParams) Unknown (hence unused) bin parameter (%s).\n",
                            paramLabel);
						ErAlert(phgBinErrString, false);
						break;
					
					default:
						sprintf(phgBinErrString, "(PhgGetBinParams) Unknown (hence unused) bin parameter (default) (%s).\n",
                            paramLabel);
						ErAlert(phgBinErrString, false);
						break;
				}
				if (!switchOkay)
					break;
		}

		/* Save dimension count */
		binParams->numDimensions = dimensionCount;
		
		/* Close the parameter file */
		LbPfClose(&phgParamFileHk);
		
		/* Verify that there are no duplicated names within the parameter files (for image files) */
		if (PhgNumBinParams > 1) {
			for (curBinParams = 0; curBinParams < PhgNumBinParams-1; curBinParams++) {
				if ((strcmp(PhgBinParams[curBinParams].countImgFilePath,
						PhgBinParams[curBinParams+1].countImgFilePath) == 0)
						&& (strlen(PhgBinParams[curBinParams].countImgFilePath) != 0)) {
					sprintf(phgBinErrString, "You have duplicated names for count images in binning parameter files:\n"
						"'%s' and '%s', the names are:\n'%s' and '%s'",
						PhgRunTimeParams.PhgBinParamsFilePath[curBinParams],
						PhgRunTimeParams.PhgBinParamsFilePath[curBinParams+1],
						PhgBinParams[curBinParams].countImgFilePath,
						PhgBinParams[curBinParams+1].countImgFilePath);
					ErStGeneric(phgBinErrString);
					goto FAIL;
				}

				if ((strcmp(PhgBinParams[curBinParams].weightImgFilePath,
						PhgBinParams[curBinParams+1].weightImgFilePath) == 0)
						&& (strlen(PhgBinParams[curBinParams].weightImgFilePath) != 0)) {
					sprintf(phgBinErrString, "You have duplicated names for weight images in binning parameter files:\n"
						"'%s' and '%s', the names are:\n'%s' and '%s'",
						PhgRunTimeParams.PhgBinParamsFilePath[curBinParams],
						PhgRunTimeParams.PhgBinParamsFilePath[curBinParams+1],
						PhgBinParams[curBinParams].weightImgFilePath,
						PhgBinParams[curBinParams+1].weightImgFilePath);
					ErStGeneric(phgBinErrString);
					goto FAIL;
				}

				if ((strcmp(PhgBinParams[curBinParams].weightSquImgFilePath,
						PhgBinParams[curBinParams+1].weightSquImgFilePath) == 0)
						&& (strlen(PhgBinParams[curBinParams].weightSquImgFilePath) != 0)) {
					sprintf(phgBinErrString, "You have duplicated names for weight squared images in binning parameter files:\n"
						"'%s' and '%s', the names are:\n'%s' and '%s'",
						PhgRunTimeParams.PhgBinParamsFilePath[curBinParams],
						PhgRunTimeParams.PhgBinParamsFilePath[curBinParams+1],
						PhgBinParams[curBinParams].weightSquImgFilePath,
						PhgBinParams[curBinParams+1].weightSquImgFilePath);
					ErStGeneric(phgBinErrString);
					goto FAIL;
				}

				if ((strcmp(PhgBinParams[curBinParams].history_path,
						PhgBinParams[curBinParams+1].history_path) == 0)
						&& (strlen(PhgBinParams[curBinParams].history_path) != 0)) {
					sprintf(phgBinErrString, "You have duplicated names for history files in binning parameter files:\n"
						"'%s' and '%s', the names are:\n'%s' and '%s'",
						PhgRunTimeParams.PhgBinParamsFilePath[curBinParams],
						PhgRunTimeParams.PhgBinParamsFilePath[curBinParams+1],
						PhgBinParams[curBinParams].history_path,
						PhgBinParams[curBinParams+1].history_path);
					ErStGeneric(phgBinErrString);
					goto FAIL;
				}
			}
		}
		
		/* See if we quit due to error */
		if (!isEOF)
			break;
		
		okay = true;
		FAIL:;
	} while (false);
	
	return (okay);
}
/*********************************************************************************
*
*			Name:			PhgLookupBinParamLabel
*
*			Summary:	Find the given label in the label table and return
*						its enumeration value.
*			Arguments:
*				char		*label	- The label to find.
*				
*			Function return: Enumerated type corresponding to label.
*
*********************************************************************************/
PhgEn_BinParamsTy  PhgLookupBinParamLabel(char *label)	
{
	PhgEn_BinParamsTy		whichParam;			/* Parameter we find */
	LbUsFourByte			labelIndex;			/* Label we are checking */

	/* Set param to null */
	whichParam = PhgBinEn_NULL;
	
	/* Search table for given label */
	for (labelIndex = 0; labelIndex < (PhgBinEn_NULL); labelIndex++){
	
		/* See if it matches */
		if (strcmp(label, phgBinParamLabels[labelIndex]) == 0) {
			whichParam = (PhgEn_BinParamsTy) labelIndex;
			break;
		}
	}
	
	return (whichParam);
}

#undef PHG_PARAMS
