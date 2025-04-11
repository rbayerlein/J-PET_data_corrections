/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1995-2005 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		PhgHdr.c
*			Revision Number:	1.2
*			Date last revised:	2007
*			Programmer:			Steven Vannoy
*			Date Originated:	14 April, 1995
*
*			Module Overview:	Manages the header file/data structure jobs.
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
*********************************************************************************
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
*						- binning by random-state and crystal number.
*						- binning PET data as SPECT
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		R Harrison
*
*			Revision date:		Jan - Feb 2005
*
*			Revision description:	Altered to read in/write out new Eight-byte
*						integer fields
*
*********************************************************************************/

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbMemory.h"
#include "LbInterface.h"
#include "LbParamFile.h"
#include "LbConvert.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhoHFile.h"
#include "PhgHdr.h"
#include "phg.h"
#include "PhgBin.h"

/* LOCAL CONSTANTS */

/* LOCAL TYPES */

/* LOCAL GLOBALS */

/* PROTOTYPES */

/* FUNCTIONS */

/**********************
*	PhgHdrGtHeaderHk
*
*	Purpose:	This routine takes a file pointer to a file known to
*				contain a PHG header and returns a header hook to it.
*
*	Arguments:
*			FILE			*headerFl	- The new file to associate with
*			LbHdrHkTy		*headerHk	- The newly created header structure.
*
*	Result:	TRUE unless an error occurs.
***********************/
Boolean  PhgHdrGtHeaderHk(FILE *headerFl, LbHdrHkTy *headerHk)
{
	return(LbHdrOpen(headerFl, PHG_HDR_HEADER_SIZE, headerHk));
}

/**********************
*	PhgHdrGtField
*
*	Purpose:	This routine retrieves data for a specific field in the
*				header.
*
*	Arguments:
*			LbHdrHkTy		*headerHk	- The newly created header structure.
*			LbUsFourByte	fieldID		- Which field to retrieve.
*			void			*fieldData	- Where to store the data.
*
*	Result:	TRUE unless an error occurs.
***********************/
Boolean  PhgHdrGtField(LbHdrHkTy *headerHk, LbUsFourByte fieldID, void *fieldData)
{
	Boolean 	okay = false;
	LbFourByte	fieldSize;
	
	
	do { /* Process Loop */
		
		/* Get the field size */
		if ((fieldSize = PhgHdrGtFieldSize(fieldID)) == -1)
			break;
			
		/* Get the field data */
		if (LbHdrGtElem(headerHk, fieldID, fieldSize, fieldData) == false) {
			
			break;
		}
		
		okay = true;
	} while (false);
	
	return (okay);	
}


/**********************
*	PhgHdrGtFieldSize
*
*	Purpose:	This routine returns the size of a field based on it's ID.
*
*	Arguments:
*			LbFourByte	fieldID	- The ID of the field in question.
*
*	Result:	The size of the field. NOTE: Returns -1 if field ID is invalid.
***********************/
LbFourByte PhgHdrGtFieldSize(LbUsFourByte fieldID)
{
	LbFourByte		fieldSize;
	PhoHFileHdrTy	params;
		
	switch (fieldID) {
	
		case HDR_PHG_PHGPARAM__ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_SUBOBJACTIVITYINDEX_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_SUBOBJACTIVITYTABLE_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_SUBOBJACTINDEXTRANS_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_SUBOBJACTIMG_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_SUBOBJATTENINDEX_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_SUBOBJATTENTABLE_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_SUBOBJATTINDEXTRANS_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_SUBOBJATTIMG_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_PRODTBLINPUTTABLE_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_PRODTBLOUTPUTTABLE_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_PHOHSTAT_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_PHOHFILEHISTORY_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_BINPARAMS_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_COLLIMATORPARAMS_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_PHG_DETECTORPARAMS_ID:
			fieldSize = PATH_LENGTH;
			break;
	
		case HDR_PHG_HEADER_SIZE_ID:
			fieldSize = sizeof(params.H.HdrSize);
			break;
			
		case HDR_PHG_HEADER_KIND_ID:
			fieldSize = sizeof(params.H.HdrKind);
			break;

		case HDR_PHG_HEADER_VERS_ID:
			fieldSize = sizeof(params.H.HdrVersion);
			break;

		case HDR_PHG_EVENTS_TO_SIM_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.Phg_EventsToSimulate);
			break;

		case HDR_PHG_IS_CALC_EVENTS_TO_SIM_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsCalcEventsToSimulate);
			break;

		case HDR_PHG_EVENTS_TO_SIM_OLD4BYTE_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.Phg_EventsToSimulateOld);
			break;

		case HDR_PHG_RAND_SEED_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgRandomSeed);
			break;

		case HDR_PHG_LEN_OF_SCAN_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.Phg_LengthOfScan);
			break;
					 
		case HDR_PHG_ACC_ANGLE_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.Phg_AcceptanceAngle);
			break;
			 
		case HDR_PHG_ACC_ANGLE_SINE_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.Phg_SineOfAcceptanceAngle);
			break;
			 
		case HDR_PHG_MIN_ENERGY_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgMinimumEnergy);
			break;
			 
		case HDR_PHG_MIN_WW_RATIO_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgMinWWRatio);
			break;
			 
		case HDR_PHG_MAX_WW_RATIO_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgMaxWWRatio);
			break;
		
		/* Due to change from two to four bytes internal representation
			we force this to be a two byte field
		*/
		case HDR_PHG_ISOTOPE_ID:
			fieldSize = sizeof(LbUsTwoByte);
			break;
			 
		case HDR_PHG_PHOTON_ENERGY_KEV_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgNuclide.photonEnergy_KEV);
			break;
			 
		case HDR_PHG_POSITRON_ENERGY_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgNuclide.maximumPositronEnergy);
			break;
			 
		case HDR_PHG_IS_FORCED_DETECTION_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsForcedDetection);
			break;
			 
		case HDR_PHG_IS_STRATIFICATION_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsStratification);
			break;
			 
		case HDR_PHG_IS_NON_ABSORPTION_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsForcedNonAbsorbtion);
			break;
			 
		case HDR_PHG_IS_SPECT_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsSPECT);
			break;
			 
		case HDR_PHG_IS_PET_COINC_ONLY_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsPETCoincidencesOnly);
			break;

		case HDR_PHG_IS_PET_COINC_PLUS_SINGLES_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsPETCoincPlusSingles);
			break;

		case HDR_PHG_IS_MULTIEMISSION_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsMultiEmission);
			break;

		case HDR_PHG_IS_PET_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsPET);
			break;

		case HDR_PHG_IS_PET_COINC_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsPETCoincidences);
			break;

		case HDR_PHG_IS_PET_SINGLES_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsPETSingles);
			break;

		case HDR_PHG_IS_HISTORY_FILE_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsHistoryFile);
			break;

		case HDR_PHG_IS_ADJ_FOR_POS_RANGE_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsAdjForPosRange);
			break;

		case HDR_PHG_IS_ADJ_FOR_COLIN_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsAdjForCollinearity);
			break;

		case HDR_PHG_IS_COMPUTED_PROD_TBL_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsComputedProductivityTbl);
			break;

		case HDR_PHG_IS_VOXEL_PTSRC_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsVoxelPointSource);
			break;

		case HDR_PHG_IS_BIN_ON_THE_FLY_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsBinOnTheFly);
			break;

		case HDR_PHG_IS_COL_ON_THE_FLY_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsCollimateOnTheFly);
			break;

		case HDR_PHG_IS_DET_ON_THE_FLY_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsDetectOnTheFly);
			break;

		case HDR_PHG_IS_VOXEL_LNSRC_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsVoxelLineSource);
			break;

		case HDR_PHG_IS_POLARIZATION_ID:
			fieldSize = sizeof(params.H.PhgRunTimeParams.PhgIsModelPolarization);
			break;

		case HDR_COL_IS_INITIALIZED_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.Initialized);
			break;

		case HDR_COL_HISTORY_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_COL_TYPE_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.ColType);
			break;

		case HDR_COL_DEPTH_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.SimplePETCol.Depth);
			break;

		case HDR_COL_UNC_HOLE_GEOM_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.UNCSPECTCol.HoleGeometry);
			break;

		case HDR_COL_UNC_RAD_OF_ROTATION_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.UNCSPECTCol.RadiusOfRotation);
			break;

		case HDR_COL_UNC_THICKNESS_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.UNCSPECTCol.Thickness);
			break;

		case HDR_COL_UNC_HOLE_RADIUS_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.UNCSPECTCol.HoleRadius);
			break;

		case HDR_COL_UNC_SEPTAL_THICKNESS_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.UNCSPECTCol.SeptalThickness);
			break;

		case HDR_COL_UNC_MIN_Z_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.UNCSPECTCol.MinZ);
			break;

		case HDR_COL_UNC_MAX_Z_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.UNCSPECTCol.MaxZ);
			break;

		case HDR_COL_UNC_START_ANGLE_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.UNCSPECTCol.StartAngle);
			break;

		case HDR_COL_UNC_STOP_ANGLE_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.UNCSPECTCol.StopAngle);
			break;

		case HDR_COL_UNC_SUM_ALL_VIEWS_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.UNCSPECTCol.SumAllViews);
			break;

		case HDR_COL_UNC_NUM_VIEWS_ID:
			fieldSize = sizeof(params.H.ColRunTimeParams.UNCSPECTCol.NumViews);
			break;

		case HDR_DET_IS_INITIALIZED_ID:
			fieldSize = sizeof(params.H.DetRunTimeParams.Initialized);
			break;

		case HDR_DET_IS_DO_HISTORY_ID:
			fieldSize = sizeof(params.H.DetRunTimeParams.DoHistory);
			break;

		case HDR_DET_HISTORY_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_DET_TYPE_ID:
			fieldSize = sizeof(params.H.DetRunTimeParams.DetectorType);
			break;

		case HDR_DET_ENERGY_RESOLUTION_PER_ID:
			fieldSize = sizeof(params.H.DetRunTimeParams.EnergyResolutionPercentage);
			break;

		case HDR_DET_REFERENCE_ENERGY_ID:
			fieldSize = sizeof(params.H.DetRunTimeParams.ReferenceEnergy);
			break;

		case HDR_DET_IS_FORCED_INTERACTION_ID:
			fieldSize = sizeof(params.H.DetRunTimeParams.DoForcedInteraction);
			break;

		case HDR_DET_PHOTON_TIME_FWHM_ID:
			fieldSize = sizeof(params.H.DetRunTimeParams.PhotonTimeFWHM);
			break;

		case HDR_DET_COINC_TIMING_WINDOW_ID:
			fieldSize = sizeof(params.H.DetRunTimeParams.CoincidenceTimingWindowNS);
			break;

		case HDR_DET_TRIPLES_METHOD_ID:
			fieldSize = sizeof(params.H.DetRunTimeParams.TriplesMethod);
			break;

		case HDR_DET_RANDOMS_HISTORY_FILE_ID:
			fieldSize = PATH_LENGTH;
			break;

		case HDR_CYL_TARGET_RADIUS_ID:
			fieldSize = sizeof(params.H.TargetCylinder.radius);
			break;

		case HDR_CYL_TARGET_MIN_Z_ID:
			fieldSize = sizeof(params.H.TargetCylinder.zMin);
			break;

		case HDR_CYL_TARGET_MAX_Z_ID:
			fieldSize = sizeof(params.H.TargetCylinder.zMax);
			break;

		case HDR_CYL_TARGET_CENTER_X_ID:
			fieldSize = sizeof(params.H.TargetCylinder.centerX);
			break;

		case HDR_CYL_TARGET_CENTER_Y_ID:
			fieldSize = sizeof(params.H.TargetCylinder.centerY);
			break;

		case HDR_CYL_CRITZONE_RADIUS_ID:
			fieldSize = sizeof(params.H.CriticalZone.radius);
			break;

		case HDR_CYL_CRITZONE_MIN_Z_ID:
			fieldSize = sizeof(params.H.CriticalZone.zMin);
			break;

		case HDR_CYL_CRITZONE_MAX_Z_ID:
			fieldSize = sizeof(params.H.CriticalZone.zMax);
			break;

		case HDR_CYL_CRITZONE_CENTER_X_ID:
			fieldSize = sizeof(params.H.CriticalZone.centerX);
			break;

		case HDR_CYL_CRITZONE_CENTER_Y_ID:
			fieldSize = sizeof(params.H.CriticalZone.centerY);
			break;
		case HDR_CYL_OBJCYL_RADIUS_ID:
			fieldSize = sizeof(params.H.ObjectCylinder.radius);
			break;

		case HDR_CYL_OBJCYL_MIN_Z_ID:
			fieldSize = sizeof(params.H.ObjectCylinder.zMin);
			break;

		case HDR_CYL_OBJCYL_MAX_Z_ID:
			fieldSize = sizeof(params.H.ObjectCylinder.zMax);
			break;

		case HDR_CYL_OBJCYL_CENTER_X_ID:
			fieldSize = sizeof(params.H.ObjectCylinder.centerX);
			break;

		case HDR_CYL_OBJCYL_CENTER_Y_ID:
			fieldSize = sizeof(params.H.ObjectCylinder.centerY);
			break;
		case HDR_CYL_LIMITCYL_RADIUS_ID:
			fieldSize = sizeof(params.H.LimitCylinder.radius);
			break;

		case HDR_CYL_LIMITCYL_MIN_Z_ID:
			fieldSize = sizeof(params.H.LimitCylinder.zMin);
			break;

		case HDR_CYL_LIMITCYL_MAX_Z_ID:
			fieldSize = sizeof(params.H.LimitCylinder.zMax);
			break;

		case HDR_CYL_LIMITCYL_CENTER_X_ID:
			fieldSize = sizeof(params.H.LimitCylinder.centerX);
			break;

		case HDR_CYL_LIMITCYL_CENTER_Y_ID:
			fieldSize = sizeof(params.H.LimitCylinder.centerY);
			break;

		case HDR_BIN1_NUM_SIMULATIONS_ID:
			fieldSize = sizeof(params.H.NumSimulations);
			break;

		case HDR_BIN1_EVENTS_TO_SIMULATION_ID:
			fieldSize = sizeof(params.H.SumEventsToSimulate);
			break;

		case HDR_BIN1_NUM_PHOTONS_ID:
			fieldSize = sizeof(params.H.NumPhotons);
			break;

		case HDR_BIN1_NUM_DECAYS_ID:
			fieldSize = sizeof(params.H.NumDecays);
			break;

		case HDR_BIN_NUM_Z_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numZBins);
			break;

		case HDR_BIN_NUM_PA_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numPABins);
			break;

		case HDR_BIN_NUM_TD_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numTDBins);
			break;

		case HDR_BIN_NUM_AA_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numAABins);
			break;

		case HDR_BIN_NUM_TOF_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numTOFBins);
			break;

		case HDR_BIN_NUM_E1_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numE1Bins);
			break;

		case HDR_BIN_NUM_E2_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numE2Bins);
			break;

		case HDR_BIN_NUM_S1_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numS1Bins);
			break;

		case HDR_BIN_NUM_S2_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numS2Bins);
			break;

		case HDR_BIN_NUM_IMAGE_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numImageBins);
			break;

		case HDR_BIN_SCATTER_PARAM_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.scatterRandomParam);
			break;

		case HDR_BIN_MIN_Z_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.minZ);
			break;
			 

		case HDR_BIN_MAX_Z_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.maxZ);
			break;

		case HDR_BIN_MIN_PA_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.minPA);
			break;

		case HDR_BIN_MAX_PA_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.maxPA);
			break;

		case HDR_BIN_MIN_TD_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.minTD);
			break;

		case HDR_BIN_MAX_TD_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.maxTD);
			break;

		case HDR_BIN_MIN_AA_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.minAA);
			break;

		case HDR_BIN_MAX_AA_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.maxAA);
			break;

		case HDR_BIN_MIN_TOF_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.minTOF);
			break;

		case HDR_BIN_MAX_TOF_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.maxTOF);
			break;

		case HDR_BIN_MIN_ENERGY_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.minE);
			break;

		case HDR_BIN_MAX_ENERGY_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.maxE);
			break;

		case HDR_BIN_MIN_SCAT_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.minS);
			break;

		case HDR_BIN_MAX_SCAT_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.maxS);
			break;

		case HDR_BIN_ADD_TO_EXISTING_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.addToExistingImg);
			break;

		case HDR_BIN_DO_COUNTS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.doCounts);
			break;

		case HDR_BIN_DO_WEIGHTS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.doWeights);
			break;

		case HDR_BIN_DO_WEIGHTS_SQU_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.doWeightsSquared);
			break;

		case HDR_BIN_Z_RANGE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.zRange);
			break;

		case HDR_BIN_E_RANGE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.eRange);
			break;

		case HDR_BIN_S_RANGE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.sRange);
			break;

		case HDR_BIN_TD_RANGE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.tdRange);
			break;

		case HDR_BIN_SCATTER2CI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.scatter2CIsize);
			break;

		case HDR_BIN_SCATTER2WI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.scatter2WIsize);
			break;

		case HDR_BIN_SCATTER2WIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.scatter2WISsize);
			break;

		case HDR_BIN_SCATTER1CI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.scatter1CIsize);
			break;

		case HDR_BIN_SCATTER1WI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.scatter1WIsize);
			break;

		case HDR_BIN_SCATTER1WIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.scatter1WISsize);
			break;

		case HDR_BIN_ENERGY2CI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.energy2CIsize);
			break;

		case HDR_BIN_ENERGY2WI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.energy2WIsize);
			break;

		case HDR_BIN_ENERGY2WIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.energy2WISsize);
			break;

		case HDR_BIN_ENERGY1CI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.energy1CIsize);
			break;

		case HDR_BIN_ENERGY1WI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.energy1WIsize);
			break;

		case HDR_BIN_ENERGY1WIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.energy1WISsize);
			break;

		case HDR_BIN_AACI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.aaCIsize);
			break;

		case HDR_BIN_AAWI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.aaWIsize);
			break;

		case HDR_BIN_AAWIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.aaWISsize);
			break;

		case HDR_BIN_TOFCI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.tofCIsize);
			break;

		case HDR_BIN_TOFWI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.tofWIsize);
			break;

		case HDR_BIN_TOFWIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.tofWISsize);
			break;

		case HDR_BIN_TDCI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.tdCIsize);
			break;

		case HDR_BIN_TDWI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.tdWIsize);
			break;

		case HDR_BIN_TDWIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.tdWISsize);
			break;

		case HDR_BIN_PACI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.paCIsize);
			break;

		case HDR_BIN_PAWI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.paWIsize);
			break;

		case HDR_BIN_PAWIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.paWISsize);
			break;

		case HDR_BIN_Z2CI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.z2CIsize);
			break;

		case HDR_BIN_Z2WI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.z2WIsize);
			break;

		case HDR_BIN_Z2WIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.z2WISsize);
			break;

		case HDR_BIN_Z1CI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.z1CIsize);
			break;

		case HDR_BIN_Z1WI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.z1WIsize);
			break;

		case HDR_BIN_Z1WIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.z1WISsize);
			break;

		case HDR_BIN_COUNT_IMG_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.countImageSize);
			break;

		case HDR_BIN_WEIGHT_IMG_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.weightImageSize);
			break;

		case HDR_BIN_WEIGHT_SQU_IMG_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.weightSquImageSize);
			break;

		case HDR_BIN_WEIGHT_IMG_TYPE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.weight_image_type);
			break;

		case HDR_BIN_COUNT_IMG_TYPE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.count_image_type);
			break;

		case HDR_BIN_DIMENSION_ORDER_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.PhgBinDimensions);
			break;

		case HDR_BIN_PHICI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.phiCIsize);
			break;

		case HDR_BIN_PHIWI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.phiWIsize);
			break;

		case HDR_BIN_PHIWIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.phiWISsize);
			break;

		case HDR_BIN_THETACI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.thetaCIsize);
			break;

		case HDR_BIN_THETAWI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.thetaWIsize);
			break;

		case HDR_BIN_THETAWIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.thetaWISsize);
			break;

		case HDR_BIN_XRCI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.xrCIsize);
			break;

		case HDR_BIN_XRWI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.xrWIsize);
			break;

		case HDR_BIN_XRWIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.xrWISsize);
			break;

		case HDR_BIN_YRCI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.yrCIsize);
			break;

		case HDR_BIN_YRWI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.yrWIsize);
			break;

		case HDR_BIN_YRWIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.yrWISsize);
			break;

		case HDR_BIN_NUM_PHI_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numPHIBins);
			break;

		case HDR_BIN_NUM_THETA_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numThetaBins);
			break;

		case HDR_BIN_NUM_XR_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numXRBins);
			break;

		case HDR_BIN_NUM_YR_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numYRBins);
			break;

		case HDR_BIN_MIN_THETA_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.minTheta);
			break;

		case HDR_BIN_MAX_THETA_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.maxTheta);
			break;

		case HDR_BIN_MIN_PHI_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.minPhi);
			break;

		case HDR_BIN_MAX_PHI_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.maxPhi);
			break;

		case HDR_BIN_MIN_XR_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.minXR);
			break;

		case HDR_BIN_MAX_XR_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.maxXR);
			break;

		case HDR_BIN_MIN_YR_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.minYR);
			break;

		case HDR_BIN_MAX_YR_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.maxYR);
			break;

		case HDR_BIN_DO_CRYSTAL_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.numCrystalBins);
			break;

		case HDR_BIN_NUM_CRYSTAL_BINS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.doCrystal);
			break;

		case HDR_BIN_CRYSTAL1CI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.crystal1CIsize);
			break;

		case HDR_BIN_CRYSTAL1WI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.crystal1WIsize);
			break;

		case HDR_BIN_CRYSTAL1WIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.crystal1WISsize);
			break;

		case HDR_BIN_CRYSTAL2CI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.crystal2CIsize);
			break;

		case HDR_BIN_CRYSTAL2WI_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.crystal2WIsize);
			break;

		case HDR_BIN_CRYSTAL2WIS_SIZE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.crystal2WISsize);
			break;

		case HDR_BIN_DOSSRB_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.doSSRB);
			break;

		case HDR_BIN_DOMSRB_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.doMSRB);
			break;

		case HDR_BIN_IMG_RADIUS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.image_radius);
			break;

		case HDR_BIN_DET_RADIUS_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.detector_radius);
			break;

		case HDR_BIN_WTIMG_PATH_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.weightImgFilePath);
			break;

		case HDR_BIN_WTSIMG_PATH_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.weightSquImgFilePath);
			break;

		case HDR_BIN_CTIMG_PATH_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.countImgFilePath);
			break;

		case HDR_BIN_THETA_RANGE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.thetaRange);
			break;

		case HDR_BIN_PHI_RANGE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.phiRange);
			break;

		case HDR_BIN_XR_RANGE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.xrRange);
			break;

		case HDR_BIN_YR_RANGE_ID:
			fieldSize = sizeof(params.H.BinRunTimeParams.yrRange);
			break;

		case HDR_BIN_SUM_WEIGHTS_ID:
			fieldSize = sizeof(params.H.weightSum);
			break;

		case HDR_BIN_SUM_WEIGHTS_SQ_ID:
			fieldSize = sizeof(params.H.weightSquSum);
			break;
		
		case HDR_HISTORY_FILE_IS_SORTED_ID:
			fieldSize = sizeof(params.H.isTimeSorted);
			break;

		case HDR_HISTORY_FILE_IS_RANDOMS_ADDED_ID:
			fieldSize = sizeof(params.H.isRandomsAdded);
			break;

		default:
			fieldSize = -1;
			break;
	}
	
	return (fieldSize);
}

/**********************
*	PhgHdrGtParams
*
*	Purpose:	Get our runtime parameters
*				from the header stored in the given file.
*	Arguments:
*			FILE			*headerFile	- The file the header will be attatched to.
*			PhoHFileHdrTy	*paramsPtr	- The run time parameters structure.
*			LbHdrHkTy		*headerHk	- The newly created header structure.
*
*	Result:	True unless an error occurs.
***********************/
Boolean PhgHdrGtParams(FILE *headerFile, PhoHFileHdrTy *paramsPtr, LbHdrHkTy *headerHk)
{
	Boolean				okay = false;			/* Process Loop */
	Boolean				cleared;
	LbUsTwoByte			isotopeConverter;		/* For backwards compatibility */
	Boolean				saveInputError;			/* save error status when checking for
												old four-byte values */
	LbUsFourByte		oldFourByteValue;		/* For reading in header values that
												used to be 4 bytes (and are now *) from
												older files */
	
	do { /* Process Loop */
			
		/* Open a new header hook */
		if (LbHdrOpen(headerFile, PHG_HDR_HEADER_SIZE, headerHk) == false){
			PhgAbort("Unable to open header", false);
			break;
		}
		
		/* One-by-one, get the header fields */

		/* Get the header size */
		if (LbHdrGtElem(headerHk, HDR_PHG_HEADER_SIZE_ID, sizeof(paramsPtr->H.HdrSize),
				(void *)&(paramsPtr->H.HdrSize)) != true) {
			
			PhgAbort("Header is missing 'size' parameter", false);
			break;
		}

		/* Get the header kind */
		if (LbHdrGtElem(headerHk, HDR_PHG_HEADER_KIND_ID, sizeof(paramsPtr->H.HdrKind),
				(void *)&(paramsPtr->H.HdrKind)) != true) {
			
			PhgAbort("Header is missing 'kind' parameter", false);
			break;
		}

		/* Get the header version */
		if (LbHdrGtElem(headerHk, HDR_PHG_HEADER_VERS_ID, sizeof(paramsPtr->H.HdrVersion),
				(void *)&(paramsPtr->H.HdrVersion)) != true) {
			
			PhgAbort("Header is missing 'version' parameter", false);
			break;
		}
		
		/* save error status before fetching values that can be either 8 or 4 byte */
		saveInputError = ErIsInError();
		
		/* get the number of events to simulate */
		if (LbHdrGtElem(headerHk, HDR_PHG_EVENTS_TO_SIM_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.Phg_EventsToSimulate),
				(void *)&paramsPtr->H.PhgRunTimeParams.Phg_EventsToSimulate) == false){
			
			/* if the above call set ErIsInError(), clear it */
			if ( ErIsInError() && (!saveInputError) ) {
			
				ErClear();
			
			}
			
			/* see if this is an old header with a four-byte number of events to simulate */	
			if (LbHdrGtElem(headerHk, HDR_PHG_EVENTS_TO_SIM_OLD4BYTE_ID,
					sizeof(paramsPtr->H.PhgRunTimeParams.Phg_EventsToSimulateOld),
					(void *)&paramsPtr->H.PhgRunTimeParams.Phg_EventsToSimulateOld) == false){
						
				PhgAbort("Header is missing 'events to simulate' parameter", false);
				break;
				
			} else {
				
				/* assign 4-byte number of events to new number of events field */
				paramsPtr->H.PhgRunTimeParams.Phg_EventsToSimulate =
							paramsPtr->H.PhgRunTimeParams.Phg_EventsToSimulateOld;
				
			}
			
		}

		if (LbHdrGtElem(headerHk, HDR_PHG_RAND_SEED_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgRandomSeed),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgRandomSeed) == false){
			
			PhgAbort("Header is missing 'random seed' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_LEN_OF_SCAN_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.Phg_LengthOfScan),
						 (void *)&paramsPtr->H.PhgRunTimeParams.Phg_LengthOfScan) == false){
			
			PhgAbort("Header is missing 'length of scan' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_ACC_ANGLE_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.Phg_AcceptanceAngle),
			 (void *)&paramsPtr->H.PhgRunTimeParams.Phg_AcceptanceAngle) == false){
			
			PhgAbort("Header is missing 'acceptance angle' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_ACC_ANGLE_SINE_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.Phg_SineOfAcceptanceAngle),
			 (void *)&paramsPtr->H.PhgRunTimeParams.Phg_SineOfAcceptanceAngle) == false){
			
			PhgAbort("Header is missing 'sine of acceptance angle' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_MIN_ENERGY_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgMinimumEnergy),
			 (void *)&paramsPtr->H.PhgRunTimeParams.PhgMinimumEnergy) == false){
			
			PhgAbort("Header is missing 'minimum energy' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_MIN_WW_RATIO_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgMinWWRatio),
			 (void *)&paramsPtr->H.PhgRunTimeParams.PhgMinWWRatio) == false){
			
			PhgAbort("Header is missing 'minimum weight window ratio' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_MAX_WW_RATIO_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgMaxWWRatio),
			 (void *)&paramsPtr->H.PhgRunTimeParams.PhgMaxWWRatio) == false){
			
			PhgAbort("Header is missing 'maximum weight window ratio' parameter", false);
			break;
		}

		
		/* Retreive as two byte and convert to enumerate */
		if (LbHdrGtElem(headerHk, HDR_PHG_ISOTOPE_ID,
				sizeof(isotopeConverter),
			 (void *)&isotopeConverter) == false){
			
			/* Set to null */
			paramsPtr->H.PhgRunTimeParams.PhgNuclide.isotope = PhgEn_IsotopType_NULL;
			
			/* Clear the error */
			ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);

			ErAlert("Warning: Header is missing 'nuclide isotope' parameter", false);
			break;
		}
		paramsPtr->H.PhgRunTimeParams.PhgNuclide.isotope = (PhgEn_IsotopeTypeTy) isotopeConverter;
		
		
		if (LbHdrGtElem(headerHk, HDR_PHG_PHOTON_ENERGY_KEV_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgNuclide.photonEnergy_KEV),
			 (void *)&paramsPtr->H.PhgRunTimeParams.PhgNuclide.photonEnergy_KEV) == false){
			
			PhgAbort("Header is missing 'photon energy' parameter", false);
			break;
		}
		
			
		if (LbHdrGtElem(headerHk, HDR_PHG_POSITRON_ENERGY_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgNuclide.maximumPositronEnergy),
			 (void *)&paramsPtr->H.PhgRunTimeParams.PhgNuclide.maximumPositronEnergy) == false){
			
			PhgAbort("Header is missing 'maximum positron energy' parameter", false);
			break;
		}

			
		
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_FORCED_DETECTION_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsForcedDetection),
			 (void *)&paramsPtr->H.PhgRunTimeParams.PhgIsForcedDetection) == false){
			
			PhgAbort("Header is missing 'forced detection on/off' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_STRATIFICATION_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsStratification),
			 (void *)&paramsPtr->H.PhgRunTimeParams.PhgIsStratification) == false){
			
			PhgAbort("Header is missing 'stratification on/off' parameter", false);
			break;
		}

		
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_NON_ABSORPTION_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsForcedNonAbsorbtion),
			 (void *)&paramsPtr->H.PhgRunTimeParams.PhgIsForcedNonAbsorbtion) == false){
			
			PhgAbort("Header is missing 'forced non-absorption on/off' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_SPECT_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsSPECT),
			 (void *)&paramsPtr->H.PhgRunTimeParams.PhgIsSPECT) == false){
			
			PhgAbort("Header is missing 'SPECT on/off' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_PET_COINC_ONLY_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincidencesOnly),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincidencesOnly) == false){
			
			PhgAbort("Header is missing 'PET coincidences only on/off' parameter", false);
			break;
		}
		
			
		/* 	old headers will not have the following fields, set them to false in this case.
			save error status and restore when this happens */
		saveInputError = ErIsInError();
		
		{
			
			if (LbHdrGtElem(headerHk, HDR_PHG_IS_PET_COINC_PLUS_SINGLES_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincPlusSingles),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincPlusSingles) == false){
				
				/* if the above call set ErIsInError(), clear it */
				if ( ErIsInError() && (!saveInputError) ) {
				
					ErClear();
				
				}
				
				/* old header - set to false */
				paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincPlusSingles = false;
				
			}
		
			if (LbHdrGtElem(headerHk, HDR_PHG_IS_MULTIEMISSION_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsMultiEmission),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsMultiEmission) == false){
				
				/* if the above call set ErIsInError(), clear it */
				if ( ErIsInError() && (!saveInputError) ) {
				
					ErClear();
				
				}
				
				/* old header - set to false */
				paramsPtr->H.PhgRunTimeParams.PhgIsMultiEmission = false;
				
			}
		
		}
		
		
		/*	the following fields are computed from the previous four */
		if ( paramsPtr->H.PhgRunTimeParams.PhgIsSPECT ) {
		
			paramsPtr->H.PhgRunTimeParams.PhgIsPET = false;
			paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincidences = false;
			paramsPtr->H.PhgRunTimeParams.PhgIsPETSingles = false;
			
		} else if ( paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincidencesOnly ) {
		
			paramsPtr->H.PhgRunTimeParams.PhgIsPET = true;
			paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincidences = true;
			paramsPtr->H.PhgRunTimeParams.PhgIsPETSingles = false;
			
		} else if ( paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincPlusSingles ) {
		
			paramsPtr->H.PhgRunTimeParams.PhgIsPET = true;
			paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincidences = true;
			paramsPtr->H.PhgRunTimeParams.PhgIsPETSingles = true;
			
		} else {
		
			PhgAbort("Header missing valid simulation-type (SPECT, or PET coinc or coinc/singles)", false);
			break;
		
		}
		
			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_HISTORY_FILE_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsHistoryFile),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsHistoryFile) == false){
			
			PhgAbort("Header is missing 'history file on/off' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_ADJ_FOR_POS_RANGE_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsAdjForPosRange),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsAdjForPosRange) == false){
			
			PhgAbort("Header is missing 'adjust for positron range' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_ADJ_FOR_COLIN_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsAdjForCollinearity),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsAdjForCollinearity) == false){
			
			PhgAbort("Header is missing 'adjust for collinearity' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_COMPUTED_PROD_TBL_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsComputedProductivityTbl),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsComputedProductivityTbl) == false){
			
			PhgAbort("Header is missing 'computed productivity' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_VOXEL_PTSRC_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsVoxelPointSource),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsVoxelPointSource) == false){
			
			PhgAbort("Header is missing 'point source voxels' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_VOXEL_LNSRC_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsVoxelLineSource),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsVoxelLineSource) == false){
			
			/* Set to false */
			paramsPtr->H.PhgRunTimeParams.PhgIsVoxelLineSource = false;

			/* Clear the error */
			ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_POLARIZATION_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsModelPolarization),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsModelPolarization) == false){
			
			/* Set to false */
			paramsPtr->H.PhgRunTimeParams.PhgIsModelPolarization = false;

			/* Clear the error */
			ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_BIN_ON_THE_FLY_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsBinOnTheFly),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsBinOnTheFly) == false){
			
			PhgAbort("Header is missing 'bin on the fly' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_COL_ON_THE_FLY_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsCollimateOnTheFly),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsCollimateOnTheFly) == false){
			
			PhgAbort("Header is missing 'collimate on the fly' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_PHG_IS_DET_ON_THE_FLY_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsDetectOnTheFly),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsDetectOnTheFly) == false){
			
			PhgAbort("Header is missing 'detect on the fly' parameter", false);
			break;
		}

					
		if (LbHdrGtElem(headerHk, HDR_PHG_PHGPARAM__ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgParamFilePath) == false){
			
			PhgAbort("Header is missing 'param file path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJACTIVITYINDEX_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjActivityIndexFilePath) == false){
			
			PhgAbort("Header is missing 'activity index path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJACTIVITYTABLE_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjActivityTableFilePath) == false){
			
			PhgAbort("Header is missing 'activity table path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJACTINDEXTRANS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjActIndexTransFilePath) == false){
			
			PhgAbort("Header is missing 'activity index translation path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJACTIMG_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjActImgFilePath) == false){
			
			PhgAbort("Header is missing 'actvity image path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJATTENINDEX_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjAttenIndexFilePath) == false){
			
			PhgAbort("Header is missing 'attenuation idex path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJATTENTABLE_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjAttenTableFilePath) == false){
			
			PhgAbort("Header is missing 'attenuation table path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJATTINDEXTRANS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjAttIndexTransFilePath) == false){
			
			PhgAbort("Header is missing 'attenuation index translation path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJATTIMG_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjAttImgFilePath) == false){
			
			PhgAbort("Header is missing 'attenuation image path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_PRODTBLINPUTTABLE_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgProdTblInputTableFilePath) == false){
			
			PhgAbort("Header is missing 'productivity input path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_PRODTBLOUTPUTTABLE_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgProdTblOutputTableFilePath) == false){
			
			PhgAbort("Header is missing 'productivity output path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_PHOHSTAT_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgPhoHStatFilePath) == false){
			
			PhgAbort("Header is missing 'statistics file path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_PHOHFILEHISTORY_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgPhoHFileHistoryFilePath) == false){
			
			PhgAbort("Header is missing 'history file path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_BINPARAMS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgBinParamsFilePath) == false){
			
			PhgAbort("Header is missing 'bin parameters file path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_COLLIMATORPARAMS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgCollimatorParamsFilePath) == false){
			
			PhgAbort("Header is missing 'collimator parameters file path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_DETECTORPARAMS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgDetectorParamsFilePath) == false){
			
			PhgAbort("Header is missing 'detector parameters file path' parameter", false);
			break;
		}

					
		if (LbHdrGtElem(headerHk, HDR_PHG_PHGPARAM__ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgParamFilePath) == false){
			
			PhgAbort("Header is missing 'param file path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJACTIVITYINDEX_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjActivityIndexFilePath) == false){
			
			PhgAbort("Header is missing 'activity index path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJACTIVITYTABLE_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjActivityTableFilePath) == false){
			
			PhgAbort("Header is missing 'activity table path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJACTINDEXTRANS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjActIndexTransFilePath) == false){
			
			PhgAbort("Header is missing 'activity index translation path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJACTIMG_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjActImgFilePath) == false){
			
			PhgAbort("Header is missing 'actvity image path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJATTENINDEX_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjAttenIndexFilePath) == false){
			
			PhgAbort("Header is missing 'attenuation idex path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJATTENTABLE_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjAttenTableFilePath) == false){
			
			PhgAbort("Header is missing 'attenuation table path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJATTINDEXTRANS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjAttIndexTransFilePath) == false){
			
			PhgAbort("Header is missing 'attenuation index translation path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_SUBOBJATTIMG_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjAttImgFilePath) == false){
			
			PhgAbort("Header is missing 'attenuation image path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_PRODTBLINPUTTABLE_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgProdTblInputTableFilePath) == false){
			
			PhgAbort("Header is missing 'productivity input path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_PRODTBLOUTPUTTABLE_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgProdTblOutputTableFilePath) == false){
			
			PhgAbort("Header is missing 'productivity output path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_PHOHSTAT_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgPhoHStatFilePath) == false){
			
			PhgAbort("Header is missing 'statistics file path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_PHOHFILEHISTORY_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgPhoHFileHistoryFilePath) == false){
			
			PhgAbort("Header is missing 'history file path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_BINPARAMS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgBinParamsFilePath) == false){
			
			PhgAbort("Header is missing 'bin parameters file path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_COLLIMATORPARAMS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgCollimatorParamsFilePath) == false){
			
			PhgAbort("Header is missing 'collimator parameters file path' parameter", false);
			break;
		}
					
		if (LbHdrGtElem(headerHk, HDR_PHG_DETECTORPARAMS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgDetectorParamsFilePath) == false){
			
			PhgAbort("Header is missing 'detector parameters file path' parameter", false);
			break;
		}
		
		/* The following were  added after early SimSET versions, so default parameters
		 are given where possible rather than aborting (for backward compatibility). */
		{
			/* get flag that tells whether the number of events to simulated was calculated or given */	
			if (LbHdrGtElem(headerHk, HDR_PHG_IS_CALC_EVENTS_TO_SIM_ID,
					sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsCalcEventsToSimulate),
					(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsCalcEventsToSimulate) == false){
				
					paramsPtr->H.PhgRunTimeParams.PhgIsCalcEventsToSimulate = false;	/* Most likely */

					/* Clear the error */
					ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
		}


		/* Read collimator params if we are doing collimation, or this is a collimation history file */
		do { /* Collimation parameters */
		if ((paramsPtr->H.PhgRunTimeParams.PhgIsCollimateOnTheFly == true) ||
				(paramsPtr->H.HdrKind == PhoHFileEn_COL)){

			if (LbHdrGtElem(headerHk, HDR_COL_IS_INITIALIZED_ID,
					sizeof(paramsPtr->H.ColRunTimeParams.Initialized),
					(void *)&paramsPtr->H.ColRunTimeParams.Initialized) == false){
				
				/* If this isn't there, and we are a PHG history file, it is okay but
					we won't get any more collimator parameters 
				*/
				if (paramsPtr->H.HdrKind == PhoHFileEn_PHG) {
					paramsPtr->H.ColRunTimeParams.Initialized = false;
					break;
				}
				else {
					PhgAbort("Header is missing 'collimator initialized flag' parameter", false);
					goto FAIL;
				}
			}


			if (LbHdrGtElem(headerHk, HDR_COL_HISTORY_ID,
				PATH_LENGTH,
				(void *)&paramsPtr->H.ColRunTimeParams.ColHistoryFilePath) == false){
				
				PhgAbort("Header is missing 'collimator params file path' parameter", false);
				goto FAIL;
			}
				
			if (LbHdrGtElem(headerHk, HDR_COL_TYPE_ID,
				sizeof(paramsPtr->H.ColRunTimeParams.ColType),
				(void *)&paramsPtr->H.ColRunTimeParams.ColType) == false){
				
				PhgAbort("Header is missing 'collimator type' parameter", false);
				goto FAIL;
			}
	
			
			/* Determine which type of collimator was used and convert its fields
			*/
			switch (paramsPtr->H.ColRunTimeParams.ColType) {
			
				case ColEn_simple_pet:
					if (LbHdrGtElem(headerHk, HDR_COL_DEPTH_ID,
						sizeof(paramsPtr->H.ColRunTimeParams.SimplePETCol.Depth),
						(void *)&paramsPtr->H.ColRunTimeParams.SimplePETCol.Depth) == false){
				
						PhgAbort("Header is missing 'simple collimator depth' parameter", false);
						goto FAIL;
					}
	
				
					break;
					
				case ColEn_monte_carlo_pet:
					
					/* Just note that the actual collimator description is not stored in the header. Hence,
						we don't have anything to do here.
					*/
					paramsPtr->H.ColRunTimeParams.MCPETCol = 0;
					break;
					
				case ColEn_simple_spect:
					break;
					
				case ColEn_unc_spect:
					
					if (LbHdrGtElem(headerHk, HDR_COL_UNC_HOLE_GEOM_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.HoleGeometry),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.HoleGeometry) == false){
				
						PhgAbort("Header is missing 'UNC collimator hole geometry' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrGtElem(headerHk, HDR_COL_UNC_RAD_OF_ROTATION_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.RadiusOfRotation),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.RadiusOfRotation) == false){
				
						PhgAbort("Header is missing 'UNC collimator radius of rotation' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrGtElem(headerHk, HDR_COL_UNC_THICKNESS_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.Thickness),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.Thickness) == false){
						
						PhgAbort("Header is missing 'UNC collimator thickness' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrGtElem(headerHk, HDR_COL_UNC_HOLE_RADIUS_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.HoleRadius),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.HoleRadius) == false){
						
						PhgAbort("Header is missing 'UNC collimator hole radius' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrGtElem(headerHk, HDR_COL_UNC_SEPTAL_THICKNESS_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.SeptalThickness),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.SeptalThickness) == false){
						
						PhgAbort("Header is missing 'UNC collimator septal thickness' parameter", false);
						goto FAIL;
					}
			
						
					if (LbHdrGtElem(headerHk, HDR_COL_UNC_MIN_Z_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.MinZ),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.MinZ) == false){
						
						PhgAbort("Header is missing 'UNC collimator minimum Z' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrGtElem(headerHk, HDR_COL_UNC_MAX_Z_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.MaxZ),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.MaxZ) == false){
						
						PhgAbort("Header is missing 'UNC collimator maximum Z' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrGtElem(headerHk, HDR_COL_UNC_START_ANGLE_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.StartAngle),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.StartAngle) == false){
						
						PhgAbort("Header is missing 'UNC collimator start angle' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrGtElem(headerHk, HDR_COL_UNC_STOP_ANGLE_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.StopAngle),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.StopAngle) == false){
				
						PhgAbort("Header is missing 'UNC collimator stop angle' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrGtElem(headerHk, HDR_COL_UNC_SUM_ALL_VIEWS_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.SumAllViews),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.SumAllViews) == false){
						
						PhgAbort("Header is missing 'UNC collimator sum all views' parameter", false);
						goto FAIL;
					}
	
				
					if (LbHdrGtElem(headerHk, HDR_COL_UNC_NUM_VIEWS_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.NumViews),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.NumViews) == false){
						
						PhgAbort("Header is missing 'UNC collimator number of views' parameter", false);
						goto FAIL;
					}
	
						
				default:
					break;
			}
		}
		} while (false);
		
		/* Read detector params if we are doing detection */
		do { /* Detector Parameters */
		if ((paramsPtr->H.PhgRunTimeParams.PhgIsDetectOnTheFly == true) ||
			(paramsPtr->H.HdrKind == PhoHFileEn_DET)) {

			if (LbHdrGtElem(headerHk, HDR_DET_IS_INITIALIZED_ID,
					sizeof(paramsPtr->H.DetRunTimeParams.Initialized),
					(void *)&paramsPtr->H.DetRunTimeParams.Initialized) == false){

				/* If this isn't there, and we are a PHG history file, it is okay but
					we won't get any more collimator parameters 
				*/
				if (paramsPtr->H.HdrKind == PhoHFileEn_PHG) {
					paramsPtr->H.DetRunTimeParams.Initialized = false;
					break;
				}
				else {
					PhgAbort("Header is missing 'detector initialized' parameter", false);
					goto FAIL;
				}
			}
			
			if (LbHdrGtElem(headerHk, HDR_DET_IS_DO_HISTORY_ID,
					sizeof(paramsPtr->H.DetRunTimeParams.DoHistory),
					(void *)&paramsPtr->H.DetRunTimeParams.DoHistory) == false){
				
				PhgAbort("Header is missing 'detector history on/off' parameter", false);
				goto FAIL;
			}
	
						
			if (LbHdrGtElem(headerHk, HDR_DET_HISTORY_ID,
				PATH_LENGTH,
				(void *)&paramsPtr->H.DetRunTimeParams.DetHistoryFilePath) == false){
				
				PhgAbort("Header is missing 'detector params file path' parameter", false);
				goto FAIL;
			}
				
			if (LbHdrGtElem(headerHk, HDR_DET_TYPE_ID,
				sizeof(paramsPtr->H.DetRunTimeParams.DetectorType),
				(void *)&paramsPtr->H.DetRunTimeParams.DetectorType) == false){
				
				PhgAbort("Header is missing 'detector type' parameter", false);
				goto FAIL;
			}
	
			
			/* Determine which type of detector was used and convert its fields
			*/
			switch (paramsPtr->H.DetRunTimeParams.DetectorType) {
			
				case DetEn_simple_pet:
					if (LbHdrGtElem(headerHk, HDR_DET_ENERGY_RESOLUTION_PER_ID,
							sizeof(paramsPtr->H.DetRunTimeParams.EnergyResolutionPercentage),
							(void *)&paramsPtr->H.DetRunTimeParams.EnergyResolutionPercentage) == false){
					
						PhgAbort("Header is missing 'simple detector energy resolution' parameter", false);
						goto FAIL;
					}
	
				
					if (LbHdrGtElem(headerHk, HDR_DET_REFERENCE_ENERGY_ID,
							sizeof(paramsPtr->H.DetRunTimeParams.ReferenceEnergy),
							(void *)&paramsPtr->H.DetRunTimeParams.ReferenceEnergy) == false){
					
						PhgAbort("Header is missing 'simple detector reference energy' parameter", false);
						goto FAIL;
					}
					break;
					
				case DetEn_simple_spect:
					break;
					
				case DetEn_unc_spect:
					break;
					
				default:
					break;
			}
			
			/* the following detector parameters were added to the header after the early version.
			 Where possible we have filled them with default values rather than aborting when they
			 are not found (for backward compatibility) */
			/* get flag that tells whether interactions were forced in the detector */	
			if (LbHdrGtElem(headerHk, HDR_DET_IS_FORCED_INTERACTION_ID,
					sizeof(paramsPtr->H.DetRunTimeParams.DoForcedInteraction),
					(void *)&paramsPtr->H.DetRunTimeParams.DoForcedInteraction) == false){
				
					paramsPtr->H.DetRunTimeParams.DoForcedInteraction = true;	/* does least harm if wrong */

					/* Clear the error */
					ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}

			
			/* get parameter for the FWHM of photon detection time */	
			if (LbHdrGtElem(headerHk, HDR_DET_PHOTON_TIME_FWHM_ID,
					sizeof(paramsPtr->H.DetRunTimeParams.PhotonTimeFWHM),
					(void *)&paramsPtr->H.DetRunTimeParams.PhotonTimeFWHM) == false){
				
					paramsPtr->H.DetRunTimeParams.PhotonTimeFWHM = 0.0;	/* default value */

					/* Clear the error */
					ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}

			
			/* get parameter for the coincidence timing window */	
			if (LbHdrGtElem(headerHk, HDR_DET_COINC_TIMING_WINDOW_ID,
					sizeof(paramsPtr->H.DetRunTimeParams.CoincidenceTimingWindowNS),
					(void *)&paramsPtr->H.DetRunTimeParams.CoincidenceTimingWindowNS) == false){
				
					paramsPtr->H.DetRunTimeParams.CoincidenceTimingWindowNS = 0.0;	/* default value */

					/* Clear the error */
					ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}

			
			/* get parameter for the triples processing method */	
			if (LbHdrGtElem(headerHk, HDR_DET_TRIPLES_METHOD_ID,
					sizeof(paramsPtr->H.DetRunTimeParams.TriplesMethod),
					(void *)&paramsPtr->H.DetRunTimeParams.TriplesMethod) == false){
				
					paramsPtr->H.DetRunTimeParams.TriplesMethod = DetEn_DeleteTriples;	/* only value! */

					/* Clear the error */
					ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}

			
			/* get file name for the randoms-added history file */	
			if (LbHdrGtElem(headerHk, HDR_DET_RANDOMS_HISTORY_FILE_ID,
					PATH_LENGTH,
					(void *)&paramsPtr->H.DetRunTimeParams.DetRandomsHistoryFilePath) == false){
				
					paramsPtr->H.DetRunTimeParams.DetRandomsHistoryFilePath[0] = '\0';	/* default value */

					/* Clear the error */
					ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
			/* set the randoms processing flag if a file was given */
			paramsPtr->H.DetRunTimeParams.DoRandomsProcessing =
						( paramsPtr->H.DetRunTimeParams.DetRandomsHistoryFilePath[0] != '\0' );

			

		}
		} while (false);
		
		/************* GENERAL COMPUTED PARAMETERS ************/
		if (LbHdrGtElem(headerHk, HDR_CYL_TARGET_RADIUS_ID,
				sizeof(paramsPtr->H.TargetCylinder.radius),
				(void *)&paramsPtr->H.TargetCylinder.radius) == false){
			
			PhgAbort("Header is missing 'target cylinder radius' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_TARGET_MIN_Z_ID,
				sizeof(paramsPtr->H.TargetCylinder.zMin),
				(void *)&paramsPtr->H.TargetCylinder.zMin) == false){
			
			PhgAbort("Header is missing 'target cylinder z min' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_TARGET_MAX_Z_ID,
				sizeof(paramsPtr->H.TargetCylinder.zMax),
				(void *)&paramsPtr->H.TargetCylinder.zMax) == false){
			
			PhgAbort("Header is missing 'target cylinder z max' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_TARGET_CENTER_X_ID,
				sizeof(paramsPtr->H.TargetCylinder.centerX),
				(void *)&paramsPtr->H.TargetCylinder.centerX) == false){
			
			PhgAbort("Header is missing 'target cylinder center X' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_TARGET_CENTER_Y_ID,
				sizeof(paramsPtr->H.TargetCylinder.centerY),
				(void *)&paramsPtr->H.TargetCylinder.centerY) == false){
			
			PhgAbort("Header is missing 'target cylinder center Y' parameter", false);
			break;
		}

		
		if (LbHdrGtElem(headerHk, HDR_CYL_CRITZONE_RADIUS_ID,
				sizeof(paramsPtr->H.CriticalZone.radius),
				(void *)&paramsPtr->H.CriticalZone.radius) == false){
			
			PhgAbort("Header is missing 'critical zone radius' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_CRITZONE_MIN_Z_ID,
				sizeof(paramsPtr->H.CriticalZone.zMin),
				(void *)&paramsPtr->H.CriticalZone.zMin) == false){
			
			PhgAbort("Header is missing 'critical zone z min' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_CRITZONE_MAX_Z_ID,
				sizeof(paramsPtr->H.CriticalZone.zMax),
				(void *)&paramsPtr->H.CriticalZone.zMax) == false){
			
			PhgAbort("Header is missing 'critical zone z max' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_CRITZONE_CENTER_X_ID,
				sizeof(paramsPtr->H.CriticalZone.centerX),
				(void *)&paramsPtr->H.CriticalZone.centerX) == false){
			
			PhgAbort("Header is missing 'critical zone center X' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_CRITZONE_CENTER_Y_ID,
				sizeof(paramsPtr->H.CriticalZone.centerY),
				(void *)&paramsPtr->H.CriticalZone.centerY) == false){
			
			PhgAbort("Header is missing 'critical zone center Y' parameter", false);
			break;
		}

			
					
		if (LbHdrGtElem(headerHk, HDR_CYL_OBJCYL_RADIUS_ID,
				sizeof(paramsPtr->H.ObjectCylinder.radius),
				(void *)&paramsPtr->H.ObjectCylinder.radius) == false){
			
			PhgAbort("Header is missing 'critical zone radius' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_OBJCYL_MIN_Z_ID,
				sizeof(paramsPtr->H.ObjectCylinder.zMin),
				(void *)&paramsPtr->H.ObjectCylinder.zMin) == false){
			
			PhgAbort("Header is missing 'object cylinder z min' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_OBJCYL_MAX_Z_ID,
				sizeof(paramsPtr->H.ObjectCylinder.zMax),
				(void *)&paramsPtr->H.ObjectCylinder.zMax) == false){
			
			PhgAbort("Header is missing 'object cylinder z max' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_OBJCYL_CENTER_X_ID,
				sizeof(paramsPtr->H.ObjectCylinder.centerX),
				(void *)&paramsPtr->H.ObjectCylinder.centerX) == false){
			
			PhgAbort("Header is missing 'object cylinder center X' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_OBJCYL_CENTER_Y_ID,
				sizeof(paramsPtr->H.ObjectCylinder.centerY),
				(void *)&paramsPtr->H.ObjectCylinder.centerY) == false){
			
			PhgAbort("Header is missing 'object cylinder center Y' parameter", false);
			break;
		}

			
	
		if (LbHdrGtElem(headerHk, HDR_CYL_LIMITCYL_RADIUS_ID,
				sizeof(paramsPtr->H.LimitCylinder.radius),
				(void *)&paramsPtr->H.LimitCylinder.radius) == false){
			
			PhgAbort("Header is missing 'object cylinder radius' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_LIMITCYL_MIN_Z_ID,
				sizeof(paramsPtr->H.LimitCylinder.zMin),
				(void *)&paramsPtr->H.LimitCylinder.zMin) == false){
			
			PhgAbort("Header is missing 'limit cylinder z min' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_LIMITCYL_MAX_Z_ID,
				sizeof(paramsPtr->H.LimitCylinder.zMax),
				(void *)&paramsPtr->H.LimitCylinder.zMax) == false){
			
			PhgAbort("Header is missing 'limit cylinder z max' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_LIMITCYL_CENTER_X_ID,
				sizeof(paramsPtr->H.LimitCylinder.centerX),
				(void *)&paramsPtr->H.LimitCylinder.centerX) == false){
			
			PhgAbort("Header is missing 'limit cylinder center X' parameter", false);
			break;
		}

			
		if (LbHdrGtElem(headerHk, HDR_CYL_LIMITCYL_CENTER_Y_ID,
				sizeof(paramsPtr->H.LimitCylinder.centerY),
				(void *)&paramsPtr->H.LimitCylinder.centerY) == false){
			
			PhgAbort("Header is missing 'limit cylinder center Y' parameter", false);
			break;
		}

	
		if (LbHdrGtElem(headerHk, HDR_BIN1_NUM_SIMULATIONS_ID,
				sizeof(paramsPtr->H.NumSimulations),
				(void *)&paramsPtr->H.NumSimulations) == false){
			
			PhgAbort("Header is missing 'number of simulations' parameter", false);
			break;
		}

	
		/* save error status before fetching values that can be either 8 or 4 byte */
		saveInputError = ErIsInError();
		
		if (LbHdrGtElem(headerHk, HDR_BIN1_EVENTS_TO_SIMULATION_ID,
				sizeof(paramsPtr->H.SumEventsToSimulate),
				(void *)&paramsPtr->H.SumEventsToSimulate) == false){
			
			/* if the above call set ErIsInError(), clear it */
			if ( ErIsInError() && (!saveInputError) ) {
			
				ErClear();
			
			}
			
			/* see if this is an old header with a four-byte number */	
			if (LbHdrGtElem(headerHk, HDR_BIN1_EVENTS_TO_SIMULATION_ID,
					sizeof(oldFourByteValue),
					(void *)&oldFourByteValue) == false){
						
				PhgAbort("Header is missing 'sum of events to simulate' parameter", false);
				break;
				
			} else {
				
				/* assign 4-byte number of events to new field */
				paramsPtr->H.SumEventsToSimulate = oldFourByteValue;
				
			}

		}

	
		/* save error status before fetching values that can be either 8 or 4 byte */
		saveInputError = ErIsInError();
		
		if (LbHdrGtElem(headerHk, HDR_BIN1_NUM_PHOTONS_ID,
				sizeof(paramsPtr->H.NumPhotons),
				(void *)&paramsPtr->H.NumPhotons) == false){
			
			/* if the above call set ErIsInError(), clear it */
			if ( ErIsInError() && (!saveInputError) ) {
			
				ErClear();
			
			}
			
			/* see if this is an old header with a four-byte number */	
			if (LbHdrGtElem(headerHk, HDR_BIN1_NUM_PHOTONS_ID,
					sizeof(oldFourByteValue),
					(void *)&oldFourByteValue) == false){
						
				PhgAbort("Header is missing 'number of photons' parameter", false);
				break;
				
			} else {
				
				/* assign 4-byte number of events to new field */
				paramsPtr->H.NumPhotons = oldFourByteValue;
				
			}

		}

	
		/* save error status before fetching values that can be either 8 or 4 byte */
		saveInputError = ErIsInError();
		
		if (LbHdrGtElem(headerHk, HDR_BIN1_NUM_DECAYS_ID,
				sizeof(paramsPtr->H.NumDecays),
				(void *)&paramsPtr->H.NumDecays) == false){
			
			/* if the above call set ErIsInError(), clear it */
			if ( ErIsInError() && (!saveInputError) ) {
			
				ErClear();
			
			}
			
			/* see if this is an old header with a four-byte number */	
			if (LbHdrGtElem(headerHk, HDR_BIN1_NUM_DECAYS_ID,
					sizeof(oldFourByteValue),
					(void *)&oldFourByteValue) == false){
						
				PhgAbort("Header is missing 'number of decays' parameter", false);
				break;
				
			} else {
				
				/* assign 4-byte number of events to new field */
				paramsPtr->H.NumDecays = oldFourByteValue;
				
			}

		}

		/******** BIN PARAMETERS **********/
		do { /* Binning parameters */
		if ((paramsPtr->H.HdrKind == PhoHFileEn_BIN_WT) ||
				(paramsPtr->H.HdrKind == PhoHFileEn_BIN_WTSQ) ||
				(paramsPtr->H.HdrKind == PhoHFileEn_BIN_CT)){
				
				
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_Z_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numZBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numZBins) == false){
				
				PhgAbort("Header is missing 'number of z bins' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_PA_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numPABins),
					(void *)&paramsPtr->H.BinRunTimeParams.numPABins) == false){
				
				PhgAbort("Header is missing 'number of PA bins' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_TD_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numTDBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numTDBins) == false){
				
				PhgAbort("Header is missing 'number of TD bins' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_AA_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numAABins),
					(void *)&paramsPtr->H.BinRunTimeParams.numAABins) == false){
				
				PhgAbort("Header is missing 'number of AA bins' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_E1_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numE1Bins),
					(void *)&paramsPtr->H.BinRunTimeParams.numE1Bins) == false){
				
				PhgAbort("Header is missing 'number of E1 bins' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_E2_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numE2Bins),
					(void *)&paramsPtr->H.BinRunTimeParams.numE2Bins) == false){
				
				PhgAbort("Header is missing 'number of E2 bins' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_S1_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numS1Bins),
					(void *)&paramsPtr->H.BinRunTimeParams.numS1Bins) == false){
				
				PhgAbort("Header is missing 'number of S1 bins' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_S2_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numS2Bins),
					(void *)&paramsPtr->H.BinRunTimeParams.numS2Bins) == false){
				
				PhgAbort("Header is missing 'number of S2 bins' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_IMAGE_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numImageBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numImageBins) == false){
				
				PhgAbort("Header is missing 'number of image bins' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_SCATTER_PARAM_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatterRandomParam),
					(void *)&paramsPtr->H.BinRunTimeParams.scatterRandomParam) == false){
				
				PhgAbort("Header is missing 'scatter binning' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_MIN_Z_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minZ),
				 	(void *)&paramsPtr->H.BinRunTimeParams.minZ) == false){
				
				PhgAbort("Header is missing 'minimum Z' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_MAX_Z_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxZ),
					(void *)&paramsPtr->H.BinRunTimeParams.maxZ) == false){
				
				PhgAbort("Header is missing 'maximum Z' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_MIN_PA_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minPA),
					(void *)&paramsPtr->H.BinRunTimeParams.minPA) == false){
				
				PhgAbort("Header is missing 'minimum PA' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_MAX_PA_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxPA),
					(void *)&paramsPtr->H.BinRunTimeParams.maxPA) == false){
				
				PhgAbort("Header is missing 'maximum PA' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_MIN_TD_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minTD),
					(void *)&paramsPtr->H.BinRunTimeParams.minTD) == false){
				
				PhgAbort("Header is missing 'minimum TD' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_MAX_TD_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxTD),
					(void *)&paramsPtr->H.BinRunTimeParams.maxTD) == false){
				
				PhgAbort("Header is missing 'maximum TD' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_MIN_AA_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minAA),
					(void *)&paramsPtr->H.BinRunTimeParams.minAA) == false){
				
				PhgAbort("Header is missing 'minimum AA' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_MAX_AA_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxAA),
					(void *)&paramsPtr->H.BinRunTimeParams.maxAA) == false){
				
				PhgAbort("Header is missing 'maximum AA' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_MIN_ENERGY_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minE),
					(void *)&paramsPtr->H.BinRunTimeParams.minE) == false){
				
				PhgAbort("Header is missing 'minimum energy' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_MAX_ENERGY_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxE),
					(void *)&paramsPtr->H.BinRunTimeParams.maxE) == false){
				
				PhgAbort("Header is missing 'maximum energy' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_MIN_SCAT_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minS),
					(void *)&paramsPtr->H.BinRunTimeParams.minS) == false){
				
				PhgAbort("Header is missing 'minimum scatters' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_MAX_SCAT_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxS),
					(void *)&paramsPtr->H.BinRunTimeParams.maxS) == false){
				
				PhgAbort("Header is missing 'maximum scatters' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_ADD_TO_EXISTING_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.addToExistingImg),
					(void *)&paramsPtr->H.BinRunTimeParams.addToExistingImg) == false){
				
				PhgAbort("Header is missing 'add to existing images' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_DO_COUNTS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.doCounts),
					(void *)&paramsPtr->H.BinRunTimeParams.doCounts) == false){
				
				PhgAbort("Header is missing 'do counts' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_DO_WEIGHTS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.doWeights),
					(void *)&paramsPtr->H.BinRunTimeParams.doWeights) == false){
				
				PhgAbort("Header is missing 'do weights' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_DO_WEIGHTS_SQU_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.doWeightsSquared),
					(void *)&paramsPtr->H.BinRunTimeParams.doWeightsSquared) == false){
				
				PhgAbort("Header is missing 'do weights squared' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_Z_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.zRange),
					(void *)&paramsPtr->H.BinRunTimeParams.zRange) == false){
				
				PhgAbort("Header is missing 'z range' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_E_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.eRange),
					(void *)&paramsPtr->H.BinRunTimeParams.eRange) == false){
				
				PhgAbort("Header is missing 'energy range' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_S_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.sRange),
					(void *)&paramsPtr->H.BinRunTimeParams.sRange) == false){
				
				PhgAbort("Header is missing 'scatter range' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_TD_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tdRange),
					(void *)&paramsPtr->H.BinRunTimeParams.tdRange) == false){
				
				PhgAbort("Header is missing 'td range' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_SCATTER2CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatter2CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.scatter2CIsize) == false){
				
				PhgAbort("Header is missing 'scatter2 count size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_SCATTER2WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatter2WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.scatter2WIsize) == false){
				
				PhgAbort("Header is missing 'scatter2 weight size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_SCATTER2WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatter2WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.scatter2WISsize) == false){
				
				PhgAbort("Header is missing 'scatter2 weight squared size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_SCATTER1CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatter1CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.scatter1CIsize) == false){
				
				PhgAbort("Header is missing 'scatter1 count size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_SCATTER1WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatter1WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.scatter1WIsize) == false){
				
				PhgAbort("Header is missing 'scatter1 weight size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_SCATTER1WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatter1WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.scatter1WISsize) == false){
				
				PhgAbort("Header is missing 'scatter1 weight squared size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_ENERGY2CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.energy2CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.energy2CIsize) == false){
				
				PhgAbort("Header is missing 'energy2 count size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_ENERGY2WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.energy2WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.energy2WIsize) == false){
				
				PhgAbort("Header is missing 'energy2 weight size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_ENERGY2WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.energy2WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.energy2WISsize) == false){
				
				PhgAbort("Header is missing 'energy2 weight squared size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_ENERGY1CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.energy1CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.energy1CIsize) == false){
				
				PhgAbort("Header is missing 'energy1 count size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_ENERGY1WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.energy1WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.energy1WIsize) == false){
				
				PhgAbort("Header is missing 'energy1 weight size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_ENERGY1WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.energy1WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.energy1WISsize) == false){
				
				PhgAbort("Header is missing 'energy1 weight squared size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_AACI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.aaCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.aaCIsize) == false){
				
				PhgAbort("Header is missing 'aa count size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_AAWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.aaWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.aaWIsize) == false){
				
				PhgAbort("Header is missing 'aa weight size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_AAWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.aaWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.aaWISsize) == false){
				
				PhgAbort("Header is missing 'aa weight squared size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_TDCI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tdCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.tdCIsize) == false){
				
				PhgAbort("Header is missing 'td count size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_TDWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tdWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.tdWIsize) == false){
				
				PhgAbort("Header is missing 'td weight size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_TDWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tdWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.tdWISsize) == false){
				
				PhgAbort("Header is missing 'td weight squared size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_PACI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.paCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.paCIsize) == false){
				
				PhgAbort("Header is missing 'pa count size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_PAWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.paWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.paWIsize) == false){
				
				PhgAbort("Header is missing 'pa weight size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_PAWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.paWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.paWISsize) == false){
				
				PhgAbort("Header is missing 'pa weight squared size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_Z2CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.z2CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.z2CIsize) == false){
				
				PhgAbort("Header is missing 'z2 count size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_Z2WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.z2WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.z2WIsize) == false){
				
				PhgAbort("Header is missing 'z2 weight size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_Z2WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.z2WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.z2WISsize) == false){
				
				PhgAbort("Header is missing 'z2 weight squared size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_Z1CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.z1CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.z1CIsize) == false){
				
				PhgAbort("Header is missing 'z1 count size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_Z1WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.z1WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.z1WIsize) == false){
				
				PhgAbort("Header is missing 'z1 weight size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_Z1WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.z1WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.z1WISsize) == false){
				
				PhgAbort("Header is missing 'z1 weight squared size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_COUNT_IMG_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.countImageSize),
					(void *)&paramsPtr->H.BinRunTimeParams.countImageSize) == false){
				
				PhgAbort("Header is missing 'count image size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_WEIGHT_IMG_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.weightImageSize),
					(void *)&paramsPtr->H.BinRunTimeParams.weightImageSize) == false){
				
				PhgAbort("Header is missing 'weight image size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_WEIGHT_SQU_IMG_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.weightSquImageSize),
					(void *)&paramsPtr->H.BinRunTimeParams.weightSquImageSize) == false){
				
				PhgAbort("Header is missing 'weight squared image size' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_WEIGHT_IMG_TYPE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.weight_image_type),
					(void *)&paramsPtr->H.BinRunTimeParams.weight_image_type) == false){
				
				PhgAbort("Header is missing 'weight image type' parameter", false);
				goto FAIL;
			}
	
		
			if (LbHdrGtElem(headerHk, HDR_BIN_COUNT_IMG_TYPE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.count_image_type),
					(void *)&paramsPtr->H.BinRunTimeParams.count_image_type) == false){
				
				PhgAbort("Header is missing 'count image type' parameter", false);
				goto FAIL;
			}
	
			/* Get the dimension ordering for the binning process */
			if (LbHdrGtElem(headerHk, HDR_BIN_DIMENSION_ORDER_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.PhgBinDimensions),
					(void *)&paramsPtr->H.BinRunTimeParams.PhgBinDimensions) == false){
				
				PhgAbort("Header is missing 'dimension ordering' parameter", false);
				goto FAIL;
			}
			
			/* THESE ARE NEW PARAMETERS ADDED TO THE HEADER BEFORE THE PUBLIC RELEASE
				BUT AFTER MANY FILES HAVE BEEN CREATED WITHIN THE IRL  For this purpose
				note that these fields are not going to be required to exist, if they
				don't we just set them to an hopefully benign value
			*/
		
			/* Get the dimension ordering for the binning process */
			if (LbHdrGtElem(headerHk, HDR_BIN_PHICI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.phiCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.phiCIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.phiCIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for the binning process */
			if (LbHdrGtElem(headerHk, HDR_BIN_PHIWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.phiWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.phiWIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.phiWIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for the binning process */
			if (LbHdrGtElem(headerHk, HDR_BIN_PHIWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.phiWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.phiWISsize) == false){
				
				paramsPtr->H.BinRunTimeParams.phiWISsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for the binning process */
			if (LbHdrGtElem(headerHk, HDR_BIN_THETACI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.thetaCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.thetaCIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.thetaCIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for the binning process */
			if (LbHdrGtElem(headerHk, HDR_BIN_THETAWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.thetaWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.thetaWIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.thetaWIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for size of the theta dimension */
			if (LbHdrGtElem(headerHk, HDR_BIN_THETAWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.thetaWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.thetaWISsize) == false){
				
				paramsPtr->H.BinRunTimeParams.thetaWISsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the size of the Xr dimension */
			if (LbHdrGtElem(headerHk, HDR_BIN_XRCI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.xrCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.xrCIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.xrCIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the size of the Xr weight dimension */
			if (LbHdrGtElem(headerHk, HDR_BIN_XRWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.xrWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.xrWIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.xrWIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the size of the Xr weight squared dimension */
			if (LbHdrGtElem(headerHk, HDR_BIN_XRWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.xrWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.xrWISsize) == false){
				
				paramsPtr->H.BinRunTimeParams.xrWISsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the size of the Yr count dimension */
			if (LbHdrGtElem(headerHk, HDR_BIN_YRCI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.yrCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.yrCIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.yrCIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the size of the Yr weight dimension */
			if (LbHdrGtElem(headerHk, HDR_BIN_YRWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.yrWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.yrWIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.yrWIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the size of the Yr weight squared dimension */
			if (LbHdrGtElem(headerHk, HDR_BIN_YRWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.yrWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.yrWISsize) == false){
				
				paramsPtr->H.BinRunTimeParams.yrWISsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the number of PHI bins */
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_PHI_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numPHIBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numPHIBins) == false){
				
				paramsPtr->H.BinRunTimeParams.numPHIBins = 0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the number of theta bins */
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_THETA_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numThetaBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numThetaBins) == false){
				
				paramsPtr->H.BinRunTimeParams.numThetaBins = 0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the number of Xr bins */
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_XR_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numXRBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numXRBins) == false){
				
				paramsPtr->H.BinRunTimeParams.numXRBins = 0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the number of Yr bins */
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_YR_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numYRBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numYRBins) == false){
				
				paramsPtr->H.BinRunTimeParams.numYRBins = 0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the min theta */
			if (LbHdrGtElem(headerHk, HDR_BIN_MIN_THETA_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minTheta),
					(void *)&paramsPtr->H.BinRunTimeParams.minTheta) == false){
				
				paramsPtr->H.BinRunTimeParams.minTheta = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get max theta */
			if (LbHdrGtElem(headerHk, HDR_BIN_MAX_THETA_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxTheta),
					(void *)&paramsPtr->H.BinRunTimeParams.maxTheta) == false){
				
				paramsPtr->H.BinRunTimeParams.maxTheta = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get min PHI */
			if (LbHdrGtElem(headerHk, HDR_BIN_MIN_PHI_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minPhi),
					(void *)&paramsPtr->H.BinRunTimeParams.minPhi) == false){
				
				paramsPtr->H.BinRunTimeParams.minPhi = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get max PHI */
			if (LbHdrGtElem(headerHk, HDR_BIN_MAX_PHI_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxPhi),
					(void *)&paramsPtr->H.BinRunTimeParams.maxPhi) == false){
				
				paramsPtr->H.BinRunTimeParams.maxPhi = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get min Xr */
			if (LbHdrGtElem(headerHk, HDR_BIN_MIN_XR_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minXR),
					(void *)&paramsPtr->H.BinRunTimeParams.minXR) == false){
				
				paramsPtr->H.BinRunTimeParams.minXR = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get max Xr */
			if (LbHdrGtElem(headerHk, HDR_BIN_MAX_XR_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxXR),
					(void *)&paramsPtr->H.BinRunTimeParams.maxXR) == false){
				
				paramsPtr->H.BinRunTimeParams.maxXR = -500.0;
	
				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get min Yr */
			if (LbHdrGtElem(headerHk, HDR_BIN_MIN_YR_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minYR),
					(void *)&paramsPtr->H.BinRunTimeParams.minYR) == false){
				
				paramsPtr->H.BinRunTimeParams.minYR = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get max Yr */
			if (LbHdrGtElem(headerHk, HDR_BIN_MAX_YR_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxYR),
					(void *)&paramsPtr->H.BinRunTimeParams.maxYR) == false){
				
				paramsPtr->H.BinRunTimeParams.maxYR = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get SSRB flag */
			if (LbHdrGtElem(headerHk, HDR_BIN_DOSSRB_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.doSSRB),
					(void *)&paramsPtr->H.BinRunTimeParams.doSSRB) == false){
				
				paramsPtr->H.BinRunTimeParams.doSSRB = false;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get MSRB flag */
			if (LbHdrGtElem(headerHk, HDR_BIN_DOMSRB_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.doMSRB),
					(void *)&paramsPtr->H.BinRunTimeParams.doMSRB) == false){
				
				paramsPtr->H.BinRunTimeParams.doMSRB = false;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the image radius */
			if (LbHdrGtElem(headerHk, HDR_BIN_IMG_RADIUS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.image_radius),
					(void *)&paramsPtr->H.BinRunTimeParams.image_radius) == false){
				
				paramsPtr->H.BinRunTimeParams.image_radius = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the detector radius */
			if (LbHdrGtElem(headerHk, HDR_BIN_DET_RADIUS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.detector_radius),
					(void *)&paramsPtr->H.BinRunTimeParams.detector_radius) == false){
				
				paramsPtr->H.BinRunTimeParams.detector_radius = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}

	
			/* Get the weight image path */
			if (LbHdrGtElem(headerHk, HDR_BIN_WTIMG_PATH_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.weightImgFilePath),
					(void *)&paramsPtr->H.BinRunTimeParams.weightImgFilePath) == false){
				
				paramsPtr->H.BinRunTimeParams.weightImgFilePath[0] = '\0';

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the weight squared image path */
			if (LbHdrGtElem(headerHk, HDR_BIN_WTSIMG_PATH_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.weightSquImgFilePath),
					(void *)&paramsPtr->H.BinRunTimeParams.weightSquImgFilePath) == false){
				
				paramsPtr->H.BinRunTimeParams.weightSquImgFilePath[0] = '\0';

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the count image path */
			if (LbHdrGtElem(headerHk, HDR_BIN_CTIMG_PATH_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.countImgFilePath),
					(void *)&paramsPtr->H.BinRunTimeParams.countImgFilePath) == false){
				
				paramsPtr->H.BinRunTimeParams.countImgFilePath[0] = '\0';

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the theta range */
			if (LbHdrGtElem(headerHk, HDR_BIN_THETA_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.thetaRange),
					(void *)&paramsPtr->H.BinRunTimeParams.thetaRange) == false){
				
				paramsPtr->H.BinRunTimeParams.thetaRange = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the phi range */
			if (LbHdrGtElem(headerHk, HDR_BIN_PHI_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.phiRange),
					(void *)&paramsPtr->H.BinRunTimeParams.phiRange) == false){
				
				paramsPtr->H.BinRunTimeParams.phiRange = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the Xr range */
			if (LbHdrGtElem(headerHk, HDR_BIN_XR_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.xrRange),
					(void *)&paramsPtr->H.BinRunTimeParams.xrRange) == false){
				
				paramsPtr->H.BinRunTimeParams.xrRange = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the Yr range */
			if (LbHdrGtElem(headerHk, HDR_BIN_YR_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.yrRange),
					(void *)&paramsPtr->H.BinRunTimeParams.yrRange) == false){
				
				paramsPtr->H.BinRunTimeParams.yrRange = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the sum weights  */
			if (LbHdrGtElem(headerHk, HDR_BIN_SUM_WEIGHTS_ID,
					sizeof(paramsPtr->H.weightSum),
					(void *)&paramsPtr->H.weightSum) == false){
				
				paramsPtr->H.weightSum = 0.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
			
	
			/* Get the sum weights squared  */
			if (LbHdrGtElem(headerHk, HDR_BIN_SUM_WEIGHTS_SQ_ID,
					sizeof(paramsPtr->H.weightSquSum),
					(void *)&paramsPtr->H.weightSquSum) == false){
				
				paramsPtr->H.weightSquSum = 0.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
			

			/* THESE ARE EVEN NEWER PARAMETERS ADDED TO THE HEADER AFTER PUBLIC RELEASE
				To help with backward compatibility,
				note that these fields are not going to be required to exist, if they
				don't we just set them to an hopefully benign value
			*/
		
			/* Get the number of TOF bins */
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_TOF_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numTOFBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numTOFBins) == false){
				
				paramsPtr->H.BinRunTimeParams.numTOFBins = 0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the min TOF */
			if (LbHdrGtElem(headerHk, HDR_BIN_MIN_TOF_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minTOF),
					(void *)&paramsPtr->H.BinRunTimeParams.minTOF) == false){
				
				paramsPtr->H.BinRunTimeParams.minTOF = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get max TOF */
			if (LbHdrGtElem(headerHk, HDR_BIN_MAX_TOF_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxTOF),
					(void *)&paramsPtr->H.BinRunTimeParams.maxTOF) == false){
				
				paramsPtr->H.BinRunTimeParams.maxTOF = -500.0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for the binning process */
			if (LbHdrGtElem(headerHk, HDR_BIN_TOFCI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tofCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.tofCIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.tofCIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for the binning process */
			if (LbHdrGtElem(headerHk, HDR_BIN_TOFWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tofWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.tofWIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.tofWIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for size of the tof dimension */
			if (LbHdrGtElem(headerHk, HDR_BIN_TOFWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tofWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.tofWISsize) == false){
				
				paramsPtr->H.BinRunTimeParams.tofWISsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get crystal binning flag */
			if (LbHdrGtElem(headerHk, HDR_BIN_DO_CRYSTAL_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.doCrystal),
					(void *)&paramsPtr->H.BinRunTimeParams.doCrystal) == false){
				
				paramsPtr->H.BinRunTimeParams.doCrystal = false;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the number of crystal bins */
			if (LbHdrGtElem(headerHk, HDR_BIN_NUM_CRYSTAL_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numCrystalBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numCrystalBins) == false){
				
				paramsPtr->H.BinRunTimeParams.numCrystalBins = 0;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for the binning process */
			if (LbHdrGtElem(headerHk, HDR_BIN_CRYSTAL1CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.crystal1CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.crystal1CIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.crystal1CIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for the binning process */
			if (LbHdrGtElem(headerHk, HDR_BIN_CRYSTAL1WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.crystal1WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.crystal1WIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.crystal1WIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for size of the crystal dimension */
			if (LbHdrGtElem(headerHk, HDR_BIN_CRYSTAL1WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.crystal1WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.crystal1WISsize) == false){
				
				paramsPtr->H.BinRunTimeParams.crystal1WISsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for the binning process */
			if (LbHdrGtElem(headerHk, HDR_BIN_CRYSTAL2CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.crystal2CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.crystal2CIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.crystal2CIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for the binning process */
			if (LbHdrGtElem(headerHk, HDR_BIN_CRYSTAL2WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.crystal2WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.crystal2WIsize) == false){
				
				paramsPtr->H.BinRunTimeParams.crystal2WIsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			/* Get the dimension ordering for size of the crystal dimension */
			if (LbHdrGtElem(headerHk, HDR_BIN_CRYSTAL2WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.crystal2WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.crystal2WISsize) == false){
				
				paramsPtr->H.BinRunTimeParams.crystal2WISsize = 1;

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
	
			
		}
		} while (false);
		
		/* History file parameters.  Added after early SimSET versions, so default parameters
		 are given where possible rather than aborting (for backward compatibility). */
		do {
			/* Get the timesorted indicator  */
			if (LbHdrGtElem(headerHk, HDR_HISTORY_FILE_IS_SORTED_ID,
					sizeof(paramsPtr->H.isTimeSorted),
					(void *)&paramsPtr->H.isTimeSorted) == false){
				
				paramsPtr->H.isTimeSorted = false;	/* Most likely */

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
			
			/* Get the randoms added indicator  */
			if (LbHdrGtElem(headerHk, HDR_HISTORY_FILE_IS_RANDOMS_ADDED_ID,
					sizeof(paramsPtr->H.isRandomsAdded),
					(void *)&paramsPtr->H.isRandomsAdded) == false){
				
				paramsPtr->H.isRandomsAdded = false;	/* Most likely */

				/* Clear the error */
				ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			}
			
		} while (false);
		okay = true;
		FAIL:;
	} while (false);

	/* Return the status */
	return (okay);
}

/*********************************************************************************
*
*			Name:			PhgHdrCreateRunTimeHeader
*
*			Summary:		Initializes a header according to current simulation
*							parameters
*
*			Arguments:
*				PhoHFileHdrKindTy	hdrKind			- The type of the header 
*				PhoHFileHdrTy 		*hdrTyPtr		- The header initialized
*				PHG_BinParamsTy		*binParams		- The binning parameters
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean PhgHdrCreateRunTimeHeader(PhoHFileHdrKindTy hdrKind, PhoHFileHdrTy *hdrTyPtr,
			PHG_BinParamsTy *binParams)	
{
	Boolean	okay = false;		/* Process flags */
	
	do { /* Process Loop */

		
		/* Initialize the header */
		hdrTyPtr->H.HdrKind = hdrKind;
		hdrTyPtr->H.HdrSize = PHG_HDR_HEADER_SIZE;
		hdrTyPtr->H.HdrVersion = PHG_HDR_HEADER_VERSION;
		hdrTyPtr->H.TargetCylinder = CylPosTargetCylinder;
		hdrTyPtr->H.CriticalZone = CylPosCriticalZone;
		hdrTyPtr->H.ObjectCylinder = CylPosObjectCylinder;
		hdrTyPtr->H.LimitCylinder = CylPosLimitCylinder;
		hdrTyPtr->H.NumSimulations = 0;
		hdrTyPtr->H.SumEventsToSimulate = PhgRunTimeParams.Phg_EventsToSimulate;
		hdrTyPtr->H.NumPhotons = 0;
		hdrTyPtr->H.NumDecays = 0;
		hdrTyPtr->H.weightSum = 0.0;
		hdrTyPtr->H.weightSquSum = 0.0;
		hdrTyPtr->H.isTimeSorted = false;
		hdrTyPtr->H.isRandomsAdded = false;
		
		/* Copy Runtime parameters */
		memcpy(&(hdrTyPtr->H.PhgRunTimeParams),
			&PhgRunTimeParams, sizeof(PhgRunTimeParamsTy));
		
		/* If we are doing collimation then copy the collimator params */
		if (PHG_IsCollimateOnTheFly()) {
			memcpy(&(hdrTyPtr->H.ColRunTimeParams),
				&ColRunTimeParams, sizeof(ColRunTimeParamsTy));
		}
		else {
			memset(&(hdrTyPtr->H.ColRunTimeParams),
				'\0', sizeof(ColRunTimeParamsTy));
		}

		/* If we are doing detection then copy the detector params */
		if (PHG_IsDetectOnTheFly()) {
			memcpy(&(hdrTyPtr->H.DetRunTimeParams),
				&DetRunTimeParams, sizeof(DetRunTimeParamsTy));
		}
		else {
			memset(&(hdrTyPtr->H.DetRunTimeParams),
				'\0', sizeof(DetRunTimeParamsTy));
		}

		/* If we are doing binning then copy the binning params */
		if (PHG_IsBinOnTheFly() && (binParams != 0)) {
			memcpy(&(hdrTyPtr->H.BinRunTimeParams),
				binParams, sizeof(PHG_BinParamsTy));
		}
		else {
			memset(&(hdrTyPtr->H.BinRunTimeParams),
				'\0', sizeof(PHG_BinParamsTy));
		}

		okay = true;
	} while (false);
	return (okay);
}

/**********************
*	PhgHdrStFile
*
*	Purpose:	This routine associates a header with a different file.
*				the previous contents of the file's header (if it has one)
*				are wiped out.
*
*	Arguments:
*			LbHdrHkTy		*headerHk	- The newly created header structure.
*			FILE			*headerFl	- The new file to associate with
*
*	Result:	TRUE unless an error occurs.
***********************/
Boolean  PhgHdrStFile(LbHdrHkTy *headerHk, FILE *headerFl)
{
	return(LbHdrStFile(headerHk, headerFl));
}

/**********************
*	PhgHdrFrHeader
*
*	Purpose:	This routine frees the memory associated with a header.
*
*	Arguments:
*			LbHdrHkTy		*headerHk	- The newly created header structure.
*
*	Result:	None.
***********************/
void  PhgHdrFrHeader(LbHdrHkTy *headerHk)
{
	LbHdrFree(headerHk);	
}

/**********************
*	PhgHdrMkHeader
*
*	Purpose:	This routine takes a file pointer, a run time params structure ptr,
*				and creates a header structure to be returned. The header is written
*				to the file, the current position within the file is left where it was
*				when the routine was called.
*
*	Arguments:
*			FILE			*headerFile	- The file the header will be attatched to.
*			PhoHFileHdrTy	*paramsPtr	- The run time parameters structure.
*			LbHdrHkTy		*headerHk	- The newly created header structure.
*
*	Result:	True unless an error occurs.
***********************/
Boolean  PhgHdrMkHeader(FILE *headerFile, PhoHFileHdrTy *paramsPtr, LbHdrHkTy *headerHk)
{
	Boolean				okay = false;			/* Process Loop */	

	do { /* Process Loop */
			
		/* Create the header within the file */
		if (LbHdrNew(headerFile, PHG_HDR_HEADER_SIZE, headerHk) == false){
			PhgAbort("Unable to create new header (PhgHdrMkHeader)", false);
			break;
		}
		
		/* Update the header */
		if (PhgHdrUpHeader(headerFile, paramsPtr, headerHk) == false) {
			break;
		}
		
		okay = true;
	} while (false);
	
	return (okay);
}

/**********************
*	PhgHdrUpHeader
*
*	Purpose:	This routine takes a file pointer, a run time params structure ptr,
*				and an existing header structure. It updates the header with the
*				current run time parameters and then writes the header to disk.
*
*	Arguments:
*			FILE			*headerFile	- The file the header will be attatched to.
*			PhoHFileHdrTy	*paramsPtr	- The run time parameters structure.
*			LbHdrHkTy		*headerHk	- The newly created header structure.
*
*	Result:	True unless an error occurs.
***********************/
Boolean  PhgHdrUpHeader(FILE *headerFile, PhoHFileHdrTy *paramsPtr, LbHdrHkTy *headerHk)
{
	Boolean				okay = false;			/* Process Loop */	
	FILE*				dummy;					/* Local copy of headerFile */
	
	dummy = headerFile;		/* Removes unused variable compiler warning */
	
	do { /* Process Loop */
			
		/* Set all of the fields */
		if (PhgHdrStFields(paramsPtr, headerHk) == false) {
			break;
		}

		/* Now write the header to the file */
		if (LbHdrWrite(headerHk) == false){
			
			PhgAbort("Unable to write header to file (PhgHdrUpHeader)", false);
			break;
		}
		
		okay = true;
	} while (false);
	
	return (okay);
}

/**********************
*	PhgHdrStFields
*
*	Purpose:	This routine sets all of the run time parameter fields into
*				the file header hook.
*
*	Arguments:
*			PhoHFileHdrTy	*paramsPtr	- The run time parameters structure.
*			LbHdrHkTy		*headerHk	- The newly created header structure.
*
*	Result:	True unless an error occurs.
***********************/
Boolean  PhgHdrStFields(PhoHFileHdrTy *paramsPtr, LbHdrHkTy *headerHk)
{
	Boolean				okay = false;			/* Process Loop */	
	LbUsTwoByte			isotopeConverter;		/* Conversion for isotope */
	
	do { /* Process Loop */
		
		/* One-by-one, get the header fields */

		/* Set the header size */
		if (LbHdrStElem(headerHk, HDR_PHG_HEADER_SIZE_ID, sizeof(paramsPtr->H.HdrSize),
				(void *)&(paramsPtr->H.HdrSize)) != true) {
			
			PhgAbort("Unable to set header 'size' parameter", false);
			break;
		}

		/* Set the header kind */
		if (LbHdrStElem(headerHk, HDR_PHG_HEADER_KIND_ID, sizeof(paramsPtr->H.HdrKind),
				(void *)&(paramsPtr->H.HdrKind)) != true) {
			
			PhgAbort("Unable to set header 'kind' parameter", false);
			break;
		}

		/* Set the header version */
		if (LbHdrStElem(headerHk, HDR_PHG_HEADER_VERS_ID, sizeof(paramsPtr->H.HdrVersion),
				(void *)&(paramsPtr->H.HdrVersion)) != true) {
			
			PhgAbort("Unable to set header 'version' parameter", false);
			break;
		}

		if (LbHdrStElem(headerHk, HDR_PHG_EVENTS_TO_SIM_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.Phg_EventsToSimulate),
				(void *)&paramsPtr->H.PhgRunTimeParams.Phg_EventsToSimulate) == false){		
			PhgAbort("Unable to set header 'events to simulate' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_CALC_EVENTS_TO_SIM_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsCalcEventsToSimulate),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsCalcEventsToSimulate) == false){		
			PhgAbort("Unable to set header 'calculate events to simulate' flag", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_RAND_SEED_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgRandomSeed),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgRandomSeed) == false){
			
			PhgAbort("Unable to set header 'random seed' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_LEN_OF_SCAN_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.Phg_LengthOfScan),
				(void *)&paramsPtr->H.PhgRunTimeParams.Phg_LengthOfScan) == false){
			
			PhgAbort("Unable to set header 'length of scan' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_ACC_ANGLE_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.Phg_AcceptanceAngle),
				(void *)&paramsPtr->H.PhgRunTimeParams.Phg_AcceptanceAngle) == false){
			
			PhgAbort("Unable to set header 'acceptance angle' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_ACC_ANGLE_SINE_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.Phg_SineOfAcceptanceAngle),
				(void *)&paramsPtr->H.PhgRunTimeParams.Phg_SineOfAcceptanceAngle) == false){
			
			PhgAbort("Unable to set header 'sine of acceptance angle' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_MIN_ENERGY_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgMinimumEnergy),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgMinimumEnergy) == false){
			
			PhgAbort("Unable to set header 'minimum energy' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_MIN_WW_RATIO_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgMinWWRatio),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgMinWWRatio) == false){
			
			PhgAbort("Unable to set header 'minimum weight window ratio' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_MAX_WW_RATIO_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgMaxWWRatio),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgMaxWWRatio) == false){
			
			PhgAbort("Unable to set header 'maximum weight window ratio' parameter", false);
			break;
		}

		/* Convert isotope enumerate to two byte for backwards compatibility */
		isotopeConverter = (LbUsTwoByte) paramsPtr->H.PhgRunTimeParams.PhgNuclide.isotope;
		
		if (LbHdrStElem(headerHk, HDR_PHG_ISOTOPE_ID,
				sizeof(isotopeConverter),
				(void *)&isotopeConverter) == false){
			
			PhgAbort("Unable to set header 'nuclide isotope' parameter", false);
			break;
		}
	
			
		if (LbHdrStElem(headerHk, HDR_PHG_PHOTON_ENERGY_KEV_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgNuclide.photonEnergy_KEV),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgNuclide.photonEnergy_KEV) == false){
			
			PhgAbort("Unable to set header 'photon energy' parameter", false);
			break;
		}
		
			
		if (LbHdrStElem(headerHk, HDR_PHG_POSITRON_ENERGY_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgNuclide.maximumPositronEnergy),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgNuclide.maximumPositronEnergy) == false){
			
			PhgAbort("Unable to set header 'maximum positron energy' parameter", false);
			break;
		}

			
		
		if (LbHdrStElem(headerHk, HDR_PHG_IS_FORCED_DETECTION_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsForcedDetection),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsForcedDetection) == false){
			
			PhgAbort("Unable to set header 'forced detection on/off' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_STRATIFICATION_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsStratification),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsStratification) == false){
			
			PhgAbort("Unable to set header 'stratification on/off' parameter", false);
			break;
		}

		
		if (LbHdrStElem(headerHk, HDR_PHG_IS_NON_ABSORPTION_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsForcedNonAbsorbtion),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsForcedNonAbsorbtion) == false){
			
			PhgAbort("Unable to set header 'forced non-absorption on/off' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_SPECT_ID,
				sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsSPECT),
				(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsSPECT) == false){
			
			PhgAbort("Unable to set header 'SPECT on/off' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_PET_COINC_ONLY_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincidencesOnly),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincidencesOnly) == false){
			
			PhgAbort("Unable to set header 'PET coincidences only on/off' parameter", false);
			break;
		}
		
			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_PET_COINC_PLUS_SINGLES_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincPlusSingles),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincPlusSingles) == false){
			
			PhgAbort("Unable to set header 'PET coincidences plus singles on/off' parameter", false);
			break;
		}
		
			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_MULTIEMISSION_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsMultiEmission),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsMultiEmission) == false){
			
			PhgAbort("Unable to set header 'Multi-emission on/off' parameter", false);
			break;
		}
		
			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_PET_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsPET),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsPET) == false){
			
			PhgAbort("Unable to set header 'PET on/off' parameter", false);
			break;
		}
		
			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_PET_COINC_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincidences),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsPETCoincidences) == false){
			
			PhgAbort("Unable to set header 'PET coincidences on/off' parameter", false);
			break;
		}
		
			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_PET_SINGLES_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsPETSingles),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsPETSingles) == false){
			
			PhgAbort("Unable to set header 'PET singles on/off' parameter", false);
			break;
		}
		
			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_HISTORY_FILE_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsHistoryFile),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsHistoryFile) == false){
			
			PhgAbort("Unable to set header 'history file on/off' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_ADJ_FOR_POS_RANGE_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsAdjForPosRange),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsAdjForPosRange) == false){
			
			PhgAbort("Unable to set header 'adjust for positron range' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_ADJ_FOR_COLIN_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsAdjForCollinearity),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsAdjForCollinearity) == false){
			
			PhgAbort("Unable to set header 'adjust for collinearity' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_COMPUTED_PROD_TBL_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsComputedProductivityTbl),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsComputedProductivityTbl) == false){
			
			PhgAbort("Unable to set header 'computed productivity' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_VOXEL_PTSRC_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsVoxelPointSource),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsVoxelPointSource) == false){
			
			PhgAbort("Unable to set header 'point source voxels' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_VOXEL_LNSRC_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsVoxelLineSource),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsVoxelLineSource) == false){
			
			PhgAbort("Unable to set header 'line source voxels' parameter", false);
			break;
		}
			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_POLARIZATION_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsModelPolarization),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsModelPolarization) == false){
			
			PhgAbort("Unable to set header 'is polarization' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_BIN_ON_THE_FLY_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsBinOnTheFly),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsBinOnTheFly) == false){
			
			PhgAbort("Unable to set header 'bin on the fly' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_COL_ON_THE_FLY_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsCollimateOnTheFly),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsCollimateOnTheFly) == false){
			
			PhgAbort("Unable to set header 'collimate on the fly' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_PHG_IS_DET_ON_THE_FLY_ID,
			sizeof(paramsPtr->H.PhgRunTimeParams.PhgIsDetectOnTheFly),
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgIsDetectOnTheFly) == false){
			
			PhgAbort("Unable to set header 'detect on the fly' parameter", false);
			break;
		}

					
		if (LbHdrStElem(headerHk, HDR_PHG_PHGPARAM__ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgParamFilePath) == false){
			
			PhgAbort("Unable to set header 'param file path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_SUBOBJACTIVITYINDEX_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjActivityIndexFilePath) == false){
			
			PhgAbort("Unable to set header 'activity index path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_SUBOBJACTIVITYTABLE_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjActivityTableFilePath) == false){
			
			PhgAbort("Unable to set header 'activity table path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_SUBOBJACTINDEXTRANS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjActIndexTransFilePath) == false){
			
			PhgAbort("Unable to set header 'activity index translation path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_SUBOBJACTIMG_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjActImgFilePath) == false){
			
			PhgAbort("Unable to set header 'actvity image path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_SUBOBJATTENINDEX_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjAttenIndexFilePath) == false){
			
			PhgAbort("Unable to set header 'attenuation idex path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_SUBOBJATTENTABLE_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjAttenTableFilePath) == false){
			
			PhgAbort("Unable to set header 'attenuation table path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_SUBOBJATTINDEXTRANS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjAttIndexTransFilePath) == false){
			
			PhgAbort("Unable to set header 'attenuation index translation path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_SUBOBJATTIMG_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgSubObjAttImgFilePath) == false){
			
			PhgAbort("Unable to set header 'attenuation image path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_PRODTBLINPUTTABLE_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgProdTblInputTableFilePath) == false){
			
			PhgAbort("Unable to set header 'productivity input path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_PRODTBLOUTPUTTABLE_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgProdTblOutputTableFilePath) == false){
			
			PhgAbort("Unable to set header 'productivity output path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_PHOHSTAT_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgPhoHStatFilePath) == false){
			
			PhgAbort("Unable to set header 'statistics file path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_PHOHFILEHISTORY_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgPhoHFileHistoryFilePath) == false){
			
			PhgAbort("Unable to set header 'history file path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_BINPARAMS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgBinParamsFilePath) == false){
			
			PhgAbort("Unable to set header 'bin parameters file path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_COLLIMATORPARAMS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgCollimatorParamsFilePath) == false){
			
			PhgAbort("Unable to set header 'collimator parameters file path' parameter", false);
			break;
		}
					
		if (LbHdrStElem(headerHk, HDR_PHG_DETECTORPARAMS_ID,
			PATH_LENGTH,
			(void *)&paramsPtr->H.PhgRunTimeParams.PhgDetectorParamsFilePath) == false){
			
			PhgAbort("Unable to set header 'detector parameters file path' parameter", false);
			break;
		}

		/* If collimator was initialized, write the rest of the elements */
		if (paramsPtr->H.PhgRunTimeParams.PhgIsCollimateOnTheFly == true) {
	
			/* Write collimator initialization flag */
			if (LbHdrStElem(headerHk, HDR_COL_IS_INITIALIZED_ID,
					sizeof(paramsPtr->H.ColRunTimeParams.Initialized),
					(void *)&paramsPtr->H.ColRunTimeParams.Initialized) == false){
				
				PhgAbort("Unable to set header 'collimator initialized' parameter", false);
				break;
			}
			if (LbHdrStElem(headerHk, HDR_COL_IS_DO_HISTORY_ID,
				sizeof(paramsPtr->H.ColRunTimeParams.DoHistory),
				(void *)&paramsPtr->H.ColRunTimeParams.DoHistory) == false){
				
				PhgAbort("Unable to set header 'collimator history on/off' parameter", false);
				break;
			}
	
						
			if (LbHdrStElem(headerHk, HDR_COL_HISTORY_ID,
				PATH_LENGTH,
				(void *)&paramsPtr->H.ColRunTimeParams.ColHistoryFilePath) == false){
				
				PhgAbort("Unable to set header 'collimator params file path' parameter", false);
				break;
			}
				
			if (LbHdrStElem(headerHk, HDR_COL_TYPE_ID,
				sizeof(paramsPtr->H.ColRunTimeParams.ColType),
				(void *)&paramsPtr->H.ColRunTimeParams.ColType) == false){
				
				PhgAbort("Unable to set header 'collimator type' parameter", false);
				break;
			}
	
			
			/* Determine which type of collimator was used and convert its fields
			*/
			switch (paramsPtr->H.ColRunTimeParams.ColType) {
			
				case ColEn_simple_pet:
					if (LbHdrStElem(headerHk, HDR_COL_DEPTH_ID,
						sizeof(paramsPtr->H.ColRunTimeParams.SimplePETCol.Depth),
						(void *)&paramsPtr->H.ColRunTimeParams.SimplePETCol.Depth) == false){
				
						PhgAbort("Unable to set header 'simple collimator depth' parameter", false);
						goto FAIL;
					}
					break;
					
				case ColEn_monte_carlo_pet:
					
					/* Just note that the actual collimator description is not stored in the header. Hence,
						we don't have anything to do here.
					*/
					paramsPtr->H.ColRunTimeParams.MCPETCol = 0;
					break;
					
				case ColEn_simple_spect:
					break;
					
				case ColEn_unc_spect:
					
					if (LbHdrStElem(headerHk, HDR_COL_UNC_HOLE_GEOM_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.HoleGeometry),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.HoleGeometry) == false){
				
						PhgAbort("Unable to set header 'UNC collimator hole geometry' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrStElem(headerHk, HDR_COL_UNC_RAD_OF_ROTATION_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.RadiusOfRotation),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.RadiusOfRotation) == false){
				
						PhgAbort("Unable to set header 'UNC collimator radius of rotation' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrStElem(headerHk, HDR_COL_UNC_THICKNESS_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.Thickness),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.Thickness) == false){
						
						PhgAbort("Unable to set header 'UNC collimator thickness' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrStElem(headerHk, HDR_COL_UNC_HOLE_RADIUS_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.HoleRadius),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.HoleRadius) == false){
						
						PhgAbort("Unable to set header 'UNC collimator hole radius' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrStElem(headerHk, HDR_COL_UNC_SEPTAL_THICKNESS_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.SeptalThickness),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.SeptalThickness) == false){
						
						PhgAbort("Unable to set header 'UNC collimator septal thickness' parameter", false);
						goto FAIL;
					}
			
						
					if (LbHdrStElem(headerHk, HDR_COL_UNC_MIN_Z_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.MinZ),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.MinZ) == false){
						
						PhgAbort("Unable to set header 'UNC collimator minimum Z' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrStElem(headerHk, HDR_COL_UNC_MAX_Z_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.MaxZ),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.MaxZ) == false){
						
						PhgAbort("Unable to set header 'UNC collimator maximum Z' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrStElem(headerHk, HDR_COL_UNC_START_ANGLE_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.StartAngle),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.StartAngle) == false){
						
						PhgAbort("Unable to set header 'UNC collimator start angle' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrStElem(headerHk, HDR_COL_UNC_STOP_ANGLE_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.StopAngle),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.StopAngle) == false){
				
						PhgAbort("Unable to set header 'UNC collimator stop angle' parameter", false);
						goto FAIL;
					}
	
						
					if (LbHdrStElem(headerHk, HDR_COL_UNC_SUM_ALL_VIEWS_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.SumAllViews),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.SumAllViews) == false){
						
						PhgAbort("Unable to set header 'UNC collimator sum all views' parameter", false);
						goto FAIL;
					}
	
				
					if (LbHdrStElem(headerHk, HDR_COL_UNC_NUM_VIEWS_ID,
							sizeof(paramsPtr->H.ColRunTimeParams.UNCSPECTCol.NumViews),
							(void *)&paramsPtr->H.ColRunTimeParams.UNCSPECTCol.NumViews) == false){
						
						PhgAbort("Unable to set header 'UNC collimator number of views' parameter", false);
						goto FAIL;
					}
	
						
				default:
					break;
			}
		}
		
		/* If the detector was initialized, write the rest of the fields */
		if (paramsPtr->H.PhgRunTimeParams.PhgIsDetectOnTheFly == true){
						
			/* Write the detector initialized flag */
			if (LbHdrStElem(headerHk, HDR_DET_IS_INITIALIZED_ID,
					sizeof(paramsPtr->H.DetRunTimeParams.Initialized),
					(void *)&paramsPtr->H.DetRunTimeParams.Initialized) == false){
				
				PhgAbort("Unable to set header 'detector initialized' parameter", false);
				break;
			}

			if (LbHdrStElem(headerHk, HDR_DET_IS_DO_HISTORY_ID,
					sizeof(paramsPtr->H.DetRunTimeParams.DoHistory),
					(void *)&paramsPtr->H.DetRunTimeParams.DoHistory) == false){
				
				PhgAbort("Unable to set header 'detector history on/off' parameter", false);
				break;
			}
	
						
			if (LbHdrStElem(headerHk, HDR_DET_HISTORY_ID,
				PATH_LENGTH,
				(void *)&paramsPtr->H.DetRunTimeParams.DetHistoryFilePath) == false){
				
				PhgAbort("Unable to set header 'detector params file path' parameter", false);
				break;
			}
				
			if (LbHdrStElem(headerHk, HDR_DET_TYPE_ID,
				sizeof(paramsPtr->H.DetRunTimeParams.DetectorType),
				(void *)&paramsPtr->H.DetRunTimeParams.DetectorType) == false){
				
				PhgAbort("Unable to set header 'detector type' parameter", false);
				break;
			}
	
			if (LbHdrStElem(headerHk, HDR_DET_IS_FORCED_INTERACTION_ID,
				sizeof(paramsPtr->H.DetRunTimeParams.DoForcedInteraction),
				(void *)&paramsPtr->H.DetRunTimeParams.DoForcedInteraction) == false){
				
				PhgAbort("Unable to set header 'do forced interaction' flag", false);
				break;
			}
	
			if (LbHdrStElem(headerHk, HDR_DET_PHOTON_TIME_FWHM_ID,
				sizeof(paramsPtr->H.DetRunTimeParams.PhotonTimeFWHM),
				(void *)&paramsPtr->H.DetRunTimeParams.PhotonTimeFWHM) == false){
				
				PhgAbort("Unable to set header 'photon detection time fwhm' parameter", false);
				break;
			}
	
			if (LbHdrStElem(headerHk, HDR_DET_COINC_TIMING_WINDOW_ID,
				sizeof(paramsPtr->H.DetRunTimeParams.CoincidenceTimingWindowNS),
				(void *)&paramsPtr->H.DetRunTimeParams.CoincidenceTimingWindowNS) == false){
				
				PhgAbort("Unable to set header 'coincidence timing window' parameter", false);
				break;
			}
	
			if (LbHdrStElem(headerHk, HDR_DET_TRIPLES_METHOD_ID,
				sizeof(paramsPtr->H.DetRunTimeParams.TriplesMethod),
				(void *)&paramsPtr->H.DetRunTimeParams.TriplesMethod) == false){
				
				PhgAbort("Unable to set header 'triples method' parameter", false);
				break;
			}
	
			if (LbHdrStElem(headerHk, HDR_DET_RANDOMS_HISTORY_FILE_ID,
				PATH_LENGTH,
				(void *)&paramsPtr->H.DetRunTimeParams.DetRandomsHistoryFilePath) == false){
				
				PhgAbort("Unable to set header 'randoms history file path' parameter", false);
				break;
			}
	
			
			/* Determine which type of detector was used and convert its fields
			*/
			switch (paramsPtr->H.DetRunTimeParams.DetectorType) {
			
				case DetEn_simple_pet:
					if (LbHdrStElem(headerHk, HDR_DET_ENERGY_RESOLUTION_PER_ID,
							sizeof(paramsPtr->H.DetRunTimeParams.EnergyResolutionPercentage),
							(void *)&paramsPtr->H.DetRunTimeParams.EnergyResolutionPercentage) == false){
					
					PhgAbort("Unable to set header 'simple detector energy resolution' parameter", false);
					goto FAIL;
				}
	
				
					if (LbHdrStElem(headerHk, HDR_DET_REFERENCE_ENERGY_ID,
							sizeof(paramsPtr->H.DetRunTimeParams.ReferenceEnergy),
							(void *)&paramsPtr->H.DetRunTimeParams.ReferenceEnergy) == false){
					
					PhgAbort("Unable to set header 'simple detector reference energy' parameter", false);
					goto FAIL;
				}
	
				
					break;
					
				case DetEn_simple_spect:
					break;
					
				case DetEn_unc_spect:
					break;
					
				default:
					break;
			}
		}
		/************* GENERAL COMPUTED PARAMETERS ************/
		if (LbHdrStElem(headerHk, HDR_CYL_TARGET_RADIUS_ID,
				sizeof(paramsPtr->H.TargetCylinder.radius),
				(void *)&paramsPtr->H.TargetCylinder.radius) == false){
			
			PhgAbort("Unable to set header 'target cylinder radius' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_TARGET_MIN_Z_ID,
				sizeof(paramsPtr->H.TargetCylinder.zMin),
				(void *)&paramsPtr->H.TargetCylinder.zMin) == false){
			
			PhgAbort("Unable to set header 'target cylinder z min' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_TARGET_MAX_Z_ID,
				sizeof(paramsPtr->H.TargetCylinder.zMax),
				(void *)&paramsPtr->H.TargetCylinder.zMax) == false){
			
			PhgAbort("Unable to set header 'target cylinder z max' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_TARGET_CENTER_X_ID,
				sizeof(paramsPtr->H.TargetCylinder.centerX),
				(void *)&paramsPtr->H.TargetCylinder.centerX) == false){
			
			PhgAbort("Unable to set header 'target cylinder center X' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_TARGET_CENTER_Y_ID,
				sizeof(paramsPtr->H.TargetCylinder.centerY),
				(void *)&paramsPtr->H.TargetCylinder.centerY) == false){
			
			PhgAbort("Unable to set header 'target cylinder center Y' parameter", false);
			break;
		}

		
		if (LbHdrStElem(headerHk, HDR_CYL_CRITZONE_RADIUS_ID,
				sizeof(paramsPtr->H.CriticalZone.radius),
				(void *)&paramsPtr->H.CriticalZone.radius) == false){
			
			PhgAbort("Unable to set header 'critical zone radius' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_CRITZONE_MIN_Z_ID,
				sizeof(paramsPtr->H.CriticalZone.zMin),
				(void *)&paramsPtr->H.CriticalZone.zMin) == false){
			
			PhgAbort("Unable to set header 'critical zone z min' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_CRITZONE_MAX_Z_ID,
				sizeof(paramsPtr->H.CriticalZone.zMax),
				(void *)&paramsPtr->H.CriticalZone.zMax) == false){
			
			PhgAbort("Unable to set header 'critical zone z max' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_CRITZONE_CENTER_X_ID,
				sizeof(paramsPtr->H.CriticalZone.centerX),
				(void *)&paramsPtr->H.CriticalZone.centerX) == false){
			
			PhgAbort("Unable to set header 'critical zone center X' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_CRITZONE_CENTER_Y_ID,
				sizeof(paramsPtr->H.CriticalZone.centerY),
				(void *)&paramsPtr->H.CriticalZone.centerY) == false){
			
			PhgAbort("Unable to set header 'critical zone center Y' parameter", false);
			break;
		}

			
					
		if (LbHdrStElem(headerHk, HDR_CYL_OBJCYL_RADIUS_ID,
				sizeof(paramsPtr->H.ObjectCylinder.radius),
				(void *)&paramsPtr->H.ObjectCylinder.radius) == false){
			
			PhgAbort("Unable to set header 'critical zone radius' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_OBJCYL_MIN_Z_ID,
				sizeof(paramsPtr->H.ObjectCylinder.zMin),
				(void *)&paramsPtr->H.ObjectCylinder.zMin) == false){
			
			PhgAbort("Unable to set header 'object cylinder z min' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_OBJCYL_MAX_Z_ID,
				sizeof(paramsPtr->H.ObjectCylinder.zMax),
				(void *)&paramsPtr->H.ObjectCylinder.zMax) == false){
			
			PhgAbort("Unable to set header 'object cylinder z max' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_OBJCYL_CENTER_X_ID,
				sizeof(paramsPtr->H.ObjectCylinder.centerX),
				(void *)&paramsPtr->H.ObjectCylinder.centerX) == false){
			
			PhgAbort("Unable to set header 'object cylinder center X' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_OBJCYL_CENTER_Y_ID,
				sizeof(paramsPtr->H.ObjectCylinder.centerY),
				(void *)&paramsPtr->H.ObjectCylinder.centerY) == false){
			
			PhgAbort("Unable to set header 'object cylinder center Y' parameter", false);
			break;
		}

			
	
		if (LbHdrStElem(headerHk, HDR_CYL_LIMITCYL_RADIUS_ID,
				sizeof(paramsPtr->H.LimitCylinder.radius),
				(void *)&paramsPtr->H.LimitCylinder.radius) == false){
			
			PhgAbort("Unable to set header 'object cylinder radius' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_LIMITCYL_MIN_Z_ID,
				sizeof(paramsPtr->H.LimitCylinder.zMin),
				(void *)&paramsPtr->H.LimitCylinder.zMin) == false){
			
			PhgAbort("Unable to set header 'limit cylinder z min' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_LIMITCYL_MAX_Z_ID,
				sizeof(paramsPtr->H.LimitCylinder.zMax),
				(void *)&paramsPtr->H.LimitCylinder.zMax) == false){
			
			PhgAbort("Unable to set header 'limit cylinder z max' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_LIMITCYL_CENTER_X_ID,
				sizeof(paramsPtr->H.LimitCylinder.centerX),
				(void *)&paramsPtr->H.LimitCylinder.centerX) == false){
			
			PhgAbort("Unable to set header 'limit cylinder center X' parameter", false);
			break;
		}

			
		if (LbHdrStElem(headerHk, HDR_CYL_LIMITCYL_CENTER_Y_ID,
				sizeof(paramsPtr->H.LimitCylinder.centerY),
				(void *)&paramsPtr->H.LimitCylinder.centerY) == false){
			
			PhgAbort("Unable to set header 'limit cylinder center Y' parameter", false);
			break;
		}

	
		if (LbHdrStElem(headerHk, HDR_BIN1_NUM_SIMULATIONS_ID,
				sizeof(paramsPtr->H.NumSimulations),
				(void *)&paramsPtr->H.NumSimulations) == false){
			
			PhgAbort("Unable to set header 'number of simulations' parameter", false);
			break;
		}

	
		if (LbHdrStElem(headerHk, HDR_BIN1_EVENTS_TO_SIMULATION_ID,
				sizeof(paramsPtr->H.SumEventsToSimulate),
				(void *)&paramsPtr->H.SumEventsToSimulate) == false){
			
			PhgAbort("Unable to set header 'sum of events to simulate' parameter", false);
			break;
		}

	
		if (LbHdrStElem(headerHk, HDR_BIN1_NUM_PHOTONS_ID,
				sizeof(paramsPtr->H.NumPhotons),
				(void *)&paramsPtr->H.NumPhotons) == false){
			
			PhgAbort("Unable to set header 'number of photons' parameter", false);
			break;
		}

	
		if (LbHdrStElem(headerHk, HDR_BIN1_NUM_DECAYS_ID,
				sizeof(paramsPtr->H.NumDecays),
				(void *)&paramsPtr->H.NumDecays) == false){
			
			PhgAbort("Unable to set header 'number of decays' parameter", false);
			break;
		}


		/******** BIN PARAMETERS **********/
		/* Set the binning parameters if we are binning */
		if ((paramsPtr->H.HdrKind == PhoHFileEn_BIN_CT) ||
				(paramsPtr->H.HdrKind == PhoHFileEn_BIN_WT) ||
				(paramsPtr->H.HdrKind == PhoHFileEn_BIN_WTSQ)) {
				
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_Z_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numZBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numZBins) == false){
				
				PhgAbort("Unable to set header 'number of z bins' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_PA_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numPABins),
					(void *)&paramsPtr->H.BinRunTimeParams.numPABins) == false){
				
				PhgAbort("Unable to set header 'number of PA bins' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_TD_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numTDBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numTDBins) == false){
				
				PhgAbort("Unable to set header 'number of TD bins' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_AA_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numAABins),
					(void *)&paramsPtr->H.BinRunTimeParams.numAABins) == false){
				
				PhgAbort("Unable to set header 'number of AA bins' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_E1_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numE1Bins),
					(void *)&paramsPtr->H.BinRunTimeParams.numE1Bins) == false){
				
				PhgAbort("Unable to set header 'number of E1 bins' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_E2_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numE2Bins),
					(void *)&paramsPtr->H.BinRunTimeParams.numE2Bins) == false){
				
				PhgAbort("Unable to set header 'number of E2 bins' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_S1_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numS1Bins),
					(void *)&paramsPtr->H.BinRunTimeParams.numS1Bins) == false){
				
				PhgAbort("Unable to set header 'number of S1 bins' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_S2_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numS2Bins),
					(void *)&paramsPtr->H.BinRunTimeParams.numS2Bins) == false){
				
				PhgAbort("Unable to set header 'number of S2 bins' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_IMAGE_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numImageBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numImageBins) == false){
				
				PhgAbort("Unable to set header 'number of image bins' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_SCATTER_PARAM_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatterRandomParam),
					(void *)&paramsPtr->H.BinRunTimeParams.scatterRandomParam) == false){
				
				PhgAbort("Unable to set header 'scatter binning' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_MIN_Z_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minZ),
				 	(void *)&paramsPtr->H.BinRunTimeParams.minZ) == false){
				
				PhgAbort("Unable to set header 'minimum Z' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_MAX_Z_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxZ),
					(void *)&paramsPtr->H.BinRunTimeParams.maxZ) == false){
				
				PhgAbort("Unable to set header 'maximum Z' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_MIN_PA_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minPA),
					(void *)&paramsPtr->H.BinRunTimeParams.minPA) == false){
				
				PhgAbort("Unable to set header 'minimum PA' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_MAX_PA_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxPA),
					(void *)&paramsPtr->H.BinRunTimeParams.maxPA) == false){
				
				PhgAbort("Unable to set header 'maximum PA' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_MIN_TD_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minTD),
					(void *)&paramsPtr->H.BinRunTimeParams.minTD) == false){
				
				PhgAbort("Unable to set header 'minimum TD' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_MAX_TD_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxTD),
					(void *)&paramsPtr->H.BinRunTimeParams.maxTD) == false){
				
				PhgAbort("Unable to set header 'maximum TD' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_MIN_AA_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minAA),
					(void *)&paramsPtr->H.BinRunTimeParams.minAA) == false){
				
				PhgAbort("Unable to set header 'minimum AA' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_MAX_AA_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxAA),
					(void *)&paramsPtr->H.BinRunTimeParams.maxAA) == false){
				
				PhgAbort("Unable to set header 'maximum AA' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_MIN_ENERGY_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minE),
					(void *)&paramsPtr->H.BinRunTimeParams.minE) == false){
				
				PhgAbort("Unable to set header 'minimum energy' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_MAX_ENERGY_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxE),
					(void *)&paramsPtr->H.BinRunTimeParams.maxE) == false){
				
				PhgAbort("Unable to set header 'maximum energy' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_MIN_SCAT_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minS),
					(void *)&paramsPtr->H.BinRunTimeParams.minS) == false){
				
				PhgAbort("Unable to set header 'minimum scatters' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_MAX_SCAT_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxS),
					(void *)&paramsPtr->H.BinRunTimeParams.maxS) == false){
				
				PhgAbort("Unable to set header 'maximum scatters' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_ADD_TO_EXISTING_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.addToExistingImg),
					(void *)&paramsPtr->H.BinRunTimeParams.addToExistingImg) == false){
				
				PhgAbort("Unable to set header 'add to existing images' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_DO_COUNTS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.doCounts),
					(void *)&paramsPtr->H.BinRunTimeParams.doCounts) == false){
				
				PhgAbort("Unable to set header 'do counts' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_DO_WEIGHTS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.doWeights),
					(void *)&paramsPtr->H.BinRunTimeParams.doWeights) == false){
				
				PhgAbort("Unable to set header 'do weights' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_DO_WEIGHTS_SQU_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.doWeightsSquared),
					(void *)&paramsPtr->H.BinRunTimeParams.doWeightsSquared) == false){
				
				PhgAbort("Unable to set header 'do weights squared' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_Z_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.zRange),
					(void *)&paramsPtr->H.BinRunTimeParams.zRange) == false){
				
				PhgAbort("Unable to set header 'z range' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_E_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.eRange),
					(void *)&paramsPtr->H.BinRunTimeParams.eRange) == false){
				
				PhgAbort("Unable to set header 'energy range' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_S_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.sRange),
					(void *)&paramsPtr->H.BinRunTimeParams.sRange) == false){
				
				PhgAbort("Unable to set header 'scatter range' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_TD_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tdRange),
					(void *)&paramsPtr->H.BinRunTimeParams.tdRange) == false){
				
				PhgAbort("Unable to set header 'td range' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_SCATTER2CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatter2CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.scatter2CIsize) == false){
				
				PhgAbort("Unable to set header 'scatter2 count size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_SCATTER2WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatter2WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.scatter2WIsize) == false){
				
				PhgAbort("Unable to set header 'scatter2 weight size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_SCATTER2WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatter2WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.scatter2WISsize) == false){
				
				PhgAbort("Unable to set header 'scatter2 weight squared size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_SCATTER1CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatter1CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.scatter1CIsize) == false){
				
				PhgAbort("Unable to set header 'scatter1 count size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_SCATTER1WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatter1WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.scatter1WIsize) == false){
				
				PhgAbort("Unable to set header 'scatter1 weight size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_SCATTER1WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.scatter1WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.scatter1WISsize) == false){
				
				PhgAbort("Unable to set header 'scatter1 weight squared size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_ENERGY2CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.energy2CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.energy2CIsize) == false){
				
				PhgAbort("Unable to set header 'energy2 count size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_ENERGY2WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.energy2WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.energy2WIsize) == false){
				
				PhgAbort("Unable to set header 'energy2 weight size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_ENERGY2WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.energy2WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.energy2WISsize) == false){
				
				PhgAbort("Unable to set header 'energy2 weight squared size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_ENERGY1CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.energy1CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.energy1CIsize) == false){
				
				PhgAbort("Unable to set header 'energy1 count size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_ENERGY1WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.energy1WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.energy1WIsize) == false){
				
				PhgAbort("Unable to set header 'energy1 weight size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_ENERGY1WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.energy1WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.energy1WISsize) == false){
				
				PhgAbort("Unable to set header 'energy1 weight squared size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_AACI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.aaCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.aaCIsize) == false){
				
				PhgAbort("Unable to set header 'aa count size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_AAWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.aaWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.aaWIsize) == false){
				
				PhgAbort("Unable to set header 'aa weight size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_AAWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.aaWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.aaWISsize) == false){
				
				PhgAbort("Unable to set header 'aa weight squared size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_TDCI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tdCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.tdCIsize) == false){
				
				PhgAbort("Unable to set header 'td count size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_TDWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tdWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.tdWIsize) == false){
				
				PhgAbort("Unable to set header 'td weight size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_TDWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tdWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.tdWISsize) == false){
				
				PhgAbort("Unable to set header 'td weight squared size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_PACI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.paCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.paCIsize) == false){
				
				PhgAbort("Unable to set header 'pa count size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_PAWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.paWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.paWIsize) == false){
				
				PhgAbort("Unable to set header 'pa weight size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_PAWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.paWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.paWISsize) == false){
				
				PhgAbort("Unable to set header 'pa weight squared size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_Z2CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.z2CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.z2CIsize) == false){
				
				PhgAbort("Unable to set header 'z2 count size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_Z2WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.z2WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.z2WIsize) == false){
				
				PhgAbort("Unable to set header 'z2 weight size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_Z2WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.z2WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.z2WISsize) == false){
				
				PhgAbort("Unable to set header 'z2 weight squared size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_Z1CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.z1CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.z1CIsize) == false){
				
				PhgAbort("Unable to set header 'z1 count size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_Z1WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.z1WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.z1WIsize) == false){
				
				PhgAbort("Unable to set header 'z1 weight size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_Z1WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.z1WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.z1WISsize) == false){
				
				PhgAbort("Unable to set header 'z1 weight squared size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_COUNT_IMG_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.countImageSize),
					(void *)&paramsPtr->H.BinRunTimeParams.countImageSize) == false){
				
				PhgAbort("Unable to set header 'count image size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_WEIGHT_IMG_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.weightImageSize),
					(void *)&paramsPtr->H.BinRunTimeParams.weightImageSize) == false){
				
				PhgAbort("Unable to set header 'weight image size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_WEIGHT_SQU_IMG_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.weightSquImageSize),
					(void *)&paramsPtr->H.BinRunTimeParams.weightSquImageSize) == false){
				
				PhgAbort("Unable to set header 'weight squared image size' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_WEIGHT_IMG_TYPE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.weight_image_type),
					(void *)&paramsPtr->H.BinRunTimeParams.weight_image_type) == false){
				
				PhgAbort("Unable to set header 'weight image type' parameter", false);
				break;
			}
	
		
			if (LbHdrStElem(headerHk, HDR_BIN_COUNT_IMG_TYPE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.count_image_type),
					(void *)&paramsPtr->H.BinRunTimeParams.count_image_type) == false){
				
				PhgAbort("Unable to set header 'count image type' parameter", false);
				break;
			}

			if (LbHdrStElem(headerHk, HDR_BIN_DIMENSION_ORDER_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.PhgBinDimensions),
					(void *)&paramsPtr->H.BinRunTimeParams.PhgBinDimensions) == false){
				
				PhgAbort("Unable to set header 'binning dimensions ordering' parameter", false);
				break;
			}

			/* THESE ARE NEW PARAMETERS ADDED TO THE HEADER BEFORE THE PUBLIC RELEASE
			*/
	
			/* Set the size of the PHI count dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_PHICI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.phiCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.phiCIsize) == false){
				
				PhgAbort("Unable to set header 'size of the PHI count dimension' parameter", false);
				break;
			}
	
			/* Set the size of the PHI weight dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_PHIWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.phiWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.phiWIsize) == false){
				
				PhgAbort("Unable to set header 'size of the PHI weight dimension' parameter", false);
				break;
			}
	
			/* Set the size of the PHI weight squared dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_PHIWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.phiWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.phiWISsize) == false){
				
				PhgAbort("Unable to set header 'size of the PHI weight squared dimension' parameter", false);
				break;
			}
	
			/* Set the size of the theta count dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_THETACI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.thetaCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.thetaCIsize) == false){
				
				PhgAbort("Unable to set header 'size of the theta count dimension' parameter", false);
				break;
			}
	
			/* Set the size of the theta weight dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_THETAWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.thetaWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.thetaWIsize) == false){
				
				PhgAbort("Unable to set header 'size of the theta weight dimension' parameter", false);
				break;
			}
	
			/* Set the size of the theta weight squared dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_THETAWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.thetaWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.thetaWISsize) == false){
				
				PhgAbort("Unable to set header 'size of the theta weight squared dimension' parameter", false);
				break;
			}
	
			/* Set the size of the Xr dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_XRCI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.xrCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.xrCIsize) == false){
				
				PhgAbort("Unable to set header 'size of the Xr dimension' parameter", false);
				break;
			}
	
			/* Set the size of the Xr weight dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_XRWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.xrWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.xrWIsize) == false){
				
				PhgAbort("Unable to set header 'size of the Xr weight dimension' parameter", false);
				break;
			}
	
			/* Set the size of the Xr weight squared dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_XRWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.xrWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.xrWISsize) == false){
				
				PhgAbort("Unable to set header 'size of the Xr weight squared dimension' parameter", false);
				break;
			}
	
			/* Set the size of the Yr count dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_YRCI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.yrCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.yrCIsize) == false){
				
				PhgAbort("Unable to set header 'size of the Yr count dimension' parameter", false);
				break;
			}
	
			/* Set the size of the Yr weight dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_YRWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.yrWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.yrWIsize) == false){
				
				PhgAbort("Unable to set header 'size of the Yr weight dimension' parameter", false);
				break;
			}
	
			/* Set the size of the Yr weight squared dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_YRWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.yrWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.yrWISsize) == false){
				
				PhgAbort("Unable to set header 'size of the Yr weight squared dimension' parameter", false);
				break;
			}
	
			/* Set the number of PHI bins */
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_PHI_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numPHIBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numPHIBins) == false){
				
				PhgAbort("Unable to set header 'number of PHI bins' parameter", false);
				break;
			}
	
			/* Set the number of theta bins */
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_THETA_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numThetaBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numThetaBins) == false){
				
				PhgAbort("Unable to set header 'number of theta bins' parameter", false);
				break;
			}
	
			/* Set the number of Xr bins */
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_XR_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numXRBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numXRBins) == false){
				
				PhgAbort("Unable to set header 'number of Xr bins' parameter", false);
				break;
			}
	
			/* Set the number of Yr bins */
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_YR_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numYRBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numYRBins) == false){
				
				PhgAbort("Unable to set header 'number of Yr bins' parameter", false);
				break;
			}
	
			/* Set the min theta */
			if (LbHdrStElem(headerHk, HDR_BIN_MIN_THETA_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minTheta),
					(void *)&paramsPtr->H.BinRunTimeParams.minTheta) == false){
				
				PhgAbort("Unable to set header 'min theta' parameter", false);
				break;
			}
	
			/* Get max theta */
			if (LbHdrStElem(headerHk, HDR_BIN_MAX_THETA_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxTheta),
					(void *)&paramsPtr->H.BinRunTimeParams.maxTheta) == false){
				
				PhgAbort("Unable to set header 'max theta' parameter", false);
				break;
			}
	
			/* Get min PHI */
			if (LbHdrStElem(headerHk, HDR_BIN_MIN_PHI_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minPhi),
					(void *)&paramsPtr->H.BinRunTimeParams.minPhi) == false){
				
				PhgAbort("Unable to set header 'min PHI' parameter", false);
				break;
			}
	
			/* Get max PHI */
			if (LbHdrStElem(headerHk, HDR_BIN_MAX_PHI_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxPhi),
					(void *)&paramsPtr->H.BinRunTimeParams.maxPhi) == false){
				
				PhgAbort("Unable to set header 'max PHI' parameter", false);
				break;
			}
	
			/* Get min Xr */
			if (LbHdrStElem(headerHk, HDR_BIN_MIN_XR_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minXR),
					(void *)&paramsPtr->H.BinRunTimeParams.minXR) == false){
				
				PhgAbort("Unable to set header 'min Xr' parameter", false);
				break;
			}
	
			/* Get max Xr */
			if (LbHdrStElem(headerHk, HDR_BIN_MAX_XR_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxXR),
					(void *)&paramsPtr->H.BinRunTimeParams.maxXR) == false){
				
				PhgAbort("Unable to set header 'max Xr' parameter", false);
				break;
			}
	
			/* Get min Yr */
			if (LbHdrStElem(headerHk, HDR_BIN_MIN_YR_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minYR),
					(void *)&paramsPtr->H.BinRunTimeParams.minYR) == false){
				
				PhgAbort("Unable to set header 'min Yr' parameter", false);
				break;
			}
	
			/* Get max Yr */
			if (LbHdrStElem(headerHk, HDR_BIN_MAX_YR_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxYR),
					(void *)&paramsPtr->H.BinRunTimeParams.maxYR) == false){
				
				PhgAbort("Unable to set header 'max Yr' parameter", false);
				break;
			}
	
			/* Get SSRB flag */
			if (LbHdrStElem(headerHk, HDR_BIN_DOSSRB_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.doSSRB),
					(void *)&paramsPtr->H.BinRunTimeParams.doSSRB) == false){
				
				PhgAbort("Unable to set header 'SSRB flag' parameter", false);
				break;
			}
	
			/* Get MSRB flag */
			if (LbHdrStElem(headerHk, HDR_BIN_DOMSRB_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.doMSRB),
					(void *)&paramsPtr->H.BinRunTimeParams.doMSRB) == false){
				
				PhgAbort("Unable to set header 'MSRB flag' parameter", false);
				break;
			}
	
			/* Set the image radius */
			if (LbHdrStElem(headerHk, HDR_BIN_IMG_RADIUS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.image_radius),
					(void *)&paramsPtr->H.BinRunTimeParams.image_radius) == false){
				
				PhgAbort("Unable to set header 'image radius' parameter", false);
				break;
			}
	
			/* Set the detector radius */
			if (LbHdrStElem(headerHk, HDR_BIN_DET_RADIUS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.detector_radius),
					(void *)&paramsPtr->H.BinRunTimeParams.detector_radius) == false){
				
				PhgAbort("Unable to set header 'detector radius' parameter", false);
				break;
			}
	
			/* Set the weight image path */
			if (LbHdrStElem(headerHk, HDR_BIN_WTIMG_PATH_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.weightImgFilePath),
					(void *)&paramsPtr->H.BinRunTimeParams.weightImgFilePath) == false){
				
				PhgAbort("Unable to set header 'weight image path' parameter", false);
				break;
			}
	
			/* Set the weight squared image path */
			if (LbHdrStElem(headerHk, HDR_BIN_WTSIMG_PATH_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.weightSquImgFilePath),
					(void *)&paramsPtr->H.BinRunTimeParams.weightSquImgFilePath) == false){
				
				PhgAbort("Unable to set header 'weight squared image path' parameter", false);
				break;
			}
	
			/* Set the count image path */
			if (LbHdrStElem(headerHk, HDR_BIN_CTIMG_PATH_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.countImgFilePath),
					(void *)&paramsPtr->H.BinRunTimeParams.countImgFilePath) == false){
				
				PhgAbort("Unable to set header 'count image path' parameter", false);
				break;
			}
	
			/* Set the theta range */
			if (LbHdrStElem(headerHk, HDR_BIN_THETA_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.thetaRange),
					(void *)&paramsPtr->H.BinRunTimeParams.thetaRange) == false){
				
				PhgAbort("Unable to set header 'theta range' parameter", false);
				break;
			}
	
			/* Set the phi range */
			if (LbHdrStElem(headerHk, HDR_BIN_PHI_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.phiRange),
					(void *)&paramsPtr->H.BinRunTimeParams.phiRange) == false){
				
				PhgAbort("Unable to set header 'phi range' parameter", false);
				break;
			}
	
			/* Set the Xr range */
			if (LbHdrStElem(headerHk, HDR_BIN_XR_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.xrRange),
					(void *)&paramsPtr->H.BinRunTimeParams.xrRange) == false){
				
				PhgAbort("Unable to set header 'Xr range' parameter", false);
				break;
			}
	
			/* Set the Yr range */
			if (LbHdrStElem(headerHk, HDR_BIN_YR_RANGE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.yrRange),
					(void *)&paramsPtr->H.BinRunTimeParams.yrRange) == false){
				
				PhgAbort("Unable to set header 'Yr range' parameter", false);
				break;
			}

	
			/* Set the sum weights  */
			if (LbHdrStElem(headerHk, HDR_BIN_SUM_WEIGHTS_ID,
					sizeof(paramsPtr->H.weightSum),
					(void *)&paramsPtr->H.weightSum) == false){
				
				PhgAbort("Unable to set header 'sum of weights' parameter", false);
				break;
			}
			
			/* Set the sum weights squared  */
			if (LbHdrStElem(headerHk, HDR_BIN_SUM_WEIGHTS_SQ_ID,
					sizeof(paramsPtr->H.weightSquSum),
					(void *)&paramsPtr->H.weightSquSum) == false){
				
				PhgAbort("Unable to set header 'sum of weights squared' parameter", false);
				break;
			}
			

			/* THESE ARE NEW PARAMETERS ADDED TO THE HEADER AFTER PUBLIC RELEASE
			*/

			/* Set the timesorted indicator  */
			if (LbHdrStElem(headerHk, HDR_HISTORY_FILE_IS_SORTED_ID,
					sizeof(paramsPtr->H.isTimeSorted),
					(void *)&paramsPtr->H.isTimeSorted) == false){
				
				PhgAbort("Unable to set header 'timesorted indicator' parameter", false);
				break;
			}
			
			/* Set the randoms added indicator  */
			if (LbHdrStElem(headerHk, HDR_HISTORY_FILE_IS_RANDOMS_ADDED_ID,
					sizeof(paramsPtr->H.isRandomsAdded),
					(void *)&paramsPtr->H.isRandomsAdded) == false){
				
				PhgAbort("Unable to set header 'randoms added indicator' parameter", false);
				break;
			}

			if (LbHdrStElem(headerHk, HDR_BIN_NUM_TOF_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numTOFBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numTOFBins) == false){
				
				PhgAbort("Unable to set header 'number of TOF bins' parameter", false);
				break;
			}
	
			/* Set the min TOF */
			if (LbHdrStElem(headerHk, HDR_BIN_MIN_TOF_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.minTOF),
					(void *)&paramsPtr->H.BinRunTimeParams.minTOF) == false){
				
				PhgAbort("Unable to set header 'min TOF' parameter", false);
				break;
			}
	
			/* Get max TOF */
			if (LbHdrStElem(headerHk, HDR_BIN_MAX_TOF_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.maxTOF),
					(void *)&paramsPtr->H.BinRunTimeParams.maxTOF) == false){
				
				PhgAbort("Unable to set header 'max TOF' parameter", false);
				break;
			}
	
			/* Set the size of the TOF count dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_TOFCI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tofCIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.tofCIsize) == false){
				
				PhgAbort("Unable to set header 'size of the TOF count dimension' parameter", false);
				break;
			}
	
			/* Set the size of the TOF weight dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_TOFWI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tofWIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.tofWIsize) == false){
				
				PhgAbort("Unable to set header 'size of the TOF weight dimension' parameter", false);
				break;
			}
	
			/* Set the size of the TOF weight squared dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_TOFWIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.tofWISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.tofWISsize) == false){
				
				PhgAbort("Unable to set header 'size of the TOF weight squared dimension' parameter", false);
				break;
			}
	
			/* Get crystal binning flag */
			if (LbHdrStElem(headerHk, HDR_BIN_DO_CRYSTAL_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.doCrystal),
					(void *)&paramsPtr->H.BinRunTimeParams.doCrystal) == false){
				
				PhgAbort("Unable to set header 'binning-by-crystal flag' parameter", false);
				break;
			}
	
			/* Set the number of crystal bins */
			if (LbHdrStElem(headerHk, HDR_BIN_NUM_CRYSTAL_BINS_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.numCrystalBins),
					(void *)&paramsPtr->H.BinRunTimeParams.numCrystalBins) == false){
				
				PhgAbort("Unable to set header 'number of crystal bins' parameter", false);
				break;
			}
	
			/* Set the size of the crystal count dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_CRYSTAL1CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.crystal1CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.crystal1CIsize) == false){
				
				PhgAbort("Unable to set header 'size of the crystal1 count dimension' parameter", false);
				break;
			}
	
			/* Set the size of the crystal1 weight dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_CRYSTAL1WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.crystal1WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.crystal1WIsize) == false){
				
				PhgAbort("Unable to set header 'size of the crystal1 weight dimension' parameter", false);
				break;
			}
	
			/* Set the size of the crystal1 weight squared dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_CRYSTAL1WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.crystal1WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.crystal1WISsize) == false){
				
				PhgAbort("Unable to set header 'size of the crystal1 weight squared dimension' parameter", false);
				break;
			}
	
			/* Set the size of the crystal2 count dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_CRYSTAL2CI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.crystal2CIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.crystal2CIsize) == false){
				
				PhgAbort("Unable to set header 'size of the crystal2 count dimension' parameter", false);
				break;
			}
	
			/* Set the size of the crystal2 weight dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_CRYSTAL2WI_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.crystal2WIsize),
					(void *)&paramsPtr->H.BinRunTimeParams.crystal2WIsize) == false){
				
				PhgAbort("Unable to set header 'size of the crystal2 weight dimension' parameter", false);
				break;
			}
	
			/* Set the size of the crystal2 weight squared dimension */
			if (LbHdrStElem(headerHk, HDR_BIN_CRYSTAL2WIS_SIZE_ID,
					sizeof(paramsPtr->H.BinRunTimeParams.crystal2WISsize),
					(void *)&paramsPtr->H.BinRunTimeParams.crystal2WISsize) == false){
				
				PhgAbort("Unable to set header 'size of the crystal2 weight squared dimension' parameter", false);
				break;
			}
	
		
	
			
		}

		okay = true;
		FAIL:;
	} while (false);

	/* Return the status */
	return (okay);
}

/**********************
*	PhgHdrStAttenuationCorrected
*
*	Purpose:	This routine puts an attenuation corrected header field
*				into the header. We assume the value is true if you are
*				setting this.
*
*	Arguments:
*			LbHdrHkTy		*headerHk	- The newly created header structure.
*
*	Result:	TRUE unless an error occurs.
***********************/
Boolean  PhgHdrStAttenuationCorrected(LbHdrHkTy *headerHk)
{
	Boolean	okay = false;			/* Process Loop */	
	Boolean	attCorrected = true;	/* Data for header field */
	
	do { /* Process Loop */
		
		/* Set the attenuation corrected flag */
		if (LbHdrStElem(headerHk, HDR_ATT_HAS_BEEN_CORRECTED_ID, sizeof(Boolean),
				(void *)&(attCorrected)) != true) {
			
			break;
		}
		
		okay = true;
	} while (false);
		
	return (okay);
}

/**********************
*	PhgHdrGtAttenuationCorrected
*
*	Purpose:	This routine checks for an attenuation corrected flag
*				in the header. If one is not found, we assume that the
*				file has not been attenuation corrected.
*
*	Arguments:
*			LbHdrHkTy	*headerHk		- The newly created header structure.
*			Boolean		*attCorrected	- The flag indicating correction.
*	Result:	TRUE unless an error occurs.
***********************/
Boolean  PhgHdrGtAttenuationCorrected(LbHdrHkTy *headerHk, Boolean *attCorrected)
{
	Boolean	okay = false;			/* Process Loop */	
	Boolean cleared = false;		/* Error handling flag */
	
	do { /* Process Loop */
		
		/* If this fails we check to verify it is because the item is not found */
		if (LbHdrGtElem(headerHk, HDR_ATT_HAS_BEEN_CORRECTED_ID, sizeof(Boolean),
				(void *)attCorrected) != true) {
			
			/* See if it wasn't retrieved because it isn't in the file.
				The alternative is that there was some failure in processing
				that should be reported.
			*/
			ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
			
			if (cleared == true)
				*attCorrected = false;
			else
				break;
		}
		
		okay = true;
	} while (false);
		
	return (okay);
}
