/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1994-2013 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			DetParams.c
*     Revision Number:		1.6
*     Date last revised:	23 July 2013
*     Programmer:			Steven Vannoy
*     Date Originated:		11 July 1994
*
*     Module Overview:	The parameter module for the detector.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*			DetLookupRunTimeParamLabel
*			DetGetRunTimeParams
*
*     Global variables defined:
*
*
*	  Global macros defined:
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		Sept 2009 - March 2010
*
*			Revision description:	Changed to allow multi-element blocks. 
*						Added computation of the block's position-in-tomograph
*						to the block record.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		July 2007
*
*			Revision description:	Changed input of block positioning 
*						within rings to cylindrical coords.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		August 2006
*
*			Revision description:	Added block detector parameters and
*						functions.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:	Added time-of-flight parameters.
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
*********************************************************************************/
#define DET_PARAMS

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
#include "LbHeader.h"
#include "LbInterface.h"
#include "Lb2DGeometry.h"

#include "Photon.h"
#include "PhgParams.h"
#include "PhgMath.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "CylPos.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "PhoHFile.h"
#include "phg.h"
#include "PhgBin.h"
#include "ProdTbl.h"
#include "SubObj.h"
#include "Detector.h"
#include "DetBlock.h"


/* LOCAL CONSTANTS */
/* LOCAL TYPES */
typedef char labelTy[LBPF_LABEL_LEN];

/* LOCAL GLOBALS */
static char			detPrmErrString[512];		/* Storage for error string */
static LbPfHkTy		detPrmFileHk;				/* Our parameter file */

/* The following table must correspond, in same order, to DetEn_RunTimeParamsTy
in DetParams.h */
static labelTy		detRunTimeParamLabels[] = {	/* Our label table */
					"detector_type",
					"photon_time_fwhm_ns",
					"energy_resolution_percentage",
					"reference_energy_keV",
					"do_forced_interaction",
					"history_file",
					"history_params_file",
					"poly_axial_cut_pos_list",
					"poly_transaxial_cut_pos_list",
					"poly_axial_cut_start",
					"poly_axial_cut_end",
					"poly_transaxial_cut_start",
					"poly_transaxial_cut_end",
					"poly_scintillator",
					"poly_radial_depth",
					"poly_num_detectors_transaxially",
					"poly_num_detectors_axially",
					"poly_layer_info_list",
					"poly_num_layers",
					"poly_cut_material",
					"poly_axial_side_cover_thickness",
					"poly_axial_side_cover_material",
					"poly_trnsax_side_cvr_thickness",
					"poly_trnsax_side_cvr_material",
					"poly_front_cover_thickness",
					"poly_front_cover_material",
					"poly_transaxial_position",
					"poly_axial_min",
					"poly_axial_max",
					"poly_block_info_list",
					"poly_transaxial_length",
					"poly_axial_length",
					"poly_num_blocks",
					"poly_inner_radius",
					"poly_num_rings",
					"poly_ring_info_list",
					"plnr_num_layers",
					"plnr_layer_material",
					"plnr_layer_depth",
					"plnr_layer_is_active",
					"plnr_layer_info_list",
					"plnr_axial_length",
					"plnr_transaxial_length",
					"plnr_inner_radius",
					"plnr_num_views",
					"plnr_min_angle",
					"plnr_max_angle",
					"cyln_num_layers",
					"cyln_layer_material",
					"cyln_layer_inner_radius",
					"cyln_layer_outer_radius",
					"cyln_layer_is_active",
					"cyln_layer_info_list",
					"cyln_gap_material",
					"cyln_transaxial_gap",
					"cyln_axial_gap",
					"cyln_num_transaxial",
					"cyln_num_axial",
					"cyln_min_z",
					"cyln_max_z",
					"cyln_ring_info_list",
					"cyln_num_rings",
					"blocktomo_position_algorithm",
					"blocktomo_num_rings",
					"blocktomo_ring_description_list",
					"blocktomo_ring_parameter_file",
					"blocktomo_ring_axial_shift",
					"blocktomo_ring_transaxial_rotation",
					"ring_x_inner_radius",
					"ring_x_outer_radius",
					"ring_y_inner_radius",
					"ring_y_outer_radius",
					"ring_z_minimum",
					"ring_z_maximum",
					"ring_num_blocks_in_ring",
					"ring_block_description_list",
					"ring_block_parameter_file",
					"ring_block_radial_position",
					"ring_block_angular_position",
					"ring_block_z_position",
					"ring_block_transaxial_orientation",
					"block_reference_x",
					"block_reference_y",
					"block_reference_z",
					"block_x_minimum",
					"block_x_maximum",
					"block_y_minimum",
					"block_y_maximum",
					"block_z_minimum",
					"block_z_maximum",
					"block_num_layers",
					"block_layer_info_list",
					"block_layer_inner_x",
					"block_layer_outer_x",
					"block_layer_num_y_changes",
					"block_layer_y_change",
					"block_layer_num_z_changes",
					"block_layer_z_change",
					"block_layer_material_elements",
					"block_material_info_list",
					"block_material_index",
					"block_material_is_active",
					"randoms_history_file",
					"coincidence_timing_window_in_ns",
					"triples_processing_method",
					""
					};

/* The following table must correspond, in same order, to DetEn_BlockDetectedPositionAlgoTy
in DetParams.h - the options for the detected photon position algorithm */
static labelTy		detBlockPosAlgoParamLabels[] = {	/* label table */
					"snap_centroid_to_crystal_center",
					"use_energy_weighted_centroid",
					""
					};

/* The following strings are used by the parameter code to determine which type
	of detector the user specified.
	Modification to these stings must be matched by mods to DetEn_DetectorTypeTy
in DetTypes.h

*/
static char *detEn_DetTypeStr[] = {
	"null",
	"simple_pet",
	"simple_spect",
	"unc_spect",
	"polygonal",
	"planar",
	"cylindrical",
	"dual_headed",
	"block"};
	
#define NUM_DET_TYPES 9

/* PROTOTYPES */
Boolean detParamsGetRingInfo(	LbFourByte	curRing,
								LbFourByte	curBlock,
								LbUsFourByte	*numActiveElementsInTomoPtr,
					DetBlockTomoRingTy		*curRingInfoPtr);

Boolean detParamsGetBlockInfo(	LbFourByte	curRing,
								LbFourByte	curBlock,
								LbUsFourByte	*numActiveElementsInTomoPtr,
					DetBlockTomoBlockTy		*curBlockInfoPtr,
					DetBlockTomoRingTy		*curRingInfoPtr);

/* FUNCTIONS */
/*********************************************************************************
*
*			Name:			DetGetRunTimeParams
*
*			Summary:	Read in the runtime parameters.
*			Arguments:
*
*			Function return: true unless an error occurs.
*
*********************************************************************************/
Boolean DetGetRunTimeParams()	
{
	double					paramBuffer[LBPF_PARAM_LEN];
	Boolean					okay = false;					/* Process flag */
	Boolean					isEOF;							/* End of file flag */
	Boolean					switchOkay;						/* Switch flag - is case statement successful? */
	Boolean					initOkay;						/* is block initialization successful? */
	LbPfEnPfTy				paramType;
	LbUsFourByte			paramSize;
	char					paramLabel[LBPF_LABEL_LEN];
	DetEn_RunTimeParamsTy	whichParam;
	LbFourByte				curRing = -1;
	LbFourByte				curBlock = -1;
	LbFourByte				curLayer = -1;
	LbFourByte				curElement = -1;
	LbFourByte				curCut = 0;				/* Yes, it is zero and not -1 */
	LbFourByte				curTransCut = 0;		/* Yes, it is zero and not -1 */
	LbFourByte				typeIndex;				/* For searching type list */
	DetBlockTomoRingTy		*curRingInfoPtr;		/* the current ring info (shortens lines in code!) */
	DetBlockTomoBlockTy		*curBlockInfoPtr;		/* the current block info (shortens lines in code!) */
	DetBlockTomoLayerTy		*curLayerInfoPtr;		/* the current layer info (shortens lines in code!) */
	double					ringOverlap;			/* amount of overlap if rings not properly layered */
	LbOneByte				magnitude = -7;			/* largest magnitude overlap acceptable in comparisons  */
	
	do { /* Process Loop */


		/* Clear out path name variables and flats in case of error */
		{
	 		memset(DetRunTimeParams[DetCurParams].DetHistoryFilePath, '\0', PATH_LENGTH);
	 		memset(DetRunTimeParams[DetCurParams].DetHistoryParamsFilePath, '\0', PATH_LENGTH);
			DetRunTimeParams[DetCurParams].DetectorType = DetEn_DetType_NULL;
		}
		
		/* Attempt to open parameter file */
		if (!LbPfOpen(PhgRunTimeParams.PhgDetectorParamsFilePath[DetCurParams], 0, &detPrmFileHk)) {
			LbInPrintf("An error occurred opening the detector parameters file '%s'\n"
				"Check your PHG parameters file, '%s' for a valid path\n", PhgRunTimeParams.PhgDetectorParamsFilePath[DetCurParams],
				PhgRunTimeParams.PhgParamFilePath);
			 
			break;
		}
		
		/* Clear our detector type for error handling */
		DetRunTimeParams[DetCurParams].DetectorType = DetEn_DetType_NULL;
		
		/* Initialize some parameter fields */
		DetRunTimeParams[DetCurParams].PlanarDetector.NumViews = 0;
		DetRunTimeParams[DetCurParams].PlanarDetector.MinAngle = 0;
		
		if (PHG_IsSPECT())
			DetRunTimeParams[DetCurParams].PlanarDetector.MaxAngle = PHGMATH_2PI;
		else
			DetRunTimeParams[DetCurParams].PlanarDetector.MaxAngle = PHGMATH_PI;

		DetRunTimeParams[DetCurParams].PlanarDetector.NumViews = 0;
		
		DetRunTimeParams[DetCurParams].BlockTomoDetector.BlockDetectedPositionAlgo = DetEn_BlockSnapCentroidToXtalCenter;
		
		/* Loop through the parameters */
		while(LbPfGetParam(&detPrmFileHk, (void *)paramBuffer,
				&paramType, &paramSize, paramLabel, &isEOF)) {
				
			/* Find the runtime parameter */
			whichParam = DetLookupRunTimeParamLabel(paramLabel);
							
			/* Verify that the type of detector is the first parameter */
			if ((DetRunTimeParams[DetCurParams].DetectorType == DetEn_DetType_NULL) &&
					(whichParam != DetEn_detector_type)){
					
				ErStGeneric("Detector type must be first parameter in detector param file!\n");
				break;
			}
							
			switchOkay = true;
			switch (whichParam) {


				case DetEn_detector_type:
				
					/* See if they gave us an enum or an int */
					if (paramType == LbPfEnEnum) {
						/* Search table for given label */
						for (typeIndex = 0; typeIndex < (NUM_DET_TYPES); typeIndex++){
						
							/* See if it matches */
							if (strcmp((char *)paramBuffer, detEn_DetTypeStr[typeIndex]) == 0) {
								DetRunTimeParams[DetCurParams].DetectorType = (DetEn_DetectorTypeTy) typeIndex;
								break;
							}
						}
						if (typeIndex == NUM_DET_TYPES) {
							LbInPrintf("Invalid detector type in detector parameters '%s', valid types are:\n",
								(char *)paramBuffer);
							for (typeIndex = 1; typeIndex < (NUM_DET_TYPES); typeIndex++){
								LbInPrintf("'%s'\n", detEn_DetTypeStr[typeIndex]);
							}
							ErStGeneric("An invalid detector type was supplied");
							switchOkay = false;
						}
					}
					else if (paramType == LbPfEnInteger) {
						if ((*(LbFourByte *)paramBuffer) < NUM_DET_TYPES) {
							DetRunTimeParams[DetCurParams].DetectorType =
								(DetEn_DetectorTypeTy) *((LbUsFourByte *) paramBuffer);
						}
						else {
							ErStGeneric("Error in detector parameter file, you specified an invalid detector type");
							switchOkay = false;
						}
					}
					else {
						ErStGeneric("Error in detector parameter file, you specified an invalid type for the detector parameter, it must be either INT or ENUM");
						switchOkay = false;
					}
					break;

				
				case DetEn_photon_tof_fwhm:
					DetRunTimeParams[DetCurParams].PhotonTimeFWHM =
						*((double *) paramBuffer);
					break;

				case DetEn_energy_resolution_percentage:
					DetRunTimeParams[DetCurParams].EnergyResolutionPercentage =
						*((double *) paramBuffer);
					break;

				case DetEn_reference_energy_keV:
					DetRunTimeParams[DetCurParams].ReferenceEnergy =
						*((double *) paramBuffer);
					break;

				case DetEn_do_forced_interaction:
					DetRunTimeParams[DetCurParams].DoForcedInteraction =
						*((Boolean *) paramBuffer);
					break;
					
				case 	DetEn_poly_inner_radius:
					DetRunTimeParams[DetCurParams].PolyDetector.InnerRadius =
						*((double *) paramBuffer);
					break;
					
				case 	DetEn_poly_num_rings:
					DetRunTimeParams[DetCurParams].PolyDetector.NumRings =
						*((LbUsFourByte *) paramBuffer);
						
					/* Now attempt to allocate the ring storage */
					if ((DetRunTimeParams[DetCurParams].PolyDetector.RingInfo =
							(DetPolRingTy *) LbMmAlloc(
							DetRunTimeParams[DetCurParams].PolyDetector.NumRings *
							sizeof(DetPolyBlockTy))) == 0) {
							
						switchOkay = false;
					}
					break;
					
				case 	DetEn_poly_ring_info_list:
				
					/*	We update the current ring each time we read this label. */
					curRing++;
					
					break;
					
				case 	DetEn_poly_num_blocks:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].NumBlocks = 
						*((LbUsFourByte *) paramBuffer);
					
						
					/* Now attempt to allocate the block storage */
					if ((DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo =
							(DetPolyBlockTy *) LbMmAlloc(
							DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].NumBlocks *
							sizeof(DetPolyBlockTy))) == 0) {
							
						switchOkay = false;
					}
					break;						
					
				case 	DetEn_poly_block_info_list:
				
					/* First verify that our block info array has been allocated */
					if (DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo == 0) {
						ErStGeneric("Error in parameter file, you probably forgot to specify"
							" the number of blocks 'poly_num_blocks', see the example files for help");
							
						switchOkay = false;
						
						break;
					}
				
					/*	We update the current block each time we read this label. */
					curBlock++;
						
					break;
					
				case 	DetEn_poly_num_layers:
				
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].NumLayers = 
						*((LbUsFourByte *) paramBuffer);
											
					/* Now attempt to allocate the layer storage */
					if ((DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].LayerInfo =
							(DetPolyLayerInfoTy *) LbMmAlloc(
							DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].NumLayers *
							sizeof(DetPolyLayerInfoTy))) == 0) {
							
						switchOkay = false;
					}
						
					break;

				case 	DetEn_poly_transaxial_position:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].TransaxialPos = 
						*((double *) paramBuffer);
					
					break;
					

				case 	DetEn_poly_axial_min:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].AxialMin = 
						*((double *) paramBuffer);
					
					break;

				case 	DetEn_poly_axial_max:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].AxialMax = 
						*((double *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_front_cover_material:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].FrontCvrMaterial = 
						*((LbUsFourByte *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_front_cover_thickness:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].FrontCvrThickness = 
						*((double *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_trnsax_side_cvr_material:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].TransSideCvrMaterial = 
						*((LbUsFourByte *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_trnsax_side_cvr_thickness:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].TransSideCvrThickness = 
						*((double *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_axial_side_cover_material:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].AxialSideCvrMaterial = 
						*((LbUsFourByte *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_axial_side_cover_thickness:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].AxialSideCvrThickness = 
						*((double *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_cut_material:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].CutMaterial = 
						*((LbUsFourByte *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_layer_info_list:
				
					/* First verify that our layer info array has been allocated */
					if (DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].LayerInfo == 0) {
						ErStGeneric("Error in parameter file, you probably forgot to specify"
							" the number of layers 'poly_num_layers', see the example files for help");
							
						switchOkay = false;
						
						break;
					}
				
					/*	We increment the current layer each time we read this label. */
					curLayer++;
					break;
					
				case 	DetEn_poly_num_detectors_axially:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].LayerInfo[curLayer].NumMatAxially = 
						*((LbUsFourByte *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_num_detectors_transaxially:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].LayerInfo[curLayer].NumMatTransaxially = 
						*((LbUsFourByte *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_radial_depth:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].LayerInfo[curLayer].RadialDepth = 
						*((double *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_axial_cut_pos_list:
				
					/*	We just reset the current cut to the first one. */
					curCut = 0;
					
					break;
					
				case 	DetEn_poly_axial_cut_start:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].LayerInfo[curLayer].AxialCutPos[curCut].CutStart = 
						*((double *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_axial_cut_end:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].LayerInfo[curLayer].AxialCutPos[curCut].CutEnd = 
						*((double *) paramBuffer);
					
					/* This ends the current cut so we advance to the next cut */
					curCut++;
					
					break;
					
				case 	DetEn_poly_transaxial_cut_pos_list:
				
					/*	We just reset the current cut to the first one. */
					curTransCut = 0;
					
					break;
					
				case 	DetEn_poly_transaxial_cut_start:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].LayerInfo[curLayer].TransaxialCutPos[curTransCut].CutStart = 
						*((double *) paramBuffer);
					
					break;
					
				case 	DetEn_poly_transaxial_cut_end:
					DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[curRing].BlockInfo[curBlock].LayerInfo[curLayer].TransaxialCutPos[curTransCut].CutEnd = 
						*((double *) paramBuffer);
					
					/* This ends the current cut so we advance to the next cut */
					curTransCut++;

					break;
					
				case 	DetEn_plnr_num_layers:
					DetRunTimeParams[DetCurParams].PlanarDetector.NumLayers =
						*((LbUsFourByte *) paramBuffer);
											
					/* Now attempt to allocate the layer storage */
					if ((DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo =
							(DetPlnrLayerInfoTy *) LbMmAlloc(
							DetRunTimeParams[DetCurParams].PlanarDetector.NumLayers *
							sizeof(DetPlnrLayerInfoTy))) == 0) {
							
						switchOkay = false;
					}
						
					break;
					
				case 	DetEn_plnr_inner_radius:
					DetRunTimeParams[DetCurParams].PlanarDetector.InnerRadius =
						*((double *) paramBuffer);
					break;
					
				case 	DetEn_plnr_transaxial_length:
					DetRunTimeParams[DetCurParams].PlanarDetector.TransaxialLength =
						*((double *) paramBuffer);
					break;
					
				case 	DetEn_plnr_axial_length:
					DetRunTimeParams[DetCurParams].PlanarDetector.AxialLength =
						*((double *) paramBuffer);
					break;
					
				case 	DetEn_plnr_layer_info_list:
					/* We update the current layer each time we hit this */
					curLayer++;
					break;
					
				case 	DetEn_plnr_layer_depth:
					DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerDepth =
						*((double *) paramBuffer);
					break;
					
				case 	DetEn_plnr_layer_material:
					DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].LayerMaterial =
						*((LbUsFourByte *) paramBuffer);
					break;
					
				case 	DetEn_plnr_layer_is_active:
					DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[curLayer].IsActive =
						*((Boolean *) paramBuffer);
					break;
					
				case 	DetEn_plnr_num_views:
					DetRunTimeParams[DetCurParams].PlanarDetector.NumViews =
						*((LbUsFourByte *) paramBuffer);
					break;
					
				case 	DetEn_plnr_min_angle:
					DetRunTimeParams[DetCurParams].PlanarDetector.MinAngle =
						PHGMATH_RadiansFromDegrees(*((double *) paramBuffer));
					break;
					
				case 	DetEn_plnr_max_angle:
					DetRunTimeParams[DetCurParams].PlanarDetector.MaxAngle =
						PHGMATH_RadiansFromDegrees(*((double *) paramBuffer));
					break;
					
				case 	DetEn_cyln_num_rings:
					DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings =
						*((LbUsFourByte *) paramBuffer);

					/* Now attempt to allocate the layer storage */
					if ((DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo =
							(DetCylnRingTy *) LbMmAlloc(
							DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings *
							sizeof(DetCylnRingTy))) == 0) {
							
						switchOkay = false;
					}
					break;
					
				case	DetEn_cyln_ring_info_list:

					/*	Increment the current ring each time this label is read. */
					curRing++;

					/* Reset the layer info each time this label is read */
					curLayer = -1;
					
					break;
					
				case 	DetEn_cyln_num_layers:
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].NumLayers =
						*((LbUsFourByte *) paramBuffer);

					/* Now attempt to allocate the layer storage */
					if ((DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].LayerInfo =
							(DetCylnLayerInfoTy *) LbMmAlloc(
							DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].NumLayers *
							sizeof(DetCylnLayerInfoTy))) == 0) {
							
						switchOkay = false;
					}
					break;
					
				case 	DetEn_cyln_min_z:
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].MinZ =
						*((double *) paramBuffer);

					break;
					
				case 	DetEn_cyln_max_z:
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].MaxZ =
						*((double *) paramBuffer);

					break;
					
				case 	DetEn_cyln_num_axial:
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].NumAxial =
						*((LbUsFourByte *) paramBuffer);

					break;
					
				case 	DetEn_cyln_num_transaxial:
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].NumTransaxial =
						*((LbUsFourByte *) paramBuffer);

					break;
					
				case 	DetEn_cyln_axial_gap:
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].AxialGap =
						*((double *) paramBuffer);

					break;
					
				case 	DetEn_cyln_transaxial_gap:
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].TransaxialGap =
						*((double *) paramBuffer);

					break;
					
				case 	DetEn_cyln_gap_material:
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].GapMaterial =
						*((LbUsFourByte *) paramBuffer);

					break;
					
				case 	DetEn_cyln_layer_info_list:
					/* We update the current layer each time we read this label */
					curLayer++;
					
					break;
					
				case 	DetEn_cyln_layer_material:
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].LayerInfo[curLayer].LayerMaterial =
						*((LbUsFourByte *) paramBuffer);

					break;
										
					
				case 	DetEn_cyln_layer_inner_radius:
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].LayerInfo[curLayer].InnerRadius =
						*((double *) paramBuffer);

					break;

				case 	DetEn_cyln_layer_outer_radius:
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].LayerInfo[curLayer].OuterRadius =
						*((double *) paramBuffer);

					break;


				case 	DetEn_cyln_layer_is_active:
					DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[curRing].LayerInfo[curLayer].IsActive =
						*((Boolean *) paramBuffer);

					break;

				case 	DetEn_blocktomo_position_algorithm:
					/* See if they gave us an enum or an int */
					if (paramType == LbPfEnEnum) {
						/* Search table for given label */
						for (typeIndex = 0; typeIndex < (LbFourByte)DetEn_BlockPosAlgoNULL; typeIndex++){
						
							/* See if it matches */
							if (strcmp((char *)paramBuffer, detBlockPosAlgoParamLabels[typeIndex]) == 0) {
								DetRunTimeParams[DetCurParams].BlockTomoDetector.BlockDetectedPositionAlgo = (DetEn_BlockDetectedPositionAlgoTy) typeIndex;
								break;
							}
						}
						if (typeIndex == (LbFourByte)DetEn_BlockPosAlgoNULL) {
							LbInPrintf("\nInvalid detected position algorithm, '%s%', in block detector parameters.  Valid algorithms:\n",
								(char *)paramBuffer);
							for (typeIndex = 0; typeIndex < (LbFourByte)DetEn_BlockPosAlgoNULL; typeIndex++){
								LbInPrintf("'%s'\n", detBlockPosAlgoParamLabels[typeIndex]);
							}
							ErStGeneric("An invalid detected position algorithm was supplied");
							switchOkay = false;
						}
					}
					else if (paramType == LbPfEnInteger) {
						if ((*(LbFourByte *)paramBuffer) < (LbFourByte)DetEn_BlockPosAlgoNULL) {
							DetRunTimeParams[DetCurParams].BlockTomoDetector.BlockDetectedPositionAlgo =
								(DetEn_BlockDetectedPositionAlgoTy) *((LbUsFourByte *) paramBuffer);
						}
						else {
							ErStGeneric("Error in block detector parameter file, you specified an invalid detected position algorithm");
							switchOkay = false;
						}
					}
					else {
						ErStGeneric("Error in block detector parameter file, you specified an invalid type for the detected position algorithm, it must be either INT or ENUM");
						switchOkay = false;
					}
					break;
					
				case 	DetEn_blocktomo_num_rings:
					DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings =
						*((LbUsFourByte *) paramBuffer);

					/* check the number of rings is valid */
					if (DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings < 1) {
						ErStGeneric("Invalid number of rings for block detector tomograph.");
						switchOkay = false;
						break;
					}

					/* Now attempt to allocate the ring storage */
					if ((DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo =
							(DetBlockTomoRingTy *) LbMmAlloc(
							DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings *
							sizeof(DetBlockTomoRingTy))) == 0) {
							
						switchOkay = false;
						break;
					}
					
					/* initialize the number of active elements in the tomograph */
					DetRunTimeParams[DetCurParams].BlockTomoDetector.NumActiveElementsInTomo = 0;
					
					/* initialize the rings - we will later check they are filled */
					for ( curRing = 0; curRing < (LbFourByte)(DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings); curRing++ ) {
						
						curRingInfoPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[curRing]);
						curRingInfoPtr->AxialShift = MAXFLOAT;
						curRingInfoPtr->TransaxialRotation = MAXFLOAT;
						curRingInfoPtr->XInnerRadius = 0.0;
						curRingInfoPtr->XOuterRadius = 0.0;
						curRingInfoPtr->YInnerRadius = 0.0;
						curRingInfoPtr->YOuterRadius = 0.0;
						curRingInfoPtr->MinZ = MAXFLOAT;
						curRingInfoPtr->MaxZ = -MAXFLOAT;
						curRingInfoPtr->NumBlocks = 0;
						curRingInfoPtr->NumActiveElementsInRing = 0;
						curRingInfoPtr->FirstActiveElementInRing = -1;
						
					}
					
					/* re-initialize curRing */
					curRing = -1;
					
					break;
					
				case	DetEn_blocktomo_ring_description_list:
					
					/*	Increment the current ring each time this label is read. */
					curRing++;
					
					/*  assign pointer to current ring */
					curRingInfoPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[curRing]);

					/* Reset the block info each time this label is read */
					curBlock = -1;
					
					break;
					
				case 	DetEn_blocktomo_ring_parameter_file:
					strcpy(curRingInfoPtr->RingParamFilePath,
						(char *) paramBuffer);
					
					/* process the ring parameter file */
					if ( !(detParamsGetRingInfo(curRing, curBlock,
							&(DetRunTimeParams[DetCurParams].BlockTomoDetector.NumActiveElementsInTomo),
							curRingInfoPtr)) )
						switchOkay = false;
					
					break;
					
				case 	DetEn_blocktomo_ring_axial_shift:
					curRingInfoPtr->AxialShift =
						*((double *) paramBuffer);

					break;
					
				case 	DetEn_blocktomo_ring_transaxial_rotation:
					curRingInfoPtr->TransaxialRotation =
						*((double *) paramBuffer);

					break;
					
				case DetEn_history_file:
						strcpy(DetRunTimeParams[DetCurParams].DetHistoryFilePath,
							(char *) paramBuffer);
						
						DetRunTimeParams[DetCurParams].DoHistory = 
							(DetRunTimeParams[DetCurParams].DetHistoryFilePath[0] != '\0');
					break;

				case DetEn_history_params_file:
						strcpy(DetRunTimeParams[DetCurParams].DetHistoryParamsFilePath,
							(char *) paramBuffer);
						
						/* !RH test for error trapping if history and custom parameters are aren't set
						DetRunTimeParams[DetCurParams].DoHistory = */
						DetRunTimeParams[DetCurParams].DoCustomHistory = 
							(DetRunTimeParams[DetCurParams].DetHistoryParamsFilePath[0] != '\0');
					break;


				case DetEn_randoms_history_file:
						strcpy(DetRunTimeParams[DetCurParams].DetRandomsHistoryFilePath,
							(char *) paramBuffer);
						
					DetRunTimeParams[DetCurParams].DoRandomsProcessing =
							(DetRunTimeParams[DetCurParams].DetRandomsHistoryFilePath[0] != '\0');
					break;


				case 	DetEn_coincidence_timing_window_in_ns:
					DetRunTimeParams[DetCurParams].CoincidenceTimingWindowNS =
						*((double *) paramBuffer);

					break;


				case DetEn_NULL:
					sprintf(detPrmErrString, "(DetGetRunTimeParams) Unknown (hence unused) parameter (%s).\n",
                        paramLabel);
					ErAlert(detPrmErrString, false);
					break;
				
				default:
					sprintf(detPrmErrString, "(DetGetRunTimeParams) Unknown (hence unused) parameter (%s).\n",
                        paramLabel);
					ErAlert(detPrmErrString, false);
					break;
			}
			
			if (!switchOkay)
				break;

		}
		
		/* Close the parameter file */
		LbPfClose(&detPrmFileHk);
		
		/* See if we quit due to error */
		if (!isEOF) {
			
			/* Print error string from above if one was given */
			/*if (detPrmErrString[0] != '\0') { */
				LbInPrintf( detPrmErrString );
			/* } */

			/* Let them know where/when we failed */
			sprintf(detPrmErrString, "Failed to process detector parameter file:\n"
				"\t%s (DetGetRunTimeParams)\n",PhgRunTimeParams.PhgDetectorParamsFilePath[DetCurParams]);
			ErAlert(detPrmErrString, false);
			
			break;
		}
		
		/* Verify that the number of parameters specified were supplied */
		if ((DetRunTimeParams[DetCurParams].DetectorType == DetEn_Planar) || (DetRunTimeParams[DetCurParams].DetectorType == DetEn_DualHeaded)){
			if ((LbFourByte)(DetRunTimeParams[DetCurParams].PlanarDetector.NumLayers) != (curLayer+1)) {
				sprintf(detPrmErrString, "(DetGetRunTimeParams) You specified %ld layers but supplied information for %ld.",
					(unsigned long)DetRunTimeParams[DetCurParams].PlanarDetector.NumLayers, 
					(long)curLayer);
					ErStGeneric(detPrmErrString);
				break;
			}
		}
			
		/* if this is a block detector, check that all necessary fields
		 were filled out - many fields were checked in detParamsGetRingInfo
		 and detParamsGetBlockInfo.  They are not checked again here.  Also
		 check that rings are contiguous and in order and 
		 call function to check for block collisions. */

		initOkay = true;	/* flag used for errors in nested do/for loops */

		if (DetRunTimeParams[DetCurParams].DetectorType == DetEn_Block) {
			
			if ( !(DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings > 0) ) {
				ErStGeneric("(DetGetRunTimeParams) No number of rings supplied for block detector.");
				initOkay = false;
				break;
			}
						
			for (curRing = 0; curRing < (LbFourByte)(DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings); curRing++) {
				
				curRingInfoPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[curRing]);
				
				if ( !(curRingInfoPtr->AxialShift < MAXFLOAT) ) {
					sprintf(detPrmErrString, "(DetGetRunTimeParams) No axial shift supplied for ring %ld.", (long)curRing);
					ErStGeneric(detPrmErrString);
					initOkay = false;
					break;
				}
				
				if ( !(curRingInfoPtr->TransaxialRotation < MAXFLOAT) ) {
					sprintf(detPrmErrString, "(DetGetRunTimeParams) No transaxial rotation supplied for ring %ld.", (long)curRing);
					ErStGeneric(detPrmErrString);
					initOkay = false;
					break;
				}
				
				if ( curRing > 0 ) {
					if ( PhgMathRealNumIsGreater(
							(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[curRing-1].MaxZ + 
							 DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[curRing-1].AxialShift), 
							(curRingInfoPtr->MinZ + curRingInfoPtr->AxialShift),
							magnitude, &ringOverlap)
						) {
						sprintf(detPrmErrString, "(DetGetRunTimeParams) Invalid minimum z supplied for ring %ld.\n"
								"Ring %ld (minimum + axial shift) must >= ring %ld (maximum + axial shift).\n",
								(long)curRing, (long)curRing, (long)(curRing-1) );
						ErStGeneric(detPrmErrString);
						initOkay = false;
						break;
					}
				}
			

				/* the following is also checked in detParamsGetRingInfo, but here serves to check that
				 detParamsGetRingInfo was run for each ring */
				if ( !(curRingInfoPtr->NumBlocks > 0) ) {
					sprintf(detPrmErrString, "(DetGetRunTimeParams) No ring parameter file supplied for ring %ld.", (long)curRing);
					ErStGeneric(detPrmErrString);
					initOkay = false;
					break;
				}
				
				/* check block min and max z are equal to the ring min and max z and perform
				 time-of-flight and energy resolution initializations */
				for (curBlock = 0; curBlock < (LbFourByte)(curRingInfoPtr->NumBlocks); curBlock++) {
					
					curBlockInfoPtr = &(curRingInfoPtr->BlockInfo[curBlock]);
					
					if ( (curBlockInfoPtr->ZMin + curBlockInfoPtr->ZPosition) != curRingInfoPtr->MinZ ) {
						sprintf(detPrmErrString, "(DetGetRunTimeParams) Block and ring minimum z must be equal, ring %ld block %ld.", (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
						initOkay = false;
						break;
					}
					
					if ( (curBlockInfoPtr->ZMax + curBlockInfoPtr->ZPosition) != curRingInfoPtr->MaxZ ) {
						sprintf(detPrmErrString, "(DetGetRunTimeParams) Block and ring minimum z must be equal, ring %ld block %ld.", (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
						initOkay = false;
						break;
					}
					
					/* compute the block position and orientation in tomograph coordinates */
					{
						
						curBlockInfoPtr->AngularPositionTomo = curBlockInfoPtr->AngularPosition + ((PHGMATH_PI/180.0) * curRingInfoPtr->TransaxialRotation);
						curBlockInfoPtr->XPositionTomo = curBlockInfoPtr->RadialPosition * cos(curBlockInfoPtr->AngularPositionTomo);
						curBlockInfoPtr->YPositionTomo = curBlockInfoPtr->RadialPosition * sin(curBlockInfoPtr->AngularPositionTomo);
						curBlockInfoPtr->ZPositionTomo = curBlockInfoPtr->ZPosition + curRingInfoPtr->AxialShift;
						curBlockInfoPtr->BlockFaceAngle = curBlockInfoPtr->AngularPositionTomo + ((PHGMATH_PI/180.0) * curBlockInfoPtr->TransaxialOrientation);
						curBlockInfoPtr->CosBlockFace = cos(curBlockInfoPtr->BlockFaceAngle);
						curBlockInfoPtr->SinBlockFace = sin(curBlockInfoPtr->BlockFaceAngle);						
						
					}
					for (curLayer = 0; curLayer < (LbFourByte)(curBlockInfoPtr->NumLayers); curLayer++ ) {
						
						/* fill in the time-of-flight and energy resolution fields of the detector elements.
						 Currently these use one global value, but having individual values will allow users
						 to model crystal-to-crystal variation if they wish and will leave us the option of
						 including it in future versions. */
						curLayerInfoPtr = &(curBlockInfoPtr->LayerInfo[curLayer]);
						
						for (curElement = 0; curElement < (LbFourByte)(curLayerInfoPtr->NumElements); curElement++ ) {
							if (curLayerInfoPtr->ElementInfo[curElement].IsActive == true) {
								curLayerInfoPtr->ElementInfo[curElement].PhotonTimeFWHM = 
										DetRunTimeParams[DetCurParams].PhotonTimeFWHM;
								curLayerInfoPtr->ElementInfo[curElement].EnergyResolutionPercentage = 
										DetRunTimeParams[DetCurParams].EnergyResolutionPercentage;
								curLayerInfoPtr->ElementInfo[curElement].ReferenceEnergy = 
										DetRunTimeParams[DetCurParams].ReferenceEnergy;
							}
						}
						
					}
					
				}
				
				if (!initOkay) {
					break;
				}
			
			}
			
			if (!initOkay) {
				break;
			}
		
			
			
			/* check that no two blocks intersect (at more than a common
			  edge or corner) */
			if (!DetBlocValidateBlocks()) {
				ErStGeneric("(DetGetRunTimeParams) Block validation failed.");
				initOkay = false;
				break;
			}
						
		}

		if (switchOkay && initOkay) {
			okay = true;
		} else {
			/* print existing error string, then report failure to break on error above */
			LbInPrintf( detPrmErrString );
			sprintf(detPrmErrString, "(DetGetRunTimeParams) Improper break from error.");
			ErStGeneric(detPrmErrString);
		}
		
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			DetLookupRunTimeParamLabel
*
*			Summary:	Find the given label in the label table and return
*						its enumeration value.
*			Arguments:
*				char		*label	- The label to find.
*			Function return: Enumerated type corresponding to label.
*
*********************************************************************************/
DetEn_RunTimeParamsTy  DetLookupRunTimeParamLabel(char *label)	
{
	DetEn_RunTimeParamsTy	whichParam;			/* Parameter we find */
	LbUsFourByte			labelIndex;			/* Label we are checking */

	/* Set param to null */
	whichParam = DetEn_NULL;
	
	/* Search table for given label */
	for (labelIndex = 0; labelIndex < (DetEn_NULL); labelIndex++){
	
		/* See if it matches */
		if (strcmp(label, detRunTimeParamLabels[labelIndex]) == 0) {
			whichParam = (DetEn_RunTimeParamsTy) labelIndex;
			break;
		}
	}
	
	return (whichParam);
}


/*********************************************************************************
*
*			Name:			detParamsGetRingInfo
*
*			Summary:	Read in the parameters for the current ring.
*
*			Arguments:
*				LbFourByte	curRing,
*				LbFourByte	curBlock,
*				LbUsFourByte		*numActiveElementsInTomoPtr,
*				DetBlockTomoRingTy	*curRingInfoPtr.
*
*			Function return: true unless an error occurs.
*
*********************************************************************************/
Boolean detParamsGetRingInfo(	LbFourByte		curRing,
								LbFourByte		curBlock,
								LbUsFourByte	*numActiveElementsInTomoPtr,
					DetBlockTomoRingTy			*curRingInfoPtr)
{
	double					paramBuffer[LBPF_PARAM_LEN];
	Boolean					okay = false;		/* Process flag */
	Boolean					isEOF;				/* End of file flag */
	Boolean					switchOkay;			/* Switch flag */
	Boolean					initOkay;			/* flag for checking initialization */
	LbFourByte				curLayer = -1;
	LbPfEnPfTy				paramType;
	LbUsFourByte			paramSize;
	char					paramLabel[LBPF_LABEL_LEN];
	DetEn_RunTimeParamsTy	whichParam;
	DetBlockTomoBlockTy		*curBlockInfoPtr;	/* the current block info (shortens lines in code!) */
	LbPfHkTy				ringParamFileHk;	/* current ring param file */
	
	do { /* Process Loop */

		/* Attempt to open parameter file */
		if (!LbPfOpen(curRingInfoPtr->RingParamFilePath, 0, &ringParamFileHk)) {
			LbInPrintf("(detParamsGetRingInfo) An error occurred opening the ring %d parameter file '%s'.\n"
				"Check your detector parameters file, '%s' for a valid path.",
				curRing,
				curRingInfoPtr->RingParamFilePath,
				PhgRunTimeParams.PhgDetectorParamsFilePath[DetCurParams]);
			 
			break;
		}

		/* Loop through the parameters */
		while ( LbPfGetParam(&ringParamFileHk, (void *)paramBuffer,
				&paramType, &paramSize, paramLabel, &isEOF) ) {
			
			/* Find the runtime parameter */
			whichParam = DetLookupRunTimeParamLabel(paramLabel);
			
			switchOkay = true;
			switch (whichParam) {
			
				case 	DetEn_ring_x_inner_radius:
					curRingInfoPtr->XInnerRadius = *((double *) paramBuffer);

					break;
						
				case 	DetEn_ring_x_outer_radius:
					curRingInfoPtr->XOuterRadius = *((double *) paramBuffer);

					break;
						
				case 	DetEn_ring_y_inner_radius:
					curRingInfoPtr->YInnerRadius = *((double *) paramBuffer);

					break;
						
				case 	DetEn_ring_y_outer_radius:
					curRingInfoPtr->YOuterRadius = *((double *) paramBuffer);

					break;
						
				case 	DetEn_ring_z_minimum:
					curRingInfoPtr->MinZ = *((double *) paramBuffer);

					break;
						
				case 	DetEn_ring_z_maximum:
					curRingInfoPtr->MaxZ = *((double *) paramBuffer);

					break;
						
				case 	DetEn_ring_num_blocks_in_ring:
					curRingInfoPtr->NumBlocks = *((LbUsFourByte *) paramBuffer);

					/* check the number of blocks is valid */
					if (curRingInfoPtr->NumBlocks < 1) {
						sprintf(detPrmErrString, "(detParamsGetRingInfo) Invalid number of blocks for ring %ld.", (long)curRing);
						ErStGeneric(detPrmErrString);
						switchOkay = false;
						break;
					}

					/* Now attempt to allocate the block storage */
					if ((curRingInfoPtr->BlockInfo =
							(DetBlockTomoBlockTy *) LbMmAlloc(
							curRingInfoPtr->NumBlocks *
							sizeof(DetBlockTomoBlockTy))) == 0) {
							
						switchOkay = false;
						break;
					}
					
					/* initialize the blocks - we will later check they are filled */
					for ( curBlock = 0; curBlock < (LbFourByte)(curRingInfoPtr->NumBlocks); curBlock++ ) {
						
						curBlockInfoPtr = &(curRingInfoPtr->BlockInfo[curBlock]);
						curBlockInfoPtr->RadialPosition = MAXFLOAT;
						curBlockInfoPtr->AngularPositionDegrees = MAXFLOAT;
						curBlockInfoPtr->AngularPosition = MAXFLOAT;
						curBlockInfoPtr->XPositionTomo = MAXFLOAT;
						curBlockInfoPtr->YPositionTomo = MAXFLOAT;
						curBlockInfoPtr->ZPosition = MAXFLOAT;
						curBlockInfoPtr->TransaxialOrientation = MAXFLOAT;
						curBlockInfoPtr->XRef = MAXFLOAT;
						curBlockInfoPtr->YRef = MAXFLOAT;
						curBlockInfoPtr->ZRef = MAXFLOAT;
						curBlockInfoPtr->XMin = MAXFLOAT;
						curBlockInfoPtr->XMax = -MAXFLOAT;
						curBlockInfoPtr->YMin = MAXFLOAT;
						curBlockInfoPtr->YMax = -MAXFLOAT;
						curBlockInfoPtr->ZMin = MAXFLOAT;
						curBlockInfoPtr->ZMax = -MAXFLOAT;
						curBlockInfoPtr->NumLayers = 0;
						curBlockInfoPtr->NumActiveElementsInBlock = 0;
						curBlockInfoPtr->FirstActiveElementInBlock = -1;
						curBlockInfoPtr->NumActiveLayers = 0;
						
					}
					
					/* re-initialize curBlock */
					curBlock = -1;
					
					break;
					
				case	DetEn_ring_block_description_list:
					
					/*	Increment the current block each time this label is read. */
					curBlock++;
					
					/* check to see that we do not exceed the defined number of blocks */
					if (curBlock >= (LbFourByte)(curRingInfoPtr->NumBlocks)) {
						
						sprintf(detPrmErrString, "(detParamsGetRingInfo) Too many blocks defined in ring %ld.", (long)curRing);
						ErStGeneric(detPrmErrString);
						switchOkay = false;
						break;
						
					}
					
					/* set a pointer to the current block */
					curBlockInfoPtr = &(curRingInfoPtr->BlockInfo[curBlock]);

					/* Reset the layer info each time this label is read */
					curLayer = -1;
					
					break;
					
				case 	DetEn_ring_block_parameter_file:
					strcpy(curBlockInfoPtr->BlockParamFilePath, (char *) paramBuffer);
					
					/* process the block parameter file */
					if ( !(detParamsGetBlockInfo(curRing, curBlock, numActiveElementsInTomoPtr, curBlockInfoPtr, curRingInfoPtr)) )
						switchOkay = false;
					
					break;
					
				case 	DetEn_ring_block_radial_position:
					curBlockInfoPtr->RadialPosition = *((double *) paramBuffer);

					break;
					
				case 	DetEn_ring_block_angular_position:
					curBlockInfoPtr->AngularPositionDegrees = *((double *) paramBuffer);
					/* Angular position input in degrees, change to radians */
					curBlockInfoPtr->AngularPosition = (PHGMATH_PI/180.0) *
							curBlockInfoPtr->AngularPositionDegrees;

					break;
					
/*				Changed input of detector position to cylindrical coordiantes, see above
				case 	DetEn_ring_block_x_position:
					curBlockInfoPtr->XPosition = *((double *) paramBuffer);

					break;
					
				case 	DetEn_ring_block_y_position:
					curBlockInfoPtr->YPosition = *((double *) paramBuffer);

					break;
*/
					
				case 	DetEn_ring_block_z_position:
					curBlockInfoPtr->ZPosition = *((double *) paramBuffer);

					break;
					
				case 	DetEn_ring_block_transaxial_orientation:
					curBlockInfoPtr->TransaxialOrientation = *((double *) paramBuffer);

					break;
					
				default:
					/* All other parameters at this point are incorrect */
					sprintf(detPrmErrString, "(detParamsGetRingInfo) Found a parameter (%s) not related to rings.", paramLabel);
					ErStGeneric(detPrmErrString);
					initOkay = false;
					
					break;
					
			}
			
			if (!switchOkay)
				break;
			
		}		

		/* Close the parameter file */
		LbPfClose(&ringParamFileHk);
		
		if (!switchOkay) {
			break;
		}
		
		/* See if we quit due to error */
		if (!isEOF) {
			
			/* Let them know where/when we failed */
			sprintf(detPrmErrString, "Failed to process ring parameter file:\n"
				"\t%s (detParamsGetRingInfo).",curRingInfoPtr->RingParamFilePath);
			ErAlert(detPrmErrString, false);
			
			break;
		}
		
		/* check that the ring was completely filled out. */
		{
			initOkay = true;
			
			if ( !(curRingInfoPtr->XInnerRadius > 0.0) ) {
				sprintf(detPrmErrString, "(detParamsGetRingInfo) No x inner radius supplied for ring %ld.", (long)curRing);
				ErStGeneric(detPrmErrString);
				initOkay = false;
				break;
			}
			
			if ( !(curRingInfoPtr->XOuterRadius > 0.0) ) {
				sprintf(detPrmErrString, "(detParamsGetRingInfo) No x outer radius supplied for ring %ld.", (long)curRing);
				ErStGeneric(detPrmErrString);
				initOkay = false;
				break;
			}
			
			if ( !(curRingInfoPtr->XOuterRadius > curRingInfoPtr->XInnerRadius) ) {
				sprintf(detPrmErrString, "(detParamsGetRingInfo) Invalid x outer radius supplied for ring %ld - must be > inner radius.", (long)curRing);
				ErStGeneric(detPrmErrString);
				initOkay = false;
				break;
			}
			
			if ( !(curRingInfoPtr->YInnerRadius > 0.0) ) {
				sprintf(detPrmErrString, "(detParamsGetRingInfo) No y inner radius supplied for ring %ld.", (long)curRing);
				ErStGeneric(detPrmErrString);
				initOkay = false;
				break;
			}
			
			if ( !(curRingInfoPtr->YOuterRadius > 0.0) ) {
				sprintf(detPrmErrString, "(detParamsGetRingInfo) No y outer radius supplied for ring %ld.", (long)curRing);
				ErStGeneric(detPrmErrString);
				initOkay = false;
				break;
			}
			
			if ( !(curRingInfoPtr->YOuterRadius > curRingInfoPtr->YInnerRadius) ) {
				sprintf(detPrmErrString, "(detParamsGetRingInfo) Invalid y outer radius supplied for ring %ld - must be > inner radius.", (long)curRing);
				ErStGeneric(detPrmErrString);
				initOkay = false;
				break;
			}
			
			if ( !(curRingInfoPtr->MinZ < MAXFLOAT) ) {
				sprintf(detPrmErrString, "(detParamsGetRingInfo) No minimum z supplied for ring %ld.", (long)curRing);
				ErStGeneric(detPrmErrString);
				initOkay = false;
				break;
			}
			
			if ( !(curRingInfoPtr->MaxZ > -MAXFLOAT) ) {
				sprintf(detPrmErrString, "(detParamsGetRingInfo) No maximum z supplied for ring %ld.", (long)curRing);
				ErStGeneric(detPrmErrString);
				initOkay = false;
				break;
			}
			
			if ( !(curRingInfoPtr->MaxZ > curRingInfoPtr->MinZ) ) {
				sprintf(detPrmErrString, "(detParamsGetRingInfo) Invalid maximum z supplied for ring %ld - must be > minimum z.", (long)curRing);
				ErStGeneric(detPrmErrString);
				initOkay = false;
				break;
			}
			
			if ( !(curRingInfoPtr->NumBlocks > 0) ) {
				sprintf(detPrmErrString, "(detParamsGetRingInfo) No number of blocks supplied for ring %ld.", (long)curRing);
				ErStGeneric(detPrmErrString);
				initOkay = false;
				break;
			}
			
			/* check that all blocks of the ring were filled in (much of this is 
			 done in detParamsGetBlockInfo - only those fields that aren't checked there
			 are checked heres).  Also compute x-y block position coords from cylindircal. */
			for (curBlock = 0; curBlock < (LbFourByte)(curRingInfoPtr->NumBlocks); curBlock++) {
			
				curBlockInfoPtr = &(curRingInfoPtr->BlockInfo[curBlock]);
				
				if ( !(curBlockInfoPtr->RadialPosition < MAXFLOAT) ) {
					sprintf(detPrmErrString, "(detParamsGetRingInfo) No radial position supplied for ring %ld, block %ld.", (long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					initOkay = false;
					break;
				}
				
				if ( !(curBlockInfoPtr->AngularPosition < MAXFLOAT) ) {
					sprintf(detPrmErrString, "(detParamsGetRingInfo) No angular position supplied for ring %ld, block %ld.", (long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					initOkay = false;
					break;
				}
				
				/* compute x-y from cylindrical coords
				{
					curBlockInfoPtr->XPosition = curBlockInfoPtr->RadialPosition * cos(curBlockInfoPtr->AngularPosition);
					curBlockInfoPtr->YPosition = curBlockInfoPtr->RadialPosition * sin(curBlockInfoPtr->AngularPosition);
				}
				
				if ( !(curBlockInfoPtr->XPosition < MAXFLOAT) ) {
					sprintf(detPrmErrString, "(detParamsGetRingInfo) No x position computed for ring %ld, block %ld.", (long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					initOkay = false;
					break;
				}
				
				if ( !(curBlockInfoPtr->YPosition < MAXFLOAT) ) {
					sprintf(detPrmErrString, "(detParamsGetRingInfo) No y position computed for ring %ld, block %ld.", (long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					initOkay = false;
					break;
				}  */
				
				if ( !(curBlockInfoPtr->ZPosition < MAXFLOAT) ) {
					sprintf(detPrmErrString, "(detParamsGetRingInfo) No z position supplied for ring %ld, block %ld.", (long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					initOkay = false;
					break;
				}
				
				if ( !(curBlockInfoPtr->TransaxialOrientation < MAXFLOAT) ) {
					sprintf(detPrmErrString, "(detParamsGetRingInfo) No transaxial orientation supplied for ring %ld, block %ld.", 
								(long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					initOkay = false;
					break;
				}
				
				/* the following is also checked in detParamsGetBlockInfo, but here serves to check that
				 detParamsGetBlockInfo was run for each block */
				if ( !(curBlockInfoPtr->NumLayers > 0) ) {
					sprintf(detPrmErrString, "(detParamsGetRingInfo) No block parameter file supplied for ring %ld, block %ld.", 
								(long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					initOkay = false;
					break;
				}

			}
			
			if (!initOkay) {
				break;
			}
			
			
		}

		if (switchOkay && initOkay) {
			okay = true;
		} else {
			/* print existing error string, then report failure to break on error above */
			LbInPrintf( detPrmErrString );
			sprintf(detPrmErrString, "(detParamsGetRingInfo) Improper break from error.");
			ErStGeneric(detPrmErrString);
		}
		
	} while (false);
	
	return (okay);

}

/*********************************************************************************
*
*			Name:			detParamsGetBlockInfo
*
*			Summary:	Read in the parameters for the current block.
*
*			Arguments:
*				LbFourByte	curRing,
*				LbFourByte	curBlock,
*				LbUsFourByte		*numActiveElementsInTomoPtr,
*				DetBlockTomoBlockTy	*curBlockInfoPtr
*				DetBlockTomoRingTy	*curRingInfoPtr.
*
*			Function return: true unless an error occurs.
*
*********************************************************************************/
Boolean detParamsGetBlockInfo(	LbFourByte	curRing,
								LbFourByte	curBlock,
								LbUsFourByte	*numActiveElementsInTomoPtr,
					DetBlockTomoBlockTy		*curBlockInfoPtr,
					DetBlockTomoRingTy		*curRingInfoPtr)
{
	double					paramBuffer[LBPF_PARAM_LEN];
	Boolean					okay = false;		/* Process flag */
	Boolean					isEOF;				/* End of file flag */
	Boolean					switchOkay;			/* Switch flag */
	Boolean					initOkay;			/* flag for checking initialization */
	LbPfEnPfTy				paramType;
	LbUsFourByte			paramSize;
	char					paramLabel[LBPF_LABEL_LEN];
	DetEn_RunTimeParamsTy	whichParam;
	LbFourByte				curLayer = -1;
	LbFourByte				curElement = -1;
	LbFourByte				curYChange = -1;
	LbFourByte				curZChange = -1;
	LbFourByte				curCrystalNum = -1;	/* number for the current crystal (active element) in block */
	LbFourByte				tempNumElements;
	DetBlockTomoLayerTy		*curLayerInfoPtr;	/* the current layer info (shortens lines in code!) */
	LbPfHkTy				blockParamFileHk;	/* current block param file */
	
	do { /* Process Loop */

		/* Attempt to open parameter file */
		if (!LbPfOpen(curBlockInfoPtr->BlockParamFilePath, 0, &blockParamFileHk)) {
			LbInPrintf("(detParamsGetBlockInfo) An error occurred opening the block %d in block %d parameter file '%s'\n"
				"Check your ring parameters file, '%s' for a valid path.",
				curBlock,
				curRing,
				curBlockInfoPtr->BlockParamFilePath,
				curRingInfoPtr->RingParamFilePath);
			 
			break;
		}

		/* Loop through the parameters */
		while ( LbPfGetParam(&blockParamFileHk, (void *)paramBuffer,
				&paramType, &paramSize, paramLabel, &isEOF) ) {
			
			/* Find the runtime parameter */
			whichParam = DetLookupRunTimeParamLabel(paramLabel);
			
			/* LbInPrintf("\nread block param\n"); */
			
			switchOkay = true;
			switch (whichParam) {
			
				case 	DetEn_block_reference_x:
					curBlockInfoPtr->XRef = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_reference_y:
					curBlockInfoPtr->YRef = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_reference_z:
					curBlockInfoPtr->ZRef = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_x_minimum:
					curBlockInfoPtr->XMin = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_x_maximum:
					curBlockInfoPtr->XMax = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_y_minimum:
					curBlockInfoPtr->YMin = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_y_maximum:
					curBlockInfoPtr->YMax = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_z_minimum:
					curBlockInfoPtr->ZMin = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_z_maximum:
					curBlockInfoPtr->ZMax = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_num_layers:
					curBlockInfoPtr->NumLayers = *((LbUsFourByte *) paramBuffer);

					/* check the number of layers is valid */
					if (curBlockInfoPtr->NumLayers < 1) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Invalid number of layers for ring %ld block %ld.", (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
						switchOkay = false;
						break;
					}

					/* Now attempt to allocate the layer storage */
					if ((curBlockInfoPtr->LayerInfo =
							(DetBlockTomoLayerTy *) LbMmAlloc(
							curBlockInfoPtr->NumLayers *
							sizeof(DetBlockTomoLayerTy))) == 0) {
							
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Unable to allocate layers for ring %ld block %ld.", (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
						
						switchOkay = false;
						break;
					}
					
					/* initialize the layers - we will later check they are filled */
					for ( curLayer = 0; curLayer < (LbFourByte)(curBlockInfoPtr->NumLayers); curLayer++ ) {
						
						curLayerInfoPtr = &(curBlockInfoPtr->LayerInfo[curLayer]);
						curLayerInfoPtr->InnerX = MAXFLOAT;
						curLayerInfoPtr->OuterX = -MAXFLOAT;
						curLayerInfoPtr->NumYChanges = -1;
						curLayerInfoPtr->NumZChanges = -1;
						curLayerInfoPtr->NumElements = 0;
						curLayerInfoPtr->NumActiveElementsInLayer = 0;
						
					}
					
					/* re-initialize curLayer */
					curLayer = -1;
					curElement = -1;
					curYChange = -1;
					curZChange = -1;
					
					break;

				case	DetEn_block_layer_info_list:
					
					/*	Increment the current layer each time this label is read. */
					curLayer++;
					
					/* check to see that we do not exceed the defined number of layers */
					if (curLayer >= (LbFourByte)(curBlockInfoPtr->NumLayers)) {
						
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Too many layers defined in block %ld ring %ld.", (long)curBlock, (long)curRing);
						ErStGeneric(detPrmErrString);
						
						switchOkay = false;
						break;
						
					}
					
					/* set a pointer to the current layer */
					curLayerInfoPtr = &(curBlockInfoPtr->LayerInfo[curLayer]);

					/* Reset the element info each time this label is read */
					curElement = -1;
					curYChange = -1;
					curZChange = -1;
					
					
					break;
					
				case 	DetEn_block_layer_inner_x:
					curLayerInfoPtr->InnerX = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_layer_outer_x:
					curLayerInfoPtr->OuterX = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_layer_num_y_changes:
					curLayerInfoPtr->NumYChanges = *((LbUsFourByte *) paramBuffer);
					
					/* check the number of y-changes is valid */
					if (curLayerInfoPtr->NumYChanges < 0) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Invalid number of y-changes for:"
									"\n    ring %ld block %ld layer %ld.", 
									(long)curRing, (long)curBlock, (long)curLayer);
						ErStGeneric(detPrmErrString);

						switchOkay = false;
						break;
					}

					/* Now attempt to allocate the y-changes storage (if there are any) */
					if (curLayerInfoPtr->NumYChanges > 0) {
					
						if ((curLayerInfoPtr->YChanges =
								(double *) LbMmAlloc(
								curLayerInfoPtr->NumYChanges *
								sizeof(double))) == 0) {
								
							sprintf(detPrmErrString, "(detParamsGetBlockInfo) Unable to allocate y-changes for:"
										"\n    ring %ld block %ld layer %ld.", 
										(long)curRing, (long)curBlock, (long)curLayer);
							ErStGeneric(detPrmErrString);
							
							switchOkay = false;
							break;
						}
						
						/* initialize the y-changes - we'll check later if they are filled */
						for ( curYChange = 0; curYChange < curLayerInfoPtr->NumYChanges; curYChange++ ) {
							curLayerInfoPtr->YChanges[curYChange] = -MAXFLOAT;
						}
						
					}
					
					/* if num-z-changes is already set, allocate and initialize the elements
					 - we will later check they are filled */
					if (curLayerInfoPtr->NumZChanges >= 0) {
						
						curLayerInfoPtr->NumElements = (curLayerInfoPtr->NumZChanges + 1) *
								(curLayerInfoPtr->NumYChanges + 1);
						
						if ((curLayerInfoPtr->ElementInfo =
								(DetBlockTomoElementTy *) LbMmAlloc(
								curLayerInfoPtr->NumElements *
								sizeof(DetBlockTomoElementTy))) == 0) {
								
							sprintf(detPrmErrString, "(detParamsGetBlockInfo) Unable to allocate elements (y) for:"
										"\n    ring %ld block %ld layer %ld.", 
										(long)curRing, (long)curBlock, (long)curLayer);
							ErStGeneric(detPrmErrString);
							
							switchOkay = false;
							break;
						}
						
						for ( curElement = 0; curElement < (LbFourByte)(curLayerInfoPtr->NumElements); curElement++ ) {
							curLayerInfoPtr->ElementInfo[curElement].MaterialIndex = SUBOBJGetNumTissues();
							curLayerInfoPtr->ElementInfo[curElement].IsActive = 2;
							curLayerInfoPtr->ElementInfo[curElement].crystalNumInBlock = -1;
							curLayerInfoPtr->ElementInfo[curElement].crystalNumInTomo = -1;
							curLayerInfoPtr->ElementInfo[curElement].PhotonTimeFWHM = 0;
							curLayerInfoPtr->ElementInfo[curElement].EnergyResolutionPercentage = 0;
							curLayerInfoPtr->ElementInfo[curElement].ReferenceEnergy = 0;
						}
						
					}
						
					/* re-initialize pointers */
					curElement = -1;
					curYChange = -1;
					
					break;

				case 	DetEn_block_layer_num_z_changes:
					curLayerInfoPtr->NumZChanges = *((LbUsFourByte *) paramBuffer);
					
					/* check the number of z-changes is valid */
					if (curLayerInfoPtr->NumZChanges < 0) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Invalid number of z-changes for ring %ld block %ld layer %ld.", 
									(long)curRing, (long)curBlock, (long)curLayer);
						ErStGeneric(detPrmErrString);
						
						switchOkay = false;
						break;
					}

					/* Now attempt to allocate the z-changes storage (if there are any) */
					if (curLayerInfoPtr->NumZChanges > 0) {
					
						if ((curLayerInfoPtr->ZChanges =
								(double *) LbMmAlloc(
								curLayerInfoPtr->NumZChanges *
								sizeof(double))) == 0) {
								
							sprintf(detPrmErrString, "(detParamsGetBlockInfo) Unable to allocate z-changes for:"
										"\n    ring %ld block %ld layer %ld.", 
										(long)curRing, (long)curBlock, (long)curLayer);
							ErStGeneric(detPrmErrString);
							
							switchOkay = false;
							break;
						}
						
						/* initialize the z-changes - we'll check later if they are filled */
						for ( curZChange = 0; curZChange < curLayerInfoPtr->NumZChanges; curZChange++ ) {
							curLayerInfoPtr->ZChanges[curZChange] = -MAXFLOAT;
						}
						
					}
					
					/* if num-y-changes is already set, allocate and initialize the elements
					 - we will later check they are filled */
					if (curLayerInfoPtr->NumYChanges >= 0) {
						
						curLayerInfoPtr->NumElements = (curLayerInfoPtr->NumZChanges + 1) *
								(curLayerInfoPtr->NumYChanges + 1);
						
						if ((curLayerInfoPtr->ElementInfo =
								(DetBlockTomoElementTy *) LbMmAlloc(
								curLayerInfoPtr->NumElements *
								sizeof(DetBlockTomoElementTy))) == 0) {
								
							sprintf(detPrmErrString, "(detParamsGetBlockInfo) Unable to allocate elements (z) for:"
										"\n    ring %ld block %ld layer %ld.", 
										(long)curRing, (long)curBlock, (long)curLayer);
							ErStGeneric(detPrmErrString);
														
							switchOkay = false;
							break;
						}
						
						for ( curElement = 0; curElement < (LbFourByte)(curLayerInfoPtr->NumElements); curElement++ ) {
							curLayerInfoPtr->ElementInfo[curElement].MaterialIndex = -1;
							curLayerInfoPtr->ElementInfo[curElement].IsActive = 2;
							curLayerInfoPtr->ElementInfo[curElement].crystalNumInBlock = -1;
							curLayerInfoPtr->ElementInfo[curElement].crystalNumInTomo = -1;
							curLayerInfoPtr->ElementInfo[curElement].PhotonTimeFWHM = 0;
							curLayerInfoPtr->ElementInfo[curElement].EnergyResolutionPercentage = 0;
							curLayerInfoPtr->ElementInfo[curElement].ReferenceEnergy = 0;
						}
						
					}
						
					/* re-initialize pointers */
					curElement = -1;
					curZChange = -1;
					
					break;

				case 	DetEn_block_layer_y_change:
					
					/* increment number of y-changes and check total */
					curYChange++;
					if (curYChange >= curLayerInfoPtr->NumYChanges) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Too many y-changes declared for:"
								"\n    layer %ld ring %ld block %ld.", 
								(long)curLayer, (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);

						switchOkay = false;
						break;
					}
					
					/* read y-change */
					curLayerInfoPtr->YChanges[curYChange] = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_layer_z_change:
					
					/* increment number of z-changes and check total */
					curZChange++;
					if (curZChange >= curLayerInfoPtr->NumZChanges) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Too many z-changes declared for:"
								"\n    layer %ld ring %ld block %ld.", 
								(long)curLayer, (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
						
						switchOkay = false;
						break;
					}
					
					/* read z-change */
					curLayerInfoPtr->ZChanges[curZChange] = *((double *) paramBuffer);

					break;
					
				case 	DetEn_block_layer_material_elements:
					
					tempNumElements = *((LbUsFourByte *) paramBuffer);
					
					if (tempNumElements != (LbFourByte)(curLayerInfoPtr->NumElements)) {
						
						if (curLayerInfoPtr->NumElements == 0) {
							sprintf(detPrmErrString, "(detParamsGetBlockInfo) Declare block_layer_num_y_changes and"
									"\n    num_z_changes before block_layer_material_elements for:\n"
									"    layer %ld ring %ld block %ld.", (long)curLayer, (long)curRing, (long)curBlock);
							ErStGeneric(detPrmErrString);
							
							switchOkay = false;
							break;
						} else {
							sprintf(detPrmErrString, "(detParamsGetBlockInfo) Declared NumElements and computed value differ for:"
									"\n    layer %ld ring %ld block %ld.",
									(long)curLayer, (long)curRing, (long)curBlock);
							ErStGeneric(detPrmErrString);
							
							switchOkay = false;
							break;
						}
					}
					
					curElement = -1;
					
					break;
					
				case 	DetEn_block_material_info_list:

					/*	Increment the current element each time this label is read. */
					curElement++;
					
					/* check to see that we do not exceed the defined number of blocks */
					if (curElement >= (LbFourByte)(curLayerInfoPtr->NumElements)) {
						
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Too many elements defined in:"
									"\n    layer %ld ring %ld block %ld.",
									(long)curLayer, (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
						
						switchOkay = false;
						break;
						
					}
					
					break;
					
				case 	DetEn_block_material_index:
					
					/* check to see if ElementInfo has been allocated */
					if (curElement >= 0) {
						
						curLayerInfoPtr->ElementInfo[curElement].MaterialIndex = *((LbFourByte *) paramBuffer);
						
					} else {
						
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) ElementInfo/material-index not allocated for:"
									"\n    layer %ld ring %ld block %ld.",
									(long)curLayer, (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
						
						switchOkay = false;
						break;
						
					}
					
					break;
					
				case 	DetEn_block_material_is_active:
					
					/* check to see if ElementInfo has been allocated */
					if (curElement >= 0) {
						
						curLayerInfoPtr->ElementInfo[curElement].IsActive = *((Boolean *) paramBuffer);
						
						/* if active, assign crystal numbers and increment active crystal numbers */
						if (curLayerInfoPtr->ElementInfo[curElement].IsActive) {
							
							/* increment local crystal counter assign crystal number in block */
							curCrystalNum++;
							curLayerInfoPtr->ElementInfo[curElement].crystalNumInBlock = curCrystalNum;
							
							/* assign crystal number in tomograph */
							curLayerInfoPtr->ElementInfo[curElement].crystalNumInTomo = *numActiveElementsInTomoPtr;
							
							/* if first active element in this block/ring, note number in block/ring info */
							if (curBlockInfoPtr->FirstActiveElementInBlock == -1) {
								curBlockInfoPtr->FirstActiveElementInBlock = *numActiveElementsInTomoPtr;
								if (curRingInfoPtr->FirstActiveElementInRing == -1) {
									curRingInfoPtr->FirstActiveElementInRing = *numActiveElementsInTomoPtr;
								}
							}
							
							/* increment the number of active elements in this layer/block/ring/tomo */
							*numActiveElementsInTomoPtr = *numActiveElementsInTomoPtr + 1;
							curBlockInfoPtr->NumActiveElementsInBlock++;
							curRingInfoPtr->NumActiveElementsInRing++;
							curLayerInfoPtr->NumActiveElementsInLayer++;
							
						}
						
					} else {
						
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) ElementInfo/material-is-active not allocated for:"
									"\n    layer %ld ring %ld block %ld.",
									(long)curLayer, (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
						
						switchOkay = false;
						break;
						
					}
					
					break;
					
				default:
					/* All other parameters at this point are incorrect */
					sprintf(detPrmErrString, "(detParamsGetBlockInfo) Found a parameter (%s) not related to blocks.", paramLabel);
					ErStGeneric(detPrmErrString);
					
					switchOkay = false;
					break;
					
			}
			
			if (!switchOkay)
				break;
			
		}		

		/* Close the parameter file */
		LbPfClose(&blockParamFileHk);
		
		if (!switchOkay) {
			break;
		}
		
		/* See if we quit switch due to error */
		if (!isEOF) {
			
			/* Print error string from above if one was given */
			/*if (detPrmErrString[0] != '\0') { */
				LbInPrintf( detPrmErrString );
			/* } */

			/* Let them know where/when we failed */
			sprintf(detPrmErrString, "Failed to process block parameter file:\n"
				"\t%s (detParamsGetBlockInfo).",curBlockInfoPtr->BlockParamFilePath);
			ErAlert(detPrmErrString, false);
			
			break;
		}
		
		/* now check that all required fields in the block were filled in */
		{
			
			initOkay = true;
			
			if ( !(curBlockInfoPtr->XRef < MAXFLOAT) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) No x reference coordinate supplied for ring %ld, block %ld.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}
			
			if ( !(curBlockInfoPtr->YRef < MAXFLOAT) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) No y reference coordinate supplied for ring %ld, block %ld.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}
			
			if ( !(curBlockInfoPtr->ZRef < MAXFLOAT) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) No z reference coordinate supplied for ring %ld, block %ld.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}
			
			if ( !(curBlockInfoPtr->XMin < MAXFLOAT) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) No minimum x supplied for ring %ld, block %ld.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}
			
			if ( !(curBlockInfoPtr->XMax > -MAXFLOAT) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) No maximum x supplied for ring %ld, block %ld.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}
			
			if ( !(curBlockInfoPtr->XMax > curBlockInfoPtr->XMin) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) Invalid maximum x supplied for:"
							"\n    ring %ld, block %ld - must be > min.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}
			
			if ( !(curBlockInfoPtr->YMin < MAXFLOAT) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) No minimum y supplied for ring %ld, block %ld.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}
			
			if ( !(curBlockInfoPtr->YMax > -MAXFLOAT) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) No maximum y supplied for ring %ld, block %ld.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}
			
			if ( !(curBlockInfoPtr->YMax > curBlockInfoPtr->YMin) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) Invalid maximum y supplied for:"
							"\n    ring %ld, block %ld - must be > min.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}
			
			if ( !(curBlockInfoPtr->ZMin < MAXFLOAT) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) No minimum z supplied for ring %ld, block %ld.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}
			
			if ( !(curBlockInfoPtr->ZMax > -MAXFLOAT) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) No maximum z supplied for ring %ld, block %ld.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}
			
			if ( !(curBlockInfoPtr->ZMax > curBlockInfoPtr->ZMin) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) Invalid maximum z supplied for:"
							"\n    ring %ld, block %ld - must be > min.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}
			
			if ( !(curBlockInfoPtr->NumLayers > 0) ) {
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) No number of layers supplied for ring %ld, block %ld.", 
							(long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				initOkay = false;
				break;
			}

			/* check that the layers of the block were filled in */
			for (curLayer = 0; curLayer < (LbFourByte)(curBlockInfoPtr->NumLayers); curLayer++) {
				
				curLayerInfoPtr = &(curBlockInfoPtr->LayerInfo[curLayer]);
				
				if ( !(curLayerInfoPtr->InnerX > -MAXFLOAT) ) {
					sprintf(detPrmErrString, "(detParamsGetBlockInfo) No inner x limit supplied for:"
								"\n   layer %ld, ring %ld, block %ld.",
								(long)curLayer, (long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					
					initOkay = false;
					break;
				}
				
				if ( !(curLayerInfoPtr->OuterX > -MAXFLOAT) ) {
					sprintf(detPrmErrString, "(detParamsGetBlockInfo) No outer x limit supplied for:"
								"\n    layer %ld, ring %ld, block %ld.",
								(long)curLayer, (long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					
					initOkay = false;
					break;
				}
				
				if ( !(curLayerInfoPtr->OuterX > curLayerInfoPtr->InnerX) ) {
					sprintf(detPrmErrString, "(detParamsGetBlockInfo) Invalid outer x limit supplied for:"
								"\n    layer %ld, ring %ld, block %ld - must be > inner.",
								(long)curLayer, (long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					
					initOkay = false;
					break;
				}
				
				if (curLayer == 0) {
					if (curLayerInfoPtr->InnerX != curBlockInfoPtr->XMin) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Invalid inner x limit supplied for layer 0:"
									"\n    ring %ld, block %ld - must be == to block x min.",
									(long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
						
						initOkay = false;
						break;
					}
				} else {
					if (curLayerInfoPtr->InnerX != curBlockInfoPtr->LayerInfo[curLayer-1].OuterX) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Invalid inner x limit supplied for:"
									"\n    layer %ld, ring %ld, block %ld - must be == outer x of previous layer.",
									(long)curLayer, (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
						
						initOkay = false;
						break;
					}
				}
				
				if (curLayer == ((LbFourByte)(curBlockInfoPtr->NumLayers) - 1)) {
					if (curLayerInfoPtr->OuterX != curBlockInfoPtr->XMax) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Invalid outer x limit supplied for last layer (%ld):"
									"\n    ring %ld, block %ld - must be == to block x min.",
									(long)curLayer, (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
						
						initOkay = false;
						break;
					}
				}
				
				if ( !(curLayerInfoPtr->NumYChanges >= 0) ) {
					sprintf(detPrmErrString, "(detParamsGetBlockInfo) Number of Y-material changes not supplied for:"
								"\n    layer %ld, ring %ld, block %ld.",
								(long)curLayer, (long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					
					initOkay = false;
					break;
				}
				
				if ( !(curLayerInfoPtr->NumZChanges >= 0) ) {
					sprintf(detPrmErrString, "(detParamsGetBlockInfo) Number of Z-material changes not supplied for:"
								"\n    layer %ld, ring %ld, block %ld.",
								(long)curLayer, (long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					
					initOkay = false;
					break;
				}

				if (curLayerInfoPtr->NumYChanges >= 1) {
					
					if (curLayerInfoPtr->YChanges[0] <= curBlockInfoPtr->YMin) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) First y-material change outside block for:"
									"\n    layer %ld, ring %ld, block %ld.",
									(long)curLayer, (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
					
						initOkay = false;
						break;
					}
					
					if (curLayerInfoPtr->YChanges[curLayerInfoPtr->NumYChanges - 1] >= curBlockInfoPtr->YMax) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Last y-material change outside block for"
									"\n    layer %ld, ring %ld, block %ld.",
									(long)curLayer, (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
					
						initOkay = false;
						break;
					}
					
				}
				
				if (curLayerInfoPtr->NumYChanges >= 2) {
					
					for (curYChange = 1; curYChange < curLayerInfoPtr->NumYChanges; curYChange++) {
						if (curLayerInfoPtr->YChanges[curYChange] <= curLayerInfoPtr->YChanges[curYChange - 1]) {
							sprintf(detPrmErrString, "(detParamsGetBlockInfo) Y-material changes must be strictly increasing:"
										"\n    layer %ld, ring %ld, block %ld.",
										(long)curLayer, (long)curRing, (long)curBlock);
							ErStGeneric(detPrmErrString);
						
							initOkay = false;
							break;
						}
					}
					
					if (!initOkay) {
						break;
					}
					
				}
				
				if (curLayerInfoPtr->NumZChanges >= 1) {
					
					if (curLayerInfoPtr->ZChanges[0] <= curBlockInfoPtr->ZMin) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) First z-material change outside block for:"
									"\n    layer %ld, ring %ld, block %ld.",
									(long)curLayer, (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
					
						initOkay = false;
						break;
					}
					
					if (curLayerInfoPtr->ZChanges[curLayerInfoPtr->NumZChanges - 1] >= curBlockInfoPtr->ZMax) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) Last z-material change outside block for"
									"\n    layer %ld, ring %ld, block %ld.",
									(long)curLayer, (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
					
						initOkay = false;
						break;
					}
					
				}
				
				if (curLayerInfoPtr->NumZChanges >= 2) {
					
					for (curZChange = 1; curZChange < curLayerInfoPtr->NumZChanges; curZChange++) {
						if (curLayerInfoPtr->ZChanges[curZChange] <= curLayerInfoPtr->ZChanges[curZChange - 1]) {
							sprintf(detPrmErrString, "(detParamsGetBlockInfo) Z-material changes must be strictly increasing:"
										"\n    layer %ld, ring %ld, block %ld.",
										(long)curLayer, (long)curRing, (long)curBlock);
							ErStGeneric(detPrmErrString);
						
							initOkay = false;
							break;
						}
					}
					
					if (!initOkay) {
						break;
					}
					
				}
				
				if (curLayerInfoPtr->NumElements != (LbUsFourByte)((curLayerInfoPtr->NumYChanges + 1) * (curLayerInfoPtr->NumZChanges + 1))) {
					sprintf(detPrmErrString, "(detParamsGetBlockInfo) Wrong number of elements for:"
								"\n    layer %ld, ring %ld, block %ld.",
								(long)curLayer, (long)curRing, (long)curBlock);
					ErStGeneric(detPrmErrString);
					
					initOkay = false;
					break;
				}
				
				for (curElement = 0; curElement < (LbFourByte)(curLayerInfoPtr->NumElements); curElement++) {
					if (!DetRunTimeParams[DetCurParams].DoRandomsProcessing) {
						/* currently the following check cannot be run when doing randoms as the SubObj module is
						not initialized.  This will change when addrandoms (and timesort) parameter files are
						incorporated into the run file */
						if ( !( (curLayerInfoPtr->ElementInfo[curElement].MaterialIndex < (LbFourByte)SUBOBJGetNumTissues()) &&
								(curLayerInfoPtr->ElementInfo[curElement].MaterialIndex >= 0) ) ) {
							sprintf(detPrmErrString, "(detParamsGetBlockInfo) Illegal MaterialIndex for:"
										"\n    element %ld layer %ld, ring %ld, block %ld.",
										(long)curElement, (long)curLayer, (long)curRing, (long)curBlock);
							ErStGeneric(detPrmErrString);
							
							initOkay = false;
							break;
						}
					}
					
					if ( !( (curLayerInfoPtr->ElementInfo[curElement].IsActive == true) ||
							(curLayerInfoPtr->ElementInfo[curElement].IsActive == false) ) ) {
						sprintf(detPrmErrString, "(detParamsGetBlockInfo) IsActive not set for:"
									"\n    element %ld layer %ld, ring %ld, block %ld.",
									(long)curElement, (long)curLayer, (long)curRing, (long)curBlock);
						ErStGeneric(detPrmErrString);
						
						initOkay = false;
						break;
					}
					
				}
				
				if (!initOkay) {
					break;
				}
			
			}
			
			if (!initOkay) {
				break;
			}

		}
		
		/* allocate and assign the block's ActiveLayers array and NumActiveLayers */
		{
			/* allocate ActiveLayers */
			if ((curBlockInfoPtr->ActiveLayers =
					(Boolean *) LbMmAlloc(
					curBlockInfoPtr->NumLayers *
					sizeof(Boolean))) == 0) {
					
				sprintf(detPrmErrString, "(detParamsGetBlockInfo) Unable to allocate active layer list for ring %ld block %ld.", (long)curRing, (long)curBlock);
				ErStGeneric(detPrmErrString);
				
				switchOkay = false;
				break;
			}

			/* assign true/false for ActiveLayers and count layers that have active elements */
			curBlockInfoPtr->NumActiveLayers = 0;
			for ( curLayer = 0; curLayer < (LbFourByte)(curBlockInfoPtr->NumLayers); curLayer++ ) {
				if ( curBlockInfoPtr->LayerInfo[curLayer].NumActiveElementsInLayer > 0 ) {
					curBlockInfoPtr->NumActiveLayers = curBlockInfoPtr->NumActiveLayers + 1;
					curBlockInfoPtr->ActiveLayers[curLayer] = true;
				} else {
					curBlockInfoPtr->ActiveLayers[curLayer] = false;
				}
			}
		}
		
		
		if (switchOkay && initOkay) {
			okay = true;
		} else {
			/* print existing error string, then report failure to break on error above */
			LbInPrintf( detPrmErrString );
			sprintf(detPrmErrString, "(detParamsGetBlockInfo) Improper break from error.");
			ErStGeneric(detPrmErrString);
		}
		
	} while (false);
	
	return (okay);

}


























#undef DET_PARAMS
