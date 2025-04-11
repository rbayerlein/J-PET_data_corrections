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
*     Module Name: 			ColParams.c
*     Revision Number:		1.4
*     Date last revised:	23 July 2013
*     Programmer:			Steven Vannoy
*     Date Originated:		11 July 1994
*
*     Module Overview:	The parameter module for the collimator.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*			ColLookupRunTimeParamLabel
*			ColGetRunTimeParams
*
*     Global variables defined:
*
*
*	  Global macros defined:
*
*********************************************************************************/
#define COL_PARAMS


#include <stdio.h>
#include <string.h>
#include <time.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
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
#include "phg.h"
#include "PhgBin.h"


/* LOCAL CONSTANTS */
/* LOCAL TYPES */
typedef char labelTy[LBPF_LABEL_LEN];

/* LOCAL GLOBALS */
static	char		colParamErrString[512];		/* Storage for error string s*/
static	LbPfHkTy	colParamFileHk;				/* Our parameter file */
static	labelTy		colRunTimeParamLabels[] = {
						"collimator_type",
						"collimator_depth",
						"layers_list",
						"segment_list",
						"segment",
						"seg_type",
						"material",
						"inner_radius",
						"outer_radius",
						"inner_min_z",
						"outer_min_z",
						"inner_max_z",
						"outer_max_z",
						"history_file",
						"history_parameters_file",
						"hole_geometry",
						"focal_length",
						"radius_of_rotation",
						"thickness",
						"hole_radius",
						"septal_thickness",
						"min_z",
						"max_z",
						"start_angle",
						"stop_angle",
						"sum_all_views",
						"num_views",
						"depth",
						"num_slats",
						"start",
						"end",
						"num_layers",
						""
						};


/* The following strings are used by the parameter code to determine which type
	of collimator the user specified.
*/
static char *colEn_ColTypeStr[] = {
	"null",
	"simple_pet",
	"monte_carlo_pet",
	"simple_spect",
	"unc_spect",
	"slat"};

#define NUM_COL_TYPES 6

/* The following strings are used by the parameter code to determine which type
	of collimator the user specified.
*/
static char *colEn_ColHoleTypeStr[] = {
	"null",
	"parallel",
	"fan_beam",
	"cone_beam"};
#define NUM_HOLE_TYPES 4

/* PROTOTYPES */

/* FUNCTIONS */
/*********************************************************************************
*
*			Name:			ColGetRunTimeParams
*
*			Summary:	Read in the runtime parameters.
*			Arguments:
*
*			Function return: true unless an error occurs.
*
*********************************************************************************/
Boolean ColGetRunTimeParams()	
{
	double					paramBuffer[LBPF_PARAM_LEN];
	Boolean					specifiedNumSlats = false;		/* Error checking flag */
	Boolean					okay = false;					/* Process flag */
	Boolean					isEOF;							/* End of file flag */
	Boolean					switchOkay;						/* Switch flag */
	LbPfEnPfTy				paramType;
	LbUsFourByte			paramSize;
	char					paramLabel[LBPF_LABEL_LEN];
	ColEn_RunTimeParamsTy	whichParam;
	LbFourByte				typeIndex;
	LbFourByte				curLayer = -1;
	LbFourByte				curSlat = -1;
	
	do { /* Process Loop */


		/* Clear out path name variables and flats in case of error */
		{
	 		memset(ColRunTimeParams[ColCurParams].ColHistoryFilePath, '\0', PATH_LENGTH);
	 		memset(ColRunTimeParams[ColCurParams].ColHistoryParamsFilePath, '\0', PATH_LENGTH);
			ColRunTimeParams[ColCurParams].ColType = ColEn_ColType_NULL;
		}
		
		/* Attempt to open parameter file */
		if (!LbPfOpen(PhgRunTimeParams.PhgCollimatorParamsFilePath[ColCurParams], 0, &colParamFileHk)) {
			LbInPrintf("An error occurred opening the collimator parameters file '%s'\n"
				"Check your PHG parameters file, '%s' for a valid path\n", PhgRunTimeParams.PhgCollimatorParamsFilePath[ColCurParams],
				PhgRunTimeParams.PhgParamFilePath);
			 
			break;
		}
		
		/* Clear our collimator type for error handling */
		ColRunTimeParams[ColCurParams].ColType = ColEn_ColType_NULL;
		
		/* Loop through the parameters */
		while(LbPfGetParam(&colParamFileHk, (void *)paramBuffer,
				&paramType, &paramSize, paramLabel, &isEOF)) {
				
				/* Find the runtime parameter */
				whichParam = ColLookupRunTimeParamLabel(paramLabel);
				
				/* Verify that the type of collimator is the first parameter */
				if ((ColRunTimeParams[ColCurParams].ColType == ColEn_NULL) &&
						(whichParam != ColEn_collimator_type)){
						
					ErStGeneric("Collimator type must be first parameter in collimator param file!\n");
					break;
				}
				
				switchOkay = true;
				switch (whichParam) {


					case ColEn_collimator_type:
					
					/* See if they gave us an enum or an int */
					if (paramType == LbPfEnEnum) {
						/* Search table for given label */
						for (typeIndex = 0; typeIndex < (NUM_COL_TYPES); typeIndex++){
						
							/* See if it matches */
							if (strcmp((char *)paramBuffer, colEn_ColTypeStr[typeIndex]) == 0) {
								ColRunTimeParams[ColCurParams].ColType = (ColEn_CollimatorTypeTy) typeIndex;
								break;
							}
						}
						if (typeIndex == NUM_COL_TYPES) {
							LbInPrintf("Invalid collimator type in collimator parameters '%s', valid types are:\n",
								(char *)paramBuffer);
							for (typeIndex = 1; typeIndex < (NUM_COL_TYPES); typeIndex++){
								LbInPrintf("'%s'\n", colEn_ColTypeStr[typeIndex]);
							}
							ErStGeneric("An invalid collimator type was supplied");
							switchOkay = false;
							}
						}
						else if (paramType == LbPfEnInteger) {
							if (((*(LbFourByte *)paramBuffer) != 0) && ((*(LbFourByte *)paramBuffer) < NUM_COL_TYPES)) {
								ColRunTimeParams[ColCurParams].ColType =
									(ColEn_CollimatorTypeTy) *((LbUsFourByte *) paramBuffer);
							}
							else {
								ErStGeneric("Error in collimator parameter file, you specified an invalid collimator type");
								switchOkay = false;
							}
						}
						else {
							ErStGeneric("Error in collimator parameter file, you specified an invalid type for the collimator parameter, it must be either INT or ENUM");
							switchOkay = false;
						}
						break;

					case ColEn_collimator_depth:
							ColRunTimeParams[ColCurParams].SimplePETCol.Depth =
								*((double *) paramBuffer);
							break;

					case ColEn_layers_list:
						switchOkay = ColGetMontCarloPETGeom(colParamFileHk,
							*((LbUsFourByte *) paramBuffer));
						break;
						

					case ColEn_history_file:
							strcpy(ColRunTimeParams[ColCurParams].ColHistoryFilePath,
								(char *) paramBuffer);
							
							ColRunTimeParams[ColCurParams].DoHistory = 
								(ColRunTimeParams[ColCurParams].ColHistoryFilePath[0] != '\0');
						break;
						

					case ColEn_history_params_file:
							strcpy(ColRunTimeParams[ColCurParams].ColHistoryParamsFilePath,
								(char *) paramBuffer);
							
							ColRunTimeParams[ColCurParams].DoCustomHistory = 
								(ColRunTimeParams[ColCurParams].ColHistoryParamsFilePath[0] != '\0');
						break;
						

					case ColEn_hole_geometry:
					
						/* See if they gave us an enum or an int */
						if (paramType == LbPfEnEnum) {
							/* Search table for given label */
							for (typeIndex = 0; typeIndex < (NUM_HOLE_TYPES); typeIndex++){
							
								/* See if it matches */
								if (strcmp((char *)paramBuffer, colEn_ColHoleTypeStr[typeIndex]) == 0) {
									ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleGeometry = (ColEn_HoleTy) typeIndex;
									break;
								}
							}
							if (typeIndex == NUM_HOLE_TYPES) {
								LbInPrintf("Invalid collimator hole type in collimator parameters '%s', valid types are:\n",
									(char *)paramBuffer);
								for (typeIndex = 1; typeIndex < (NUM_HOLE_TYPES); typeIndex++){
									LbInPrintf("'%s'\n", colEn_ColHoleTypeStr[typeIndex]);
								}
								ErStGeneric("An invalid collimator hole type was supplied");
								switchOkay = false;
								}
							}
							else if (paramType == LbPfEnInteger) {
								if (((*(LbFourByte *)paramBuffer) != 0) && ((*(LbFourByte *)paramBuffer) < NUM_HOLE_TYPES)) {
								ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleGeometry =
									(ColEn_HoleTy) *((LbUsFourByte *) paramBuffer);
								}
								else {
									ErStGeneric("Error in collimator parameter file, you specified an invalid collimator hole type");
									switchOkay = false;
								}
							}
							else {
								ErStGeneric("Error in collimator parameter file, you specified an invalid type for the collimator parameter, it must be either INT or ENUM");
								switchOkay = false;
							}
							break;

							ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleGeometry =
								(ColEn_HoleTy) *((LbUsFourByte *) paramBuffer);
						break;
						

					case ColEn_focal_length:
							ColRunTimeParams[ColCurParams].UNCSPECTCol.FocalLength =
								*((double *) paramBuffer);
						break;
						

					case ColEn_radius_of_rotation:
							ColRunTimeParams[ColCurParams].UNCSPECTCol.RadiusOfRotation =
								*((double *) paramBuffer);
						break;
						

					case ColEn_thickness:
							ColRunTimeParams[ColCurParams].UNCSPECTCol.Thickness =
								*((double *) paramBuffer);
						break;
						

					case ColEn_hole_radius:
							ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleRadius =
								*((double *) paramBuffer);
						break;
						

					case ColEn_septal_thickness:
							ColRunTimeParams[ColCurParams].UNCSPECTCol.SeptalThickness =
								*((double *) paramBuffer);
						break;
						

					case ColEn_min_z:
							ColRunTimeParams[ColCurParams].UNCSPECTCol.MinZ =
								*((double *) paramBuffer);
						break;
						

					case ColEn_max_z:
							ColRunTimeParams[ColCurParams].UNCSPECTCol.MaxZ =
								*((double *) paramBuffer);
						break;
						

					case ColEn_start_angle:
							ColRunTimeParams[ColCurParams].UNCSPECTCol.StartAngle =
								*((double *) paramBuffer);
						break;
						

					case ColEn_stop_angle:
							ColRunTimeParams[ColCurParams].UNCSPECTCol.StopAngle =
								*((double *) paramBuffer);
						break;
						

					case ColEn_sum_all_views:
						ErAlert("Sum_all_views is no longer used in the collimator module."
							"  The simultion will ignore this value, you may want to remove it"
							" from the collimator parameter file to avoid this message", false);
						break;
						

					case ColEn_num_views:
							ColRunTimeParams[ColCurParams].UNCSPECTCol.NumViews =
								*((LbUsFourByte *) paramBuffer);
							
						break;

					case ColEn_num_slats:
							ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].NumSlats =
								*((LbUsFourByte *) paramBuffer);
							
							/* Allocate memory for slats */
							if ((ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].Slats =
									(Col_Slat_Slat_Ty *) LbMmAlloc(
									ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].NumSlats *
									sizeof(Col_Slat_Slat_Ty))) == 0) {
								
								switchOkay = false;
							}
							
							curSlat = 0;
							specifiedNumSlats = true;
							
						break;

					case ColEn_num_layers:
							ColRunTimeParams[ColCurParams].SlatCol.NumLayers =
								*((LbUsFourByte *) paramBuffer);
							
							/* Allocate memory for layers */
							if ((ColRunTimeParams[ColCurParams].SlatCol.Layers =
									(Col_Slat_Layer_Ty *) LbMmAlloc(
									ColRunTimeParams[ColCurParams].SlatCol.NumLayers *
									sizeof(Col_Slat_Layer_Ty))) == 0) {
								
								switchOkay = false;
							}
						break;

					case ColEn_depth:
							ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].Depth =
								*((double *) paramBuffer);
							
						break;

					case ColEn_start:
					
							/* Verify that the number of slats for this layer has been specified */
							if (specifiedNumSlats == false) {
								sprintf(colParamErrString, "You have not specified the number of slats for layer %ld in your collimator parameters file. (ColGetRunTimeParams) \n",
		                            (long)curLayer);
		                 
								ErStGeneric(colParamErrString);
								switchOkay = false;
								break;
							}
							
							ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].Slats[curSlat].Start =
								*((double *) paramBuffer);
							
						break;

					case ColEn_end:
							ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].Slats[curSlat].End =
								*((double *) paramBuffer);
						break;

					case ColEn_material:
							ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].Slats[curSlat].Material =
								*((LbFourByte *) paramBuffer);
							curSlat += 1;
						break;
					
					case ColEn_inner_radius:
							/* Clear specified number of slats flag to be sure it is done
							on this layer */
							specifiedNumSlats = false;
							
							/* Increment layer count */
							curLayer += 1;
							
							/* Record inner radius */
							ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].InnerRadius =
								*((double *) paramBuffer);
							
						break;
					
					case ColEn_outer_radius:
							ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].Depth =
							*((double *) paramBuffer) -
							ColRunTimeParams[ColCurParams].SlatCol.Layers[curLayer].InnerRadius;
							
						break;


					case ColEn_NULL:
						sprintf(colParamErrString, "(ColGetRunTimeParams) Unknown (hence unused) parameter (%s).\n",
                            paramLabel);
						ErAlert(colParamErrString, false);
						break;
					
					default:
						sprintf(colParamErrString, "(ColGetRunTimeParams) Unknown (hence unused) parameter (%s).\n",
                            paramLabel);
						ErAlert(colParamErrString, false);
						break;
								}
				if (!switchOkay)
					break;
		}
		/* Close the parameter file */
		LbPfClose(&colParamFileHk);
		
		/* See if we quit due to error */
		if (!isEOF || !switchOkay) {
			/* Let them know where/when we failed */
			sprintf(colParamErrString, "Failed to process collimator parameter file:\n"
				"\t%s (ColGetRunTimeParams)\n",PhgRunTimeParams.PhgCollimatorParamsFilePath[ColCurParams]);
			ErAlert(colParamErrString, false);
			
			break;
		}
		okay = true;
	} while (false);
	
	return (okay);
}
/*********************************************************************************
*
*			Name:			ColLookupRunTimeParamLabel
*
*			Summary:	Find the given label in the label table and return
*						its enumeration value.
*			Arguments:
*				char		*label	- The label to find.
*			Function return: Enumerated type corresponding to label.
*
*********************************************************************************/
ColEn_RunTimeParamsTy  ColLookupRunTimeParamLabel(char *label)	
{
	ColEn_RunTimeParamsTy	whichParam;			/* Parameter we find */
	LbUsFourByte			labelIndex;			/* Label we are checking */

	/* Set param to null */
	whichParam = ColEn_NULL;
	
	/* Search table for given label */
	for (labelIndex = 0; labelIndex < (ColEn_NULL); labelIndex++){
	
		/* See if it matches */
		if (strcmp(label, colRunTimeParamLabels[labelIndex]) == 0) {
			whichParam = (ColEn_RunTimeParamsTy) labelIndex;
			break;
		}
	}
	
	return (whichParam);
}
#undef COL_PARAMS
