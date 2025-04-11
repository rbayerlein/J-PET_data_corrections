/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1994-2012 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/

/*********************************************************************************
*
*			Module Name:		Collimator.c
*			Revision Number:	2.2
*			Date last revised:	17 September 2012
*			Programmer:			Steven Vannoy
*			Date Originated:	11 July, 1994
*
*			Module Overview:	Simulates collimator functionality.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:
*				ColInitialize
*				ColGtOutsideRadius
*				ColPETPhotons
*				ColSPECTPhotons
*				ColTerminate
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
*
*			Revision description:
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		23 January 2012
*
*			Revision description:	Changed form of user functions to pointers
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*
*			Revision description:
*						- Support for simulate_PET_coincidences_only option
*						- Modified for compatibility with block detectors.
*						- Added initialization for more variables.
*						- Eight-byte integer for number of decays support.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		17 July 2004
*
*			Revision description:
*						- Changed ColPETPhotons to zero the number of photons
*						passing the collimator.  The lack of this was causing
*						errors when processing list mode data (history files).
*						- Changed colMCPETFindInnerDist to use the quadratic
*						solver in PhgMath.  It is more robust than the one that
*						was implemented in this routine.
*
*********************************************************************************/

#define COLLIMATOR


#include <stdio.h>
#include <memory.h>

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


/* Local types */
typedef enum { colEnAc_Null,
				colEnAc_Detect,
				colEnAc_Discard,
				colEnAc_Interact,
				colEnAc_Absorb,
				colEnAc_AxialCross,
				colEnAc_LayerCross} ColEn_ActionTy;

/* The following strings are used by the parameter code to determine which type
	of collimator segment is being used.
*/
typedef enum { null, parallel, tapered} colEn_SegTypeTy;

static char *colEn_ColSegTypeStr[] = {
	"null",
	"parallel",
	"tapered"};
#define NUM_SEG_TYPES 3


/* Prototypes */
void			colSimplePET(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					CollimatedPhotonsTy *colPhotonsPtr);
void			colMonteCarloPET(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					CollimatedPhotonsTy *colPhotonsPtr);
void 			colSimpleSPECT(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
					CollimatedPhotonsTy *colPhotonsPtr);
LbFourByte		colMCPETFindAxialSeg(double zPos, LbFourByte curLayer);
ColEn_ActionTy	colMCPETTrack(PHG_TrackingPhoton *photonPtr);
double			colMCPETFindInnerDist(PHG_Position *posPtr,
					PHG_Direction *dirPtr);
ColEn_ActionTy colMCPETProject(PHG_TrackingPhoton *photonPtr,
					Col_Monte_Carlo_PET_SegTy *curSegInfo, Boolean firstTime, 
					PHG_Position *newPosPtr, LbFourByte *curSegPtr,
					LbFourByte *curLayerPtr, double *distPtr);
double			colMCPETFindWallIntersect(PHG_Position *posPtr, PHG_Direction *dirPtr,
					double innerR, double outerR, double innerZ, double outerZ);
void			colMonteCarloSPECT(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
					CollimatedPhotonsTy *colPhotonsPtr);
void			colInitCylinders(LbUsFourByte curLayer);
					
		

/* Local Global variables */
static char							colErrStr[1024];					/* Storage for creating error strings */

/*********************************************************************************
*
*			Name:		ColGetMontCarloPETGeom
*
*			Summary:	Get the monte carlo PET collimator geometry from the 
*						parameter file. This routine assumes the parameter file
*						contains a description of a monte carlo PET collimator,
*						error handling is minimal.
*						NOTE THIS ROUTINE IS RECURSIVE.
*
*			Arguments:
*				LbPfHkTy		paramFlHk	- The parameter file hook.
*				LbUsFourByte	numParams	- The number of parameters.
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean	ColGetMontCarloPETGeom(LbPfHkTy paramFlHk, LbUsFourByte numParams)
{
	double					paramBuffer[LBPF_PARAM_LEN];	/* Parameter buffer */
	Boolean					okay = false;					/* Process flag */
	Boolean					isEOF;							/* End of file flag */
	Boolean					switchOkay;						/* Switch flag */
	LbPfEnPfTy				paramType;						/* Parameter type */
	LbUsFourByte			paramSize;						/* Parameter size */
	LbUsFourByte			typeIndex;						/* Loop control */
	char					paramLabel[LBPF_LABEL_LEN];		/* Parameter label space */
	char					errString[256];					/* Storage for error strings */
	ColEn_RunTimeParamsTy	whichParam;						/* Also current parameter */
	LbUsFourByte			paramIndex;						/* Current parameter */
	
	do { /* Process Loop */
	
		/* Allocate the collimator if this is the first time in */
		if (colData[ColCurParams].colCurLayer == -1 ) {
			if ((ColRunTimeParams[ColCurParams].MCPETCol = 
					(Col_Monte_Carlo_PET_Ty *)
					LbMmAlloc(sizeof(Col_Monte_Carlo_PET_Ty) *
					numParams))
					== 0) {
				
				break;
			}
			
			/* Store the number of layers it has */
			ColRunTimeParams[ColCurParams].MCPETCol[0].NumLayers = numParams;
		
			/* Set our current layer to zero */
			colData[ColCurParams].colCurLayer = 0;
		}
		
		/* Loop through the parameters at this level */
		for (paramIndex = 0; paramIndex < numParams; paramIndex++) {
			
			/* Clear process flag */
			okay = false;
		
			/* Get the current parameter */
			if (!LbPfGetParam(&paramFlHk, (void *)paramBuffer,
					&paramType, &paramSize, paramLabel, &isEOF)) {
				
				break;
			}
			
			/* Process the switch */
			whichParam = ColLookupRunTimeParamLabel(paramLabel);
			switchOkay = true;
			switch (whichParam) {
				
				case ColEn_segment_list:
				
				
					ColRunTimeParams[ColCurParams].MCPETCol[colData[ColCurParams].colCurLayer].NumSegments
						= *((LbUsFourByte *) paramBuffer);
					
			
					/* Allocate the collimator segments */
					if ((ColRunTimeParams[ColCurParams].MCPETCol[colData[ColCurParams].colCurLayer].Segments
							= (Col_Monte_Carlo_PET_SegTy *)
							LbMmAlloc(sizeof(Col_Monte_Carlo_PET_SegTy) *
							ColRunTimeParams[ColCurParams].MCPETCol[colData[ColCurParams].colCurLayer].NumSegments))
							== 0) {
						
						switchOkay = false;
					}
					
					/* Clear the current segment */
					colData[ColCurParams].colCurSeg = 0;
					
					/* Call this routine to process the segments */
					switchOkay = ColGetMontCarloPETGeom(paramFlHk,
							*((LbUsFourByte *) paramBuffer));
						
					/* Increment to Layer */
					colData[ColCurParams].colCurLayer++;
					
					break;
					
				case ColEn_segment:

					/* Set local current segment, this cuts down on the long data type
						references
					*/
					colData[ColCurParams].colCurMCPETSeg =
						&(ColRunTimeParams[ColCurParams].MCPETCol[colData[ColCurParams].colCurLayer].Segments[colData[ColCurParams].colCurSeg]);
						
					/* Call this routine to process the segments */
					switchOkay = ColGetMontCarloPETGeom(paramFlHk,
							*((LbUsFourByte *) paramBuffer));

					/* Set local current segment, this cuts down on the long data type
						references */
						
					/* Increment to next segment */
					colData[ColCurParams].colCurSeg++;
					
					break;

				case ColEn_seg_type:
					/* See if they gave us an enum or an int */
					if (paramType == LbPfEnEnum) {
						/* Search table for given label */
						for (typeIndex = 0; typeIndex < (NUM_SEG_TYPES); typeIndex++){
						
							/* See if it matches */
							if (strcmp((char *)paramBuffer, colEn_ColSegTypeStr[typeIndex]) == 0) {
								colData[ColCurParams].colCurMCPETSeg->SegType = (ColEn_SegTy) typeIndex;
								break;
							}
						}
						if (typeIndex == NUM_SEG_TYPES) {
							LbInPrintf("Invalid collimator segment type in collimator parameters '%s', valid types are:\n",
								(char *)paramBuffer);
							for (typeIndex = 1; typeIndex < (NUM_SEG_TYPES); typeIndex++){
								LbInPrintf("'%s'\n", colEn_ColSegTypeStr[typeIndex]);
							}
							ErStGeneric("An invalid collimator segment type was supplied");
							switchOkay = false;
							}
						}
						else if (paramType == LbPfEnInteger) {
							if (((*(LbFourByte *)paramBuffer) != 0) && ((*(LbFourByte *)paramBuffer) < NUM_SEG_TYPES)) {
							colData[ColCurParams].colCurMCPETSeg->SegType = 
									(ColEn_SegTy) *((LbUsFourByte *) paramBuffer);
							}
							else {
								ErStGeneric("Error in collimator parameter file, you specified an invalid collimator segment type");
								switchOkay = false;
							}
						}
						else {
							ErStGeneric("Error in collimator parameter file, you specified an invalid type for the collimator segment parameter, it must be either INT or ENUM");
							switchOkay = false;
						}
						break;

					break;

				case ColEn_material:
					colData[ColCurParams].colCurMCPETSeg->MaterialIndex = 
						*((LbUsFourByte *) paramBuffer);
					break;
					
				case ColEn_inner_radius:
					colData[ColCurParams].colCurMCPETSeg->InnerRadius = 
						*((double *) paramBuffer);
					if ( colData[ColCurParams].colCurMCPETSeg->InnerRadius < colData[ColCurParams].colInnermostRadius ) {
						colData[ColCurParams].colInnermostRadius = colData[ColCurParams].colCurMCPETSeg->InnerRadius;
					}
					break;
					
				case ColEn_outer_radius:
					colData[ColCurParams].colCurMCPETSeg->OuterRadius = 
						*((double *) paramBuffer);
					if ( colData[ColCurParams].colCurMCPETSeg->OuterRadius > colData[ColCurParams].colOutermostRadius ) {
						colData[ColCurParams].colOutermostRadius = colData[ColCurParams].colCurMCPETSeg->OuterRadius;
					}
					break;

				case ColEn_inner_min_z:
					colData[ColCurParams].colCurMCPETSeg->InnerMinZ = 
						*((double *) paramBuffer);
						
					/* If seg type is 1 then OuterMinZ the same */
					if (colData[ColCurParams].colCurMCPETSeg->SegType == parallel) {
						colData[ColCurParams].colCurMCPETSeg->OuterMinZ = colData[ColCurParams].colCurMCPETSeg->InnerMinZ;
					}
					break;

				case ColEn_outer_min_z:

					colData[ColCurParams].colCurMCPETSeg->OuterMinZ = 
						*((double *) paramBuffer);

					/* If seg type = 1 then this must be the same as inner min Z */
					if ((colData[ColCurParams].colCurMCPETSeg->SegType == parallel) &&
						(colData[ColCurParams].colCurMCPETSeg->InnerMinZ != colData[ColCurParams].colCurMCPETSeg->OuterMinZ)) {
						sprintf(errString, "When the collimator segment type is parallel, the inner and outer Z dimensions must be the same!");
						switchOkay = false;
						break;
					}
					break;

				case ColEn_inner_max_z:
					colData[ColCurParams].colCurMCPETSeg->InnerMaxZ = 
						*((double *) paramBuffer);
	
						
					/* If seg type is 1 then OuterMinZ the same */
					if (colData[ColCurParams].colCurMCPETSeg->SegType == parallel) {
						colData[ColCurParams].colCurMCPETSeg->OuterMaxZ = colData[ColCurParams].colCurMCPETSeg->InnerMaxZ;
					}
				break;

				case ColEn_outer_max_z:
					colData[ColCurParams].colCurMCPETSeg->OuterMaxZ = 
						*((double *) paramBuffer);
						
					/* If seg type = 1 then this must be the same as inner min Z */
					if ((colData[ColCurParams].colCurMCPETSeg->SegType == parallel) &&
						(colData[ColCurParams].colCurMCPETSeg->InnerMaxZ != colData[ColCurParams].colCurMCPETSeg->OuterMaxZ)) {
						sprintf(errString, "When the collimator segment type is parallel, the inner and outer Z dimensions must be the same!");
						switchOkay = false;
					}
					break;

				
				default:
					sprintf(errString, "(ColGetMontCarloPETGeom) Unknown (hence unused) parameter (%s).\n",
						paramLabel);
					ErAlert(errString, false);
					break;
			}
			if (!switchOkay)
				break;
				
			okay = true;
		}
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			ColPrintParams
*
*			Summary:		Prints the collimation parameters.
*
*			Arguments:
*			Function return: None.
*
*********************************************************************************/
void ColPrintParams()
{
	LbUsFourByte	i;				/* Index for looping */
	LbFourByte		j;				/* Index for looping */
	LbFourByte		currentColNum;	/* Index for looping through the collimator parameters */
	LbFourByte		ColCurParamsSave;	/* for saving the current collimator parameters so
											that this function leaves it unchanged */


	ColCurParamsSave = ColCurParams;
	
	for (currentColNum = 0; currentColNum < ColNumParams; currentColNum++) {
	
		ColCurParams = currentColNum;	/* this global must be set correctly for some macros
							and functions to work correctly */
	
		LbInPrintf("\n Collimator Parameters for '%s'\n", 
					PhgRunTimeParams.PhgCollimatorParamsFilePath[ColCurParams]);

		/* Print out parameters based on type of collimator */
		switch(ColRunTimeParams[ColCurParams].ColType) {
		
			case ColEn_simple_pet:
				LbInPrintf("\n\tCollimator is a simple PET model");
				break;
				
			case  ColEn_monte_carlo_pet:
				LbInPrintf("\n\tCollimator is a Monte Carlo PET model");
				LbInPrintf("\n\t\tThe number of layers is %d",
					ColRunTimeParams[ColCurParams].MCPETCol[0].NumLayers);
				for (i = 0; i < ColRunTimeParams[ColCurParams].MCPETCol[0].NumLayers; i++){
					LbInPrintf("\n\t\t\tFor layer %d", i);
					for (j = 0; (LbUsFourByte)j < ColRunTimeParams[ColCurParams].MCPETCol[i].NumSegments; j++){
						switch (ColRunTimeParams[ColCurParams].MCPETCol[i].Segments[j].SegType) {
						
							case parallel:
								LbInPrintf("\n\t\t\t\tSegment type is parallel");
								LbInPrintf("\n\t\t\t\tMax Z = %3.2f cm Min Z %3.2f",
										ColRunTimeParams[ColCurParams].MCPETCol[i].Segments[j].InnerMaxZ,
										ColRunTimeParams[ColCurParams].MCPETCol[i].Segments[j].InnerMinZ);
							break;
							
							case tapered:
							LbInPrintf("\n\t\t\t\tSegment type is tapered");
							LbInPrintf("\n\t\t\t\tInner Max Z = %3.2f cm Inner Min Z %3.2f",
								ColRunTimeParams[ColCurParams].MCPETCol[i].Segments[j].InnerMaxZ,
								ColRunTimeParams[ColCurParams].MCPETCol[i].Segments[j].InnerMinZ);
							LbInPrintf("\n\t\t\t\tOuter Max Z = %3.2f cm Outer Min Z %3.2f",
								ColRunTimeParams[ColCurParams].MCPETCol[i].Segments[j].OuterMaxZ,
								ColRunTimeParams[ColCurParams].MCPETCol[i].Segments[j].OuterMinZ);
							break;				
						
						default:
							break;
						}
						LbInPrintf("\n\t\t\t\tInner Radius = %3.2f cm Outer Radius %3.2f",
								ColRunTimeParams[ColCurParams].MCPETCol[i].Segments[j].InnerRadius,
								ColRunTimeParams[ColCurParams].MCPETCol[i].Segments[j].OuterRadius);
						LbInPrintf("\n\t\t\t\tMaterial is %s",
							SubObjGtAttenuationMaterialName(ColRunTimeParams[ColCurParams].MCPETCol[i].Segments[j].MaterialIndex));
					}
				}
					
				break;
			
			case ColEn_unc_spect:
				UNCColPrintParams( );
				break;
			
			case ColEn_slat:
				LbInPrintf("\n\tCollimator is a slat model");
				LbInPrintf("\n\t\tThe number of layers is %d",
					ColRunTimeParams[ColCurParams].SlatCol.NumLayers);
				for (i = 0; i < ColRunTimeParams[ColCurParams].SlatCol.NumLayers; i++){
					LbInPrintf("\n\t\t\tFor layer %d", i);
						LbInPrintf("\n\t\t\tInner radius = %3.2f cm  outer radius is %3.2f cm",
							ColRunTimeParams[ColCurParams].SlatCol.Layers[i].InnerRadius,
							ColRunTimeParams[ColCurParams].SlatCol.Layers[i].InnerRadius + ColRunTimeParams[ColCurParams].SlatCol.Layers[i].Depth);
					LbInPrintf("\n\t\t\tThere are %d slats:", ColRunTimeParams[ColCurParams].SlatCol.Layers[i].NumSlats);
					for (j = 0; j < colData[ColCurParams].colSlatNumSegs[i]; j++){
						LbInPrintf("\n\t\t\t\tStart = %3.2f cm End is %3.2f cm Material = %s",
							colData[ColCurParams].colSlatSegs[i][j].Start,
							colData[ColCurParams].colSlatSegs[i][j].End,
							SubObjGtAttenuationMaterialName(colData[ColCurParams].colSlatSegs[i][j].Material));
					}
					LbInPrintf("\n");				
				}
				break;
				
			default:
				PhgAbort("Invalid collimator type in ColRunTimeParams[ColCurParams] (ColPrintParams)",false);
				break;
				
			/* Print out history parameters */
			if (colData[ColCurParams].colHistFileHk.doCustom == true) {
				LbInPrintf("\nHistory file parameters for collimator module");
				PhoHFilePrintParams(&(colData[ColCurParams].colHistFileHk));
			}
		}
	}
	
	ColCurParams = ColCurParamsSave;
	
}

/*********************************************************************************
*
*			Name:			ColPrintReport
*
*			Summary:		Prints a report of the final statistics.
*
*			Arguments:
*			Function return: None.
*
*********************************************************************************/
void ColPrintReport()
{

	for (ColCurParams = 0; ColCurParams < ColNumParams; ColCurParams++) {
		/* Print out some collimator statistics */
		if (!ErIsInError()) {
			LbInPrintf("\n\n\n***************** Collimated Photon Statistics *************\n");
			
			/* Print out our statistics */

			if ( PHG_IsPET() ) {
				LbInPrintf("\nTotal blue photons reaching collimator in coincidence = %lld.", colData[ColCurParams].colTotBluePhotons);
				LbInPrintf("\nTotal pink photons reaching collimator in coincidence = %lld.", colData[ColCurParams].colTotPinkPhotons);
				LbInPrintf("\n\nTotal blue photons passing collimator = %lld.", colData[ColCurParams].colTotAccBluePhotons);
				LbInPrintf("\nTotal pink photons passing collimator = %lld.", colData[ColCurParams].colTotAccPinkPhotons);
				LbInPrintf("\n\nTotal blue photon weight reaching collimator = %3.4e.", colData[ColCurParams].colTotBluePhotonWt);
				LbInPrintf("\nTotal pink photon weight reaching collimator = %3.4e.", colData[ColCurParams].colTotPinkPhotonWt);
				LbInPrintf("\nTotal coincidence weight reaching collimator = %3.4e.", colData[ColCurParams].colTotInCoincWt);
				LbInPrintf("\n\nTotal blue photon weight passing collimator = %3.4e.", colData[ColCurParams].colTotAccBluePhotonWt);
				LbInPrintf("\nTotal pink photon weight passing collimator = %3.4e.", colData[ColCurParams].colTotAccPinkPhotonWt);
				LbInPrintf("\nTotal coincidence weight passing collimator = %3.4e.", colData[ColCurParams].colTotAccCoincWt);
				LbInPrintf("\n\nTotal blue photon weight rejected by collimator = %3.4e.", colData[ColCurParams].colTotRejBluePhotonWt);
				LbInPrintf("\nTotal pink photon weight rejected by collimator = %3.4e.", colData[ColCurParams].colTotRejPinkPhotonWt);
			}
			else {
				LbInPrintf("\nTotal photons reaching collimator in  = %lld.", colData[ColCurParams].colTotBluePhotons);
				LbInPrintf("\n\nTotal photons passing collimator = %lld.", colData[ColCurParams].colTotAccBluePhotons);
			}
			
			/* Print out parameters based on type of collimator */
			switch(ColRunTimeParams[ColCurParams].ColType) {
			
				case ColEn_simple_pet:
				case  ColEn_monte_carlo_pet:
					
					break;
				
				case ColEn_unc_spect:
					UNCColPrintReport();
					break;
					
				case ColEn_slat:
					LbInPrintf("\nTotal weight of object-scattered photons leaving collimator without interacting = %3.6e.",
								colData[ColCurParams].colTotScatWtPassedThroughCollimator);
					LbInPrintf("\nTotal weight of primary photons leaving collimator without interacting = %3.6e.",
								colData[ColCurParams].colTotPrimWtPassedThroughCollimator);
					break;
					
				default:
					PhgAbort("Invalid collimator type in ColRunTimeParams[ColCurParams] (ColPrintReport)",false);
					break;
					
			}
				
			/* Print out history parameters */
			if (colData[ColCurParams].colHistFileHk.doCustom == true) {
				LbInPrintf("\nHistory file report for collimator module");
				PhoHFilePrintReport(&(colData[ColCurParams].colHistFileHk));
			}

			LbInPrintf("\n\n****************************************************************\n");
		}

	}
	/* Make sure things go out */
	fflush(stdout);
}

/*********************************************************************************
*
*			Name:			ColInitialize
*
*			Summary:		Initialize the binning module.
*
*			Arguments:
*				Boolean		doHistory	- Do we create a history file (overides all)
*			Function return: TRUE unless an error occurs.
*
*********************************************************************************/
Boolean ColInitialize(Boolean doHistory)
{
	Boolean 		okay = false;	/* Process Flag */
	double			minZ;			/* Minimum z on collimator */
	double			maxZ;			/* Maximum z on collimator */
	do {
		
		/* Clear processing  variables */
		colData[ColCurParams].colCurLayer = -1;
		colData[ColCurParams].colTotBluePhotons = 0;
		colData[ColCurParams].colTotPinkPhotons = 0;
		colData[ColCurParams].colTotReachingCollimator = 0;
		colData[ColCurParams].colTotAccBluePhotons = 0;
		colData[ColCurParams].colTotAccPinkPhotons = 0;
		colData[ColCurParams].colTotScatWtPassedThroughCollimator = 0.0;
		colData[ColCurParams].colTotPrimWtPassedThroughCollimator = 0.0;
		colData[ColCurParams].colTotBluePhotonWt = 0.0;
		colData[ColCurParams].colTotPinkPhotonWt = 0.0;
		colData[ColCurParams].colTotInCoincWt = 0.0;
		colData[ColCurParams].colTotAccBluePhotonWt = 0.0;
		colData[ColCurParams].colTotAccPinkPhotonWt = 0.0;
		colData[ColCurParams].colTotAccCoincWt = 0.0;
		colData[ColCurParams].colTotRejBluePhotonWt = 0.0;
		colData[ColCurParams].colTotRejPinkPhotonWt = 0.0;
		
		/* initialize innermost and outermost radii */
		colData[ColCurParams].colInnermostRadius = LBFLOAT_MAX;
		colData[ColCurParams].colOutermostRadius = 0.0;
		
		/* Get the collimator parameters */
		if (ColGetRunTimeParams() == false) {
			break;
		}
	
		/* Compute the inner and outer most radii (and other info) */
		switch(ColRunTimeParams[ColCurParams].ColType) {
		
			case ColEn_simple_pet:
	
				/* Get the current target cylinder */
				colData[ColCurParams].colInnermostRadius = CylPosGetTargetRadius();
				
				colData[ColCurParams].colOutermostRadius = colData[ColCurParams].colInnermostRadius +
					ColRunTimeParams[ColCurParams].SimplePETCol.Depth;
				
				/* You are done */
				break;
				
			case  ColEn_monte_carlo_pet:
			
				/* Start by setting inner and outer radius to that of first segment */
				colData[ColCurParams].colCurInnerRadius = ColRunTimeParams[ColCurParams].MCPETCol[0].Segments[0].InnerRadius;
				colData[ColCurParams].colCurOuterRadius = ColRunTimeParams[ColCurParams].MCPETCol[0].Segments[0].OuterRadius;
				minZ = ColRunTimeParams[ColCurParams].MCPETCol[0].Segments[0].InnerMinZ;
				maxZ = ColRunTimeParams[ColCurParams].MCPETCol[0].Segments[0].InnerMaxZ;
				
				if (ColRunTimeParams[ColCurParams].MCPETCol[0].Segments[0].OuterMinZ < minZ)
					minZ = ColRunTimeParams[ColCurParams].MCPETCol[0].Segments[0].OuterMinZ;
					
				if (ColRunTimeParams[ColCurParams].MCPETCol[0].Segments[0].OuterMaxZ > maxZ)
					maxZ = ColRunTimeParams[ColCurParams].MCPETCol[0].Segments[0].OuterMinZ;
								
				/* Now initialize the cylinders */
				colData[ColCurParams].colInBoundCyl.radius = colData[ColCurParams].colCurInnerRadius;
				colData[ColCurParams].colInBoundCyl.zMin = minZ;
				colData[ColCurParams].colInBoundCyl.zMax = maxZ;
				colData[ColCurParams].colInBoundCyl.centerX = 0.0;
				colData[ColCurParams].colInBoundCyl.centerY = 0.0;

				colData[ColCurParams].colOutBoundCyl.radius = colData[ColCurParams].colCurOuterRadius;
				colData[ColCurParams].colOutBoundCyl.zMin = minZ;
				colData[ColCurParams].colOutBoundCyl.zMax = maxZ;
				colData[ColCurParams].colOutBoundCyl.centerX = 0.0;
				colData[ColCurParams].colOutBoundCyl.centerY = 0.0;
				
				break;
			
			case ColEn_unc_spect:
				if (!UNCColInitialize())
					goto FAIL;
				break;
			
			case ColEn_slat:
					/* None necessary here, this is primarily done from a routine named
						ColSlatSetParamsFromDetector which is called from the detector
						initialization
					 */
					break;
				
			default:
				PhgAbort("Invalid collimator type in ColRunTimeParams[ColCurParams] (ColInitialize)",false);
				break;
				
		}
		
		/* Override do history flag if they want it turned off */
		if (doHistory == false)
			ColRunTimeParams[ColCurParams].DoHistory = doHistory;
		
		/* Set our initialization flag */
		ColRunTimeParams[ColCurParams].Initialized = true;
		
		/* Call the user initialization routine */
		if (ColUsrInitializeFPtr) {
			(*ColUsrInitializeFPtr)(&ColRunTimeParams[ColCurParams]);
		}


		/* Create history file if were are supposed to */
		if (COL_IsDoHistory()) {
			if (PhoHFileCreate(ColRunTimeParams[ColCurParams].ColHistoryFilePath,
					ColRunTimeParams[ColCurParams].ColHistoryParamsFilePath,
					PhoHFileEn_COL, &(colData[ColCurParams].colHistFileHk)) == false) {
					
				sprintf(colErrStr,"Failed to create history file specified in collimator parameters file named:\n"
					"'%s'\n"
					"The custom parameters file name is: '%s'",
					ColRunTimeParams[ColCurParams].ColHistoryFilePath, ColRunTimeParams[ColCurParams].ColHistoryParamsFilePath);
				PhgAbort(colErrStr, false);
			}
		}
		
		okay = true;
		FAIL:;
	} while (false);
	
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			colInitCylinders
*
*			Summary:		Initialize our cylinders based on current
*							segment.
*
*			Arguments:
*				LbUsFourByte	curLayer			- The current layer
*			Function return: None.
*
*********************************************************************************/
void colInitCylinders(LbUsFourByte curLayer)
{
	double			minZ;			/* Minimum z on collimator */
	double			maxZ;			/* Maximum z on collimator */
	LbUsFourByte	segIndex;		/* Index for traversing segments */
	
	/* Start by setting inner and outer radius to that of first segment */
	colData[ColCurParams].colCurInnerRadius = ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[0].InnerRadius;
	colData[ColCurParams].colCurOuterRadius = ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[0].OuterRadius;
	minZ = ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[0].InnerMinZ;
	maxZ = ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[0].InnerMaxZ;
	
	if (ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[0].OuterMinZ < minZ)
		minZ = ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[0].OuterMinZ;
		
	if (ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[0].OuterMaxZ > maxZ)
		maxZ = ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[0].OuterMinZ;
	
	/* Loop through remaining segments to get largest and smallest values */
	for (segIndex = 1; segIndex < ColRunTimeParams[ColCurParams].MCPETCol[curLayer].NumSegments; segIndex++) {
	
		/* Get minimum/maximum inner/outer radii */
		if (ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].InnerRadius < 
				colData[ColCurParams].colCurInnerRadius) {
		
			colData[ColCurParams].colCurInnerRadius = ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].InnerRadius;
		}
		if (ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].OuterRadius > 
				colData[ColCurParams].colCurOuterRadius) {
		
			colData[ColCurParams].colCurOuterRadius = ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].OuterRadius;
		}
		
		/* Get minimum/maximum inner/outer Z ranges */
		if (ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].OuterMinZ < minZ)
			minZ = ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].OuterMinZ;
			
		if (ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].OuterMaxZ > maxZ)
			maxZ = ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].OuterMaxZ;

		if (ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].InnerMinZ < minZ)
			minZ = ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].InnerMinZ;
			
		if (ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].InnerMaxZ > maxZ)
			maxZ = ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].InnerMaxZ;
	}
	
	/* Now initialize the cylinders */
	colData[ColCurParams].colInBoundCyl.radius = colData[ColCurParams].colCurInnerRadius;
	colData[ColCurParams].colInBoundCyl.zMin = minZ;
	colData[ColCurParams].colInBoundCyl.zMax = maxZ;
	colData[ColCurParams].colInBoundCyl.centerX = 0.0;
	colData[ColCurParams].colInBoundCyl.centerY = 0.0;
	
	colData[ColCurParams].colOutBoundCyl.radius = colData[ColCurParams].colCurOuterRadius;
	colData[ColCurParams].colOutBoundCyl.zMin = minZ;
	colData[ColCurParams].colOutBoundCyl.zMax = maxZ;
	colData[ColCurParams].colOutBoundCyl.centerX = 0.0;
	colData[ColCurParams].colOutBoundCyl.centerY = 0.0;
}

/*********************************************************************************
*
*			Name:			ColGtOutsideRadius
*
*			Summary:		Returns the outside radius of the collimator. Originally
*							added for the binning module.
*
*			Arguments:
*			Function return: Outside radius of the collimator.
*
*********************************************************************************/
double	ColGtOutsideRadius()
{
	double outsideRadius;

	/* If we aren't doing run-time collimation just return the target radius */
	if (PHG_IsCollimateOnTheFly() == false) {
	
			outsideRadius = CylPosGetTargetRadius();
	}
	else {
		/* Access the appropriate data structure based on collimator type */
		switch(ColRunTimeParams[ColCurParams].ColType) {
		
			case ColEn_simple_spect:
			case ColEn_simple_pet:
	
				/* Get the current target cylinder */
				outsideRadius = CylPosGetTargetRadius();
				
				/* Extend the cylinder the requested amount */
				outsideRadius += ColRunTimeParams[ColCurParams].SimplePETCol.Depth;
				
				/* You are done */
				break;
			
			case ColEn_monte_carlo_pet:
				outsideRadius = colData[ColCurParams].colOutermostRadius;
				break;
				
			case ColEn_unc_spect:
				outsideRadius = (ColRunTimeParams[ColCurParams].UNCSPECTCol.RadiusOfRotation +
					ColRunTimeParams[ColCurParams].UNCSPECTCol.Thickness);
			break;
				
			case ColEn_slat:
				outsideRadius = (ColRunTimeParams[ColCurParams].SlatCol.Layers[ColRunTimeParams[ColCurParams].SlatCol.NumLayers-1].InnerRadius +
					ColRunTimeParams[ColCurParams].SlatCol.Layers[ColRunTimeParams[ColCurParams].SlatCol.NumLayers-1].Depth);
			break;
			
			default:
				PhgAbort("Invalid collimator type in ColRunTimeParams[ColCurParams] (ColGtOutsideRadius)",false);
				break;
				
		}
	}
	return (outsideRadius);
}

/*********************************************************************************
*
*			Name:			ColGtInsideRadius
*
*			Summary:		Returns the inside radius of the collimator. 
*			Arguments:
*			Function return: Inside radius of the collimator.
*
*********************************************************************************/
double	ColGtInsideRadius()
{
	double insideRadius;

	/* If we aren't doing run-time collimation just return the target radius */
	if (PHG_IsCollimateOnTheFly() == false) {
	
			insideRadius = CylPosGetTargetRadius();
	}
	else {
		/* Access the appropriate data structure based on collimator type */
		switch(ColRunTimeParams[ColCurParams].ColType) {
		
			case ColEn_simple_spect:
			case ColEn_simple_pet:
	
				/* Get the current target cylinder */
				insideRadius = CylPosGetTargetRadius();
				break;
			
			case ColEn_monte_carlo_pet:
				insideRadius = colData[ColCurParams].colInnermostRadius;
				break;
				
			case ColEn_unc_spect:
				insideRadius = ColRunTimeParams[ColCurParams].UNCSPECTCol.RadiusOfRotation;
			break;
					
			case ColEn_slat:
				insideRadius = ColRunTimeParams[ColCurParams].SlatCol.Layers[0].InnerRadius;
				break;
			
			default:
				PhgAbort("Invalid collimator type in ColRunTimeParams[ColCurParams] (ColGtInsideRadius)",false);
				break;
				
		}
	}
	return (insideRadius);
}

/*********************************************************************************
*
*			Name:			ColIsSlat
*
*			Summary:		Returns true if slat collimation is being done.
*
*			Arguments:
*			Function return: true if slat collimation is being done.
*
*********************************************************************************/
Boolean	ColIsSlat()
{
return((ColRunTimeParams[ColCurParams].ColType == ColEn_slat) ? true : false);
}

/*********************************************************************************
*
*			Name:			ColIsUNC
*
*			Summary:		Returns true if slat collimation is being done.
*
*			Arguments:
*			Function return: true if UNC collimation is being done.
*
*********************************************************************************/
Boolean	ColIsUNC()
{
return((ColRunTimeParams[ColCurParams].ColType == ColEn_unc_spect) ? true : false);
}

/*********************************************************************************
*
*			Name:			ColGtZLimits
*
*			Summary:		Returns the Z limits for the collimator.
*
*			Arguments:
*					double	*zMin	- Minimum Z limit
*					double	*zMax	- Maximum Z limit
*			Function return: None.
*
*********************************************************************************/
void	ColGtZLimits(double *zMin, double *zMax)
{

	#ifdef PHG_DEBUG
	/* If we aren't doing run-time collimation this is a bogus call */
	if (PHG_IsCollimateOnTheFly() == false) {
	
		PhgAbort("Called 'ColGtZLimits' when not doing collimation (ColGtZLimits)", false);
	}
	#endif

	/* Access the appropriate data structure based on collimator type */
	switch(ColRunTimeParams[ColCurParams].ColType) {
	
		case ColEn_simple_pet:

			PhgAbort("Call to ColGtZLimits for unsupported collimator type, ColEn_simple_pet (ColGtOutsideRadius)",false);
			break;
		
		case ColEn_monte_carlo_pet:
			PhgAbort("Call to ColGtZLimits for unsupported collimator type, ColEn_monte_carlo_pet (ColGtOutsideRadius)",false);
			break;
			
		case ColEn_unc_spect:
			*zMin = ColRunTimeParams[ColCurParams].UNCSPECTCol.MinZ;
			*zMax = ColRunTimeParams[ColCurParams].UNCSPECTCol.MaxZ;
		break;
		
		default:
			PhgAbort("Invalid collimator type in ColRunTimeParams[ColCurParams] (ColGtZLimits)",false);
			break;
			
	}
}

/*********************************************************************************
*
*			Name:			ColGtNumViews
*
*			Summary:		Returns the number of views specified for the SPECT
*							UNC Collimator.
*
*			Arguments:
*			Function return: Number of views.
*
*********************************************************************************/
LbUsFourByte	ColGtNumViews()
{
	LbUsFourByte	numViews;
	
	#ifdef PHG_DEBUG
		/* If we aren't doing run-time collimation bolt */
		if (PHG_IsCollimateOnTheFly() == false) {
			PhgAbort("You can't call this routine if you aren't collimating! (ColGtNumViews)",false);
		}
		
	#endif
	
		/* Access the appropriate data structure based on collimator type */
		switch(ColRunTimeParams[ColCurParams].ColType) {
			
			case ColEn_monte_carlo_pet:
			case ColEn_simple_spect:
			case ColEn_simple_pet:
			case ColEn_slat:
				/* We don't expect to get the call for these types of collimators */
				PhgAbort("Unexpected call given the type of collimator and detector (ColGtNumViews)",false);
				break;
				
			case ColEn_unc_spect:
				numViews = ColRunTimeParams[ColCurParams].UNCSPECTCol.NumViews;
			break;
				
			default:
				PhgAbort("Invalid collimator type in ColRunTimeParams[ColCurParams] (ColGtNumViews)",false);
				break;
				
		}
	
	/* Return number of views */
	return (numViews);
}

/*********************************************************************************
*
*			Name:			ColPETPhotons
*
*			Summary:		Update the binning images with the current batch of
*							detected photons.
*
*			Arguments:
*				PHG_Decay			decayPtr			- The decayPtr that started the process.
*				PHG_TrackingPhoton *bluePhotons		- The blue photons detected.
*				LbUsFourByte 		numBluePhotons	- The number of blue photons.
*				PHG_TrackingPhoton *pinkPhotons		- The  pink photons detected.
*				LbUsFourByte		numPinkPhotons	- The number of pink photons.
*				CollimatedPhotonsTy *colPhotonsPtr	- The accepted collimated photons
*			Function return: None.
*
*********************************************************************************/
void ColPETPhotons(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		CollimatedPhotonsTy *colPhotonsPtr)
{
	
	/* Make sure that the number of photons clearing collimator gets initialized to 0.  */
	colPhotonsPtr->NumCollimatedBluePhotons = 0;
	colPhotonsPtr->NumCollimatedPinkPhotons = 0;
	colData[ColCurParams].colDetectedPinkPhotonIndex = 0;
	colData[ColCurParams].colDetectedBluePhotonIndex = 0;

	/* Bolt if there is nothing to do */
	if ((numBluePhotons == 0) && (numPinkPhotons == 0)) {
		return;
	}
		

	/* Compute statistics */
	colData[ColCurParams].colTotBluePhotons += numBluePhotons;
	colData[ColCurParams].colTotPinkPhotons += numPinkPhotons;

	/* Call appropriate collimation routine based on user specified type */
	switch(ColRunTimeParams[ColCurParams].ColType) {
	
		case ColEn_simple_pet:
		
			/* Call simple pet collimator */
			colSimplePET(decayPtr, bluePhotons, numBluePhotons,
				pinkPhotons, numPinkPhotons,
				colPhotonsPtr);
				
			/* You are done */
			break;
			
		case ColEn_monte_carlo_pet:
			colMonteCarloPET(decayPtr,
				bluePhotons,  numBluePhotons,
				pinkPhotons,  numPinkPhotons,
				colPhotonsPtr);
			break;
				
		case ColEn_slat:
			ColSlatDualHeaded(decayPtr,
				bluePhotons,  numBluePhotons,
				pinkPhotons,  numPinkPhotons,
				colPhotonsPtr);
			break;
	
		default:
			PhgAbort("Invalid collimator type in ColRunTimeParams[ColCurParams] (ColPETPhotons)",false);
			break;
			
	}
	
	/* Update statistics */
	colData[ColCurParams].colTotAccBluePhotons += colPhotonsPtr->NumCollimatedBluePhotons; 
	colData[ColCurParams].colTotAccPinkPhotons += colPhotonsPtr->NumCollimatedPinkPhotons; 
	
	/* Write collimated photons to history file, if requested */
	if (COL_IsDoHistory()) {

		/* Write the photons */
		if (PhoHFileWriteDetections(&(colData[ColCurParams].colHistFileHk), decayPtr,
				colData[ColCurParams].colDetectedBluePhotons,
				colData[ColCurParams].colDetectedBluePhotonIndex,
				colData[ColCurParams].colDetectedPinkPhotons,
				colData[ColCurParams].colDetectedPinkPhotonIndex) == false) {
			
			/* Abort Program execution */
			PhgAbort("Got failure from PhoHFileWriteDetections trapped in ColPETPhotons.",false);
		}
	}

}

/*********************************************************************************
*
*			Name:			colSimplePET
*
*			Summary:		Perform "simple" collimation. This routine simply
*							projects the photon the distance of the collimator
*							depth. If the photon is outside of the axial limits
*							it is rejected.
*
*			Arguments:
*				PHG_Decay			decayPtr			- The decayPtr that started the process.
*				PHG_TrackingPhoton *bluePhotons		- The blue photons detected.
*				LbUsFourByte 		numBluePhotons	- The number of blue photons.
*				PHG_TrackingPhoton *pinkPhotons		- The  pink photons detected.
*				LbUsFourByte		numPinkPhotons	- The number of pink photons.
*				CollimatedPhotonsTy *colPhotonsPtr	- The accepted collimated photons
*			Function return: None.
*
*********************************************************************************/
void colSimplePET(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		CollimatedPhotonsTy *colPhotonsPtr)
{
	PHG_Position		newPosition;				/* Projected location */
	double				distance;					/* Distance traveled to extended cylinder */
	CylPosCylinderTy	extendedCylinder;			/* The cylinder we will project to */	
	LbUsFourByte		photonIndex;				/* The current photon */
	LbUsFourByte		bluePhotonIndex;			/* The current blue photon when processing coincidences */

	do { /* Process Loop */
	
		/* Clear the counters */
		colPhotonsPtr->NumCollimatedBluePhotons = 0;
		colPhotonsPtr->NumCollimatedPinkPhotons = 0;
		colData[ColCurParams].colDetectedPinkPhotonIndex = 0;
		colData[ColCurParams].colDetectedBluePhotonIndex = 0;

		/* Get the current target cylinder */
		extendedCylinder = CylPosTargetCylinder;
		
		/* Extend the cylinder the requested amount */
		extendedCylinder.radius += ColRunTimeParams[ColCurParams].SimplePETCol.Depth;
		
		/* Collimate the blue photons */
		for (photonIndex = 0; photonIndex < numBluePhotons; photonIndex++){
			
			/* keep track of the incoming weight */
			colData[ColCurParams].colTotBluePhotonWt += 
					bluePhotons[photonIndex].photon_current_weight * decayPtr->startWeight;

			/* Let user modify and/or reject photons */
			if (ColUsrModPETPhotonsFPtr && 
					(*ColUsrModPETPhotonsFPtr)(&ColRunTimeParams[ColCurParams], decayPtr,
						&(bluePhotons[photonIndex])) == false) {
				
				/* keep track of rejected weight */
				colData[ColCurParams].colTotRejBluePhotonWt += 
						bluePhotons[photonIndex].photon_current_weight * decayPtr->startWeight;
				
				/* They rejected it so go to next blue */
				continue;
			}
			
			/* Project blue photon to cylinder, bolt if unsuccessful */
			if (CylPosProjectToCylinder(&(bluePhotons[photonIndex].location),
					&(bluePhotons[photonIndex].angle), &extendedCylinder,
					&newPosition, &distance) == false) {
			
				/* keep track of rejected weight */
				colData[ColCurParams].colTotRejBluePhotonWt += 
						bluePhotons[photonIndex].photon_current_weight * decayPtr->startWeight;
				
				/* photon is rejected by collimator, go to next blue */
				continue;
			}
				
					
			/* See if we fall outside of axial limits */
			if ((newPosition.z_position > extendedCylinder.zMax) ||
					(newPosition.z_position < extendedCylinder.zMin)) {

				/* keep track of rejected weight */
				colData[ColCurParams].colTotRejBluePhotonWt += 
						bluePhotons[photonIndex].photon_current_weight * decayPtr->startWeight;
				
				/* photon is rejected by collimator, go to next blue */
				continue;
			}
			else {
				/* Save the new location */
				bluePhotons[photonIndex].location = newPosition;
				bluePhotons[photonIndex].travel_distance += distance;

				/* add this weight to the total accepted blue photon weight */
				colData[ColCurParams].colTotAccBluePhotonWt += 
						bluePhotons[photonIndex].photon_current_weight * decayPtr->startWeight;
			}

			/* We made it to here so save the collimated photons */
			{
				colPhotonsPtr->CollimatedTrkngBluePhotons[colPhotonsPtr->NumCollimatedBluePhotons]
					= bluePhotons[photonIndex];

				/* Increment the counters */
				colPhotonsPtr->NumCollimatedBluePhotons++;

				/* Update our detected photon block if doing history file */
				if  (COL_IsDoHistory()) {
					ColUpdateCollimatedPhotonBlock(&bluePhotons[photonIndex]);
				}
			}
		}
		
		/* Collimate the pink photons */
		for (photonIndex = 0; photonIndex < numPinkPhotons; photonIndex++){

			/* keep track of incoming photon weight */
			colData[ColCurParams].colTotPinkPhotonWt += 
					pinkPhotons[photonIndex].photon_current_weight * decayPtr->startWeight;

			/* also keep track of the total incoming coincidence weight */
			for (bluePhotonIndex = 0; bluePhotonIndex < numBluePhotons; bluePhotonIndex++){
				colData[ColCurParams].colTotInCoincWt += pinkPhotons[photonIndex].photon_current_weight *
					bluePhotons[bluePhotonIndex].photon_current_weight * decayPtr->startWeight;
			}
			
			/* Let user modify and/or reject photons */
			if (ColUsrModPETPhotonsFPtr && 
					(*ColUsrModPETPhotonsFPtr)(&ColRunTimeParams[ColCurParams], decayPtr,
						&(pinkPhotons[photonIndex])) == false) {
				
				/* keep track of rejected weight */
				colData[ColCurParams].colTotRejPinkPhotonWt += 
						pinkPhotons[photonIndex].photon_current_weight * decayPtr->startWeight;

				/* They rejected it so go to next pink */
				continue;
			}

			/* Project pink photon to cylinder, bolt if unsuccessful */
			if (CylPosProjectToCylinder(&(pinkPhotons[photonIndex].location),
					&(pinkPhotons[photonIndex].angle), &extendedCylinder,
					&newPosition, &distance) == false) {

				/* keep track of rejected weight */
				colData[ColCurParams].colTotRejPinkPhotonWt += 
						pinkPhotons[photonIndex].photon_current_weight * decayPtr->startWeight;

				/* Collimator rejected it so go to next pink */
				continue;
			}
					
			/* See if we fall outside of axial limits */
			if ((newPosition.z_position > extendedCylinder.zMax) ||
					(newPosition.z_position < extendedCylinder.zMin)) {
				
				/* keep track of rejected weight */
				colData[ColCurParams].colTotRejPinkPhotonWt += 
						pinkPhotons[photonIndex].photon_current_weight * decayPtr->startWeight;

				/* Collimator rejected it so go to next pink */
				continue;
			}
			else {
				/* Save the new location */
				pinkPhotons[photonIndex].location = newPosition;
				pinkPhotons[photonIndex].travel_distance += distance;

				/* track total pink weight that clears collimator */
				colData[ColCurParams].colTotAccPinkPhotonWt += 
						pinkPhotons[photonIndex].photon_current_weight * decayPtr->startWeight;

				/* also keep track of the total coincidence weight clearing collimator */
				for (bluePhotonIndex = 0; bluePhotonIndex < numBluePhotons; bluePhotonIndex++){
					colData[ColCurParams].colTotAccCoincWt += pinkPhotons[photonIndex].photon_current_weight *
						bluePhotons[bluePhotonIndex].photon_current_weight * decayPtr->startWeight;
				}
			}

			/* We made it to here so save the collimated photons */
			{
				colPhotonsPtr->CollimatedTrkngPinkPhotons[colPhotonsPtr->NumCollimatedPinkPhotons]
					= pinkPhotons[photonIndex];

				/* Increment the counters */
				colPhotonsPtr->NumCollimatedPinkPhotons++;

				/* Update our detected photon block if doing history file */
				if  (COL_IsDoHistory()) {
					ColUpdateCollimatedPhotonBlock(&pinkPhotons[photonIndex]);
				}
			}
		}
	} while (false);
	
}

/*********************************************************************************
*
*			Name:			colMonteCarloPET
*
*			Summary:		Perform Monte Carlo collimation.
*
*			Arguments:
*				PHG_Decay			decayPtr			- The decayPtr that started the process.
*				PHG_TrackingPhoton *bluePhotons		- The blue photons detected.
*				LbUsFourByte 		numBluePhotons	- The number of blue photons.
*				PHG_TrackingPhoton *pinkPhotons		- The  pink photons detected.
*				LbUsFourByte		numPinkPhotons	- The number of pink photons.
*				CollimatedPhotonsTy *colPhotonsPtr	- The accepted collimated photons
*			Function return: None.
*
*********************************************************************************/
void colMonteCarloPET(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		CollimatedPhotonsTy *colPhotonsPtr)
{
	Boolean				doingBlue = true;	/* Flag for doing blue or pink photons */
	LbUsFourByte		photonIndex;		/* The current photon */
	LbUsFourByte		bluePhotonIndex;	/* The current blue photon for coincidence weight sub-loops */
	LbUsFourByte		loopIndex;			/* Loop index for blue/pink photons */
	PHG_TrackingPhoton	*photons;			/* Generic pointer for accessing both blue and pink photons */
	LbUsFourByte		numPhotons;			/* Generic count of blue or pink photons */
	ColEn_ActionTy		action;				/* What action to take during tracking */
	
	do { /* Process Loop */
	
		/* Clear the counters */
		colPhotonsPtr->NumCollimatedBluePhotons = 0;
		colPhotonsPtr->NumCollimatedPinkPhotons = 0;
		colData[ColCurParams].colDetectedPinkPhotonIndex = 0;
		colData[ColCurParams].colDetectedBluePhotonIndex = 0;
		#ifdef PHG_DEBUG
		ColDiscardBlue[ColCurParams] = true;
		ColDiscardPink[ColCurParams] = true;
		#endif
		/* Set our generic pointers */
		photons = bluePhotons;
		numPhotons = numBluePhotons;
		
		
		/* We use a loop to process the blues and then the pinks. This
			minimizes code duplication.
		*/
		for (loopIndex = 0; loopIndex < 2; loopIndex++) {
			/* Collimate the photons */
			for (photonIndex = 0; photonIndex < numPhotons; photonIndex++){
			
				/* keep track of the total incoming blue and pink weight */
				if (doingBlue == true) {
					colData[ColCurParams].colTotBluePhotonWt += 
						photons[photonIndex].photon_current_weight * decayPtr->startWeight;
				}
				else {
					colData[ColCurParams].colTotPinkPhotonWt += 
						photons[photonIndex].photon_current_weight * decayPtr->startWeight;
					
					/* also keep track of the total incoming coincidence weight */
					for (bluePhotonIndex = 0; bluePhotonIndex < numBluePhotons; bluePhotonIndex++){
						colData[ColCurParams].colTotInCoincWt += photons[photonIndex].photon_current_weight *
							bluePhotons[bluePhotonIndex].photon_current_weight * decayPtr->startWeight;
					}
				}
				
				/* Let user modify and/or reject photons */
				if (ColUsrModPETPhotonsFPtr && 
						(*ColUsrModPETPhotonsFPtr)(&ColRunTimeParams[ColCurParams], decayPtr,
							&(photons[photonIndex])) == false) {
					
					/* They rejected it so go to next photon */
					continue;
				}
				
				/* Track the photon */
				action = colMCPETTrack(&(photons[photonIndex]));
								
				/* Store the photon if the action is to detect */
				if (action == colEnAc_Detect) {	 

					/* We made it to here so save the collimated photons */
					if (doingBlue == true) {
					
						/* add this weight to the total accepted blue photon weight */
						colData[ColCurParams].colTotAccBluePhotonWt += 
								photons[photonIndex].photon_current_weight * decayPtr->startWeight;

#ifdef PHG_DEBUG
ColDiscardBlue[ColCurParams] = false;
#endif
						colPhotonsPtr->CollimatedTrkngBluePhotons[colPhotonsPtr->NumCollimatedBluePhotons]
							= bluePhotons[photonIndex];
		
						/* Increment the counters */
						colPhotonsPtr->NumCollimatedBluePhotons++;
		
						/* Update our detected photon block if doing history file */
						if  (COL_IsDoHistory()) {
							ColUpdateCollimatedPhotonBlock(&bluePhotons[photonIndex]);
						}
					}
					else {
						
						/* increment the counter that tracks total accepted pink weight */
						colData[ColCurParams].colTotAccPinkPhotonWt += 
								photons[photonIndex].photon_current_weight * decayPtr->startWeight;

						/* also keep track of the total coincidence weight clearing collimator */
						for (bluePhotonIndex = 0; bluePhotonIndex < numBluePhotons; bluePhotonIndex++){
							colData[ColCurParams].colTotAccCoincWt += photons[photonIndex].photon_current_weight *
								bluePhotons[bluePhotonIndex].photon_current_weight * decayPtr->startWeight;
						}

#ifdef PHG_DEBUG
ColDiscardPink[ColCurParams] = false;
#endif
						colPhotonsPtr->CollimatedTrkngPinkPhotons[colPhotonsPtr->NumCollimatedPinkPhotons]
							= pinkPhotons[photonIndex];
		
						/* Increment the counters */
						colPhotonsPtr->NumCollimatedPinkPhotons++;
		
						/* Update our detected photon block if doing history file */
						if  (COL_IsDoHistory()) {
							ColUpdateCollimatedPhotonBlock(&pinkPhotons[photonIndex]);
						}
					}
				}
				else {
				
					if (doingBlue == true) {
						colData[ColCurParams].colTotRejBluePhotonWt += 
								photons[photonIndex].photon_current_weight * decayPtr->startWeight;
					}
					else {
						colData[ColCurParams].colTotRejPinkPhotonWt += 
								photons[photonIndex].photon_current_weight * decayPtr->startWeight;
					}
				}				
			}
			
			/* if no blue photons were found and this is a coincidence-only
			simulation, break out of the loop */
			if ( PHG_IsPETCoincidencesOnly() && (colPhotonsPtr->NumCollimatedBluePhotons == 0) ) {
				
				break;
				
			}
	
			/* We are at the bottom of the loop so set our pointers to the
				pink photons
			*/
			photons = pinkPhotons;
			numPhotons = numPinkPhotons;
			doingBlue = false;
		}
	} while (false);
	
}

/*********************************************************************************
*
*			Name:			colMCPETTrack
*
*			Summary:		Perform the Monte Carlo tracking of photons
*							through the PET collimator.
*
*			Arguments:
*				PHG_TrackingPhoton	*photonPtr	- The photon.
*
*			Function return: The action to take based on the results of tracking.
*
*********************************************************************************/
ColEn_ActionTy colMCPETTrack(PHG_TrackingPhoton *photonPtr)
{
	Boolean						firstTime = true;		/* Flag indicating first leg of projection */
	double						distance;				/* Distance traveled to extended cylinder */
	double						innerRsquared;			/* Inner radius of collimator squared */
	double						rSquared;				/* Current "radial" position of photon */
	double						fpToGo;					/* Free paths to travel before interacting */
	double						attenuation;			/* The attenuation of the collimator material */
	double						comptonToScatterProbability;		/* Ratio of prob of compton to prob of scatter */
	double						interactionProbability;	/* Cumulative probability of compton and coherent scatter */
	double						scatterProbability;	/* Probability of absorption */
	LbUsFourByte				materialIndex;			/* Index into material for collimator */
	LbFourByte					curLayer;				/* Which layer we are in */
	LbFourByte					newLayer;				/* Which layer we are in */
	LbFourByte					curSeg;					/* The photon's segment */
	ColEn_ActionTy				action;					/* What action will result? */
	PHG_Position				newPosition;			/* The photon's new location */
	Col_Monte_Carlo_PET_SegTy	*curSegInfo;			/* The current MCPET collimator segment (if we are doing that type of collimator) */
	
	/* Assume action is null */
	action = colEnAc_Null;
	 
	do { /* Process Loop */
	
		/*	Initialize our cylinders, this must be done each time because it is
			dependent on the current layer
		*/
		colInitCylinders(0);
		
		/* Compute optimizing parameters */
		innerRsquared = PHGMATH_Square(colData[ColCurParams].colInBoundCyl.radius);

		/* Copy the current location */
		newPosition = photonPtr->location;
		
		/* Clear the distance traveled */
		distance = 0.0;
	
		/* The target cylinder of the PHG and the inner-most radius of the
			collimator may be different. Hence, we will first check to
			see if we are on the inner-most surface. If not, we will
			project to there
		*/
		{
			rSquared = PHGMATH_Square(photonPtr->location.x_position) +
				PHGMATH_Square(photonPtr->location.y_position);
				
		
			if (rSquared < innerRsquared) {
				
				/* Project photon to cylinder */
				if (CylPosProjectToCylinder(&(photonPtr->location),
						&(photonPtr->angle), &(colData[ColCurParams].colInBoundCyl),
						&newPosition, &distance) == false) {
				
					action = colEnAc_Discard;
					break;
				}
				
				/* See if we fall outside of axial limits */
				if ((newPosition.z_position > colData[ColCurParams].colInBoundCyl.zMax) ||
						(newPosition.z_position < colData[ColCurParams].colInBoundCyl.zMin)) {
					
					action = colEnAc_Discard;
					break;
				}
				
				/* We made it here so update our position on the cylinder */
				photonPtr->location = newPosition;
				photonPtr->travel_distance += distance;
			}
		}
		
		
		/* See if we fall outside of axial limits */
		if ((newPosition.z_position > colData[ColCurParams].colInBoundCyl.zMax) ||
				(newPosition.z_position < colData[ColCurParams].colInBoundCyl.zMin)) {
			
			action = colEnAc_Discard;
			break;
		}
		
		/* We know we are starting in the first layer */
		curLayer = 0;
		
		/* Determine which axial segment we are starting in */
		curSeg = colMCPETFindAxialSeg(newPosition.z_position, curLayer);
		
		
		/* Loop until photon escapes, or is absorbed */
		do {
		
			/* Set local current segment, this cuts down on the long data type
				references
			*/
			curSegInfo =
				&(ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[curSeg]);
				
			/* Get the material of this collimator segment */
			materialIndex = curSegInfo->MaterialIndex;
			
			/* Get the attenuation for the segment material */
			SubObjGetAttenuationInTomo(materialIndex, photonPtr->energy, &attenuation);
			
			/* Sample free paths to travel */
			PhgMathGetTotalFreePaths(&fpToGo);
			
			/* Compute the distance to travel, considering
			   no boundary intersections
			 */
			 distance = fpToGo/attenuation;
			
			/* Copy current layer */
			newLayer = curLayer;
			
			/* Project the photon through the current segment */
			action = colMCPETProject(photonPtr, curSegInfo, firstTime,
				&newPosition, &curSeg, &newLayer, &distance);
				
			/* Check results */
			{
				/* See if photon crossed a layer */
				if (action == colEnAc_LayerCross) {

					/* Check layer crossing for boundaries. If outer crossing
						is outer layer, it is detected. If inner crossing is
						minimum layer it is discarded. Otherwise the loop
						continues.
					*/
					if (newLayer == (LbFourByte)(ColRunTimeParams[ColCurParams].MCPETCol[0].NumLayers)) {
						action = colEnAc_Detect;
						break;
					}
					else if (newLayer < 0) {
						action = colEnAc_Discard;
						break;
					}
					else {
						/* Set new layer */
						curLayer = newLayer;
						
						/* Update our cylinders */
						colInitCylinders(curLayer);

						/* Determine which axial segment we are starting in */
						curSeg = colMCPETFindAxialSeg(newPosition.z_position, curLayer);
					}
				}
				
				/* See if photon crossed a segment */
				if (action == colEnAc_AxialCross) {
					/* Check segment crossing for boundaries. If crossing
						is boundary segment, then discard it.
					*/
					if (curSeg == (LbFourByte)(ColRunTimeParams[ColCurParams].MCPETCol[curLayer].NumSegments)) {
						action = colEnAc_Discard;
						break;
					}
					
					/* Check for minimum boundary */
					if (curSeg < 0) {
						action = colEnAc_Discard;
						break;
					}
				}
				
				/* If we made it to here, we want to update the 
					position and distance traveled.
				*/
				photonPtr->location = newPosition;
				photonPtr->travel_distance += distance;
				
				
				/* If photon should interact, pick a scatter angle and
					adjust the angle
				*/
				if (action == colEnAc_Interact) {
									
					/* Get probability of compton scatter */
					comptonToScatterProbability = SubObjGetProbComptToScatInTomo2(materialIndex,
						photonPtr->energy);
					
					/* Get probability of  scatter */
					scatterProbability = SubObjGetProbScatterInTomo2(materialIndex,
						photonPtr->energy);
					
					/* Sample for probability of interaction type */
					interactionProbability = PhgMathGetRandomNumber();
					
					/* See if absorption should occur */
					if (interactionProbability > scatterProbability) {
	
						/* Mark the photon for discarding (it's being absorbed) */
						action = colEnAc_Discard;
						
						/* Terminate the switch */
						break;
						
					} /* See if compton scatter should occur */
					else if (interactionProbability > (scatterProbability*comptonToScatterProbability)) {
						if (PHG_IsModelCoherentInTomo()) {
							EmisListDoCoherent(photonPtr, materialIndex);
						}
						else {
							continue;
						}
					}
					else {					
						/* Model Compton scatter (Using Klein-Nishina) */
						EmisListDoComptonInteraction(photonPtr);
					}
					
					/* Check energy for minimum threshold */
					if (photonPtr->energy < PhgRunTimeParams.PhgMinimumEnergy) {
						
						/* Account for the photon loss */
						PhoHStatIncTotLowEnPhotons();
						
						/* Mark the photon for discard */
						action = colEnAc_Discard;
						
						/* Terminate the tracking loop */
						break;
					}
					
					/* Increment the photon's scatter count */
					photonPtr->scatters_in_col += 1;

				}
				
				
			}
			
			/* Clear "first time" flag */
			firstTime = false;
				
		} while ((action != colEnAc_Detect) && (action != colEnAc_Discard) &&
			(action != colEnAc_Absorb));
		
		/* If we detect the photons, then update it's information */
		if (action == colEnAc_Detect) {
			/* Store the new location and distance traveled for the photon */
			{
				/* Save the new location */
				photonPtr->location = newPosition;
				photonPtr->travel_distance += distance;
			}
		}
		
	} while (false);
	
	#ifdef PHG_DEBUG
		/* Make sure we set some action */
		if (action == colEnAc_Null)
			PhgAbort("Failed to set an action in tracking through collimator (colMCPETTrack).", true);
	#endif
			
	return(action);
}

/*********************************************************************************
*
*			Name:			colMCPETProject
*
*			Summary:		Project a photon through the Monte Carlo PET collimator.
*							Truncate the projection based on the free paths to travel,
*							the inner/outer cylinder boundaries, the current
*							segment boundary.
*
*			Arguments:
*				PHG_TrackingPhoton				*photonPtr		- The photon.
*				Col_Monte_Carlo_PET_SegTy		*curSegInfo		- The current segment.
*				Boolean							firstTime		- Is this the first projection?
*				PHG_Position					*newPosPtr		- The new photon position.
*				LbFourByte						*curSegPtr		- The current segment index.
*				LbFourByte						*curLayerPtr	- The current layer index.
*				double							*distPtr		- The distance traveled.
*
*			Function return: The action to take based on the results of tracking.
*
*********************************************************************************/
ColEn_ActionTy colMCPETProject(PHG_TrackingPhoton *photonPtr,
					Col_Monte_Carlo_PET_SegTy *curSegInfo, Boolean firstTime, 
					PHG_Position *newPosPtr, LbFourByte *curSegPtr,
					LbFourByte *curLayerPtr, double *distPtr)
{
	double				distToInner = 0.0;			/* Distance to inner cylinder */
	double				distToOuter = 0.0;			/* Distance to inner cylinder */
	double				distToWall = 0.0;			/* Distance to donut wall */
	double				distToWall1 = 0.0;			/* Distance to donut wall */
	double				distToWall2 = 0.0;			/* Distance to donut wall */
	ColEn_ActionTy		action = colEnAc_Null;		/* What action will result? */
	
	 
	do { /* Process Loop */
	
		
		/*	If this is the first projection, we know the photon is
			traveling outwards so we don't worry about intersections
			with the inner cylinder.
			
			However, if not, if the photon is going to hit the inner cylinder, it will
			do it before hitting the outer cylinder. So we will start by
			getting the distance to the inner cylinder. If it is non-
			negative, and not "equal" to zero,
			then the photon is traveling towards the inner
			cylinder.
		*/
		if ((firstTime == false) && ((distToInner = colMCPETFindInnerDist((&photonPtr->location),
				&(photonPtr->angle))) > 0) && (
				PhgMathRealNumAreEqual(distToInner, 0.0, -7, 0, 0, 0) == false)){
		
			/* Find intersection with donut walls */
			distToWall1 = colMCPETFindWallIntersect(&(photonPtr->location),
				&(photonPtr->angle), colData[ColCurParams].colInBoundCyl.radius,
				colData[ColCurParams].colOutBoundCyl.radius,
				curSegInfo->InnerMaxZ,
				curSegInfo->OuterMaxZ);

			distToWall2 = colMCPETFindWallIntersect(&(photonPtr->location),
				&(photonPtr->angle), colData[ColCurParams].colInBoundCyl.radius,
				colData[ColCurParams].colOutBoundCyl.radius,
				curSegInfo->InnerMinZ,
				curSegInfo->OuterMinZ);
				
			/* Save distToWall as maximum positive (non-zero) of 1 & 2 */
			/* NOTE: if both are negative it is weeded out below */
			if ((distToWall1 > 0.0) || (distToWall2 > 0.0)){
				if (distToWall1 <= 0.0)
					distToWall = distToWall2;
				else if (distToWall2 <= 0.0)
					distToWall = distToWall1;
				else
					distToWall =  ((distToWall1 > distToWall2) 
						? distToWall1 : distToWall2);
			}
							
			/* Pick the smallest positive solution (we know distToInner is > 0) */
			if ((distToWall > 0) && (distToWall < distToInner)) {
				
				/* Compare distance to wall with given distance to travel,
					determines if photon will interact
				 */
				if (distToWall > *distPtr) {
				
					/* Photon will interact at distance given */
					action = colEnAc_Interact;
					break;
				
				}
				else {
				
					/* Photon will interact at cell boundary */
					*distPtr = distToWall;
					
					/* Inc/dec segment index depending on which wall was intersected */
					if (distToWall == distToWall1)
						*curSegPtr += 1;
					else
						*curSegPtr -= 1;
						
					/* Set our action */
					action = colEnAc_AxialCross;
					break;
				}
			}
			else {
				
				/* Compare distance to inner cylinder with given distance to travel,
					determines if photon will interact
				 */
				if (distToInner > *distPtr) {
				
					/* Photon will interact at given distance */
					action = colEnAc_Interact;
					break;
				
				}
				else {
					
					/* The photon will reach the inner cylinder surface
						without interacting. 
					*/
					*distPtr = distToInner;
					action = colEnAc_LayerCross;
					*curLayerPtr -= 1;
					break;
				}
			}
			
		}
		else {
		
			/* The photon is going to intersect with the outer cylinder */
			if (CylPosProjectToCylinder(&(photonPtr->location),
					&(photonPtr->angle), &(colData[ColCurParams].colOutBoundCyl), newPosPtr, &distToOuter)
					== false) {

				/* Distance is too large */
				distToOuter = MAXFLOAT;
			}
			
		
			/* Find intersection with donut wall */
			distToWall1 = colMCPETFindWallIntersect(&(photonPtr->location),
				&(photonPtr->angle), colData[ColCurParams].colInBoundCyl.radius,
				colData[ColCurParams].colOutBoundCyl.radius,
				curSegInfo->InnerMaxZ,
				curSegInfo->OuterMaxZ);
			
			distToWall2 = colMCPETFindWallIntersect(&(photonPtr->location),
				&(photonPtr->angle), colData[ColCurParams].colInBoundCyl.radius,
				colData[ColCurParams].colOutBoundCyl.radius,
				curSegInfo->InnerMinZ,
				curSegInfo->OuterMinZ);
			
			/* Save distToWall as minimum positive of 1 & 2 */
			/* NOTE: if both are negative it is weeded out below */
			if ((distToWall1 > 0.0) || (distToWall2 > 0.0)){
				if (distToWall1 <= 0.0)
					distToWall = distToWall2;
				else if (distToWall2 <= 0.0)
					distToWall = distToWall1;
				else
					distToWall =  ((distToWall1 > distToWall2) 
						? distToWall1 : distToWall2);
			}

			/* If distToWall is less than distToOuter, that is the maximum distance */
			if ((distToWall > 0) && (distToWall < distToOuter)) {
				
				/* Compare distance to wall with given distance to travel,
					determines if photon will interact
				 */
				if (distToWall > *distPtr) {
				
					/* Photon will interact at given distance */
					action = colEnAc_Interact;
					break;
				
				}
				else {
				
					/* Photon will interact at wall */
					*distPtr = distToWall;
					
					/* Inc/dec segment index depending on which wall was intersected */
					if (distToWall == distToWall1)
						*curSegPtr += 1;
					else
						*curSegPtr -= 1;
						
					/* Set our action */
					action = colEnAc_AxialCross;
					break;
				}
			}
			else {
				
				/* Compare distance to inner cylinder with given distance to travel,
					determines if photon will interact
				 */
				if (distToOuter > *distPtr) {
				
					/* Photon will interact at given distance */
					action = colEnAc_Interact;
					break;
				
				}
				else {
					
					/* The photon will reach the outer cylinder surface
						without interacting. 
					*/
					*distPtr = distToOuter;
					action = colEnAc_LayerCross;
					*curLayerPtr += 1;
					
					break;
				}
			}
		}
				
	} while (false);
	
	#ifdef PHG_DEBUG
		/* Make sure we set some action */
		if (action == colEnAc_Null)
			PhgAbort("Failed to set an action in tracking through collimator (colMCPETProject).", true);
	#endif
	
	/* Project to new location */
	newPosPtr->x_position = photonPtr->location.x_position +
		(photonPtr->angle.cosine_x * *distPtr);
			
	newPosPtr->y_position = 	photonPtr->location.y_position +
		(photonPtr->angle.cosine_y * *distPtr);
			
	newPosPtr->z_position = 	photonPtr->location.z_position +
		(photonPtr->angle.cosine_z * *distPtr);
			
	return(action);
}

/*********************************************************************************
*
*			Name:			colMCPETFindInnerDist
*
*			Summary:		Determine the distance to the inner cylinder.
*
*			Arguments:
*				PHG_Position	*posPtr	- The position of the photon.
*				PHG_Direction	*dirPtr	- The direction of the photon.
*
*			Function return: Distance to inner cylinder, negative if
*								it doesn't intersect.
*
*********************************************************************************/
double colMCPETFindInnerDist(PHG_Position *posPtr,
				PHG_Direction *dirPtr)
{
	double		quadA = 0.0;				/* "a" in the quadratic formula */
	double		quadB = 0.0;				/* "b" in the quadratic formula */
	double		quadC = 0.0;				/* "c" in the quadratic formula */
	double		xCord = 0.0;				/* X coordinate normalized to center of cylinder */
	LbUsFourByte	numRoots;				/* number of roots returned from quadratic formula */	
	double		yCord = 0.0;				/* Y coordinate normalized to center of cylinder */
	double		dist1 = 0.0;				/* Temp distance to inner cylinder */
	double		dist2 = 0.0;				/* Temp distance to inner cylinder */
	double		distToInner = -1.0;			/* Distance to inner cylinder */
	
	/* Compute shortest positive distance to inner cylinder.
		If the intersection occurs,
		which it may not, intersectInner gets set to true
		and innerDist is initialized
	 */
	 do { /* Process Loop */
	 
		/* Compute normalized x/y coordinates */
		xCord = posPtr->x_position - colData[ColCurParams].colInBoundCyl.centerX;
		yCord = posPtr->y_position - colData[ColCurParams].colInBoundCyl.centerY;
		
		/* Compute a */
		quadA = 1.0 - PHGMATH_Square(dirPtr->cosine_z);
		
		/* Compute b */
		quadB = 2.0 * ( (dirPtr->cosine_x * xCord) +
			(dirPtr->cosine_y * yCord));
		
		/* Compute c */
		quadC = 2.0 * PHGMATH_Square(xCord) + PHGMATH_Square(yCord) -
				PHGMATH_Square(colData[ColCurParams].colInBoundCyl.radius);
		
		numRoots = PhgMathSolveQuadratic(quadA, quadB, quadC, &dist1, &dist2);
		
		/* if numRoots is 0 we don't intersect the inner cylinder */
		if (numRoots == 0)
			break;

		/* We want the shortest positive distance.
			We'll start by checking the case of
			one of the distances being essentially
			zero. This is necessary because a distance
			might be negative, but so small that we want to
			consider it to be zero.
			
			Also note, that if it is considered zero
			then there can not be any smaller positive
			value.
		*/
		
		if (numRoots == 1) {
			
			if (PhgMathRealNumAreEqual(dist1, 0.0, -5, 0, 0, 0) == true) {
			
				/* We will set distance to be zero */
				distToInner = 0.0;
				
			} else if ( dist1 > 0 ) {
			
				distToInner = dist1;
				
			} else {
			
				/* leave distToInner = -1--there is no positive intersection */
				break;
				
			}
			
		/* otherwise there are two roots */
		} else if (PhgMathRealNumAreEqual(dist1, 0.0, -5, 0, 0, 0) == true
				|| PhgMathRealNumAreEqual(dist2, 0.0, -5, 0, 0, 0) == true ) {
			
			/* We will set distance to be zero */
			distToInner = 0.0;
			
		} else if (dist1 > 0) {
		
			distToInner = dist1;
			
		} else if (dist2 > 0) {
		
			distToInner = dist2;
		
		} else {
		
			break;
			
		}
	} while (false);			

	return(distToInner);
}

/*********************************************************************************
*
*			Name:			colMCPETFindWallIntersect
*
*			Summary:		Determine the distance to the intersection of the
*							collimator axial wall.
*
*			Arguments:
*				PHG_Position	*posPtr	- The position of the photon.
*				PHG_Direction	*dirPtr	- The direction of the photon.
*				double			innerR	- The inner adius of the collimator.
*				double			outerR	- The outer radius of the collimator.
*				double			innerZ	- See formula.
*				double			outerZ	- See formula.
*
*			Function return: Distance to intersection (negative if it doesn't).
*
*********************************************************************************/
double colMCPETFindWallIntersect(PHG_Position *posPtr, PHG_Direction *dirPtr,
			double innerR, double outerR, double innerZ, double outerZ)
{
	double 	distToWall = -1.0;	/* Distance to wall intersection */
	double	s;					/* Substitution value for simplification */
	double	d;					/* Substitution value for simplification */
	double	a, b, c;			/* Parameters for quadratic formula */
	double	bSq_4ac;			/* Value under radical of quadratic formula */
	double	Sqrt_bSq_4ac;		/* Square root of value under radical of quadratic formula */
	double	t1, t2;				/* Two solutions to quadratic formula */
	
	do { /* Process Loop */
	
		/* NOTE: There could be an optimization here for rectanguler (non tapered)
			collimators, but not enough information is provided. Depending on the
			cosine_z parameter, we need the inner Min/Max and outer Min/Max Z values.
			
			For now we'll just use the general case.
		*/
		/* Check z values for rectangular collimator */
		if (innerZ == outerZ) {
			
			/* Compute the distance */
			distToWall = (outerZ - posPtr->z_position)/dirPtr->cosine_z;

			/*	Break out of process loop.
				 NOTE: if cosine_z == 0, no intersection distToWall remains negative
			*/	
			break;
		}
		
		/* IF we made it to here then we use the general formula for tapered
			collimators.
		*/
		
		/* Compute S */
		s = (outerR-innerR)/(outerZ-innerZ);
		
		/* Compute D */
		d = (s * posPtr->z_position) - (s * innerZ) + innerR;
		
		/* Compute a */
		a = PHGMATH_Square(dirPtr->cosine_x) + PHGMATH_Square(dirPtr->cosine_y)
			- (PHGMATH_Square(s) * PHGMATH_Square(dirPtr->cosine_z));
		
		/* Compute b */
		b = (2 * dirPtr->cosine_x * posPtr->x_position) +
			(2 * dirPtr->cosine_y * posPtr->y_position) -
			(2 * s * dirPtr->cosine_z * d);	
		
		/* Compute c */
		c = PHGMATH_Square(posPtr->x_position) +
			PHGMATH_Square(posPtr->y_position) -
			PHGMATH_Square(d);
		
		/* Compute value under the radical in quadratic formula */	
		bSq_4ac = PHGMATH_Square(b) - (4 * a * c);
		
		/* If value under radical is negative, there is no solution */
		if (bSq_4ac < 0) {
			break;
		}
		
		/* Take the square root */
		Sqrt_bSq_4ac = PHGMATH_SquareRoot(bSq_4ac);
		
		/* Compute two solutions to quadratic */
		t1 = (-b + Sqrt_bSq_4ac)/(2*a);
		t2 = (-b - Sqrt_bSq_4ac)/(2*a);
		
		/* If t1 or t2 is essentially zero, we don't want it */
		if (PhgMathRealNumAreEqual(t1, 0.0, -7, 0, 0, 0)
				== true) {
			
			t1 = -1;
			
		}
		else  if (PhgMathRealNumAreEqual(t2, 0.0, -7, 0, 0, 0)
				== true) {
			
			t2 = -1;
			
		} 
		
		if ((t1 > 0) && (t2 > 0)) {
			distToWall = ((t1 < t2) ? t1 : t2);
		}
		else  if ((t1 > 0) || (t2 > 0)){
			distToWall = ((t1 > 0) ? t1 : t2);
		}

	} while (false);
	
	return(distToWall);
}

/*********************************************************************************
*
*			Name:			colMCPETFindAxialSeg
*
*			Summary:		Figure out which axial segment the photon is in.
*							
*			Arguments:
*				double			zPos		- The z position of the photon.
*				LbUsFourByte	curLayer	- The current layer.
*
*			Function return: None.
*
*********************************************************************************/
LbFourByte colMCPETFindAxialSeg(double zPos, LbFourByte curLayer)
{
	LbUsFourByte	segIndex;		/* The photon's segment */
	
	/* Loop  segments to find the one that contains the photon */
	for (segIndex = 0; segIndex < ColRunTimeParams[ColCurParams].MCPETCol[curLayer].NumSegments; segIndex++) {
	
		if ((ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].InnerMinZ <= zPos)
				&&
				(ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex].InnerMaxZ >= zPos)) {
			break;
		}
	}
	
	#ifdef PHG_DEBUG
		/* Verify we found a segment */
		if (segIndex == ColRunTimeParams[ColCurParams].MCPETCol[curLayer].NumSegments) {
			sprintf(colErrStr, "Photon starting out on collimator layer out of z-bounds (colMCPETFindAxialSeg)"
				"\nzPos = %3.5lf, zMin = %3.5lf, zMax = %3.5lf\n",
				zPos, ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[0].InnerMinZ,
				ColRunTimeParams[ColCurParams].MCPETCol[curLayer].Segments[segIndex-1].InnerMaxZ);
				
			PhgAbort(colErrStr, true);
		}
		
	#endif
	
	
	return(segIndex);
}

/*********************************************************************************
*
*			Name:			ColSPECTPhotons
*
*			Summary:		Update the binning images with the current batch of
*							detected photons.
*
*			Arguments:
*				PHG_Decay			decayPtr			- The decayPtr that started the process.
*				PHG_TrackingPhoton *photons			- The photons detected.
*				LbUsFourByte 		numPhotons		- The number of blue photons.
*				CollimatedPhotonsTy *colPhotonsPtr	- The accepted collimated photons
*			Function return: None.
*
*********************************************************************************/
void ColSPECTPhotons(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
		CollimatedPhotonsTy *colPhotonsPtr)
{
			
	/* Clear the counters, note you must clear the PINKs also */			
	colPhotonsPtr->NumCollimatedBluePhotons = 0;
	colPhotonsPtr->NumCollimatedPinkPhotons = 0;
	colData[ColCurParams].colDetectedPinkPhotonIndex = 0;
	colData[ColCurParams].colDetectedBluePhotonIndex = 0;

	/* Compute statistics */
	colData[ColCurParams].colTotBluePhotons += numPhotons;

	/* Call appropriate collimation routine based on user specified type */
	switch(ColRunTimeParams[ColCurParams].ColType) {
	
		case ColEn_simple_spect:
		
			/* Call simple pet collimator */
			colSimpleSPECT(decayPtr, photons, numPhotons,
				colPhotonsPtr);
				
			/* You are done */
			break;
			
		case ColEn_monte_carlo_pet:
			colMonteCarloSPECT(decayPtr,
				photons,  numPhotons,
				colPhotonsPtr);
			break;
		
		case ColEn_unc_spect:
			UNCCollimate(decayPtr, photons, numPhotons, colPhotonsPtr);
			break;
		
		case ColEn_slat:
			ColSlatSPECT(decayPtr, photons, numPhotons, colPhotonsPtr);
			break;
			
		default:
			PhgAbort("Invalid collimator type in ColRunTimeParams[ColCurParams] (ColSPECTPhotons)",false);
			break;
			
	}
	
	/* Update statistics */
	colData[ColCurParams].colTotAccBluePhotons += colPhotonsPtr->NumCollimatedBluePhotons; 
	
	/* Write collimated photons to history file, if requested */
	if (COL_IsDoHistory()) {

		/* Write the photons */
		if (PhoHFileWriteDetections(&(colData[ColCurParams].colHistFileHk), decayPtr,
				colData[ColCurParams].colDetectedBluePhotons,
				colData[ColCurParams].colDetectedBluePhotonIndex,
				colData[ColCurParams].colDetectedPinkPhotons,
				colData[ColCurParams].colDetectedPinkPhotonIndex) == false) {
			
			/* Abort Program execution */
			PhgAbort("Got failure from PhoHFileWriteDetections trapped in ColSPECTPhotons.",false);
		}
	}

}

/*********************************************************************************
*
*			Name:			colSimpleSPECT
*
*			Summary:		Perform "simple" collimation. This routine simply
*							projects the photon the distance of the collimator
*							depth. If the photon is outside of the axial limits
*							it is rejected.
*
*			Arguments:
*				PHG_Decay			decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton *photons			- The blue photons detected.
*				LbUsFourByte 		numPhotons		- The number of blue photons.
*				CollimatedPhotonsTy *colPhotonsPtr	- The accepted collimated photons
*			Function return: None.
*
*********************************************************************************/
void colSimpleSPECT(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
		CollimatedPhotonsTy *colPhotonsPtr)
{
	LbUsFourByte		photonIndex;				/* The current photon */

	/* Bolt if there is nothing to do */
	if (numPhotons == 0) {
		return;
	}
	
	do { /* Process Loop */
	
		/* Clear the counters */
		colPhotonsPtr->NumCollimatedBluePhotons = 0;
		colPhotonsPtr->NumCollimatedPinkPhotons = 0;
		colData[ColCurParams].colDetectedPinkPhotonIndex = 0;
		colData[ColCurParams].colDetectedBluePhotonIndex = 0;

		/* Collimate the photons */
		for (photonIndex = 0; photonIndex < numPhotons; photonIndex++){
			
			/* Let user modify and/or reject photons */
			if (ColUsrModSPECTPhotonsFPtr && 
					(*ColUsrModSPECTPhotonsFPtr)(&ColRunTimeParams[ColCurParams], decayPtr,
						&(photons[photonIndex])) == false) {
				
				/* They rejected it so go to next pink */
				continue;
			}
			

			/* We made it to here so save the collimated photons */
			{
				colPhotonsPtr->CollimatedTrkngBluePhotons[colPhotonsPtr->NumCollimatedBluePhotons]
					= photons[photonIndex];

				/* Increment the counters */
				colPhotonsPtr->NumCollimatedBluePhotons++;

				/* Update our detected photon block if doing history file */
				if  (COL_IsDoHistory()) {
					ColUpdateCollimatedPhotonBlock(&photons[photonIndex]);
				}
			}
		}

	} while (false);
	
}

/*********************************************************************************
*
*			Name:			colMonteCarloSPECT
*
*			Summary:		Perform "PET" collimation for singles studies.
*
*			Arguments:
*				PHG_Decay			decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton *photons			- The photons detected.
*				LbUsFourByte 		numPhotons		- The number of photons.
*				CollimatedPhotonsTy *colPhotonsPtr	- The accepted collimated photons
*			Function return: None.
*
*********************************************************************************/
void colMonteCarloSPECT(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
		CollimatedPhotonsTy *colPhotonsPtr)
{
	LbUsFourByte		photonIndex;		/* The current photon */
	ColEn_ActionTy		action;				/* What action to take during tracking */
	
	do { /* Process Loop */
	
		/* Clear the counters */
		colPhotonsPtr->NumCollimatedBluePhotons = 0;
		colPhotonsPtr->NumCollimatedPinkPhotons = 0;
		colData[ColCurParams].colDetectedPinkPhotonIndex = 0;
		colData[ColCurParams].colDetectedBluePhotonIndex = 0;
		
		/* We use a loop to process the blues and then the pinks. This
			minimizes code duplication.
		*/
		/* Collimate the photons */
		for (photonIndex = 0; photonIndex < numPhotons; photonIndex++){
			
			/* Let user modify and/or reject photons */
			if (ColUsrModSPECTPhotonsFPtr && 
					(*ColUsrModSPECTPhotonsFPtr)(&ColRunTimeParams[ColCurParams], decayPtr,
						&(photons[photonIndex])) == false) {
				
				/* They rejected it so go to next photon */
				continue;
			}
			
			/* Track the photon */
			action = colMCPETTrack(&(photons[photonIndex]));
							
			/* Store the photon if the action is to detect */
			if (action == colEnAc_Detect) {	 

				/* We made it to here so save the collimated photons */
					colPhotonsPtr->CollimatedTrkngBluePhotons[colPhotonsPtr->NumCollimatedBluePhotons]
						= photons[photonIndex];
	
					/* Increment the counters */
					colPhotonsPtr->NumCollimatedBluePhotons++;
	
					/* Update our detected photon block if doing history file */
					if  (COL_IsDoHistory()) {
						ColUpdateCollimatedPhotonBlock(&photons[photonIndex]);
					}
			}
		}
	} while (false);
	
}

/*********************************************************************************
*
*			Name:			ColTerminate
*
*			Summary:		Process the bin buffers and terminate the module.
*
*			Arguments:
*			Function return: None.
*
*********************************************************************************/
void ColTerminate()
{
	Boolean			okay = false;	/* Process Flag */
	LbUsFourByte	layers;			/* Loop Control */
	LbUsFourByte	index;			/* Loop Control */
	
	
	do { /* Process Loop */
		
		/* Call the user termination routine */
		if (ColUsrTerminateFPtr) {
			(*ColUsrTerminateFPtr)(&ColRunTimeParams[ColCurParams]);
		}
		
		/* Call UNC termination if it is being used */
		if (ColRunTimeParams[ColCurParams].ColType == ColEn_unc_spect)
			UNCColTerminate();
				
		/* Close the history file if we created one */
		if (COL_IsDoHistory()) {
			if (!PhoHFileClose(&(colData[ColCurParams].colHistFileHk))) {
				break;
			}
		}
		
		/* Free up allocated memory */
		if (ColRunTimeParams[ColCurParams].MCPETCol != 0) {
			for (layers = 0; layers < ColRunTimeParams[ColCurParams].MCPETCol[0].NumLayers; layers++) {
				for (index = 0; index < ColRunTimeParams[ColCurParams].MCPETCol[0].NumSegments;index++) {
					if (ColRunTimeParams[ColCurParams].MCPETCol[layers].Segments != 0)
						LbMmFree((void **) &(ColRunTimeParams[ColCurParams].MCPETCol[layers].Segments));
				}
			}
			
			LbMmFree((void **) &ColRunTimeParams[ColCurParams].MCPETCol);
		}
		okay = true;
	} while (false);
	
	/* Do error handling here */
	if (!okay) {
		ErHandle("Failed to process binning files.", false);
	}
}

/*********************************************************************************
*
*			Name:		ColUpdateDetectedPhotonBlock
*
*			Summary:	Make a new detected photon (To write to the history file).
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*				
*			Function return: None.
*
*********************************************************************************/
void ColUpdateCollimatedPhotonBlock(PHG_TrackingPhoton	*trackingPhotonPtr)	
{
	
	/* Store deteced photon in bin */
	if (PHG_IsBlue(trackingPhotonPtr)) {
	
		colData[ColCurParams].colDetectedBluePhotons[colData[ColCurParams].colDetectedBluePhotonIndex] = *trackingPhotonPtr;
		colData[ColCurParams].colDetectedBluePhotonIndex++;
	}
	else {
	
		colData[ColCurParams].colDetectedPinkPhotons[colData[ColCurParams].colDetectedPinkPhotonIndex] = *trackingPhotonPtr;
		colData[ColCurParams].colDetectedPinkPhotonIndex++;
	}

}

#undef COLLIMATOR
