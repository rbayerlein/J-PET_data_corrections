/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1992-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		SubObj.c
*			Revision Number:	1.10
*			Date last revised:	23 July 2012
*			Programmer:			Steven Vannoy
*			Date Originated:	20 August 1992
*
*			Module Overview:	Sub-Object Access routines.
*
*			References:			'Sub-Object Access Processes' PHG design.
*
**********************************************************************************
*
*			Global functions defined:
*				SubObjAdjPosDecayLocation
*				SubObjCalcTimeBinDecays
*				SubObjCreate
*				SubObjGenVoxAngCellDecay
*				SubObjGetAttCellIndexes
*				SubObjGetInnerCellDistance
*				SubObjGetObjCylinder
*				SubObjGetProbComptToScatter
*				SubObjGetStartingProdValues
*				SubObjDumpObjects
*
*			Global variables defined:		none
*
*			Macros defined:
*				SUBOBJGetSliceMinX
*				SUBOBJGetSliceMinY
*				SUBOBJGetSliceMinZ
*				SUBOBJGetSliceManX
*				SUBOBJGetSliceManY
*				SUBOBJGetSliceManZ
*				SUBOBJGetSliceMaxXIndex
*				SUBOBJGetSliceMaxYIndex
*				SUBOBJGetSliceVoxelWidth
*				SUBOBJGetSliceVoxelHeight
*				SUBOBJGetSliceVoxelDepth
*				SUBOBJGetSliceWidth
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
*			Revision date:		19 December 2012
*
*			Revision description:	Converted timing actions to use LbTiming
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		10 October 2012
*
*			Revision description:	Reduced timing options by system to just two
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		R Harrison
*
*			Revision date:		Oct 2009
*
*			Revision description:  Fixed SubObjGenVoxAngCellDecay to set
*							decayType.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):	R Harrison	
*
*			Revision date:	April 28, 2008	
*
*			Revision description:  Changed call to SubObjGtPositionIndexes at
*				end of SubObjGenVoxAngCellDecay to occur only if a new decay
*				was created.  (Previous failure to do so was a bug.)
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		R Harrison
*
*			Revision date:		2007
*
*			Revision description:	Support for:
*				- eight byte number of decays and real scan time.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		R Harrison
*
*			Revision date:		August 2006
*
*			Revision description:	Moved SubObjNumTissues to SubObj.h
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		22 December 2005
*
*			Revision description:	Added parameterized GetProb functions.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		January 2005
*
*			Revision description:	Added decay time to decay record.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		21 November 2003
*
*			Revision description:	Added changes for Mac CW 9 compatibility:
*										Changed PATH_SEPARATOR for Mac to '/'
*
*********************************************************************************/

#define	SUB_OBJECT

#include <stdio.h>
#include <string.h>
#include <ctype.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbFile.h"
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
#include "PhoHStat.h"
#include "ProdTbl.h"
#include "SubObj.h"
#include "ColUsr.h"
#include "phg.h"
#include "PhgBin.h"

/* LOCAL CONSTANTS */
#define SUBOBJ_ACC_STRAT_BINS			48								/* Count of stratification bins within acceptance angle */
#define SUBOBJ_OBJECT_FILEPREFIX		"phgobj."						/* Prefix to object description file */	
#define SUBOBJ_ATTEN_FILE_NAME			"phg.attenuation"				/* Prefix to attenuation description file */	
#define SUBOBJ_ACT_FILE_NAME			"phg.activity"					/* Prefix to activity description file */
#define SUBOBJ_CONST_ACT				-1								/* Activity is constant */

#define SUBOBJ_VERYSMALLATTENUATION		0.0000000001
#define SUBOBJ_NUM_ENERGY_BINS			1000								/* Number of energy bins */
#define SUBOBJ_NUM_TISSUE_TRANSLATIONS	2301								/* Entries in the translation tables */
#define SUBOBJ_NUM_COH_THETAS			180								/* Number of angles in coherent scatter angle table */

#define SUBOBJ_NUM_1KEV_ENERGIES	150
#define SUBOBJ_NUM_10KEV_ENERGIES	15
#define SUBOBJ_NUM_100KEV_ENERGIES	7
#define SUBOBJ_MAX_10KEV_ENERGIES	(SUBOBJ_NUM_1KEV_ENERGIES + (SUBOBJ_NUM_10KEV_ENERGIES * 10))
#define SUBOBJ_MAX_100KEV_ENERGIES	(SUBOBJ_NUM_1KEV_ENERGIES + SUBOBJ_MAX_10KEV_ENERGIES + (SUBOBJ_NUM_100KEV_ENERGIES * 100))
#define	SUBOBJ_NUM_COH_ENERGIES		172	/* Energies total */
#define SUBOBJ_NUM_COH_ANGLES			100	/* Cosine theta split into 100 bins */

#define SUBOBJ_NUM_ANGLES		100
#define MAX_NUM_MATERIALS		2301

#ifdef __MWERKS__
/* Changed to Unix format for easier portability of parameter files */
/*#define	PATH_SEPARATOR			':'*/
#define	PATH_SEPARATOR			'/'
#else
#define PATH_SEPARATOR			'/'
#endif

/* LOCAL TYPES */	
typedef LbUsEightByte					subObjDecaySliceTy;

typedef double							subObjDecayWeightSliceTy;		/* Decay weight per angle */

typedef struct {
	double	D;		/* Density */
	double	A;		/* Effective Atomic Weight */
	double	Z;		/* Effective Atomic Number */
} subObjMaterialDAZTy;

typedef struct {
	double	energy;														/* Energy level */
	double	attenuation;												/* Attenuation coefficient */
	double	probComptonToScatter;										/* Probability of compton/scatter */
	double	probScatter;												/* Probability of compton + coherent */
} subObjAttenTy;

typedef subObjAttenTy					subObjTissueAttenTy[SUBOBJ_NUM_ENERGY_BINS];		/* Attenuation values for a given tissue */
typedef subObjTissueAttenTy				*subObjTissueAttenTblTy;		/* Attenuation table */

#ifdef ABANDONED_WAY
typedef struct {
	double	theta;		/* The scatter angle */
	double	cumProb;	/* Cumulative probability of that angle */
} subObjCohScatAngleTy;
#endif


typedef	struct {

	LbUsFourByte	energy;
	LbUsFourByte	materialIndex;
	double			angleProbabilities[SUBOBJ_NUM_COH_ANGLES];
	} subObjCoScatAngleTy;

typedef subObjCoScatAngleTy subObjCoScatAngleTblTy[SUBOBJ_NUM_COH_ENERGIES];
typedef char				subObjNameTy[80];
typedef char				subObjLabelTy[1024];

/* LOCAL GLOBALS */
static	Boolean						subObjIsInitialized = false;		/* Initialization flag */
static	char						subObjErrStr[1024];
static	double						subObjSimulatedDivRealDetected;		/* Ratio used to determine number of decays per voxel */
static	LbUsFourByte				subObjCurTimeBin;					/* Our current time bin */
static	LbUsFourByte				SubObjNumActIndexes;				/* Number of activity indexes */
static	double						SubObjCurTimeBinDuration;			/* Duration of current time bin */
static	LbUsFourByte				SubObjAngleRoundUpCount;			/* Count of number of times an angle was rounded up */
static	subObjDecaySliceTy			*SubObjDecaySlice = 0;				/* Number of decays for the slice */
static	subObjDecayWeightSliceTy	*SubObjDecayWeightSlice = 0;			/* Weight of decays for the slice */
static	SubObjTissueTableTy			SubObjTissueTable;					/* Tissue activity table */
static	subObjTissueAttenTblTy		SubObjTissueAttenTableNoCoh = 0;	/* Tissue attenuation table */
static	subObjTissueAttenTblTy		SubObjTissueAttenTableCoh = 0;		/* Tissue attenuation table */
static	subObjNameTy				*subObjMaterialNames;				/* Names of attenuation materials */
static 	subObjMaterialDAZTy			*subObjMaterialDAZ;					/* Density, weight, and number of material */
static	subObjCoScatAngleTblTy		*subObjCohScatAngles;
static	LbUsFourByte				subObjNumCohMaterials;
#ifdef old_way
static subObjLabelTy	convertCohMaterialLabels[] = {
				"air",
				"water",
				"blood",
				"bone",
				"brain",
				"heart",
				"lung",
				"muscle",
				"lead",
				"NaI",
				"BGO",
				"iron",
				"graphite",
				"tin",
				"GI_tract",
				"con_tissue",
				"copper",
				"perfect_absorber",
				"LSO",
				"GSO",
				"aluminum",
				"tungsten"
				};
#endif
static subObjLabelTy	convertCohMaterialLabels[MAX_NUM_MATERIALS];

static	LbTmTimingType		SubObjTrackPhotonsStartTime;		/* Time we started tracking photons */

#ifdef PHG_DEBUG
static	Boolean	subObjDoWriteRandSeed = false;
static	Boolean	subObjDoReadRandSeed = false;
#endif
/* LOCAL MACROS */

/*********************************************************************************
*
*			Name:			SubObjGetNumActVoxels
*
*			Summary:		Return the number of activity voxels in a given slice.
*
*			Arguments:
*				LbUsFourByte	sliceIndex	- Current slice.
*
*			Function return: Number of voxels.
*
*********************************************************************************/
#define SubObjGetNumActVoxels(sliceIndex) \
	(SubObjObject[(sliceIndex)].actNumYBins * SubObjObject[(sliceIndex)].actNumXBins)

/*********************************************************************************
*
*			Name:			SubObjGetTissueActivity
*
*			Summary:		Return the activity for a given tissue/time.
*
*			Arguments:
*				LbUsFourByte	tissueIndex	- Which tissue.
*				LbUsFourByte	simTime		- Time within simulation.
*
*			Function return: Activity level.
*
*********************************************************************************/
#define SubObjGetTissueActivity(tissueIndex, simTime) \
	(SubObjTissueTable.tissueValues[(tissueIndex)].activityValues[(simTime)].activityLevel)

/*********************************************************************************
*
*			Name:			SubObjGetTissueIndex
*
*			Summary:		Return the tissue index for given x/y voxel.
*
*			Arguments:
*				LbUsFourByte	sliceIndex	- Current slice.
*				LbUsFourByte	yIndex		- Which y value.
*				LbUsFourByte	xIndex		- Which x value.
*
*			Function return: Tissue index.
*
*********************************************************************************/
#define SubObjGetTissueIndex(sliceIndex, yIndex, xIndex) \
	(SubObjObject[(sliceIndex)].activityArray[(yIndex) * SubObjObject[(sliceIndex)].actNumXBins + (xIndex)])


/*********************************************************************************
*
*			Name:		SUBOBJGetSliceActivity
*
*			Summary:	Return activity for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: LbUsFourByte, the slice's activity.
*
*********************************************************************************/
#define SUBOBJGetSliceActivity(sliceIndex) (SubObjObject[sliceIndex].sliceActivity)

/* PROTOTYPES */
double			subObjCellGetNumReal(LbUsFourByte sliceIndex, LbUsFourByte angleIndex,
					double sliceActivity);	
Boolean			subObjCreateObject(void);
void			subObjCalcSliceActivity(LbUsFourByte sliceIndex);
Boolean			subObjGetNextVoxAngCellDecay(void);
double			subObjVoxelCellGetNumReal(LbUsFourByte sliceIndex, LbUsFourByte angleIndex,
					LbUsFourByte voxelIndex, LbUsFourByte currentTime);	
void			subObjCalcSliceTimeBinDecays(LbUsFourByte currentTimeBin,
					LbUsFourByte sliceIndex);
Boolean			subObjAllocDecaySlice(LbUsFourByte sliceIndex);
Boolean			subObjInitCoherentTbl(void);
Boolean			subObjInitCoherentBinaryTbl(void);
void			subObjPrintStatus(void);
void			subObjComputePosRangeConstants(	double atomicNumber, double atomicWeight,
					double *b1Ptr, double *b2Ptr );

/*********************************************************************************
*
*			Name:		SubObjGetObjGeometry
*
*			Summary:	Get the subobject geometry from the parameter file.
*						NOTE THIS ROUTINE IS RECURSIVE.
*
*			Arguments:
*				LbPfHkTy		paramFlHk	- The parameter file hook.
*				LbUsFourByte	numParams	- The number of parameters.
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean	SubObjGetObjGeometry(LbPfHkTy paramFlHk, LbUsFourByte numParams)
{
	double					paramBuffer[LBPF_PARAM_LEN];	/* Parameter buffer */
	Boolean					okay = false;					/* Process flag */
	Boolean					isEOF;							/* End of file flag */
	Boolean					switchOkay;						/* Switch flag */
	LbPfEnPfTy				paramType;						/* Parameter type */
	LbUsFourByte			paramSize;						/* Parameter size */
	char					paramLabel[LBPF_LABEL_LEN];		/* Parameter label space */
	char					errString[256];					/* Storage for error strings */
	PhgEn_RunTimeParamsTy	whichParam;						/* Also current parameter */
	LbUsFourByte			paramIndex;						/* Current parameter */
	
	do { /* Process Loop */
	
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
			whichParam = PhgLookupRunTimeParamLabel(paramLabel);
			switchOkay = true;
			switch (whichParam) {
				
				case SubObjEn_num_slices:
						SubObjNumSlices = *((LbUsFourByte *) paramBuffer);

						/* Allocate the slice info */
						if ((SubObjObject = (SubObjSliceInfoArrayTy)
								LbMmAlloc(sizeof(SubObjSliceInfoTy) * SubObjNumSlices)) == 0) {
							
							switchOkay = false;
						}
						SubObjCurSliceIndex = 0;
					break;
					
				case SubObjEn_slice:
				
					/* Call this routine to process the slice */
					switchOkay = SubObjGetObjGeometry(paramFlHk,
							*((LbUsFourByte *) paramBuffer));
						
					/* Increment to next slice */
					SubObjCurSliceIndex++;
					break;
					
				case SubObjEn_slice_number:
					SubObjObject[SubObjCurSliceIndex].sliceNum = 
						*((LbUsFourByte *) paramBuffer);
					break;

				case SubObjEn_zMin:
					SubObjObject[SubObjCurSliceIndex].zMin = 
						*((double *) paramBuffer);
					break;

				case SubObjEn_zMax:
					SubObjObject[SubObjCurSliceIndex].zMax = 
						*((double *) paramBuffer);
					break;
					
				case SubObjEn_xMin:
					SubObjObject[SubObjCurSliceIndex].xMin = 
						*((double *) paramBuffer);
					break;
					
				case SubObjEn_xMax:
					SubObjObject[SubObjCurSliceIndex].xMax = 
						*((double *) paramBuffer);
					break;

				case SubObjEn_yMin:
					SubObjObject[SubObjCurSliceIndex].yMin = 
						*((double *) paramBuffer);
					break;

				case SubObjEn_yMax:
					SubObjObject[SubObjCurSliceIndex].yMax = 
						*((double *) paramBuffer);
					break;

				case SubObjEn_num_act_X_bins:
					SubObjObject[SubObjCurSliceIndex].actNumXBins = 
						*((LbUsFourByte *) paramBuffer);
					break;

				case SubObjEn_num_att_X_bins:
					SubObjObject[SubObjCurSliceIndex].attNumXBins = 
						*((LbUsFourByte *) paramBuffer);
					break;

				case SubObjEn_num_act_Y_bins:
					SubObjObject[SubObjCurSliceIndex].actNumYBins = 
						*((LbUsFourByte *) paramBuffer);
					break;

				case SubObjEn_num_att_Y_bins:
					SubObjObject[SubObjCurSliceIndex].attNumYBins = 
						*((LbUsFourByte *) paramBuffer);
					break;
				
				/*	Note, if this is an old style parameter file then
					both counts are the same
				*/
				case SubObjEn_num_X_bins:
					SubObjObject[SubObjCurSliceIndex].actNumXBins = 
						*((LbUsFourByte *) paramBuffer);
					SubObjObject[SubObjCurSliceIndex].attNumXBins = 
						*((LbUsFourByte *) paramBuffer);
					break;

				
				/*	Note, if this is an old style parameter file then
					both counts are the same
				*/
				case SubObjEn_num_Y_bins:
					SubObjObject[SubObjCurSliceIndex].actNumYBins = 
						*((LbUsFourByte *) paramBuffer);
					SubObjObject[SubObjCurSliceIndex].attNumYBins = 
						*((LbUsFourByte *) paramBuffer);
					break;
				default:
					sprintf(errString, "(SubObjGetObjGeometry) Unknown (hence unused) parameter (%s).\n",
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
*			Name:		SubObjGetProbScatterInObj2
*
*			Summary:	Get probability of scatter given material
*						and photon energy.
*			Arguments:
*			LbUsFourByte	index	- The tissue index.
*			double			energy	- The tissue energy.
*
*			Function return: Probability of compton scatter.
*
*********************************************************************************/
double	SubObjGetProbScatterInObj2(LbUsFourByte index, double energy)
{
	double			probScatter;		/* Probability of coherent scatter */
	LbUsFourByte	energyIndex;		/* Energy index */
	
	do { /* Process Loop  */
		
		/* Convert energy to tissue index */
		energyIndex = (LbUsFourByte) energy;
		
		/* Round up */
		if ((energy - (LbUsFourByte) energy) > .5) {
			energyIndex++;
		}				
		/* Subtract for lowest supported energy */
		energyIndex -= PHG_MIN_PHOTON_ENERGY;
		
		/* Get attenuation for given energy */
		if (PHG_IsModelCoherentInObj()){
			probScatter = SubObjTissueAttenTableCoh[index][energyIndex].probScatter;
		}
		else {
			probScatter = SubObjTissueAttenTableNoCoh[index][energyIndex].probScatter;
		}
		
	} while (false);
	
	return (probScatter);
}

/*********************************************************************************
*
*			Name:		SubObjGetProbScatterInObj
*
*			Summary:	Get probability of scatter for photon's
*						current state.  Note this is actually the cumulative
*						probability of compton + coherent.
*			Arguments:
*			PHG_TrackingPhoton	*trackingPhotonPtr	- Our photon
*
*			Function return: Probability of compton scatter.
*
*********************************************************************************/
double	SubObjGetProbScatterInObj(PHG_TrackingPhoton *trackingPhotonPtr)
{
	double			probScatter;		/* Probability of coherent scatter */
	LbUsFourByte	index;				/* The tissue index */

	do { /* Process Loop  */
		
		/* Get tissue index */
		index = SubObjObject[trackingPhotonPtr->sliceIndex].attenuationArray[(trackingPhotonPtr->yIndex*SubObjObject[trackingPhotonPtr->sliceIndex].attNumXBins)+trackingPhotonPtr->xIndex];
		if (index >= SubObjNumTissues)
			PhgAbort("Index out of bounds for photon position (SubObjGetProbScatterInObj)", true);
			
		
		probScatter = SubObjGetProbScatterInObj2(index, trackingPhotonPtr->energy);
		
	} while (false);
	
	return (probScatter);
}

/*********************************************************************************
*
*			Name:		SubObjGetProbScatterInTomo2
*
*			Summary:	Get probability of scatter given material
*						and photon energy.
*			Arguments:
*			LbUsFourByte	index	- The tissue index.
*			double			energy	- The tissue energy.
*
*			Function return: Probability of compton scatter.
*
*********************************************************************************/
double	SubObjGetProbScatterInTomo2(LbUsFourByte index, double energy)
{
	double			probScatter;		/* Probability of coherent scatter */
	LbUsFourByte	energyIndex;		/* Energy index */
	
	do { /* Process Loop  */
		
		/* Convert energy to tissue index */
		energyIndex = (LbUsFourByte) energy;
		
		/* Round up */
		if ((energy - (LbUsFourByte) energy) > .5) {
			energyIndex++;
		}				
		/* Subtract for lowest supported energy */
		energyIndex -= PHG_MIN_PHOTON_ENERGY;
		
		/* Get attenuation for given energy */
		if (PHG_IsModelCoherentInTomo()){
			probScatter = SubObjTissueAttenTableCoh[index][energyIndex].probScatter;
		}
		else {
			probScatter = SubObjTissueAttenTableNoCoh[index][energyIndex].probScatter;
		}
		
	} while (false);
	
	return (probScatter);
}

/*********************************************************************************
*
*			Name:		SubObjGetProbScatterInTomo
*
*			Summary:	Get probability of scatter for photon's
*						current state.  Note this is actually the cumulative
*						probability of compton + coherent.
*			Arguments:
*			PHG_TrackingPhoton	*trackingPhotonPtr	- Our photon
*
*			Function return: Probability of compton scatter.
*
*********************************************************************************/
double	SubObjGetProbScatterInTomo(PHG_TrackingPhoton *trackingPhotonPtr)
{
	double			probScatter;		/* Probability of coherent scatter */
	LbUsFourByte	index;				/* The tissue index */

	do { /* Process Loop  */
		
		/* Get tissue index */
		index = SubObjObject[trackingPhotonPtr->sliceIndex].attenuationArray[(trackingPhotonPtr->yIndex*SubObjObject[trackingPhotonPtr->sliceIndex].attNumXBins)+trackingPhotonPtr->xIndex];
		if (index >= SubObjNumTissues)
			PhgAbort("Index out of bounds for photon position (SubObjGetProbScatterInTomo)", true);
			
		
		probScatter = SubObjGetProbScatterInTomo2(index, trackingPhotonPtr->energy);
		
	} while (false);
	
	return (probScatter);
}

/*********************************************************************************
*
*			Name:		SubObjGetProbComptToScatInTomo2
*
*			Summary:	Get probability of compton scatter given material
*						and photon energy.
*			Arguments:
*			LbUsFourByte	index	- The tissue index.
*			double			energy	- The tissue energy.
*
*			Function return: Probability of compton scatter.
*
*********************************************************************************/
double	SubObjGetProbComptToScatInTomo2(LbUsFourByte index, double energy)
{
	double			probComptonToScatter;	/* Probability of compton scatter */
	LbUsFourByte	energyIndex;			/* Energy index */

		
	/* Convert energy to tissue index */
	energyIndex = (LbUsFourByte) energy;
	
	/* Round up */
	if ((energy - (LbUsFourByte) energy) > .5) {
		energyIndex++;
	}				
	/* Subtract for lowest supported energy */
	energyIndex -= PHG_MIN_PHOTON_ENERGY;
	
	/* Get attenuation for given energy */
	if (PHG_IsModelCoherentInTomo()) {
		probComptonToScatter = SubObjTissueAttenTableCoh[index][energyIndex].probComptonToScatter;
	}
	else {
		probComptonToScatter = SubObjTissueAttenTableNoCoh[index][energyIndex].probComptonToScatter;
	}
	
	return (probComptonToScatter);
}

/*********************************************************************************
*
*			Name:		SubObjGetProbComptToScatInObj2
*
*			Summary:	Get probability of compton scatter given material
*						and photon energy.
*			Arguments:
*			LbUsFourByte	index	- The tissue index.
*			double			energy	- The tissue energy.
*
*			Function return: Probability of compton scatter.
*
*********************************************************************************/
double	SubObjGetProbComptToScatInObj2(LbUsFourByte index, double energy)
{
	double			probComptonToScatter;	/* Probability of compton scatter */
	LbUsFourByte	energyIndex;			/* Energy index */

		
	/* Convert energy to tissue index */
	energyIndex = (LbUsFourByte) energy;
	
	/* Round up */
	if ((energy - (LbUsFourByte) energy) > .5) {
		energyIndex++;
	}				
	/* Subtract for lowest supported energy */
	energyIndex -= PHG_MIN_PHOTON_ENERGY;
	
	/* Get attenuation for given energy */
	if (PHG_IsModelCoherentInObj()) {
		probComptonToScatter = SubObjTissueAttenTableCoh[index][energyIndex].probComptonToScatter;
	}
	else {
		probComptonToScatter = SubObjTissueAttenTableNoCoh[index][energyIndex].probComptonToScatter;
	}
	
	return (probComptonToScatter);
}


/*********************************************************************************
*
*			Name:		SubObjGetProbComptToScatter
*
*			Summary:	Get ratio of probability of compton scatter to probability of any scatter
*						for photon'current state.
*			Arguments:
*			PHG_TrackingPhoton	*trackingPhotonPtr	- Our photon
*
*			Function return: Probability of compton scatter.
*
*********************************************************************************/
double	SubObjGetProbComptToScatter(PHG_TrackingPhoton *trackingPhotonPtr)
{
	double			probComptonToScatter;	/* Probability of compton scatter */
	LbUsFourByte	index;					/* The tissue index */
	
	
	/* Get tissue index */
	index = SubObjObject[trackingPhotonPtr->sliceIndex].attenuationArray[(trackingPhotonPtr->yIndex*SubObjObject[trackingPhotonPtr->sliceIndex].attNumXBins)+trackingPhotonPtr->xIndex];
	
	#ifdef PHG_DEBUG
	if (index >= SubObjNumTissues)
		PhgAbort("Index out of bounds for photon position (SubObjGetProbComptToScatter)", true);
	#endif
	
	probComptonToScatter = SubObjGetProbComptToScatInObj2(index, trackingPhotonPtr->energy);
	
	
	return (probComptonToScatter);
}

/*********************************************************************************
*
*		Name:		SubObjGetProbScatter
*
*		Summary:	Get probability of scatter, depending on modelingCohScatter, 
*						given material and photon energy.
*
*		Arguments:
*			LbUsFourByte	index				- The tissue index.
*			double			energy				- The photon energy.
*			Boolean			modelingCohScatter	- If coherent scatter is simulated.
*
*		Function return: Probability of a scatter.
*
*********************************************************************************/

double	SubObjGetProbScatter(LbUsFourByte index, double energy, 
								Boolean modelingCohScatter)

{
	LbUsFourByte	energyIndex;		/* Energy index */
	double			probScatter;		/* Probability of scatter */
	
	
	/* Convert energy to tissue index (rounding up) */
	if (energy <= PHG_MIN_PHOTON_ENERGY) {
		energyIndex = 0;
	}
	else {
		energyIndex = (LbUsFourByte)(energy - PHG_MIN_PHOTON_ENERGY + 0.5);
	}
	
	/* Get scatter probability for the given energy */
	if (modelingCohScatter){
		probScatter = SubObjTissueAttenTableCoh[index][energyIndex].probScatter;
	}
	else {
		probScatter = SubObjTissueAttenTableNoCoh[index][energyIndex].probScatter;
	}
	
	return (probScatter);
}

/*********************************************************************************
*
*		Name:		SubObjGetProbComptonCondnl
*
*		Summary:	Get ratio of probability of Compton scatter to probability of any scatter, 
*						depending on modelingCohScatter, given material and photon energy.
*					This is the conditional probability of a Compton scatter, given a scatter.
*
*		Arguments:
*			LbUsFourByte	index				- The tissue index.
*			double			energy				- The photon energy.
*			Boolean			modelingCohScatter	- If coherent scatter is simulated.
*
*		Function return: Conditional probability of a Compton scatter.
*
*********************************************************************************/

double	SubObjGetProbComptonCondnl(LbUsFourByte index, double energy, 
									Boolean modelingCohScatter)

{
	LbUsFourByte	energyIndex;		/* Energy index */
	double			probScatter;		/* Conditional probability of scatter */
	
	
	/* Convert energy to tissue index (rounding up) */
	if (energy <= PHG_MIN_PHOTON_ENERGY) {
		energyIndex = 0;
	}
	else {
		energyIndex = (LbUsFourByte)(energy - PHG_MIN_PHOTON_ENERGY + 0.5);
	}
	
	/* Get conditional scatter probability for the given energy */
	if (modelingCohScatter){
		probScatter = SubObjTissueAttenTableCoh[index][energyIndex].probComptonToScatter;
	}
	else {
		probScatter = SubObjTissueAttenTableNoCoh[index][energyIndex].probComptonToScatter;
	}
	
	return (probScatter);
}

/*********************************************************************************
*
*			Name:			SubObjCalcTimeBinDecays
*
*			Summary:		Calculate decay distribution for current simulation
*							time bin.
*
*			Arguments:
*				LbUsFourByte	currentTimeBin.
*				LbUsFourByte	numDecaysToSimulate;
*
*			Function return: None.
*
*********************************************************************************/
void SubObjCalcTimeBinDecays(LbUsFourByte currentTimeBin,
		LbUsEightByte numDecaysToSimulate)	
{
    double			sliceProductivity;			/* Productivity of a given slice */
    double			estimatedRealDetected;			/* Number of events that would be detected */
    LbUsFourByte	sliceIndex;					/* Current slice */
    LbUsFourByte	angleIndex;					/* Current angle */
    
    
    /* Clear counters */
    sliceProductivity = 0.0;
    estimatedRealDetected = 0.0;
    subObjSimulatedDivRealDetected = 0.0;
    
    /* Save our current time bin */
    subObjCurTimeBin = currentTimeBin;
    
    /* Calculate number of real events that would be detected */
    for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
	
		/* Clear slice productivity accumulator */
		sliceProductivity = 0.0;
	
		/* Calculate slice productivity  */
		for (angleIndex = 0; angleIndex < PRODTBLGetNumAngleCells(); angleIndex++) {
	    
	    	/* Add to sum, productivity of current angle * size of angle */
	    	sliceProductivity +=
				PRODTBLGetMaxCellProductivity(sliceIndex, angleIndex) *
		   		(PRODTBLGetProdTblAngleSize(sliceIndex, angleIndex)/2);
		}
	
		/* Add to sum, (slice activity * decays per curie * 
	   		slice productivity) 
	   	*/
		estimatedRealDetected += (SUBOBJ_DECAYS_PER_CURIE * 
			SUBOBJGetSliceActivity(sliceIndex) 
			* sliceProductivity);
    }
    
   	/* Verify we have a valid productivity situation */
   	if (estimatedRealDetected == 0) {
		PhgAbort("Computed zero for estimated real decays, please check your object (SubObjCalcTimeBinDecays).",
			false);
   	}
    
    /* Calculate subObjSimulatedDivRealDetected (Ratio of desired/real) */
    subObjSimulatedDivRealDetected = numDecaysToSimulate/estimatedRealDetected;
    
	/* Clear count of "round ups" this is for statistical info and error checking */
    SubObjAngleRoundUpCount = 0;

	/* Compute the slice time bin decays for slice zero */
	subObjCalcSliceTimeBinDecays(subObjCurTimeBin, 0);    

    /* Now reset current indexes because we'll be retreiving these decays next */
    SubObjCurSliceIndex = 0;
    SubObjCurVoxelIndex = 0;
    SubObjCurAngleIndex = 0;
	
	/* Get timing information */
	LbTmStartTiming(&SubObjTrackPhotonsStartTime);
}

/*********************************************************************************
*
*			Name:			subObjCalcSliceTimeBinDecays
*
*			Summary:		Calculate decay distribution for current simulation
*							time bin and slice.
*
*			Arguments:
*				LbUsFourByte	currentTimeBin				- The current time bin
*				LbUsFourByte	sliceIndex					- The current slice
*
*			Function return: None.
*
*********************************************************************************/
void subObjCalcSliceTimeBinDecays(LbUsFourByte currentTimeBin,
		LbUsFourByte sliceIndex)
{
    double			sliceActivity;				/* Total slice activity */
    double			voxelActivity;				/* Total voxel activity */
    double			realDecaysInCell;			/* Number of real events that the cell would produce */
    double			realDecaysInVoxelAng;		/* Number of real events that the cell would produce */
    double			simDetectedInCell;			/* Number of simulated events that the cell will produce */
	double			numDecaysForCell;			/* Number of decays to generate for the cell */
	double			voxelSize;					/* Size of the voxel */
    LbUsFourByte	angleIndex;					/* Current angle */
    LbUsFourByte	voxelIndex;					/* Current voxel */
    
    
    #ifdef PHG_DEBUG
    	/* Write random seed value if requested, always resetting to beggining */
    	if (subObjDoWriteRandSeed == true)
    		PhgMathWriteSeed(-1);
    		
    	/* Read random seed value if requested */
    	if (subObjDoReadRandSeed == true)
    		PhgMathReadSeed(-1);
    		
    #endif
    	    
	/* Get slice activity */
	sliceActivity = SUBOBJGetSliceActivity(sliceIndex);
	
	/* Calculate voxel size */
	voxelSize = SubObjObject[sliceIndex].actVoxelWidth *
		SubObjObject[sliceIndex].actVoxelHeight *
		SubObjObject[sliceIndex].sliceDepth;
	
	/* Compute decays for each voxel */
	for (voxelIndex = 0; voxelIndex < SubObjGetNumActVoxels(sliceIndex); voxelIndex++) {
    
    	/* Compute total voxel activity */
    	voxelActivity = 
			SubObjGetTissueActivity(SubObjObject[(sliceIndex)].activityArray[voxelIndex], 
			currentTimeBin)
	    	* voxelSize * SubObjCurTimeBinDuration;
    
    	/* Distribute voxel activity throughout angles */
    	for (angleIndex = 0; angleIndex < PRODTBLGetNumAngleCells(); angleIndex++) {
	
			/* Calc voxel/cell number of decays */
			if ((sliceActivity > 0) && (voxelActivity > 0.0)) {
	    
	   	 		/* Get real decays for angle cell */
	    		realDecaysInCell = subObjCellGetNumReal(sliceIndex,
					angleIndex, SUBOBJGetSliceActivity(sliceIndex));
	    
	    		/* Calculate simulated decays for angle cell */
	    		simDetectedInCell =  subObjSimulatedDivRealDetected * realDecaysInCell *
					PRODTBLGetMaxCellProductivity(sliceIndex, angleIndex);
	    
	    		/* Calculate number to simulate for this angle cell of the voxel */
				numDecaysForCell = ((voxelActivity/sliceActivity) * 
					simDetectedInCell);
				
				/* Get real decays for the current voxel/angle */	
				realDecaysInVoxelAng = subObjVoxelCellGetNumReal(sliceIndex,  angleIndex,
						voxelIndex, 0);
						
				/* Assign if greater than zero */
				if (numDecaysForCell > 0.0){
				
					/* Truncate number of decays for the angle cell */
					SubObjDecaySlice[(voxelIndex*PRODTBLGetNumAngleCells()) + angleIndex] =
						(LbUsEightByte) numDecaysForCell;

					/* Russian Roulette any fractional decays 
					 (leftover from numDecaysForCell) */
					if ( (numDecaysForCell - SubObjDecaySlice[(voxelIndex*PRODTBLGetNumAngleCells()) + angleIndex])
							 > PhgMathGetRandomNumber() ) {

						/* Increment number of decays for the angle cell by one */
						SubObjDecaySlice[(voxelIndex*PRODTBLGetNumAngleCells()) + angleIndex] +=
							1;

						SubObjAngleRoundUpCount++;
					}
	
					/* Set decay weight to expected real decays divided by simulated decays */
					SubObjDecayWeightSlice[(voxelIndex*PRODTBLGetNumAngleCells()) + angleIndex] = 
						realDecaysInVoxelAng / numDecaysForCell;
	    		}

	    		/* the following branch should never be used, but is here as safety precaution */
	    		else {
		    		SubObjDecaySlice[(voxelIndex*PRODTBLGetNumAngleCells()) + angleIndex] = 0;
					SubObjDecayWeightSlice[(voxelIndex*PRODTBLGetNumAngleCells()) + angleIndex] = 0.0;
			    }

   			}
			else {
	    		SubObjDecaySlice[(voxelIndex*PRODTBLGetNumAngleCells()) + angleIndex] = 0;
				SubObjDecayWeightSlice[(voxelIndex*PRODTBLGetNumAngleCells()) + angleIndex] = 0.0;
			}
		}
    
	} /* Loop to next voxel */
}

/*********************************************************************************
*
*			Name:			subObjCellGetNumReal
*
*			Summary:		Get number of real decays for a given angle cell
*							and time.
*
*			Arguments:
*				LbUsFourByte	sliceIndex		- Current slice.
*				LbUsFourByte	angleIndex		- Current angle.
*				double			sliceActivity	- Activity of slice.
*
*			Function return: Number of real events to occur at the given location.
*
*********************************************************************************/
double subObjCellGetNumReal(LbUsFourByte sliceIndex, LbUsFourByte angleIndex,
		double sliceActivity)	
{
	double			numRealDecays;	/* Number of real decays for a given voxel */
	double			tempCount;		/* Temporary double value for calculation */
	
	/* Perform 1st calculation */
	tempCount = ((SUBOBJ_DECAYS_PER_CURIE * sliceActivity) *
		PRODTBLGetProdTblAngleSize(sliceIndex, angleIndex));

	/* Angles range from [-1,1], so angle sizes sum to 2 */ 
	numRealDecays = tempCount/2;
		
	return (numRealDecays);
}

/*********************************************************************************
*
*			Name:		subObjAllocDecaySlice
*
*			Summary:	Allocate memory for decays of the current slice.
*
*			Arguments:
*				LbUsFourByte	sliceIndex	- The current slice
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean	subObjAllocDecaySlice(LbUsFourByte sliceIndex)
{
	Boolean			okay = false;			/* Success flag */
	LbUsFourByte	numVoxels;				/* Number of voxels in current slice */
	LbUsFourByte	numPrevVoxels;			/* Number of voxels in previous slice */
	
	do { /* Process Loop */
		
		/* Get number of voxels in this slice */
		numVoxels = SubObjGetNumActVoxels(sliceIndex);

		/* If sliceIndex is not zero, then we will free the slice
			if the number of voxels is different then previous
			slice, otherwise we can skip this whole step
		*/
		
		if (sliceIndex != 0) {

			/* Get voxels in previous slice */
			numPrevVoxels = SubObjGetNumActVoxels(sliceIndex-1);
			
		}
		else {
			
			/* Set number of voxels in previous slice to zero, forcing
				us to allocate memory for the first slice
			*/
			numPrevVoxels = 0;
			
		}
		
		/* if number of voxels in previous slice is not equal to number in this slice
			then process the slice */
		if (numPrevVoxels != numVoxels) {
		
			/* Free voxels from previous slice if there was one */
			if (numPrevVoxels != 0) {
							
				/* Free the decay slice and decay weight slice */
				LbMmFree((void **)&(SubObjDecaySlice));
				LbMmFree((void **)&(SubObjDecayWeightSlice));

			}			
			
			/* Allocate decay slice array */
			if ((SubObjDecaySlice = (subObjDecaySliceTy *) LbMmAlloc(
					(sizeof(subObjDecaySliceTy) * numVoxels * 
					PRODTBLGetNumAngleCells()))) == 0) {
				
				goto FAIL;
			}
			
			/* Allocate decay weight slice array */
			if ((SubObjDecayWeightSlice = (subObjDecayWeightSliceTy *) LbMmAlloc(
					(sizeof(subObjDecayWeightSliceTy) * numVoxels * 
					PRODTBLGetNumAngleCells()))) == 0) {
				
				goto FAIL;
			}
	
		}
		okay = true;
		FAIL:;
	} while (false);
	
	/* If we failed, free any memory allocated because it won't get used */
	if (!okay) {
		if (SubObjDecaySlice != 0) {

			/* Free the decay slice and decay weight slice */
			LbMmFree((void **)&(SubObjDecaySlice));
			
		}
		if (SubObjDecayWeightSlice != 0) {
			
			/* Free the decay slice and decay weight slice */
			LbMmFree((void **)&(SubObjDecayWeightSlice));
		}
	}
	
	return(okay);
}

/*********************************************************************************
*
*			Name:		SubObjGetCohTheta2
*
*			Summary:	Return scatter angle for coherent scatter.
*
*			Arguments:
*				LbFourByte	materialIndex	- What material the photon is traveling in
*				double		energy			- The photon energy
*
*			Function return: The cosine of the scatter angle.
*
*********************************************************************************/
double	SubObjGetCohTheta2(LbUsFourByte materialIndex, double energy)
{
	double cosTheta;
	LbFourByte	adIndex;
	LbFourByte	energyIndex = floor(energy);
	LbFourByte	energy1, energy2;
	double		cosThetaE1, cosThetaE2;
	double		adValue;
	
	/* Verify our material is supported for coherent scatter */
	if (materialIndex >= subObjNumCohMaterials) {
		sprintf(subObjErrStr, "Attempted to compute coherent scatter for unsupported material, index = %ld",
			(unsigned long)materialIndex);
		
		PhgAbort(subObjErrStr, true);
	}
	
	/* Convert energy index into table index */
	if (energyIndex <= SUBOBJ_NUM_1KEV_ENERGIES) {

		energy1 = energyIndex-1;
		energy2 = energyIndex;

	}	/* See if energy is in 10 keV bin range */
	else if ((energyIndex > SUBOBJ_NUM_1KEV_ENERGIES) && (energyIndex < SUBOBJ_MAX_10KEV_ENERGIES)) {
	
		energy1 = ((energyIndex - SUBOBJ_NUM_1KEV_ENERGIES) * SUBOBJ_NUM_10KEV_ENERGIES)/
			(SUBOBJ_NUM_10KEV_ENERGIES * 10);
		energy1 += (SUBOBJ_NUM_1KEV_ENERGIES-1);
		
		energy2 = energy1+1;
	}
	else { /* Energy must  be in 100 keV bin range */
		energy1 = ((energyIndex - (SUBOBJ_MAX_10KEV_ENERGIES))*SUBOBJ_NUM_100KEV_ENERGIES)/
			(SUBOBJ_NUM_100KEV_ENERGIES*100);
		energy1 += ( SUBOBJ_NUM_1KEV_ENERGIES + SUBOBJ_NUM_10KEV_ENERGIES) - 1;
		
		energy2 = energy1+1;
	}
	
	/* Now compute theta */
	
	/* The coherent scattering angle tables are indexed by cumulative probability*/
	/* Get a random cumulative probability index */
	adValue = SUBOBJ_NUM_COH_ANGLES * PhgMathGetRandomNumber();
	adIndex = (LbFourByte) floor(adValue);

	/* Use linear interpolation to choose scatter angles for the two bracketing energy bins */
	if (adIndex == 0) {
		/* For index 0 we interpolate between 1, which is the value for cumulative
			probability=0, and */
		/* the first value in the table, which has cumulative probability of
			1/SUBOBJ_NUM_COH_ANGLES */
		cosThetaE1 =
			adValue * subObjCohScatAngles[materialIndex][energy1].angleProbabilities[0]
				+
				(1 - adValue);
			
			cosThetaE2 =
				adValue * subObjCohScatAngles[materialIndex][energy2].angleProbabilities[0]
				+
				(1 - adValue);
	}
	else if (adIndex < SUBOBJ_NUM_COH_ANGLES) {
		/* This is the more typical interpolation between an entry and the previous
			entry */
		cosThetaE1 =
			(adValue -
			adIndex)*subObjCohScatAngles[materialIndex][energy1].angleProbabilities[adIndex]
			+
			(adIndex+1 -
			adValue)*subObjCohScatAngles[materialIndex][energy1].angleProbabilities[adIndex-1];
	
		cosThetaE2 =
			(adValue -
			adIndex)*subObjCohScatAngles[materialIndex][energy2].angleProbabilities[adIndex]
			+
			(adIndex+1 -
			adValue)*subObjCohScatAngles[materialIndex][energy2].angleProbabilities[adIndex-1];
	}
	else {
		/* exceptional case where adIndex points past end of table--set cosTheta to
			the last value, -1 */ 
		cosThetaE1 = -1;
		cosThetaE2 = -1;
	}
		
	/* Now use linear interpolation between the two energy bins */
	cosTheta = (((energy -
		subObjCohScatAngles[materialIndex][energy1].energy)*cosThetaE2) +
		((subObjCohScatAngles[materialIndex][energy2].energy-energy)*cosThetaE1))/
		(subObjCohScatAngles[materialIndex][energy2].energy -
		subObjCohScatAngles[materialIndex][energy1].energy);

	#ifdef PHG_DEBUG
	if ((cosTheta < -1.0) || (cosTheta > 1.0)) {
		ErAbort("Computed invalid cosTheta in (SubObjGetCohTheta2)");
	}
	#endif
	
	return(cosTheta);
}

/*********************************************************************************
*
*			Name:		SubObjGetCohTheta
*
*			Summary:	Return scatter angle for coherent scatter.
*
*			Arguments:
*
*			Function return: Cosine of the scatter angle.
*
*********************************************************************************/
double	SubObjGetCohTheta(PHG_TrackingPhoton *trPhoton)
{
	LbUsFourByte	index;				/* The tissue index */
	double			cosTheta;
	
	/* Get tissue index */
	index = SubObjObject[trPhoton->sliceIndex].attenuationArray[(trPhoton->yIndex*SubObjObject[trPhoton->sliceIndex].attNumXBins)+trPhoton->xIndex];
	
	#ifdef PHG_DEBUG
	if (index >= subObjNumCohMaterials)
		PhgAbort("Index out of bounds for photon position (SubObjGetCohTheta)", true);
	#endif
	
	cosTheta = SubObjGetCohTheta2(index, trPhoton->energy);
	
	return(cosTheta);
}

/*********************************************************************************
*
*			Name:		subObjInitCoherentTbl
*
*			Summary:	Reads in the data for the coherent scatter angle table.
*
*			Arguments:
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean	subObjInitCoherentTbl()
{
	Boolean			okay = false;			/* Process flat */
	char			inputBuffer[256];		/* Enough to read in lines of material info */
	char			fileNameBuffer[256];	/* Enough to read in lines of material info */
	char			kevStr[1024];
	char			materialStr[1024];
	char			theMaterialStr[1024];
	FILE			*cohFile = 0;		/* Lists of angular distribution files */
	FILE			*adFile = 0;		/* Angular distribution file */
	LbFourByte		materialIndex;		/* LCV for materials */
	LbUsFourByte	lineCount = 0;		/* LCV for angle probabilities */
	LbFourByte		energyIndex = -1;	/* LCV for angle probabilities */
	LbFourByte		angleIndex = 0;	/* LCV for angle probabilities */
	double			angleProbability;	/* Data for angle probabilities */
	LbUsFourByte	energy;				/* Data for angle probabilities */
	unsigned long	energyUSL;
	LbUsFourByte	numItems;			/* Error checking for text conversion */
	char			*materialName;		/* Name of material within path */
	do { /* Process Loop */
						
		/* Open the coherent scatter angle file */
		if ((strlen(PhgRunTimeParams.PhgSubObjCohScatTableFilePath) == 0) ||
				(cohFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjCohScatTableFilePath, "r")) == 0) {
			
			sprintf(subObjErrStr, "Unable to open coherent scatter angle file list '%s' (SubObjInitCoherentTbl).",
					PhgRunTimeParams.PhgSubObjCohScatTableFilePath);
			ErStFileError(subObjErrStr);
			break;
		}

		/* Clear counter */
		subObjNumCohMaterials = 0;
		
		/* Loop through and count the names */
		while (LbFlFGetS(inputBuffer, sizeof(inputBuffer)-1, cohFile) != NULL) {
		
			/* Strip off the new line from LbFlFGetS */
			inputBuffer[strlen(inputBuffer)-1] = '\0';
		
			/* If this is a blank line, assume we are at the end of file and have extra line feeds */
			if (inputBuffer[0] == '\0')
				break;
				
			/* Peel off the material name */
			if ((materialName = strrchr(inputBuffer, PATH_SEPARATOR)) == 0) {
				*(strstr(inputBuffer, ".ad")) = '\0';
				strcpy(convertCohMaterialLabels[subObjNumCohMaterials], inputBuffer);
			}
			else {
				/* Skip past separator and copy name */
				materialName += 1;
				*(strstr(materialName, ".ad")) = '\0';
				strcpy(convertCohMaterialLabels[subObjNumCohMaterials], materialName);
			}
			
			subObjNumCohMaterials++;
			
			if (subObjNumCohMaterials == MAX_NUM_MATERIALS) {
				ErStGeneric("You have too many materials in your table, edit SubObj.c and increase MAX_NUM_MATERIALS to support more materials.(SubObjInitCoherentTbl)");
				break;
			}
		}
		
		/* Allocated the table */
		if ((subObjCohScatAngles = (subObjCoScatAngleTblTy *) LbMmAlloc(sizeof(subObjCoScatAngleTblTy) * subObjNumCohMaterials)) == 0) {
			ErStGeneric("Unable to allocate memory for coherent scatter angle table (SubObjInitCoherentTbl)");
			break;
		}

		/* Go back to beginning of file */
		if (fseek(cohFile, 0, SEEK_SET) != 0){
			ErStFileError("Unable to seek to beginning of coherent scatter angle file list  (SubObjInitCoherentTbl).");
			goto FAIL;
		}
		
		materialIndex = 0;
		
		/* Loop through and process the files */
		while (LbFlFGetS(fileNameBuffer, sizeof(inputBuffer)-1, cohFile) != NULL) {

			/* If there is a blank line, skip over it */
			if (strlen(fileNameBuffer) <= 1)
				continue;
				
			/* Strip the new line */
			fileNameBuffer[strlen(fileNameBuffer)-1] = '\0';
			
			/* Open the material file */
			if ((adFile = LbFlFileOpen(fileNameBuffer, "rb")) == 0) {
				sprintf(subObjErrStr, "Unable to create/open angular distribution file\n'%s'\n"
				"  Edit file '%s' to specify a valid path (SubObjInitCoherentTbl).", fileNameBuffer, PhgRunTimeParams.PhgSubObjCohScatTableFilePath);
				ErStFileError(subObjErrStr);
				goto FAIL;
			}	
			
			/* Read through material file and convert text to binary table */
			while (LbFlFGetS(inputBuffer,sizeof(inputBuffer)-1, adFile) != NULL) {

				/* If there is a blank line, assume it is extra line feeds at the end of file and terminate */
				if (strlen(inputBuffer) <= 1) {
					sprintf(subObjErrStr, "Found unexpected blank line in '%s'"
						"\nIt may be at the end of the file, regardless please edit and remove the blank lines  (SubObjInitCoherentTbl)",
						fileNameBuffer);
					ErStGeneric(subObjErrStr);
					goto FAIL;
				}
			
				/* If it's the first line, then convert text to header info */
				if ((lineCount == 0) || ((lineCount) % (SUBOBJ_NUM_ANGLES+1) == 0)) {
					
					/* Reset energy counter */
					energyIndex++;
					angleIndex = 0;
					/* Convert values */
					if ((numItems = sscanf(inputBuffer, "%ld %s %s = %s",
							&energyUSL, kevStr, materialStr, theMaterialStr)) != 4){
							
							sprintf(subObjErrStr,"Error from sscanf on header info, string = '%s'\n"
								" expected something like 1 keV, material = air. (SubObjInitCoherentTbl)",
								inputBuffer);
							ErStGeneric(subObjErrStr);
							goto FAIL;
					}
					energy = (LbUsFourByte)energyUSL;
					
					/* Get the material index */
					if ((materialIndex = SubObjGtMaterialIndex(theMaterialStr)) == -1) {
						goto FAIL;
					}
					
					/* Check indexes */
					if (materialIndex >= (LbFourByte)subObjNumCohMaterials) {
						sprintf(subObjErrStr,"Got too large of material index, %ld (SubObjInitCoherentTbl)",
							(unsigned long)materialIndex);
						ErStGeneric(subObjErrStr);
						goto FAIL;
					}
					
					/* Check indexes */
					if (energyIndex >= SUBOBJ_NUM_COH_ENERGIES) {
						sprintf(subObjErrStr,"Got too large of energyIndex, %ld (SubObjInitCoherentTbl)",
							(unsigned long)energyIndex);
						ErStGeneric(subObjErrStr);
						goto FAIL;
					}
					subObjCohScatAngles[materialIndex][energyIndex].energy = energy;
					subObjCohScatAngles[materialIndex][energyIndex].materialIndex = materialIndex;
					
				}
				else {

					/* Convert values */
					if ((sscanf(inputBuffer, "%lf",
							&angleProbability)) != 1){
							
							sprintf(subObjErrStr,"Error from sscanf on angle probability, string = '%s'\n"
								" expected something like .00874. (SubObjInitCoherentTbl)",
								inputBuffer);
							ErStGeneric(subObjErrStr);
							goto FAIL;
					}
					
					/* Check indexes */
					if (angleIndex >= SUBOBJ_NUM_ANGLES) {
						sprintf(subObjErrStr,"Got too large of angle index, %ld (SubObjInitCoherentTbl)",
							(unsigned long)angleIndex);
						ErStGeneric(subObjErrStr);
						goto FAIL;
					}
					subObjCohScatAngles[materialIndex][energyIndex].angleProbabilities[angleIndex] = angleProbability;
					angleIndex++;
					
				}
				lineCount++;
			}

		
			fclose(adFile);
			adFile = 0;
			energyIndex = -1;
			materialIndex++;
		}
		
		okay = true;
		FAIL:;
	} while (false);
	
	if (adFile != 0) {
		fclose(adFile);
		adFile = 0;
	}
		
	return (okay);
}
/**********************
*	SubObjGtMaterialIndex
*
*	Purpose:	Convert text string to material index.
*	Arguments:
*		char	* theMaterialStr - text label of material string
*	Result:	The index of the material.
***********************/

LbFourByte	SubObjGtMaterialIndex(char *theMaterialStr)
{
	LbUsFourByte	index;
	LbFourByte		materialIndex = -1;
	
	#ifdef PHG_DEBUG
		if ((PHG_IsModelCoherentInObj() == false) && (PHG_IsModelCoherentInTomo() == false)) {
			PhgAbort("\nYou are trying to call SubObjGtMaterialIndex when coherent scatter is not being modelled\n", false);
		}
	#endif
	
	/* Search table for given label */
	for (index = 0; index < subObjNumCohMaterials; index++){
	
		/* See if it matches */
		if (strcmp(theMaterialStr, convertCohMaterialLabels[index]) == 0) {
			materialIndex = index;
			break;
		}
	}

	/* If not found, register an error */
	if (materialIndex == -1){
		sprintf(subObjErrStr, "Unable to find '%s' in material table (SubObjGtMaterialIndex)",
			theMaterialStr);
			
		ErStGeneric(subObjErrStr);
		ErAbort("Failed to find attenuation material in material list");
	}
	return(materialIndex);
}


/*********************************************************************************
*
*			Name:		subObjInitCoherentBinaryTbl
*
*			Summary:	Reads in the data for the coherent scatter angle table.
*
*			Arguments:
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean	subObjInitCoherentBinaryTbl(void)
{
	Boolean			okay = false;		/* Process flat */
	char			inputBuffer[256];	/* Enough to read in lines of material info */
	FILE			*cohFile = 0;		/* Lists of angular distribution files */
	FILE			*adFile = 0;		/* Angular distribution file */
	LbUsFourByte	materialIndex;		/* LCV for materials */
	
	do { /* Process Loop */
						
		/* Open the coherent scatter angle file */
		if ((strlen(PhgRunTimeParams.PhgSubObjCohScatTableFilePath) == 0) ||
				(cohFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjCohScatTableFilePath, "r")) == 0) {
			
			sprintf(subObjErrStr, "Unable to open coherent scatter angle file list '%s' (subObjInitCoherentBinaryTbl).",
					PhgRunTimeParams.PhgSubObjCohScatTableFilePath);
			ErStFileError(subObjErrStr);
			break;
		}

		/* Clear counter */
		subObjNumCohMaterials = 0;
		
		/* Loop through and count the names */
		while (LbFlFGetS(inputBuffer, sizeof(inputBuffer), cohFile) != NULL) {
			subObjNumCohMaterials++;
		}
		
		/* Allocated the table */
		if ((subObjCohScatAngles = (subObjCoScatAngleTblTy *) LbMmAlloc(sizeof(subObjCoScatAngleTblTy) * subObjNumCohMaterials)) == 0) {
			ErStGeneric("Unable to allocate memory for coherent scatter angle table (SubObjInitCoherentTbl)");
			break;
		}

		/* Go back to beginning of file */
		if (fseek(cohFile, 0, SEEK_SET) != 0){
			ErStFileError("Unable to seek to beginning of coherent scatter angle file list  (SubObjInitCoherentTbl).");
			goto FAIL;
		}
		
		materialIndex = 0;
		
		/* Loop through and count the names */
		while (LbFlFGetS(inputBuffer, sizeof(inputBuffer), cohFile) != NULL) {

			/* If there is a blank line, skip over it */
			if (strlen(inputBuffer) <= 1)
				continue;
				
			/* Strip the new line */
			inputBuffer[strlen(inputBuffer)-1] = '\0';
			
			/* Open the material file */
			if ((adFile = LbFlFileOpen(inputBuffer, "rb")) == 0) {
				sprintf(subObjErrStr, "Unable to create/open angular distribution file\n'%s'\n"
				"  Edit file '%s' to specify a valid path (subObjInitCoherentBinaryTbl).", inputBuffer, PhgRunTimeParams.PhgSubObjCohScatTableFilePath);
				ErStFileError(subObjErrStr);
				goto FAIL;
			}		

			/* Read the angular distribution file */
			if (fread(subObjCohScatAngles[materialIndex], sizeof(subObjCoScatAngleTblTy), 1, adFile) != 1) {
				sprintf(subObjErrStr, "Error reading angular distribution data from file '%s' (subObjInitCoherentBinaryTbl)",
					inputBuffer);
				ErStFileError(subObjErrStr);
				goto FAIL;
			}
		
			fclose( adFile );
			adFile = 0;
			
			materialIndex++;
		}
		
		okay = true;
		FAIL:;
	} while (false);
	
	if (adFile != 0) {
		fclose(adFile);
		adFile = 0;
	}
		
	return (okay);
}

/*********************************************************************************
*
*			Name:		SubObjCreate
*
*			Summary:	Convert the scan_object from the Object Editor into forms usable
*						by the photon history generator.
*
*			Arguments:
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean	SubObjCreate()
{
	Boolean			okay = false;				/* Process flag */
	char			inputBuffer[256];			/* Enough to read in lines of attenuation info */
	char			errString[1024];				/* Error information */
	char			materialString[256];
	char			materialD[32];
	char			materialA[32];
	char			materialZ[32];
	double			voxelAttenuation;			/* Voxel attenuation */
	FILE			*attenuationFile = 0;		/* Attenuation info file */
	FILE            *activityFile= 0;     	 	/* Activity info file */
	FILE			*activityTransFile = 0;		/* Activity index translation file */
	FILE			*attenuationTransFile = 0;	/* Attenuation index translation file */
	FILE			*actImgFile = 0;			/* Activity image file */
	FILE			*attImgFile = 0;			/* Attenuation image file */
	float			voxelActivity;				/* Activity of current voxel. */
	float			voxelAttAsFloat;			/* Attenuation converted to float */
	LbUsFourByte	sliceIndex;					/* Current slice */
	LbUsFourByte	xIndex;						/* Current x index */
	LbUsFourByte	yIndex;						/* Current y index */
	LbUsFourByte	attenIndex;					/* Attenuation coefficient index */
	LbUsFourByte	tissueIndex;				/* Attenuation tissue index */
	unsigned long	tissueIndexUSL;
	LbUsFourByte	loopIndex;					/* Loop control variable */
	LbUsFourByte	*attenuationTransTbl = 0;	/* Attenuation index translation table */
	LbUsFourByte	*activityTransTbl = 0;		/* Activity index translation table */
	LbUsFourByte	tissueIndexTranslation;		/* Translated value of tissue index */
	unsigned long	tissueIndexTranslationUSL;
	LbUsFourByte	numColumns;					/* Used for testing attenuation table */
	LbUsFourByte	attInTranI;					/* Used for parsing index translations */
	LbUsFourByte	materialDAZIndex;			/* Used to parse material string */
	LbUsFourByte	materialStringLen;
	LbUsFourByte	numberI;
	
	do { /* Process Loop */

		/* Initialize coherent scatter table  if it's being modelled */
		if (PHG_IsModelCoherentInObj() || PHG_IsModelCoherentInTomo())
			if (!subObjInitCoherentTbl())
				break;
		
		/* Create the  object */
		{
			/* Clear current slice counter */
	    	SubObjCurSliceIndex = 0;
			
			/* Call subroutine that reads activity and attenuation indexes */
			if (!subObjCreateObject())
				break;
	    	
			/* Restore current slice counter */
			SubObjCurSliceIndex = 0;
		}
		
		/* Allocate decay tables for first slice */
		if (subObjAllocDecaySlice(SubObjCurSliceIndex) == false) {
			break;
		}

		/* Create the tissue table */
		{
							
			/* Open the tissue activity file */
			if ((activityFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjActivityTableFilePath, "r")) == 0) {
				sprintf(errString, "Unable to open activity description file '%s'.",
						PhgRunTimeParams.PhgSubObjActivityTableFilePath);
				ErStFileError(errString);
				goto FAIL;
			}

			/* Get the first line of the file, this is the number of activity indexes */
			if (LbFlFGetS(inputBuffer, sizeof(inputBuffer), activityFile) == 0) {
				ErStGeneric("Error from LbFlFGetS while reading number of tissue indexes.");
				goto FAIL;
			}
		
			/* Set table size to supported count */
			SubObjNumActIndexes = atol(inputBuffer);
			
			/* Allocate space for the table */
			if ((SubObjTissueTable.tissueValues = (SubObjTissueArrayTy) LbMmAlloc(
					sizeof(SubObjTissueElemTy) * SubObjNumActIndexes)) == 0){
			
				goto FAIL;
			}
			
			/* Initialize the table */			
			/* Get input a line at a time and convert it */
			for (tissueIndex = 0; tissueIndex < SubObjNumActIndexes; tissueIndex++) {
					
				/* Get tissue name */
				if (LbFlFGetS(inputBuffer, sizeof(inputBuffer), activityFile) == 0) {
					sprintf(errString, "Error from LbFlFGetS while reading line number %ld"
						"\nfrom tissue activity file named: %s", (unsigned long)(tissueIndex+1),
						PhgRunTimeParams.PhgSubObjActivityTableFilePath);
						
					ErStGeneric(errString);
					goto FAIL;
				}

				SubObjTissueTable.tissueValues[tissueIndex].totalDuration = SUBOBJ_CONST_ACT;
				SubObjTissueTable.tissueValues[tissueIndex].numberOfBins = 1;
				if ((SubObjTissueTable.tissueValues[tissueIndex].activityValues = 
						(SubObjTACArrayTy) LbMmAlloc(sizeof(SubObjTACBinTy))) == 0) {
						
					goto FAIL;
				}
				
				/* Initialize individual activity arrays */
				SubObjTissueTable.tissueValues[tissueIndex].activityValues[0].activityLevel
					= atof(inputBuffer);
				SubObjTissueTable.tissueValues[tissueIndex].activityValues[0].duration
					= SUBOBJ_CONST_ACT;
			}			
		}

		/* Create the tissue attenuation table */
		{			
			/* Open the file */
			if ((attenuationFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjAttenTableFilePath, "r")) == 0) {
				sprintf(errString, "Unable to open attenuation description file '%s'.",
						PhgRunTimeParams.PhgSubObjAttenTableFilePath);
				ErStFileError(errString);
				goto FAIL;
			}
		
			/* Get the number of tissue types */
			if (LbFlFGetS(inputBuffer, sizeof(inputBuffer), attenuationFile) == 0) {
				ErStGeneric("Error from LbFlFGetS while reading number of tissue indexes.");
				goto FAIL;
			}
		
			/* Set table size to supported count */
			SubObjNumTissues = atol(inputBuffer);
			
			/* Allocate space for the table */
			if ((SubObjTissueAttenTableNoCoh = (subObjTissueAttenTblTy) LbMmAlloc(
					sizeof(subObjTissueAttenTy) * SubObjNumTissues)) == 0){
			
				break;
			}
			if ((SubObjTissueAttenTableCoh = (subObjTissueAttenTblTy) LbMmAlloc(
					sizeof(subObjTissueAttenTy) * SubObjNumTissues)) == 0){
			
				break;
			}

				/* Allocate space for tissue name table */
				if ((subObjMaterialNames = (subObjNameTy *) LbMmAlloc(
						sizeof(subObjNameTy) * SubObjNumTissues)) == 0){
				
					break;
				}

				/* Allocate space for tissue density, weight, and number table */
				if ((subObjMaterialDAZ = (subObjMaterialDAZTy *) LbMmAlloc(
						sizeof(subObjMaterialDAZTy) * SubObjNumTissues)) == 0){
				
					break;
				}

			/* Get input a line at a time and convert it */
			for (tissueIndex = 0; tissueIndex < SubObjNumTissues; tissueIndex++) {
	
				/* Get tissue name, density, weight, and number */
				if (LbFlFGetS(materialString, sizeof(materialString), attenuationFile) == 0) {
					ErStGeneric("Error from LbFlFGetS while reading tissue seperator in tissue attenuation file.");
					goto FAIL;
				}
				
				/* Get length of string for boundary checking */
				materialStringLen = strlen(materialString);
				
				/* Copy up to first space into name portion */
				materialDAZIndex = 0;
				do {
					subObjMaterialNames[tissueIndex][materialDAZIndex] = materialString[materialDAZIndex];
					materialDAZIndex++;
				} while ((materialDAZIndex < materialStringLen) && (materialString[materialDAZIndex] != ' '));
				subObjMaterialNames[tissueIndex][materialDAZIndex] = '\0';
				
				/* If we have more on the line then try to convert to D,Z,A values */
				/* Do error checking */
				if (materialDAZIndex != materialStringLen) {
					
					/* Go to the equal sign */
					do {
						materialDAZIndex++;
					} while ((materialDAZIndex < materialStringLen) && (materialString[materialDAZIndex] != '='));
					
					/* Do error checking */
					if (materialDAZIndex == materialStringLen) {
						sprintf(errString,"Your attenuation data file doesn't contain right information (format) for the material, it may be outdated, it needs to contain Density, Atomic Weight and Atomic Number.\n it currently contains '%s'\n",  materialString);
						ErStGeneric(errString);
						goto FAIL;
					}
					
					/* Go to the number */
					do {
						materialDAZIndex++;
					} while ((materialDAZIndex < materialStringLen) && (materialString[materialDAZIndex] == ' '));
					
					/* Do error checking */
					if (materialDAZIndex == materialStringLen) {
						sprintf(errString,"Your attenuation data file doesn't contain right information (format) for the material, it may be outdated, it needs to contain Density, Atomic Weight and Atomic Number.\n it currently contains '%s'\n",  materialString);
						ErStGeneric(errString);
						goto FAIL;
					}
					
					/* Copy the number */
					numberI = 0;
					do {
						materialD[numberI] = materialString[materialDAZIndex];
						materialDAZIndex++;
						numberI++;
					} while ((materialDAZIndex < materialStringLen-1) && (materialString[materialDAZIndex] != ' '));
					materialD[numberI] = '\0';
					subObjMaterialDAZ[tissueIndex].D = atof(materialD);
					
					/* Go to the equal sign */
					do {
						materialDAZIndex++;
					} while ((materialDAZIndex < materialStringLen) && (materialString[materialDAZIndex] != '='));
					
					/* Do error checking */
					if (materialDAZIndex == materialStringLen) {
						sprintf(errString,"Your attenuation data file doesn't contain right information (format) for the material, it may be outdated, it needs to contain Density, Atomic Weight and Atomic Number.\n it currently contains '%s'\n",  materialString);
						ErStGeneric(errString);
						goto FAIL;
					}
					
					/* Go to the number */
					do {
						materialDAZIndex++;
					} while ((materialDAZIndex < materialStringLen) && (materialString[materialDAZIndex] == ' '));
					
					/* Do error checking */
					if (materialDAZIndex == materialStringLen) {
						sprintf(errString,"Your attenuation data file doesn't contain right information (format) for the material, it may be outdated, it needs to contain Density, Atomic Weight and Atomic Number.\n it currently contains '%s'\n",  materialString);
						ErStGeneric(errString);
						goto FAIL;
					}
					
					/* Copy the number */
					numberI = 0;
					do {
						materialA[numberI] = materialString[materialDAZIndex];
						materialDAZIndex++;
						numberI++;
					} while ((materialDAZIndex < materialStringLen-1) && (materialString[materialDAZIndex] != ' '));
					materialA[numberI] = '\0';
					subObjMaterialDAZ[tissueIndex].A = atof(materialA);
					
					/* Do error checking */
					if (materialDAZIndex == materialStringLen) {
						sprintf(errString,"Your attenuation data file doesn't contain right information (format) for the material, it may be outdated, it needs to contain Density, Atomic Weight and Atomic Number.\n it currently contains '%s'\n",  materialString);
						ErStGeneric(errString);
						goto FAIL;
					}
					
					/* Go to the equal sign */
					do {
						materialDAZIndex++;
					} while ((materialDAZIndex < materialStringLen) && (materialString[materialDAZIndex] != '='));
					
					/* Do error checking */
					if (materialDAZIndex == materialStringLen) {
						sprintf(errString,"Your attenuation data file doesn't contain right information (format) for the material, it may be outdated, it needs to contain Density, Atomic Weight and Atomic Number.\n it currently contains '%s'\n",  materialString);
						ErStGeneric(errString);
						goto FAIL;
					}
					
					/* Go to the number */
					do {
						materialDAZIndex++;
					} while ((materialDAZIndex < materialStringLen) && (materialString[materialDAZIndex] == ' '));
					
					/* Copy the number */
					numberI = 0;
					do {
						materialZ[numberI] = materialString[materialDAZIndex];
						materialDAZIndex++;
						numberI++;
					} while ((materialDAZIndex < materialStringLen-1) && (materialString[materialDAZIndex] != ' '));
					materialZ[numberI] = '\0';
					subObjMaterialDAZ[tissueIndex].Z = atof(materialZ);
				}
				else {
					subObjMaterialDAZ[tissueIndex].D = -1;
					subObjMaterialDAZ[tissueIndex].A = -1;
					subObjMaterialDAZ[tissueIndex].Z = -1;
				}
				
				/* Now get the numbers */
				for (attenIndex = 0; attenIndex < SUBOBJ_NUM_ENERGY_BINS; attenIndex++) {
					
					/* Get energy specific line */
					if (LbFlFGetS(inputBuffer, sizeof(inputBuffer), attenuationFile) == 0) {
						ErStGeneric("Error from LbFlFGetS while reading tissue values in tissue attenuation file.");
						goto FAIL;
					}
					
					/* Parse the string */
					numColumns = sscanf(inputBuffer, " %lf %lf %lf %lf", &(SubObjTissueAttenTableCoh[tissueIndex][attenIndex].energy),
								&(SubObjTissueAttenTableCoh[tissueIndex][attenIndex].attenuation),
								&(SubObjTissueAttenTableCoh[tissueIndex][attenIndex].probScatter),
								&(SubObjTissueAttenTableCoh[tissueIndex][attenIndex].probComptonToScatter));
								
					/* See if this is an old style file (just compton probability) */
					if (numColumns != 4){		
						sprintf(errString,"Error from sscanf while reading attenuation table '%s',\n got unexpected number of columns %ld from tissue #%ld on energy bin #%ld (SubObjCreate).\nYou are probably using a pre 2.5 attenuation table",
							PhgRunTimeParams.PhgSubObjAttenTableFilePath, (unsigned long)numColumns, 
							(unsigned long)tissueIndex, (unsigned long)attenIndex); 
						ErStGeneric(errString);
						goto FAIL;
					}
					
					/* Copy  coherent data into non coherent table */
					SubObjTissueAttenTableNoCoh[tissueIndex][attenIndex].energy = SubObjTissueAttenTableCoh[tissueIndex][attenIndex].energy;
					SubObjTissueAttenTableNoCoh[tissueIndex][attenIndex].attenuation = SubObjTissueAttenTableCoh[tissueIndex][attenIndex].attenuation;
					SubObjTissueAttenTableNoCoh[tissueIndex][attenIndex].probScatter = SubObjTissueAttenTableCoh[tissueIndex][attenIndex].probScatter;
					SubObjTissueAttenTableNoCoh[tissueIndex][attenIndex].probComptonToScatter = SubObjTissueAttenTableCoh[tissueIndex][attenIndex].probComptonToScatter;

					/* If we are not doing coherent, modify the table entries */
					if ((PHG_IsModelCoherentInObj() == false) || (PHG_IsModelCoherentInTomo() == false)) {
						double probPhotoelectric, probCompton;
						
						probPhotoelectric = 1 - SubObjTissueAttenTableNoCoh[tissueIndex][attenIndex].probScatter;
						probCompton = SubObjTissueAttenTableNoCoh[tissueIndex][attenIndex].probScatter 
							* SubObjTissueAttenTableNoCoh[tissueIndex][attenIndex].probComptonToScatter;
						SubObjTissueAttenTableNoCoh[tissueIndex][attenIndex].attenuation =
								(probPhotoelectric + probCompton) *
						 		SubObjTissueAttenTableNoCoh[tissueIndex][attenIndex].attenuation;

						SubObjTissueAttenTableNoCoh[tissueIndex][attenIndex].probScatter = probCompton / (probCompton + probPhotoelectric);
						SubObjTissueAttenTableNoCoh[tissueIndex][attenIndex].probComptonToScatter = 1.0;
 					}
				}
			}
		}

		/* Create the attenuation translation table */
		{			
			/* Open the file */
			if ((attenuationTransFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjAttIndexTransFilePath, "r")) == 0) {
				sprintf(errString, "Unable to open attenuation translation file '%s'.",
					PhgRunTimeParams.PhgSubObjAttIndexTransFilePath);
				ErStFileError(errString);
				goto FAIL;
			}
		
			/* Allocate space for the table */
			if ((attenuationTransTbl = (LbUsFourByte *) LbMmAlloc(
					SUBOBJ_NUM_TISSUE_TRANSLATIONS * sizeof(LbUsFourByte))) == 0){
			
				break;
			}
			
			/* Clear the table */
			for (tissueIndex = 0; tissueIndex < SUBOBJ_NUM_TISSUE_TRANSLATIONS; tissueIndex++)
				attenuationTransTbl[tissueIndex] = 0;
				
			/* Read in the file header */
			do {
				
				/* Get header lines */
				if (LbFlFGetS(inputBuffer, sizeof(inputBuffer), attenuationTransFile) == 0) {
					ErStFileError("Error while reading header from attenuation index translation file.");
					goto FAIL;
				}
			} while (inputBuffer[0] == '#');
			
			/* Read in lines and convert themm (first line already read from above, its the first non-comment line */
			loopIndex = 0;
			
			do {
				/* Reset parsing string */
				attInTranI = 0;
				
				/* Skip over white space */
				while (!isdigit(inputBuffer[attInTranI])) {
				
					/* See if we hit the end of line, probably a blank one */
					if (inputBuffer[attInTranI] == '\0')
						goto NEXT_LINE;
					attInTranI++;
				}
				/* Convert first item to integer */
				tissueIndex = atol(inputBuffer+attInTranI);
				
				/* Skip over non-white space */
				while (isdigit(inputBuffer[attInTranI])) {
				
					/* See if we hit the end of line, probably a blank one */
					if (inputBuffer[attInTranI] == '\0')
						goto NEXT_LINE;

					attInTranI++;
				}
				
				/* Skip over white space */
				while (!isdigit(inputBuffer[attInTranI])) {
				
					/* See if we hit the end of line, probably a blank one */
					if (inputBuffer[attInTranI] == '\0')
						goto NEXT_LINE;

					attInTranI++;
				}
				
				/* Convert first item to integer */
				tissueIndexTranslation = atol(inputBuffer+attInTranI);
				/* Verify that value is not too large */
				if (tissueIndex > SUBOBJ_NUM_TISSUE_TRANSLATIONS) {
					sprintf(errString, "Attenuation index out of range, value is %ld.\n"
						"This occurred on translation pair #%ld of the file named %s.\n"
						"The out of range value is the number in the first column.\n",
						(unsigned long)tissueIndex, (unsigned long)loopIndex+1, 
						PhgRunTimeParams.PhgSubObjAttIndexTransFilePath);
						
					ErStGeneric(errString);
					goto FAIL;
				}
				
				/* Verify that neither value is too large */
				if (tissueIndexTranslation > SubObjNumTissues) {
					sprintf(errString, "Attenuation translation index out of range, value is %ld.\n"
						"This occurred on translation pair #%ld of the file named %s.\n"
						"The out of range value is the number in the second column.\n",
						(unsigned long)tissueIndexTranslation, (unsigned long)loopIndex+1, 
						PhgRunTimeParams.PhgSubObjAttIndexTransFilePath);
					ErStGeneric("Translation tissue index out of range in attenuation index translation table.\n");
					goto FAIL;
				}
				
				/* Store translation in tranlsation table, this affectively performs the translation */
				attenuationTransTbl[tissueIndex] = tissueIndexTranslation;
				
				/* Get next pair of translation values */
				NEXT_LINE:;
				if (tissueIndex < SUBOBJ_NUM_TISSUE_TRANSLATIONS-1) {
					if (LbFlFGetS(inputBuffer, sizeof(inputBuffer), attenuationTransFile) == 0) {
						ErStFileError("Error while reading attenuation index translation value.");
						goto FAIL;
					}
				}
				
				loopIndex++;
			} while (loopIndex <= SUBOBJ_NUM_TISSUE_TRANSLATIONS);
			
			/* Now run through attenuation indexes and perform the translation */
			for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
				for (yIndex = 0; yIndex < SubObjObject[sliceIndex].attNumYBins; yIndex++) {
					for (xIndex = 0; xIndex < SubObjObject[sliceIndex].attNumXBins; xIndex++) {
																	
						/* Translate the voxel index */
						tissueIndex = SubObjObject[sliceIndex].attenuationArray[
							(yIndex * SubObjObject[sliceIndex].attNumXBins)
							+ xIndex];
							
						/* Verify that value is not too great */
						if (tissueIndex > SUBOBJ_NUM_TISSUE_TRANSLATIONS) {
							sprintf(errString, "Tissue index out of range in attenuation object.\n"
								"The tissue value is %ld, the maximum allowed value is %d",
								(unsigned long)tissueIndex, SUBOBJ_NUM_TISSUE_TRANSLATIONS);
							ErStGeneric(errString);
							goto FAIL;
						}
						
						/* Stuff the translation value in the attenuation array */
						SubObjObject[sliceIndex].attenuationArray[
							(yIndex * SubObjObject[sliceIndex].attNumXBins)
							+ xIndex] = 
							attenuationTransTbl[tissueIndex];				
					}
				}
			}
			
			/* Dump the attenuation image if they requested it. */
			if (strlen(PhgRunTimeParams.PhgSubObjAttImgFilePath) != 0) {
				
				/* Attempt to open the file for re-write */
				if ((attImgFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjAttImgFilePath, "wb")) == 0){
					sprintf(errString, "\nUnable to create attenuation image file named '%s'\n",
						PhgRunTimeParams.PhgSubObjAttImgFilePath);
						
					ErAlert(errString, false);
				}
				else {
					/* Dump the image file */
					do {
						for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
							for (yIndex = 0; yIndex < SubObjObject[sliceIndex].attNumYBins; yIndex++) {
								for (xIndex = 0; xIndex < SubObjObject[sliceIndex].attNumXBins; xIndex++) {
									
									/* Get attenutation of current cell, fails if we are out of the object */
									if (!SubObjGetCellAttenuation(sliceIndex,
											xIndex, yIndex,
											PhgRunTimeParams.PhgNuclide.photonEnergy_KEV,
											&voxelAttenuation)) {
										PhgAbort("Attempt to get attenuation for cell not in subobject (SubObjCreate).",
											true);
									}
									
									/* Conver to float */
									voxelAttAsFloat = (float) voxelAttenuation;
									if (fwrite(&voxelAttAsFloat, sizeof(float), 1, attImgFile) != 1){
										
										ErAlert("\nUnable to write to attenuation image file.\n", false);
										goto FAIL_WRITE_ATT;				
									}
								}
							}
						}
						FAIL_WRITE_ATT:;
					}
					while (false);
					
					/* Close the file */
					fclose(attImgFile);
					attImgFile = 0;
				}
				
			}
			/* Close the file */
			fclose(attenuationTransFile);
		
			/* Free the memory */
			LbMmFree((void **)&attenuationTransTbl);
			
		}


		/* Create the activity translation table */
		{			
			/* Open the file */
			if ((activityTransFile = 
					LbFlFileOpen(PhgRunTimeParams.PhgSubObjActIndexTransFilePath, "r"))
					== 0) {
					
				sprintf(errString, "Unable to open activity translation file '%s'.",
					PhgRunTimeParams.PhgSubObjActIndexTransFilePath);
				ErStFileError(errString);
				goto FAIL;
			}
		
			/* Allocate space for the table */
			if ((activityTransTbl = (LbUsFourByte *) LbMmAlloc(
					SUBOBJ_NUM_TISSUE_TRANSLATIONS * sizeof(LbUsFourByte))) == 0){
			
				break;
			}
			
			/* Clear the table */
			for (tissueIndex = 0; tissueIndex < SUBOBJ_NUM_TISSUE_TRANSLATIONS; tissueIndex++)
				activityTransTbl[tissueIndex] = 0;
				
			/* Read in the file header */
			do {
				
				/* Get header lines */
				if (LbFlFGetS(inputBuffer, sizeof(inputBuffer), activityTransFile) == 0) {
					ErStFileError("Error while reading header from activity index translation file.");
					goto FAIL;
				}
			} while (inputBuffer[0] == '#');
			
			/* Read in lines and convert them */
			loopIndex = 0;
			do {
				/* Convert line to index values */
				if (sscanf(inputBuffer, " %ld %ld", 
						&tissueIndexUSL, &tissueIndexTranslationUSL) != 2){
						
						ErStGeneric("Error from sscanf on activity index translation value.\n");
						goto FAIL;
				}
				tissueIndex = (LbUsFourByte)tissueIndexUSL;
				tissueIndexTranslation = (LbUsFourByte)tissueIndexTranslationUSL;
				
				/* Verify index will not be too large */
				if (tissueIndex >= SUBOBJ_NUM_TISSUE_TRANSLATIONS) {
						
					ErStFileError("Invalid index in left hand column of activity translation table.");
					goto FAIL;
				}
				
				/* Store translation */
				activityTransTbl[tissueIndex] = tissueIndexTranslation;
				
				/* Get tissue name */
				if (tissueIndex < SUBOBJ_NUM_TISSUE_TRANSLATIONS-1) {
					if (LbFlFGetS(inputBuffer, sizeof(inputBuffer), activityTransFile) == 0) {
						ErStFileError("Error while reading activity index translation value.");
						goto FAIL;
					}
				}
				
				loopIndex++;
			} while (loopIndex <= SUBOBJ_NUM_TISSUE_TRANSLATIONS);
			
			/* Initialize global that keeps track of total activity * time (* decays/Ci/sec) */
			SubObjTotalRealDecays = 0;
			
			/* Now run through activity indexes and perform the translation */
			for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
				for (yIndex = 0; yIndex < SubObjObject[sliceIndex].actNumYBins; yIndex++) {
					for (xIndex = 0; xIndex < SubObjObject[sliceIndex].actNumXBins; xIndex++) {
						
						/* Verify index will not be too large */
						if (SubObjObject[sliceIndex].activityArray[
								(yIndex * SubObjObject[sliceIndex].actNumXBins)
								+ xIndex] > SUBOBJ_NUM_TISSUE_TRANSLATIONS) {
								
							ErStFileError("Invalid index in activity index file.");
							goto FAIL;
						}
						
						/* Translate the voxel index */
						SubObjObject[sliceIndex].activityArray[
							(yIndex * SubObjObject[sliceIndex].actNumXBins)
							+ xIndex] = 
							activityTransTbl[SubObjObject[sliceIndex].activityArray[
							(yIndex * SubObjObject[sliceIndex].actNumXBins)
							+ xIndex]];				
					}
				}
				
				/* Compute total slice activity */
				subObjCalcSliceActivity(sliceIndex);
				
				/* keep track of total activity */
				SubObjTotalRealDecays += SubObjObject[sliceIndex].sliceActivity;
				
			}
			
			/* complete calcuation of SubObjTotalRealDecays */
			SubObjTotalRealDecays *= SUBOBJ_DECAYS_PER_CURIE;
			
			/* Compute the total number of decays expected in a real-world scan */
			
		
			/* Dump the activity image if they requested it. */
			if (strlen(PhgRunTimeParams.PhgSubObjActImgFilePath) != 0) {
				
				/* Attempt to open the file for re-write */
				if ((actImgFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjActImgFilePath, "wb")) == 0){
					sprintf(errString, "\nUnable to create activity image file named '%s'\n",
						PhgRunTimeParams.PhgSubObjActImgFilePath);
						
					ErAlert(errString, false);
				}
				else {
					/* Dump the image file */
					do {
						for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
							for (yIndex = 0; yIndex < SubObjObject[sliceIndex].actNumYBins; yIndex++) {
								for (xIndex = 0; xIndex < SubObjObject[sliceIndex].actNumXBins; xIndex++) {
									
									/* Retreive tissue index for voxel */
									tissueIndex = SubObjGetTissueIndex(sliceIndex, yIndex, xIndex);
						
									/* Look up activity and add to running total */
									voxelActivity = (float) SubObjGetTissueActivity(tissueIndex, 0);
									if (fwrite(&voxelActivity, sizeof(float), 1, actImgFile) != 1){
										
										ErAlert("\nUnable to write to activity image file.\n",
											false);
										goto FAIL_WRITE_ACT;				
									}
								}
							}
						}
						FAIL_WRITE_ACT:;
					}
					while (false);
					
					/* Close the file */
					fclose(actImgFile);
					actImgFile = 0;
				}
				
			}

			/* Close the file */
			fclose(activityTransFile);
			
			/* Free the memory */
			LbMmFree((void **)&activityTransTbl);
		}

		okay = true;
		FAIL:;
	} while (false);
	
	/* Free up any memory we can upon failure */
	if (!okay) {
				
		if (SubObjTissueTable.tissueValues != 0) {
			for (tissueIndex = 0; tissueIndex < SubObjNumActIndexes; tissueIndex++) {
				if (SubObjTissueTable.tissueValues[tissueIndex].activityValues != 0) {
					LbMmFree((void **) &SubObjTissueTable.tissueValues[tissueIndex].activityValues);
				}
			}
			LbMmFree((void **) &(SubObjTissueTable.tissueValues));
		}

		if (SubObjTissueAttenTableNoCoh != 0) {
			LbMmFree((void **) &SubObjTissueAttenTableNoCoh);
		}

		if (SubObjTissueAttenTableCoh != 0) {
			LbMmFree((void **) &SubObjTissueAttenTableCoh);
		}

		if (attenuationTransTbl != 0) {
			LbMmFree((void **) &attenuationTransTbl);
		}

		if (activityTransTbl != 0) {
			LbMmFree((void **) &attenuationTransTbl);
		}

		for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
		
			if (SubObjObject[sliceIndex].activityArray != 0) {
				LbMmFree((void **) &(SubObjObject[sliceIndex].activityArray));
			}
		
			if (SubObjObject[sliceIndex].attenuationArray != 0) {
				LbMmFree((void **) &(SubObjObject[sliceIndex].attenuationArray));
			}
		}
	}
	
	return (okay);
}

/*********************************************************************************
*
*			Name:		subObjCreateObject
*
*			Summary:	Create the object from the input file.
*
*			Arguments:
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean	subObjCreateObject()
{
	Boolean					okay = false;			/* Process Flag */
	char					errString[256];			/* Space for building error strings */
	FILE					*actIndexFile = 0;		/* File containing voxel indexes */
	FILE					*attIndexFile = 0;		/* File containing voxel indexes */
	LbUsFourByte			voxelIndex;				/* Index for the voxel */
	LbUsTwoByte				sliceIndex;				/* Current slice */
	LbUsTwoByte				xIndex;					/* Current x index */
	LbUsTwoByte				yIndex;					/* Current y index */
	
	do { /* Process Loop */

		/* Open the activity index file */
		if ((actIndexFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjActivityIndexFilePath, "rb")) == 0) {
			sprintf(errString, "Unable to open activity index file '%s'.", 
				PhgRunTimeParams.PhgSubObjActivityIndexFilePath);
			ErStFileError(errString);
			goto FAIL;
		}

		/* Open the attenuation index file */
		if ((attIndexFile = LbFlFileOpen(PhgRunTimeParams.PhgSubObjAttenIndexFilePath, "rb")) == 0) {
			sprintf(errString, "Unable to open attenuation index file '%s'.",
				PhgRunTimeParams.PhgSubObjAttenIndexFilePath);
			ErStFileError(errString);
			goto FAIL;
		}
		
		/* Calculate the computed values from those initialized by the parameter file */
		for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {


			/* Compute slice width */
			SubObjObject[sliceIndex].sliceWidth = 
				SubObjObject[sliceIndex].xMax - 
				SubObjObject[sliceIndex].xMin;
				
			/* Verify that slice has positive width */
			if (SubObjObject[sliceIndex].sliceWidth <= 0){;
				sprintf(errString, "\nSlice %d has invalid xMin/xMax.",
					sliceIndex);
					
				ErStGeneric(errString);
				goto FAIL;
			}
			
			/* Compute slice height */
			SubObjObject[sliceIndex].sliceHeight = 
				SubObjObject[sliceIndex].yMax - 
				SubObjObject[sliceIndex].yMin;
			
				
			/* Verify that slice has positive height */
			if (SubObjObject[sliceIndex].sliceHeight <= 0){;
				sprintf(errString, "\nSlice %d has invalid yMin/yMax.",
					sliceIndex);
					
				ErStGeneric(errString);
				goto FAIL;
			}

			/* Compute slice depth */
			SubObjObject[sliceIndex].sliceDepth = fabs(
				SubObjObject[sliceIndex].zMax - 
				SubObjObject[sliceIndex].zMin);
				
			/* Verify that slice has positive depth */
			if (SubObjObject[sliceIndex].sliceDepth <= 0){;
				sprintf(errString, "\nSlice %d has invalid zMin/zMax.",
					sliceIndex);
					
				ErStGeneric(errString);
				goto FAIL;
			}

			SubObjObject[sliceIndex].actVoxelWidth = SubObjObject[sliceIndex].sliceWidth/
				SubObjObject[sliceIndex].actNumXBins;

			SubObjObject[sliceIndex].attVoxelWidth = SubObjObject[sliceIndex].sliceWidth/
				SubObjObject[sliceIndex].attNumXBins;

			SubObjObject[sliceIndex].actVoxelHeight = SubObjObject[sliceIndex].sliceHeight/
				SubObjObject[sliceIndex].actNumYBins;

			SubObjObject[sliceIndex].attVoxelHeight = SubObjObject[sliceIndex].sliceHeight/
				SubObjObject[sliceIndex].attNumYBins;

			/* Allocate slices for activity slice */
			if ((SubObjObject[sliceIndex].activityArray = (SubObjActVoxelTy *) LbMmAlloc(
					(sizeof(SubObjActVoxelTy) * (SubObjObject[sliceIndex].actNumXBins * 
					SubObjObject[sliceIndex].actNumYBins)))) == 0) {
				
				goto FAIL;
			}
	
			/* Allocate slices for attenuation slice */
			if ((SubObjObject[sliceIndex].attenuationArray = (LbUsFourByte *) LbMmAlloc(
					(sizeof(LbUsFourByte) * (SubObjObject[sliceIndex].attNumXBins * 
					SubObjObject[sliceIndex].attNumYBins)))) == 0) {
				
				goto FAIL;
			}
			
			/* Read in activity index data */
			for (yIndex = 0; yIndex < SubObjObject[sliceIndex].actNumYBins; yIndex++) {
				
				for (xIndex = 0; xIndex < SubObjObject[sliceIndex].actNumXBins; xIndex++) {
				
					/* Read from the activity index file */
					if (fread(&voxelIndex, sizeof(LbUsFourByte), 1, actIndexFile) != 1) {
						sprintf(errString, "Failed trying to read slice - voxel activity index (%d, (%d, %d)), from file %s.\n",
							sliceIndex, xIndex, yIndex, PhgRunTimeParams.PhgSubObjActivityIndexFilePath);
						ErStFileError(errString);
						goto FAIL;
					}
					
					/* Save the voxel index */
					SubObjObject[sliceIndex].activityArray[
					    (yIndex * SubObjObject[sliceIndex].actNumXBins)
						+ xIndex] = voxelIndex;
				}
			}

			/* Read  in attenuation index data */
			for (yIndex = 0; yIndex < SubObjObject[sliceIndex].attNumYBins; yIndex++) {
				
				for (xIndex = 0; xIndex < SubObjObject[sliceIndex].attNumXBins; xIndex++) {
										
					/* Read from the attenuation index file */
					if (fread(&voxelIndex, sizeof(LbUsFourByte), 1, attIndexFile) != 1) {
						sprintf(errString, "Failed trying to read slice - voxel attenuation index (%d, (%d, %d)), from file %s.\n",
							sliceIndex, xIndex, yIndex, PhgRunTimeParams.PhgSubObjAttenIndexFilePath);
						ErStFileError(errString);
						goto FAIL;
					}
					
					/* Save the voxel index */
					SubObjObject[sliceIndex].attenuationArray[
						(yIndex * SubObjObject[sliceIndex].attNumXBins)
						+ xIndex] = voxelIndex;				
				}
			}
		 }

		/* Verify that the slices are contiguous and of the same x/y dimension */
		for (sliceIndex = 1; sliceIndex < (SubObjNumSlices); sliceIndex++) {
	
			/* Verify xMax of current slice is the same as previous slice */
			if (SubObjObject[sliceIndex-1].xMax != SubObjObject[sliceIndex].xMax){
				ErStGeneric("Slices must have same x/y dimensions!");
				goto FAIL;
			}
			/* Verify xMin of current slice is the same as previous slice */
			if (SubObjObject[sliceIndex-1].xMin != SubObjObject[sliceIndex].xMin){
				ErStGeneric("Slices must have same x/y dimensions!");
				goto FAIL;
			}
			/* Verify yMax of current slice is the same as previous slice */
			if (SubObjObject[sliceIndex-1].yMax != SubObjObject[sliceIndex].yMax){
				ErStGeneric("Slices must have same x/y dimensions!");
				goto FAIL;
			}
			/* Verify yMin of current slice is the same as previous slice */
			if (SubObjObject[sliceIndex-1].yMin != SubObjObject[sliceIndex].yMin){
				ErStGeneric("Slices must have same x/y dimensions!");
				goto FAIL;
			}
				
			/* Check that slice-1 == slice */
			if (SubObjObject[sliceIndex-1].zMax != SubObjObject[sliceIndex].zMin){;
				sprintf(errString, "\nSlice %d and %d are not contiguous in the axial direction.",
					sliceIndex-1, sliceIndex);
					
				ErStGeneric(errString);
				goto FAIL;
			}
		 }

		/* Verify xMax of current slice is the same as yMax */
		if (SubObjObject[0].xMax != SubObjObject[0].yMax){
			sprintf(errString,"\nCurrently, the object cylinder must be round.\n\t"
				"In slice %d, you have specified x-max as %3.2f and y-max as %3.2f\n\t"
				"These values must be equal.\n", 0, SubObjObject[0].xMax,
				SubObjObject[0].yMax);
			ErStGeneric(errString);
			goto FAIL;
		}

		/* Verify xMin of current slice is the same as yMin */
		if (SubObjObject[0].xMin != SubObjObject[0].yMin){
			sprintf(errString,"\nCurrently, the object cylinder must be round.\n\t"
				"In slice %d, you have specified x-min as %3.2f and y-min as %3.2f\n\t"
				"These values must be equal.\n", 0, SubObjObject[0].xMin,
				SubObjObject[0].yMin);
			ErStGeneric(errString);
			goto FAIL;
		}

		/* Verify object cylinder is round */
		if ((SubObjObject[0].xMin != -SubObjObject[0].xMax)){
			sprintf(errString,"\nCurrently, the object cylinder must be round.\n\t"
				"In slice %d, you have specified x-min as %3.2f and x-max as %3.2f\n\t"
				"These values must be symmetric.\n", 0, SubObjObject[0].xMin,
				SubObjObject[0].xMax);
			ErStGeneric(errString);
			goto FAIL;
		}

		/* Verify object cylinder is centered */
		if ((SubObjObject[0].xMin + SubObjObject[0].yMax) != 0){
			sprintf(errString,"\nCurrently, the object cylinder must be centered.\n\t"
				"In slice %d, you have specified x-min as %3.2f and x-max as %3.2f\n\t"
				"These values must be symmetric.\n", 0, SubObjObject[0].xMin,
				SubObjObject[0].xMax);
			ErStGeneric(errString);
			goto FAIL;
		}

		okay = true;
		FAIL:;
	} while (false);
	
	/* Close the files */
	if (actIndexFile != 0)
		fclose(actIndexFile);

	if (attIndexFile != 0)
		fclose(attIndexFile);
	
	/* If we failed free any memory that was allocated */
	if (!okay) {
		for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
			if (SubObjObject[sliceIndex].activityArray != 0) {
				LbMmFree((void **) &(SubObjObject[sliceIndex].activityArray));
			}
			/* Allocate slices for activity slice */
			if (SubObjObject[sliceIndex].attenuationArray != 0) {
				LbMmFree((void **) &(SubObjObject[sliceIndex].attenuationArray));
			}
		}
	}
	
	return (okay);
}

/*********************************************************************************
*
*			Name:		SubObjGtPositionIndexes
*
*			Summary:	Given a position, get the indexes into the attenuation
*						object.
*
*			Arguments:
*				PHG_Position	*posPtr 					- Storage for the new decay.
*				LbUsFourByte	*sliceIndexPtr				- Storage for slice the decay came from.
*				LbUsFourByte	*xIndexPtr					- Storage for x index into object.
*				LbUsFourByte	*yIndexPtr					- Storage for y index into object.
*	
*			Function return: None.
*
*********************************************************************************/
void	SubObjGtPositionIndexes(PHG_Position *posPtr,
			LbFourByte *sliceIndexPtr,
			LbFourByte *xIndexPtr,
			LbFourByte *yIndexPtr)
{
	
	/* Find which slice we are in */
	*sliceIndexPtr = 0;
	do {
	
		/* Break if we are in this slice, note that this loop should
			always terminate via this condition
		*/
		if (posPtr->z_position <= SubObjObject[*sliceIndexPtr].zMax)
			break;
			
		(*sliceIndexPtr)++;
	} while (*sliceIndexPtr < (LbFourByte)SubObjNumSlices);
	
	#ifdef PHG_DEBUG
		/* This should never happen */
		if (*sliceIndexPtr >= (LbFourByte)SubObjNumSlices)
			PhgAbort("Invalid calculation of current slice (SubObjGtPositionIndexes)", false);
	#endif
		
	/* Find which x bin we are in */
	*xIndexPtr = 0;
	do {
	
		/* Break if we are in this X voxel, note that this loop should
			always terminate via this condition
		*/
	    if (posPtr->x_position <= (SubObjObject[*sliceIndexPtr].xMin +
				(((*xIndexPtr)+1) * SubObjObject[*sliceIndexPtr].attVoxelWidth)))
			break;
			
		(*xIndexPtr)++;
	} while (*xIndexPtr < (LbFourByte)(SubObjObject[*sliceIndexPtr].attNumXBins));
		
	/* Find which y bin we are in */
	*yIndexPtr = 0;
	do {
	
		/* Break if we are in this Y voxel, note that this loop should
			always terminate via this condition
		*/
		if (posPtr->y_position >= (SubObjObject[*sliceIndexPtr].yMax -
				(((*yIndexPtr)+1) * SubObjObject[*sliceIndexPtr].attVoxelHeight)))
			break;
			
		(*yIndexPtr)++;
	} while (*yIndexPtr < (LbFourByte)(SubObjObject[*sliceIndexPtr].attNumYBins));
	
}


/*********************************************************************************
*
*			Name:		subObjPrintStatus
*
*			Summary:	Prints CPU time and percentage complete status.
*
*			Arguments:
*	
*			Function return: none.
*
*********************************************************************************/
void subObjPrintStatus()
{
	char				timeStr[32];				/* String for printing the time */
	time_t				curTime;					/* The current wall clock time */
	double				perDone;					/* Percentage completed */
	Boolean				timingValid;				/* Whether timing exists on this system */
	double				trackPhotonsTime;			/* User time to track photons */
	double				trackPhotonsCPUTime;		/* System time to track photons */
	
	
	/* Get timing info */
	timingValid = LbTmStopTiming(&SubObjTrackPhotonsStartTime, 
									&trackPhotonsTime, &trackPhotonsCPUTime);
	
	/* Get current time */
	time(&curTime);
	
	/* Convert time to string format */
	strftime(timeStr, 31, "(%I:%M %p - %m/%d/%y)", localtime(&curTime));
	
	/* Compute percentage completed, two steps eliminates overflow */
	perDone = (double)SUBOBJGetDecaysProcessed()/(double)PhgRunTimeParams.Phg_EventsToSimulate;
	
	perDone *= 100;
	
	if (timingValid) {
		LbInPrintf(" %3.0f%% tracked,\tCPU Time = %7.0f seconds\t%s.\n",
			perDone, trackPhotonsCPUTime, timeStr);
	}
	else {
		LbInPrintf(" %3.0f%% tracked,\t%s.\n",
			perDone, timeStr);
	}
}

/*********************************************************************************
*
*			Name:		SubObjGenVoxAngCellDecay
*
*			Summary:	Generate a new decay.
*
*			Arguments:
*				PHG_Decay		*newDecayPtr 				- Storage for the new decay.
*				PHG_Direction	*newDecayEmissionAnglePtr	- Storage for the new decay's emmision angle.
*				LbUsFourByte	*sliceIndexPtr				- Storage for slice the decay came from.
*				LbUsFourByte	*angleIndexPtr				- Storage for angle the decay came from.
*				LbUsFourByte	*xIndexPtr					- Storage for x index into object.
*				LbUsFourByte	*yIndexPtr					- Storage for y index into object.
*	
*			Function return: true if there is a decay to generate.
*
*********************************************************************************/
Boolean	SubObjGenVoxAngCellDecay(PHG_Decay *newDecayPtr,
			PHG_Direction *newDecayEmissionAnglePtr, LbFourByte *sliceIndexPtr,
			LbFourByte *angleIndexPtr,  LbFourByte *xIndexPtr,
			LbFourByte *yIndexPtr)
{
	Boolean 		gotOne = false;			/* Storage for whether we got one or not */
	Boolean			validLocation = false;	/* Is our location valid */
	double			cosAlpha;				/* Cosine of alpha angle */
	double			dTemp;					/* Temporary value called d */
	double			theta;					/* Theta angle */
	double			rand1;					/* First random number */
	double			rand2;					/* Second random number */
	LbUsFourByte	rowNumber;				/* The current row */
	LbFourByte		dummySlice;				/* Used in getting attenuation x/y indexes */
	
	/* Loop until valid location or run out of decays */
	do {
    
	    #ifdef PHG_DEBUG
	    	/* Write random seed value if requested, always resetting to beggining */
	    	if (subObjDoWriteRandSeed == true)
	    		PhgMathWriteSeed(1);
	    		
	    	/* Read random seed value if requested */
	    	if (subObjDoReadRandSeed == true)
	    		PhgMathReadSeed(1);
	    		
	    #endif
		
		/* Setup next voxel decay */
		if ((gotOne = subObjGetNextVoxAngCellDecay()) == false) {
		
			/* Note, to get a screen output of 100% always test to see if we need to print a message */
			if (SubObjDecaysProcessed < PhgRunTimeParams.Phg_EventsToSimulate) {
				SubObjDecaysProcessed = PhgRunTimeParams.Phg_EventsToSimulate;
				subObjPrintStatus();
			}
			
			/* Since there are no more decays, bust out of here */	
			break;
		}
		
		/* Increment counter */
		SubObjDecaysProcessed++;

		/* Print status report if we're at a 10% increment */
		if ((SubObjDecaysProcessed % (PhgRunTimeParams.Phg_EventsToSimulate/10)) == 0)
			subObjPrintStatus();
			
			
		/* Save the starting location */
		*sliceIndexPtr = SubObjCurSliceIndex;
		*angleIndexPtr = SubObjCurAngleIndex;

		/* Calculate the current row we are on (starting from top of object going down, 0 to n) */
		rowNumber = (SubObjCurVoxelIndex/
					 SubObjObject[SubObjCurSliceIndex].actNumXBins);
		
		/* Save the x/y indexes, note these are currently computed from the activity object,
		   at the end of this routine we'll set them according to the attenuation object
		   for tracking
		 */
		*xIndexPtr = SubObjCurVoxelIndex % SubObjObject[SubObjCurSliceIndex].actNumXBins;
		*yIndexPtr = rowNumber;
		
		/* If we have point source voxels, center decay in voxel */
		if (PHG_IsVoxelPointSource()) {
	
			/* Choose the X position within the voxel */
			newDecayPtr->location.x_position = (SubObjObject[SubObjCurSliceIndex].xMin +
				(SubObjObject[SubObjCurSliceIndex].actVoxelWidth *
				*xIndexPtr) +
				(SubObjObject[SubObjCurSliceIndex].actVoxelWidth/2));
			
	
			/* Choose the Y position within the voxel, noting y indexes go opposite to y position */
			newDecayPtr->location.y_position = (SubObjObject[SubObjCurSliceIndex].yMax -
				(SubObjObject[SubObjCurSliceIndex].actVoxelHeight * *yIndexPtr) -
				(SubObjObject[SubObjCurSliceIndex].actVoxelHeight/2));
	
			/* Choose the Z position within the voxel */
			newDecayPtr->location.z_position = SubObjObject[SubObjCurSliceIndex].zMin +
				(SubObjObject[SubObjCurSliceIndex].sliceDepth/2);
		}
		else if (PHG_IsVoxelLineSource()) {
	
			/* Choose the X position within the voxel */
			newDecayPtr->location.x_position = (SubObjObject[SubObjCurSliceIndex].xMin +
				(SubObjObject[SubObjCurSliceIndex].actVoxelWidth *
				*xIndexPtr) +
				(SubObjObject[SubObjCurSliceIndex].actVoxelWidth/2));
			
	
			/* Choose the Y position within the voxel, noting y indexes go opposite to y position */
			newDecayPtr->location.y_position = (SubObjObject[SubObjCurSliceIndex].yMax -
				(SubObjObject[SubObjCurSliceIndex].actVoxelHeight * *yIndexPtr) -
				(SubObjObject[SubObjCurSliceIndex].actVoxelHeight/2));
	
			/* Choose the Z position within the voxel */
			newDecayPtr->location.z_position = SubObjObject[SubObjCurSliceIndex].zMin +
				(PhgMathGetRandomNumber() * SubObjObject[SubObjCurSliceIndex].sliceDepth);
		}
		else {			
			/*	The location within the voxel is chosen randomly with uniform 
				distribution over the voxel area.
			*/
	
			/* Choose the X position within the voxel */
			newDecayPtr->location.x_position = (SubObjObject[SubObjCurSliceIndex].xMin +
				(SubObjObject[SubObjCurSliceIndex].actVoxelWidth *
				*xIndexPtr) +
				(SubObjObject[SubObjCurSliceIndex].actVoxelWidth * PhgMathGetRandomNumber()));
			
	
			/* Choose the Y position within the voxel */
			newDecayPtr->location.y_position = (SubObjObject[SubObjCurSliceIndex].yMax -
				(SubObjObject[SubObjCurSliceIndex].actVoxelHeight * *yIndexPtr) -
				(SubObjObject[SubObjCurSliceIndex].actVoxelHeight * PhgMathGetRandomNumber()));
			
			/* Choose the Z position within the voxel */
			newDecayPtr->location.z_position = SubObjObject[SubObjCurSliceIndex].zMin +
				(SubObjObject[SubObjCurSliceIndex].sliceDepth * PhgMathGetRandomNumber());
	
		}
		
		/* IF location is outside object, start this loop again */
		validLocation = (CylPosIsOutsideObjCylinder(newDecayPtr->location) == false);
		if (validLocation == false) {
		
			/* Increment count of decays outside object cylinder */
			PhoHStatIncTotInvalPhoLocations();
			
			/* Start loop over */
			continue;
		}

		/* Get the decay weight */
		newDecayPtr->startWeight =
			SubObjDecayWeightSlice[(SubObjCurVoxelIndex*PRODTBLGetNumAngleCells())+SubObjCurAngleIndex];
		
		/* Set the decay time randomly within the scan */
		newDecayPtr->decayTime = PhgMathGetDPRandomNumber() * 
							(double)(PHGGetLengthOfScan());
		
		/* Set the decay type; note that PhgEn_Complex is not yet supported, PhgEn_PETRandom is
		added in addrandom.c */
		if ( PHG_IsSPECT() ) {
			newDecayPtr->decayType = PhgEn_SinglePhoton;
		} else if ( PHG_IsPETCoincPlusSingles() || PHG_IsPETCoincidencesOnly() ) {
			newDecayPtr->decayType = PhgEn_Positron;
		} else {
			newDecayPtr->decayType = PhgEn_Unknown;
		}
						
		/*	The location within the angle bin of the axial emission angle is also
			chosen randomly with uniform distribution over the
			angle bin range.
		*/

		/* Get random numbers */
		rand1 = PhgMathGetRandomNumber();
		rand2 = PhgMathGetRandomNumber();

		/* Calculate theta */
		theta = PHGMATH_2PI * rand1;

		/* Get cosine alpha from random pick within current stratification cell */
		cosAlpha = PRODTBLGetProdTblAngleStart(SubObjCurSliceIndex, SubObjCurAngleIndex) +
			((PRODTBLGetProdTblAngleEnd(SubObjCurSliceIndex, SubObjCurAngleIndex) -
			PRODTBLGetProdTblAngleStart(SubObjCurSliceIndex, SubObjCurAngleIndex)) * rand2);

		/* Choose the emission axial angle within the angle bin */
		{
			/* Calculate temporary value d */
			dTemp = PHGMATH_SquareRoot(1 - PHGMATH_Square(cosAlpha));

			/* X cosine is just d(cos(theta)) */
			newDecayEmissionAnglePtr->cosine_x = dTemp * PHGMATH_Cosine(theta);
	
			/* Cosine y is just d(sine(theta)) */
			newDecayEmissionAnglePtr->cosine_y = dTemp * PHGMATH_Sine(theta);
		}

		/* Assign Z cosine as cosAlpha */
		newDecayEmissionAnglePtr->cosine_z = cosAlpha;	

        #ifdef SUBOBJ_DEBUG
		{
			char 	debugStr[150];
			double	absDiff;		/* Absolute difference */
	
			/* Make sure we are close to 1.0 */
			if (PhgMathRealNumAreEqual(PHGMATH_Square(newDecayEmissionAnglePtr->cosine_x) +
					PHGMATH_Square(newDecayEmissionAnglePtr->cosine_y) +
					PHGMATH_Square(newDecayEmissionAnglePtr->cosine_z), 1, -7, 0, &absDiff, 0) == false) {
	
				sprintf(debugStr, "Invalid calculation of direction vector (SubObjGenVoxAngCellDecay). \nSum of squares = %3.2e, 1.0 - Sum of squares = %3.2e",
					PHGMATH_Square(newDecayEmissionAnglePtr->cosine_x) +
					PHGMATH_Square(newDecayEmissionAnglePtr->cosine_y) +
					PHGMATH_Square(newDecayEmissionAnglePtr->cosine_z),
					absDiff);
			
				PhgAbort(debugStr, true);
		    }
        }
        #endif

		/* If they have requested a fixed direction vector, slam over the previous calculations */
		#ifdef PHG_DEBUG
			if (PHGDEBUG_FixedDirection()) {
			     
			   /* Set from #defines in SubObj.h  */ 
				newDecayEmissionAnglePtr->cosine_x = SUBOBJ_FIXED_DIR_COSINE_X;
				newDecayEmissionAnglePtr->cosine_y = SUBOBJ_FIXED_DIR_COSINE_Y;
				newDecayEmissionAnglePtr->cosine_z = SUBOBJ_FIXED_DIR_COSINE_Z;
			  
 
				/* 45 degree angle x/y 
				newDecayEmissionAnglePtr->cosine_x = 0.707;
				newDecayEmissionAnglePtr->cosine_y = 0.707;
				newDecayEmissionAnglePtr->cosine_z = 0.000001;
			    */
			    
				/* 45 degree angle x/z  
				newDecayEmissionAnglePtr->cosine_x = 0.707;
				newDecayEmissionAnglePtr->cosine_y = 0.000001;
				newDecayEmissionAnglePtr->cosine_z = 0.707;
			   */
				/* 30 degree angle y/x
				newDecayEmissionAnglePtr->cosine_x = 0.866;
				newDecayEmissionAnglePtr->cosine_y = 0.5;
				newDecayEmissionAnglePtr->cosine_z = 0.000001;
			 */

                                /* odd angle--fill in your own here
                                newDecayEmissionAnglePtr->cosine_x = 0.9018099;
                                newDecayEmissionAnglePtr->cosine_y = 0.4205208;
                                newDecayEmissionAnglePtr->cosine_z = 0.0995037;
                        */

				/* pencilBeam3 
				newDecayEmissionAnglePtr->cosine_x = 0.9805806;
				newDecayEmissionAnglePtr->cosine_y = 0.1961161;
				newDecayEmissionAnglePtr->cosine_z = 0.0;
			 */
			 
                                /* pencilBeam4 
                                newDecayEmissionAnglePtr->cosine_x = 0.8939764;
                                newDecayEmissionAnglePtr->cosine_y = 0.416868;
                                newDecayEmissionAnglePtr->cosine_z = -0.1643989;
			*/
                                /* pencilBeam5, 7, 8 
                                newDecayEmissionAnglePtr->cosine_x = 0.9018099;
                                newDecayEmissionAnglePtr->cosine_y = 0.4205208;
                                newDecayEmissionAnglePtr->cosine_z = 0.0995037;
			*/

			}
		#endif	
	} while (!validLocation);
	
	
	/* Now set x and y index according to attenuation object voxelization
	(if we got a decay) */
	if (gotOne) {
		SubObjGtPositionIndexes(&newDecayPtr->location,
			&dummySlice,
			xIndexPtr,
			yIndexPtr);
	}
			
			
	return (gotOne);
}

/*********************************************************************************
*
*			Name:		subObjGetNextVoxAngCellDecay
*
*			Summary:	Get next voxel angle cell decay.
*
*			Arguments:
*
*			Function return: True if there is another decay to get.
*
*********************************************************************************/
Boolean	subObjGetNextVoxAngCellDecay()
{
	Boolean	gotOne = true;	/* Assume there is another one to get */
	
	do {
		/* Set flag true, contrary condition sets it false */
		gotOne = true;

		/* See if decay left in current voxel */
		if (SubObjDecaySlice[(SubObjCurVoxelIndex*PRODTBLGetNumAngleCells())+SubObjCurAngleIndex] != 0){
			
			/* Decrement the count */
			(SubObjDecaySlice[(SubObjCurVoxelIndex*PRODTBLGetNumAngleCells())+SubObjCurAngleIndex])--;
		}
		else {
			
			/* Go to next angle */
			SubObjCurAngleIndex++;
			
			/* See if through all angles */
			if (SubObjCurAngleIndex == PRODTBLGetNumAngleCells()) {
				SubObjCurAngleIndex = 0;
				SubObjCurVoxelIndex++;
				
				/* See if through all voxels */
				if (SubObjCurVoxelIndex == SubObjGetNumActVoxels(SubObjCurSliceIndex)) {
					SubObjCurVoxelIndex = 0;
					
					/* Increment index to next slice */
					SubObjCurSliceIndex++;
					
					/* See if through all slices */
					if (SubObjCurSliceIndex == SubObjNumSlices) {
						gotOne = false;
					}
					else {
						/* Setup memory for next slice */
						if (subObjAllocDecaySlice(SubObjCurSliceIndex) == false) {
							PhgAbort("Unable to allocate memory for next slice!", true);
						}
							
						/* Calculate time bin decays for new slice */
						subObjCalcSliceTimeBinDecays(subObjCurTimeBin, SubObjCurSliceIndex);
						
						/* Decrement the count */
						if (SubObjDecaySlice[(SubObjCurVoxelIndex*PRODTBLGetNumAngleCells())+SubObjCurAngleIndex]
							!= 0) {
							(SubObjDecaySlice[(SubObjCurVoxelIndex*PRODTBLGetNumAngleCells())+SubObjCurAngleIndex])--;
						}
						else
							gotOne = false;
					}
				}
				else {
							
					/* Decrement the count */
					if (SubObjDecaySlice[(SubObjCurVoxelIndex*PRODTBLGetNumAngleCells())+SubObjCurAngleIndex]
						!= 0) {
						(SubObjDecaySlice[(SubObjCurVoxelIndex*PRODTBLGetNumAngleCells())+SubObjCurAngleIndex])--;
					}
					else
						gotOne = false;
				}
			}
			else {
						
				/* Decrement the count */
				if (SubObjDecaySlice[(SubObjCurVoxelIndex*PRODTBLGetNumAngleCells())+SubObjCurAngleIndex]
					!= 0) {
					(SubObjDecaySlice[(SubObjCurVoxelIndex*PRODTBLGetNumAngleCells())+SubObjCurAngleIndex])--;
				}
				else
					gotOne = false;
			}
		}
	} while ((SubObjCurSliceIndex < SubObjNumSlices) && (gotOne == false));
	
	return (gotOne);
}

/*********************************************************************************
*
*			Name:		SubObjGetCellAttenuation
*
*			Summary:	Get attenuation of given cell.
*			Arguments:
*				LbUsFourByte		sliceIndex	- The slice index of the position.
*				LbUsFourByte		xIndex		- The x index of the position.
*				LbUsFourByte		yIndex		- The y index of the position.
*				double				energy		- The energy of the photon.
*				double				*attenPtr	- Storage for the attenuation.
*
*			Function return: True unless given indexes are outside the object.
*
*********************************************************************************/
Boolean	SubObjGetCellAttenuation(LbFourByte sliceIndex,
			LbFourByte xIndex, LbFourByte yIndex, double energy,
			double *attenPtr)
{
	Boolean			isInside = false;	/* Are the indeces outside */
	LbUsFourByte	index;				/* The tissue index */
	
	/* Set the default value */
	*attenPtr = SUBOBJ_VERYSMALLATTENUATION;

	do { /* Process Loop  */

		#ifdef PHG_DEBUG
		/* Check indexes for object boundaries, then assign attenuation */
			if ((sliceIndex < 0) || ((LbUsFourByte)sliceIndex >= SubObjNumSlices)) {
				LbInPrintf("sliceIndex out of range %d, [%d, %d]\n", sliceIndex, 0, SubObjNumSlices-1);
				break;
			}
			
			else if ((xIndex < 0) || ((LbUsFourByte)xIndex >= SubObjObject[sliceIndex].attNumXBins)) {
				LbInPrintf("xIndex out of range %d, [%d, %d]\n", xIndex, 0, SubObjObject[sliceIndex].attNumXBins-1);
				break;
			}
	
			else if ((yIndex < 0) || ((LbUsFourByte)yIndex >= SubObjObject[sliceIndex].attNumYBins)) {
				LbInPrintf("yIndex out of range %d, [%d, %d]\n", yIndex, 0, SubObjObject[sliceIndex].attNumYBins-1);
				break;
			}
		#endif	
		
 
	   	/* Get tissue index */
		index = SubObjObject[sliceIndex].attenuationArray[(yIndex*SubObjObject[sliceIndex].attNumXBins)+xIndex];
		if (index >= SubObjNumTissues) {
			sprintf(subObjErrStr, "Tissue index '%ld' retrieved for sliceIndex = %ld"
				" xIndex = %ld and yIndex = %ld is greater than number of supported tissues = %ld (SubObjGetCellAttenuation)",
				(unsigned long)index, (unsigned long)sliceIndex, 
				(unsigned long)xIndex,(unsigned long) yIndex, (unsigned long)SubObjNumTissues);
			PhgAbort(subObjErrStr, true);
		}
		
		/* Look up the attenuation */
		SubObjGetAttenuationInObj(index, energy, attenPtr);

		
		/* We made it here so we are inside the object */
		isInside = true;
	} while (false);
	
	#ifdef PHG_DEBUG
		/* If we were outside, this is an error */
		if (isInside == false)
			PhgAbort("SubObjGetCellAttenuation called for voxel indexes with position outside object!.",
				true);
	#endif
	
	return (isInside);
}

/*********************************************************************************
*
*			Name:		SubObjGtAttenuationMaterialName
*
*			Summary:	Get name of attenuation material.
*			Arguments:
*				LbUsFourByte		materialIndex	- The material index.
*				char				*namePtr		
*
*			Function return: char * - Pointer to the name, not to be changed!
*
*********************************************************************************/
char *	SubObjGtAttenuationMaterialName(LbFourByte materialIndex)
{
	return(subObjMaterialNames[materialIndex]);
}

/*********************************************************************************
*
*			Name:		SubObjGetAttenuationInObj
*
*			Summary:	Get attenuation in object of given value.
*			Arguments:
*				LbUsFourByte		materialIndex	- The material index.
*				double				energy			- The energy of the photon.
*				double				*attenPtr		- Storage for the attenuation.
*
*			Function return: None.
*
*********************************************************************************/
void	SubObjGetAttenuationInObj(LbFourByte materialIndex,
			double energy, double *attenPtr)
{
	LbUsFourByte	energyIndex;		/* Energy index */
	
	/* Convert energy to tissue index */
	energyIndex = ((LbUsFourByte) (energy + .499));
		
	#ifdef PHG_DEBUG
		if (subObjIsInitialized == false) {
			PhgAbort("\nYou are trying to call a SubObj routine before it is initialized\n", false);
		}
		
		/* Verify we are not below the minimum supported energy */
		if (energyIndex < PHG_MIN_PHOTON_ENERGY) {
		
			sprintf(subObjErrStr, "Incoming energy below threshold.\t"
				" energy = %f, integer index = %ld (SubObjGetAttenuation).",
				energy, (unsigned long)energyIndex);
				
			PhgAbort(subObjErrStr, false);
		}
	#endif
	
	/* Subtract for lowest supported energy */
	energyIndex -= PHG_MIN_PHOTON_ENERGY;
	
	#ifdef PHG_DEBUG
		/* Prevent array boundary error */
		if (materialIndex >= (LbFourByte)SubObjNumTissues) {
			
			sprintf(subObjErrStr, "Invalid material index for attenuation retrieval.\t"
				" matierialIndex = %ld, max material index = %ld (SubObjGetAttenuation).",
				(unsigned long)materialIndex, (unsigned long)SubObjNumTissues);
				
			PhgAbort(subObjErrStr, false);
		}
				
		if (energyIndex >= SUBOBJ_NUM_ENERGY_BINS) {
			
			sprintf(subObjErrStr, "Invalid energy index for attenuation retrieval.\t"
				" energyIndex = %ld, max energy index = %d (SubObjGetAttenuation).",
				(unsigned long)energyIndex, SUBOBJ_NUM_ENERGY_BINS);

			PhgAbort(subObjErrStr, false);
		}
	#endif
	
	
	/* Get attenuation for given energy */
	if (PHG_IsModelCoherentInObj()){
		*attenPtr = SubObjTissueAttenTableCoh[materialIndex][energyIndex].attenuation;	
	}
	else {
		*attenPtr = SubObjTissueAttenTableNoCoh[materialIndex][energyIndex].attenuation;	
	}
}

/*********************************************************************************
*
*			Name:		SubObjGetAttenuationInTomo
*
*			Summary:	Get attenuation in object of given value.
*			Arguments:
*				LbUsFourByte		materialIndex	- The material index.
*				double				energy			- The energy of the photon.
*				double				*attenPtr		- Storage for the attenuation.
*
*			Function return: None.
*
*********************************************************************************/
void	SubObjGetAttenuationInTomo(LbFourByte materialIndex,
			double energy, double *attenPtr)
{
	LbUsFourByte	energyIndex;		/* Energy index */
	
	/* Convert energy to tissue index */
	energyIndex = ((LbUsFourByte) (energy + .499));
		
	#ifdef PHG_DEBUG
		if (subObjIsInitialized == false) {
			PhgAbort("\nYou are trying to call a SubObj routine before it is initialized\n", false);
		}
		
		/* Verify we are not below the minimum supported energy */
		if (energyIndex < PHG_MIN_PHOTON_ENERGY) {
		
			sprintf(subObjErrStr, "Incoming energy below threshold.\t"
				" energy = %f, integer index = %ld (SubObjGetAttenuation).",
				energy, (unsigned long)energyIndex);
				
			PhgAbort(subObjErrStr, false);
		}
	#endif
	
	/* Subtract for lowest supported energy */
	energyIndex -= PHG_MIN_PHOTON_ENERGY;
	
	#ifdef PHG_DEBUG
		/* Prevent array boundary error */
		if (materialIndex >= (LbFourByte)SubObjNumTissues) {
			
			sprintf(subObjErrStr, "Invalid material index for attenuation retrieval.\t"
				" matierialIndex = %ld, max material index = %ld (SubObjGetAttenuation).",
				(unsigned long)materialIndex, (unsigned long)SubObjNumTissues);
				
			PhgAbort(subObjErrStr, false);
		}
				
		if (energyIndex >= SUBOBJ_NUM_ENERGY_BINS) {
			
			sprintf(subObjErrStr, "Invalid energy index for attenuation retrieval.\t"
				" energyIndex = %ld, max energy index = %d (SubObjGetAttenuation).",
				(unsigned long)energyIndex, SUBOBJ_NUM_ENERGY_BINS);

			PhgAbort(subObjErrStr, false);
		}
	#endif
	
	
	/* Get attenuation for given energy */
	if (PHG_IsModelCoherentInTomo()){
		*attenPtr = SubObjTissueAttenTableCoh[materialIndex][energyIndex].attenuation;	
	}
	else {
		*attenPtr = SubObjTissueAttenTableNoCoh[materialIndex][energyIndex].attenuation;	
	}
}

/*********************************************************************************
*
*			Name:		SubObjGetInnerCellDistance
*
*			Summary:	Compute signed distance to voxel edge(s) in X, Y, & Z 
*						direction.
*			Arguments:
*				PHG_Position		*posPtr			- The position of interest.
*				PHG_Direction		*dirPtr			- Direction photon is traveling.
*				LbUsFourByte		sliceIndex		- The slice index of the position.
*				LbUsFourByte		xIndex			- The x index of the position.
*				LbUsFourByte		*yIndex			- The y index of the position.
*				double				*zDistPtr		- Distance to slice edge.
*				double				*xDistPtr		- Distance to cell edge.
*				double				*yDistPtr		- Distance to cell edge.
*
*			Function return: none.
*
*********************************************************************************/
void	SubObjGetInnerCellDistance(PHG_Position *posPtr, PHG_Direction *dirPtr,
			LbFourByte sliceIndex, LbFourByte xIndex, LbFourByte yIndex,
			double *xDistPtr, double *yDistPtr, double *zDistPtr)
{
		
	/* Calculate distance to edge of slice */
	if (dirPtr->cosine_z >= 0) {
		*zDistPtr = SubObjObject[sliceIndex].zMax - posPtr->z_position;
	}
	else {
		*zDistPtr = SubObjObject[sliceIndex].zMin - posPtr->z_position;
	}
	
	/* Calculate the x distance to voxel edge */
	if (dirPtr->cosine_x >= 0) {
		*xDistPtr = (SubObjObject[sliceIndex].xMin + 
		((xIndex + 1) * SubObjObject[sliceIndex].attVoxelWidth)) - posPtr->x_position;
	}
	else {
		*xDistPtr = (SubObjObject[sliceIndex].xMin + 
		(xIndex * SubObjObject[sliceIndex].attVoxelWidth)) - posPtr->x_position;
	}
		
	/* Calculate the y distance to voxel edge NOTE inverse coordinate direction from X */
	if (dirPtr->cosine_y >= 0) {
		*yDistPtr = (SubObjObject[sliceIndex].yMax - 
		(yIndex * SubObjObject[sliceIndex].attVoxelHeight)) - posPtr->y_position;
	}
	else {
		*yDistPtr = (SubObjObject[sliceIndex].yMax - 
		((yIndex+1) * SubObjObject[sliceIndex].attVoxelHeight)) - posPtr->y_position;
	}

	/* Double check for correct interpretation above */
	#ifdef PHG_DEBUG
		if (dirPtr->cosine_z >= 0) {
			if ((*zDistPtr < 0) || (*zDistPtr > SubObjObject[sliceIndex].sliceDepth)){
				if ((*zDistPtr > SubObjObject[sliceIndex].sliceDepth) ||
						(PhgMathRealNumAreEqual(*zDistPtr, 0, -7, 0,0,0) == false)){
					PhgAbort("Invalid computation for z distance (SubObjGetInnerCellDistance)", true);
				}
			}
		}
		else {
			if ((*zDistPtr > 0) || (*zDistPtr < -SubObjObject[sliceIndex].sliceDepth)){
				if ((*zDistPtr < -SubObjObject[sliceIndex].sliceDepth) ||
						(PhgMathRealNumAreEqual(*zDistPtr, 0, -7, 0,0,0) == false)){
					PhgAbort("Invalid computation for z distance (SubObjGetInnerCellDistance)", true);
				}
			}
		}
		
		if (dirPtr->cosine_x >= 0) {
			if ((*xDistPtr < 0) || (*xDistPtr > SubObjObject[sliceIndex].attVoxelWidth)){
				if ((*xDistPtr > SubObjObject[sliceIndex].attVoxelWidth) ||
						(PhgMathRealNumAreEqual(*xDistPtr, 0, -7, 0,0,0) == false)){
					PhgAbort("Invalid computation for x distance (SubObjGetInnerCellDistance)", true);
				}
			}
		}
		else {
			if ((*xDistPtr > 0) || (*xDistPtr < -SubObjObject[sliceIndex].attVoxelWidth)){
				if ((*xDistPtr < -SubObjObject[sliceIndex].attVoxelWidth) ||
						(PhgMathRealNumAreEqual(*xDistPtr, 0, -7, 0,0,0) == false)){
					PhgAbort("Invalid computation for x distance (SubObjGetInnerCellDistance)", true);
				}
			}
		}
		
		if (dirPtr->cosine_y >= 0) {
			if ((*yDistPtr < 0) || (*yDistPtr > SubObjObject[sliceIndex].attVoxelHeight)){
				if ((*yDistPtr > SubObjObject[sliceIndex].attVoxelHeight) ||
						(PhgMathRealNumAreEqual(*yDistPtr, 0, -7, 0,0,0) == false)){
					PhgAbort("Invalid computation for y distance (SubObjGetInnerCellDistance)", true);
				}
			}
		}
		else {
			if ((*yDistPtr > 0) || (*yDistPtr < -SubObjObject[sliceIndex].attVoxelHeight)){
				if ((*yDistPtr < -SubObjObject[sliceIndex].attVoxelHeight) ||
						(PhgMathRealNumAreEqual(*yDistPtr, 0, -7, 0,0,0) == false)){
					PhgAbort("Invalid computation for y distance (SubObjGetInnerCellDistance)", true);
				}
			}
		}
	#endif

}

/*********************************************************************************
*
*			Name:		SubObjGetObjCylinder
*
*			Summary:	Return parameters of object that define it's containing
*						cylinder.
*			Arguments:
*				double	*radiusPtr	- The radius.
*				double	*zMinPtr	- The minimum z value.
*				double	*zMaxPtr	- The maximum z value.
*				double	*centerX	- The center x value.
*				double	*centerY	- The center y value.
*
*			Function return: none.
*
*********************************************************************************/
void	SubObjGetObjCylinder(double	*radiusPtr, double *zMinPtr,
			double *zMaxPtr, double *centerX, double *centerY)
{
	LbUsFourByte	sliceIndex;	/* Current slice */
	double			minX;		/* Minimum x coordinate */
	double			maxX;		/* Maximum x coordinate */
	double			minY;		/* Minimum y coordinate */
	double			maxY;		/* Maximum y coordinate */
	
	/* Get zMin */
	*zMinPtr = SubObjObject[0].zMin;
	
	/* Get zMax */
	*zMaxPtr = SubObjObject[SubObjNumSlices-1].zMax;
	
	/* Find the other parameters */
	minX = MAXFLOAT;
	maxX = -MAXFLOAT;
	minY = MAXFLOAT;
	maxY = -MAXFLOAT;

	for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
				
		/* Calculate minimum/maximum x values */
		{
			if (minX > SubObjObject[sliceIndex].xMin)
				minX = SubObjObject[sliceIndex].xMin;
			
			if (maxX < SubObjObject[sliceIndex].xMax)
				maxX = SubObjObject[sliceIndex].xMax;
		}
		
		/* Calculate minimum/maximum y values */
		{
			if (minY > SubObjObject[sliceIndex].yMin)
				minY = SubObjObject[sliceIndex].yMin;
			
			if (maxY < SubObjObject[sliceIndex].yMax)
				maxY = SubObjObject[sliceIndex].yMax;
		}
	}

	*centerX = maxX - fabs((maxX - minX)/2);
	*centerY = maxY - fabs((maxY - minY)/2);
	
	*radiusPtr = (maxX - minX)/2;
}

/*********************************************************************************
*
*			Name:		subObjCalcSliceActivity
*
*			Summary:	Calculate the activity for the entire slice.
*
*			Arguments:
*				LbUsFourByte	sliceIndex - The slice to calculate.
*
*			Function return: None.
*
*********************************************************************************/
void	subObjCalcSliceActivity(LbUsFourByte sliceIndex)
{
	double 			totalActivity = 0.0;	/* Cumulative activity for slice */
	double			voxelSize;				/* Size of voxel */
	LbUsFourByte	tissueIndex;			/* Current tissue index */
	LbUsFourByte	xIndex;					/* X axis voxel bin */
	LbUsFourByte	yIndex;					/* Y axis voxel bin */
	
	
	/* Calculate voxel size */
	voxelSize = SubObjObject[sliceIndex].actVoxelWidth *
		SubObjObject[sliceIndex].actVoxelHeight * 
		SubObjObject[sliceIndex].sliceDepth;
		
	/* Sum activity over all voxels */
	for (yIndex = 0; yIndex < SubObjObject[sliceIndex].actNumYBins; yIndex++) {
		for (xIndex = 0; xIndex < SubObjObject[sliceIndex].actNumXBins; xIndex++) {
		
			/* Retreive tissue index for voxel */
			tissueIndex = SubObjGetTissueIndex(sliceIndex, yIndex, xIndex);

			/* Look up activity and add to running total */
			totalActivity += (SubObjGetTissueActivity(tissueIndex, 0) *
					voxelSize * SubObjCurTimeBinDuration);
		}
	}	
	
	/* Set the slice's activity */
	SubObjObject[sliceIndex].sliceActivity = totalActivity;
}

/*********************************************************************************
*
*			Name:		SubObjGetStartingProdValues
*
*			Summary:	Create the productivity table.
*
*			Arguments:
*				ProdTblProdTblInfoTy	*prodTableInfoPtr - Information for object creation.
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean	SubObjGetStartingProdValues(ProdTblProdTblInfoTy *prodTableInfoPtr)
{
	Boolean 		okay = false;	/* Process Flag */
	double			zMax;			/* Max z value */
	double			zMin;			/* Min z value */
	LbUsFourByte	sliceIndex;	/* Current slice */
	
	do {	/* Process Loop */
	
		/* Initialize our values of the productivity table info */
		prodTableInfoPtr->numSlices = SubObjNumSlices;
		
		/* Create the productivity table */
		if (!ProdTblCreateTable(prodTableInfoPtr)) {
			break;
		}
	
		/* Verify that z value's match */
		for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
			
			/* Get productivity table z values */
			ProdTblGetZmaxZmin(sliceIndex, &zMax, &zMin);
			
			/* Make sure they match if stratification is on */
			if (PHG_IsStratification()) {
				if ((SubObjObject[sliceIndex].zMax != zMax) ||
						(SubObjObject[sliceIndex].zMin != zMin)) {
					
					ErStGeneric("Non-matching z values between object and productivity table.");
					goto FAIL;
				}
			}
		}
		
		okay = true;
		FAIL:;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			SubObjInitialize
*
*			Summary:		Initialize the sub object code.
*
*			Arguments:
*
*			Function return: TRUE unless an error occurs.
*
*********************************************************************************/
Boolean SubObjInitialize()	
{
	do { /* Process Loop */
	
		/* Set the length of the current time bin */
		SubObjCurTimeBinDuration = (double)PHGGetLengthOfScan();
		
		/* Clear this pointer for memory recovery procedures */
		SubObjTissueTable.tissueValues = 0;
		
		SubObjDecaysProcessed = 0;
		
		
		/* Do some verification on subobject related parameters */
		if (PHG_IsVoxelPointSource() && PHG_IsVoxelLineSource()) {
			ErStGeneric("You have specified point source voxels AND line source voxels,\n"
				"please set one or both of these to false in the parameter file.");
			break;
		}
		
		subObjIsInitialized = true;
	} while (false);
		
	return (subObjIsInitialized);
}

/*********************************************************************************
*
*			Name:		SubObjGetCellPosRangeConstants
*
*			Summary:	Get positron range constants of given cell.
*			Arguments:
*				LbUsFourByte		sliceIndex	- The slice index of the position.
*				LbUsFourByte		xIndex		- The x index of the position.
*				LbUsFourByte		yIndex		- The y index of the position.
*				double				*b1Ptr		- Storage for positron range constant.
*				double				*b2Ptr		- Storage for positron range constant.
*				double				*densityPtr	- Storage for the density.
*
*			Function return: True unless given indexes are outside the object.
*
*********************************************************************************/
Boolean	SubObjGetCellPosRangeConstants(LbFourByte sliceIndex,
			LbFourByte xIndex, LbFourByte yIndex, double *b1Ptr,
			double *b2Ptr, double *densityPtr)

{
	Boolean			isInside = false;	/* Are the indeces outside */
	LbUsFourByte	index;				/* The tissue index */
	double			atomicNumber;		/* effective atomic number */
	double			atomicWeight;		/* effective atomic weight */
	
	/* Set the default values */
	*b1Ptr = -1;
	*b2Ptr = -1;
	*densityPtr = -1;

	do { /* Process Loop  */

		#ifdef PHG_DEBUG
		/* Check indexes for object boundaries, then assign attenuation */
			if ((sliceIndex < 0) || ((LbUsFourByte)sliceIndex >= SubObjNumSlices)) {
				LbInPrintf("sliceIndex out of range %d, [%d, %d]\n", sliceIndex, 0, SubObjNumSlices-1);
				break;
			}
			
			else if ((xIndex < 0) || ((LbUsFourByte)xIndex >= SubObjObject[sliceIndex].attNumXBins)) {
				LbInPrintf("xIndex out of range %d, [%d, %d]\n", xIndex, 0, SubObjObject[sliceIndex].attNumXBins-1);
				break;
			}
	
			else if ((yIndex < 0) || ((LbUsFourByte)yIndex >= SubObjObject[sliceIndex].attNumYBins)) {
				LbInPrintf("yIndex out of range %d, [%d, %d]\n", yIndex, 0, SubObjObject[sliceIndex].attNumYBins-1);
				break;
			}
		#endif	
		
 
	   	/* Get tissue index */
		index = SubObjObject[sliceIndex].attenuationArray[(yIndex*SubObjObject[sliceIndex].attNumXBins)+xIndex];

		#ifdef PHG_DEBUG
		if (index >= SubObjNumTissues) {
			sprintf(subObjErrStr, "Tissue index '%ld' retrieved for sliceIndex = %ld"
				" xIndex = %ld and yIndex = %ld is greater than number of supported tissues = %ld (SubObjGetCellAttenuation)",
				(unsigned long)index, (unsigned long)sliceIndex, 
				(unsigned long)xIndex, (unsigned long)yIndex, (unsigned long)SubObjNumTissues);
			PhgAbort(subObjErrStr, true);
		}
		#endif	
		
		/* Get the density, atomic number and weight */
		*densityPtr = subObjMaterialDAZ[index].D;
		atomicNumber = subObjMaterialDAZ[index].Z;
		atomicWeight = subObjMaterialDAZ[index].A;
		
		if (*densityPtr < 0.0) {
			sprintf(subObjErrStr, "Attempt to compute positron range constants for material not supported for positron range.");
			PhgAbort(subObjErrStr, true);
		}

		/* Compute the positron range constants */
		subObjComputePosRangeConstants( atomicNumber, atomicWeight, b1Ptr, b2Ptr );

		/* We made it here so we are inside the object */
		isInside = true;
	} while (false);
	
	#ifdef PHG_DEBUG
		/* If we were outside, this is an error */
		if (isInside == false)
			PhgAbort("SubObjGetCellAttenuation called for voxel indexes with position outside object!.",
				true);
	#endif
	
	return (isInside);
}

/*********************************************************************************
*
*			Name:		SubObjGetWaterPosRangeConstants
*
*			Summary:	Get positron range constants of water.
*			Arguments:
*				double				*b1Ptr		- Storage for positron range constant.
*				double				*b2Ptr		- Storage for positron range constant.
*				double				*densityPtr	- Storage for the density.
*
*			Function return: None.
*
*********************************************************************************/
void SubObjGetWaterPosRangeConstants(double *b1Ptr, double *b2Ptr, double *densityPtr)

{
	LbUsFourByte	index;				/* The tissue index */
	double			atomicNumber;		/* effective atomic number */
	double			atomicWeight;		/* effective atomic weight */
	
	/* Set the default values */
	atomicNumber = -1;
	atomicWeight = -1;
	*densityPtr = -1;
	
	/* We assume that water is at index 1 */
	index = 1;
	/* And complain if it isn't */
	if (strcmp(SubObjGtAttenuationMaterialName(index), "water") != 0) {
		PhgAbort("You must have material '1' set to water in your attenuation table to model positron range", false);
	}


	do { /* Process Loop  */
		
		/* Get the density, atomic number and weight */
		*densityPtr = subObjMaterialDAZ[index].D;
		atomicNumber = subObjMaterialDAZ[index].Z;
		atomicWeight = subObjMaterialDAZ[index].A;
		
		if (*densityPtr < 0.0) {
			sprintf(subObjErrStr, "Water positron range constants not provided.");
			PhgAbort(subObjErrStr, true);
		}

		/* Compute the positron range constants */
		subObjComputePosRangeConstants( atomicNumber, atomicWeight, b1Ptr, b2Ptr );

		/* We made it here so we are inside the object */
	} while (false);
	
	return;
}

/*********************************************************************************
*
*			Name:		subObjComputePosRangeConstants
*
*			Summary:	Compute b1 and b2 using formulaes 4 from Palmer and Brownell,
*						IEEE TMI 11:3:373-378, 1992.
*			Arguments:
*				double				atomicNumber	- effective atomic number.
*				double				atomicWeight	- effective atomic weight.
*				double				*b1Ptr			- Storage for positron range constant.
*				double				*b2Ptr			- Storage for positron range constant.
*
*			Function return: True unless given indexes are outside the object.
*
*********************************************************************************/
void	subObjComputePosRangeConstants(	double atomicNumber, double atomicWeight,
										double *b1Ptr, double *b2Ptr )

{
	*b1Ptr = (4.569 * atomicWeight) / pow( atomicNumber,1.209 );
	*b2Ptr = 1 / (2.873 - 0.02309 * atomicNumber);
}


/*********************************************************************************
*
*			Name:			SubObjTerminate
*
*			Summary:		Terminate the sub object code.
*
*			Arguments:
*
*			Function return: None.
*
*********************************************************************************/
void SubObjTerminate()	
{
	LbUsFourByte	sliceIndex;				/* Current slice */
	LbUsFourByte	tissueIndex;			/* The number of voxels */
	
	/* Only do this if we were initialized */
	if (subObjIsInitialized) {
		/* Free the activity arrays */
		if (SubObjObject != 0) {
			for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {
			
				if (SubObjObject[sliceIndex].activityArray != 0) {
					LbMmFree((void **) &(SubObjObject[sliceIndex].activityArray));
				}
			
				if (SubObjObject[sliceIndex].attenuationArray != 0) {
					LbMmFree((void **) &(SubObjObject[sliceIndex].attenuationArray));
				}
			}
	
			LbMmFree((void **) &SubObjObject);
		}
		
		/* Free the final slice of the decay count and weight arrays */
		if (SubObjDecaySlice != 0){
			LbMmFree((void **)&(SubObjDecaySlice));
		}
		if (SubObjDecayWeightSlice != 0){
			LbMmFree((void **)&(SubObjDecayWeightSlice));
		}
	
		if (SubObjTissueTable.tissueValues != 0) {
			for (tissueIndex = 0; tissueIndex < SubObjNumActIndexes; tissueIndex++) {
				if (SubObjTissueTable.tissueValues[tissueIndex].activityValues != 0)
					LbMmFree((void **)&(SubObjTissueTable.tissueValues[tissueIndex].activityValues));
			}
			LbMmFree((void **)&(SubObjTissueTable.tissueValues));
		}
		
		if (SubObjTissueAttenTableNoCoh != 0)
			LbMmFree((void **)&(SubObjTissueAttenTableNoCoh));
		
		if (SubObjTissueAttenTableCoh != 0)
			LbMmFree((void **)&(SubObjTissueAttenTableCoh));
		
		if (subObjCohScatAngles != 0)
			LbMmFree((void **)&(subObjCohScatAngles));
			
		/* Clear our initialization flag */
		subObjIsInitialized = false;
	}
	
}

/*********************************************************************************
*
*			Name:			subObjVoxelCellGetNumReal
*
*			Summary:		Get number of real decays for a given voxel/angle cell
*							and time.
*
*			Arguments:
*				LbUsFourByte	sliceIndex		- Current slice.
*				LbUsFourByte	voxelIndex		- Current voxel.
*				LbUsFourByte	angleIndex		- Current angle.
*				LbUsFourByte	currentTime		- Current time.
*
*			Function return: Number of real events to occur at the given location.
*
*********************************************************************************/
double subObjVoxelCellGetNumReal(LbUsFourByte sliceIndex, LbUsFourByte angleIndex,
		LbUsFourByte voxelIndex, LbUsFourByte currentTime)	
{
	double	numRealDecays;	/* Number of real decays for a given voxel */
	double	voxelSize;		/* Size of voxels in this slice */
	double	voxelActivity;	/* Activity value for the voxel */
	
	/* Calculate voxel size */
	voxelSize = SubObjObject[sliceIndex].actVoxelWidth *
		SubObjObject[sliceIndex].actVoxelHeight *
		SubObjObject[sliceIndex].sliceDepth;

	/* Calculate the activity value for the voxel */
	voxelActivity = (SUBOBJ_DECAYS_PER_CURIE * SubObjGetTissueActivity(
		SubObjObject[sliceIndex].activityArray[voxelIndex], currentTime)
		* voxelSize * SubObjCurTimeBinDuration);
		
	/* Calculate the number of real decays for the voxel/cell */
	numRealDecays = (voxelActivity *
		PRODTBLGetProdTblAngleSize(sliceIndex, angleIndex))/2;
		
	return (numRealDecays);
}

#ifdef PHG_DEBUG
/*********************************************************************************
*
*			Name:			SubObjDumpObjects
*
*			Summary:		Dump the objects for debuging purposes.
*
*			Arguments:
*
*			Function return: None.
*
*********************************************************************************/
void SubObjDumpObjects()	
{
/* #ifdef FOR_SEVER_DEBUGGING_ONLY */
	Boolean					okay = false;	/* Process Flag */
	LbUsTwoByte				sliceIndex;		/* Current slice */
	LbUsTwoByte				xIndex;			/* Current x index */
	LbUsTwoByte				yIndex;			/* Current y index */
	
	do { /* Process Loop */


		/* Print out the slice information */
		for (sliceIndex = 0; sliceIndex < SubObjNumSlices; sliceIndex++) {

			/* Print slice info */
			LbInPrintf("\nSlice # %d", SubObjObject[sliceIndex].sliceNum);
			
			LbInPrintf("\n xMin\t= %03.2lf \t\t\txMax\t= %03.2lf \tsliceWidth\t= %03.2lf\n",
				SubObjObject[sliceIndex].xMin, SubObjObject[sliceIndex].xMax, SubObjObject[sliceIndex].sliceWidth);

			LbInPrintf(" yMin\t= %03.2lf \t\t\tyMax\t= %03.2lf \tsliceHeight\t= %03.2lf\n zMin\t= %03.2lf \t\t\tzMax\t= %03.2lf",
				SubObjObject[sliceIndex].yMin, SubObjObject[sliceIndex].yMax, SubObjObject[sliceIndex].sliceHeight,
				SubObjObject[sliceIndex].zMin, SubObjObject[sliceIndex].zMax);
				
			LbInPrintf("\tsliceDepth\t= %03.2lf\n-\n num horizontal activity voxels\t\t= %ld\n num vertical activity voxels\t\t\t= %ld\n",
				SubObjObject[sliceIndex].sliceDepth, 
				(unsigned long)(SubObjObject[sliceIndex].actNumXBins), 
				(unsigned long)(SubObjObject[sliceIndex].actNumYBins));
				
			LbInPrintf("\n-\n num horizontal attenuation voxels\t\t= %ld\n num vertical attenuation voxels\t\t\t= %ld\n",
				(unsigned long)(SubObjObject[sliceIndex].attNumXBins), 
				(unsigned long)(SubObjObject[sliceIndex].attNumYBins));
				
			LbInPrintf(" activity voxel width\t\t\t\t\t= %03.2lf\n activity voxel height\t\t\t\t\t= %03.2lf\n\n",
				SubObjObject[sliceIndex].actVoxelWidth, SubObjObject[sliceIndex].actVoxelHeight);
			
			LbInPrintf(" attenuation voxel width\t\t\t\t\t= %03.2lf\n attenuation voxel height\t\t\t\t\t= %03.2lf\n\n",
				SubObjObject[sliceIndex].attVoxelWidth, SubObjObject[sliceIndex].attVoxelHeight);

			/* Print seperator */		
			LbInPrintf("\nTissue indexes by voxel.\n");
			
			/* Print out tissue indexes; Loop through each y index */
			for (yIndex = 0; yIndex < SubObjObject[sliceIndex].actNumYBins; yIndex++) {
				for (xIndex = 0; xIndex < SubObjObject[sliceIndex].actNumXBins; xIndex++) {
		
					/* Print out tissue index */
					LbInPrintf("%d ",
						SubObjObject[sliceIndex].activityArray[
						(yIndex * SubObjObject[sliceIndex].actNumXBins) + xIndex]);
				}
				LbInPrintf("\n");
			}

			/* Print seperator */		
			LbInPrintf("\n\nVoxel attenuation coefficients and scatter to cross section ratios.\n");
			
			/* Print out  attenuation; Loop through each y index */
			for (yIndex = 0; yIndex < SubObjObject[sliceIndex].attNumYBins; yIndex++) {				
				for (xIndex = 0; xIndex < SubObjObject[sliceIndex].attNumXBins; xIndex++) {
					
					/* Print out attenuation index */
					LbInPrintf("(%ld) ",
						(unsigned long)(SubObjObject[sliceIndex].attenuationArray[
						(yIndex * SubObjObject[sliceIndex].attNumXBins) + xIndex]));
				}
				LbInPrintf("\n");
			}
						
			/* Print the blank line */
			LbInPrintf("\n\n");
			
		} /* Loop to next slice */
		

		/* Print out round up ratio */
		LbInPrintf("\nNumber of voxel/cell counts rounded up %3.2e.\n", SubObjAngleRoundUpCount);
		
		okay = true;
	} while (false);
	
	if (!okay) {
		ErHandle("Dump failed.", false);
		PhgAbort("Program terminated (SubObjDumpObjects).", true);
	}
/* #endif */
}
#endif /* PHG_DEBUG */


#undef SUB_OBJECT
