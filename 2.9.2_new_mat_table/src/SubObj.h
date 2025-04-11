/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1992-2006 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		SubObj.h
*			Revision Number:	1.1
*			Date last revised:	2007
*			Programmer:			Steven Vannoy
*			Date Originated:	20 August 1992
*
*			Module Overview:	Definitions for SubObj.c.
*
*			References:			'Sub Object Processes' PHG design.
*
**********************************************************************************
*
*			Global functions defined:	
*
*			Global variables defined:		none
*
*			Global Macros
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
*			Revision description:	Moved SubObjNumTissues here from SubObj.c
*				and added SUBOBJGetNumTissues.
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
*********************************************************************************/

#ifndef SUB_OBJ
#define SUB_OBJ

#ifdef SUB_OBJECT
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */
#define SUBOBJ_DECAYS_PER_CURIE			37000000000.0		/* Number of decays per second per Curie */
#ifdef PHG_DEBUG
	#define SUBOBJ_FIXED_DIR_COSINE_X	1.0					/* For fixed direction debug option */
	#define SUBOBJ_FIXED_DIR_COSINE_Y	0.0					/* For fixed direction debug option */
	#define SUBOBJ_FIXED_DIR_COSINE_Z	0.0					/* For fixed direction debug option */
#endif	

/* TYPES */
	/* Activity Object */
typedef	LbUsFourByte				SubObjActVoxelTy;		/* A voxel is a tissue type, is an index */


	/* Time activity curve */
typedef struct {
	double	activityLevel;							/* Activity in curries */
	int		duration;								/* Number of seconds for this activity level */
} SubObjTACBinTy;
typedef SubObjTACBinTy		*SubObjTACArrayTy;		/* Array of TACs */

	/* Tissue */
typedef struct {
	LbFourByte			totalDuration;				/* Total duration of TAC */
	LbFourByte			numberOfBins;				/* Number of time bins */
	SubObjTACArrayTy	activityValues;				/* Values for the activity */
} SubObjTissueElemTy, *SubObjTissueArrayTy;

	/* Tissue Table */
typedef struct {
	SubObjTissueArrayTy	tissueValues;				/* Array of tissues */
} SubObjTissueTableTy;


	/* Slice Information */
typedef struct {
	LbUsFourByte		sliceNum;			/* Number of slice */
	double				zMin;				/* Minimum z coordinate */
	double				zMax;				/* Maximum z coordinate */
	double				xMin;				/* Minimum x coordinate */
	double				xMax;				/* Maximum x coordinate */
	double				yMin;				/* Minimum y coordinate */
	double				yMax;				/* Maximum y coordinate */
	double				sliceWidth;			/* Width of the slice */
	double				sliceHeight;		/* Height of the slice */
	double				sliceDepth;			/* Depth of slice */
	double				actVoxelWidth;		/* Width of a voxel in activity object */
	double				actVoxelHeight;		/* Height of a voxel (in activity object */
	LbUsFourByte		actNumXBins;		/* Number of x bins in activity object */
	LbUsFourByte		actNumYBins;		/* Number of y bins in activity object */
	double				attVoxelWidth;		/* Width of a voxel in attenuation object */
	double				attVoxelHeight;		/* Height of a voxel in attenuation object */
	LbUsFourByte		attNumXBins;		/* Number of x bins in attenuation object */
	LbUsFourByte		attNumYBins;		/* Number of y bins in attenuation object */
	double				sliceActivity;		/* Total activity in the slice */
	SubObjActVoxelTy	*activityArray;		/* The voxel activity information */
	LbUsFourByte		*attenuationArray;	/* The voxel attenuation information */
} SubObjSliceInfoTy;
typedef SubObjSliceInfoTy	*SubObjSliceInfoArrayTy;

/* GLOBALS */
LOCALE	LbUsFourByte				SubObjNumSlices;			/* Number of slices in object */
LOCALE	SubObjSliceInfoArrayTy		SubObjObject;				/* Our Object */
LOCALE	LbUsFourByte				SubObjCurSliceIndex;		/* Current slice index */
LOCALE	LbUsFourByte				SubObjCurVoxelIndex;		/* Current voxel index */
LOCALE	LbUsFourByte				SubObjCurAngleIndex;		/* Current angle index */
LOCALE	LbUsEightByte				SubObjDecaysProcessed;		/* Number of decays processed */
LOCALE	LbUsFourByte				SubObjNumTissues;			/* Number of tissues in attenuation table */
LOCALE	double						SubObjTotalRealDecays;		/* Total expected real-world decays */
/* MACROS */		

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceMinX
*
*			Summary:	Return minimum x value for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: double, the slice's minimum x value.
*
*********************************************************************************/
#define SUBOBJGetSliceMinX(sliceIndex) (SubObjObject[sliceIndex].xMin)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceMinY
*
*			Summary:	Return minimum y value for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: double, the slice's minimum y value.
*
*********************************************************************************/
#define SUBOBJGetSliceMinY(sliceIndex) (SubObjObject[sliceIndex].yMin)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceMinZ
*
*			Summary:	Return minimum z value for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: double, the slice's minimum z value.
*
*********************************************************************************/
#define SUBOBJGetSliceMinZ(sliceIndex) (SubObjObject[sliceIndex].zMin)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceMaxX
*
*			Summary:	Return maximum x value for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: double, the slice's maximum x value.
*
*********************************************************************************/
#define SUBOBJGetSliceMaxX(sliceIndex) (SubObjObject[sliceIndex].xMax)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceMaxY
*
*			Summary:	Return maximum y value for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: double, the slice's maximum y value.
*
*********************************************************************************/
#define SUBOBJGetSliceMaxY(sliceIndex) (SubObjObject[sliceIndex].yMax)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceMaxZ
*
*			Summary:	Return maximum x value for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: double, the slice's maximum z value.
*
*********************************************************************************/
#define SUBOBJGetSliceMaxZ(sliceIndex) (SubObjObject[sliceIndex].zMax)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceMaxXIndex
*
*			Summary:	Return maximum x index for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: LbUsFourByte, the slice's maximum x index.
*
*********************************************************************************/
#define SUBOBJGetSliceMaxXIndex(sliceIndex) (SubObjObject[sliceIndex].numXBins -1)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceMaxYIndex
*
*			Summary:	Return maximum y index for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: LbUsFourByte, the slice's maximum y index.
*
*********************************************************************************/
#define SUBOBJGetSliceMaxYIndex(sliceIndex) (SubObjObject[sliceIndex].numYBins -1)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceNumAttXBins
*
*			Summary:	Return number of x bins.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: LbUsFourByte, the number of x bins.
*
*********************************************************************************/
#define SUBOBJGetSliceAttNumXBins(sliceIndex) (SubObjObject[sliceIndex].attNumXBins)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceAttNumYBins
*
*			Summary:	Return number of y bins.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: LbUsFourByte, the number of y bins.
*
*********************************************************************************/
#define SUBOBJGetSliceAttNumYBins(sliceIndex) (SubObjObject[sliceIndex].attNumYBins)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceAttVoxelWidth
*
*			Summary:	Return voxel width for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: double, the slice's voxel width.
*
*********************************************************************************/
#define SUBOBJGetSliceAttVoxelWidth(sliceIndex) (SubObjObject[sliceIndex].attVoxelWidth)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceAttVoxelHeight
*
*			Summary:	Return voxel height for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: double, the slice's voxel height.
*
*********************************************************************************/
#define SUBOBJGetSliceAttVoxelHeight(sliceIndex) (SubObjObject[sliceIndex].attVoxelHeight)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceVoxelDepth
*
*			Summary:	Return voxel depth for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: double, the slice's voxel depth.
*
*********************************************************************************/
#define SUBOBJGetSliceVoxelDepth(sliceIndex) (SubObjObject[sliceIndex].sliceDepth)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceWidth
*
*			Summary:	Return width for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: double, the slice's width.
*
*********************************************************************************/
#define SUBOBJGetSliceWidth(sliceIndex) (SubObjObject[sliceIndex].sliceWidth)

/*********************************************************************************
*
*			Name:		SUBOBJGetSliceSliceHeight
*
*			Summary:	Return height for given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex - The current slice.
*				
*			Function return: double, the slice's height.
*
*********************************************************************************/
#define SUBOBJGetSliceSliceHeight(sliceIndex) (SubObjObject[sliceIndex].sliceHeight)

/*********************************************************************************
*
*			Name:		SUBOBJGetTissueIndex
*
*			Summary:	Return tissue index for photon's current location.
*			Arguments:
*				PHG_TrackingPhotonPtr	trPhoton.
*				
*			Function return: double, the slice's height.
*
*********************************************************************************/
#define SUBOBJGetTissueIndex(trPhoton) (SubObjObject[trPhoton->sliceIndex].attenuationArray[(trPhoton->yIndex*SubObjObject[trPhoton->sliceIndex].attNumXBins)+trPhoton->xIndex])

/*********************************************************************************
*
*			Name:		SUBOBJGetDecaysProcessed
*
*			Summary:	Return number of decays processed.
*			Arguments:
*				
*			Function return: double, the slice's height.
*
*********************************************************************************/
#define SUBOBJGetDecaysProcessed() (SubObjDecaysProcessed)

/*********************************************************************************
*
*			Name:		SUBOBJGetNumTissues
*
*			Summary:	Return number of different attenuation materials.
*			Arguments:
*				
*			Function return: double, the slice's height.
*
*********************************************************************************/
#define SUBOBJGetNumTissues() (SubObjNumTissues)

/*********************************************************************************
*
*			Name:		SUBOBJGetTotalRealDecays
*
*			Summary:	Return expected number of decays in a real world scan
*				using simulation setup.
*			Arguments:
*				
*			Function return: double, the slice's height.
*
*********************************************************************************/
#define SUBOBJGetTotalRealDecays() (SubObjTotalRealDecays)

/* PROTOTYPES */
void	SubObjGtPositionIndexes(PHG_Position *posPtr,
			LbFourByte *sliceIndexPtr,
			LbFourByte *xIndexPtr,
			LbFourByte *yIndexPtr);
void	SubObjCalcTimeBinDecays(LbUsFourByte currentTimeBin,
			LbUsEightByte numberOfDecaysToSimulate);
Boolean	SubObjCreate(void);
Boolean	SubObjGenVoxAngCellDecay(PHG_Decay *newDecayPtr,
			PHG_Direction *newDecayEmissionAnglePtr, LbFourByte *sliceIndexPtr,
			LbFourByte *angleIndexPtr,  LbFourByte *xIndexPtr,
			LbFourByte *yIndexPtr);
Boolean	SubObjGetCellAttenuation(LbFourByte sliceIndex, LbFourByte xIndex,
			LbFourByte yIndex, double energy, double *attenPtr);
void	SubObjGetAttenuationInObj(LbFourByte materialIndex,
			double energy, double *attenPtr);
void	SubObjGetAttenuationInTomo(LbFourByte materialIndex,
			double energy, double *attenPtr);
void	SubObjGetInnerCellDistance(PHG_Position *posPtr, PHG_Direction *dirPtr,
			LbFourByte	 sliceIndex, LbFourByte xIndex, LbFourByte yIndex,
			double *xDistPtr, double *yDistPtr, double *zDistPtr);
void	SubObjGetObjCylinder(double	*radiusPtr, double *zMinPtr,
			double *zMaxPtr, double *centerX, double *centerY);
Boolean	SubObjGetObjGeometry(LbPfHkTy paramFlHk, LbUsFourByte numParams);
double	SubObjGetProbComptToScatter(PHG_TrackingPhoton *trackingPhotonPtr);
double	SubObjGetProbComptToScatInTomo2(LbUsFourByte index, double energy);
double	SubObjGetProbComptToScatInObj2(LbUsFourByte index, double energy);
double	SubObjGetProbScatterInObj(PHG_TrackingPhoton *trackingPhotonPtr);
double	SubObjGetProbScatterInObj2(LbUsFourByte index, double energy);
double	SubObjGetProbScatterInTomo(PHG_TrackingPhoton *trackingPhotonPtr);
double	SubObjGetProbScatterInTomo2(LbUsFourByte index, double energy);
double	SubObjGetProbScatter(LbUsFourByte index, double energy, 
			Boolean modelingCohScatter);
double	SubObjGetProbComptonCondnl(LbUsFourByte index, double energy, 
			Boolean modelingCohScatter);
char *	SubObjGtAttenuationMaterialName(LbFourByte materialIndex);
Boolean	SubObjGetStartingProdValues(ProdTblProdTblInfoTy *prodTableInfoPtr);
Boolean	SubObjInitialize(void);
void	SubObjTerminate(void);
double	SubObjGetCohTheta(PHG_TrackingPhoton *trPhoton);
double	SubObjGetCohTheta2(LbUsFourByte materialIndex, double energy);
LbFourByte		SubObjGtMaterialIndex(char *theMaterialStr);

Boolean			SubObjGetCellPosRangeConstants(LbFourByte sliceIndex,
					LbFourByte xIndex, LbFourByte yIndex, double *b1Ptr,
					double *b2Ptr, double *densityPtr);
void 			SubObjGetWaterPosRangeConstants(double *b1Ptr, double *b2Ptr,
					double *densityPtr);

#ifdef PHG_DEBUG
void	SubObjDumpObjects(void);
#endif
#undef LOCALE
#endif /* SUB_OBJ */

