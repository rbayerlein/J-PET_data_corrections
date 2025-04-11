/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 2002-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/

/*********************************************************************************
*
*			Module Name:		DetBlock.c
*			Revision Number:	3.27
*			Date last revised:	23 July 2013
*			Programmer:			Steven Gillispie
*			Date Originated:	9 August 2002
*
*			Module Overview:	Simulates block detector functionality
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:
*
*			Global variables defined:	None
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
*			Revision date:		4 June 2012
*
*			Revision description:	Modified DetBlocFindDetPosition to respond to 
*										new BlockDetectedPositionAlgo parameter
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		11 February 2010
*
*			Revision description:	Changed to new block parameter fields
*									Corrected 3D dir cos misuses in detBlocCrossZoneBounds
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		7 October 2008 - 20 May 2009
*
*			Revision description:	Added changes for non-homogeneous blocks
*
*********************************************************************************/

#define DETECTOR_BLOCK


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
#include "Lb2DGeometry.h"

#include "Photon.h"
#include "PhgParams.h"
#include "PhgMath.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "CylPos.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "ProdTbl.h"
#include "SubObj.h"
#include "EmisList.h"
#include "PhoHStat.h"
#include "PhoHFile.h"
#include "ColUsr.h"
#include "UNCCollimator.h"
#include "Collimator.h"
#include "ColSlat.h"
#include "DetGeometric.h"
#include "DetCylinder.h"
#include "Detector.h"
#include "DetUsr.h"
#include "phg.h"
#include "PhgBin.h"

#include "DetBlock.h"


/* TYPES */
typedef enum {							/* Possible photon destinations */
	detBlocEvt_Null, 						/* Never set value; should never occur */
	detBlocEvt_OutEnd, 						/* Exited cylinder end */
	detBlocEvt_OuterCyl, 					/* Hit outer cylinder */
	detBlocEvt_InnerCyl, 					/* Hit inner cylinder */
	detBlocEvt_NextRing, 					/* Went into next ring */
	detBlocEvt_PrevRing, 					/* Went into previous ring */
	detBlocEvt_NextZone, 					/* Went into next zone */
	detBlocEvt_PrevZone, 					/* Went into previous zone */
	detBlocEvt_Block						/* Hit a detector block */
} DetBlocEvents;
typedef struct  {						/* Angle specifier */
	double	xDirCos;						/* X direction cosine of angle line */
	double	yDirCos;						/* Y direction cosine of angle line */
} DetAngleSpec;
typedef struct  {						/* Arc range specifier */
	DetAngleSpec		minAngle;			/* Minimum angle of subtended arc */
	DetAngleSpec		maxAngle;			/* Maximum angle of subtended arc */
} ArcRangeSpec;
typedef struct  {						/* Internal record for a detector block */
	LbUsFourByte		itsRing;			/* Ring number of block */
	LbUsFourByte		itsBlock;			/* Block number of block */
	Lb2D_Rect			itsRect;			/* 2D rect of the block in (x,y) plane */
	double				zMin;				/* Minimum z coordinate of block */
	double				zMax;				/* Maximum z coordinate of block */
	DetAngleSpec		minAngle;			/* Minimum angle of subtended arc */
	DetAngleSpec		maxAngle;			/* Maximum angle of subtended arc */
	DetBlockTomoBlockTy	*userBlockData;		/* Block data specified by user */
} DetectorBlock;
typedef DetectorBlock*	DetectorBlockPtr;

/* GLOBALS */
static char				detBlocErrStr[1024];		/* Storage for creating error strings */
static LbFourByte		detBlocNumRings = 0;		/* Number of rings */
static LbFourByte		detBlocMaxRingBlocks = 0;	/* Maximum number of blocks in a ring */
static DetectorBlockPtr detBlocBlocksDatabase = NULL;/* Internal database by (r,b) coords */
													/*	= [totalRings][maxRingBlocks] */
static LbFourByte		detBlocNumRingZones = 0;	/* Number of zones per ring */
static DetAngleSpec		*detBlocZoneBounds = NULL;	/* List of ring zone boundaries */
static LbFourByte		detBlocNumZoneBlocks = 0;	/* Maximum number of blocks per zone */
static DetectorBlockPtr* detBlocIndexDatabase = NULL;/* Internal database by (r,z,i) coords */
													/*	= [totalRings][#zones][zoneBlocks] */
static CylPosCylinderTy	detBlocInnerCylinder;		/* Inner cylindrical boundary */
static double			detBlocOutRadius = 0.0;		/* Radius of outer cylinder */
static double			detBlocOutRadSqrd = 0.0;	/* Radius of outer cylinder (squared) */

/* PROTOTYPES */
void	detBlocGetZoneRange(	LbFourByte			ringNum, 
								LbFourByte			zoneNum, 
								ArcRangeSpec		*zoneArcRange);
Boolean detBlocGetZone( PHG_Position	*thePosition, 
						LbFourByte		*ringNum, 
						LbFourByte		*zoneNum);
Boolean detBlocGetBlock(LbFourByte			ringNum, 
						LbFourByte			zoneNum, 
						LbFourByte			blockIndex, 
						DetectorBlockPtr	*blockPtr);
Boolean detBlocGetNumberedBlock(LbFourByte			ringNum, 
								LbFourByte			blockNum, 
								DetectorBlockPtr	*blockPtr);
Boolean	detBlocGetElementIndex(	PHG_Position *thePosition, 
								DetElementPosIndexTy *posIndex);
Boolean detBlocGetElementCorners(	DetElementPosIndexTy	*elementDataPtr, 
									PHG_Position			*corner1, 
									PHG_Position			*corner2);
Boolean detBlocGetElementCenter(	DetElementPosIndexTy	*elementDataPtr, 
									PHG_Position			*itsCenterPtr);
void detBlocTomoToBlockCoord(	DetBlockTomoBlockTy	*curBlockInfoPtr, 
								PHG_Position		*tomoPos, 
								PHG_Position		*blockPos);
void detBlocBlockToTomoCoord(	DetBlockTomoBlockTy	*curBlockInfoPtr, 
								PHG_Position		*blockPos, 
								PHG_Position		*tomoPos);
void detBlocTomoToBlockDirection(	DetBlockTomoBlockTy	*curBlockInfoPtr, 
									PHG_Direction		*tomoDir, 
									PHG_Direction		*blockDir);
void detBlocBlockToTomoDirection(	DetBlockTomoBlockTy	*curBlockInfoPtr, 
									PHG_Direction		*blockDir, 
									PHG_Direction		*tomoDir);
Boolean detBlocProjAxially(	PHG_TrackingPhoton	*thePhotonPtr, 
							double				zDest);
Boolean detBlocProjAcrossGap(PHG_TrackingPhoton	*photonPtr);
void detBlocGetInnerCylinder(CylPosCylinderTy	*theCylinder);
Boolean detBlocPtInInnerCylinder(	CylPosCylinderTy	*innerCyl, 
									Lb2D_Point			*testPoint);
Boolean detBlocXInnerCylinder(	CylPosCylinderTy	*innerCyl, 
								Lb2D_Point			*point1, 
								Lb2D_Point			*point2);
Boolean detBlocProjectToCylinder(	PHG_Position *positionPtr,
									PHG_Direction *directionPtr, CylPosCylinderTy *cylinderPtr, 
									PHG_Position *newPosPtr, double *distPtr);
Boolean detBlocHitInnerCylinder(	CylPosCylinderTy	*innerCyl, 
									PHG_TrackingPhoton	*photonPtr, 
									PHG_Position		*intersectPt, 
									double				*distance);
Boolean detBlocPtInOuterCylinder(Lb2D_Point		*testPoint);
Boolean detBlocHitOuterCylinder(	PHG_TrackingPhoton	*photonPtr, 
									PHG_Position		*intersectPt, 
									double				*distance);
DetBlocEvents detBlocCrossZoneBounds(	PHG_TrackingPhoton	*thePhotonPtr, 
										LbFourByte			itsRing, 
										LbFourByte			itsZone, 
										PHG_Position		*intersectPt, 
										double				*distance);
Boolean detBlocCalcDistanceToBlock(	PHG_Position		*photonPosition, 
									PHG_Direction		*photonDirection, 
									DetectorBlockPtr	theBlockPtr, 
									double				*theDistance, 
									PHG_Position		*hitPoint);
Boolean detBlocMakeRBDatabase(void);
double detBlocIntraFreePaths(	PHG_TrackingPhoton	*thePhotonPtr, 
								DetectorBlockPtr	theBlockPtr, 
								double				travelDistance);
Boolean detBlocSelfConsCheck(void);
Boolean detBlocOverlapCheck(void);
void detBlocIntraDistance(	PHG_TrackingPhoton	*thePhotonPtr, 
							DetectorBlockPtr	theBlockPtr, 
							double				freePaths, 
							double				maxTravelDist, 
							double				*travelDistance, 
							double				*freePathsUsed, 
							LbUsFourByte		*detMaterial, 
							Boolean				*isActive);
DetBlocEvents detBlocGetDistToExit(	PHG_TrackingPhoton	*thePhotonPtr, 
									DetectorBlockPtr	theBlockPtr, 
									PHG_Position		*blockExitPoint, 
									double				*travelDist);
double	detBlocCompFreePathsToExit(PHG_TrackingPhoton *photonPtr, 
								LbUsFourByte curRingIndex, LbUsFourByte curDetIndex);
Boolean detBlocNextBlock(	PHG_TrackingPhoton	*thePhotonPtr, 
							LbFourByte			itsRing, 
							LbFourByte			itsZone, 
							DetectorBlockPtr	itsBlockPtr, 
							DetectorBlockPtr	*nextBlockPtr, 
							LbFourByte			*blockRing, 
							LbFourByte			*blockZone, 
							LbFourByte			*blockIndex);
Boolean detBlocIntraSearchCentroidLayer(	PHG_Position			*fixedPositionPtr, 
											DetectorBlockPtr		theBlockPtr, 
											DetElementPosIndexTy	*closestElemDataPtr, 
											double					*closestDistance);
Boolean detBlocIntraFindCentroid(	PHG_Position		*origCentroid, 
									LbFourByte			itsRing, 
									LbFourByte			itsBlock, 
									PHG_Position		*newCentroid, 
									LbFourByte			*crystalNum);


/*********************************************************************************
*
*		Name:			detBlocGetZoneRange
*
*		Summary:		Return arc range of a specified zone.
*
*		Arguments:
*			LbFourByte		   ringNum			- Ring number of the zone.
*			LbFourByte		   zoneNum			- Zone number of the zone.
*			ArcRangeSpec	   *zoneArcRange	- Arc range of the zone.
*
*		Function return: None.
*
*********************************************************************************/

void detBlocGetZoneRange(	LbFourByte			ringNum, 
							LbFourByte			zoneNum, 
							ArcRangeSpec		*zoneArcRange)

{
	/* NOTE:  At this time, all rings use the same zones, but this could change 
		in the future, so the ring is still passed as a parameter */
	if (ringNum) {};			/* Remove unused variable compiler warning */
	
	zoneArcRange->minAngle = detBlocZoneBounds[zoneNum];
	zoneArcRange->maxAngle = detBlocZoneBounds[(zoneNum+1) % detBlocNumRingZones];
}


/*********************************************************************************
*
*		Name:			detBlocGetZone
*
*		Summary:		Return the ring and zone of the given point.
*						Unless (r,z) are (-1,-1), assume (r,z) are already close.
*
*		Arguments:
*			PHG_Position		*thePosition	- Given position.
*			LbFourByte			*ringNum		- Ring number of the position.
*			LbFourByte			*zoneNum		- Zone number of the position.
*
*		Function return: True if determined, false if not inside a ring.
*
*********************************************************************************/

Boolean detBlocGetZone( PHG_Position	*thePosition, 
						LbFourByte		*ringNum, 
						LbFourByte		*zoneNum)

{
	Boolean				validRing;			/* Function return */
	LbFourByte			r;					/* Index through checked rings */
	DetBlockTomoRingTy	*ringPtr;			/* Pointer to current ring info */
	double				ringMin;			/* Current ring minimum axial edge */
	double				ringMax;			/* Current ring maximum axial edge */
	LbFourByte			itsRing;			/* Returned ring number */
	LbFourByte			delta;				/* Increment in ring checking */
	LbFourByte			itsZone;			/* Returned zone number */
	
	
	validRing = false;
	if (*ringNum >= 0) {
		/* Has supplied value; check it */
		r = *ringNum;
		if (r < detBlocNumRings) {
			ringPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[r]);
			ringMin = ringPtr->MinZ + ringPtr->AxialShift;
			ringMax = ringPtr->MaxZ + ringPtr->AxialShift;
			if (thePosition->z_position >= ringMin) {
				if (thePosition->z_position <= ringMax) {
					/* The point is in this ring */
					validRing = true;
					itsRing = r;
				}
				else {
					/* Check higher */
					r++;
					delta = 1;
				}
			}
			else {
				/* Check lower */
				r--;
				delta = -1;
			}
		}
		else {
			/* Invalid ring number */
			r = 0;
			delta = 1;
		}
	}
	else {
		/* No starting guess */
		r = 0;
		delta = 1;
	}
	
	if (! validRing) {
		/* Determine the ring by comparing the z with each ring */
		while ((0<=r) && (r<detBlocNumRings)) {
			ringPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[r]);
			ringMin = ringPtr->MinZ + ringPtr->AxialShift;
			ringMax = ringPtr->MaxZ + ringPtr->AxialShift;
			if ((ringMin <= thePosition->z_position) && 
					(thePosition->z_position <= ringMax)) {
				/* The point is in this ring */
				validRing = true;
				itsRing = r;
				break;
			}
			r += delta;
		}
	}
	
	
	if (validRing) {
		/* Determine the zone */
		
		DetAngleSpec		itsAngle;			/* Angle of the position */
		LbFourByte			z;					/* Index through checked zones */
		ArcRangeSpec		zoneArcRange;		/* Arc range of checked zones */
		
		
		/* Get the 2D arc of the position */
		{
			Lb2D_Point			origin;				/* (0,0) */
			
			origin.x_position = 0.0;
			origin.y_position = 0.0;
			Lb2DDirCosines(	&origin, (Lb2D_Point*)thePosition, 
								&itsAngle.xDirCos, &itsAngle.yDirCos);
		}
		
		itsZone = -1;
		if (*zoneNum >= 0) {
			/* Has supplied value; check it */
			detBlocGetZoneRange(itsRing, *zoneNum, &zoneArcRange);
			if (Lb2DDirCosComp(itsAngle.xDirCos, itsAngle.yDirCos, 
								zoneArcRange.minAngle.xDirCos, zoneArcRange.minAngle.yDirCos) 
					>= 0) {
				if (Lb2DDirCosComp(itsAngle.xDirCos, itsAngle.yDirCos, 
								zoneArcRange.maxAngle.xDirCos, zoneArcRange.maxAngle.yDirCos) 
						<= 0) {
					/* Position is within this zone */
					itsZone = *zoneNum;
				}
				else {
					/* Check higher */
					z = *zoneNum + 1;
					if (z >= detBlocNumRingZones) {
						z = 0;
					}
					delta = 1;
				}
			}
			else {
				/* Check lower */
				z = *zoneNum - 1;
				if (z < 0) {
					z = detBlocNumRingZones-1;
				}
				delta = -1;
			}
		}
		else {
			/* No starting guess */
			z = 0;
			delta = 1;
		}
		
		if (itsZone == -1) {
			/* Check the zones within the ring */
			do {
				detBlocGetZoneRange(itsRing, z, &zoneArcRange);
				if ((Lb2DDirCosComp(itsAngle.xDirCos, itsAngle.yDirCos, 
								zoneArcRange.minAngle.xDirCos, zoneArcRange.minAngle.yDirCos) 
						>= 0) && 
					(Lb2DDirCosComp(itsAngle.xDirCos, itsAngle.yDirCos, 
								zoneArcRange.maxAngle.xDirCos, zoneArcRange.maxAngle.yDirCos) 
						<= 0)) {
					
					/* Position is within this zone */
					itsZone = z;
					break;
				}
				
				z += delta;
				if (z < 0) {
					z = detBlocNumRingZones-1;
				}
				else if (z >= detBlocNumRingZones) {
					z = 0;
				}
			} while (true);
		}
	}
	
	
	/* Return the results */
	if (validRing) {
		*ringNum = itsRing;
		*zoneNum = itsZone;
	}
	return (validRing);
}


/*********************************************************************************
*
*		Name:			detBlocGetBlock
*
*		Summary:		Return a pointer to the requested block, by (r,z,i).
*
*		Arguments:
*			LbFourByte			ringNum			- Ring number of the block.
*			LbFourByte			zoneNum			- Zone number of the block.
*			LbFourByte			blockIndex		- Zone block index of the block.
*			DetectorBlockPtr	*blockPtr		- Pointer to the detector block;
*													Null if not found.
*
*		Function return: Boolean; true if found, false if not.
*
*********************************************************************************/

Boolean detBlocGetBlock(LbFourByte			ringNum, 
						LbFourByte			zoneNum, 
						LbFourByte			blockIndex, 
						DetectorBlockPtr	*blockPtr)

{
	LbUsFourByte		arrayIndex;			/* Index into the index blocks database */
	Boolean				blockFound;			/* Return value */
	
	
	/* Access the index database for the (r,z,i) requested block */
	if ((0<=ringNum) && (ringNum<detBlocNumRings) && 
			(0<=zoneNum) && (zoneNum<detBlocNumRingZones) && 
			(0<=blockIndex) && (blockIndex<detBlocNumZoneBlocks)) {
		arrayIndex = (ringNum * detBlocNumRingZones  +  zoneNum) * detBlocNumZoneBlocks  +  
						blockIndex;
		*blockPtr = detBlocIndexDatabase[arrayIndex];
		if (*blockPtr) {
			if ((*blockPtr)->userBlockData) {
				blockFound = true;
			}
			else {
				*blockPtr = NULL;
				blockFound = false;
			}
		}
		else {
			blockFound = false;
		}
	}
	else {
		*blockPtr = NULL;
		blockFound = false;
	}
	
	return (blockFound);
}


/*********************************************************************************
*
*		Name:			detBlocGetNumberedBlock
*
*		Summary:		Return a pointer to the requested block by (r,b).
*
*		Arguments:
*			LbFourByte			ringNum			- Ring number of the block.
*			LbFourByte			blockNum		- User block number of the block.
*			DetectorBlockPtr	*blockPtr		- Pointer to the detector block; 
*													Null if not found.
*
*		Function return: True if found, false if not.
*
*********************************************************************************/

Boolean detBlocGetNumberedBlock(LbFourByte			ringNum, 
								LbFourByte			blockNum, 
								DetectorBlockPtr	*blockPtr)

{
	Boolean				blockFound;			/* Return value */
	LbFourByte			theBlockIndex;		/* Index of the block in the database */
	DetectorBlockPtr	theBlockPtr;		/* Pointer to the block */
	
	
	/* Access the blocks database for the (r,b) requested block */
	blockFound = false;
	if ((0<=ringNum) && (ringNum<detBlocNumRings) && 
			(0<=blockNum) && (blockNum<detBlocMaxRingBlocks)) {
		theBlockIndex = ringNum * detBlocMaxRingBlocks  +  blockNum;
		theBlockPtr = &(detBlocBlocksDatabase[theBlockIndex]);
		if (theBlockPtr->userBlockData) {
			/* Valid database block */
			blockFound = true;
			*blockPtr = theBlockPtr;
		}
	}
	if (! blockFound) {
		*blockPtr = NULL;
	}
	
	return (blockFound);
}


/*********************************************************************************
*
*		Name:			detBlocGetElementIndex
*
*		Summary:		Fill in the layer and element in the position index for 
*							the supplied position.
*
*		Arguments:
*			PHG_Position			*thePosition	- Pointer to supplied position.
*			DetElementPosIndexTy	*posIndex		- Returned position index, with 
*														ring and block supplied.
*
*		Function return: True if found, false if not.
*
*********************************************************************************/

Boolean detBlocGetElementIndex(	PHG_Position			*thePosition, 
								DetElementPosIndexTy	*posIndex)

{
	const double		kTolerance = 1E-12;	/* Block edges tolerance of thePosition */
	
	Boolean				valid;				/* Return value */
	DetectorBlockPtr	theBlockPtr;		/* Pointer to thePosition's block */
	DetBlockTomoLayerTy	*layerPtr;			/* Current layer */
	LbFourByte			maxLayer;			/* Index of outermost layer in block */
	LbFourByte			i;					/* Index through layers and elements */
	LbUsFourByte		yElement;			/* Index of element in y direction */
	LbUsFourByte		zElement;			/* Index of element in z direction */
	LbUsFourByte		elementNum;			/* Calculated element number */
	
	
	/* NOTE:  thePosition is assumed to be in local block coordinates;
				of course, this is irrelevant for homogeneous blocks.  
			  It is assumed that thePosition is meant to be inside a block. */
	
	
	/* Get the supplied block */
	valid = detBlocGetNumberedBlock(posIndex->ringNum, posIndex->blockNum, 
										&theBlockPtr);
	
	if (valid) {
		/* Find the layer corresponding to the position */
		if (theBlockPtr->userBlockData->NumLayers == 1) {
			/* Only one layer */
			posIndex->layerNum = 0;
		}
		else {
			/* Multiple layers in the block */
			valid = false;
			layerPtr = &(theBlockPtr->userBlockData->LayerInfo[0]);
			if (fabs(thePosition->x_position - layerPtr->InnerX) < kTolerance) {
				posIndex->layerNum = 0;
				valid = true;
			}
			else {
				maxLayer = theBlockPtr->userBlockData->NumLayers-1;
				layerPtr = &(theBlockPtr->userBlockData->LayerInfo[maxLayer]);
				if (fabs(thePosition->x_position - layerPtr->OuterX) < kTolerance) {
					posIndex->layerNum = maxLayer;
					valid = true;
				}
				else {
					for (i=0; i<=maxLayer; i++) {
						layerPtr = &(theBlockPtr->userBlockData->LayerInfo[i]);
						if ((layerPtr->InnerX <= thePosition->x_position) && 
								(thePosition->x_position < layerPtr->OuterX)) {
							/* Found the layer */
							posIndex->layerNum = i;
							valid = true;
							break;
						}
					}
				}
			}
		}
		
		if (valid) {
			/* Find the element corresponding to the position */
			if (theBlockPtr->userBlockData->LayerInfo[posIndex->layerNum].NumElements == 1) {
				/* Only one element */
				posIndex->elementNum = 0;
			}
			else {
				/* Multiple elements in the layer */
				valid = false;
				if (((theBlockPtr->userBlockData->YMin - kTolerance) <= thePosition->y_position) && 
						(thePosition->y_position <= (theBlockPtr->userBlockData->YMax + kTolerance))) {
					layerPtr = &(theBlockPtr->userBlockData->LayerInfo[posIndex->layerNum]);
					yElement = layerPtr->NumYChanges;	/* Default */
					for (i=0; i<layerPtr->NumYChanges; i++) {
						if (thePosition->y_position < layerPtr->YChanges[i]) {
							/* Found the correct y element */
							yElement = i;
							break;
						}
					}
					
					if (((theBlockPtr->userBlockData->ZMin - kTolerance) <= thePosition->z_position) && 
							(thePosition->z_position <= (theBlockPtr->userBlockData->ZMax + kTolerance))) {
						zElement = layerPtr->NumZChanges;	/* Default */
						for (i=0; i<layerPtr->NumZChanges; i++) {
							if (thePosition->z_position < layerPtr->ZChanges[i]) {
								/* Found the correct z element */
								zElement = i;
								break;
							}
						}
						
						/* The (y,z) position was within the block */
						valid = true;
						
						/* Calculate and set the element number */
						elementNum = zElement * (layerPtr->NumYChanges+1)  +  yElement;
						if (elementNum >= layerPtr->NumElements) {
							/* Element number exceeds number of elements */
							valid = false;
						}
						else {
							posIndex->elementNum = elementNum;
						}
					}
				}
			}
		}
	}
	
	return (valid);
}


/*********************************************************************************
*
*		Name:			detBlocGetElementCorners
*
*		Summary:		Return two diagonal corners of the supplied element.
*
*		Arguments:
*			DetElementPosIndexTy	*elementDataPtr	- Supplied element data pointer. 
*			PHG_Position			*corner1		- Pointer to first returned corner.
*			PHG_Position			*corner2		- Pointer to second returned corner.
*
*		Function return: True if found, false if not (invalid element).
*
*********************************************************************************/

Boolean detBlocGetElementCorners(	DetElementPosIndexTy	*elementDataPtr, 
									PHG_Position			*corner1, 
									PHG_Position			*corner2)

{
	Boolean				valid;				/* Return value */
	DetectorBlockPtr	itsBlockPtr;		/* Pointer to the element's block */
	DetBlockTomoBlockTy	*blockPtr;			/* Pointer to the block's user data */
	DetBlockTomoLayerTy	*layerPtr;			/* Pointer to the block layer data */
	LbFourByte			yIndex;				/* Index of element in y direction */
	LbFourByte			zIndex;				/* Index of element in z direction */
	
	
	/* Get pointers to the element's data */
	valid = detBlocGetNumberedBlock(elementDataPtr->ringNum, elementDataPtr->blockNum, 
										&itsBlockPtr);
	if (valid) {
		blockPtr = itsBlockPtr->userBlockData;
		layerPtr = &(blockPtr->LayerInfo[elementDataPtr->layerNum]);
		
		/* Set the x coordinates */
		corner1->x_position = layerPtr->InnerX;
		corner2->x_position = layerPtr->OuterX;
		
		/* Calculate the y and z indices */
		yIndex = elementDataPtr->elementNum % (layerPtr->NumYChanges+1);
		zIndex = elementDataPtr->elementNum / (layerPtr->NumYChanges+1);
		
		/* Find the coordinate of the lesser y surface */
		if ( yIndex == 0 ) {
			corner1->y_position = blockPtr->YMin;
		}
		else {
			corner1->y_position = layerPtr->YChanges[yIndex-1];
		}
		
		/* Find the coordinate of the greater y surface */
		if ( yIndex == layerPtr->NumYChanges ) {
			corner2->y_position = blockPtr->YMax;
		}
		else {
			corner2->y_position = layerPtr->YChanges[yIndex];
		}
		
		/* Find the coordinate of the lesser z surface */
		if ( zIndex == 0 ) {
			corner1->z_position = blockPtr->ZMin;
		}
		else {
			corner1->z_position = layerPtr->ZChanges[zIndex-1];
		}
		
		/* Find the coordinate of the greater z surface */
		if ( zIndex == layerPtr->NumZChanges ) {
			corner2->z_position = blockPtr->ZMax;
		}
		else {
			corner2->z_position = layerPtr->ZChanges[zIndex];
		}
	}
	
	return (valid);
}


/*********************************************************************************
*
*		Name:			detBlocGetElementCenter
*
*		Summary:		Find the center of the supplied element.
*
*		Arguments:
*			DetElementPosIndexTy	*elementDataPtr	- Supplied element data pointer. 
*			PHG_Position			*itsCenterPtr	- Pointer to returned center.
*
*		Function return: True if found, false if not (invalid element).
*
*********************************************************************************/

Boolean detBlocGetElementCenter(	DetElementPosIndexTy	*elementDataPtr, 
									PHG_Position			*itsCenterPtr)

{
	Boolean				valid;				/* Return value */
	PHG_Position		corner1;			/* First corner of the element */
	PHG_Position		corner2;			/* Second corner of the element */
	
	
	/* Get the corners of the element */
	valid = detBlocGetElementCorners(elementDataPtr, &corner1, &corner2);
	
	if (valid) {
		/* Calculate the center coordinates */
		itsCenterPtr->x_position = (corner1.x_position + corner2.x_position) / 2.0;
		itsCenterPtr->y_position = (corner1.y_position + corner2.y_position) / 2.0;
		itsCenterPtr->z_position = (corner1.z_position + corner2.z_position) / 2.0;
	}
	
	return (valid);
}


/*********************************************************************************
*
*		Name:			detBlocTomoToBlockCoord
*
*		Summary:		Transform tomograph (x,y,z) coordinates of a point
*						to block (x,y,z) coordinates.
*
*		Arguments:
*			DetBlockTomoBlockTy	*curBlockInfoPtr	- Current block info.
*			PHG_Position		*tomoPos			- Point in tomo coords.
*			PHG_Position		*blockPos			- Point in block coords.
*
*		Function return: None.
*
*********************************************************************************/

void detBlocTomoToBlockCoord(	DetBlockTomoBlockTy	*curBlockInfoPtr, 
								PHG_Position		*tomoPos, 
								PHG_Position		*blockPos)

{
	double		cosAlpha;	/* Cosine of block face angle */ 
	double		sinAlpha;	/* Sine of block face angle */ 
	double		xTomo;		/* Reference point offset values of blockPos */
	double		yTomo;
	double		zTomo;
	
	
	/* Get rotation angle values from tomograph to block detector coords */
	cosAlpha = curBlockInfoPtr->CosBlockFace;
	sinAlpha = curBlockInfoPtr->SinBlockFace;
	
	/* Shift the tomograph coordinates from the block position */
	xTomo = tomoPos->x_position - curBlockInfoPtr->XPositionTomo;
	yTomo = tomoPos->y_position - curBlockInfoPtr->YPositionTomo;
	zTomo = tomoPos->z_position - curBlockInfoPtr->ZPositionTomo;
	
	/* Rotate and translate into block detector coords */
	blockPos->x_position = (cosAlpha * xTomo) + (sinAlpha * yTomo) + curBlockInfoPtr->XRef;
	blockPos->y_position = -(sinAlpha * xTomo) + (cosAlpha * yTomo) + curBlockInfoPtr->YRef;
	blockPos->z_position = zTomo + curBlockInfoPtr->ZRef;
}


/*********************************************************************************
*
*		Name:			detBlocBlockToTomoCoord
*
*		Summary:		Transform block (x,y,z) coordinates of a point
*						to tomograph (x,y,z) coordinates.
*
*		Arguments:
*			DetBlockTomoBlockTy	*curBlockInfoPtr	- Current block info.
*			PHG_Position		*blockPos			- Point in block coords.
*			PHG_Position		*tomoPos			- Point in tomo coords.
*
*		Function return: None.
*
*********************************************************************************/

void detBlocBlockToTomoCoord(	DetBlockTomoBlockTy	*curBlockInfoPtr, 
								PHG_Position		*blockPos, 
								PHG_Position		*tomoPos)

{
	double		cosAlpha;	/* Cosine of block face angle */ 
	double		sinAlpha;	/* Sine of block face angle */ 
	double		xBlock;		/* Reference point offset values of blockPos */
	double		yBlock;
	double		zBlock;
	
	
	/* Get rotation angle values from block detector to tomograph coords */
	cosAlpha = curBlockInfoPtr->CosBlockFace;
	sinAlpha = - curBlockInfoPtr->SinBlockFace;
	
	/* Shift the block detector coordinates from the block reference point */
	xBlock = blockPos->x_position - curBlockInfoPtr->XRef;
	yBlock = blockPos->y_position - curBlockInfoPtr->YRef;
	zBlock = blockPos->z_position - curBlockInfoPtr->ZRef;
	
	/* Rotate and translate into tomo coords */
	tomoPos->x_position = (cosAlpha * xBlock) + (sinAlpha * yBlock) + curBlockInfoPtr->XPositionTomo;
	tomoPos->y_position = -(sinAlpha * xBlock) + (cosAlpha * yBlock) + curBlockInfoPtr->YPositionTomo;
	tomoPos->z_position = zBlock + curBlockInfoPtr->ZPositionTomo;
}


/*********************************************************************************
*
*		Name:			detBlocTomoToBlockDirection
*
*		Summary:		Transform tomograph cosine_(x,y,z) directions of an angle 
*						to block cosine_(x,y,z) directions.
*
*		Arguments:
*			DetBlockTomoBlockTy	*curBlockInfoPtr	- Current block info.
*			PHG_Direction		*tomoDir			- Direction in tomo coordinates.
*			PHG_Direction		*blockDir			- Direction in block coordinates.
*
*		Function return: None.
*
*********************************************************************************/

void detBlocTomoToBlockDirection(	DetBlockTomoBlockTy	*curBlockInfoPtr, 
									PHG_Direction		*tomoDir, 
									PHG_Direction		*blockDir)

{
	double		cosAlpha;	/* Cosine of block face angle */ 
	double		sinAlpha;	/* Sine of block face angle */ 
	
	
	/* Get rotation angle values from tomograph to block detector coords */
	cosAlpha = curBlockInfoPtr->CosBlockFace;
	sinAlpha = curBlockInfoPtr->SinBlockFace;
	
	/* Convert the tomograph directions to block detector directions:  
		cos(a - b) = cos(a) cos(b)  +  sin(a) sin(b) 
		sin(a - b) = sin(a) cos(b)  -  cos(a) sin(b)  */
	blockDir->cosine_x = tomoDir->cosine_x * cosAlpha  +  tomoDir->cosine_y * sinAlpha;
	blockDir->cosine_y = tomoDir->cosine_y * cosAlpha  -  tomoDir->cosine_x * sinAlpha;
	blockDir->cosine_z = tomoDir->cosine_z;		/* No change in z-direction */
}


/*********************************************************************************
*
*		Name:			detBlocBlockToTomoDirection
*
*		Summary:		Transform block cosine_(x,y,z) directions of an angle 
*						to tomograph cosine_(x,y,z) directions.
*
*		Arguments:
*			DetBlockTomoBlockTy	*curBlockInfoPtr	- Current block info.
*			PHG_Direction		*blockDir			- Direction in block coordinates.
*			PHG_Direction		*tomoDir			- Direction in tomo coordinates.
*
*		Function return: None.
*
*********************************************************************************/

void detBlocBlockToTomoDirection(	DetBlockTomoBlockTy	*curBlockInfoPtr, 
									PHG_Direction		*blockDir, 
									PHG_Direction		*tomoDir)

{
	double		cosAlpha;	/* Cosine of block face angle */ 
	double		sinAlpha;	/* Sine of block face angle */ 
	
	
	/* Get rotation angle values from tomograph to block detector coords */
	cosAlpha = curBlockInfoPtr->CosBlockFace;
	sinAlpha = curBlockInfoPtr->SinBlockFace;
	
	/* Convert the block detector directions to tomograph directions:  
		cos(a + b) = cos(a) cos(b)  -  sin(a) sin(b) 
		sin(a + b) = sin(a) cos(b)  +  cos(a) sin(b)  */
	blockDir->cosine_x = tomoDir->cosine_x * cosAlpha  -  tomoDir->cosine_y * sinAlpha;
	blockDir->cosine_y = tomoDir->cosine_y * cosAlpha  +  tomoDir->cosine_x * sinAlpha;
	blockDir->cosine_z = tomoDir->cosine_z;		/* No change in z-direction */
}


/*********************************************************************************
*
*		Name:			detBlocProjAxially
*
*		Summary:		Project the photon to the supplied axial destination.
*
*		Arguments:
*			PHG_TrackingPhoton	*thePhotonPtr		- The photon.
*			double				zDest				- Axial (z) destination.
*
*		Function return: Boolean--true if moved, false if not moving axially.
*
*********************************************************************************/

Boolean detBlocProjAxially(	PHG_TrackingPhoton	*thePhotonPtr, 
							double				zDest)

{
	PHG_Direction		direction;		/* Local copy of photon direction */
	PHG_Position		*locationPtr;	/* Pointer to photon position */
	double				distance;		/* Projected z distance */
	Boolean				moved;			/* Function return value */
	double				travelDist;		/* Distance traveled in (x,y,z) space */
	
	
	direction = thePhotonPtr->angle;
	locationPtr = &(thePhotonPtr->location);
	distance = zDest - locationPtr->z_position;
	moved = false;
	/* If cosine_z == 0, will never intersect */
	if (fabs(direction.cosine_z) > 1E-12) {
		if (distance != 0.0) {
			travelDist = distance / direction.cosine_z;
			locationPtr->x_position += travelDist * direction.cosine_x;
			locationPtr->y_position += travelDist * direction.cosine_y;
			locationPtr->z_position = zDest;
			thePhotonPtr->travel_distance += travelDist;
			moved = true;
		}
	}
	
	return (moved);
}


/*********************************************************************************
*
*		Name:			detBlocProjAcrossGap
*
*		Summary:		Project a photon across any gap between rings.
*
*		Arguments:
*			PHG_TrackingPhoton	*photonPtr			- The photon.
*
*		Function return: Boolean--true if moved, false if can't enter new ring.
*
*********************************************************************************/

Boolean detBlocProjAcrossGap(PHG_TrackingPhoton	*photonPtr)

{
	double				zPos;				/* Photon z position */
	double				zDistance;			/* Axial distance to nearest ring */
	LbFourByte			r;					/* Index through rings */
	DetBlockTomoRingTy	*ringPtr;			/* Pointer to current ring info */
	double				zRing;				/* Axial boundary of nearest ring */
	Boolean				valid;				/* True if photon moved inside a ring */
	
	
	/* NOTE:  It is assumed that the photon is not currently in a ring */
	
	zPos = photonPtr->location.z_position;
	if (photonPtr->angle.cosine_z > 0.0) {
		/* Moving axially positive */
		/* Determine next ring to enter */
		zDistance = -1.0;
		for (r=0; r<detBlocNumRings; r++) {
			ringPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[r]);
			zRing = ringPtr->MinZ + ringPtr->AxialShift;
			if (zPos < zRing) {
				zDistance = zRing - zPos;
				break;
			}
		}
		if (zDistance < 0.0) {
			/* Leaving the detector axially at the positive end */
			valid = false;
		}
		else {
			/* Project the photon to its closest ring */
			detBlocProjAxially(photonPtr, zRing);
			
			valid = true;
		}
	}
	else if (photonPtr->angle.cosine_z < 0.0) {
		/* Moving axially negative */
		/* Determine next ring to enter */
		zDistance = 1.0;
		for (r=detBlocNumRings-1; r>=0; r--) {
			ringPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[r]);
			zRing = ringPtr->MaxZ + ringPtr->AxialShift;
			if (zPos > zRing) {
				zDistance = zRing - zPos;
				break;
			}
		}
		if (zDistance > 0.0) {
			/* Leaving the detector axially at the negative end */
			valid = false;
		}
		else {
			/* Project the photon to its closest ring */
			detBlocProjAxially(photonPtr, zRing);
			
			valid = true;
		}
	}
	else {
		/* Only moving radially between rings, leaving the detector */
		valid = false;
	}
	
	if (valid) {
		/* Make sure the photon is still inside the outer detector cylinder */
		if (! detBlocPtInOuterCylinder((Lb2D_Point*)&(photonPtr->location))) {
			valid = false;
		}
		else {
			/* Make sure the photon is not inside the inner detector cylinder */
			if (detBlocPtInInnerCylinder(&detBlocInnerCylinder, 
											(Lb2D_Point*)&(photonPtr->location))) {
				valid = false;
			}
		}
	}
	
	return (valid);
}


/*********************************************************************************
*
*		Name:			detBlocGetInnerCylinder
*
*		Summary:		Return the (outermost) collimator or object cylinder.
*
*		Arguments:
*			CylPosCylinderTy	*theCylinder	- Outer (forbidden zone) cylinder.
*
*		Function return: None.
*
*********************************************************************************/

void detBlocGetInnerCylinder(CylPosCylinderTy		*theCylinder)

{
	/* ### RH will be writing new function for this */
	
	/* NOTE:  This is for circular cylinders only */
	
	theCylinder->radius = ColGtOutsideRadius();		/* ### Should check object, also */
	theCylinder->zMin = -10000.0;					/* ### ??  Arbitrary for now */
	theCylinder->zMax = 10000.0;					/* ### ??  Arbitrary for now */
	theCylinder->centerX = 0.0;						/* ### ??  Always correct?? */
	theCylinder->centerY = 0.0;						/* ### ??  Always correct?? */
}


/*********************************************************************************
*
*		Name:			detBlocPtInInnerCylinder
*
*		Summary:		Test whether the point is inside the collimator cylinder.
*
*		Arguments:
*			CylPosCylinderTy	*innerCyl		- The inner cylinder.
*			Lb2D_Point			*testPoint		- Position to be tested.
*
*		Function return: True if inside, false otherwise.
*
*********************************************************************************/

Boolean detBlocPtInInnerCylinder(	CylPosCylinderTy	*innerCyl, 
									Lb2D_Point			*testPoint)

{
	double			colmatrRadius;		/* Collimator radius */
	double			colmatrRSqrd;		/* Collimator radius (squared) */
	double			rSqd;				/* Point radial distance (squared) */
	Boolean			result;				/* Returned function result */
	
	
	/* ### RH will be writing new function for this */
	
	/* NOTE:  This is for circular cylinders only */
	/* NOTE:  Can be simplified by eliminating intermediate variables */
	
	/* Get the outer edge of the collimator (inner bound of the detectors) */
	colmatrRadius = innerCyl->radius;
	colmatrRSqrd = PHGMATH_Square(colmatrRadius);	/* Will work with R^2 (faster) */
	
	/* Get the point distance (squared) */
	rSqd = PHGMATH_Square(testPoint->x_position) + PHGMATH_Square(testPoint->y_position);
	
	result = (rSqd < colmatrRSqrd);
	
	return (result);
}


/*********************************************************************************
*
*		Name:			detBlocXInnerCylinder
*
*		Summary:		Test whether the segment will cross the collimator cylinder.
*
*		Arguments:
*			CylPosCylinderTy	*innerCyl		- The inner cylinder.
*			Lb2D_Point			*point1			- First point (outside collimator).
*			Lb2D_Point			*point2			- Second point of line segment.
*
*		Function return: True if crosses, false otherwise.
*
*********************************************************************************/

Boolean detBlocXInnerCylinder(	CylPosCylinderTy	*innerCyl, 
								Lb2D_Point			*point1, 
								Lb2D_Point			*point2)

{
	double			colmatrRadius;		/* Collimator radius */
	double			xLength;			/* Line segment x length */
	double			yLength;			/* Line segment y length */
	double			segLengthSqrd;		/* Length of point1<->point2 segment (squared) */
	double			segLength;			/* Length of point1<->point2 segment */
	PHG_Direction	segDir;				/* Direction of point2 from point1 */
	Boolean			result;				/* Returned function result */
	double			d1, d2;				/* Distances to collimator cylinder */
	
	
	/* ### RH will be writing new function for this */
	
	/* NOTE:  This should be a collimator module function */
	/* NOTE:  This is for circular cylinders only */
	
	/* Get the outer edge of the collimator (inner bound of the detectors) */
	colmatrRadius = innerCyl->radius;
	
	/* Calculate line segment length and direction */
	xLength = point2->x_position - point1->x_position;
	yLength = point2->y_position - point1->y_position;
	segLengthSqrd = PHGMATH_Square(xLength) + PHGMATH_Square(yLength);
	segLength = PHGMATH_SquareRoot(segLengthSqrd);
	segDir.cosine_x = xLength / segLength;
	segDir.cosine_y = yLength / segLength;
	
	result = false;
	if (CylPosFind2dIntersection(point1->x_position, point1->y_position, 
									segDir.cosine_x, segDir.cosine_y, colmatrRadius, 
									&d1, &d2)) {
		if (d2 > 0) {
			/* The radius-relative distance from a chord to its circle is ~ c^2/8R^2.
				Use the square root instead for speed. */
			if (PhgMathRealNumAreEqual(0.353553391*(d2-d1)/colmatrRadius, 0.0, -4, 0, 0, 0)) {
				/* Edge just or barely touches (tangent) */
				/* Allow this */
			}
			else {
				/* Part of ray is inside the collimator */
				if (d1 < segLength) {
					/* Part of line segment is inside the collimator */
					result = true;
				}
			}
		}
	}
	
	return (result);
}


/*********************************************************************************
*
*		Name:		detBlocProjectToCylinder
*
*		Summary:	Given a location and angle, project to cylinder.
*
*		Arguments:
*			PHG_Position		*positionPtr	- The current position.
*			PHG_Direction		*directionPtr	- The current direction.
*			CylPosCylinderTy	*cylinderPtr	- The cylinder to project to.
*			PHG_Position		*newPosPtr		- The new position.
*			double				*distPtr		- Distance to projection.
*
*		Function return: TRUE unless cosine_z == +/-1.
*
*********************************************************************************/

Boolean detBlocProjectToCylinder(	PHG_Position *positionPtr,
									PHG_Direction *directionPtr, CylPosCylinderTy *cylinderPtr,
									PHG_Position *newPosPtr, double *distPtr)

{
	/* CylPosProjectToCylinder refuses to calculate intersections for flat 
		(in z-direction) trajectories, so it has to be duplicated here 
		(but fixed) for detBlocHitInnerCylinder to work correctly. */
	
	
	Boolean	intersects = false;	/* Do we intersect */
	
	do {
		/* If cosine_z == +/- 1, we will never intersect */
		if (fabs(directionPtr->cosine_z) == 1.0) {
			break;
		}
		
		/* Compute distance to surface */
		CylPosCalcDistanceToCylSurface(positionPtr,
			directionPtr, cylinderPtr, distPtr);
		
		if (*distPtr > 0.0) {
			/* Calculate z intersection */
			newPosPtr->z_position = positionPtr->z_position +
				(*distPtr * directionPtr->cosine_z);
			
			/* Calculate x intersection */
			newPosPtr->x_position = positionPtr->x_position + 
				(*distPtr * directionPtr->cosine_x);
			
			/* Calculate y intersection */
			newPosPtr->y_position = positionPtr->y_position +
				(*distPtr * directionPtr->cosine_y);
		}
		
		intersects = true;
	} while (false);
	
	return (intersects);
}


/*********************************************************************************
*
*		Name:			detBlocHitInnerCylinder
*
*		Summary:		Test whether the photon will intersect the inner cylinder.
*						If so, return the intersection point and the distance to it.
*
*		Arguments:
*			CylPosCylinderTy	*innerCyl		- The inner cylinder.
*			PHG_TrackingPhoton	*photonPtr		- The photon.
*			PHG_Position		*intersectPt	- Position of the intersection.
*			double				*distance		- Distance to the intersection.
*
*		Function return: True if will intersect, false otherwise.
*
*********************************************************************************/

Boolean detBlocHitInnerCylinder(	CylPosCylinderTy	*innerCyl, 
									PHG_TrackingPhoton	*photonPtr, 
									PHG_Position		*intersectPt, 
									double				*distance)

{
	double				radProjection;		/* Used to determine photon radial direction */
	double				radDirection;		/* Photon radial direction indicator */
	Boolean				hitInnerCyl;		/* Whether photon will hit inner cylinder */
	double				innerLineCos;		/* Normal form photon line cosine */
	double				innerLineSin;		/* Normal form photon line sine */
	double				innerLineDist;		/* Normal form photon line distance */
	double				innerLineMinDist;	/* Closest distance of photon path to center */
	
	
	/* ### RH will be writing a new one of this */
	
	/* NOTE:  This is for circular cylinders only */
	
	/* NOTE:  This assumes the current photon position is outside of the inner cylinder */
	
	/* Determine the photon's radial direction */
	/* Use dot product of photon's x-y direction and a vector from the origin */
	/* +1 outward, -1 inward, 0 perp */
	radProjection = 
		photonPtr->angle.cosine_x * photonPtr->location.x_position  +  
		photonPtr->angle.cosine_y * photonPtr->location.y_position;
	radDirection = copysign(1.0, radProjection);
	
	if (radDirection >= 0.0) {
		/* Won't hit inner cylinder (if is circular) */
		hitInnerCyl = false;
	}
	else {
		/* Has some motion radially inward */
		
		/* Determine shortest distance on potential path to center of cylinders */
		/* Is just abs(distance) in normal form of line starting at current position 
			and extending through a projected point on the line */
		Lb2D_Point		projPt;			/* Projected point on photon line */
		
		projPt.x_position = photonPtr->location.x_position + 
							photonPtr->angle.cosine_x;
		projPt.y_position = photonPtr->location.y_position + 
							photonPtr->angle.cosine_y;
		Lb2DNormalLine((Lb2D_Point*)&(photonPtr->location), &projPt, 
						&innerLineCos, &innerLineSin, &innerLineDist);
		innerLineMinDist = fabs(innerLineDist);
		
		if (innerLineMinDist > innerCyl->radius /* Major elliptical axis */) {
			/* Won't hit inner cylinder */
			hitInnerCyl = false;
		}
		else if (innerLineMinDist <= innerCyl->radius /* Minor elliptical axis */) {
			/* Will hit inner cylinder */
			hitInnerCyl = true;
		}
		else {
			/* Must test for inner cylinder intersection */
			/* For circular cylinders will never occur */
			hitInnerCyl = false;
		}
	}
	
	if (hitInnerCyl) {
		/* Will hit inner cylinder */
		
		PHG_Position		closestInnerPt;		/* Closest point to center of photon path */
		double				photonXYTravDist;	/* 2D photon travel distance to closest pt */
		double				photonSinZ;			/* Sine (non-zero) of photon acos(cosine_z) */
		PHG_Position		projPosition;		/* New photon position */
		PHG_Direction		reverseAngle;		/* Reverse direction of photon */
		double				cylDistance;		/* Distance to inner cylinder */
		
		/* Find point on potential photon path closest to center axis */
		/* Use normal form of path computed above:
			x and y follow easily from 2D line normal form definition;
			z is just z-distance traveled on path projected to closest point, 
				added to current position
		*/
		closestInnerPt.x_position =  - innerLineDist * innerLineCos;
		closestInnerPt.y_position =  - innerLineDist * innerLineSin;
		photonXYTravDist = PHGMATH_RadialPos( 
			photonPtr->location.x_position - closestInnerPt.x_position, 
			photonPtr->location.y_position - closestInnerPt.y_position);
		photonSinZ = PHGMATH_SquareRoot(1.0 - PHGMATH_Square(photonPtr->angle.cosine_z));
		closestInnerPt.z_position = photonPtr->location.z_position  +  
			photonXYTravDist / photonSinZ * photonPtr->angle.cosine_z;
		
		/* Compute distance to first inner cylinder intersection */
		if (innerLineMinDist == innerCyl->radius) {	/* Implies circular cylinder */
			/* Just touches inner cylinder */
			projPosition = closestInnerPt;
		}
		else {
			/* Potential photon path passes through inner cylinder */
			/* Project from inner point back towards current point */
			/* Must intersect because is not on inner cylinder already */
			reverseAngle.cosine_x =  - photonPtr->angle.cosine_x;
			reverseAngle.cosine_y =  - photonPtr->angle.cosine_y;
			reverseAngle.cosine_z =  - photonPtr->angle.cosine_z;
			detBlocProjectToCylinder(&closestInnerPt, 
										&reverseAngle, 
										innerCyl, 
										&projPosition, 
										&cylDistance /* throwaway */ );
		}
		cylDistance = PHGMATH_SquareRoot( 
				PHGMATH_Square(photonPtr->location.x_position - 
									projPosition.x_position) + 
				PHGMATH_Square(photonPtr->location.y_position - 
									projPosition.y_position) + 
				PHGMATH_Square(photonPtr->location.z_position - 
									projPosition.z_position) );
		
		/* Set return values */
		*intersectPt = projPosition;
		*distance = cylDistance;
	}
	
	return (hitInnerCyl);
}


/*********************************************************************************
*
*		Name:			detBlocPtInOuterCylinder
*
*		Summary:		Test whether the point is inside the detector cylinder.
*
*		Arguments:
*			Lb2D_Point			*testPoint		- Position to be tested.
*
*		Function return: True if inside, false otherwise.
*
*********************************************************************************/

Boolean detBlocPtInOuterCylinder(Lb2D_Point		*testPoint)

{
	double			rSqd;				/* Point radial distance (squared) */
	Boolean			result;				/* Returned function result */
	
	
	/* Get the point distance (squared) */
	rSqd = PHGMATH_Square(testPoint->x_position) + PHGMATH_Square(testPoint->y_position);
	
	result = (rSqd < detBlocOutRadSqrd);
	
	return (result);
}


/*********************************************************************************
*
*		Name:			detBlocHitOuterCylinder
*
*		Summary:		Test whether the photon will intersect the outer cylinder.
*						If so, return the intersection point and the distance to it.
*
*		Arguments:
*			PHG_TrackingPhoton	*photonPtr		- The photon.
*			PHG_Position		*intersectPt	- Position of the intersection.
*			double				*distance		- Distance to the intersection.
*
*		Function return: True if will intersect, false otherwise.
*
*********************************************************************************/

Boolean detBlocHitOuterCylinder(	PHG_TrackingPhoton	*photonPtr, 
									PHG_Position		*intersectPt, 
									double				*distance)

{
	CylPosCylinderTy	outerCylinder;		/* Outer cylindrical boundary */
	Boolean				intersects;			/* Function return value */
	
	
	/* Exiting an elliptical outer cylinder and a circular outer cylinder are 
		equivalent as long as nothing is in the space between the two cylinder spaces */
	
	if (! detBlocPtInOuterCylinder((Lb2D_Point*)&(photonPtr->location))) {
		/* Photon is already outside outer cylinder */
		intersects = false;
	}
	else {
		/* Create an outer cylinder */
		outerCylinder.radius = detBlocOutRadius;
		outerCylinder.zMin = -10000.0;
		outerCylinder.zMax = 10000.0;
		outerCylinder.centerX = 0.0;
		outerCylinder.centerY = 0.0;
		
		intersects = CylPosProjectToCylinder(	&(photonPtr->location), 
												&(photonPtr->angle), 
												&outerCylinder, 
												intersectPt, 
												distance);
	}
	
	return (intersects);
}


/*********************************************************************************
*
*		Name:			detBlocCrossZoneBounds
*
*		Summary:		Test whether the photon will cross a zone boundary.
*						If so, return the intersection point and the distance to it.
*
*		Arguments:
*			PHG_TrackingPhoton	*photonPtr		- The photon.
*			LbFourByte			itsRing			- Current ring.
*			LbFourByte			itsZone			- Current ring zone.
*			PHG_Position		*intersectPt	- Position of the zone intersection.
*			double				*distance		- Distance to the zone intersection.
*
*		Function return: Type of zone crossing event, or null.
*
*********************************************************************************/

DetBlocEvents detBlocCrossZoneBounds(	PHG_TrackingPhoton	*thePhotonPtr, 
										LbFourByte			itsRing, 
										LbFourByte			itsZone, 
										PHG_Position		*intersectPt, 
										double				*distance)

{
	DetBlocEvents		photonEvent;		/* Result of the function */
	ArcRangeSpec		zoneArcRange;		/* Arc range of the current zone */
	double				photon2dRadius;		/* 2D distance of photon from center */
	double				photonPosCosX;		/* Dir cos (x) to photon position */
	double				photonPosCosY;		/* Dir cos (y) to photon position */
	double				dirCos2DNormalizer;	/* 2D "distance" of 3D x-y dir cosines */
	double				photonDirCosX;		/* Dir cos (x) of photon direction */
	double				photonDirCosY;		/* Dir cos (y) of photon direction */
	int					angDirection;		/* Photon relative angular direction indicator */
	double				photonSinZ;			/* Photon dir sin */
	double				photonPosTheta;		/* Photon position angle */
	double				photonDirTheta;		/* Photon direction angle */
	double				photonPosDirTheta;	/* Angle between photon position and direction */
	double				zoneBoundTheta;		/* Angle of zone boundary */
	double				zonePhotonTheta;	/* Angle between photon and zone boundary */
	double				zoneProjTheta;		/* Angle to photon/zone intersection */
	double				zoneProjXYDist;		/* 2D distance to photon/zone intersection */
	double				zoneProjDistance;	/* 3D distance to photon/zone intersection */
	PHG_Position		zoneBoundPos;		/* Photon/zone intersection */
	
	
	/* Assume default of no zone boundary crossing */
	photonEvent = detBlocEvt_Null;
	
	/* Get arc range of current photon zone */
	detBlocGetZoneRange(itsRing, itsZone, &zoneArcRange);
	
	/* Determine the photon's angular direction */
	/* Compare with vector from origin to photon position */
	/* +1 counterclockwise, -1 clockwise, 0 radial only */
	photon2dRadius = PHGMATH_RadialPos( 
		thePhotonPtr->location.x_position, thePhotonPtr->location.y_position);
	photonPosCosX = thePhotonPtr->location.x_position / photon2dRadius;
	photonPosCosY = thePhotonPtr->location.y_position / photon2dRadius;
	dirCos2DNormalizer = PHGMATH_RadialPos( 
		thePhotonPtr->angle.cosine_x, thePhotonPtr->angle.cosine_y);
	photonDirCosX = thePhotonPtr->angle.cosine_x / dirCos2DNormalizer;
	photonDirCosY = thePhotonPtr->angle.cosine_y / dirCos2DNormalizer;
	angDirection = Lb2DDirCosComp(photonDirCosX, photonDirCosY, photonPosCosX, photonPosCosY);
	
	/* Calculate photon sin(z) where cos(z)=thePhotonPtr->angle.cosine_z */
	photonSinZ = PHGMATH_SquareRoot( 
					1.0 - PHGMATH_Square(thePhotonPtr->angle.cosine_z));
	
	/* Check distance to zone boundaries */
	if (angDirection == 0) {
		/* Moving radially:  can't cross a zone boundary */
	}
	else {
		/* Use law of sines; need to compute triangle angles and sines */
		photonPosTheta = acos(photonPosCosX);
		if (photonPosCosY < 0.0) {
			photonPosTheta = PHGMATH_2PI - photonPosTheta;
		}
		photonDirTheta = acos(photonDirCosX);
		if (photonDirCosY < 0.0) {
			photonDirTheta = PHGMATH_2PI - photonDirTheta;
		}
		if (angDirection > 0) {
			/* Moving counterclockwise */
			if (photonDirTheta < photonPosTheta) {
				/* Lines split across x-axis */
				photonPosDirTheta = -PHGMATH_PI - (photonDirTheta - photonPosTheta);
			}
			else {
				photonPosDirTheta = PHGMATH_PI - (photonDirTheta - photonPosTheta);
			}
			zoneBoundTheta = acos(zoneArcRange.maxAngle.xDirCos);
			if (zoneArcRange.maxAngle.yDirCos < 0.0) {
				zoneBoundTheta = PHGMATH_2PI - zoneBoundTheta;
			}
			if (zoneBoundTheta < photonPosTheta) {
				/* Lines split across x-axis */
				zonePhotonTheta = PHGMATH_2PI + zoneBoundTheta - photonPosTheta;
			}
			else {
				zonePhotonTheta = zoneBoundTheta - photonPosTheta;
			}
		}
		else {
			/* Moving clockwise */
			if (photonPosTheta < photonDirTheta) {
				/* Lines split across x-axis */
				photonPosDirTheta = -PHGMATH_PI - (photonPosTheta - photonDirTheta);
			}
			else {
				photonPosDirTheta = PHGMATH_PI - (photonPosTheta - photonDirTheta);
			}
			zoneBoundTheta = acos(zoneArcRange.minAngle.xDirCos);
			if (zoneArcRange.minAngle.yDirCos < 0.0) {
				zoneBoundTheta = PHGMATH_2PI - zoneBoundTheta;
			}
			if (photonPosTheta < zoneBoundTheta) {
				/* Lines split across x-axis */
				zonePhotonTheta = PHGMATH_2PI + photonPosTheta - zoneBoundTheta;
			}
			else {
				zonePhotonTheta = photonPosTheta - zoneBoundTheta;
			}
		}
		zoneProjTheta = PHGMATH_PI - (zonePhotonTheta + photonPosDirTheta);
		if (zoneProjTheta <= 0.0) {
			/* Triangle cannot be formed; photon can't cross zone boundary */
		}
		else {
			zoneProjXYDist = sin(zonePhotonTheta) * 
								photon2dRadius / sin(zoneProjTheta);
			zoneProjDistance = zoneProjXYDist / photonSinZ;
			zoneBoundPos.x_position = thePhotonPtr->location.x_position  +  
				thePhotonPtr->angle.cosine_x * zoneProjDistance;
			zoneBoundPos.y_position = thePhotonPtr->location.y_position  +  
				thePhotonPtr->angle.cosine_y * zoneProjDistance;
			zoneBoundPos.z_position = thePhotonPtr->location.z_position  +  
				thePhotonPtr->angle.cosine_z * zoneProjDistance;
			
			if (angDirection > 0) {
				/* Hit the greater zone boundary */
				photonEvent = detBlocEvt_NextZone;
			}
			else {
				/* Hit the lesser zone boundary */
				photonEvent = detBlocEvt_PrevZone;
			}
		}
	}
	
	/* Return the results */
	if (photonEvent != detBlocEvt_Null) {
		*intersectPt = zoneBoundPos;
		*distance = zoneProjDistance;
	}
	return (photonEvent);
}


/*********************************************************************************
*
*		Name:			detBlocCalcDistanceToBlock
*
*		Summary:		Return the distance to the supplied block and its 
*							intersection point starting from the 
*							specified photon position and direction.
*
*		Arguments:
*			PHG_Position		*photonPosition		- Photon position.
*			PHG_Direction		*photonDirection	- Photon direction.
*			DetectorBlockPtr	theBlockPtr			- Pointer to the block.
*			double				*theDistance		- Returned block distance.
*			PHG_Position		*hitPoint			- Intersection with block.
*
*		Function return: Boolean: true if intersected, false if not.
*
*********************************************************************************/

Boolean detBlocCalcDistanceToBlock(	PHG_Position		*photonPosition, 
									PHG_Direction		*photonDirection, 
									DetectorBlockPtr	theBlockPtr, 
									double				*theDistance, 
									PHG_Position		*hitPoint)

{
	DetBlockTomoRingTy*	curRingPtr;		/* Ring info for the block's ring */
	double				ringMinus,		/* Lesser and greater detector Z-boundaries */
						ringPlus;
	double				minDistance;	/* Shortest distance to hit point */
	Lb2D_Rect			theBlockRect;	/* Block boundary 2D rectangle */
	Lb2D_PosDesc		phoBlockPos;	/* Relation of photon to block rectangle */
	
	
	/* Get the ring boundaries of the detector */
	curRingPtr = 
		&(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[theBlockPtr->itsRing]);
	ringMinus = curRingPtr->MinZ + curRingPtr->AxialShift;
	ringPlus = curRingPtr->MaxZ + curRingPtr->AxialShift;
	
	/* Verify that the photon and block are in the same ring */
	if ((photonPosition->z_position < ringMinus) || (ringPlus < photonPosition->z_position)) {
		/* Invalid parameters; indicated with a negative distance */
		minDistance = -2.0;
	}
	else {
		/* Get the given block's 2D rectangle boundary */
		theBlockRect = theBlockPtr->itsRect;
		
		/* Compare the photon's position to the block rectangle */
		phoBlockPos = Lb2DPointRectIntersect((Lb2D_Point*)photonPosition, &theBlockRect);
		
		switch (phoBlockPos) {
			case Lb2D_Inside:
			{
				/* Already inside block, so calculation is done */
				*hitPoint = *photonPosition;
				minDistance = 0.0;
			}
			break;
			
			
			case Lb2D_OnBound:
			{
				/* On edge of block; determine whether heading inward or outward */
				
				Lb2D_Point			point1;		/* Photon position as 2D point */
				Lb2D_Point			point2;		/* Other end of path as line segment */
				Lb2D_Point			*edgeCornerPtrs[5];/* Corners (wrapped) of theBlockRect */
				Boolean				inward;		/* Photon heading into the block or not */
				int					e;			/* Index through rect edges */
				Lb2D_PosDesc		sideIntersect;/* Relation of path line segment to 
													rect edge */
				
				/* Create a long line segment with the current point at one end */
				point1.x_position = photonPosition->x_position;
				point1.y_position = photonPosition->y_position;
				point2.x_position = point1.x_position + 8192.0;
				point2.y_position = point1.y_position + 8192.0;
				
				/* Set up a point array for looping through the block's sides */
				edgeCornerPtrs[0] = &(theBlockRect.corner_1);
				edgeCornerPtrs[1] = &(theBlockRect.corner_2);
				edgeCornerPtrs[2] = &(theBlockRect.corner_3);
				edgeCornerPtrs[3] = &(theBlockRect.corner_4);
				edgeCornerPtrs[4] = &(theBlockRect.corner_1);
				
				inward = true;			/* Intersection indicator */
				
				/* Check for intersection with each of the block's sides */
				for (e=0; e<4; e++) {
					sideIntersect = Lb2DLineSegsIntersect(&point1, &point2, 
										edgeCornerPtrs[e], edgeCornerPtrs[e+1]);
					if (sideIntersect == Lb2D_OnBound) {
						/* There was an intersection; 
							see if the photon passes across the block */
						double		edgeCos, 		/* Normal line descriptors of rect edge */
									edgeSin, 
									edgeDist;
						int			edgeOpp;		/* Index of opposite block edge */
						double		edgeCosOpp, 	/* Normal line descriptors of opposite edge */
									edgeSinOpp, 
									edgeDistOpp;
						double		halfDist;		/* Half of two edge dists */
						Lb2D_Point	point2Opp;		/* Other end of path as line segment */
						
						/* Get the normal lines of the intersecting and opposite block edges */
						Lb2DNormalLine(edgeCornerPtrs[e], edgeCornerPtrs[e+1], 
										&edgeCos, &edgeSin, &edgeDist);
						edgeOpp = (e+2) % 4;
						Lb2DNormalLine(edgeCornerPtrs[edgeOpp], edgeCornerPtrs[edgeOpp+1], 
										&edgeCosOpp, &edgeSinOpp, &edgeDistOpp);
						
						/* Project to a point that won't cross the opposite edge */
						halfDist = fabs((edgeDist - edgeDistOpp)/2);
						point2Opp.x_position = 
							point1.x_position + halfDist * photonDirection->cosine_x;
						point2Opp.y_position = 
							point1.y_position + halfDist * photonDirection->cosine_y;
						
						/* Check if the halfway point is between the edges */
						if (Lb2DPointParLinesIntersect(&point2Opp, 
								edgeCos, edgeSin, edgeDist, 
								edgeCosOpp, edgeSinOpp, edgeDistOpp) != Lb2D_Inside) {
							/* Photon can't enter the block */
							inward = false;
							break;
						}
					}
				}
				
				if (inward) {
					/* Photon is on the block boundary and headed inward; 
						already at block, so calculation is done */
					*hitPoint = *photonPosition;
					minDistance = 0.0;
				}
				else {
					/* Photon is on the block boundary but headed outward; 
						intersection cannot occur */
					minDistance = -3.0;
				}
			}
			break;
			
			
			case Lb2D_Outside:
			{
				Lb2D_Point			point1;		/* Photon position as 2D point */
				Lb2D_Point			point2;		/* Other end of path as line segment */
				double				segCos, 	/* Normal line descriptors of path segment */
									segSin, 
									segDist;
				Lb2D_Point			*edgeCornerPtrs[5];/* Corners (wrapped) of theBlockRect */
				int					e;			/* Index through rect edges */
				Lb2D_PosDesc		sideIntersect;/* Relation of path line segment to 
													rect edge */
				
				minDistance = -1.0;		/* Non-intersection indicator */
				
				/* Convert the photon's future path to a long line segment */
				point1.x_position = photonPosition->x_position;
				point1.y_position = photonPosition->y_position;
				point2.x_position = point1.x_position + 8192.0 * photonDirection->cosine_x;
				point2.y_position = point1.y_position + 8192.0 * photonDirection->cosine_y;
				if ((point1.x_position == point2.x_position) && 
						(point1.y_position == point2.y_position)) {
					/* Photon is moving axially, so can't possibly hit the block */
					break;
				}
				Lb2DNormalLine(&point1, &point2, &segCos, &segSin, &segDist);
				
				/* Set up a point array for looping through the block's sides */
				edgeCornerPtrs[0] = &(theBlockRect.corner_1);
				edgeCornerPtrs[1] = &(theBlockRect.corner_2);
				edgeCornerPtrs[2] = &(theBlockRect.corner_3);
				edgeCornerPtrs[3] = &(theBlockRect.corner_4);
				edgeCornerPtrs[4] = &(theBlockRect.corner_1);
				
				/* Check for intersection with each of the block's sides */
				for (e=0; e<4; e++) {
					sideIntersect = Lb2DLineSegsIntersect(&point1, &point2, 
										edgeCornerPtrs[e], edgeCornerPtrs[e+1]);
					if (sideIntersect != Lb2D_Outside) {
						/* There was an intersection; 
							find the intersection point and the distance to it */
						double		edgeCos, 		/* Normal line descriptors of rect edge */
									edgeSin, 
									edgeDist;
						double		det;			/* Determinant of non-parallel line coeffs */
						double		hitX,			/* Coordinates of possible hit point */
									hitY,
									hitZ;
						double		sideDistance2;	/* 2D distance^2 to hit point */
						double		sineZ2;			/* Sine^2 of photon Z-direction */
						double		hitDistance;	/* 3D distance to hit point */
						
						
						/* Calculate the 2D and 3D distances to the side */
						Lb2DNormalLine(edgeCornerPtrs[e], edgeCornerPtrs[e+1], 
										&edgeCos, &edgeSin, &edgeDist);
						det = segCos*edgeSin - edgeCos*segSin;
						if (fabs(det) < 1E-15) {
							/* Consider the lines to be overlapping */
							/* But this means the photon never enters the block */
							hitDistance = -1.0;
						}
						else {
							hitX = (segSin*edgeDist - edgeSin*segDist) / det;
							hitY = (edgeCos*segDist - segCos*edgeDist) / det;
							if (sideIntersect == Lb2D_OnBound) {
								/* The photon crosses through a corner */
								/* Move the hit point exactly to the corner; 
									use 1-norm to find closest corner */
								if (((hitX-edgeCornerPtrs[e]->x_position) + 
										(hitY-edgeCornerPtrs[e]->y_position)) < 
									((hitX-edgeCornerPtrs[e+1]->x_position) + 
										(hitY-edgeCornerPtrs[e+1]->y_position))) {
									/* Closest to [e] corner */
									hitX = edgeCornerPtrs[e]->x_position;
									hitY = edgeCornerPtrs[e]->y_position;
								}
								else {
									/* Closest to [e+1] corner */
									hitX = edgeCornerPtrs[e+1]->x_position;
									hitY = edgeCornerPtrs[e+1]->y_position;
								}
							}
							sideDistance2 = PHGMATH_Square(hitX - point1.x_position) + 
											PHGMATH_Square(hitY - point1.y_position);
							sineZ2 = 									/* Won't be 0 */
								(1.0 - PHGMATH_Square(photonDirection->cosine_z));
							hitDistance = PHGMATH_SquareRoot(sideDistance2 / sineZ2);
							
							/* Determine the z position for this hit point */
							hitZ = photonPosition->z_position + 
								hitDistance * photonDirection->cosine_z;
							
							/* Compare the hit point to the ring boundaries */
							if (hitZ < ringMinus) {
								/* Photon first exited the lesser side of the ring; no hit */
								hitDistance = -1.0;
							}
							else if (hitZ > ringPlus) {
								/* Photon first exited the greater side of the ring; no hit */
								hitDistance = -1.0;
							}
							else {
								/* Photon hit inside the ring; a hit point has been found */
							}
						}
						
						/* If a valid minimum distance... */
						if (hitDistance >= 0.0) {
							/* Check for it being a shorter one */
							if ((minDistance < 0.0) || (hitDistance < minDistance)) {
								/* Shorter distance found */
								minDistance = hitDistance;
								hitPoint->x_position = hitX;
								hitPoint->y_position = hitY;
								hitPoint->z_position = hitZ;
							}
						}
					}
				}
			}
			break;
			
			
			default:
				/* Shouldn't happen */
			break;
		}
	}
	
	*theDistance = minDistance;
	return (minDistance >= 0.0);
}


/*********************************************************************************
*
*		Name:			detBlocMakeRBDatabase
*
*		Summary:		Make the (r,b) database of the block detector blocks.
*
*		Arguments:		None.
*
*		Function return: True if successful, false if memory failure.
*
*********************************************************************************/

Boolean detBlocMakeRBDatabase(void)

{
	LbUsFourByte		numRings;			/* Number of rings */
	LbUsFourByte		maxRingBlocks;		/* Maximum number of blocks in a ring */
	LbUsFourByte		r;					/* Index through rings */
	DetBlockTomoRingTy	*curRingPtr;		/* Current ring info */
	LbUsFourByte		curBlocks;			/* Number of blocks in current ring */
	DetectorBlockPtr	curBlockPtr;		/* Current database block */
	LbUsFourByte		b;					/* Index through ring blocks */
	DetBlockTomoBlockTy	*curBlockInfoPtr;	/* Current block info */
	PHG_Position		blockPos;			/* Rect corner in block coordinates */
	PHG_Position		tomoPos;			/* Rect corner in tomo coordinates */
	
	
	/* Count the blocks in the detector */
	numRings = DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings;
	maxRingBlocks = 0;
	for (r=0; r<numRings; r++) {
		curRingPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[r]);
		curBlocks = curRingPtr->NumBlocks;
		if (curBlocks > maxRingBlocks) {
			maxRingBlocks = curBlocks;
		}
	}
	
	/* Record the number of rings and maximum blocks per ring */
	detBlocNumRings = numRings;
	detBlocMaxRingBlocks = maxRingBlocks;
	
	/* Allocate memory for the (r,b) blocks database */
	detBlocBlocksDatabase = ( DetectorBlockPtr )LbMmAlloc( 
								numRings * maxRingBlocks * sizeof(DetectorBlock));
	if (!detBlocBlocksDatabase) {
		/* Couldn't allocate the memory */
		return(false);
	}
		
	/* Fill the database */
	curBlockPtr = detBlocBlocksDatabase;
	for (r=0; r<numRings; r++) {
		curRingPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[r]);
		curBlocks = curRingPtr->NumBlocks;
		for (b=0; b<curBlocks; b++) {
			curBlockPtr->itsRing = r;
			curBlockPtr->itsBlock = b;
			
			/* Compute the four corners of itsRect */
			{
				curBlockInfoPtr = &(curRingPtr->BlockInfo[b]);
				blockPos.x_position = curBlockInfoPtr->XMin;
				blockPos.y_position = curBlockInfoPtr->YMin;
				blockPos.z_position = curBlockInfoPtr->ZMin;
				detBlocBlockToTomoCoord(curBlockInfoPtr, &blockPos, &tomoPos);
				curBlockPtr->itsRect.corner_1.x_position = tomoPos.x_position;
				curBlockPtr->itsRect.corner_1.y_position = tomoPos.y_position;
				blockPos.x_position = curBlockInfoPtr->XMax;
				detBlocBlockToTomoCoord(curBlockInfoPtr, &blockPos, &tomoPos);
				curBlockPtr->itsRect.corner_2.x_position = tomoPos.x_position;
				curBlockPtr->itsRect.corner_2.y_position = tomoPos.y_position;
				blockPos.y_position = curBlockInfoPtr->YMax;
				detBlocBlockToTomoCoord(curBlockInfoPtr, &blockPos, &tomoPos);
				curBlockPtr->itsRect.corner_3.x_position = tomoPos.x_position;
				curBlockPtr->itsRect.corner_3.y_position = tomoPos.y_position;
				blockPos.x_position = curBlockInfoPtr->XMin;
				detBlocBlockToTomoCoord(curBlockInfoPtr, &blockPos, &tomoPos);
				curBlockPtr->itsRect.corner_4.x_position = tomoPos.x_position;
				curBlockPtr->itsRect.corner_4.y_position = tomoPos.y_position;
			}
			
			/* Set the z bounds */
			{
				curBlockPtr->zMin = tomoPos.z_position;
				blockPos.z_position = curBlockInfoPtr->ZMax;
				detBlocBlockToTomoCoord(curBlockInfoPtr, &blockPos, &tomoPos);
				curBlockPtr->zMax = tomoPos.z_position;
			}
			
			/* Compute and record the subtended arcs for each block rectangle */
			/* Since rectangles are finite and at positive distance from the origin, they will 
				subtend less than pi radians of arc (so Lb2DDirCosComp results are valid) */
			{
				Lb2D_Point		origin;				/* (0,0) */
				DetAngleSpec	curRectMinAngle;	/* Rect minimum angle of subtended arc */
				DetAngleSpec	curRectMaxAngle;	/* Rect maximum angle of subtended arc */
				double			xCos;				/* Direction cosines for rect corners */
				double			yCos;
				
				
				/* Initialize origin point */
				origin.x_position = 0.0;
				origin.y_position = 0.0;
				
				/* Compute the subtended arc (using the four rectangle corners) */
				Lb2DDirCosines(&origin, &(curBlockPtr->itsRect.corner_1), 
								&(curRectMinAngle.xDirCos), 
								&(curRectMinAngle.yDirCos));
				curRectMaxAngle = curRectMinAngle;
				
				Lb2DDirCosines(&origin, &(curBlockPtr->itsRect.corner_2), &xCos, &yCos);
				if (Lb2DDirCosComp(xCos, yCos, 
							curRectMinAngle.xDirCos, curRectMinAngle.yDirCos) < 0) {
					curRectMinAngle.xDirCos = xCos;
					curRectMinAngle.yDirCos = yCos;
				}
				else if (Lb2DDirCosComp(xCos, yCos, 
							curRectMaxAngle.xDirCos, curRectMaxAngle.yDirCos) > 0) {
					curRectMaxAngle.xDirCos = xCos;
					curRectMaxAngle.yDirCos = yCos;
				}
				
				Lb2DDirCosines(&origin, &(curBlockPtr->itsRect.corner_3), &xCos, &yCos);
				if (Lb2DDirCosComp(xCos, yCos, 
							curRectMinAngle.xDirCos, curRectMinAngle.yDirCos) < 0) {
					curRectMinAngle.xDirCos = xCos;
					curRectMinAngle.yDirCos = yCos;
				}
				else if (Lb2DDirCosComp(xCos, yCos, 
							curRectMaxAngle.xDirCos, curRectMaxAngle.yDirCos) > 0) {
					curRectMaxAngle.xDirCos = xCos;
					curRectMaxAngle.yDirCos = yCos;
				}
				
				Lb2DDirCosines(&origin, &(curBlockPtr->itsRect.corner_4), &xCos, &yCos);
				if (Lb2DDirCosComp(xCos, yCos, 
							curRectMinAngle.xDirCos, curRectMinAngle.yDirCos) < 0) {
					curRectMinAngle.xDirCos = xCos;
					curRectMinAngle.yDirCos = yCos;
				}
				else if (Lb2DDirCosComp(xCos, yCos, 
							curRectMaxAngle.xDirCos, curRectMaxAngle.yDirCos) > 0) {
					curRectMaxAngle.xDirCos = xCos;
					curRectMaxAngle.yDirCos = yCos;
				}
				curBlockPtr->minAngle = curRectMinAngle;
				curBlockPtr->maxAngle = curRectMaxAngle;
			}
			
			curBlockPtr->userBlockData = curBlockInfoPtr;
			
			curBlockPtr++;
		}
		
		/* Mark any remaining ring blocks as unused */
		for (b=curBlocks; b<maxRingBlocks; b++) {
			curBlockPtr->userBlockData = NULL;
			curBlockPtr++;
		}
	}
	
	return(true);
}


/*********************************************************************************
*
*		Name:			DetBlocValidateBlocks
*
*		Summary:		Validate that the user-specified block positions are 
*							self-consistent and do not overlap each other or 
*							the collimator cylinder.
*						Also initialize some global variables.
*
*		Arguments:		None.
*
*		Function return: Boolean; true if all blocks are valid.
*
*********************************************************************************/

Boolean DetBlocValidateBlocks(void)

{
	Boolean				blocksValid;		/* Flag for all blocks okay */
	
	
	blocksValid = true;
	
	/* Create the (r,b) blocks database */
	blocksValid = detBlocMakeRBDatabase();
	
	if (blocksValid) {
		/* Get info on inner cylinder */
		detBlocGetInnerCylinder(&detBlocInnerCylinder);
	}
	
	if (blocksValid) {
		/* Test the blocks for self-consistency and boundary conditions */
		/* Also find and record the outer block cylinder radius */
		blocksValid = detBlocSelfConsCheck();
	}
	
	if (blocksValid) {
		/* Test the blocks for overlap */
		blocksValid = detBlocOverlapCheck();
		/*LbInPrintf( "\n\nWARNING!!!  Block overlap check not being executed!\n");*/
	}
	
	return (blocksValid);
}


/*********************************************************************************
*
*		Name:			detBlocSelfConsCheck
*
*		Summary:		Validate the user-specified block positions 
*						  for self-consistency and boundary conditions.
*						Also find and record the outer block cylinder radius.
*
*		Arguments:		None.
*
*		Function return: Boolean; true if all blocks are consistent and valid.
*
*********************************************************************************/

Boolean detBlocSelfConsCheck(void)

{
	Boolean				blocksValid;		/* Flag for all blocks okay */
	double				maxDistSqrd;		/* Outer radius of detector (^2) */
	DetectorBlockPtr	curRingBlockPtr;	/* Pointer to start of (r,b) ring blocks */
	DetectorBlockPtr	curBlockPtr;		/* Pointer to current (r,b) ring block */
	DetBlockTomoRingTy	*curRingPtr;		/* Current ring info */
	CylPosCylinderTy	inRingCyl;			/* Inner boundary cylinder for current ring */
	CylPosCylinderTy	outRingCyl;			/* Outer boundary cylinder for current ring */
	double				ringMinZ;			/* Minimum axial boundary of current ring */
	double				ringMaxZ;			/* Maximum axial boundary of current ring */
	LbFourByte			r;					/* Index through rings of current ring */
	LbFourByte			b;					/* Index through block in current ring */
	Lb2D_Rect			curRect;			/* Current rectangle being validated */
	double				minCornerDistSqrd;	/* Distance^2 of closest corner of curRect */
	double				maxCornerDistSqrd;	/* Distance^2 of farthest corner of curRect */
	int					minCorner;			/* Index of closest corner of curRect */
	int					maxCorner;			/* Index of farthest corner of curRect */
	double				distSqrd;			/* Working distance^2 of current corner */
	Lb2D_Point			point1;				/* Closest corner of curRect */
	int					adjCorner;			/* Index of corner adjacent to minCorner */
	Lb2D_Point			point2;				/* Corner of curRect adjacent to point1 */
	
	
	/* Initializations */
	blocksValid = true;
	maxDistSqrd = 0.0;
	
	/* Proceed through the user blocks, ring by ring */
	for (r=0; r<detBlocNumRings; r++) {
		curRingBlockPtr = &(detBlocBlocksDatabase[r*detBlocMaxRingBlocks]);
		curBlockPtr = curRingBlockPtr;
		
		/* Get boundaries for this ring */
		curRingPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[r]);
		inRingCyl = detBlocInnerCylinder;
		inRingCyl.radius = PHGMATH_Min(curRingPtr->XInnerRadius, curRingPtr->YInnerRadius);
		outRingCyl = detBlocInnerCylinder;
		outRingCyl.radius = PHGMATH_Max(curRingPtr->XOuterRadius, curRingPtr->YOuterRadius);
		ringMinZ = curRingPtr->MinZ;
		ringMaxZ = curRingPtr->MaxZ;
		
		for (b=0; b<detBlocMaxRingBlocks; b++) {
			if (curBlockPtr->userBlockData) {
				/* Check that the block is within its axial bounds */
				if ((curBlockPtr->userBlockData->ZMin < ringMinZ) || 
						(curBlockPtr->userBlockData->ZMax > ringMaxZ)) {
					/* Part of block extends beyond its axial bounds */
					blocksValid = false;
					sprintf(detBlocErrStr, 
								"Block %ld in ring %ld is outside its axial bounds.", 
								(long)b, (long)r);
					PhgAbort(detBlocErrStr, false);
				}
				
				/* Get the block rectangle */
				curRect = curBlockPtr->itsRect;
				
				/* Find the innermost and outermost corners of the rectangle; 
					in case there are two (there will be no more than two), 
					either one is fine */
				/* NOTE:  Assumes center is at (0,0) */
				minCornerDistSqrd = PHGMATH_Square(curRect.corner_1.x_position) + 
									PHGMATH_Square(curRect.corner_1.y_position);
				maxCornerDistSqrd = minCornerDistSqrd;
				minCorner = 1;
				maxCorner = 1;
				distSqrd = PHGMATH_Square(curRect.corner_2.x_position) + 
							PHGMATH_Square(curRect.corner_2.y_position);
				if (distSqrd < minCornerDistSqrd) {
					minCornerDistSqrd = distSqrd;
					minCorner = 2;
				}
				else if (distSqrd > maxCornerDistSqrd) {
					maxCornerDistSqrd = distSqrd;
					maxCorner = 2;
				}
				distSqrd = PHGMATH_Square(curRect.corner_3.x_position) + 
							PHGMATH_Square(curRect.corner_3.y_position);
				if (distSqrd < minCornerDistSqrd) {
					minCornerDistSqrd = distSqrd;
					minCorner = 3;
				}
				else if (distSqrd > maxCornerDistSqrd) {
					maxCornerDistSqrd = distSqrd;
					maxCorner = 3;
				}
				distSqrd = PHGMATH_Square(curRect.corner_4.x_position) + 
							PHGMATH_Square(curRect.corner_4.y_position);
				if (distSqrd < minCornerDistSqrd) {
					minCornerDistSqrd = distSqrd;
					minCorner = 4;
				}
				else if (distSqrd > maxCornerDistSqrd) {
					maxCornerDistSqrd = distSqrd;
					maxCorner = 4;
				}
				
				/* Check that each rectangle corner is outside of its inner ring 
					boundary and the collimator boundary
					(just need to check the minimum corner) */
				switch (minCorner) {
					case 1:		point1 = curRect.corner_1;  break;
					case 2:		point1 = curRect.corner_2;  break;
					case 3:		point1 = curRect.corner_3;  break;
					case 4:		point1 = curRect.corner_4;  break;
					default:	blocksValid = false;  break;
				}
				if (detBlocPtInInnerCylinder(&inRingCyl, &point1)) {
					/* Block is inside the inner ring boundary -- Error!! */
					sprintf(detBlocErrStr, 
								"Block %ld in ring %ld has a corner inside its inner bound.", 
								(long)b, (long)r);
					PhgAbort(detBlocErrStr, false);
				}
				if (detBlocPtInInnerCylinder(&detBlocInnerCylinder, &point1)) {
					/* Block is inside the collimator -- Error!! */
					sprintf(detBlocErrStr, 
								"Block %ld in ring %ld has a corner inside the collimator.", 
								(long)b, (long)r);
					PhgAbort(detBlocErrStr, false);
				}
				
				if (blocksValid) {
					/* Verify that no interior point of the rectangle is inside the 
						same boundaries, as well */
					
					/* First adjacent edge */
					adjCorner = (minCorner % 4) + 1;
					switch (adjCorner) {
						case 1:		point2 = curRect.corner_1;  break;
						case 2:		point2 = curRect.corner_2;  break;
						case 3:		point2 = curRect.corner_3;  break;
						case 4:		point2 = curRect.corner_4;  break;
						default:	blocksValid = false;  break;
					}
					if (detBlocXInnerCylinder(&inRingCyl, &point1, &point2)) {
						/* Part of block is inside the inner ring bound */
						blocksValid = false;
						sprintf(detBlocErrStr, 
									"Block %ld in ring %ld has leading edge inside its inner bound.", 
									(long)b, (long)r);
						PhgAbort(detBlocErrStr, false);
					}
					if (detBlocXInnerCylinder(&detBlocInnerCylinder, &point1, &point2)) {
						/* Part of block is inside the collimator */
						blocksValid = false;
						sprintf(detBlocErrStr, 
									"Block %ld in ring %ld has leading edge inside the collimator.", 
									(long)b, (long)r);
						PhgAbort(detBlocErrStr, false);
					}
					
					/* Second adjacent edge */
					if (blocksValid) {
						adjCorner = minCorner - 1;
						if (adjCorner == 0)
							adjCorner = 4;
						switch (adjCorner) {
							case 1:		point2 = curRect.corner_1;  break;
							case 2:		point2 = curRect.corner_2;  break;
							case 3:		point2 = curRect.corner_3;  break;
							case 4:		point2 = curRect.corner_4;  break;
							default:	blocksValid = false;  break;
						}
						if (detBlocXInnerCylinder(&inRingCyl, &point1, &point2)) {
							/* Part of block is inside the inner ring bound */
							blocksValid = false;
							sprintf(detBlocErrStr, 
										"Block %ld in ring %ld has trailing edge inside its inner bound.", 
										(long)b, (long)r);
							PhgAbort(detBlocErrStr, false);
						}
						if (detBlocXInnerCylinder(&detBlocInnerCylinder, &point1, &point2)) {
							/* Part of block is inside the collimator */
							blocksValid = false;
							sprintf(detBlocErrStr, 
										"Block %ld in ring %ld has trailing edge inside the collimator.", 
										(long)b, (long)r);
							PhgAbort(detBlocErrStr, false);
						}
					}
				}
				
				if (blocksValid) {
					/* Check that each rectangle corner is inside of its outer ring 
						boundary
						(just need to check the maximum corner) */
					switch (maxCorner) {
						case 1:		point1 = curRect.corner_1;  break;
						case 2:		point1 = curRect.corner_2;  break;
						case 3:		point1 = curRect.corner_3;  break;
						case 4:		point1 = curRect.corner_4;  break;
						default:	blocksValid = false;  break;
					}
					if (! detBlocPtInInnerCylinder(&outRingCyl, &point1)) {
						/* Block is outside the outer ring boundary -- Error!! */
						sprintf(detBlocErrStr, 
									"Block %ld in ring %ld has a corner outside its outer bound.", 
									(long)b, (long)r);
						PhgAbort(detBlocErrStr, false);
					}
				}
				
				if (blocksValid) {
					/* Update the outer radius^2 */
					if (maxCornerDistSqrd > maxDistSqrd) {
						maxDistSqrd = maxCornerDistSqrd;
					}
				}
				
				if (! blocksValid) {
					/* Something went wrong; time to quit */
					sprintf(detBlocErrStr, 
						"An unknown error occurred while validating detector blocks in ring %ld.", 
						(long)r);
					PhgAbort(detBlocErrStr, false);
				}
			}
			else {
				/* No more valid blocks in this ring */
			}
			
			/* Get next block */
			curBlockPtr++;
		}
	}
	
	/* Set the radii of the block detector area (guaranteed within) */
	detBlocOutRadSqrd = maxDistSqrd + 1.0;
	detBlocOutRadius = PHGMATH_SquareRoot(detBlocOutRadSqrd);
	
	return (blocksValid);
}


/*********************************************************************************
*
*		Name:			detBlocOverlapCheck
*
*		Summary:		Validate that the user-specified block positions 
*						  do not overlap each other.
*
*		Arguments:		None.
*
*		Function return: Boolean; true if all blocks are valid.
*
*********************************************************************************/

Boolean detBlocOverlapCheck(void)

{
	Boolean				blocksValid;		/* Flag for all blocks okay */
	DetectorBlockPtr	curRingBlockPtr;	/* Pointer to start of (r,b) ring blocks */
	DetectorBlockPtr	curBlockPtr;		/* Pointer to current (r,b) ring block */
	LbFourByte			r;					/* Index through rings of current ring */
	LbFourByte			b;					/* Index through block in current ring */
	Lb2D_Rect			curRect;			/* Current rectangle being validated */
	LbFourByte			i;					/* Index through validated blocks */
	Lb2D_PosDesc		intersectPos;		/* Rectangle intersection result */
	
	
	/* Proceed through the user blocks, ring by ring */
	blocksValid = true;
	for (r=0; r<detBlocNumRings; r++) {
		curRingBlockPtr = &(detBlocBlocksDatabase[r*detBlocMaxRingBlocks]);
		curBlockPtr = curRingBlockPtr;
		for (b=0; b<detBlocMaxRingBlocks; b++) {
			if (curBlockPtr->userBlockData) {
				/* Get the block rectangle */
				curRect = curBlockPtr->itsRect;
				
				/* Verify that the rectangle does not intersect any previous one */
				for (i=0; i<b; i++) {
					intersectPos = Lb2DRectsIntersect(
										&curRect, &(curRingBlockPtr[i].itsRect));
					if (intersectPos == Lb2D_Inside) {
						/* Interior rectangle intersection -- Error!! */
						sprintf(detBlocErrStr, 
							"Block %ld intersects with block %ld in ring %ld.", 
							(long)b, (long)i, (long)r);
						PhgAbort(detBlocErrStr, false);
					}
				}
				
				if (! blocksValid) {
					/* Something went wrong; time to quit */
					sprintf(detBlocErrStr, 
						"An unknown error occurred while checking detector blocks in ring %ld.", 
						(long)r);
					PhgAbort(detBlocErrStr, false);
				}
			}
			else {
				/* No more valid blocks in this ring */
			}
			
			/* Get next block */
			curBlockPtr++;
		}
	}
	
	return (blocksValid);
}


/*********************************************************************************
*
*		Name:			DetBlocDivideZones
*
*		Summary:		Divide the ring spaces into zones (cheesecake slices) 
*							so no single zone contains too many detector blocks.
*						All rings use the same zone boundaries.
*						Create the database of blocks within each zone;
*							blocks are indexed in order, starting with 0.
*						WARNING:  Blocks can be in multiple zones.
*
*		Arguments:
*			LbUsFourByte	maxZBlocks		- Maximum allowable blocks in a zone.
*
*		Function return: Number of zones; 0 if memory failure.
*
*********************************************************************************/

LbUsFourByte DetBlocDivideZones(LbUsFourByte		maxZBlocks)

{
	#define kMaxZones 300					/* Maximum number of zones per ring */
	
	DetAngleSpec		zoneBounds[kMaxZones];/* Zone boundaries (lesser edges) */
	LbFourByte			numZones;			/* Current number of ring zones */
	Boolean				done;				/* Flag for loop completion */
	Boolean				exceedsMaxBlocks;	/* maxZBlocks exceeded to keep zones manageable */
	LbUsFourByte		maxZoneBlocks;		/* Actual maximum of blocks per zone */
	Boolean				overcrowded;		/* Flag indicating too many blocks per zone */
	DetectorBlockPtr	curBlockPtr;		/* The current blocks database block */
	LbFourByte			r;					/* Index through rings */
	LbFourByte			z;					/* Index through zones */
	LbUsFourByte		zoneCounts[kMaxZones];/* Blocks per zone */
	LbFourByte			b;					/* Index through ring blocks */
	DetAngleSpec		zoneBoundLow;		/* Lesser angle boundary of current zone */
	DetAngleSpec		zoneBoundHigh;		/* Greater angle boundary of current zone */
	LbFourByte			numIndices;			/* Size of index database */
	DetectorBlockPtr*	blockPtrPtr;		/* Blocks pointer into detBlocIndexDatabase */
	LbUsFourByte		ringIndex;			/* Index to ring start of index database */
	
	
	/* Initialize the zone boundary list to the 4 quadrants */
	zoneBounds[0].xDirCos = 1.0;
	zoneBounds[0].yDirCos = 0.0;
	zoneBounds[1].xDirCos = 0.0;
	zoneBounds[1].yDirCos = 1.0;
	zoneBounds[2].xDirCos = -1.0;
	zoneBounds[2].yDirCos = 0.0;
	zoneBounds[3].xDirCos = 0.0;
	zoneBounds[3].yDirCos = -1.0;
	numZones = 4;
	
	/* Check the zones for overcrowding, subdividing if so, until everything fits */
	done = false;
	exceedsMaxBlocks = false;
	while (! done) {
		/* Proceed through the user blocks, ring by ring, checking for overcrowded zones */
		maxZoneBlocks = 0;
		overcrowded = false;
		curBlockPtr = detBlocBlocksDatabase;
		for (r=0; r<detBlocNumRings; r++) {
			for (z=0; z<numZones; z++) {
				zoneCounts[z] = 0;
			}
			for (b=0; b<detBlocMaxRingBlocks; b++) {
				if (curBlockPtr->userBlockData) {
					/* Check each zone for intersection with the block */
					for (z=0; z<numZones; z++) {
						zoneBoundLow = zoneBounds[z];
						if ((z+1) < numZones) {
							zoneBoundHigh = zoneBounds[z+1];
						}
						else {
							zoneBoundHigh = zoneBounds[0];
						}
						if ((Lb2DDirCosComp(curBlockPtr->maxAngle.xDirCos, 
											curBlockPtr->maxAngle.yDirCos, 
											zoneBoundLow.xDirCos, 
											zoneBoundLow.yDirCos) > 0)  && 
							(Lb2DDirCosComp(curBlockPtr->minAngle.xDirCos, 
											curBlockPtr->minAngle.yDirCos, 
											zoneBoundHigh.xDirCos, 
											zoneBoundHigh.yDirCos) < 0)) {
							/* Block intersects with zone */
							zoneCounts[z]++;
							if (zoneCounts[z] > maxZoneBlocks) {
								maxZoneBlocks = zoneCounts[z];
							}
							
							if (zoneCounts[z] > maxZBlocks) {
								/* Too many blocks in a zone */
								if (exceedsMaxBlocks) {
									/* But can't create more zones, so continue on */
								}
								else {
									overcrowded = true;
									break;
								}
							}
						}
					}
				}
				curBlockPtr++;
				
				if (overcrowded) {
					/* Stop checking the ring blocks */
					break;
				}
			}
			
			if (overcrowded) {
				/* Stop checking the rings */
				break;
			}
		}
		
		if (overcrowded) {
			if ((numZones*2) > kMaxZones) {
				/* Splitting further would create too many zones */
				exceedsMaxBlocks = true;
			}
			else {
				/* Binarily split the zones and try again */
				double			thetaLow;			/* Angle of zoneBoundLow */
				double			thetaHigh;			/* Angle of zoneBoundHigh */
				double			thetaMid;			/* Average of thetaLow and thetaHigh */
				
				for (z=numZones-1; z>=0; z--) {
					/* Calculate the intermediate zone boundary */
					zoneBoundLow = zoneBounds[z];
					if ((z+1) < numZones) {
						zoneBoundHigh = zoneBounds[z+1];
					}
					else {
						zoneBoundHigh = zoneBounds[0];
					}
					thetaLow = acos(zoneBoundLow.xDirCos);
					if (zoneBoundLow.yDirCos < 0.0) {
						thetaLow = PHGMATH_2PI - thetaLow;
					}
					thetaHigh = acos(zoneBoundHigh.xDirCos);
					if (zoneBoundHigh.yDirCos <= 0.0) {
						thetaHigh = PHGMATH_2PI - thetaHigh;
					}
					thetaMid = (thetaLow + thetaHigh) / 2.0;
					
					/* Shift up the current zone */
					zoneBounds[2*z] = zoneBounds[z];
					
					/* Set the intermediate zone boundary */
					zoneBounds[2*z+1].xDirCos = cos(thetaMid);
					zoneBounds[2*z+1].yDirCos = sin(thetaMid);
				}
				numZones *= 2;
			}
		}
		else {
			/* Checked all blocks without any zone exceeding the limit */
			done = true;
		}
	}
	
	/* Record the final number of zones */
	detBlocNumRingZones = numZones;
	
	/* Allocate memory for the ring zones list */
	detBlocZoneBounds = ( DetAngleSpec* )LbMmAlloc(numZones * sizeof(DetAngleSpec));
	if (!detBlocZoneBounds) {
		/* Couldn't allocate the memory */
		return(0);
	}
	
	/* Copy zoneBounds to the global ring zones list */
	for (z=0; z<numZones; z++) {
		detBlocZoneBounds[z] = zoneBounds[z];
	}
	
	/* Record the maximum block count per zone */
	detBlocNumZoneBlocks = maxZoneBlocks;
	
	/* Allocate memory for the block index database */
	numIndices = detBlocNumRings * numZones * maxZoneBlocks;
	detBlocIndexDatabase = ( DetectorBlockPtr* )LbMmAlloc(
										numIndices * sizeof(DetectorBlockPtr));
	if (!detBlocIndexDatabase) {
		/* Couldn't allocate the memory */
		LbMmFree((void **)&detBlocZoneBounds);
		return(0);
	}
	
	/* Initialize every block index to unused (= NULL) */
	blockPtrPtr = detBlocIndexDatabase;
	for (b=0; b<numIndices; b++) {
		*blockPtrPtr = NULL;
		blockPtrPtr++;
	}
	
	/* Proceed through the user blocks, ring by ring, adding each block to the database */
	curBlockPtr = detBlocBlocksDatabase;
	for (r=0; r<detBlocNumRings; r++) {
		ringIndex = r * numZones * maxZoneBlocks;
		for (z=0; z<numZones; z++) {
			zoneCounts[z] = 0;
		}
		for (b=0; b<detBlocMaxRingBlocks; b++) {
			if (curBlockPtr->userBlockData) {
				/* Check each ring zone for intersection with the block */
				for (z=0; z<numZones; z++) {
					zoneBoundLow = zoneBounds[z];
					if ((z+1) < numZones) {
						zoneBoundHigh = zoneBounds[z+1];
					}
					else {
						zoneBoundHigh = zoneBounds[0];
					}
					if ((Lb2DDirCosComp(curBlockPtr->maxAngle.xDirCos, 
										curBlockPtr->maxAngle.yDirCos, 
										zoneBoundLow.xDirCos, 
										zoneBoundLow.yDirCos) > 0)  && 
						(Lb2DDirCosComp(curBlockPtr->minAngle.xDirCos, 
										curBlockPtr->minAngle.yDirCos, 
										zoneBoundHigh.xDirCos, 
										zoneBoundHigh.yDirCos) < 0)) {
						
						/* Block intersects with zone; add it to the (r,z,i) index database */
						detBlocIndexDatabase[ringIndex + z*maxZoneBlocks + zoneCounts[z]] = 
							curBlockPtr;
						
						zoneCounts[z]++;
					}
				}
			}
			
			curBlockPtr++;
		}
	}
	
	return(numZones);
}


/*********************************************************************************
*
*		Name:			detBlocIntraFreePaths
*
*		Summary:		Given a photon with its location and direction and the 
*							block that it is in, determine the total number of 
*							free paths consumed in traveling within the block the 
*							supplied travel distance to the block exit point.
*
*		Arguments:
*			PHG_TrackingPhoton	*thePhotonPtr	- Info about the photon.
*			DetectorBlockPtr	theBlockPtr		- Pointer to the block.
*			double				travelDistance	- Distance to block exit point.
*
*		Function return: 	Free paths consumed along the path.
*
*********************************************************************************/

double detBlocIntraFreePaths(	PHG_TrackingPhoton	*thePhotonPtr, 
								DetectorBlockPtr	theBlockPtr, 
								double				travelDistance)

{
	double			freePaths;			/* Returned function value */
	double			curAtten;			/* Current block material attenuation */
	
	
	if (travelDistance <= 0.0) {
		return(0.0);		/* NOTE!!!  Alternate return point */
	}
	
	
	if ((theBlockPtr->userBlockData->NumLayers == 1) && 
			(theBlockPtr->userBlockData->LayerInfo[0].NumElements == 1)) {
		/* Only one element in the block */
		
		/* Get the attenuation of the block */
		SubObjGetAttenuationInTomo(
			theBlockPtr->userBlockData->LayerInfo[0].ElementInfo[0].MaterialIndex, 
			thePhotonPtr->energy, &curAtten);
		
		/* Accumulate the free paths */
		freePaths = curAtten * travelDistance;
	}
	else {
		/* Multiple elements in the block */
		
		PHG_Position	curPos;			/* Current position along photon path through block */
		PHG_Direction	curDir;			/* Direction of photon path through block */
		Boolean			xPosDir;		/* If photon traveling in positive x direction */
		Boolean			yPosDir;		/* If photon traveling in positive y direction */
		Boolean			zPosDir;		/* If photon traveling in positive z direction */
		Boolean			xNegDir;		/* If photon traveling in negative x direction */
		Boolean			yNegDir;		/* If photon traveling in negative y direction */
		Boolean			zNegDir;		/* If photon traveling in negative z direction */
		DetElementPosIndexTy
						thePosIndex;	/* Position index of the current element */
		Boolean			validElem;		/* Whether current element was found in block */
		DetBlockTomoElementTy
						*curElemDataPtr;/* Pointer to current element */
		PHG_Position	elemCorner1;	/* Position of "smallest" element corner */
		PHG_Position	elemCorner2;	/* Position of "largest" element corner */
		double			xElemTravDist;	/* Element travel distance in x direction */
		double			yElemTravDist;	/* Element travel distance in y direction */
		double			zElemTravDist;	/* Element travel distance in z direction */
		double			curElemTravDist;/* Total element travel distance */
		
		
		/* Transform the photon position and direction from tomo coords to block coords */
		detBlocTomoToBlockCoord(theBlockPtr->userBlockData, &(thePhotonPtr->location), &curPos);
		detBlocTomoToBlockDirection(theBlockPtr->userBlockData, &(thePhotonPtr->angle), &curDir);
		
		/* Assuming the photon really is in the block, make sure roundoff didn't contradict that */
		if (curPos.x_position < theBlockPtr->userBlockData->XMin ) {
			curPos.x_position = theBlockPtr->userBlockData->XMin;
		}
		else if (curPos.x_position > theBlockPtr->userBlockData->XMax ) {
			curPos.x_position = theBlockPtr->userBlockData->XMax;
		}
		if (curPos.y_position < theBlockPtr->userBlockData->YMin ) {
			curPos.y_position = theBlockPtr->userBlockData->YMin;
		}
		else if (curPos.y_position > theBlockPtr->userBlockData->YMax ) {
			curPos.y_position = theBlockPtr->userBlockData->YMax;
		}
		if (curPos.z_position < theBlockPtr->userBlockData->ZMin ) {
			curPos.z_position = theBlockPtr->userBlockData->ZMin;
		}
		else if (curPos.z_position > theBlockPtr->userBlockData->ZMax ) {
			curPos.z_position = theBlockPtr->userBlockData->ZMax;
		}
		
		/* Determine the +/- travel directions */
		xPosDir = (curDir.cosine_x > 1E-12);
		yPosDir = (curDir.cosine_y > 1E-12);
		zPosDir = (curDir.cosine_z > 1E-12);
		xNegDir = (curDir.cosine_x < 1E-12);
		yNegDir = (curDir.cosine_y < 1E-12);
		zNegDir = (curDir.cosine_z < 1E-12);
		
		
		/* Proceed through the block elements */
		thePosIndex.ringNum = theBlockPtr->itsRing;
		thePosIndex.blockNum = theBlockPtr->itsBlock;
		freePaths = 0.0;
		do {
			/* Get the current position's element */
			validElem = detBlocGetElementIndex(&curPos, &thePosIndex);
			
			if (validElem) {
				/* Get the current position's element data */
				curElemDataPtr = &(theBlockPtr->userBlockData->LayerInfo[thePosIndex.layerNum
														].ElementInfo[thePosIndex.elementNum]);
				detBlocGetElementCorners(&thePosIndex, &elemCorner1, &elemCorner2);
				
				/* Get the current element (x,y,z) travel distances */
				xElemTravDist = 100000000.0;
				if (xPosDir) {
					xElemTravDist = (elemCorner2.x_position - curPos.x_position) / curDir.cosine_x;
				}
				else if (xNegDir) {
					xElemTravDist = (elemCorner1.x_position - curPos.x_position) / curDir.cosine_x;
				}
				yElemTravDist = 100000000.0;
				if (yPosDir) {
					yElemTravDist = (elemCorner2.y_position - curPos.y_position) / curDir.cosine_y;
				}
				else if (yNegDir) {
					yElemTravDist = (elemCorner1.y_position - curPos.y_position) / curDir.cosine_y;
				}
				zElemTravDist = 100000000.0;
				if (zPosDir) {
					zElemTravDist = (elemCorner2.z_position - curPos.z_position) / curDir.cosine_z;
				}
				else if (zNegDir) {
					zElemTravDist = (elemCorner1.z_position - curPos.z_position) / curDir.cosine_z;
				}
				
				/* Find the current element travel distance */
				if (xElemTravDist <= yElemTravDist) {
					if (xElemTravDist <= zElemTravDist) {
						/* X distance is shortest */
						curElemTravDist = xElemTravDist;
					}
					else {
						/* Z distance is shortest */
						curElemTravDist = zElemTravDist;
					}
				}
				else {
					if (yElemTravDist <= zElemTravDist) {
						/* Y distance is shortest */
						curElemTravDist = yElemTravDist;
					}
					else {
						/* Z distance is shortest */
						curElemTravDist = zElemTravDist;
					}
				}
				
				/* Advance to the element boundary */
				{
					/* Proceed slightly beyond boundary to avoid element membership quandary */
					curElemTravDist += 1E-9;
					
					curPos.x_position += curElemTravDist * curDir.cosine_x;
					curPos.y_position += curElemTravDist * curDir.cosine_y;
					curPos.z_position += curElemTravDist * curDir.cosine_z;
					travelDistance -= curElemTravDist;
					SubObjGetAttenuationInTomo(
						curElemDataPtr->MaterialIndex, thePhotonPtr->energy, &curAtten);
					freePaths += curElemTravDist * curAtten;
				}
			}
			else {
				/* Photon previously exited the block, so done */
				travelDistance = 0.0;
			}
		} while (fabs(travelDistance) > 1E-7);
	}
	
	/* Return the value */
	return (freePaths);
}


/*********************************************************************************
*
*		Name:			detBlocIntraDistance
*
*		Summary:		Given a photon with its location and direction and the 
*							block that it is in, determine the distance to travel 
*							that consumes the supplied number of free paths.
*							The maximum travel distance within the block is also 
*							supplied.  If the returned distance to travel exceeds 
*							the maximum travel distance, instead return the 
*							distance that would be traveled if the block extended 
*							far enough to consume all of the free paths.  In either 
*							case, return the total number of free paths actually 
*							consumed in traveling only within the block along the 
*							photon's path.  Finally, return the material index and 
*							active status of the last material traveled through.
*
*		Arguments:
*			PHG_TrackingPhoton	*thePhotonPtr	- Info about the photon.
*			DetectorBlockPtr	theBlockPtr		- Pointer to the block.
*			double				freePaths		- Free paths to travel.
*			double				maxTravelDist	- Maximum block travel distance.
*			double				*travelDistance	- Free paths travel distance.
*			double				*freePathsUsed	- Actual consumed free paths.
*			LbUsFourByte		*detMaterial	- The last used detector material.
*			Boolean				*isActive		- Active status of detMaterial.
*
*		Function return: 	None.
*
*********************************************************************************/

void detBlocIntraDistance(	PHG_TrackingPhoton	*thePhotonPtr, 
							DetectorBlockPtr	theBlockPtr, 
							double				freePaths, 
							double				maxTravelDist, 
							double				*travelDistance, 
							double				*freePathsUsed, 
							LbUsFourByte		*detMaterial, 
							Boolean				*isActive)

{
	double			curAtten;			/* Current block material attenuation */
	double			distance;			/* Virtual distance traveled */
	
	
	if ((theBlockPtr->userBlockData->NumLayers == 1) && 
			(theBlockPtr->userBlockData->LayerInfo[0].NumElements == 1)) {
		/* Only one element in the block */
		
		/* Set the block material */
		*detMaterial = theBlockPtr->userBlockData->LayerInfo[0].ElementInfo[0].MaterialIndex;
		
		/* Get the attenuation of the block */
		SubObjGetAttenuationInTomo(*detMaterial, thePhotonPtr->energy, &curAtten);
		
		/* Compute the distance to travel, if no boundary intersections */
		distance = freePaths / curAtten;
		
		if (distance <= maxTravelDist) {
			/* Photon consumed all free paths and stayed within the block */
			*freePathsUsed = freePaths;
		}
		else {
			/* Photon left the block without consuming all of the free paths */
			*freePathsUsed = maxTravelDist * curAtten;
		}
		
		/* Return the remaining values */
		*travelDistance = distance;
		*isActive = theBlockPtr->userBlockData->LayerInfo[0].ElementInfo[0].IsActive;
	}
	else {
		/* Multiple elements in the block */
		
		PHG_Position	curPos;			/* Current position along photon path through block */
		PHG_Direction	curDir;			/* Direction of photon path through block */
		Boolean			xPosDir;		/* If photon traveling in positive x direction */
		Boolean			yPosDir;		/* If photon traveling in positive y direction */
		Boolean			zPosDir;		/* If photon traveling in positive z direction */
		Boolean			xNegDir;		/* If photon traveling in negative x direction */
		Boolean			yNegDir;		/* If photon traveling in negative y direction */
		Boolean			zNegDir;		/* If photon traveling in negative z direction */
		DetElementPosIndexTy
						thePosIndex;	/* Position index of the current element */
		double			curFreePaths;	/* Current remaining free paths to consume */
		Boolean			done;			/* Whether element searching is completed */
		Boolean			validElem;		/* Whether current element was found in block */
		DetBlockTomoElementTy
						*curElemDataPtr;/* Pointer to current element */
		PHG_Position	elemCorner1;	/* Position of "smallest" element corner */
		PHG_Position	elemCorner2;	/* Position of "largest" element corner */
		double			xElemTravDist;	/* Element travel distance in x direction */
		double			yElemTravDist;	/* Element travel distance in y direction */
		double			zElemTravDist;	/* Element travel distance in z direction */
		double			curElemTravDist;/* Total element travel distance */
		LbFourByte		elemMaterial;	/* Material of current element */
		double			curDistance;	/* Current element travel distance */
		
		
		/* Transform the photon position and direction from tomo coords to block coords */
		detBlocTomoToBlockCoord(theBlockPtr->userBlockData, &(thePhotonPtr->location), &curPos);
		detBlocTomoToBlockDirection(theBlockPtr->userBlockData, &(thePhotonPtr->angle), &curDir);
		
		/* Assuming the photon really is in the block, make sure roundoff didn't contradict that */
		if (curPos.x_position < theBlockPtr->userBlockData->XMin ) {
			curPos.x_position = theBlockPtr->userBlockData->XMin;
		}
		else if (curPos.x_position > theBlockPtr->userBlockData->XMax ) {
			curPos.x_position = theBlockPtr->userBlockData->XMax;
		}
		if (curPos.y_position < theBlockPtr->userBlockData->YMin ) {
			curPos.y_position = theBlockPtr->userBlockData->YMin;
		}
		else if (curPos.y_position > theBlockPtr->userBlockData->YMax ) {
			curPos.y_position = theBlockPtr->userBlockData->YMax;
		}
		if (curPos.z_position < theBlockPtr->userBlockData->ZMin ) {
			curPos.z_position = theBlockPtr->userBlockData->ZMin;
		}
		else if (curPos.z_position > theBlockPtr->userBlockData->ZMax ) {
			curPos.z_position = theBlockPtr->userBlockData->ZMax;
		}
		
		/* Determine the +/- travel directions */
		xPosDir = (curDir.cosine_x > 1E-12);
		yPosDir = (curDir.cosine_y > 1E-12);
		zPosDir = (curDir.cosine_z > 1E-12);
		xNegDir = (curDir.cosine_x < -1E-12);
		yNegDir = (curDir.cosine_y < -1E-12);
		zNegDir = (curDir.cosine_z < -1E-12);
		
		
		/* Proceed through the block elements */
		thePosIndex.ringNum = theBlockPtr->itsRing;
		thePosIndex.blockNum = theBlockPtr->itsBlock;
		curFreePaths = freePaths;
		distance = 0.0;
		done = false;
		do {
			if (distance > maxTravelDist) {
				/* Photon has exited the block */
				validElem = false;
			}
			else {
				/* Get the current position's element */
				validElem = detBlocGetElementIndex(&curPos, &thePosIndex);
			}
			
			if (validElem) {
				/* Get the current position's element data */
				curElemDataPtr = &(theBlockPtr->userBlockData->LayerInfo[thePosIndex.layerNum
														].ElementInfo[thePosIndex.elementNum]);
				detBlocGetElementCorners(&thePosIndex, &elemCorner1, &elemCorner2);
				
				/* Get the current element (x,y,z) travel distances */
				xElemTravDist = 100000000.0;
				if (xPosDir) {
					xElemTravDist = (elemCorner2.x_position - curPos.x_position) / curDir.cosine_x;
				}
				else if (xNegDir) {
					xElemTravDist = (elemCorner1.x_position - curPos.x_position) / curDir.cosine_x;
				}
				yElemTravDist = 100000000.0;
				if (yPosDir) {
					yElemTravDist = (elemCorner2.y_position - curPos.y_position) / curDir.cosine_y;
				}
				else if (yNegDir) {
					yElemTravDist = (elemCorner1.y_position - curPos.y_position) / curDir.cosine_y;
				}
				zElemTravDist = 100000000.0;
				if (zPosDir) {
					zElemTravDist = (elemCorner2.z_position - curPos.z_position) / curDir.cosine_z;
				}
				else if (zNegDir) {
					zElemTravDist = (elemCorner1.z_position - curPos.z_position) / curDir.cosine_z;
				}
				
				/* Find the current element travel distance */
				if (xElemTravDist <= yElemTravDist) {
					if (xElemTravDist <= zElemTravDist) {
						/* X distance is shortest */
						curElemTravDist = xElemTravDist;
					}
					else {
						/* Z distance is shortest */
						curElemTravDist = zElemTravDist;
					}
				}
				else {
					if (yElemTravDist <= zElemTravDist) {
						/* Y distance is shortest */
						curElemTravDist = yElemTravDist;
					}
					else {
						/* Z distance is shortest */
						curElemTravDist = zElemTravDist;
					}
				}
				
				/* Get the element material */
				elemMaterial = curElemDataPtr->MaterialIndex;
				
				/* Get the attenuation of the element */
				SubObjGetAttenuationInTomo(elemMaterial, thePhotonPtr->energy, &curAtten);
				
				/* Compute the distance to travel, if no boundary intersections */
				curDistance = curFreePaths / curAtten;
				
				if (curDistance > curElemTravDist) {
					/* Photon exited the element; advance to the element boundary */
					{
						double		nudgeDist = 1E-9;			/* Nudge distance */
						
						
						/* Proceed slightly beyond boundary to avoid element membership quandary */
						curElemTravDist += nudgeDist;
						
						curPos.x_position += curElemTravDist * curDir.cosine_x;
						curPos.y_position += curElemTravDist * curDir.cosine_y;
						curPos.z_position += curElemTravDist * curDir.cosine_z;
						
						curFreePaths -= (curElemTravDist * curAtten);
						if ( curFreePaths < 0.0 ) {
							/* The distance nudging has pushed it beyond the allowed free paths; 
								correct it using the average error */
							curFreePaths = nudgeDist / 2 * curAtten;
						}
						distance += curElemTravDist;
					}
				}
				else {
					/* Photon stopped within the element, so done */
					distance += curDistance;
					*travelDistance = distance;
					*freePathsUsed = freePaths;
					*detMaterial = elemMaterial;
					*isActive = curElemDataPtr->IsActive;
					
					done = true;
				}
			}
			else {
				/* Photon previously exited the block, so done */
				
				if (curFreePaths > 0.0) {
					/* Increase to extended travel distance */
					distance += curFreePaths / curAtten;
				}
				*travelDistance = distance;
				*freePathsUsed = freePaths - curFreePaths;
				*detMaterial = elemMaterial;
				*isActive = curElemDataPtr->IsActive;
				
				done = true;
			}
		} while ( ! done );
	}
}


/*********************************************************************************
*
*		Name:			detBlocGetDistToExit
*
*		Summary:		Given a photon with its location and direction, 
*							determine the exit point from the given block and 
*							the travel distance to that point.
*
*		Arguments:
*			PHG_TrackingPhoton	*thePhotonPtr	- Info about the photon.
*			DetectorBlockPtr	theBlockPtr		- Pointer to the block.
*			PHG_Position		*blockExitPoint	- Photon exit point from block.
*			double				*travelDist		- Distance to exit point.
*
*		Function return: 	Indicator of how photon exited the block.
*
*********************************************************************************/

DetBlocEvents detBlocGetDistToExit(	PHG_TrackingPhoton	*thePhotonPtr, 
									DetectorBlockPtr	theBlockPtr, 
									PHG_Position		*blockExitPoint, 
									double				*travelDist)

{
	PHG_Position		exitPoint;		/* Local block exit point */
	double				distanceTraveled;/* Local distance to exit point */
	DetBlocEvents		exitManner;		/* Function return value */
	PHG_Position		phoOrigPos;		/* Starting photon position */
	Lb2D_Point			point1;			/* Photon position as 2D point */
	DetBlockTomoRingTy	*curRingPtr;	/* Ring info for the block's ring */
	double				ringMinus,		/* Lesser and greater detector Z-boundaries */
						ringPlus;
	Lb2D_Rect			theBlockRect;	/* Block boundary 2D rectangle */
	Lb2D_PosDesc		phoBlockPos;	/* Relation of photon to block rectangle */
	int					validEdges;		/* Indicates valid rect intersection edges */
	Lb2D_Point			*edgeCornerPtrs[5];/* Corners (wrapped) of theBlockRect */
	int					e;				/* Index through rect edges */
	double				edgeCos, 		/* Normal line descriptors of rect edge */
						edgeSin, 
						edgeDist;
	
	
	/* Failure indicators */
	exitPoint.x_position = 0.0;
	exitPoint.y_position = 0.0;
	exitPoint.z_position = 0.0;
	distanceTraveled = 0.0;
	exitManner = detBlocEvt_InnerCyl;
	
	/* Save the photon's current position */
	phoOrigPos = thePhotonPtr->location;
	
	/* Local 2D variable of starting point */
	point1.x_position = phoOrigPos.x_position;
	point1.y_position = phoOrigPos.y_position;
	
	/* Get the ring boundaries of the detector */
	curRingPtr = 
		&(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[theBlockPtr->itsRing]);
	ringMinus = curRingPtr->MinZ + curRingPtr->AxialShift;
	ringPlus = curRingPtr->MaxZ + curRingPtr->AxialShift;
	
	/* Get the given block's 2D rectangle boundary */
	theBlockRect = theBlockPtr->itsRect;
	
	/* Compare the photon's position to the block rectangle */
	phoBlockPos = Lb2DPointRectIntersect(&point1, &theBlockRect);
	
	validEdges = 15;	/* = all edges */
	switch (phoBlockPos) {
		/* NOTE:  case fallthroughs are >>>intentional<<<  */
		
		case Lb2D_Outside:
		{
			/* NOTE:  This function will normally be called only when the photon is 
						on or within the block boundaries */
			
			LbFourByte			itsRing;		/* Current ring */
			LbFourByte			itsZone;		/* Current zone */
			DetectorBlockPtr	nextBlockPtr;	/* Next block */
			LbFourByte			blockRing;		/* Ignored */
			LbFourByte			blockZone;		/* Ignored */
			LbFourByte			blockIndex;		/* Ignored */
			
			/* Project the photon to the block boundary */
			itsRing = -1;
			itsZone = -1;
			if (! detBlocGetZone(&(thePhotonPtr->location), &itsRing, &itsZone)) {
				/* Not within a ring; must project until is */
				if (detBlocProjAcrossGap(thePhotonPtr)) {
					/* Get its new ring and zone */
					if (! detBlocGetZone(&(thePhotonPtr->location), &itsRing, &itsZone)) {
						/* Possible roundoff problem? */
						PhgAbort("Wasn't in projected ring (detBlocGetDistToExit).", false);
					}
				}
				else {
					/* Photon couldn't enter a ring; probably because it exited */
					exitManner = detBlocEvt_OuterCyl;
					break;
				}
			}
			if (detBlocNextBlock(thePhotonPtr, itsRing, itsZone, NULL, 
								&nextBlockPtr, &blockRing, &blockZone, &blockIndex)) {
				if (nextBlockPtr == theBlockPtr) {
					/* Update the first point of the line segment */
					point1.x_position = thePhotonPtr->location.x_position;
					point1.y_position = thePhotonPtr->location.y_position;
				}
				else {
					/* Could not project to the requested block!!! */
					break;
				}
			}
			else {
				/* Could not project to the requested block!!! */
				break;
			}
		}
		/* Fall through to the next case... */
		
		case Lb2D_OnBound:
		{
			/* Determine which edges the photon is on */
			
			double			sideDistance;			/* Distance to rect edge */
			
			/* Set up a point array for looping through the block's sides */
			edgeCornerPtrs[0] = &(theBlockRect.corner_1);
			edgeCornerPtrs[1] = &(theBlockRect.corner_2);
			edgeCornerPtrs[2] = &(theBlockRect.corner_3);
			edgeCornerPtrs[3] = &(theBlockRect.corner_4);
			edgeCornerPtrs[4] = &(theBlockRect.corner_1);
			
			/* Check each side for intersection */
			for (e=0; e<4; e++) {
				/* Calculate the 2D distance to the side */
				Lb2DNormalLine(edgeCornerPtrs[e], edgeCornerPtrs[e+1], 
								&edgeCos, &edgeSin, &edgeDist);
				sideDistance = fabs((edgeCos*point1.x_position + 
										edgeSin*point1.y_position + 
										edgeDist));
				
				if (PhgMathRealNumAreEqual(sideDistance, 0.0, -7, 0, 0, 0)) {
					/* Remove the side from the valid exit edges */
					validEdges ^= (1<<e);
				}
			}
		}
		/* Fall through to the next case... */
		
		case Lb2D_Inside:
		{
			/* Photon is within the block */
			
			double			kAxialZ = 1E-8;	/* Value axial direction considered pure axial; 
												Carefully balanced between precision of doubles 
												and cos, sqrt consequences */
			Lb2D_Point		point2;			/* Endpoint of path as line segment */
			double			segCos, 		/* Normal line descriptors of path segment */
							segSin, 
							segDist;
			double			sineZ2;			/* Sine^2 of photon Z-direction */
			double			hitDistance;	/* 3D distance to hit point */
			Boolean			edgeBounce;		/* True if photon bounced off block surface */
			Lb2D_PosDesc	sideIntersect;	/* Relation of path line segment to rect edge */
			double			det;			/* Determinant of non-parallel line coeffs */
			double			hitX,			/* Coordinates of possible hit point */
							hitY;
			double			sideDistance2;	/* 2D distance^2 to hit point */
			
			
			/* Set up a point array for looping through the block's sides */
			edgeCornerPtrs[0] = &(theBlockRect.corner_1);
			edgeCornerPtrs[1] = &(theBlockRect.corner_2);
			edgeCornerPtrs[2] = &(theBlockRect.corner_3);
			edgeCornerPtrs[3] = &(theBlockRect.corner_4);
			edgeCornerPtrs[4] = &(theBlockRect.corner_1);
			
			if ((1.0-fabs(thePhotonPtr->angle.cosine_z)) <= kAxialZ) {
				/* Very axially directed photon; immediately project to the ring boundary */
				exitPoint.x_position = point1.x_position;
				exitPoint.y_position = point1.y_position;
				if (thePhotonPtr->angle.cosine_z > 0.0) {
					exitPoint.z_position = ringPlus;
					distanceTraveled = ringPlus - phoOrigPos.z_position;
				}
				else {
					exitPoint.z_position = ringMinus;
					distanceTraveled = phoOrigPos.z_position - ringMinus;
				}
				exitManner = detBlocEvt_Block;
				
				break;		/* All done here */
			}
			
			/* Convert the photon's future path to a very long line segment */
			point2.x_position = point1.x_position + 8388608.0 * thePhotonPtr->angle.cosine_x;
			point2.y_position = point1.y_position + 8388608.0 * thePhotonPtr->angle.cosine_y;
			Lb2DNormalLine(&point1, &point2, &segCos, &segSin, &segDist);
			
			/* Compute photon z angle sine^2 */
			sineZ2 = 									/* Won't be 0 */
				(1.0 - PHGMATH_Square(thePhotonPtr->angle.cosine_z));
			
			/* Check for intersection with each of the block's sides */
			hitDistance = -10.0;
			edgeBounce = true;
			for (e=0; e<4; e++) {
				sideIntersect = Lb2DLineSegsIntersect(&point1, &point2, 
									edgeCornerPtrs[e], edgeCornerPtrs[e+1]);
				if (sideIntersect != Lb2D_Outside) {
					/* There was an intersection */
					
					/* Make sure it is a valid one */
					if ((1<<e) & validEdges) {
						edgeBounce = false;
						
						/* Calculate the 2D and 3D distances to the side */
						Lb2DNormalLine(edgeCornerPtrs[e], edgeCornerPtrs[e+1], 
										&edgeCos, &edgeSin, &edgeDist);
						det = segCos*edgeSin - edgeCos*segSin;	/* Won't be 0 */
						hitX = (segSin*edgeDist - edgeSin*segDist) / det;
						hitY = (edgeCos*segDist - segCos*edgeDist) / det;
						sideDistance2 = PHGMATH_Square(hitX - point1.x_position) + 
										PHGMATH_Square(hitY - point1.y_position);
						hitDistance = PHGMATH_SquareRoot(sideDistance2 / sineZ2);
						
						/* Only one side intersection needs to be considered */
						break;
					}
				}
			}
			
			if (edgeBounce) {
				/* The photon was inside the block but crossed none of the edges it 
					wasn't on; therefore it must have bounced outwards off of an edge */
				
				/* Take the maximum distance to each edge the photon is on */
				for (e=0; e<4; e++) {
					if (!((1<<e) & validEdges)) {
						double		tempX,			/* Coordinates of possible hit point */
									tempY;
						double		tempDistance;	/* 3D distance to hit point */
						
						
						/* Calculate the 2D and 3D distances to the side */
						Lb2DNormalLine(edgeCornerPtrs[e], edgeCornerPtrs[e+1], 
										&edgeCos, &edgeSin, &edgeDist);
						det = segCos*edgeSin - edgeCos*segSin;	/* Won't be 0 */
						tempX = (segSin*edgeDist - edgeSin*segDist) / det;
						tempY = (edgeCos*segDist - segCos*edgeDist) / det;
						sideDistance2 = PHGMATH_Square(tempX - point1.x_position) + 
										PHGMATH_Square(tempY - point1.y_position);
						tempDistance = PHGMATH_SquareRoot(sideDistance2 / sineZ2);
						
						/* Take the largest one */
						if (tempDistance > hitDistance) {
							hitX = tempX;
							hitY = tempY;
							hitDistance = tempDistance;
						}
					}
				}
			}
			
			if (hitDistance >= 0.0) {
				/* Find the intersection point and the distance to it */
				
				/* Determine the z position for this exit point */
				exitPoint.z_position = thePhotonPtr->location.z_position + 
					hitDistance * thePhotonPtr->angle.cosine_z;
				
				/* Compare the exit point to the ring boundaries */
				if (exitPoint.z_position < ringMinus) {
					/* Photon exited the lesser side of the ring; adjust the exit point */
					exitPoint.z_position = ringMinus;
					distanceTraveled = (ringMinus - thePhotonPtr->location.z_position) / 
										thePhotonPtr->angle.cosine_z;  /* cos < 0 */
					exitPoint.x_position = point1.x_position + 
						distanceTraveled * thePhotonPtr->angle.cosine_x;
					exitPoint.y_position = point1.y_position + 
						distanceTraveled * thePhotonPtr->angle.cosine_y;
					exitManner = detBlocEvt_PrevRing;
				}
				else if (exitPoint.z_position > ringPlus) {
					/* Photon exited the greater side of the ring; adjust the exit point */
					exitPoint.z_position = ringPlus;
					distanceTraveled = (ringPlus - thePhotonPtr->location.z_position) / 
										thePhotonPtr->angle.cosine_z;
					exitPoint.x_position = point1.x_position + 
						distanceTraveled * thePhotonPtr->angle.cosine_x;
					exitPoint.y_position = point1.y_position + 
						distanceTraveled * thePhotonPtr->angle.cosine_y;
					exitManner = detBlocEvt_NextRing;
				}
				else {
					/* Photon exited inside the ring; the exit point has been found */
					exitPoint.x_position = hitX;
					exitPoint.y_position = hitY;
					distanceTraveled = hitDistance;
					exitManner = detBlocEvt_Block;
				}
			}
		}
		/* Fall through to the next case... */
		
		default:
			/* Nothing else left */
			break;
	}
	
	/* Restore the photon's position */
	thePhotonPtr->location = phoOrigPos;
	
	/* Return the results */
	*blockExitPoint = exitPoint;
	*travelDist = distanceTraveled;
	return (exitManner);
}


/*********************************************************************************
*
*		Name:			detBlocCompFreePathsToExit
*
*		Summary:		Compute the free-paths necessary for the photon
*							to exit all of the detector blocks.
*
*		Arguments:
*			PHG_TrackingPhoton		*photonPtr		- The photon.
*			LbUsFourByte			curRingIndex	- The current ring.
*			LbUsFourByte			curDetIndex		- The current detector block.
*
*		Function return: The free paths to exit.
*
*********************************************************************************/

double detBlocCompFreePathsToExit(PHG_TrackingPhoton *photonPtr, 
								LbUsFourByte curRingIndex, LbUsFourByte curDetIndex)

{
	double				fpToGo;			/* Returned number of free paths */
	PHG_Position		phoOrigPos;		/* Starting photon position */
	double				phoOrigDist;	/* Starting photon travel distance */
	LbFourByte			curRing;		/* Current ring */
	DetectorBlockPtr	curBlockPtr;	/* Current block */
	Lb2D_Rect			theBlockRect;	/* Block boundary 2D rectangle */
	PHG_Position		blockExitPoint;	/* Exit point of current block */
	double				travelDistance;	/* Travel distance in current block */
	LbFourByte			curZone;		/* Current zone */
	LbFourByte			curIndex;		/* Current zone block index */
	Boolean				newBlock;		/* If another block is on the path */
	
	
	/* Compute the free paths on the photon's current path, intersecting all 
		possible blocks on the path line */
	
	fpToGo = 0.0;
	
	/* Save the photon's current position and travel distance */
	phoOrigPos = photonPtr->location;
	phoOrigDist = photonPtr->travel_distance;
	
	curRing = curRingIndex;
	
	/* Get the current block */
	curBlockPtr = NULL;		/* Just to be safe */
	if (detBlocGetNumberedBlock(curRingIndex, curDetIndex, &curBlockPtr)) {
		
		/* Get the given block's 2D rectangle boundary */
		theBlockRect = curBlockPtr->itsRect;
		
		if (Lb2DPointRectIntersect((Lb2D_Point*)&(photonPtr->location), &theBlockRect) != 
									Lb2D_Outside) {
			/* Compute the exit point and the distance traveled to it */
			detBlocGetDistToExit(photonPtr, curBlockPtr, 
									&blockExitPoint, &travelDistance);
			
			/* Accumulate the free paths along the exit path */
			fpToGo += detBlocIntraFreePaths(photonPtr, curBlockPtr, travelDistance);
			
			/* Move the photon to the exit point */
			photonPtr->location = blockExitPoint;
		}
		else {
			/* Photon wasn't in a block yet */
			curBlockPtr = NULL;
		}
	}
	else {
		/* Bad block; go ahead and re-project to detectors */
	}
	
	do {
		/* Find and move to the next block intersected by the photon's path, if any */
		curZone = -1;
		if (! detBlocGetZone(&(photonPtr->location), &curRing, &curZone)) {
			/* Not within a ring; must project until is */
			if (detBlocProjAcrossGap(photonPtr)) {
				/* Get its new ring and zone */
				if (! detBlocGetZone(&(photonPtr->location), &curRing, &curZone)) {
					/* Possible roundoff problem? */
					PhgAbort("Wasn't in projected ring (detBlocCompFreePathsToExit).", false);
				}
			}
			else {
				/* Photon couldn't enter a ring; probably because it exited */
				break;
			}
			curBlockPtr = NULL;
		}
		newBlock = detBlocNextBlock(photonPtr, curRing, curZone, curBlockPtr, 
									&curBlockPtr, 
									&curRing, &curZone, &curIndex);
		
		if (newBlock) {
			/* Get the given block's 2D rectangle boundary */
			theBlockRect = curBlockPtr->itsRect;
			
			if (Lb2DPointRectIntersect((Lb2D_Point*)&(photonPtr->location), &theBlockRect) != 
										Lb2D_Outside) {
				/* Compute the exit point and the distance traveled to it */
				detBlocGetDistToExit(photonPtr, curBlockPtr, 
										&blockExitPoint, &travelDistance);
				
				/* Accumulate the free paths along the exit path */
				fpToGo += detBlocIntraFreePaths(photonPtr, curBlockPtr, travelDistance);
				
				/* Move the photon to the exit point */
				photonPtr->location = blockExitPoint;
			}
			else {
				/* Photon wasn't in a block yet */
			}
		}
	} while (newBlock);
	
	/* Restore the photon's position and travel distance */
	photonPtr->location = phoOrigPos;
	photonPtr->travel_distance = phoOrigDist;
	
	return (fpToGo);
}


/*********************************************************************************
*
*		Name:			DetBlocGtTruncatedFreePaths
*
*		Summary:		Return the free paths to travel for a forced first interaction.
*
*		Arguments:		
*			PHG_TrackingPhotonPtr	photonPtr		- The photon.
*			LbUsFourByte			ringNum			- The current ring.
*			LbUsFourByte			detNum			- The current detector.
*			double					*weight			- The photon weight.
*
*		Function return: Free-paths to travel.
*
*********************************************************************************/

double DetBlocGtTruncatedFreePaths(PHG_TrackingPhoton *photonPtr, 
						LbUsFourByte ringNum, LbUsFourByte detNum, double *weight)

{
	double 			freePathsToExit;		/* Unadjusted free paths to exit */
	double 			newWeight;
	double			fpToGo;
	double			randFromExp;
	LbUsFourByte	wholeNum;				/* Whole number temp variable */
	
	
	/* Compute the free paths on the photon's current path, intersecting all 
		possible blocks on the path line */
	/* Copied from SV -- can be improved via global function change */
	
	
	/* Compute the free paths to exit */
	freePathsToExit = detBlocCompFreePathsToExit(photonPtr, ringNum, detNum);
	
	/* Adjust the photon's weight by the probability that an interaction would occur */
	newWeight = (*weight) * (1 - exp(-freePathsToExit));
	
	/* Increment global counter for reporting */
	detData[DetCurParams].detWeightAdjusted += ((*weight) - newWeight);
	
	/* Update weight */
	*weight = newWeight;
	
	
	/* Pick a free path to go based on a truncated exponential distribution */
	{
		/* Start with a random from the exponential distribution */
		PhgMathGetTotalFreePaths(&randFromExp);
		
		/* Truncate to desired range */
		{
			/* Now compute the free paths to go */
			/* Truncate to desired range */
			{
				/* Compute whole number part NOTE I AM IGNORING THE
				   POSSIBILITY OF INTEGER OVERFLOW
				*/
				wholeNum = (LbUsFourByte) (randFromExp/
					freePathsToExit);

				fpToGo = ((randFromExp/
			   		freePathsToExit) -
					wholeNum)
					* freePathsToExit;
					
				if (fpToGo > freePathsToExit) {
					#ifdef HEAVY_PHG_DEBUG
						/* send alert, but don't terminate simulation */
						ErAlert("Invalid calculation of fpToGo for forced interaction (detBlock)--PLEASE REPORT", false);
					#endif
					/* this should be an infrequent numerical condition, set to boundary value */
					fpToGo = freePathsToExit;
				}
			}

		}	
	}				
	
	return (fpToGo);
}


/*********************************************************************************
*
*			Name:			detBlocNextBlock
*
*			Summary:		Given a photon with its location and direction, 
*								determine the next detector block it will enter 
*								and move the photon to its edge.
*
*			Arguments:
*				PHG_TrackingPhoton	*thePhotonPtr	- The given photon.
*				LbFourByte			itsRing			- Starting ring.
*				LbFourByte			itsZone			- Starting ring zone.
*				DetectorBlockPtr	itsBlockPtr		- Starting block (or NULL).
*				DetectorBlockPtr	*nextBlockPtr	- Next block.
*				LbFourByte			*blockRing		- Next block's ring.
*				LbFourByte			*blockZone		- Next block's ring zone.
*				LbFourByte			*blockIndex		- Next block's zone block index.
*
*			Function return: Boolean; true if block was found, false if exits detector.
*
*********************************************************************************/

Boolean detBlocNextBlock(	PHG_TrackingPhoton	*thePhotonPtr, 
							LbFourByte			itsRing, 
							LbFourByte			itsZone, 
							DetectorBlockPtr	itsBlockPtr, 
							DetectorBlockPtr	*nextBlockPtr, 
							LbFourByte			*blockRing, 
							LbFourByte			*blockZone, 
							LbFourByte			*blockIndex)

{
	Boolean				blockResult;		/* Function result */
	DetBlockTomoRingTy*	curRingPtr;			/* Ring info for the block's ring */
	double				posRingZ;			/* Position of ring boundary (to positive) */
	double				negRingZ;			/* Position of ring boundary (to negative) */
	DetBlocEvents		photonEvent;		/* Destination event of photon */
	Boolean				hitInnerCyl;		/* Whether photon could hit inner cylinder */
	Boolean				outEnd;				/* Whether photon exits end of a cylinder */
	PHG_Position		projPosition;		/* New photon position */
	double				shortestDistance;	/* Distance to nearest boundary */
	DetBlocEvents		zoneEvent;			/* Result of testing zone boundary crossing */
	PHG_Position		zoneBoundPos;		/* Photon/zone intersection */
	double				zoneProjDistance;	/* 3D distance to photon/zone intersection */
	LbFourByte			b;					/* Index through list of zone blocks */
	DetectorBlockPtr	curBlockPtr;		/* Pointer to current block descriptor */
	double				blkDistance;		/* Distance to current block */
	PHG_Position		blkPoint;			/* Intersection with current block */
	LbFourByte			newRing;			/* Ring photon went into */
	DetBlockTomoRingTy*	nextRingPtr;		/* Ring info for the adjacent ring */
	double				nextRingZ;			/* Position of adjacent ring boundary */
	LbFourByte			newZone;			/* Zone photon went into */
	
	
	/* NOTE:  It is assumed that the photon is in a ring and not in a block */
	
	photonEvent = detBlocEvt_Null;
	
	/* Verify that the ring number is valid */
	if ((itsRing < 0) || (itsRing >= detBlocNumRings)) {
		/* Photon is not in detector */
		photonEvent = detBlocEvt_OutEnd;
		blockResult = false;
	}
	else {
		/* Get ring boundaries */
		curRingPtr = 
			&(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[itsRing]);
		negRingZ = curRingPtr->MinZ + curRingPtr->AxialShift;
		posRingZ = curRingPtr->MaxZ + curRingPtr->AxialShift;
		
		
		/* Determine distances to possible boundaries; the shortest is the one */
		
		/* Check for intersection with inner cylinder */
		hitInnerCyl = detBlocHitInnerCylinder(&detBlocInnerCylinder, thePhotonPtr, 
												&projPosition, &shortestDistance);
		
		if ( ! hitInnerCyl) {
			/* Won't hit inner cylinder */
			
			/* Check against the outer cylinder */
			outEnd = ! detBlocHitOuterCylinder(thePhotonPtr, &projPosition, &shortestDistance);
			if (outEnd) {
				/* Photon went out a cylinder end */
				if (thePhotonPtr->angle.cosine_z > 0.0) {
					/* Went into the next ring */
					photonEvent = detBlocEvt_NextRing;
				}
				else {
					/* Went into the previous ring */
					photonEvent = detBlocEvt_PrevRing;
				}
			}
			else {
				/* Hit the cylinder; check for ring exit before hitting the cylinder */
				if (projPosition.z_position > posRingZ) {
					/* Went into the next ring */
					outEnd = true;
					photonEvent = detBlocEvt_NextRing;
				}
				else if (projPosition.z_position < negRingZ) {
					/* Went into the previous ring */
					outEnd = true;
					photonEvent = detBlocEvt_PrevRing;
				}
				else {
					/* Hit the outer cylinder in the ring */
					photonEvent = detBlocEvt_OuterCyl;
				}
			}
			if (outEnd) {
				/* Project the photon to the ring boundary */
				if (photonEvent == detBlocEvt_NextRing) {
					/* Positive ring boundary */
					shortestDistance = 
						(posRingZ - thePhotonPtr->location.z_position) / 
						thePhotonPtr->angle.cosine_z;
					projPosition.x_position = thePhotonPtr->location.x_position  +  
						thePhotonPtr->angle.cosine_x * shortestDistance;
					projPosition.y_position = thePhotonPtr->location.y_position  +  
						thePhotonPtr->angle.cosine_y * shortestDistance;
					projPosition.z_position = posRingZ;
				}
				else {
					/* Negative ring boundary */
					shortestDistance = 
						(negRingZ - thePhotonPtr->location.z_position) / 
						thePhotonPtr->angle.cosine_z;	/* cos(z) < 0 */
					projPosition.x_position = thePhotonPtr->location.x_position  +  
						thePhotonPtr->angle.cosine_x * shortestDistance;
					projPosition.y_position = thePhotonPtr->location.y_position  +  
						thePhotonPtr->angle.cosine_y * shortestDistance;
					projPosition.z_position = negRingZ;
				}
			}
			
			/* Check distance to zone boundaries */
			zoneEvent = detBlocCrossZoneBounds(thePhotonPtr, itsRing, itsZone, 
												&zoneBoundPos, &zoneProjDistance);
			if (zoneEvent == detBlocEvt_Null) {
				/* Won't hit a zone boundary */
				/* So shortest projection already known */
			}
			else {
				/* Check for shorter projection distance */
				if (zoneProjDistance < shortestDistance) {
					/* Hit zone boundary first */
					photonEvent = zoneEvent;
					shortestDistance = zoneProjDistance;
					projPosition = zoneBoundPos;
				}
			}
		}
		else {
			/* Could intersect inner cylinder */
			
			/* Start with inner cylinder intersection as default event */
			photonEvent = detBlocEvt_InnerCyl;
			
			/* Check for shorter ring exit distance */
			if (projPosition.z_position > posRingZ) {
				/* Went into the next ring */
				photonEvent = detBlocEvt_NextRing;
				shortestDistance = 
					(posRingZ - thePhotonPtr->location.z_position) / 
					thePhotonPtr->angle.cosine_z;
				projPosition.x_position = thePhotonPtr->location.x_position  +  
					thePhotonPtr->angle.cosine_x * shortestDistance;
				projPosition.y_position = thePhotonPtr->location.y_position  +  
					thePhotonPtr->angle.cosine_y * shortestDistance;
				projPosition.z_position = posRingZ;
			}
			else if (projPosition.z_position < negRingZ) {
				/* Went into the previous ring */
				photonEvent = detBlocEvt_PrevRing;
				shortestDistance = 
					(negRingZ - thePhotonPtr->location.z_position) / 
					thePhotonPtr->angle.cosine_z;	/* cos(z) < 0 */
				projPosition.x_position = thePhotonPtr->location.x_position  +  
					thePhotonPtr->angle.cosine_x * shortestDistance;
				projPosition.y_position = thePhotonPtr->location.y_position  +  
					thePhotonPtr->angle.cosine_y * shortestDistance;
				projPosition.z_position = negRingZ;
			}
			
			/* Check distance to zone boundaries */
			zoneEvent = detBlocCrossZoneBounds(thePhotonPtr, itsRing, itsZone, 
												&zoneBoundPos, &zoneProjDistance);
			if (zoneEvent == detBlocEvt_Null) {
				/* Won't hit a zone boundary */
				/* So shortest projection already known */
			}
			else {
				/* Check for shorter projection distance */
				if (zoneProjDistance < shortestDistance) {
					/* Hit zone boundary first */
					photonEvent = zoneEvent;
					shortestDistance = zoneProjDistance;
					projPosition = zoneBoundPos;
				}
			}
		}
		
		/* Check each detector block in the starting ring zone */
		for (b=0; b<detBlocNumZoneBlocks; b++) {
			if (detBlocGetBlock(itsRing, itsZone, b, &curBlockPtr)) {
				/* Avoid checking the block being exited */
				if ((itsBlockPtr == NULL) || (itsBlockPtr != curBlockPtr)) {
					/* Check the distance to the block */
					if (detBlocCalcDistanceToBlock(	&(thePhotonPtr->location), 
													&(thePhotonPtr->angle), 
													curBlockPtr, 
													&blkDistance, &blkPoint)) {
						/* The path does cross the block */
						if (blkDistance < shortestDistance) {
							/* New shorter distance */
							photonEvent = detBlocEvt_Block;
							shortestDistance = blkDistance;
							
							/* Record the location */
							projPosition = blkPoint;
							*nextBlockPtr = curBlockPtr;
							*blockIndex = b;
						}
					}
				}
			}
			else {
				/* No more blocks in this zone */
				break;
			}
		}
		
		
		/* For photons leaving the zone, project to the zone boundary */
		switch (photonEvent) {
			case detBlocEvt_NextRing:
				newRing = itsRing + 1;
				if (newRing >= detBlocNumRings) {
					/* Leaving the detector axially at the positive end */
					photonEvent = detBlocEvt_OutEnd;
				}
				else {
					/* Project across the ring gap */
					nextRingPtr = 
						&(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[newRing]);
					nextRingZ = nextRingPtr->MinZ + nextRingPtr->AxialShift;
					detBlocProjAxially(thePhotonPtr, nextRingZ);
					if (detBlocPtInOuterCylinder((Lb2D_Point*)&(thePhotonPtr->location))) {
						if (! detBlocPtInInnerCylinder(&detBlocInnerCylinder, 
													(Lb2D_Point*)&(thePhotonPtr->location))) {
							newZone = itsZone;
							if (detBlocGetZone(&(thePhotonPtr->location), &newRing, &newZone)) {
								/* Photon has been projected to its new (r,z) location */
							}
							else {
								/* Possible roundoff problem? */
								PhgAbort("Wasn't in next ring (detBlocNextBlock).", false);
							}
						}
						else {
							/* Photon entered the inner cylinder between rings */
							photonEvent = detBlocEvt_InnerCyl;
						}
					}
					else {
						/* Photon exited the outer cylinder between rings */
						photonEvent = detBlocEvt_OuterCyl;
					}
				}
				break;
			
			case detBlocEvt_PrevRing:
				newRing = itsRing - 1;
				if (newRing < 0) {
					/* Leaving the detector axially at the negative end */
					photonEvent = detBlocEvt_OutEnd;
				}
				else {
					/* Project across the ring gap */
					nextRingPtr = 
						&(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[newRing]);
					nextRingZ = nextRingPtr->MaxZ + nextRingPtr->AxialShift;
					detBlocProjAxially(thePhotonPtr, nextRingZ);
					if (detBlocPtInOuterCylinder((Lb2D_Point*)&(thePhotonPtr->location))) {
						if (! detBlocPtInInnerCylinder(&detBlocInnerCylinder, 
													(Lb2D_Point*)&(thePhotonPtr->location))) {
							newZone = itsZone;
							if (detBlocGetZone(&(thePhotonPtr->location), &newRing, &newZone)) {
								/* Photon has been projected to its new (r,z) location */
							}
							else {
								/* Possible roundoff problem? */
								PhgAbort("Wasn't in previous ring (detBlocNextBlock).", false);
							}
						}
						else {
							/* Photon entered the inner cylinder between rings */
							photonEvent = detBlocEvt_InnerCyl;
						}
					}
					else {
						/* Photon exited the outer cylinder between rings */
						photonEvent = detBlocEvt_OuterCyl;
					}
				}
				break;
			
			case detBlocEvt_NextZone:
				newRing = itsRing;
				newZone = itsZone + 1;
				if (newZone >= detBlocNumRingZones) {
					newZone = 0;
				}
				thePhotonPtr->location = projPosition;
				thePhotonPtr->travel_distance += shortestDistance;
				break;
			
			case detBlocEvt_PrevZone:
				newRing = itsRing;
				newZone = itsZone - 1;
				if (newZone < 0) {
					newZone = detBlocNumRingZones-1;
				}
				thePhotonPtr->location = projPosition;
				thePhotonPtr->travel_distance += shortestDistance;
				break;
			
			default:
				/* Other events are not affected */
				break;
		}
		
		
		/* Handle the photon event that occurred */
		switch (photonEvent) {
			case detBlocEvt_OutEnd:
			case detBlocEvt_OuterCyl:
			case detBlocEvt_InnerCyl:
				/* Photon exited the detector */
				blockResult = false;
				break;
			
			case detBlocEvt_NextRing:
			case detBlocEvt_PrevRing:
			case detBlocEvt_NextZone:
			case detBlocEvt_PrevZone:
				/* Photon went into another zone; continue the tracking there */
				blockResult = detBlocNextBlock( thePhotonPtr, 
												newRing, newZone, NULL, 
												nextBlockPtr, 
												blockRing, blockZone, blockIndex);
				break;
			
			case detBlocEvt_Block:
				/* Photon hit a block within the current zone */
				blockResult = true;
				
				/* Project to the block boundary */
				thePhotonPtr->location = projPosition;
				thePhotonPtr->travel_distance += shortestDistance;
				
				*blockRing = itsRing;
				*blockZone = itsZone;
				break;
			
			default:
				/* Some kind of inappropriate error occurred */
				PhgAbort("Invalid photon event in detBlocNextBlock.", false);
				break;
		}
	}
	
	return (blockResult);
}


/*********************************************************************************
*
*		Name:			DetBlocProjectToDetector
*
*		Summary:		Project the photon to the nearest detector block edge.
*
*		Arguments:
*			PHG_TrackingPhoton 	*photonPtr	- The photon to track and detect.
*			LbFourByte			*ringNum	- The ring projected to.
*			LbFourByte			*detNum		- The detector projected to.
*
*		Function return: 	True if photon not rejected (still valid).
*
*********************************************************************************/

Boolean DetBlocProjectToDetector(PHG_TrackingPhoton *photonPtr, 
									LbFourByte *ringNum, LbFourByte *detNum)

{
	CylPosCylinderTy	innerCylinder;		/* Collimator outer cylinder */
	DetBlockTomoRingTy	*ringPtr;			/* Pointer to current ring info */
	Boolean				valid;				/* True if photon not rejected */
	LbFourByte			itsRing;			/* Photon starting ring */
	LbFourByte			itsZone;			/* Photon starting zone */
	
	
	/* Create a circular cylinder for the inner cylinder */
	detBlocGetInnerCylinder(&innerCylinder);
	ringPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[0]);
	innerCylinder.zMin = ringPtr->MinZ + ringPtr->AxialShift;
	ringPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[detBlocNumRings-1]);
	innerCylinder.zMax = ringPtr->MaxZ + ringPtr->AxialShift;
	
	/* Project the photon, validly, to the inner cylinder */
	if (detBlocPtInInnerCylinder(&detBlocInnerCylinder, 
									(Lb2D_Point*)(&(photonPtr->location)))) {
		/* Project the photon to the inner cylinder */
		
		PHG_Position	newPos;					/* Projected position */
		double			distance;				/* Projected distance */
		
		
		if (CylPosProjectToCylinder(&(photonPtr->location),
				&(photonPtr->angle), &innerCylinder,
				&newPos, &distance)) {
			
			if ((newPos.z_position < innerCylinder.zMin) || 
					(newPos.z_position > innerCylinder.zMax)) {
				/* Photon left detector axially */
				valid = false;
			}
			else {
				/* Update the position */
				photonPtr->location = newPos;
				photonPtr->travel_distance += distance;
				
				valid = true;
			}
		}
		else {
			/* Photon is only moving axially, leaving the detector */
			valid = false;
		}
	}
	else {
		/* Radially outside of inner cylinder already */
		if ((photonPtr->location.z_position < innerCylinder.zMin) || 
				(photonPtr->location.z_position > innerCylinder.zMax)) {
			/* But outside detector axial limits */
			valid = false;
		}
		else {
			/* Everything is okay so far */
			valid = true;
		}
	}
	
	
	if (valid) {
		/* Make sure the photon is inside a ring */
		
		Boolean			inRings;			/* True if photon inside a ring */
		
		itsRing = -1;
		itsZone = -1;
		inRings = detBlocGetZone(&(photonPtr->location), &itsRing, &itsZone);
		if (! inRings) {
			/* Not within a ring; must project until is */
			valid = detBlocProjAcrossGap(photonPtr);
			if (valid) {
				/* Get its new ring and zone */
				if (! detBlocGetZone(&(photonPtr->location), &itsRing, &itsZone)) {
					/* Possible roundoff problem? */
					PhgAbort("Wasn't in projected ring (DetBlocProjectToDetector).", false);
				}
			}
		}
	}
	
	if (valid) {
		/* Make sure the photon is still inside the outer detector cylinder */
		if (! detBlocPtInOuterCylinder((Lb2D_Point*)&(photonPtr->location))) {
			valid = false;
		}
	}
	
	
	if (valid) {
		/* Project the photon to the nearest block, if any */
		
		DetectorBlockPtr	nearestBlockPtr;	/* Nearest detector block */
		LbFourByte			blockRing;			/* Nearest detector block ring */
		LbFourByte			blockZone;			/* Ignored */
		LbFourByte			blockIndex;			/* Ignored */
		
		
		valid = detBlocNextBlock(photonPtr, itsRing, itsZone, NULL, 
									&nearestBlockPtr, &blockRing, &blockZone, &blockIndex);
		
		if (valid) {
			/* Return the nearest block's (r,b) coordinates */
			*ringNum = blockRing;
			*detNum = nearestBlockPtr->itsBlock;
		}
	}
	
	
	return (valid);
}


/*********************************************************************************
*
*		Name:			DetBlocFindNextInteraction
*
*		Summary:		Find and move the photon to the next interaction point.
*
*		Arguments:
*			PHG_TrackingPhoton 		*photonPtr		- The photon to track and detect.
*			LbFourByte				*curRingNum		- The current ring number.
*			LbFourByte				*curDetNum		- The current detector number.
*			detEn_ActionTy			*actionType		- The type of interaction that occurred.
*			double					*fpToGo			- Remaining free paths to travel.
*			LbUsFourByte			*detMaterial	- The interaction material of the detector.
*			Boolean					*isActive		- Active status of the interaction material.
*
*		Function return:	None.
*
*********************************************************************************/

void DetBlocFindNextInteraction(PHG_TrackingPhoton *photonPtr, 
								LbFourByte *curRingNum, LbFourByte *curDetNum, 
								detEn_ActionTy *actionType, double *fpToGo, 
								LbUsFourByte *detMaterial, Boolean *isActive)

{
	LbFourByte			curRing;		/* Current ring */
	DetectorBlockPtr	curBlockPtr;	/* Current block */
	DetBlocEvents		exitSide;		/* How photon exited the block */
	PHG_Position		blockExitPoint;	/* Where photon exited the block */
	double				maxTravelDist;	/* Distance to block exit point */
	double				travelDist;		/* Distance actually traveled */
	double				fpUsed;			/* Free paths used in travelDist */
	LbUsFourByte		lastMaterial;	/* Material at last point */
	Boolean				lastActive;		/* Active status at last point */
	detEn_ActionTy		detAction;		/* Photon action type */
	LbFourByte			curZone;		/* Current zone */
	DetectorBlockPtr	nextBlockPtr;	/* New block */
	LbFourByte			nextRing;		/* New block ring */
	LbFourByte			nextZone;		/* Ignored */
	LbFourByte			nextBlock;		/* Ignored */
	
	
	/* Get the current photon position's block */
	curRing = *curRingNum;
	detBlocGetNumberedBlock(curRing, *curDetNum, &curBlockPtr);
	
	if (! curBlockPtr) {
		/* The block was not found--shouldn't happen */
		exitSide = detBlocEvt_InnerCyl;
		blockExitPoint = photonPtr->location;
		maxTravelDist = 0.0;
		travelDist = 1.0;
		fpUsed = 0.0;
	}
	else {
		/* Compute the maximum travel distance within this block */
		exitSide = detBlocGetDistToExit(photonPtr, curBlockPtr, 
											&blockExitPoint, &maxTravelDist);
		
		if (exitSide == detBlocEvt_InnerCyl) {
			/* There was a failure */
			travelDist = 1.0;
			fpUsed = 0.0;
		}
		else {
			/* Determine the block distance to travel for the given free paths */
			detBlocIntraDistance(photonPtr, curBlockPtr, *fpToGo, maxTravelDist, 
									&travelDist, &fpUsed, &lastMaterial, &lastActive);
		}
	}
	
	if (travelDist <= maxTravelDist) {
		/* Photon interacted within this block */
		detAction = detEnAc_Interact;
		photonPtr->location.x_position += photonPtr->angle.cosine_x * travelDist;
		photonPtr->location.y_position += photonPtr->angle.cosine_y * travelDist;
		photonPtr->location.z_position += photonPtr->angle.cosine_z * travelDist;
		photonPtr->travel_distance += travelDist;
		*fpToGo = 0.0;
		*detMaterial = lastMaterial;
		*isActive = lastActive;
	}
	else {
		/* Photon exited the block without interaction */
		
		photonPtr->location = blockExitPoint;
		photonPtr->travel_distance += maxTravelDist;
		*fpToGo -= fpUsed;		
		
		/* Check the type of the block exit */
		if (exitSide == detBlocEvt_PrevRing) {
			/* Photon exited the lesser side of the ring */
			
			if (curRing == 0) {
				/* Leaving first ring:  out of the detector, so discard the photon */
				detAction = detEnAc_Discard;
			}
			else {
				detAction = detEnAc_AxialCross;
				curRing--;
				
				/* Project the photon across the ring gap */
				{
					DetBlockTomoRingTy*	nextRingPtr;	/* Ring info for the adjacent ring */
					double				nextRingZ;		/* Adjacent ring axial boundary */
					
					nextRingPtr = 
						&(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[curRing]);
					nextRingZ = nextRingPtr->MaxZ + nextRingPtr->AxialShift;
					detBlocProjAxially(photonPtr, nextRingZ);
				}
			}
		}
		else if (exitSide == detBlocEvt_NextRing) {
			/* Photon exited the greater side of the ring */
			
			if (curRing >= (LbFourByte)(DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings-1)) {
				/* Leaving last ring:  out of the detector, so discard the photon */
				detAction = detEnAc_Discard;
			}
			else {
				detAction = detEnAc_AxialCross;
				curRing++;
				
				/* Project the photon across the ring gap */
				{
					DetBlockTomoRingTy*	nextRingPtr;	/* Ring info for the adjacent ring */
					double				nextRingZ;		/* Adjacent ring axial boundary */
					
					nextRingPtr = 
						&(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[curRing]);
					nextRingZ = nextRingPtr->MinZ + nextRingPtr->AxialShift;
					detBlocProjAxially(photonPtr, nextRingZ);
				}
			}
		}
		else if (exitSide == detBlocEvt_Block) {
			/* Photon left the block within the current ring */
			detAction = detEnAc_LayerCross;
		}
		else {
			/* The photon was not in the block and could not enter and exit it */
			/* NOTE:  Shouldn't ever occur */
			detAction = detEnAc_Discard;
			#ifdef PHG_DEBUG
				PhgAbort("Photon failed to exit block properly (DetBlocFindNextInteraction).", true);
			#endif
		}
		
		if (detAction != detEnAc_Discard) {
			/* Make sure the photon is still inside the outer detector cylinder */
			if (! detBlocPtInOuterCylinder((Lb2D_Point*)&(photonPtr->location))) {
				detAction = detEnAc_Discard;
			}
			else {
				/* Make sure the photon is not inside the inner detector cylinder */
				if (detBlocPtInInnerCylinder(&detBlocInnerCylinder, 
												(Lb2D_Point*)&(photonPtr->location))) {
					detAction = detEnAc_Discard;
				}
			}
		}
		
		if (detAction != detEnAc_Discard) {
			/* Get the photon's new (r,z) coordinates */
			curZone = -1;
			if (! detBlocGetZone(&blockExitPoint, &curRing, &curZone)) {
				/* Possible roundoff problem? */
				PhgAbort("Wasn't in projected ring (DetBlocFindNextInteraction).", false);
			}
			
			/* Don't try to re-enter the current block, if exiting it */
			if (detAction != detEnAc_LayerCross) {
				curBlockPtr = NULL;
			}
			
			/* Determine next block */
			if (detBlocNextBlock(photonPtr, curRing, curZone, curBlockPtr, 
									&nextBlockPtr, 
									&nextRing, &nextZone, &nextBlock)) {
				*curRingNum = nextRing;
				*curDetNum = nextBlockPtr->itsBlock;
			}
			else {
				/* Exited the detector without hitting a block */
				detAction = detEnAc_Discard;
			}
		}
	}
	
	/* Return the action */
	*actionType = detAction;
}


/*********************************************************************************
*
*		Name:			DetBlocFindDetPosition
*
*		Summary:		Find the photon's detected position.
*
*		Arguments:
*			PHG_TrackingPhoton		*photonPtr			- The photon.
*			LbUsFourByte 			interactionIndex	- The number of interactions.
*
*		Function return: None.
*
*********************************************************************************/

void DetBlocFindDetPosition(PHG_TrackingPhoton *photonPtr, LbUsFourByte interactionIndex)

{
	LbUsFourByte	i;								/* Index through interactions */
	Boolean			processed[MAX_DET_INTERACTIONS];/* Which blocks have been handled */
	LbUsFourByte	curInt;							/* Current interaction */
	Boolean			interactsSet[MAX_DET_INTERACTIONS];	/* Which interactions to use */
	detInteractionInfoTy	
					*curInteractPtr;				/* Pointer to current interaction */
	LbUsFourByte	curRing;						/* Ring of current interaction */
	LbUsFourByte	curBlock;						/* Block of current interaction */
	PHG_Position	centroidPos;					/* The computed centroid position */
	double			depositedEnergy;				/* The amount of energy deposited */
	LbUsFourByte	detRing;						/* Ring of detected position */
	LbUsFourByte	detBlock;						/* Block of detected position */
	
	
	for (i=0; i<=interactionIndex; i++) {
		processed[i] = false;
	}
	photonPtr->energy = 0.0;
	
	/* Loop through the different sets of interactions */
	for (curInt=0; curInt<=interactionIndex; curInt++) {
		if (! processed[curInt]) {
			/* Indicate the requested set of interactions */
			for (i=0; i<curInt; i++) {
				interactsSet[i] = false;
			}
			curInteractPtr = &(photonPtr->det_interactions[curInt]);
			curRing = curInteractPtr->posIndices.ringNum;
			curBlock = curInteractPtr->posIndices.blockNum;
			interactsSet[curInt] = true;
			processed[curInt] = true;
			for (i=curInt+1; i<=interactionIndex; i++) {
				curInteractPtr = &(photonPtr->det_interactions[i]);
				if ((curInteractPtr->posIndices.ringNum == curRing) && 
						(curInteractPtr->posIndices.blockNum == curBlock)) {
					interactsSet[i] = true;
					processed[i] = true;
				}
				else {
					interactsSet[i] = false;
				}
			}
			
			/* Compute the centroid and deposited energy for the requested set */
			DetGeomCompCentroid(photonPtr, 
									&centroidPos, &depositedEnergy, 
									interactionIndex, interactsSet);
			
			/* Save the block data having the greatest deposited energy */
			if (depositedEnergy >= photonPtr->energy) {
				photonPtr->energy = depositedEnergy;
				photonPtr->detLocation = centroidPos;
				detRing = curRing;
				detBlock = curBlock;
			}
		}
	}
	
	/* Find the adjusted active crystal centroid position and crystal number */
	centroidPos = photonPtr->detLocation;
	if (! detBlocIntraFindCentroid(&centroidPos, detRing, detBlock, 
									&(photonPtr->detLocation), &(photonPtr->detCrystal))) {
		/* There was a problem */
		PhgAbort("Could not adjust the centroid to an active crystal.", false);
	}
	if (DetRunTimeParams[DetCurParams].BlockTomoDetector.BlockDetectedPositionAlgo == 
			DetEn_BlockUseCentroid) {
		/* Restore the original centroid position (but still record the nearest crystal) */
		photonPtr->detLocation = centroidPos;
	}
}


/*********************************************************************************
*
*		Name:			detBlocIntraSearchCentroidLayer
*
*		Summary:		Given a fixed position, search a given (ring, block, layer) 
*							for an active element in the layer whose center is 
*							closest to the given fixed position.
*						Fill in the element number in the element index.
*						Return the distance.
*						All parameter positions are in block coordinates.
*
*		Arguments:
*			PHG_Position			*fixedPositionPtr	- Fixed position.
*			DetectorBlockPtr		theBlockPtr			- Pointer to the block.
*			DetElementPosIndexTy	*closestElemDataPtr	- Element index data.
*			double					*closestDistance	- Closest distance.
*
*		Function return: 	True if found, false for layer with no active element.
*
*********************************************************************************/

Boolean detBlocIntraSearchCentroidLayer(	PHG_Position			*fixedPositionPtr, 
											DetectorBlockPtr		theBlockPtr, 
											DetElementPosIndexTy	*closestElemDataPtr, 
											double					*closestDistance)

{
	Boolean					valid;				/* Return value */
	DetElementPosIndexTy	testElemIndex;		/* Position index of a test element */
	DetBlockTomoLayerTy		*layerPtr;			/* Pointer to the given layer */
	LbUsFourByte			numLayerElems;		/* Number of layer elements */
	double					closestDistSqrd;	/* Squared distance to the closest element center */
	LbUsFourByte			elem;				/* Index through layer elements */
	PHG_Position			testCenter;			/* Current test element center */
	double					testDistSqrd;		/* Squared distance to test element center */
	LbUsFourByte			closestElement;		/* Index of the closest element */
	
	
	/* NOTE:  The layer is assumed to exist. */
	
	/* Check all active elements in the supplied (ring, block, layer) 
		for the element center that is closest to the given position */
	testElemIndex.ringNum = closestElemDataPtr->ringNum;
	testElemIndex.blockNum = closestElemDataPtr->blockNum;
	testElemIndex.layerNum = closestElemDataPtr->layerNum;
	layerPtr = &(theBlockPtr->userBlockData->LayerInfo[testElemIndex.layerNum]);
	numLayerElems = layerPtr->NumElements;
	valid = false;
	closestDistSqrd = 1000000000.0;		/* Really big distance (squared) */
	for (elem=0; elem<numLayerElems; elem++) {
		testElemIndex.elementNum = elem;
		if (layerPtr->ElementInfo[elem].IsActive) {
			if (detBlocGetElementCenter(&testElemIndex, &testCenter)) {
				/* Determine the distance (squared) */
				testDistSqrd = 	PHGMATH_Square(testCenter.x_position - 
													fixedPositionPtr->x_position) + 
								PHGMATH_Square(testCenter.y_position - 
													fixedPositionPtr->y_position) + 
								PHGMATH_Square(testCenter.z_position - 
													fixedPositionPtr->z_position);
				
				if (testDistSqrd < closestDistSqrd) {
					/* Closer element center */
					closestDistSqrd = testDistSqrd;
					closestElement = elem;
					valid = true;
				}
			}
			else {
				/* No such element; shouldn't actually happen */
				int		i;
				i = 10;		/* Breakpoint location */
			}
		}
	}
	
	if (valid) {
		/* Return the closest element */
		closestElemDataPtr->elementNum = closestElement;
		*closestDistance = PHGMATH_SquareRoot(closestDistSqrd);
	}
	
	return (valid);
}


/*********************************************************************************
*
*		Name:			detBlocIntraFindCentroid
*
*		Summary:		Given a centroid position in a given ring and block, 
*							find the new position for the centroid centered within 
*							the closest active element in the block and its 
*							active crystal identification number.
*						All parameter positions are in tomograph coordinates.
*
*		Arguments:
*			PHG_Position		*origCentroidPtr	- Original centroid position.
*			LbFourByte			itsRing				- Ring containing the position.
*			LbFourByte			itsBlock			- Block containing the position.
*			PHG_Position		*newCentroidPtr		- New centroid position.
*			LbFourByte			*crystalNum			- Active crystal number.
*
*		Function return: 	True if found, false for block with no active element.
*
*********************************************************************************/

Boolean detBlocIntraFindCentroid(	PHG_Position		*origCentroidPtr, 
									LbFourByte			itsRing, 
									LbFourByte			itsBlock, 
									PHG_Position		*newCentroidPtr, 
									LbFourByte			*crystalNum)

{
	Boolean					valid;				/* Return value */
	DetectorBlockPtr		theBlockPtr;		/* Block containing the centroid */
	DetElementPosIndexTy	thePosIndex;		/* Position index of the centroid */
	DetBlockTomoElementTy	*centElemDataPtr;	/* Data for the centroid's element */
	
	
	/* Get the block containing the centroid */
	valid = detBlocGetNumberedBlock(itsRing, itsBlock, &theBlockPtr);
	if (valid) {
		thePosIndex.ringNum = itsRing;
		thePosIndex.blockNum = itsBlock;
	}
	
	if (valid) {
		/* Check for single or multiple elements */
		if ((theBlockPtr->userBlockData->NumLayers == 1) && 
				(theBlockPtr->userBlockData->LayerInfo[0].NumElements == 1)) {
			/* Only one element in the block */
			
			/* Get the element containing the centroid */
			valid = detBlocGetElementIndex(origCentroidPtr, &thePosIndex);
			
			if (valid) {
				centElemDataPtr = &(theBlockPtr->userBlockData->LayerInfo[thePosIndex.layerNum
														].ElementInfo[thePosIndex.elementNum]);
				if (centElemDataPtr->IsActive) {
					/* Find its center, which is the center of the block
						(only two diagonal points are needed, but using all four is more robust) */
					newCentroidPtr->x_position = (theBlockPtr->itsRect.corner_1.x_position + 
												theBlockPtr->itsRect.corner_2.x_position + 
												theBlockPtr->itsRect.corner_3.x_position + 
												theBlockPtr->itsRect.corner_4.x_position) / 4.0;
					newCentroidPtr->y_position = (theBlockPtr->itsRect.corner_1.y_position + 
												theBlockPtr->itsRect.corner_2.y_position + 
												theBlockPtr->itsRect.corner_3.y_position + 
												theBlockPtr->itsRect.corner_4.y_position) / 4.0;
					newCentroidPtr->z_position = (theBlockPtr->zMin + theBlockPtr->zMax) / 2.0;
					
					*crystalNum = centElemDataPtr->crystalNumInTomo;
				}
				else {
					valid = false;
				}
			}
		}
		else {
			/* Multiple elements in the block */
			
			PHG_Position			origBlockCentroid;	/* Original centroid in block coords */
			DetBlockTomoElementTy	*centroidElementPtr;/* Pointer to chosen active element */
			DetElementPosIndexTy	centroidIndex;		/* Position index of the chosen element */
			double					centroidDist;		/* Distance to the chosen element center */
			LbFourByte				curLayer;			/* Current layer to search through */
			DetElementPosIndexTy	testIndex;			/* Position index of a test layer/element */
			double					centroidDistPlus;	/* Distance to the closest upper center */
			double					centroidDistMinus;	/* Distance to the closest lower center */
			PHG_Position			newBlockCentroid;	/* New centroid in block coords */
			
			
			/* Convert the centroid to block coordinates */
			detBlocTomoToBlockCoord(theBlockPtr->userBlockData, origCentroidPtr, &origBlockCentroid);
			
			/* Since the centroid really is in the block, make sure roundoff didn't violate that */
			if (origBlockCentroid.x_position < theBlockPtr->userBlockData->XMin ) {
				origBlockCentroid.x_position = theBlockPtr->userBlockData->XMin;
			}
			else if (origBlockCentroid.x_position > theBlockPtr->userBlockData->XMax ) {
				origBlockCentroid.x_position = theBlockPtr->userBlockData->XMax;
			}
			if (origBlockCentroid.y_position < theBlockPtr->userBlockData->YMin ) {
				origBlockCentroid.y_position = theBlockPtr->userBlockData->YMin;
			}
			else if (origBlockCentroid.y_position > theBlockPtr->userBlockData->YMax ) {
				origBlockCentroid.y_position = theBlockPtr->userBlockData->YMax;
			}
			if (origBlockCentroid.z_position < theBlockPtr->userBlockData->ZMin ) {
				origBlockCentroid.z_position = theBlockPtr->userBlockData->ZMin;
			}
			else if (origBlockCentroid.z_position > theBlockPtr->userBlockData->ZMax ) {
				origBlockCentroid.z_position = theBlockPtr->userBlockData->ZMax;
			}
			
			/* Get the element containing the centroid */
			valid = detBlocGetElementIndex(&origBlockCentroid, &thePosIndex);
			
			if (valid) {
				centroidElementPtr = NULL;
				
				centElemDataPtr = &(theBlockPtr->userBlockData->LayerInfo[thePosIndex.layerNum
														].ElementInfo[thePosIndex.elementNum]);
				if (centElemDataPtr->IsActive) {
					/* The original centroid is already in the desired active element */
					centroidElementPtr = centElemDataPtr;
					centroidIndex = thePosIndex;
				}
				else {
					if ( theBlockPtr->userBlockData->ActiveLayers[thePosIndex.layerNum] ) {
						/* Get the closest active element within the centroid's layer */
						if (detBlocIntraSearchCentroidLayer(&origBlockCentroid, theBlockPtr, 
																&thePosIndex, &centroidDist)) {
							centroidElementPtr = &(theBlockPtr->userBlockData->LayerInfo[
													thePosIndex.layerNum].ElementInfo[
													thePosIndex.elementNum]);
							centroidIndex = thePosIndex;
						}
						else {
							/* There were no active elements in the layer; shouldn't happen */
							int		i;
							i = 10;		/* Breakpoint location */
						}
					}
					else {
						/* Check the elements in the two closest layers on either side (in x)
							for the closest active centroid element center */
						
						/* Layer above (in x) */
						curLayer = thePosIndex.layerNum;
						do {
							curLayer++;
							if ( theBlockPtr->userBlockData->ActiveLayers[curLayer] ) {
								/* The current layer has active elements */
								testIndex = thePosIndex;
								testIndex.layerNum = curLayer;
								if (detBlocIntraSearchCentroidLayer(&origBlockCentroid, theBlockPtr, 
																	&testIndex, &centroidDistPlus)) {
									centroidElementPtr = &(theBlockPtr->userBlockData->LayerInfo[
														curLayer].ElementInfo[testIndex.elementNum]);
									centroidIndex = testIndex;
								}
								else {
									/* There were no active elements in the layer; shouldn't happen */
									int		i;
									i = 10;		/* Breakpoint location */
								}
								break;
							}
						} while (true);
						
						/* Layer below (in x) */
						curLayer = thePosIndex.layerNum;
						do {
							curLayer--;
							if ( theBlockPtr->userBlockData->ActiveLayers[curLayer] ) {
								/* The current layer has active elements */
								testIndex = thePosIndex;
								testIndex.layerNum = curLayer;
								if (detBlocIntraSearchCentroidLayer(&origBlockCentroid, theBlockPtr, 
																	&testIndex, &centroidDistMinus)) {
									if (centroidDistMinus < centroidDistPlus) {
										centroidElementPtr = &(theBlockPtr->userBlockData->LayerInfo[
															curLayer].ElementInfo[testIndex.elementNum]);
										centroidIndex = testIndex;
									}
								}
								else {
									/* There were no active elements in the layer; shouldn't happen */
									int		i;
									i = 10;		/* Breakpoint location */
								}
								break;
							}
						} while (true);
					}
				}
				
				if (centroidElementPtr) {
					/* Get the center of the chosen active element */
					detBlocGetElementCenter(&centroidIndex, &newBlockCentroid);
					
					/* Convert the block centroid back to tomograph coordinates and save it */
					detBlocBlockToTomoCoord(theBlockPtr->userBlockData, 
												&newBlockCentroid, newCentroidPtr);
					
					*crystalNum = centroidElementPtr->crystalNumInTomo;
				}
				else {
					/* No active element was found */
					valid = false;
				}
			}
		}
	}
	
	return (valid);
}


/*********************************************************************************
*
*		Name:			DetBlocFreeData
*
*		Summary:		Free any dynamically allocated memory.
*
*		Arguments:		None.
*
*		Function return: None.
*
*********************************************************************************/

void DetBlocFreeData(void)

{
	LbMmFree((void **)&detBlocBlocksDatabase);
	LbMmFree((void **)&detBlocZoneBounds);
	LbMmFree((void **)&detBlocIndexDatabase);
}
