/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1992-2011 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		CylPos.c
*			Revision Number:	1.1
*			Date last revised:	8 December 2011
*
*			Programmer:			Steven Vannoy
*			Date Originated:	Wednesday, September 9, 1992
*
*			Module Overview:	Cylinder Positions routines.
*
*			References:			'Cylinder Position Processes' PHG design.
*
**********************************************************************************
*
*			Global functions defined:
*				CylPosCalcDistanceToCylinder
*				CylPosCalcDistanceToLimitSurface
*				CylPosCalcDistanceToObjectSurface
*				CylPosCalcForcedAngAndWeight
*				CylPosClipToLimitCylinder
*				CylPosGetLimitZRange
*				CylPosGetTargetZRange
*				CylPosGetObjectRadius
*				CylPosGetTargetRadius
*				CylPosInitCriticalZone
*				CylPosInitLimitCylinder
*				CylPosInitObjectCylinder
*				CylPosInitTargetCylinder
*				CylPosIsOutsideObjCylinder
*				CylPosIsOutsideTarget
*				CylPosProjectToTargetCylinder
*				CylPosProjectToCylinder
*				CylPosWillIntersectCritZone
*				CylPosDumpObjects
*
*			Global Macros
*				CYLPOSGetObjZMin
*				CYLPOSGetObjZMax
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
*			Programmer(s):	R Harrison	
*
*			Revision date:	21 Nov 2006
*
*			Revision description:	changed CylPosFind2dIntersection to
*					use PhgMathSolveQuadratic.
*
*********************************************************************************/
#define CYLINDER_POSITION

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbFile.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ProdTbl.h"
#include "SubObj.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "UNCCollimator.h"
#include "PhgMath.h"
#include "ColUsr.h"
#include "CylPos.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "PhoHFile.h"
#include "phg.h"
#include "PhgBin.h"

						/*		#define DUMP_INTERSECTION */
						/*		#define DUMP_CLIPS */
								
								
/* LOCAL CONSTANTS */
#define	CYLPOS_TARGETFILE	"phgcylpos.target"			/* Name of target cylinder description file */

/* LOCAL TYPES */

/* LOCAL GLOBALS */
char	cylPosErrString[1024];

#ifdef CYLPOS_DEBUG
Boolean				CylPosTargetInitialized = false;	/* Have we initialized? */
#endif

#ifdef DUMP_INTERSECTION
FILE				*CylPosIntersectFile;		/* x-y intersections of critical zone */
#endif
#ifdef DUMP_CLIPS
FILE				*CylPosClipFile;			/* clip positions */
#endif

/* PROTOTYPES */
/* LOCAL MACROS */
/*********************************************************************************
*
*			Name:			CYLPOSGetTargetCylAxialLength
*
*			Summary:		Returns the axial length of the target cylinder.
*			Arguments:
*
*			Function return: double: the axial length of the target cylinder.
*
*********************************************************************************/
#define CYLPOSGetTargetCylAxialLength() (fabs(CylPosTargetCylinder.zMax - CylPosTargetCylinder.zMin))

/*********************************************************************************
*
*			Name:			CYLPOSGetCriticalZoneAxialLength
*
*			Summary:		Returns the axial length of the critical zone.
*			Arguments:
*
*			Function return: double: the axial length of the critical zone.
*
*********************************************************************************/
#define CYLPOSGetCriticalZoneAxialLength() (fabs(CylPosCriticalZone.zMax - CylPosCriticalZone.zMin))

/*********************************************************************************
*
*			Name:			CYLPOSGetObjectCylAxialLength
*
*			Summary:		Returns the axial length of the object cylinder.
*			Arguments:
*
*			Function return: double: the axial length of the object cylinder.
*
*********************************************************************************/
#define CYLPOSGetObjectCylAxialLength() (fabs(CylPosObjectCylinder.zMax - CylPosObjectCylinder.zMin))

/*********************************************************************************
*
*			Name:			CYLPOSGetLimitCylAxialLength
*
*			Summary:		Returns the axial length of the limit cylinder.
*			Arguments:
*
*			Function return: double: the axial length of the limit cylinder.
*
*********************************************************************************/
#define CYLPOSGetLimitCylAxialLength() (fabs(CylPosLimitCylinder.zMax - CylPosLimitCylinder.zMin))



/* FUNCTIONS */

/*********************************************************************************
*
*			Name:			CylPosInitTargetCylinder
*
*			Summary:	Get the user-chosen values for the target cylinder.
*
*			Arguments:
*				LbPfHkTy		paramFlHk	- The parameter file hook.
*				LbUsFourByte	numParams	- The number of parameters.
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean CylPosInitTargetCylinder(LbPfHkTy paramFlHk, LbUsFourByte numParams)	
{
	double					paramBuffer[LBPF_PARAM_LEN];	/* Storage for params */
	Boolean					okay = false;					/* Process flag */
	Boolean					isEOF;							/* End of file flag */
	Boolean					switchOkay;						/* Switch flag */
	LbPfEnPfTy				paramType;						/* Type of parameter read */
	LbUsFourByte			paramSize;						/* Size of parameter read */
	char					paramLabel[LBPF_LABEL_LEN];		/* Label for parameter */
	char					errString[256];					/* Storage for error strings */
	PhgEn_RunTimeParamsTy	whichParam;						/* Enum of parameter read */
	LbUsFourByte			paramIndex;						/* Current parameter */

	do {	/* Process Loop */
	
	
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

				case CylPosEn_zMin:
					CylPosTargetCylinder.zMin = 
						*((double *) paramBuffer);
					break;

				case CylPosEn_zMax:
					CylPosTargetCylinder.zMax = 
						*((double *) paramBuffer);
					break;

				case CylPosEn_centerX:
					CylPosTargetCylinder.centerX = 
						*((double *) paramBuffer);
					break;

				case CylPosEn_centerY:
					CylPosTargetCylinder.centerY = 
						*((double *) paramBuffer);
					break;

				case CylPosEn_radius:
					CylPosTargetCylinder.radius = 
						*((double *) paramBuffer);
					break;

				
				default:
					sprintf(errString, "(CylPosInitTargetCylinder) Unknown (hence unused) parameter (%s).\n",
						paramLabel);
					ErAlert(errString, false);
					break;
			}
			if (!switchOkay)
				break;
		}
		#ifdef CYLPOS_DEBUG
			CylPosTargetInitialized = true;
		#endif
		
		okay = true;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			CylPosFind2dIntersection
*			Summary:	Calculate the 2 dimensional intersection of a line with
*						a cylinder. All coordinates are assumed to be normalized
*						to coordinate system origin.
*			Arguments:
*				double	x0		- Initial position on X axis.
*				double	y0		- Initial position on Y axis.
*				double	cosX	- Direction component in X direction.
*				double	cosY	- Direction component in Y direction.
*				double	radius	- Radius of target cylinder.
*				double	*d1		- Shortest distance to intercept.
*				double	*d2		- Longest distance intercept.
*
*			Function return: True if there is an intersection.
*
*********************************************************************************/
Boolean CylPosFind2dIntersection(double x0, double y0, double cosX, double cosY,
			double radius, double *d1, double *d2)	
{
	Boolean	hasIntersection = false;	/* Intersection flag */
	double	a;							/* Temp for quadratic */
	double	b;							/* Temp for quadratic */
	double	c;							/* Temp for quadratic */
	LbUsFourByte numRoots;				/* PhgMathSolveQuadratic number roots */

	do { /* Process Loop */
	
		/* Compute a, b, c */
		a = PHGMATH_Square(cosX) + PHGMATH_Square(cosY);
		b = 2 * ((x0*cosX) + (y0*cosY));
		c = PHGMATH_Square(x0) + PHGMATH_Square(y0) - PHGMATH_Square(radius);
		numRoots = PhgMathSolveQuadratic( a, b, c, d1, d2 );
		
		/* See if there is no intersection */
		if (numRoots == 0)
			break;
		
		/* if only one root was found, assign d2
		 (some functions using this function look for d1==d2) */
		if (numRoots == 1) {
			*d2 = *d1;
		}
				
		hasIntersection = true;
	} while (false);
	
	return (hasIntersection);
}


void OldCylPosCalcDistanceToCylSurface(PHG_Position *posPtr, PHG_Direction *dirPtr,
		CylPosCylinderTy *cylinderPtr, double *distPtr);
/*********************************************************************************
*
*			Name:			CylPosCalcDistanceToCylSurface
*			Summary:	Calculate the distance from a starting location/direction
*						to the surface of the cylinder.
*			Arguments:
*				PHG_Position		posPtr			- The position.
*				PHG_Direction		dirPtr			- The direction.
*				CylPosCylinderTy	*cylinderPtr	- Radius of the cylinder.
*				double				*distPtr		- The distance.
*
*			Function return: None.
*
*********************************************************************************/
void CylPosCalcDistanceToCylSurface(PHG_Position *posPtr, PHG_Direction *dirPtr,
		CylPosCylinderTy *cylinderPtr, double *distPtr)	
{
	double			xCord;					/* X coordinate adjusted for center of cylinder */
	double			yCord;					/* Y coordinate adjusted for center of cylinder */
	double			a, b, c;				/* Quadratic formula variables */
	double			minRoot, maxRoot;		/* Roots returned by quadratic formula */
	LbUsFourByte	numRoots;				/* number of roots returned by quadratic formula */
	


	/* Calculate "distance" using PhgMath quadratic routine */
	{
		/* Compute x/y coordinates adjusted for center of cylinder */
		xCord = posPtr->x_position - cylinderPtr->centerX;
		yCord = posPtr->y_position - cylinderPtr->centerY;
		
		/* Calculate a, b, and c of the quadratic formula */
		a = 1 - PHGMATH_Square(dirPtr->cosine_z);
		b = 2 * ( (xCord * dirPtr->cosine_x) + (yCord * dirPtr->cosine_y) );
		c = PHGMATH_Square(xCord) + PHGMATH_Square(yCord) - PHGMATH_Square(cylinderPtr->radius);
		
		numRoots = PhgMathSolveQuadratic(a, b, c, &minRoot, &maxRoot);
		

		#ifdef CYLPOS_DEBUG
		/*  These debug tests look to see if this routine is being called at inappropriate times */
		{
			if (numRoots == 0){
				PhgAbort("CylPosCalcDistanceToCylSurface:  No roots in distance calculation.",
								true);
			} else if ( (numRoots == 1) && (fabs(minRoot) > 0.0000001) ) {
				PhgAbort("CylPosCalcDistanceToCylSurface:  Photon outside cylinder in distance calculation.",
								true);
			} else if ( (numRoots == 2) && (minRoot > 0.0000001) ) {
				PhgAbort("CylPosCalcDistanceToCylSurface:  Photon outside cylinder in distance calculation.",
								true);
			}
		}
		#endif
		
		/* Now distance */
		if (numRoots == 2){
			*distPtr = maxRoot;
		} else {
			*distPtr = minRoot;
		}

/* #define TEMP_DEBUG */
		#ifdef TEMP_DEBUG
		/* This debug test compares the distance obtained by this version of CylPosCalcDistanceToCylSurface
		with the original version */
		{
			PHG_Position TEMPpos;
			PHG_Direction TEMPdir;
			CylPosCylinderTy TEMPCylinder;
			double TEMPdist;
			double projX, projY, projOldX, projOldY;	/* x,y coordinates projected the distance
														given by CylPosCalcDistanceToCylSurface */
			double newDifSquare, oldDifSquare;	/* difference between the radius of the cylinder
												and the radius of the photon projected distance
												given by CylPosCalcDistanceToCylSurface */
			
			TEMPpos = *posPtr;
			TEMPdir = *dirPtr;
			TEMPCylinder = *cylinderPtr;
			TEMPdist = *distPtr;
			
			OldCylPosCalcDistanceToCylSurface(&TEMPpos, &TEMPdir, &TEMPCylinder, &TEMPdist);
			
			if (!PhgMathRealNumAreEqual(*distPtr, TEMPdist, -7, 0, 0, 0)) {
				projX = xCord + (*distPtr * dirPtr->cosine_x);
				projY = yCord + (*distPtr * dirPtr->cosine_y);
				projOldX = xCord + (TEMPdist * dirPtr->cosine_x);
				projOldY = yCord + (TEMPdist * dirPtr->cosine_y);
				newDifSquare = fabs( PHGMATH_Square(cylinderPtr->radius) -
						( PHGMATH_Square(projX) + PHGMATH_Square(projY) ) );
				oldDifSquare = fabs( PHGMATH_Square(cylinderPtr->radius) -
						( PHGMATH_Square(projOldX) + PHGMATH_Square(projOldY) ) );
				if ( newDifSquare > oldDifSquare ) {
					PhgAbort("CylPosCalcDistanceToCylSurface:  Old algorithm works better.",
								true);
				}
			}
		}
		#endif

						
	}	
}

/*********************************************************************************
*
*			Name:			CylPosCalcDistanceToCylSurface
*			Summary:	Calculate the distance from a starting location/direction
*						to the surface of the cylinder.
*			Arguments:
*				PHG_Position		posPtr			- The position.
*				PHG_Direction		dirPtr			- The direction.
*				CylPosCylinderTy	*cylinderPtr	- Radius of the cylinder.
*				double				*distPtr		- The distance.
*
*			Function return: None.
*
*********************************************************************************/
		
void OldCylPosCalcDistanceToCylSurface(PHG_Position *posPtr, PHG_Direction *dirPtr,
		CylPosCylinderTy *cylinderPtr, double *distPtr)	
{
	double			temp1;					/* Temporary variable */
	double			oneMinusCosSquGamma;	/* Temporary variable */
	double			temp3;					/* Temporary variable */
	double			xCord;					/* X coordinate adjusted for center of cylinder */
	double			yCord;					/* Y coordinate adjusted for center of cylinder */


	/* Calculate "distance" from Design documents*/
	{
		/* Compute x/y coordinates adjusted for center of cylinder */
		xCord = posPtr->x_position - cylinderPtr->centerX;
		yCord = posPtr->y_position - cylinderPtr->centerY;
		
		/* Calculate first temp */
		temp1 = -((xCord * dirPtr->cosine_x) +
			(yCord * dirPtr->cosine_y));
			
		/* Calculate second temp */
		oneMinusCosSquGamma = 1 - PHGMATH_Square(dirPtr->cosine_z);

		#ifdef CYLPOS_DEBUG
		{
			char	errStr[256];			/* Error string */
			if (oneMinusCosSquGamma == 0){
				sprintf(errStr,
					"oneMinusCosSquGamma should not be zero, directionPtr->cosine_z == %3.2e  (CylPosCalcDistanceToCylSurface).",
					dirPtr->cosine_z);
				PhgAbort(errStr, true);
			}
		}
		#endif
		
		/* Calculate value under the square root in quadratic formula */
		temp3 = (oneMinusCosSquGamma * PHGMATH_Square(cylinderPtr->radius)) -
			PHGMATH_Square((dirPtr->cosine_x * yCord) -
			(dirPtr->cosine_y * xCord));
						
		/* Now distance */
		*distPtr = (temp1 + PHGMATH_SquareRoot(temp3))/oneMinusCosSquGamma;
		
	}	
}
/*********************************************************************************
*
*			Name:			CylPosGetTargetZRange
*
*			Summary:	Return the range limiters for the target cylinder.
*			Arguments:
*				double	*zMinPtr	- Minimum z value.
*				double	*zMaxPtr	- Maximum z value.
*
*			Function return: None.
*
*********************************************************************************/
void CylPosGetTargetZRange(double *zMinPtr, double *zMaxPtr)	
{
	*zMinPtr = CylPosTargetCylinder.zMin;
	*zMaxPtr = CylPosTargetCylinder.zMax;
}

/*********************************************************************************
*
*			Name:			CylPosGetTargetRadius
*
*			Summary:	Return the radius of the target cylinder.
*			Arguments:
*				.
*
*			Function return: None.
*
*********************************************************************************/
double CylPosGetTargetRadius()	
{
	return (CylPosTargetCylinder.radius);
}


/*********************************************************************************
*
*			Name:			CylPosProjectToCylinder
*
*			Summary:	Given a location and angle, project to cylinder.
*			Arguments:
*				PHG_Position		*positionPtr	- The current position.
*				PHG_Direction		*directionPtr	- The current direction.
*				CylPosCylinderTy	*cylinderPtr	- The cylinder to project to.
*				PHG_Position		*newPosPtr		- The new position.
*				double				*distPtr		- Distance to projection.
*
*			Function return: TRUE unless cosine_z == +/-1.
*
*********************************************************************************/
Boolean CylPosProjectToCylinder(PHG_Position *positionPtr,
			PHG_Direction *directionPtr, CylPosCylinderTy *cylinderPtr,
			PHG_Position *newPosPtr, double *distPtr)	
{
	Boolean	intersects = false;	/* Do we intersect */
	
	do {
	
		/* If cosine_z == +/- 1, we will never intersect */
		if (PhgMathRealNumAreEqual(directionPtr->cosine_z, 1.0, -7, 0, 0, 0)) {
			break;
		}
		if (PhgMathRealNumAreEqual(directionPtr->cosine_z, -1.0, -7, 0, 0, 0)) {
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
			
		#ifdef CYLPOS_DEBUG_HEAVY
		{
			char	errString[1024];
			/* Make sure we are on the surface */
			if (!PhgMathRealNumAreEqual((PHGMATH_Square(newPosPtr->x_position-cylinderPtr->centerX) +
					PHGMATH_Square(newPosPtr->y_position-cylinderPtr->centerY)) -
					PHGMATH_Square(cylinderPtr->radius), 0.0, -7, 0, 0, 0)) {
				
				/* Build error string message */
				sprintf(errString, "\nInvalid projection from CylPosProjectToCylinder.\n\tradius = %3.2e, *distPtr= %3.2e\n\tx = %3.2e, y = %3.2e, z = %3.2e, x^2 + y ^2 = %3.2e,\n\torig_x = %3.2e, orig_y = %3.2e, orig_z = %3.2e\ncos(x) = %3.2e, cos(y) = %3.2e, cos(z) = %3.2e.",
				cylinderPtr->radius, *distPtr,
				newPosPtr->x_position, newPosPtr->y_position, newPosPtr->z_position,
				(PHGMATH_Square(newPosPtr->x_position-cylinderPtr->centerX) +
				PHGMATH_Square(newPosPtr->y_position-cylinderPtr->centerY)),
				positionPtr->x_position, positionPtr->y_position, positionPtr->z_position,
				directionPtr->cosine_x, directionPtr->cosine_y, directionPtr->cosine_z);
				
				/* Alert the user */
				PhgAbort(errString, true);
			}
		}
		#endif
		
		intersects = true;
	} while (false);
	
	return (intersects);
}

/*********************************************************************************
*
*			Name:			CylPosIsOutsideObjCylinder
*
*			Summary:	Return a boolean value indicating whether the photon is 
*						outside the object cylinder or not.
*			Arguments:
*				PHG_Position	photon_position	- The current position.
*
*			Function return: True if photon position is outside of object cylinder.
*
*********************************************************************************/
Boolean CylPosIsOutsideObjCylinder(PHG_Position	photon_position)	
{
	Boolean isOutside = true;	/* Assume photon is outside */
	
	do {	/* Process Loop */
	
		/* See if outside x,y  boundary */
		if ((PHGMATH_Square(photon_position.x_position - CylPosObjectCylinder.centerX) +
				PHGMATH_Square(photon_position.y_position - CylPosObjectCylinder.centerY)) >
				PHGMATH_Square(CylPosObjectCylinder.radius)) {
			break;
		}
		
		/* If we made it here, the photon is not outside */
		isOutside = false;
	} while (false);
		
	return (isOutside);
}


/*********************************************************************************
*
*			Name:			CylPosGetObjectRadius
*
*			Summary:	Return the radius of the object cylinder.
*			Arguments:
*				.
*
*			Function return: None.
*
*********************************************************************************/
double CylPosGetObjectRadius()	
{
	return (CylPosObjectCylinder.radius);
}

/*********************************************************************************
*
*			Name:			CylPosCalcDistanceToLimitSurface
*			Summary:	Calculate the distance from a starting location/direction
*						to the surface of the limit cylinder.
*			Arguments:
*				PHG_Position	posPtr		- The position.
*				PHG_Direction	dirPtr		- The direction.
*				double			*distPtr	- The distance.
*
*			Function return: None.
*
*********************************************************************************/
void CylPosCalcDistanceToLimitSurface(PHG_Position *posPtr, PHG_Direction *dirPtr,
		double *distPtr)	
{
	CylPosCalcDistanceToCylSurface(posPtr, dirPtr, 
	    &CylPosLimitCylinder, distPtr);
}

/*********************************************************************************
*
*			Name:			CylPosCalcDistanceToObjectSurface
*			Summary:	Calculate the distance from a starting location/direction
*						to the surface of the object cylinder.
*			Arguments:
*				PHG_Position	posPtr		- The position.
*				PHG_Direction	dirPtr		- The direction.
*				double			*distPtr	- The distance.
*
*			Function return: None.
*
*********************************************************************************/
void CylPosCalcDistanceToObjectSurface(PHG_Position *posPtr, PHG_Direction *dirPtr,
		double *distPtr)	
{
	CylPosCalcDistanceToCylSurface(posPtr, dirPtr, &CylPosObjectCylinder,
		distPtr);
		
}

/*********************************************************************************
*
*			Name:			CylPosCalcDistanceToTargetSurface
*			Summary:	Calculate the distance from a starting location/direction
*						to the surface of the target cylinder.
*			Arguments:
*				PHG_Position	posPtr		- The position.
*				PHG_Direction	dirPtr		- The direction.
*				double			*distPtr	- The distance.
*
*			Function return: None.
*
*********************************************************************************/
void CylPosCalcDistanceToTargetSurface(PHG_Position *posPtr, PHG_Direction *dirPtr,
		double *distPtr)	
{
	CylPosCalcDistanceToCylSurface(posPtr, dirPtr, &CylPosTargetCylinder,
		distPtr);
}

/*********************************************************************************
*
*			Name:			CylPosClipToLimitCylinder
*
*			Summary:	If the given position has exited the limit cylinder, 
*						then clip the position along its direction 
*						to where it leaves the limit cylinder, 
*						otherwise change nothing.
*			Arguments:
*				PHG_Position	*posPtr			- The position to clip.
*				PHG_Direction	*directionPtr	- The direction it is traveling at.
*				Boolean			*wasClipped		- Did the position need to be clipped?
*
*			Function return: None.
*
*********************************************************************************/
void CylPosClipToLimitCylinder(PHG_Position *posPtr, PHG_Direction *directionPtr)	
{
	double			distToSurface;			/* Distance to surface of cylinder */
	double			adjustmentRatio;		/* Ratio used to adjust distance to travel */
	
	/* If we are out of the radius, do the clipping */
	if (((PHGMATH_Square(posPtr->x_position) + PHGMATH_Square(posPtr->y_position)) >
			PHGMATH_Square(CylPosLimitCylinder.radius))) {
		
		/* Calculate distance to the surface */
		CylPosCalcDistanceToLimitSurface(posPtr, directionPtr, &distToSurface);
		
		/* If projection to surface puts us past "end" of cylinder, adjust distance to travel */
		if ((posPtr->z_position + (distToSurface * directionPtr->cosine_z)) >
				CylPosLimitCylinder.zMax) {
		
			/* Calculate ratio of delta z to determine new distance to travel */
			adjustmentRatio = (CylPosLimitCylinder.zMax - posPtr->z_position)/
				(distToSurface * directionPtr->cosine_z);
				
			distToSurface = distToSurface * adjustmentRatio;
		}
		else if ((posPtr->z_position + (distToSurface * directionPtr->cosine_z)) <
				CylPosLimitCylinder.zMin) {
		
			/* Calculate ratio of delta z to determine new distance to travel */
			adjustmentRatio = (CylPosLimitCylinder.zMin - posPtr->z_position)/
				(distToSurface * directionPtr->cosine_z);
				
			distToSurface = distToSurface * adjustmentRatio;
		}
	
		/* Save new position */
		posPtr->x_position = posPtr->x_position + (distToSurface * directionPtr->cosine_x);
		posPtr->y_position = posPtr->y_position + (distToSurface * directionPtr->cosine_y);
		posPtr->z_position = posPtr->z_position + (distToSurface * directionPtr->cosine_z);

	}
}

/*********************************************************************************
*
*			Name:			CylPosGetLimitZRange
*
*			Summary:	Return the range limiters for the limit cylinder.
*			Arguments:
*				double	*zMinPtr	- Minimum z value.
*				double	*zMaxPtr	- Maximum z value.
*
*			Function return: None.
*
*********************************************************************************/
void CylPosGetLimitZRange(double *zMinPtr, double *zMaxPtr)	
{
	*zMinPtr = CylPosLimitCylinder.zMin;
	*zMaxPtr = CylPosLimitCylinder.zMax;
}


/*********************************************************************************
*
*			Name:			CylPosInitCriticalZone
*
*			Summary:	Get the user-chosen values for the critical zone.
*
*			Arguments:
*				double	acceptanceAngle	- The acceptance angle.
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean CylPosInitCriticalZone(double acceptanceAngle)	
{
	Boolean 		okay = false;		/* Process success flag */
	double			minAngle, maxAngle; /* Temps for cone beam formula */
	double			radiusOfFocalCircle;
	
	do {	/* Process Loop */
	
		if (	PHG_IsPET() ||
				( PHG_IsSPECT() && (UNCCOLIsConeBeam() == false) )	){
			/* Find the appropriate radius */
			CylPosCriticalZone.radius = CylPosObjectCylinder.radius;

			/* Critical zone is centered over the object cylinder */
			CylPosCriticalZone.centerX = CylPosObjectCylinder.centerX;
			CylPosCriticalZone.centerY = CylPosObjectCylinder.centerY;
			
			/* The axial length is extended as a function of the target size/objectsize/acceptance angle */
			{
				/* Calculate amount to extend each end */
				if(acceptanceAngle != 90.0) {
					CylPosCriticalZone.zMax = CylPosTargetCylinder.zMax + (
						PHGMATH_TanFromDegrees(acceptanceAngle) *
						(CylPosTargetCylinder.radius + CylPosObjectCylinder.radius));
				}
				else {
					CylPosCriticalZone.zMax = CylPosObjectCylinder.zMax;
				}
				
				/* Now clip to object cylinder */
				if (CylPosCriticalZone.zMax > CylPosObjectCylinder.zMax)
					CylPosCriticalZone.zMax = CylPosObjectCylinder.zMax;

				if(acceptanceAngle != 90.0) {
					CylPosCriticalZone.zMin = CylPosTargetCylinder.zMin - (
						PHGMATH_TanFromDegrees(acceptanceAngle) *
						(CylPosTargetCylinder.radius + CylPosObjectCylinder.radius));
				}
				else {
					CylPosCriticalZone.zMin = CylPosObjectCylinder.zMin;
				}
				
				/* Now clip to object cylinder */
				if (CylPosCriticalZone.zMin < CylPosObjectCylinder.zMin)
					CylPosCriticalZone.zMin = CylPosObjectCylinder.zMin;

			}
		}
		else {
			/* Set the appropriate radius */
			CylPosCriticalZone.radius = CylPosObjectCylinder.radius;

			/* Critical zone is centered over the object cylinder */
			CylPosCriticalZone.centerX = CylPosObjectCylinder.centerX;
			CylPosCriticalZone.centerY = CylPosObjectCylinder.centerY;
			
			if (acceptanceAngle < 90.0) {
			
				maxAngle = atan(CylPosTargetCylinder.zMax/UNCCOLGetFocalLength())
					- PHGMATH_RadiansFromDegrees(acceptanceAngle);
					
				minAngle = atan(CylPosTargetCylinder.zMin/UNCCOLGetFocalLength())
					+ PHGMATH_RadiansFromDegrees(acceptanceAngle);
	
				radiusOfFocalCircle = UNCCOLGetFocalLength() - ColRunTimeParams[ColCurParams].UNCSPECTCol.RadiusOfRotation;
				
				CylPosCriticalZone.zMax = CylPosTargetCylinder.zMax -
					(CylPosTargetCylinder.radius-CylPosObjectCylinder.radius) * tan(maxAngle);

				CylPosCriticalZone.zMin =  CylPosTargetCylinder.zMin -
					(CylPosTargetCylinder.radius-CylPosObjectCylinder.radius) * tan(minAngle);
				
				if (CylPosCriticalZone.zMin < CylPosObjectCylinder.zMin) {
					CylPosCriticalZone.zMin = CylPosObjectCylinder.zMin;
				}
				if (CylPosCriticalZone.zMax > CylPosObjectCylinder.zMax) {
					CylPosCriticalZone.zMax = CylPosObjectCylinder.zMax;
				}
					
			}
			else {
	
				CylPosCriticalZone.zMin = CylPosObjectCylinder.zMin;
				CylPosCriticalZone.zMax = CylPosObjectCylinder.zMax;
				
			}
		}	
		#ifdef DUMP_INTERSECTION
			/* Create the file */
			if ((CylPosIntersectFile = LbFlFileOpen("cylpos.intersections", "w")) == 0) {
				
				ErStGeneric("Unable to open intersection file.");
				break;
			}
		#endif
		
		okay = true;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			CylPosInitLimitCylinder
*
*			Summary:	Determine the size of the limit cylinder.
*
*			Arguments:
*
*			Function return: None.
*
*********************************************************************************/
void CylPosInitLimitCylinder(void)	
{
	/* Just calculate the smallest enclosing cylinder */
	CylPosLimitCylinder = CylPosObjectCylinder;
	
	if (CylPosLimitCylinder.radius < CylPosTargetCylinder.radius)
		CylPosLimitCylinder.radius = CylPosTargetCylinder.radius;
		
	if (CylPosLimitCylinder.zMin > CylPosTargetCylinder.zMin)
		CylPosLimitCylinder.zMin = CylPosTargetCylinder.zMin;
		
	if (CylPosLimitCylinder.zMax < CylPosTargetCylinder.zMax)
		CylPosLimitCylinder.zMax = CylPosTargetCylinder.zMax;
		
			
	#ifdef DUMP_CLIPS
		/* Create the file */
		if ((CylPosClipFile = LbFlFileOpen("cylpos.clipping", "w")) == 0) {
			
			PhgAbort("Unable to open clip file (CylPosInitLimitCylinder).", false);
		}
	#endif
}

/*********************************************************************************
*
*			Name:			CylPosInitObjectCylinder
*
*			Summary:	Get the user-chosen values for the object cylinder.
*
*			Arguments:
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean CylPosInitObjectCylinder(void)	
{
	Boolean okay = false;
	
	do {	/* Process Loop */
		#ifdef CYLPOS_DEBUG
			if (!CylPosTargetInitialized == true)
				PhgAbort("You must initialize the target before the object cylinder (CylPosInitObjectCylinder).",
					false);
		#endif
	
		/* The Object cylinder is defined by the sub-object module */
		SubObjGetObjCylinder(&(CylPosObjectCylinder.radius),
			&(CylPosObjectCylinder.zMin), &(CylPosObjectCylinder.zMax),
			&(CylPosObjectCylinder.centerX), &(CylPosObjectCylinder.centerY));
		
		/* Verify object is within target cylinder */
       if (CylPosObjectCylinder.radius > CylPosTargetCylinder.radius) {
            sprintf(cylPosErrString,
                    "You can't have an object bigger than the target!\n"
                    "Object Cylinder Radius = %3.2f\n"
                    "Target Cylinder Radius = %3.2f\n",
                    CylPosObjectCylinder.radius,
                    CylPosTargetCylinder.radius);
            ErStGeneric(cylPosErrString);
            break;
        }
		
		okay = true;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			CylPosProjectToTargetCylinder
*
*			Summary:	Given a location and angle, project to target cylinder.
*			Arguments:
*				PHG_Position		*positionPtr	- The current position.
*				PHG_Direction		*directionPtr	- The current direction.
*				double				*distPtr		- The distance to the cylinder.
*
*			Function return: True if projection is within z limits.
*
*********************************************************************************/
Boolean CylPosProjectToTargetCylinder(PHG_Position *positionPtr,
			PHG_Direction *directionPtr, double *distPtr)	
{
	Boolean			isWithin = false;		/* Is the photon within the z limits? */
	PHG_Position	projectedPosition;		/* Position on cylinder */

	/* Assume we are not within boundary */
	isWithin = false;
	
	do { /* Process Loop */
		
		/* Project it to the cylinder */
		if (CylPosProjectToCylinder(positionPtr,
				directionPtr, &CylPosTargetCylinder,
				&projectedPosition,
				distPtr) == false) {
			break;
		}
				
			
		/* See of outside z limit */
		if ((projectedPosition.z_position > CylPosTargetCylinder.zMax) ||
				(projectedPosition.z_position < CylPosTargetCylinder.zMin)) {
	
			break;
		}
		else {
	
			/* Calculate z intersection */
			positionPtr->z_position = projectedPosition.z_position;

			/* Calculate x intersection */
			positionPtr->x_position = projectedPosition.x_position;
				
			/* Calculate y intersection */
			positionPtr->y_position = projectedPosition.y_position;
		}
		
		/* We made it here, we are within */
		isWithin = true;
	} while (false);
	
	return (isWithin);
}

/*********************************************************************************
*
*			Name:			CylPosWillIntersectCritZone
*
*			Summary:	If the photon's path extended from its position along its 
*						direction will intersect the critical zone, 
*						then indicate that and return the intersection, 
*						otherwise indicate that there is no intersection.
*
*			Arguments:
*				PHG_Position		photon_position				- The current position.
*				PHG_Direction		photon_direction			- The current direction.
*				PHG_Intersection	*critZoneIntersectionPtr	- The intersection information.
*
*			Function return: True if photon will intersect the critical zone.
*
*********************************************************************************/
Boolean CylPosWillIntersectCritZone(PHG_Position photon_position,
			PHG_Direction photon_direction, PHG_Intersection *critZoneIntersectionPtr)	
{
	Boolean 		willIntersect;			/* Intersection flag */
	
	do {	/* Process Loop */
	
		
		/* Assume no intersection */
		willIntersect = false;

		/* If cosine_z == +/- 1 we will not intersect */
		if (PhgMathRealNumAreEqual(fabs(photon_direction.cosine_z), 1.0, -7, 0,
						0, 0)) {
			break;
		}
		
		/* We will not be changing the direciton, so save it */
		critZoneIntersectionPtr->photonsDirection = photon_direction;
		
		/* See if photon will intersect the critical zone on its current path */
		{
			/* See if we are already within the critical zone */
			if ((photon_position.z_position >= CylPosCriticalZone.zMin)
					&& (photon_position.z_position <=
					CylPosCriticalZone.zMax)) {
					
				/* Photon has already intersected */
				willIntersect = true;
				critZoneIntersectionPtr->distToEnter = 0;
				critZoneIntersectionPtr->startingPosition = photon_position;
			}
			else if ((photon_position.z_position < CylPosCriticalZone.zMin)
					&& (photon_direction.cosine_z > 0)) {
					
				/* Photon is below the zone, and traveling up */
				
				/* Calculate distance to enter zone */
				critZoneIntersectionPtr->distToEnter = (CylPosCriticalZone.zMin -
					photon_position.z_position)/photon_direction.cosine_z;
					
				/* Calculate intersection point */
				critZoneIntersectionPtr->startingPosition.x_position = photon_position.x_position
					+ (critZoneIntersectionPtr->distToEnter * photon_direction.cosine_x);
					
				critZoneIntersectionPtr->startingPosition.y_position = photon_position.y_position
					+ (critZoneIntersectionPtr->distToEnter * photon_direction.cosine_y);
					
				critZoneIntersectionPtr->startingPosition.z_position = CylPosCriticalZone.zMin;
				
				/* Verify that after traveling this far, we are within the critical zone */
				if ((PHGMATH_Square((critZoneIntersectionPtr->startingPosition.x_position-
						CylPosCriticalZone.centerX)/
						CylPosCriticalZone.radius) +
						PHGMATH_Square((critZoneIntersectionPtr->startingPosition.y_position-
							CylPosCriticalZone.centerY)/
						CylPosCriticalZone.radius)) <= 1) {
						
					/* We will be intersecting */
					willIntersect = true;
				}
			}
			else if ((photon_position.z_position > CylPosCriticalZone.zMax)
					&& (photon_direction.cosine_z < 0)) {
	
				/* Photon is above, and traveling down */
				
				/* Calculate distance to enter zone */
				critZoneIntersectionPtr->distToEnter = (CylPosCriticalZone.zMax -
					photon_position.z_position)/photon_direction.cosine_z;
					
				/* Calculate intersection point */
				critZoneIntersectionPtr->startingPosition.x_position = 
					photon_position.x_position
					+ (critZoneIntersectionPtr->distToEnter * photon_direction.cosine_x);
					
				critZoneIntersectionPtr->startingPosition.y_position = 
					photon_position.y_position
					+ (critZoneIntersectionPtr->distToEnter * photon_direction.cosine_y);
					
				critZoneIntersectionPtr->startingPosition.z_position = CylPosCriticalZone.zMax;
				
				/* Verify that after traveling this far, we are within the critical zone */
				if ((PHGMATH_Square((critZoneIntersectionPtr->startingPosition.x_position-
						CylPosCriticalZone.centerX)/
						CylPosCriticalZone.radius) +
						PHGMATH_Square((critZoneIntersectionPtr->startingPosition.y_position-
						CylPosCriticalZone.centerY)/
						CylPosCriticalZone.radius)) <= 1) {
						
					/* We will be intersecting */
					willIntersect = true;
				}
				
			}
		}
		
		/* If we won't intersect, bolt */
		if (willIntersect == false)
			break;
			
		/* See if the photon is moving down */
		if (photon_direction.cosine_z < 0) {
			
			/* Calculate distance to exit */
			critZoneIntersectionPtr->distToExit = (CylPosCriticalZone.zMin - 
			   	photon_position.z_position)/photon_direction.cosine_z;
				
			/* Calculate exit point */
			critZoneIntersectionPtr->finalPosition.x_position = 
				photon_position.x_position			
				+ (critZoneIntersectionPtr->distToExit * 
			   photon_direction.cosine_x);
				
			critZoneIntersectionPtr->finalPosition.y_position = 
				photon_position.y_position			
				+ (critZoneIntersectionPtr->distToExit * 
				photon_direction.cosine_y);

			/* See if we go out of the cylinder */
			if ((PHGMATH_Square((critZoneIntersectionPtr->finalPosition.x_position-
					CylPosCriticalZone.centerX)/
					CylPosCriticalZone.radius) +
					PHGMATH_Square((critZoneIntersectionPtr->finalPosition.y_position-
					CylPosCriticalZone.centerY)/
					CylPosCriticalZone.radius)) > 1) {
				
				/* Project to the critical zone */
				(void) CylPosProjectToCylinder(&photon_position,
					&(photon_direction), &CylPosCriticalZone,
					&(critZoneIntersectionPtr->finalPosition),
					&critZoneIntersectionPtr->distToExit);
							
			}
			else {
				
				critZoneIntersectionPtr->finalPosition.z_position = photon_position.z_position			
					+ (critZoneIntersectionPtr->distToExit * photon_direction.cosine_z);

			}
		}
		else if (photon_direction.cosine_z > 0) {
		
			/* Photon is moving upward */

			/* Calculate distance to exit */
			critZoneIntersectionPtr->distToExit = (CylPosCriticalZone.zMax - 
				photon_position.z_position)/photon_direction.cosine_z;
				
			/* Calculate exit point */
			critZoneIntersectionPtr->finalPosition.x_position = photon_position.x_position			
				+ (critZoneIntersectionPtr->distToExit * photon_direction.cosine_x);
				
			critZoneIntersectionPtr->finalPosition.y_position = photon_position.y_position			
				+ (critZoneIntersectionPtr->distToExit * photon_direction.cosine_y);

			/* See if we go out of the cylinder */
			if ((PHGMATH_Square((critZoneIntersectionPtr->finalPosition.x_position-
					CylPosCriticalZone.centerX)/
					CylPosCriticalZone.radius) +
					PHGMATH_Square((critZoneIntersectionPtr->finalPosition.y_position-
					CylPosCriticalZone.centerY)/
					CylPosCriticalZone.radius)) > 1) {
					
				/* Project to the critical zone */
				(void) CylPosProjectToCylinder(
 						&photon_position,
						&(photon_direction), &CylPosCriticalZone,
						&(critZoneIntersectionPtr->finalPosition),
 						&critZoneIntersectionPtr->distToExit);
							
			}
			else {
				
				critZoneIntersectionPtr->finalPosition.z_position = 
					photon_position.z_position			
					+ (critZoneIntersectionPtr->distToExit * 
					photon_direction.cosine_z);

			}
			
		}
		else {
			/* Photon is moving parallel to the x-y plane */
						
			/* Project to the critical zone */
			(void) CylPosProjectToCylinder(
				&photon_position,
				&(photon_direction), &CylPosCriticalZone,
				&(critZoneIntersectionPtr->finalPosition), 
				&critZoneIntersectionPtr->distToExit);
		}
			
		#ifdef PHG_DEBUG
			if ((willIntersect == true) && (
					(critZoneIntersectionPtr->distToEnter < 0) ||
					(critZoneIntersectionPtr->distToExit < 0))) {
				
				PhgAbort("Invalid distance calculations (CylPosWillIntersectCritZone)", true);
			}

			if ((willIntersect == true) && (
					(critZoneIntersectionPtr->distToExit < 
					 critZoneIntersectionPtr->distToEnter))) {
				PhgAbort("Invalid calculation (CylPosWillIntersectCritZone: 2)", true);
			}
		#endif
			
		#ifdef DUMP_INTERSECTION
		{
			char	debugStr[1024];
			double	absDif;
			
			/* Verify position is valid */
			if ((critZoneIntersectionPtr->finalPosition.z_position > CylPosCriticalZone.zMin)
					&& (critZoneIntersectionPtr->finalPosition.z_position < CylPosCriticalZone.zMin)) {
				
				if (PhgMathRealNumAreEqual(
						(PHGMATH_Square(critZoneIntersectionPtr->finalPosition.x_position-
						CylPosCriticalZone.centerX)) +
						PHGMATH_Square(critZoneIntersectionPtr->finalPosition.y_position-
						CylPosCriticalZone.centerY)),
						PHGMATH_Square(CylPosCriticalZone.radius), -7, 0,
						&absDif, 0) == false) {
						
					sprintf(debugStr, "\nInvalid photon/critical zone intersection calculate\n\t x = %3.2e, y = %3.2e, z = %3.2e, r = %3.2e\n x^2 = %3.2e, y^2 = %3.2e, z^2 = %3.2e, sum = %3.2e, r^2 = %3.2e absDif = %3.2e\n",
					critZoneIntersectionPtr->finalPosition.x_position, critZoneIntersectionPtr->finalPosition.y_position, critZoneIntersectionPtr->finalPosition.z_position,
					CylPosCriticalZone.radius,
					PHGMATH_Square(critZoneIntersectionPtr->finalPosition.x_position),
					PHGMATH_Square(critZoneIntersectionPtr->finalPosition.y_position),
					PHGMATH_Square(critZoneIntersectionPtr->finalPosition.z_position),
					PHGMATH_Square(critZoneIntersectionPtr->finalPosition.x_position) + 
					PHGMATH_Square(critZoneIntersectionPtr->finalPosition.y_position),
					PHGMATH_Square(CylPosCriticalZone.radius),
					absDif);
					
					ErBreak(debugStr);
				}
			}
			
			fprintf(CylPosIntersectFile, "%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\t%3.2e\n",
				photon_position.x_position, photon_position.y_position, photon_position.x_position,
				photon_direction.cosine_x, photon_direction.cosine_y, photon_direction.cosine_z,
				critZoneIntersectionPtr->startingPosition.x_position, critZoneIntersectionPtr->startingPosition.y_position,  critZoneIntersectionPtr->startingPosition.z_position,
				critZoneIntersectionPtr->finalPosition.x_position, critZoneIntersectionPtr->finalPosition.y_position,  critZoneIntersectionPtr->finalPosition.z_position,
				critZoneIntersectionPtr->distToEnter, critZoneIntersectionPtr->distToExit);
		}
		#endif
	}	while (false);
	
	
	#ifdef NOT_NOW
	/* Do a little checking for proper inclusion/exclusion */
	{
		#include "EmisList.h"
		
		double	gamma;
		double	temp1, temp2;
		double	cosGamma;
		
		temp1 = ((CylPosObjectCylinder.zMax - CylPosTargetCylinder.zMax)/2) - 
			CylPosCriticalZone.zMax;
			
		temp2 = (CylPosTargetCylinder.zMax + temp1)/CylPosObjectCylinder.radius;
		
		gamma = (PHGMATH_PI_DIV2) - atan(temp2);
		
		cosGamma = cos(gamma);
		
		if (willIntersect) {
			
			if (cosGamma > photon_direction.cosine_z) {
				LbInPrintf("\nInvalid computation from critical zone intersection.");
			}
				
		}
		else {

			
			if (cosGamma < photon_direction.cosine_z) {
				LbInPrintf("\nInvalid computation from critical zone intersection.");
			}
		}
		temp2 = gamma;
	}
	#endif
	return (willIntersect);
}

#ifdef PHG_DEBUG
/*********************************************************************************
*
*			Name:			CylPosDumpObjects
*
*			Summary:		Dump the objects for debuging purposes.
*
*			Arguments:
*
*			Function return: None.
*
*********************************************************************************/
void CylPosDumpObjects()	
{
	if (PhgDebugDumpFile == 0)
		return;
		
	/* Print out target cylinder */
	fprintf(PhgDebugDumpFile, "\nTarget Cylinder\n\tradius\t= % 3.2f\n\tzMin\t= % 3.2f\n\tzMax\t= % 3.2f\n",
		CylPosTargetCylinder.radius, CylPosTargetCylinder.zMin, CylPosTargetCylinder.zMax);

	fprintf(PhgDebugDumpFile, "\tcenterX\t= % 3.2f\tcenterY\t= % 3.2f\n",
		CylPosTargetCylinder.centerX, CylPosTargetCylinder.centerY);

	/* Print out object cylinder */
	fprintf(PhgDebugDumpFile, "\nObject Cylinder\n\tradius\t= % 3.2f\n\tzMin\t\t= % 3.2f\tzMax\t= % 3.2f\n",
		CylPosObjectCylinder.radius, CylPosObjectCylinder.zMin, CylPosObjectCylinder.zMax);

	fprintf(PhgDebugDumpFile, "\tcenterX\t= % 3.2f\tcenterY\t= % 3.2f\n",
		CylPosObjectCylinder.centerX, CylPosObjectCylinder.centerY);

	/* Print out limit cylinder */
	fprintf(PhgDebugDumpFile, "\nLimit Cylinder\n\tradius\t= % 3.2f\n\tzMin\t\t= % 3.2f\tzMax\t= % 3.2f\n",
		CylPosLimitCylinder.radius, CylPosLimitCylinder.zMin, CylPosLimitCylinder.zMax);

	fprintf(PhgDebugDumpFile, "\tcenterX\t= % 3.2f\tcenterY\t= % 3.2f\n",
		CylPosLimitCylinder.centerX, CylPosLimitCylinder.centerY);

	/* Print out object cylinder */
	fprintf(PhgDebugDumpFile, "\nCritical Zone\n\tradius\t= % 3.2f\n\tzMin\t= % 3.2f\tzMax\t= % 3.2f\n",
		CylPosCriticalZone.radius, CylPosCriticalZone.zMin, CylPosCriticalZone.zMax);

	fprintf(PhgDebugDumpFile, "\tcenterX\t= % 3.2f\tcenterY\t= % 3.2f\n",
		CylPosCriticalZone.centerX, CylPosCriticalZone.centerY);
		
	fprintf(PhgDebugDumpFile, "\n\n");
}
#endif /* PHG_DEBUG */
#undef CYLINDER_POSITION
