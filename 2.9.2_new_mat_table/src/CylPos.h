/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1992-1997 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		CylPos.h
*			Revision Number:	1.0
*			Date last revised:	Wednesday, September 28, 1993
*			Programmer:			Steven Vannoy
*			Date Originated:	Wednesday, September 9, 1992
*
*			Module Overview:	Definitions for CylPos.c.
*
*			References:			'Cylinder Position Processes' PHG design.
*
**********************************************************************************
*
*			Global functions defined:	
*
*			Global variables defined:		none
*
*			Global Macros
*				CYLPOSGetObjZMin
*				CYLPOSGetObjZMax
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
*********************************************************************************/
#ifndef CYL_POS_HDR
#define CYL_POS_HDR

#ifdef CYLINDER_POSITION
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */
/* TYPES */
typedef struct {
	double			radius;		/* Radius of cylinder */
	double			zMin;		/* Minimum z coordinate of cylinder */
	double			zMax;		/* Maximum z coordinate of cylinder */
	double			centerX;	/* Center point on x axis */
	double			centerY;	/* Center point on y axis */
} CylPosCylinderTy;

/* MACROS */

/* GLOBALS */
LOCALE	CylPosCylinderTy	CylPosTargetCylinder;		/* The target cylinder */
LOCALE	CylPosCylinderTy	CylPosCriticalZone;			/* The critical zone */
LOCALE	CylPosCylinderTy	CylPosObjectCylinder;		/* The object cylinder */
LOCALE	CylPosCylinderTy	CylPosLimitCylinder;		/* The limit cylinder */

/* MACROS */

/*********************************************************************************
*
*			Name:		CYLPOSGetObjZMin
*
*			Summary:	Return min Z value for object cylinder.
*			Arguments:
*				
*			Function return: double, the objects z min.
*
*********************************************************************************/
#define CYLPOSGetObjZMin() (CylPosObjectCylinder.zMin)

/*********************************************************************************
*
*			Name:		CYLPOSGetObjZMax
*
*			Summary:	Return max Z value for object cylinder.
*			Arguments:
*				
*			Function return: double, the objects z max.
*
*********************************************************************************/
#define CYLPOSGetObjZMax() (CylPosObjectCylinder.zMax)

/*********************************************************************************
*
*			Name:		CYLPOSGetCritZMin
*
*			Summary:	Return min Z value for critical zone.
*			Arguments:
*				
*			Function return: double, the critical zone's z min.
*
*********************************************************************************/
#define CYLPOSGetCritZMin() (CylPosCriticalZone.zMin)

/*********************************************************************************
*
*			Name:		CYLPOSGetCritZMax
*
*			Summary:	Return max Z value for critical zone.
*			Arguments:
*				
*			Function return: double, the critical zone's z max.
*
*********************************************************************************/
#define CYLPOSGetCritZMax() (CylPosCriticalZone.zMax)

/* PROTOTYPES */
Boolean CylPosFind2dIntersection(double x0, double y0, double cosX, double cosY,
			double radius, double *d1, double *d2);
void 	CylPosCalcDistanceToLimitSurface(PHG_Position *posPtr, PHG_Direction *dirPtr,
			double *distPtr);	
void 	CylPosCalcDistanceToObjectSurface(PHG_Position *posPtr, PHG_Direction *dirPtr,
			double *distPtr);	
void	CylPosCalcDistanceToTargetSurface(PHG_Position *posPtr, PHG_Direction *dirPtr,
		double *distPtr);
void 	CylPosClipToLimitCylinder(PHG_Position *posPtr, PHG_Direction *directionPtr);
void	CylPosGetLimitZRange(double *zMinPtr, double *zMaxPtr);
void	CylPosGetTargetZRange(double *zMinPtr, double *zMaxPtr);
double	CylPosGetTargetRadius(void);
Boolean CylPosInitCriticalZone(double acceptanceAngle);
void 	CylPosInitLimitCylinder(void);
double	CylPosGetObjectRadius(void);
double	CylPosGetTargetRadius(void);
Boolean CylPosInitObjectCylinder(void);
Boolean CylPosInitTargetCylinder(LbPfHkTy paramFlHk, LbUsFourByte numParams);
Boolean CylPosIsOutsideObjCylinder(PHG_Position	photon_position);
Boolean	CylPosProjectToCylinder(PHG_Position *positionPtr,
			PHG_Direction *directionPtr, CylPosCylinderTy *cylinderPtr,
			PHG_Position *newPosPtr, double *distPtr);
	
Boolean CylPosProjectToTargetCylinder(PHG_Position *positionPtr,
			PHG_Direction *directionPtr, double *distPtr);
Boolean CylPosWillIntersectCritZone(PHG_Position photon_position,
			PHG_Direction photon_direction, PHG_Intersection *critZone_Intersection);	
void CylPosCalcDistanceToCylSurface(PHG_Position *posPtr, PHG_Direction *dirPtr,
		CylPosCylinderTy *cylinderPtr, double *distPtr);
		

#ifdef PHG_DEBUG
void	CylPosDumpObjects(void);
#endif
#undef LOCALE
#endif /* CYL_POS_HDR */
