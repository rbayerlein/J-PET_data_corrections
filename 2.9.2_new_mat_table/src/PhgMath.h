/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992 - 2012 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			PhgMath.h
*     Revision Number:		1.2
*     Date last revised:	1 October 2012
*     Programmer:			Steven Vannoy
*     Date Originated:		Friday September 18, 1992
*
*     Module Overview:	Definitions for math module.
*
*     References:       'PHG Physics/Math Fuctions' PHG design.  
*
**********************************************************************************
*
*     Global functions defined:
*			PHGMATH_Cosine
*			PHGMATH_CosFromDegrees
*			PHGMATH_ElectronRadius
*			PHGMATH_GetRandCosine
*			PHGMATH_Log
*			PHGMATH_Max
*			PHGMATH_RadiansFromDegrees
*			PHGMATH_Square
*			PHGMATH_SquareRoot
*			PHGMATH_Sine
*			PHGMATH_Tangent
*			PHGMATH_TanFromDegrees
*
*     Global variables defined:   None.
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
*			Revision date:		1 October 2012
*
*			Revision description:	Moved some definitions to LbMath.h.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		July 2006
*
*			Revision description:	Changed PI and other PI constants to
*					40 decimal places (enough even if we end up using
*					quad precision), changed the speed of light to its
*					exact value.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		Sept 1, 2006
*
*			Revision description:	Added PhgMathRealNumIsGreater prototype
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		20 May 2005
*
*			Revision description:	Revised random number generation functions.
*
*********************************************************************************/

#ifndef PHG_MATH
#define PHG_MATH

#include "LbMath.h"


/* CONSTANTS */
#define PHGMATH_PI					LBMATH_PI
#define PHGMATH_2PI					LBMATH_2PI
#define PHGMATH_PI_DIV2				1.5707963267948966192313216916397514420986
									/* What we'll use for  pi/2 */
#define PHGMATH_SPEED_OF_LIGHT		2.99792458E10
									/* Centimeters per second */

/* TYPES */

/* PROTOTYPES */
Boolean			PhgMathInit(LbFourByte *randSeed);
void			PhgMathTerminate(void);
void			PhgMathInitRNGFromSeed(LbFourByte randSeed);
void			PhgMathReadSeed(LbFourByte	resetCount);
void			PhgMathWriteSeed(LbFourByte	resetCount);
double 			PhgMathGetRandomNumberOld(void);
double 			PhgMathGetRandomNumber(void);
double 			PhgMathGetDPRandomNumber(void);
void			PhgMathGetTotalFreePaths(double *totalFreePathLength);
Boolean			PhgMathRealNumAreEqual(double r1, double r2, LbOneByte absMag,
					LbOneByte perMag, double *absDifPtr, double *perDifPtr);
Boolean PhgMathRealNumIsGreater(double r1, 			/* First real value */
								double r2, 			/* Second real value */
								LbOneByte Mag,	/* Magnitude of allowable absolute tolerance */
								double *DifPtr	/* Absolute difference (computed here) */
								);
double			PhgMathSampleFromGauss(double mean,
					double standDev);
LbUsFourByte	PhgMathSolveQuadratic(double a, double b, double c,
					double *minRoot, double *maxRoot);

/* MACROS */

/*********************************************************************************
*
*			Name:		PHGMATH_ArcCosine
*
*			Summary:	Calculate an arc cosine.
*			Arguments:
*				double	cosine - Cosine specified in radians.
*				
*			Function return: double, the cosine's angle.
*
*********************************************************************************/
#define PHGMATH_ArcCosine(cosine) ((double) acos(cosine))

/*********************************************************************************
*
*			Name:		PHGMATH_Cosine
*
*			Summary:	Calculate a cosine.
*			Arguments:
*				double	angle - Angle specified in radians.
*				
*			Function return: double, the angle's cosine.
*
*********************************************************************************/
#define PHGMATH_Cosine(angle) ((double) cos(angle))

/*********************************************************************************
*
*			Name:		PHGMATH_CosFromDegrees
*
*			Summary:	Calculate a cosine from degrees.
*			Arguments:
*				double	angle - Angle specified in degrees.
*				
*			Function return: double, the angle's cosine.
*
*********************************************************************************/
#define PHGMATH_CosFromDegrees(angle) ((double) (cos(PHGMATH_RadiansFromDegrees(angle))))

/*********************************************************************************
*
*			Name:		PHGMATH_ElectronRadius
*
*			Summary:	Return the classic electron radius
*			Arguments:
*				
*			Function return: double, the a randomly generated cosine.
*
*********************************************************************************/
#define PHGMATH_ElectronRadius() ((double) (2.817938 * pow(10.0, -13.0)))

/*********************************************************************************
*
*			Name:		PHGMATH_GetRandCosine
*
*			Summary:	Create a random cosine.
*			Arguments:
*				
*			Function return: double, the a randomly generated cosine.
*
*********************************************************************************/
#define PHGMATH_GetRandCosine() ((double) cos(PhgMathGetRandomNumber() * (PHGMATH_2PI)))

/*********************************************************************************
*
*			Name:			PHGMATH_Log
*
*			Summary:		Return the natural log of x.
*
*			Arguments:
*				void			x.
*
*			Function return: Log of x.
*
*********************************************************************************/
#define PHGMATH_Log(x) log(x)

/*********************************************************************************
*
*			Name:			PHGMATH_Max
*
*			Summary:		Return the maximum of two values.
*
*			Arguments:
*				void			value1.
*				void			value2.
*
*			Function return: Maximum of two values.
*
*********************************************************************************/
#define PHGMATH_Max(value1, value2) (((value1) > (value2)) ? (value1) : (value2))

/*********************************************************************************
*
*			Name:			PHGMATH_Min
*
*			Summary:		Return the minimum of two values.
*
*			Arguments:
*				void			value1.
*				void			value2.
*
*			Function return: Minimum of two values.
*
*********************************************************************************/
#define PHGMATH_Min(value1, value2) (((value1) < (value2)) ? (value1) : (value2))

/*********************************************************************************
*
*			Name:		PHGMATH_RadiansFromDegrees
*
*			Summary:	Convert degrees to radians.
*			Arguments:
*				double	angle - Angle specified in degrees.
*				
*			Function return: double, the angle in radians.
*
*********************************************************************************/
#define PHGMATH_RadiansFromDegrees(angle) ((double) (angle) * (PHGMATH_PI/180))

/*********************************************************************************
*
*			Name:		PHGMATH_DegreesFromRadians
*
*			Summary:	Convert radians to degrees.
*			Arguments:
*				double	angle - Angle specified in radians.
*				
*			Function return: double, the angle in degrees.
*
*********************************************************************************/
#define PHGMATH_DegreesFromRadians(angle) ((double) (angle) * (180/PHGMATH_PI))

/*********************************************************************************
*
*			Name:			PHGMATH_Square
*
*			Summary:		Return the square of a value.
*
*			Arguments:
*				void			value1.
*
*			Function return: Square of value.
*
*********************************************************************************/
#define PHGMATH_Square(value1)	LBMATH_Square(value1)

/*********************************************************************************
*
*			Name:			PHGMATH_SquareRoot
*
*			Summary:		Return the square root of a value.
*
*			Arguments:
*				void		value1.
*
*			Function return: Square root of value.
*
*********************************************************************************/
#define PHGMATH_SquareRoot(value1)	LBMATH_SquareRoot(value1)

/*********************************************************************************
*
*			Name:			PHGMATH_RadialPos
*
*			Summary:		Return the square root of the sum of squares.
*
*			Arguments:
*				void		value1.
*				void		value2
*
*			Function return: Square root of value1 squared + value2 squared.
*
*********************************************************************************/
#define PHGMATH_RadialPos(value1, value2) PHGMATH_SquareRoot(PHGMATH_Square(value1)+PHGMATH_Square(value2))

/*********************************************************************************
*
*			Name:		PHGMATH_Sine
*
*			Summary:	Calculate a sine.
*			Arguments:
*				double	angle - Angle specified in radians.
*				
*			Function return: double, the angle's sine.
*
*********************************************************************************/
#define PHGMATH_Sine(angle) ((double) sin(angle))

/*********************************************************************************
*
*			Name:		PHGMATH_SinFromDegrees
*
*			Summary:	Calculate a sine.
*			Arguments:
*				double	angle - Angle specified in degrees.
*				
*			Function return: double, the angle's sine.
*
*********************************************************************************/
#define PHGMATH_SinFromDegrees(angle) ((double) PHGMATH_Sine(PHGMATH_RadiansFromDegrees(angle)))

/*********************************************************************************
*
*			Name:		PHGMATH_Tangent
*
*			Summary:	Calculate a tangent.
*			Arguments:
*				double	angle - Angle specified in radians.
*				
*			Function return: double, the angle's tangent.
*
*********************************************************************************/
#define PHGMATH_Tangent(angle) ((double) tan(angle))

/*********************************************************************************
*
*			Name:		PhgMathTanFromDegrees
*
*			Summary:	Calculate a tangent from degrees.
*			Arguments:
*				double	angle - Angle specified in degrees.
*				
*			Function return: double, the angle's tangent.
*
*********************************************************************************/
#define PHGMATH_TanFromDegrees(angle) ((double) (tan(PHGMATH_RadiansFromDegrees(angle))))

#endif /* PHG_MATH */
