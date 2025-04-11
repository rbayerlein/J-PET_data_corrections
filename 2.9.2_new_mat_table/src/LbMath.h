/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                 (C) Copyright 2012 Department of Radiology                    *
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			LbMath.h
*     Revision Number:		1.0
*     Date last revised:	1 October 2012
*     Programmer:			Steven Gillispie
*     Date Originated:		1 October 2012
*
*     Module Overview:	Mathematical definitions needed by libraries and programs.
*
*     References:       See PhgMath.h for original definitions.  
*
**********************************************************************************
*
*     Global functions defined:
*			LBMATH_Square
*			LBMATH_SquareRoot
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
*********************************************************************************/

#ifndef LB_MATH
#define LB_MATH


/* CONSTANTS */
#define LBMATH_PI					3.1415926535897932384626433832795028841971
									/* What we'll use for PI */
#define LBMATH_2PI					6.2831853071795864769252867665590057683942
									/* What we'll use for 2 pi */

/* TYPES */

/* PROTOTYPES */

/* MACROS */

/*********************************************************************************
*
*			Name:			LBMATH_Square
*
*			Summary:		Return the square of a value.
*
*			Arguments:
*				void			value1.
*
*			Function return: Square of value.
*
*********************************************************************************/
#define LBMATH_Square(value1) ((value1) * (value1))


/*********************************************************************************
*
*			Name:			LBMATH_SquareRoot
*
*			Summary:		Return the square root of a value.
*
*			Arguments:
*				void		value1.
*
*			Function return: Square root of value.
*
*********************************************************************************/
#define LBMATH_SquareRoot(value1) sqrt(value1)


#endif /* LB_MATH */
