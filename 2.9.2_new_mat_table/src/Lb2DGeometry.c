/*********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 2002-2012 Department of Radiology	   			*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/*********************************************************************************
*
*		Module Name:		Lb2DGeometry.c
*		Revision Number:	1.17
*		Date last revised:	1 October 2012
*		Programmer:			Steven Gillispie 
*		Date Originated:	9 July 2002
*
*		Module Overview:	Functions for 2D geometry functions.
*
*		References:			None.
*
**********************************************************************************
*
*		Global functions defined:	
*			Lb2DDirCosines
*			Lb2DDirCosComp
*			Lb2DNormalLine
*    	 	Lb2DLineSegsIntersect
*			Lb2DPointParLinesIntersect
*			Lb2DPointRectIntersect
*			Lb2DRectsIntersect
*
*		Global variables defined:	None.
*
**********************************************************************************
*
*		Revision Section (Also update version number, if relevant)
*
*		Programmer(s):		
*
*		Revision date:		
*
*		Revision description:
*
**********************************************************************************
*
*		Revision Section (Also update version number, if relevant)
*
*		Programmer(s):		Steven Gillispie
*
*		Revision date:		1 October 2012
*
*		Revision description:	Changed to use LbMath.h instead of PhgMath.h.
*
**********************************************************************************
*
*		Revision Section (Also update version number, if relevant)
*
*		Programmer(s):		Steven Gillispie
*
*		Revision date:		24 March 2010
*
*		Revision description:	Tightened tolerances from 10^-7 to 10^-10 in 
*							Lb2DLineSegsIntersect and Lb2DPointParLinesIntersect.
*
*********************************************************************************/

#include "Lb2DGeometry.h"


#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMath.h"


/*	CONSTANTS */

/*  LOCAL GLOBALS */

/*	LOCAL MACROS */

/*********************************************************************************
*
*	Name:			LB2D_AreEqual
*
*	Summary:		Compares two values for fuzzy equality:  |v1 - v2| < 10^(-10).
*
*	Arguments:		value1 and value2  (numeric values).
*
*	Function return: C Boolean (zero = false, non-zero = true).
*
*********************************************************************************/

#define LB2D_AreEqual(value1, value2) (fabs((value1) - (value2)) < 1.0E-10)


/*	LOCAL FUNCTIONS */

/*	FUNCTIONS */


/*********************************************************************************
*
*	Name:			Lb2DDirCosines
*
*	Summary:		Compute and return the directon cosines of the line from 
*						point1 to point2.
*					NOTE:  For two identical points, the line points along the x-axis.
*
*	Arguments:		
*		Lb2D_Point*			point1	- The first point.
*		Lb2D_Point*			point2	- The second point.
*		double*				xCos	- X direction cosine.
*		double*				yCos	- Y direction cosine.
*
*	Function return:		None.
*
*********************************************************************************/

void Lb2DDirCosines(	Lb2D_Point	*point1, 
						Lb2D_Point	*point2, 
						double		*xCos, 
						double		*yCos )

{
	double			xOffset;			/* X distance from point1 to point2 */
	double			yOffset;			/* Y distance from point1 to point2 */
	double			segLength;			/* Distance from point1 to point2 */
	
	
	/* Compute the relative distances */
	xOffset = point2->x_position - point1->x_position;
	yOffset = point2->y_position - point1->y_position;
	
	/* Compute the segment length */
	segLength = LBMATH_SquareRoot( 
					LBMATH_Square( xOffset ) + LBMATH_Square( yOffset ) );
	
	/* Compute the direction cosines */
	if ( segLength == 0.0 ) {
		/* Identical points */
		*xCos = 1.0;
		*yCos = 0.0;
	}
	else {
		*xCos = xOffset / segLength;
		*yCos = yOffset / segLength;
	}
}


/*********************************************************************************
*
*	Name:			Lb2DDirCosComp
*
*	Summary:		Compares the lines (xCos1,yCos1) and (xCos2,yCos2) to see which 
*						precedes the other.
*					Order is counterclockwise, when the angle between the lines is 
*						less than pi; for pi separation, the result is arbitrary.
*
*	Arguments:		
*		double				xCos1	- Line 1 x-direction cosine.
*		double				yCos1	- Line 1 y-direction cosine.
*		double				xCos2	- Line 2 x-direction cosine.
*		double				yCos2	- Line 2 y-direction cosine.
*
*	Function return:		int:  -1 for 1<2, 0 for 1=2, 1 for 1>2.
*
*********************************************************************************/

int Lb2DDirCosComp(	double		xCos1, 
					double		yCos1, 
					double		xCos2, 
					double		yCos2 )

{
	int				result;				/* Result of function */
	double			theta1;				/* Angle of line 1 */
	double			theta2;				/* Angle of line 2 */
	
	
	if ( ( xCos1 == xCos2 ) && ( yCos1 == yCos2 ) ) {
		/* Identical lines */
		result = 0;
	}
	else {
		if ( copysign( 1.0, yCos1 ) == copysign( 1.0, yCos2 ) ) {
			/* Both lines on same side of x-axis */
			if ( copysign( 1.0, yCos1 ) > 0.0 ) {
				/* Upper half-plane */
				result = ( xCos1 > xCos2 ) ? -1 : 1;
			}
			else {
				/* Below x-axis */
				result = ( xCos1 > xCos2 ) ? 1 : -1;
			}
		}
		else {
			/* Only one line below x-axis */
			if ( copysign( 1.0, xCos1 ) == copysign( 1.0, xCos2 ) ) {
				/* Both lines on same side of y-axis */
				if ( copysign( 1.0, xCos1 ) > 0.0 ) {
					/* Right half-plane */
					result = ( yCos1 < yCos2 ) ? -1 : 1;
				}
				else {
					/* Left of y-axis */
					result = ( yCos1 < yCos2 ) ? 1 : -1;
				}
			}
			else {
				/* Lines in diagonal quadrants */
				if ( xCos1 == -xCos2 ) {
					/* Lines separated by pi */
					result = 1;
				}
				else {
					/* Determine angles of each line */
					theta1 = acos( xCos1 );
					if ( yCos1 < 0.0 ) {
						theta1 = LBMATH_2PI - theta1;
					}
					theta2 = acos( xCos2 );
					if ( yCos2 < 0.0 ) {
						theta2 = LBMATH_2PI - theta2;
					}
					
					/* Compare the angles */
					if ( theta1 < theta2 ) {
						if ( theta2 < ( theta1 + LBMATH_PI ) ) {
							/* Line 1 precedes line 2 */
							result = -1;
						}
						else {
							/* Line 2 precedes line 1 */
							result = 1;
						}
					}
					else {
						if ( theta1 < ( theta2 + LBMATH_PI ) ) {
							/* Line 2 precedes line 1 */
							result = 1;
						}
						else {
							/* Line 1 precedes line 2 */
							result = -1;
						}
					}
				}
			}
		}
	}
	
	return ( result );
}


/*********************************************************************************
*
*	Name:			Lb2DNormalLine
*
*	Summary:		Computes the coefficients for the (modified) normal form equation 
*						of the line passing through point1 and point2:
*						lineCos * x  +  lineSin * y  +  lineDist = 0.
*					NOTE:  For efficiency, the points are assumed to be *distinct*.
*					NOTE:  "Modified" means lineDist is negative distance.
*
*	Arguments:		
*		Lb2D_Point*			point1	- First point.
*		Lb2D_Point*			point2	- Second point.
*		double*				lineCos	- Returned direction cosine.
*		double*				lineCos	- Returned direction cosine.
*		double*				lineDist- Returned line distance.
*
*	Function return:		None.
*
*********************************************************************************/

void Lb2DNormalLine( 	Lb2D_Point	*point1, 
						Lb2D_Point	*point2, 
						double		*lineCos, 
						double		*lineSin, 
						double		*lineDist )

{
	double			lineA;			/* Parameter A in Ax + By + C = 0 */
	double			lineB;			/* Parameter B in Ax + By + C = 0 */
	double			lineC;			/* Parameter C in Ax + By + C = 0 */
	double			segLength;		/* Length of line segment [point1, point2] */
	
	
	/* Compute the general form equation of the line */
	lineA = point2->y_position - point1->y_position;
	lineB = point1->x_position - point2->x_position;
	lineC = point2->x_position * point1->y_position  -  
			point1->x_position * point2->y_position;
	
	/* Compute the segment length */
	segLength = LBMATH_SquareRoot( LBMATH_Square(lineA) + LBMATH_Square(lineB) );
	
	/* Adjust the length to get the proper sign */
	if ( lineC == 0.0 ) {
		segLength = copysign( segLength, lineB );
	}
	else {
		segLength = - copysign( segLength, lineC );
	}
	
	/* Convert to the modified normal form equation */
	*lineCos = lineA / segLength;
	*lineSin = lineB / segLength;
	*lineDist = lineC / segLength;
}


/*********************************************************************************
*
*	Name:			Lb2DLineSegsIntersect
*
*	Summary:		Compares two line segments for intersection.
*					NOTE:  For efficiency, it is assumed that line segment points are 
*						*distinct*.
*
*	Arguments:		
*		Lb2D_Point*			seg1_Pt1	- First segment, first point.
*		Lb2D_Point*			seg1_Pt2	- First segment, second point.
*		Lb2D_Point*			seg2_Pt1	- Second segment, first point.
*		Lb2D_Point*			seg2_Pt2	- Second segment, second point.
*
*	Function return:		Lb2D_PosDesc:  don't intersect, 
*								intersect at endpoint(s), intersect.
*
*********************************************************************************/

Lb2D_PosDesc Lb2DLineSegsIntersect( Lb2D_Point	*seg1_Pt1, Lb2D_Point	*seg1_Pt2, 
									Lb2D_Point	*seg2_Pt1, Lb2D_Point	*seg2_Pt2 )

{
	Lb2D_PosDesc	result;			/* Result of function */
	double			line1Cos;		/* Cos w for line 1 */
	double			line1Sin;		/* Sin w for line 1 */
	double			line1Dist;		/* Origin distance for line 1 */
	double			line2Cos;		/* Cos w for line 2 */
	double			line2Sin;		/* Sin w for line 2 */
	double			line2Dist;		/* Origin distance for line 2 */
	double			seg1Min;		/* Endpoint of segment 1 closest to origin */
	double			seg1Max;		/* Endpoint of segment 1 farthest from origin */
	double			seg2Min;		/* Endpoint of segment 2 closest to origin */
	double			seg2Max;		/* Endpoint of segment 2 farthest from origin */
	double			p1DistL2;		/* Distance of seg1_Pt1 to line 2 */
	double			p2DistL2;		/* Distance of seg1_Pt2 to line 2 */
	double			p1DistL1;		/* Distance of seg2_Pt1 to line 1 */
	double			p2DistL1;		/* Distance of seg2_Pt2 to line 1 */
	double			signL1s;		/* Carrier of sign of p1DistL1 * p2DistL1 */
	double			signL2s;		/* Carrier of sign of p1DistL2 * p2DistL2 */
	
	
	/* Compute normal forms of lines through each line segment */
	Lb2DNormalLine(seg1_Pt1, seg1_Pt2, 
					&line1Cos, &line1Sin, &line1Dist);
	Lb2DNormalLine(seg2_Pt1, seg2_Pt2, 
					&line2Cos, &line2Sin, &line2Dist);
	
	/* Check for parallel lines */
	if (LB2D_AreEqual(line1Cos, line2Cos) && LB2D_AreEqual(line1Sin, line2Sin)) {
		/* Parallel -- check for collinearity */
		if (LB2D_AreEqual(line1Dist, line2Dist)) {
			/* Collinear */
			if (line1Cos > line1Sin) {
				/* Changing most in X direction */
				if (seg1_Pt1->x_position > seg1_Pt2->x_position) {
					seg1Min = seg1_Pt2->x_position;
					seg1Max = seg1_Pt1->x_position;
				}
				else {
					seg1Min = seg1_Pt1->x_position;
					seg1Max = seg1_Pt2->x_position;
				}
				if (seg2_Pt1->x_position > seg2_Pt2->x_position) {
					seg2Min = seg2_Pt2->x_position;
					seg2Max = seg2_Pt1->x_position;
				}
				else {
					seg2Min = seg2_Pt1->x_position;
					seg2Max = seg2_Pt2->x_position;
				}
			}
			else {
				/* Changing most in Y direction */
				if (seg1_Pt1->y_position > seg1_Pt2->y_position) {
					seg1Min = seg1_Pt2->y_position;
					seg1Max = seg1_Pt1->y_position;
				}
				else {
					seg1Min = seg1_Pt1->y_position;
					seg1Max = seg1_Pt2->y_position;
				}
				if (seg2_Pt1->y_position > seg2_Pt2->y_position) {
					seg2Min = seg2_Pt2->y_position;
					seg2Max = seg2_Pt1->y_position;
				}
				else {
					seg2Min = seg2_Pt1->y_position;
					seg2Max = seg2_Pt2->y_position;
				}
			}
			
			if ((seg2Min > seg1Max) || (seg1Min > seg2Max)) {
				/* Don't overlap */
				result = Lb2D_Outside;
			}
			else if (LB2D_AreEqual(seg2Min, seg1Max) || LB2D_AreEqual(seg1Min, seg2Max)) {
				/* Intersect only at endpoints */
				result = Lb2D_OnBound;
			}
			else {
				/* Intersect in a line segment */
				result = Lb2D_Inside;
			}
		}
		else {
			/* Parallel but not collinear */
			result = Lb2D_Outside;
		}
	}
	else {
		/* Lines are not parallel */
		
		/* Evaluate distances to lines */
		p1DistL2 = line2Cos * seg1_Pt1->x_position  +  
					line2Sin * seg1_Pt1->y_position  +  
					line2Dist;
		if (LB2D_AreEqual(p1DistL2, 0.0)) {
			p1DistL2 = 0.0;
		}
		p2DistL2 = line2Cos * seg1_Pt2->x_position  +  
					line2Sin * seg1_Pt2->y_position  +  
					line2Dist;
		if (LB2D_AreEqual(p2DistL2, 0.0)) {
			p2DistL2 = 0.0;
		}
		p1DistL1 = line1Cos * seg2_Pt1->x_position  +  
					line1Sin * seg2_Pt1->y_position  +  
					line1Dist;
		if (LB2D_AreEqual(p1DistL1, 0.0)) {
			p1DistL1 = 0.0;
		}
		p2DistL1 = line1Cos * seg2_Pt2->x_position  +  
					line1Sin * seg2_Pt2->y_position  +  
					line1Dist;
		if (LB2D_AreEqual(p2DistL1, 0.0)) {
			p2DistL1 = 0.0;
		}
		
		/* Compare signs */
		signL1s = p1DistL1 * p2DistL1;
		signL2s = p1DistL2 * p2DistL2;
		
		if (signL1s > 0) {
			/* Both segment 2 points on same side of line 1 -- can't intersect */
			result = Lb2D_Outside;
		}
		else if (signL2s > 0) {
			/* Both segment 1 points on same side of line 2 -- can't intersect */
			result = Lb2D_Outside;
		}
		else if ((signL1s < 0) && (signL2s < 0)) {
			/* Both endpoints are on opposite sides of the other line segment -- 
				must intersect */
			result = Lb2D_Inside;
		}
		else {
			/* At least one line segment intersects only at an endpoint */
			result = Lb2D_OnBound;
		}
	}
	
	return ( result );
}


/*********************************************************************************
*
*	Name:			Lb2DPointParLinesIntersect
*
*	Summary:		Compares a point and two parallel lines for intersection.
*					NOTE:  For efficiency, the lines are assumed to be parallel.
*
*	Arguments:		
*		Lb2D_Point*			thePoint	- The point.
*		double				line1Cos	- First line, normal form direction cosine.
*		double				line1Sin	- First line, normal form direction sine.
*		double				line1Dist	- First line, normal form line distance.
*		double				line2Cos	- Second line, normal form direction cosine.
*		double				line2Sin	- Second line, normal form direction sine.
*		double				line2Dist	- Second line, normal form line distance.
*
*	Function return:		Lb2D_PosDesc:  point outside lines, 
*								on a line, between the lines.
*
*********************************************************************************/

Lb2D_PosDesc Lb2DPointParLinesIntersect( 	Lb2D_Point	*thePoint, 
											double		line1Cos, 
											double		line1Sin, 
											double		line1Dist, 
											double		line2Cos, 
											double		line2Sin, 
											double		line2Dist )

{
	Lb2D_PosDesc	result;			/* Result of function */
	double			dist1;			/* Distance of thePoint to line 1 */
	double			dist2;			/* Distance of thePoint to line 2 */
	double			prod;			/* Product of dist1 and dist2 */
	
	
	/* Evaluate thePoint in relation to the two lines */
	dist1 = line1Cos * thePoint->x_position  +  
			line1Sin * thePoint->y_position  +  
			line1Dist;
	if (LB2D_AreEqual(dist1, 0.0)) {
		/* Consider < 10^-10 = 0.0 */
		dist1 = 0.0;
	}
	dist2 = line2Cos * thePoint->x_position  +  
			line2Sin * thePoint->y_position  +  
			line2Dist;
	if (LB2D_AreEqual(dist2, 0.0)) {
		/* Consider < 10^-10 = 0.0 */
		dist2 = 0.0;
	}
	
	/* Calculate their distance product (only sign of product is relevant) */
	prod = dist1 * dist2;
	
	/* Determine the line positions relative to the origin */
	if ( fabs( line1Cos ) > fabs( line1Sin ) ) {
		/* Lines more vertical than horizontal */
		if ( copysign( 1.0, line1Cos ) == copysign( 1.0, line2Cos ) ) {
			/* Lines on same side of origin:  no change */
		}
		else {
			/* Lines on opposite side of origin:  reverse sign */
			prod = -prod;
		}
	}
	else {
		/* Lines more horizontal than vertical */
		if ( copysign( 1.0, line1Sin ) == copysign( 1.0, line2Sin ) ) {
			/* Lines on same side of origin:  no change */
		}
		else {
			/* Lines on opposite side of origin:  reverse sign */
			prod = -prod;
		}
	}
	
	/* Set the result */
	if ( prod > 0 ) {
		/* Point is outside lines */
		result = Lb2D_Outside;
	}
	else if ( prod < 0 ) {
		/* Point is inside lines */
		result = Lb2D_Inside;
	}
	else {
		/* Point is on a line */
		result = Lb2D_OnBound;
	}
	
	return ( result );
}


/*********************************************************************************
*
*	Name:			Lb2DPointRectIntersect
*
*	Summary:		Compares a point and a rectangle for intersection.
*
*	Arguments:		
*		Lb2D_Point*			thePoint	- The point.
*		Lb2D_Rect*			theRect		- The rectangle.
*
*	Function return:		Lb2D_PosDesc:  don't intersect, 
*								intersect only on edge, intersect.
*
*********************************************************************************/

Lb2D_PosDesc Lb2DPointRectIntersect( 	Lb2D_Point	*thePoint, 
										Lb2D_Rect	*theRect )

{
	Lb2D_PosDesc	result1;		/* Result of first edge pair comparison */
	Lb2D_PosDesc	result2;		/* Result of second edge pair comparison */
	Lb2D_PosDesc	result;			/* Result of function */
	double			line1Cos;		/* Cos w for line 1 */
	double			line1Sin;		/* Sin w for line 1 */
	double			line1Dist;		/* Origin distance for line 1 */
	double			line2Cos;		/* Cos w for line 2 */
	double			line2Sin;		/* Sin w for line 2 */
	double			line2Dist;		/* Origin distance for line 2 */
	
	
	/* Compute normal forms of lines through one pair of opposite rectangle edges */
	Lb2DNormalLine( &theRect->corner_1, &theRect->corner_2, 
					&line1Cos, &line1Sin, &line1Dist );
	Lb2DNormalLine( &theRect->corner_3, &theRect->corner_4, 
					&line2Cos, &line2Sin, &line2Dist );
	
	/* Evaluate thePoint in relation to the two lines */
	result1 = Lb2DPointParLinesIntersect( thePoint, 
											line1Cos, line1Sin, line1Dist, 
											line2Cos, line2Sin, line2Dist );
	
	if ( result1 == Lb2D_Outside ) {
		/* Point is outside the rectangle */
		result = Lb2D_Outside;
	}
	else {
		/* Compute normal forms of lines through other pair of opposite rectangle edges */
		Lb2DNormalLine( &theRect->corner_1, &theRect->corner_4, 
						&line1Cos, &line1Sin, &line1Dist );
		Lb2DNormalLine( &theRect->corner_2, &theRect->corner_3, 
						&line2Cos, &line2Sin, &line2Dist );
		
		/* Evaluate thePoint in relation to the two lines */
		result2 = Lb2DPointParLinesIntersect( thePoint, 
												line1Cos, line1Sin, line1Dist, 
												line2Cos, line2Sin, line2Dist );
		
		if ( result2 == Lb2D_Outside ) {
			/* Point is outside the rectangle */
			result = Lb2D_Outside;
		}
		else {
			if ( ( result1 == Lb2D_Inside ) && ( result2 == Lb2D_Inside ) ) {
				/* Point is inside the rectangle */
				result = Lb2D_Inside;
			}
			else {
				/* Point is on the edges of the rectangle */
				result = Lb2D_OnBound;
			}
		}
	}
	
	return ( result );
}


/*********************************************************************************
*
*	Name:			Lb2DRectsIntersect
*
*	Summary:		Compares two rectangles for intersection.
*
*	Arguments:		
*		Lb2D_Rect*			rect_1		- First rectangle.
*		Lb2D_Rect*			rect_2		- Second rectangle.
*
*	Function return:		Lb2D_PosDesc:  don't intersect, 
*								intersect only at edges, intersect.
*
*********************************************************************************/

Lb2D_PosDesc Lb2DRectsIntersect( 	Lb2D_Rect	*rect_1, 
									Lb2D_Rect	*rect_2 )

{
	Lb2D_PosDesc	result;			/* Result of function */
	Boolean			onBoundary;		/* Whether boundary intersection occurs */
	
	
	/* Check the rectangle diagonals for intersection */
	result = Lb2DLineSegsIntersect( &rect_1->corner_1, &rect_1->corner_3, 
									&rect_2->corner_1, &rect_2->corner_3 );
	if ( result != Lb2D_Inside ) {
		result = Lb2DLineSegsIntersect( &rect_1->corner_1, &rect_1->corner_3, 
										&rect_2->corner_2, &rect_2->corner_4 );
	}
	if ( result != Lb2D_Inside ) {
		result = Lb2DLineSegsIntersect( &rect_1->corner_2, &rect_1->corner_4, 
										&rect_2->corner_1, &rect_2->corner_3 );
	}
	if ( result != Lb2D_Inside ) {
		result = Lb2DLineSegsIntersect( &rect_1->corner_2, &rect_1->corner_4, 
										&rect_2->corner_2, &rect_2->corner_4 );
	}
	
	if ( result != Lb2D_Inside ) {
		/* Check if any rectangle corners are interior */
		result = Lb2DPointRectIntersect( &rect_1->corner_1, rect_2 );
		onBoundary = ( result == Lb2D_OnBound );
		if ( result != Lb2D_Inside ) {
			result = Lb2DPointRectIntersect( &rect_1->corner_2, rect_2 );
			onBoundary |= ( result == Lb2D_OnBound );
		}
		if ( result != Lb2D_Inside ) {
			result = Lb2DPointRectIntersect( &rect_1->corner_3, rect_2 );
			onBoundary |= ( result == Lb2D_OnBound );
		}
		if ( result != Lb2D_Inside ) {
			result = Lb2DPointRectIntersect( &rect_1->corner_4, rect_2 );
			onBoundary |= ( result == Lb2D_OnBound );
		}
		if ( result != Lb2D_Inside ) {
			result = Lb2DPointRectIntersect( &rect_2->corner_1, rect_1 );
			onBoundary |= ( result == Lb2D_OnBound );
		}
		if ( result != Lb2D_Inside ) {
			result = Lb2DPointRectIntersect( &rect_2->corner_2, rect_1 );
			onBoundary |= ( result == Lb2D_OnBound );
		}
		if ( result != Lb2D_Inside ) {
			result = Lb2DPointRectIntersect( &rect_2->corner_3, rect_1 );
			onBoundary |= ( result == Lb2D_OnBound );
		}
		if ( result != Lb2D_Inside ) {
			result = Lb2DPointRectIntersect( &rect_2->corner_4, rect_1 );
			onBoundary |= ( result == Lb2D_OnBound );
		}
		
		if ( result != Lb2D_Inside ) {
			/* Check for boundary intersection */
			if ( onBoundary ) {
				/* Rectangles never intersected internally, but did on edges */
				result = Lb2D_OnBound;
			}
			else {
				/* Rectangles never intersected in any way */
				result = Lb2D_Outside;
			}
		}
		else {
			/* Rectangles intersect */
			/* (result already set) */
		}
	}
	
	return ( result );
}
