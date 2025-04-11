/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 2002-2006 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		Lb2DGeometry.h
*			Revision Number:	1.8
*			Date last revised:	13 April 2006
*			Programmer:			Steven Gillispie 
*			Date Originated:	9 July 2002
*
*			Module Overview:	Declarations for 2D geometry functions
*
*			References:			none
*
**********************************************************************************
*
*			Global functions defined:	none
*
*			Global variables defined:	none
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

#ifndef LB_2DGEOMETRY
#define LB_2DGEOMETRY


/* Position descriptions for intersecting geometric objects */
typedef enum {
	Lb2D_Outside = -1,
	Lb2D_OnBound = 0,
	Lb2D_Inside = 1
} Lb2D_PosDesc;	

/* 2-dimensional point for 2D geometry calculations */
/*	Note:  Make sure these fields are similar to PHG_Position, then 
		PHG_Position* can be used for Lb2D_Point* */
typedef struct  {	
	double x_position;							/* X position */
	double y_position;							/* Y position */
} Lb2D_Point;

/* 2-dimensional rectangle */
/*	NOTE:  The points are assumed to be described adjacently, but the rotational
			order is irrelevant; thus 1 and 3 are diagonal, as are 2 and 4.  */
typedef struct  {	
	Lb2D_Point corner_1;						/* Corner 1 */
	Lb2D_Point corner_2;						/* Corner 2 */
	Lb2D_Point corner_3;						/* Corner 3 */
	Lb2D_Point corner_4;						/* Corner 4 */
} Lb2D_Rect;


void Lb2DDirCosines(	Lb2D_Point	*point1, 
						Lb2D_Point	*point2, 
						double		*xCos, 
						double		*yCos );
int Lb2DDirCosComp(	double		xCos1, 
					double		yCos1, 
					double		xCos2, 
					double		yCos2 );
void Lb2DNormalLine( 	Lb2D_Point	*point1, 
						Lb2D_Point	*point2, 
						double		*lineCos, 
						double		*lineSin, 
						double		*lineDist );
Lb2D_PosDesc Lb2DLineSegsIntersect( Lb2D_Point	*seg1_Pt1, Lb2D_Point	*seg1_Pt2, 
									Lb2D_Point	*seg2_Pt1, Lb2D_Point	*seg2_Pt2 );
Lb2D_PosDesc Lb2DPointParLinesIntersect( 	Lb2D_Point	*thePoint, 
											double		line1Cos, 
											double		line1Sin, 
											double		line1Dist, 
											double		line2Cos, 
											double		line2Sin, 
											double		line2Dist );
Lb2D_PosDesc Lb2DPointRectIntersect( 	Lb2D_Point	*thePoint, 
										Lb2D_Rect	*theRect );
Lb2D_PosDesc Lb2DRectsIntersect( 	Lb2D_Rect	*rect_1, 
									Lb2D_Rect	*rect_2 );


#endif /* LB_2DGEOMETRY */
