/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                 (C) Copyright 2012 Department of Radiology	    			*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			PhgIsotopes.c
*     Revision Number:		1.0
*     Date last revised:	10 September 2012
*     Programmer:			Steven Gillispie
*     Date Originated:		10 September 2012
*
*     Module Overview:		This file brings all the isotope variables together.
*
*     References:			Items moved from PhgParams.h and Photon.h
*
**********************************************************************************
*
*     Global functions defined:
*
*	  Global variables defined:   
*					phgEn_IsotopeStr
*
*	  Global macros defined:
*
*********************************************************************************
*
*		Revision Section (Also update version number, if relevant)
*
*		Programmer(s):		
*
*		Revision date:		
*
*		Revision description:	 
*
*********************************************************************************/


/* CONSTANTS */


/* TYPES */


/* VARIABLES */

/* The following list needs to have the same order as the isotope file and as the
	PhgEn_IsotopeTypeTy definition in PhgIsotopes.h.  NUM_ISOTOPE_TYPES should also be 
	changed when adding to the list. */
char *phgEn_IsotopeStr[] = {
	"null",
	"c11",
	"n13",
	"o15",
	"f18",
	"na22",
	"ga68",
	"rb82"
};
