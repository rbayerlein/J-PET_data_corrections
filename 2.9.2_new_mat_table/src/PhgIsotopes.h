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
*     Module Name: 			PhgIsotopes.h
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

#ifndef PHG_ISOTOPES
#define PHG_ISOTOPES


/* CONSTANTS */

/* The following constant must equal the number of enum types in PhgEn_IsotopeTypeTy and 
	the number of strings in phgEn_IsotopeStr below. */
#define NUM_ISOTOPE_TYPES	8


/* TYPES */

/* The following list needs to be in the same order as the isotope file and as the
	phgEn_IsotopeStr variable below. */
typedef enum {
	PhgEn_IsotopType_NULL,
	PhgEn_c11,
	PhgEn_n13,
	PhgEn_o15,
	PhgEn_f18,
	PhgEn_na22,
	PhgEn_ga68,
	PhgEn_rb82
} PhgEn_IsotopeTypeTy;

/* Nuclide */
typedef struct {
	PhgEn_IsotopeTypeTy		isotope;					/* The isotope */
	double					photonEnergy_KEV;			/* The photon's energy */
	double					maximumPositronEnergy;		/* Maximum positron energy */
} PHG_Nuclide;			/* Our nuclide type */


/* VARIABLES */

/* The list of strings in phgEn_IsotopeStr below (defined in PhgIsotopes.c) needs to have 
	the same order as the isotope file and as the PhgEn_IsotopeTypeTy definition above.  
	NUM_ISOTOPE_TYPES (above) should also be changed when adding to the list. */
extern char *phgEn_IsotopeStr[];


#endif /* PHG_ISOTOPES */
