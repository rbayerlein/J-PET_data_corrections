/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1994-2006 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			DetTypes.h
*     Revision Number:		1.0
*     Date last revised:	Aug 2006
*     Programmer:			Steven Vannoy
*     Date Originated:		11 July, 1994
*
*     Module Overview:	This is the global types include file for the detector.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*     Global variables defined:   
*				
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		
*
*			Revision date:		
*
*			Revision description:	
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		August 2006
*
*			Revision description:	Added DetEn_Block to DetEn_DetectorTypeTy.
*
*********************************************************************************/
#ifndef DET_TYPES
#define DET_TYPES


/* CONSTANTS */
/* TYPES */
#define	 DET_MAX_SCATTERS			9		/* Maximum number of scatters tracked */
/* TYPES */

/* The following enumerate indicates the various detector modules supported */
/* Modification to this enumerate must be matched by mods to detEn_DetTypeStr
in DetParams.c */
typedef enum {
	DetEn_DetType_NULL,
	DetEn_simple_pet,
	DetEn_simple_spect,
	DetEn_unc_spect,
	DetEn_Polygonal,
	DetEn_Planar,
	DetEn_Cylindrical,
	DetEn_DualHeaded,
	DetEn_Block
} DetEn_DetectorTypeTy;


/* This data structure is used to create a block of collimated photons.
	The EmisList module will pass an empty block to the collimator module,
	photons that successfully pass through the collimator will be put into
	the block. These photons will then be passed onto the detector module.
*/
typedef struct {
	PHG_TrackingPhoton	DetectedTrkngBluePhotons[PHG_MAX_DETECTED_PHOTONS];	/* Collimated blue photons */
	PHG_TrackingPhoton	DetectedTrkngPinkPhotons[PHG_MAX_DETECTED_PHOTONS];	/* Collimated pink photons */
	LbUsFourByte		NumDetectedBluePhotons;							/* Number of collimated blue photons */
	LbUsFourByte		NumDetectedPinkPhotons;							/* Number of collimated pink photons */
} DetectedPhotonsTy;

#endif /* DET_TYPES */


