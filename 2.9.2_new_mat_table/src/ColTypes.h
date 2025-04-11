/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1994 - 1999 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			ColTypes.h
*     Revision Number:		1.0
*     Date last revised:	13 September, 1999
*     Programmer:			Steven Vannoy
*     Date Originated:		11 July, 1994
*
*     Module Overview:	This is the global types include file for the collimator.
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
*********************************************************************************/

#ifndef COL_TYPES
#define COL_TYPES

#include "LbTypes.h"

#include "Photon.h"

/* CONSTANTS */
#define	 COL_MAX_SCATTERS			9		/* Maximum number of scatters tracked */
/* TYPES */
	
/* The following enumeration flags what type of collimator we have */
typedef enum {
	ColEn_ColType_NULL,				/* 0 */
	ColEn_simple_pet,				/* 1 */
	ColEn_monte_carlo_pet,			/* 2 */
	ColEn_simple_spect,				/* 3 */
	ColEn_unc_spect,				/* 4 */
	ColEn_slat						/* 5 */
} ColEn_CollimatorTypeTy;


/* This data structure is used to create a block of collimated photons.
	The EmisList module will pass an empty block to the collimator module,
	photons that successfully pass through the collimator will be put into
	the block. These photons will then be passed onto the detector module.
*/
typedef struct {
	PHG_TrackingPhoton	CollimatedTrkngBluePhotons[PHG_MAX_DETECTED_PHOTONS];	/* Collimated blue photons */
	PHG_TrackingPhoton	CollimatedTrkngPinkPhotons[PHG_MAX_DETECTED_PHOTONS];	/* Collimated pink photons */
	LbUsFourByte		NumCollimatedBluePhotons;							/* Number of collimated blue photons */
	LbUsFourByte		NumCollimatedPinkPhotons;							/* Number of collimated pink photons */
} CollimatedPhotonsTy;


/* The following data structures are used to define the various data types
	necessary to support the variety of collimators used in SimSET. This includes
	both PET and SPECT. So, there are a lot of types defined here. 
	
	The run-time parameters for the collimator module will have a field for
	each possible collimator type even though only one type will be used
	for any given simulation. This is somewhat undesireable from an software
	engineering viewpoint but is done none the less to simplify the data
	structures.
	
*/

/* The simplest type of collimation we support for PET is simply end-plate
	truncation. The photons are projected from the face of the target cylinder
	a user-specified distance. If they are within the axial limits of the
	Target cylinder then they are accepted, otherwise they are rejected.
*/
typedef struct {
	double	Depth;	/* Depth of the collimator "end-plate" */
} Col_SimplePETTy;

/* For performing Monte-Carlo simulation of PET collimators the user has
	the choice of specifying cylindrical collimators or tapered collimators.
	For tapered, the collimators may actually be "double" tapered. The following
	struct definse the cylindrical and tappered collimator segments. For 
	double tappering the user defines two "rings" or "layers" of collimators,
	the first is defined with inward tappers and the second with outward
	tapers. NOTE: for cylindrical collimators OuterMinZ and OuterMaxZ are
	the same.
*/

/* The following enumeration flags what type of collimator segment we have */
typedef	enum {SegTy_NULL, Cone, Cylinder} ColEn_SegTy;

/* The following struct defines the parameters for Monte Carlo PET collimators.
	These are basically a list of collimator segments. Each segment can be either
	a cone or a cylinder. It will be defined by its type (cone/cylinder), and its
	material (lead, tungsten, air, etc).
*/
typedef struct {

	ColEn_SegTy		SegType;			/* What type of segment is this */
	LbUsFourByte	MaterialIndex;		/* What type of material is the segment made of? */
	double			InnerRadius;		/* Inner radius of collimator */
	double			OuterRadius;		/* Outer radius of collimator */
	double			InnerMinZ;			/* Minimum axial coordinate for inner radius*/
	double			InnerMaxZ;			/* Maximum axial coordinate for inner radius */
	double			OuterMinZ;			/* Minimum axial coordinate for outer radius */
	double			OuterMaxZ;			/* Maximum axial coordinate for outer radius */
	
} Col_Monte_Carlo_PET_SegTy;

/* The following struct defines the parameters for describing a collection
	of the monte carlo pet segments, comprising of an entire monte carlo
	collimator.
*/
typedef struct {
	LbUsFourByte				NumLayers;		/* How many layers do we have */
	LbUsFourByte				NumSegments;	/* How many segments per layer */
	Col_Monte_Carlo_PET_SegTy	*Segments;		/* The segment list */
} Col_Monte_Carlo_PET_Ty;

/* The following enumerated type is for defining hole geometry for the UNC
	style SPECT colimators
*/
typedef enum {HoleTy_NULL, PARALLEL, FAN, CONE} ColEn_HoleTy;

/* The following struct defines the parameters for describing the collimators used
	in the UNC style SPECT Collimator
*/
typedef struct {
	ColEn_HoleTy	HoleGeometry;
	double			FocalLength;
	double			RadiusOfRotation;
	double			Thickness;
	double			HoleRadius;
	double			SeptalThickness;
	double			MinZ;
	double			MaxZ;
	double			StartAngle;
	double			StopAngle;
	Boolean			SumAllViews;
	LbUsFourByte	NumViews;
} Col_UNC_SPECT_Ty;


/* The following struct defines the parameters for describing a slat hole collimator.
	NOTE: it is assumed that a detector module will define the overall geometry of the
	collimator so those parameters are not specified here.  The user must specify a
	width and gap that works out to be an integral of the detector size
*/
typedef struct {
double 				Start;		/* Where the slat starts */
double 				End;		/* Where the slat ends */
LbFourByte			Material;	/* What the slat is made of */
} Col_Slat_Slat_Ty;

typedef struct {
	double 						InnerRadius;		/* Inner radius */
	double						Depth;				/* Depth of collimator */
	LbUsFourByte				NumSlats;			/* Number of slats */
	Col_Slat_Slat_Ty			*Slats;				/* Where and what the slats are */
} Col_Slat_Layer_Ty;

typedef struct {
	LbUsFourByte				NumLayers;			/* Number of layers */
	Col_Slat_Layer_Ty			*Layers;			/* List of layers */
	double						MinAngle;			/* Specified by DETECTOR */
	double						MinZ;				/* Specified by DETECTOR */
	double						MaxZ;				/* Specified by DETECTOR */
} Col_Slat_Ty;


#endif /* COL_TYPES */


