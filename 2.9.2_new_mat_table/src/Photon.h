/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1990-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		Photon.h
*			Revision Number:	1.5
*			Date last revised:	4 June 2013
*			Programmer:			Steven Vannoy, Tom Keenan, Steven Gillispie 
*			Date Originated:	14 March 1990
*
*			Module Overview:	Photon Header
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
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		10 September 2012
*
*			Revision description:	Moved PhgEn_IsotopeTypeTy and PHG_Nuclide to 
*										PhgIsotopes.h.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		January 2008
*
*			Revision description:	 Added element indexing data to detInteractionInfoTy.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*
*			Revision description:  Modification of decay and photon types to:
*				- support eight-byte number of decays
*				- random coincidences
*				- binning by detector crystal
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		March 2005
*
*			Revision description:	 Added Na22 to the supported isotopes for
*								positron range.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		January 2005
*
*			Revision description:	Put a time-stamp on the decays.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:	Put a new field in the detected photon
*						to fix a history file bug.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):			Steven Vannoy		
*		
*			Revision date:			19 August, 1992
*
*			Revision description:	Changed names of types for new PHG
*									implementation.
*
*********************************************************************************/

#ifndef PHOTON
#define PHOTON


/* These are tracking photon modifier flags. Note maximum of 8 because field is a one byte. */
#define	PHGFg_PhotonBlue				LBFlag0
#define PHGFg_TrackAsScatter			LBFlag1
#define PHGFg_TrackAsPrimary			LBFlag2

#define	PHG_IsBlue(trPhoton)			LbFgIsSet((trPhoton)->flags, PHGFg_PhotonBlue)		/* Is this a blue photon? (if not it is pink) */
#define PHG_IsTrackAsScatter(trPhoton)	LbFgIsSet((trPhoton)->flags, PHGFg_TrackAsScatter)	/* Are we tracking as scatter */
#define PHG_IsTrackAsPrimary(trPhoton)	LbFgIsSet((trPhoton)->flags, PHGFg_TrackAsPrimary)	/* Are we tracking as primary */
#define PHGSetTrackAsScatter(trPhoton)	LbFgSet((trPhoton)->flags, PHGFg_TrackAsScatter)	/* Are we tracking as scatter */
#define PHGSetTrackAsPrimary(trPhoton)	LbFgSet((trPhoton)->flags, PHGFg_TrackAsPrimary)	/* Are we tracking as primary */


/* These are flags for the history file decays and detected photon structures */
#define PHGFg_IsADecay					LBFlag0
#define PHGFg_IsAPhoton					LBFlag1
	/* in standard history files, when flagged as photon, flags 2-7 currently used for number of scatters. */

#define PHG_IsADecay(flags)			LbFgIsSet(flags, PHGFg_IsADecay)
#define PHG_IsAPhoton(flags)		LbFgIsSet(flags, PHGFg_IsAPhoton)
#define PHG_SetIsADecay(flags)		LbFgSet(flags, PHGFg_IsADecay)
#define PHG_SetIsAPhoton(flags)		LbFgSet(flags, PHGFg_IsAPhoton)


#define PHG_MAXIMUM_STARTS				20		/* Maximum number of "starts" allowed */
/* Position */
typedef struct  {	
	double x_position;							/* X position */
	double y_position;							/* Y position */
	double z_position;							/* Z position */
}PHG_Position;

/* Direction */
typedef struct {
	double cosine_x;							/* Cosine of the angle along the x direction */
	double cosine_y;							/* Cosine of the angle along the y direction */
	double cosine_z;							/* Cosine of the angle along the z direction */
} PHG_Direction;

/* Detected Position  (reduced to float) */
typedef struct  {	
	float x_position;							/* X position */
	float y_position;							/* Y position */
	float z_position;							/* Z position */
} PHG_DetPosition;

/* Detected Direction (reduced to float) */
typedef struct {
	float cosine_x;							/* Cosine of the angle along the x direction */
	float cosine_y;							/* Cosine of the angle along the y direction */
	float cosine_z;							/* Cosine of the angle along the z direction */
} PHG_DetDirection;

/* Generated decay type */
typedef enum {
	PhgEn_SinglePhoton,
	PhgEn_Positron,
	PhgEn_PETRandom,		/* Artificial random coincidence events */
	PhgEn_Complex,		/* For isotopes with multiple possible decay products - not implemented */
	PhgEn_Unknown		/* for unassigned or error situations */
} PhgEn_DecayTypeTy;	


/* Decay */
typedef struct {
	PHG_Position		location;		/* Origination of decay */
	double				startWeight;	/* starting weight for this decay */
	/* double			eventWeight;	old decay weight for possible changes during tracking */
	double				decayTime;		/* time, in seconds, between scan start and the decay */
	PhgEn_DecayTypeTy	decayType;		/* single photon/positron/multi-emission/random etc. */
} PHG_Decay;			

/* Old Decay - for reading old history files */
typedef struct {
	PHG_Position	location;			/* Origination of decay */
	double			startWeight;		/* Decay weight */
	double			eventWeight;		/* Decay weight */
} PHG_OldDecay;			

/* Decay -- appears unused, delete?
typedef struct {*/
	/*PHG_DetPosition	location;*/			/* Origination of decay */
	/*double			weight;*/				/* Decay weight */
/*} PHG_DetDecay;
*/


/* Interaction Information */
typedef struct {
	LbFourByte	angleIndex;			/* Stratification angle at time of interaction */
	LbFourByte	sliceIndex;			/* Slice it originated from */
} PHG_InteractionInfo;

/* Intersection */
typedef struct {
	double			distToEnter;				/* Distance to enter zone */
	double			distToExit;					/* Distance to exit zone */
	PHG_Position	startingPosition;			/* Starting position for intersection */
	PHG_Position	finalPosition;				/* Final position for intersection */
	PHG_Direction	photonsDirection;			/* Direction at time of intersection */
} PHG_Intersection;

/*	
	The following two structs define a type for storing the position of
	an interaction along with the energy deposited there.
	
*/
typedef struct {	/* Ring, block, layer, and element of a position */
	LbUsFourByte	ringNum;
	LbUsFourByte	blockNum;
	LbUsFourByte	layerNum;
	LbUsFourByte	elementNum;
} DetElementPosIndexTy;

typedef struct {
	PHG_Position			pos;
	double					energy_deposited;
	Boolean					isActive;	/* Is the layer the interaction took place in active */
	DetElementPosIndexTy	posIndices;	/* Usually used only by the block detector */
} detInteractionInfoTy;

#define MAX_DET_INTERACTIONS	30

/* Tracking Photon */
typedef struct  {	
	LbUsOneByte					flags;									/* Modifier flags for the photon - see PHGFg_... above */
	PHG_Position				location;								/* Photon location */
	PHG_Direction				angle;									/* Photon direction */
	LbFourByte					sliceIndex;								/* Slice it is in */
	LbFourByte					angleIndex;								/* Stratification angle it is in */
	LbFourByte					origSliceIndex;							/* Slice it originated from */
	LbFourByte					origAngleIndex;							/* Stratification angle it originated from */
	LbFourByte					xIndex;									/* X index into object */
	LbFourByte					yIndex;									/* Y index into object */
	LbUsFourByte				num_of_scatters;						/* Number of photon scatters */
	LbUsFourByte				scatters_in_col;						/* Number of photon scatters in the collimator */
	double						photon_scatter_weight;					/* Photon scatter weight in PHG */
	double						photon_primary_weight;					/* Photon primary weight in PHG */
	double						photon_current_weight;					/* Photon weight in modules after the PHG */
	double						scatter_target_weight;					/* Used to for weight windowing */
	double						decay_weight;							/* Decay weight */
	double						energy;									/* Photon energy */
	double						travel_distance;						/* Photon travel distance */
	LbUsEightByte				number;									/* Just to help debug */
	double						transaxialPosition;						/* For SPECT, transaxial position on "back" of collimator/detector */
	LbFourByte					azimuthalAngleIndex;					/* For SPECT/DHCI, index of collimator/detector angle */
	double						axialPosition;							/* For SPECT, axial position on "back" of collimator */
	double						detectorAngle;							/* For SPECT, DHCI, angle of detector */
	LbUsFourByte				numStarts;								/* Number of starts used */
	PHG_InteractionInfo			starts_list[PHG_MAXIMUM_STARTS];		/* Photon starts list */
	detInteractionInfoTy		det_interactions[MAX_DET_INTERACTIONS];	/* Interactions in the detector */
	LbUsFourByte				num_det_interactions;					/* The number of detector interactions */
	PHG_Position				detLocation;							/* Centroid location in detector coordinates */
	LbFourByte					detCrystal;								/* for block detectors, the crystal number for detection */
}PHG_TrackingPhoton;

/* time */
typedef struct {
	int hour;
	int min;
	int sec;
} PHG_Time;

/* date */
typedef struct  {
	int month;
	int day;
	int year;
}PHG_Date;

/* detected_photon used for history files */
typedef struct {
	PHG_DetPosition		location;				/* photon current location or, in detector list mode, detection position */
	PHG_DetDirection	angle;					/* photon durrent direction.  perhaps undefined in detector list mode. */
	LbUsOneByte			flags;					/* Photon flags */
	double				photon_weight;			/* Photon's weight */
	float				energy;					/* Photon's energy */
	double				time_since_creation;	/* In seconds */
	float				transaxialPosition;		/* For SPECT, transaxial position on "back" of collimator/detector */
	LbTwoByte			azimuthalAngleIndex;	/* For SPECT/DHCI, index of collimator/detector angle */
	float				detectorAngle;			/* For SPECT, DHCI, angle of detector */
	LbFourByte			detCrystal;				/* for block detectors, the crystal number for detection */
} PHG_DetectedPhoton;

/* detected_photon used for reading version 2.6.2.5 and 2.6.2.6 history file photons */
typedef struct {
	PHG_DetPosition		location;				/* photon current location or, in detector list mode, detection position */
	PHG_DetDirection	angle;					/* photon durrent direction.  perhaps undefined in detector list mode. */
	LbUsOneByte			flags;					/* Photon flags */
	double				photon_weight;			/* Photon's weight */
	float				energy;					/* Photon's energy */
	double				time_since_creation;	/* In seconds */
	float				transaxialPosition;		/* For SPECT, transaxial position on "back" of collimator/detector */
	LbTwoByte			azimuthalAngleIndex;	/* For SPECT/DHCI, index of collimator/detector angle */
	float				detectorAngle;			/* For SPECT, DHCI, angle of detector */
} PHG_OldDetectedPhoton2;

/* old detected_photon, used for reading in old history file photons */
typedef struct {
	PHG_DetPosition		location;				/* Location? */
	PHG_DetDirection	angle;					/* Angle? */
	LbUsOneByte			flags;					/* Photon flags */
	double				photon_weight;			/* Photon's weight */
	float				energy;					/* Photon's energy */
	double				time_since_creation;	/* In seconds */
	float				transaxialPosition;		/* For SPECT, transaxial position on "back" of collimator/detector */
	LbTwoByte			azimuthalAngleIndex;	/* For SPECT/DHCI, index of collimator/detector angle */
} PHG_OldDetectedPhoton1;

/* header - Photon history header template */
typedef struct  {
	PHG_Date	update_date;				/* Current date of header  */
	PHG_Time	update_time;				/* Current time of header  */
	char 		*string_header;				/* Current text of header  */
}PHG_Header;

#endif /* PHOTON */
