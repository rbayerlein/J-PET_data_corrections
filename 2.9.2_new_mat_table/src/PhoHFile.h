/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1992-2005 Department of Radiology       		     *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		PhoHFile.h
*			Revision Number:	1.3
*			Date last revised:	September 2005
*			Programmer:			Steven Vannoy
*			Date Originated:	24 August 1992
*
*			Module Overview:	Definitions for PhoHFile.c.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:	
*
*			Global variables defined:		none
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
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*
*			Revision description:	
*				- New custom history file options:  decay position, decay time,
*				detector crystal, detector angle.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		September 2005
*
*			Revision description:	Added isRandomsAdded field to header.
*									Added prototype for new function
*									PhoHFileRewriteDetections.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		24 August 2005
*
*			Revision description:	Added isTimeSorted field to header.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		Jan-Feb 2005
*
*			Revision description:	
*					- Added PhoHFileEn_PHG2625, PhoHFileEn_COL2625,
*						PhoHFileEn_DET2625
*					- added custom header fields for decay location and time
*					- changed some counting fields to LbUsEightByte to allow
*						longer simulations.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		5 January 2004
*
*			Revision description:	Removed unused PhoHFileStHeader;
*									Added PhoHFileEventType and PhoHFileReadEvent.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		15 July 2004
*
*			Revision description:	
*					- modified PhoHFileHdrKindTy to include old types,
*					- added PhoHFileOldReadEvent
*
*********************************************************************************/

#ifndef HIST_FILE_HDR
#define HIST_FILE_HDR

#ifdef PHOTON_HIST_FILE
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/*	CONSTANTS	*/
#define	PHOHFILE_OLD_HEADER_VERSION	1.00
#define	PHOHFILE_CUR_HEADER_VERSION	1.00		/* This is the current header version.
												It is also the first version of our
												headers that contains a version number.
												As changes in the header occur so
												too will the header version change.
												All attempts will be made to
												provide backward compatibility.
												*/
										
/* GLOBAL TYPES */

/* There is one header for all history/histogram files, however the information
	that has actually been initialized may be different depending on
	whether or not collimation and detection were performed.  The
	following enumeration indicates what type of information has been
	initialized in the header.
	
	When history/histogram files are combined together the runtime parameters
	will be set to the initial file that contributed to the combined one.  There
	are several statistics that maintain a running total and reflect the combined
	contribution of all files.  These parameters are declared below the field
	"NumSimulations".
*/
/* This enum gives the possible types of history files.
	They should remain in this order so that old files can be read by
	new software.
	The "OLD" types (e.g. PhoHFileEn_PHGOLD) were renamed when the detectedPhoton type
	was changed to fix a bug (with the addition of detectorAngle in version 2.6.2.5).
	
	The "2625" types (e.g. PhoHFileEn_PHG2625) were renamed when, in version 2.9,
	(1) the PHG_Decay type was changed to include decay time and decay type, and
	(2) the PHG_DetectedPhoton type was changed to include detCrystal. (2625 stands
	for versions 2.6.2.5 - 2.6.2.6)
	
	Version 2.9 and above use the PhoHFileEn_PHG/COL/DET enums and include all the
	above additions to the decays and history file detected photons.)
*/
typedef enum { PhoHFileEn_PHGOLD,
				PhoHFileEn_COLOLD,
				PhoHFileEn_DETOLD,
				PhoHFileEn_BIN_WT,
				PhoHFileEn_BIN_WTSQ,
				PhoHFileEn_BIN_CT,
				PhoHFileEn_BIN,
				PhoHFileEn_PHG2625,
				PhoHFileEn_COL2625,
				PhoHFileEn_DET2625,
				PhoHFileEn_PHG,
				PhoHFileEn_COL,
				PhoHFileEn_DET
} PhoHFileHdrKindTy;

/* This is the header that used to get written to the history file */
typedef union {

	struct {
		LbUsFourByte			HdrSize;				/* The size of the header */
		PhoHFileHdrKindTy		HdrKind;				/* What type of header (see above) */
		float					HdrVersion;				/* What version is this header */
		PhgRunTimeParamsTy		PhgRunTimeParams;		/* The run time parameters for the PHG */
		ColRunTimeParamsTy		ColRunTimeParams;		/* The run time parameters for the collimator */
		DetRunTimeParamsTy		DetRunTimeParams;		/* The run time parameters for the collimator */
		PHG_BinParamsTy			BinRunTimeParams;		/* The run time parameters for the binning module */
		CylPosCylinderTy		TargetCylinder;			/* The target cylinder */
		CylPosCylinderTy		CriticalZone;			/* The critical zone */
		CylPosCylinderTy		ObjectCylinder;			/* The object cylinder */
		CylPosCylinderTy		LimitCylinder;			/* The limit cylinder */
		LbUsFourByte			NumSimulations;			/* The number of simulations contributing to 
															the file. In order to support combining
															multiple history and histogram files we
															need to keep track of the number of decays
															requested as well as the number of photons
															detected. This flag indicates that this
															process must be done to analyse the file
														*/
		LbUsEightByte		SumEventsToSimulate;		/*	Total number of decays requested, this 
															value is summed when multiple history or
															histogram files are combined into one.
														*/
		LbUsEightByte		NumPhotons;					/* The number of photons in the file */
		LbUsEightByte		NumDecays;					/* The number of decays in the file */
		double				weightSum;					/* Sum of weights in image */
		double				weightSquSum;				/* Sum of weights squared in image */
		Boolean				isTimeSorted;				/* True if history file has been timesorted */
		Boolean				isRandomsAdded;				/* True if history file has been timesorted */
	} H;
	char Buffer[8192];
} PhoHFileHdrTy;

/*	This type is returned by the PhoHFileReadEvent function to indicate what type of 
		event was read */
typedef enum { PhoHFileNullEvent, PhoHFileDecayEvent, PhoHFilePhotonEvent } PhoHFileEventType;


/* PROGRAM TYPES */

/*	This list of enumerates specifies the possible custom history file fields.
	This allows one monster data structure to house a custom history file definition.
	When adding to this list, always add a label in the same position in PhoHFileRunTimeParamLabels
	(in PhoHFile.c) and a field in PhoHFileRunTimeParamsTy (below).
	Never change the order of items already in the list as this will destroy backward compatibility
	Always leave PhoHFileEn_NULL as the last enum in this list: it is the enum that ends
	loops.
*/
typedef enum {
		/* Our label table */
		PhoHFileEn_x_position,
		PhoHFileEn_x_pos_min,
		PhoHFileEn_x_pos_max,
		PhoHFileEn_y_position,
		PhoHFileEn_y_pos_min,
		PhoHFileEn_y_pos_max,
		PhoHFileEn_z_position,
		PhoHFileEn_z_pos_min,
		PhoHFileEn_z_pos_max,
		PhoHFileEn_x_cosine,
		PhoHFileEn_x_cos_min,
		PhoHFileEn_x_cos_max,
		PhoHFileEn_y_cosine,
		PhoHFileEn_y_cos_min,
		PhoHFileEn_y_cos_max,
		PhoHFileEn_z_cosine,
		PhoHFileEn_z_cos_min,
		PhoHFileEn_z_cos_max,
		PhoHFileEn_scatters_in_object,
		PhoHFileEn_obj_scatters_min,
		PhoHFileEn_obj_scatters_max,
		PhoHFileEn_scatters_in_collimator,
		PhoHFileEn_col_scatters_min,
		PhoHFileEn_col_scatters_max,
		PhoHFileEn_decay_weight,
		PhoHFileEn_weight,
		PhoHFileEn_energy,
		PhoHFileEn_energy_min,
		PhoHFileEn_energy_max,
		PhoHFileEn_travel_distance,
		PhoHFileEn_travel_distance_min,
		PhoHFileEn_travel_distance_max,
		PhoHFileEn_transaxial_distance,
		PhoHFileEn_transaxial_distance_min,
		PhoHFileEn_transaxial_distance_max,
		PhoHFileEn_azimuthal_angle_index,
		PhoHFileEn_aa_index_min,
		PhoHFileEn_aa_index_max,
		PhoHFileEn_axial_position,
		PhoHFileEn_axial_position_min,
		PhoHFileEn_axial_position_max,
		PhoHFileEn_detector_x_position,
		PhoHFileEn_detector_x_position_min,
		PhoHFileEn_detector_x_position_max,
		PhoHFileEn_detector_y_position,
		PhoHFileEn_detector_y_position_min,
		PhoHFileEn_detector_y_position_max,
		PhoHFileEn_detector_z_position,
		PhoHFileEn_detector_z_position_min,
		PhoHFileEn_detector_z_position_max,
		PhoHFileEn_num_detector_interactions,
		PhoHFileEn_num_detector_interactions_min,
		PhoHFileEn_num_detector_interactions_max,
		PhoHFileEn_det_interaction_positions,
		PhoHFileEn_detector_angle,
		PhoHFileEn_detector_angle_min,
		PhoHFileEn_detector_angle_max,
		PhoHFileEn_decay_x_position,
		PhoHFileEn_decay_x_pos_min,
		PhoHFileEn_decay_x_pos_max,
		PhoHFileEn_decay_y_position,
		PhoHFileEn_decay_y_pos_min,
		PhoHFileEn_decay_y_pos_max,
		PhoHFileEn_decay_z_position,
		PhoHFileEn_decay_z_pos_min,
		PhoHFileEn_decay_z_pos_max,
		PhoHFileEn_decay_time,
		PhoHFileEn_decay_type,
		PhoHFileEn_detector_crystal,
		PhoHFileEn_NULL
}PhoHFileEn_RunTimeParamsTy;

/* see comment above PhoHFileEn_RunTimeParamsTy when altering this struct */
typedef struct {
	Boolean		doXposition;
	double		xPosMin;
	double		xPosMax;
	Boolean		doXrange;
	Boolean		doYposition;
	double		yPosMin;
	double		yPosMax;
	Boolean		doYrange;
	Boolean		doZposition;
	double		zPosMin;
	double		zPosMax;
	Boolean		doZrange;
	Boolean		doXcosine;
	double		xCosMin;
	double		xCosMax;
	Boolean		doXCosRange;
	Boolean		doYcosine;
	double		yCosMin;
	double		yCosMax;
	Boolean		doYCosRange;
	Boolean		doZcosine;
	double		zCosMin;
	double		zCosMax;
	Boolean		doZCosRange;
	Boolean		doScattersInObject;
	LbFourByte	objScattersMin;
	LbFourByte	objScattersMax;
	Boolean		doObjScatRange;
	Boolean		doScattersInCollimator;
	LbFourByte	colScattersMin;
	LbFourByte	colScattersMax;
	Boolean		doColScatRange;
	Boolean		doDecayWeight;
	Boolean		doWeight;
	Boolean		doEnergy;
	double		energyMin;
	double		energyMax;
	Boolean		doEnergyRange;
	Boolean		doTravelDistance;
	double		travelDistanceMin;
	double		travelDistanceMax;
	Boolean		doTravDistRange;
	Boolean		doTransaxialDistance;
	double		transaxialDistanceMin;
	double		transaxialDistanceMax;
	Boolean		doTransDistRange;
	Boolean		doAzimuthalAngleIndex;
	LbFourByte	aaIndexMin;
	LbFourByte	aaIndexMax;
	Boolean		doAziAngRange;
	Boolean		doAxialPosition;
	double		axialPositionMin;
	double		axialPositionMax;
	Boolean		doAxiPosRange;
	Boolean		doDetectorXposition;
	double		detectorXpositionMin;
	double		detectorXpositionMax;
	Boolean		doDetectorXposRange;
	Boolean		doDetectorYposition;
	double		detectorYpositionMin;
	double		detectorYpositionMax;
	Boolean		doDetectorYposRange;
	Boolean		doDetectorZposition;
	double		detectorZpositionMin;
	double		detectorZpositionMax;
	Boolean		doDetectorZposRange;
	Boolean		doNumDetectorInteractions;
	LbFourByte	numDetInteractionsMin;
	LbFourByte	numDetInteractionsMax;
	Boolean		doNumDetInteractionsRange;
	Boolean		doDetInteractionPos;
	Boolean		doDetectorAngle;
	Boolean		doDetectorAngleRange;
	double		detectorAngleMin;
	double		detectorAngleMax;
	Boolean		doDecayXposition;
	double		decayXPosMin;
	double		decayXPosMax;
	Boolean		doDecayXrange;
	Boolean		doDecayYposition;
	double		decayYPosMin;
	double		decayYPosMax;
	Boolean		doDecayYrange;
	Boolean		doDecayZposition;
	double		decayZPosMin;
	double		decayZPosMax;
	Boolean		doDecayZrange;
	Boolean		doDecayTime;
	Boolean		doDecayType;
	Boolean		doDetectorCrystal;
} PhoHFileRunTimeParamsTy;

/* This is the data structure that gets used by the history file to know what
	header, what file, etc. to access on calls.
*/
typedef struct {
	PhoHFileHdrTy				header;					/* The header */
	Boolean						doCustom;				/* Do custom writing? */
	PhoHFileRunTimeParamsTy		customParams;			/* For custom fields */	
	char						customParamsName[256];	/* Name of custom parameter file */
	LbUsEightByte				bluesReceived;			/* Blues received */
	LbUsEightByte				bluesAccepted;			/* Blues accepted */
	LbUsEightByte				pinksReceived;			/* Pinks received */
	LbUsEightByte				pinksAccepted;			/* Pinks accepted */
	FILE						*histFile;				/* The history file */
	LbHdrHkTy					headerHk;				/* Hook to the file header */
} PhoHFileHkTy;


/* PROTOTYPES */
Boolean	PhoHFileClose(PhoHFileHkTy *hdrHkTyPtr);
Boolean PhoHFileCreate(char *histFilePath, char *histParamsFilePath, PhoHFileHdrKindTy hdrType,
			PhoHFileHkTy *hdrHkTyPtr);	
Boolean PhoHFileWriteDetections(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay,
			PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
			PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons);
Boolean	PhoHFileRewriteDetections(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay,
			PHG_DetectedPhoton *bluePhotons, LbUsFourByte numBluePhotons,
			PHG_DetectedPhoton *pinkPhotons, LbUsFourByte numPinkPhotons);
Boolean PhoHFileGetRunTimeParams(char *paramPath, PhoHFileRunTimeParamsTy *customParams);	
PhoHFileEn_RunTimeParamsTy  PhoHFileLookupRunTimeParamLabel(char *label);	
void	PhoHFilePrintParams(PhoHFileHkTy *hdrHkTyPtr);
void	PhoHFilePrintReport(PhoHFileHkTy *hdrHkTyPtr);
Boolean PhoHFileReadFields(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay,
			PHG_TrackingPhoton *photon, Boolean *photonAccepted);
void	PhoHFilePrintFields(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay,
			PHG_TrackingPhoton *photon);
PhoHFileEventType PhoHFileReadEvent(FILE *historyFile, 
			PHG_Decay *decayPtr, PHG_DetectedPhoton *photonPtr);
PhoHFileEventType PhoHFileOldReadEvent(FILE *historyFile, 
			PHG_OldDecay *oldDecayPtr, PHG_DetectedPhoton *photonPtr, Boolean isOldPhotons1, Boolean isOldPhotons2 );

#undef LOCALE
#endif /* HIST_FILE_HDR */

