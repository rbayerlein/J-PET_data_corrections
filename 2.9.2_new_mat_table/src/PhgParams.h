/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1993-2012 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			PhgParams.h
*     Revision Number:		1.1
*     Date last revised:	10 September 2012
*     Programmer:			Steven Vannoy
*     Date Originated:		28, July, 1993
*
*     Module Overview:	This is the global include file for the PhgParams.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*	  Global variables defined:   
*				PhgRunTimeParams
*
*	  Global macros defined:
*				PHGGetSineOfAccAngle
*				PHGGetLengthOfScan
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		10 September 2012
*
*			Revision description:	Moved phgEn_IsotopeStr and NUM_ISOTOPE_TYPES to 
*										PhgIsotopes.h
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*			Revision description:
*						- support for randoms and eight-byte number of decays
*						- binning by random-state and crystal number.
*						- binning PET data as SPECT
*						- moved binning options here from PhgBin.h
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:	 Added support for time-of-flight.
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
*********************************************************************************/

#ifndef PHG_PARAMS_HDR_HDR
#define PHG_PARAMS_HDR_HDR

#ifdef PHG_PARAMS
	#define LOCALE
#else	
	#define LOCALE extern
#endif


#include "PhgIsotopes.h"



/* PROGRAM CONSTANTS */

/* PROCESSING CONSTANTS */

#define PATH_LENGTH		256					/* Length of file name paths */
#ifdef __MWERKS__
#define PHG_MAX_DETECTED_PHOTONS	50		/* Maximum number of detected photons for a single decay */
#else
#define PHG_MAX_DETECTED_PHOTONS	50		/* Maximum number of detected photons for a single decay */
#endif
#define	PHG_MIN_PHOTON_ENERGY		1.0		/* This is the absolute minimum energy supported */
#ifdef __MWERKS__
#define PHG_MAX_PARAM_FILES	10				/* Maximum number of parameter files available */
#else
#define PHG_MAX_PARAM_FILES	100				/* Maximum number of parameter files available */
#endif

/* OPTION MACROS */

	/* Calculate the number of decays from activity and scan time? */
#define Phg_IsCalcEventsToSimulate()	PhgRunTimeParams.PhgIsCalcEventsToSimulate

	/* Stratification */
#define	PHG_IsStratification()			PhgRunTimeParams.PhgIsStratification

	/* Forced detection on */
#define PHG_IsForcedDetection()			PhgRunTimeParams.PhgIsForcedDetection

	/* No forcing of non absorption */
#define PHG_IsNoForcedNonAbsorbtion()	!PhgRunTimeParams.PhgIsForcedNonAbsorbtion

	/* No generation of primary events */
#define PHG_IsNoPrimaryEvents()			false

	/* No scatter */
#define PHG_IsNoScatter()				false

	/* Did user request SPECT simulation */
#define PHG_IsSPECT()					PhgRunTimeParams.PhgIsSPECT

	/* Did user request PET coincidence-only (no randoms) simulation */
#define PHG_IsPETCoincidencesOnly()		PhgRunTimeParams.PhgIsPETCoincidencesOnly

	/* Did user request PET coincidence plus singles simulation */
#define PHG_IsPETCoincPlusSingles()		PhgRunTimeParams.PhgIsPETCoincPlusSingles

	/* Did user request multi-emission isotope simulation */
#define PHG_IsMultiEmission()			PhgRunTimeParams.PhgIsMultiEmission

	/* Did user request one of the PET options? */
#define PHG_IsPET()						PhgRunTimeParams.PhgIsPET

	/* Did user request one of the PET coincidence options */
#define PHG_IsPETCoincidences()			PhgRunTimeParams.PhgIsPETCoincidences

	/* Did user request one of the PET singles options */
#define PHG_IsPETSingles()				PhgRunTimeParams.PhgIsPETSingles

	/* Did user request no history file */
#define PHG_IsNoHist()					!PhgRunTimeParams.PhgIsHistoryFile

	/* Did user request history file */
#define PHG_IsHist()					PhgRunTimeParams.PhgIsHistoryFile

	/* Did user request custom history file parameters */
#define PHG_IsHistParams()				PhgRunTimeParams.PhgIsHistoryParamsFile

	/* Did user request adjustment for positron range */
#define PHG_IsRangeAdjust()				PhgRunTimeParams.PhgIsAdjForPosRange

	/* Did user request adjustment for non-collinearity */
#define PHG_IsNonCollinearityAdjust()	PhgRunTimeParams.PhgIsAdjForCollinearity

	/* Did user specify an input productivity table */
#define PHG_IsProductivityComputed()	PhgRunTimeParams.PhgIsComputedProductivityTbl

	/* Our we simulating voxels as point sources */
#define PHG_IsVoxelPointSource()		PhgRunTimeParams.PhgIsVoxelPointSource

	/* Our we simulating voxels as line sources */
#define PHG_IsVoxelLineSource()			PhgRunTimeParams.PhgIsVoxelLineSource

	/* Do we bin on the fly */
#define PHG_IsBinOnTheFly()				PhgRunTimeParams.PhgIsBinOnTheFly

	/* Do we collimate on the fly */
#define PHG_IsCollimateOnTheFly()		PhgRunTimeParams.PhgIsCollimateOnTheFly

	/* Do we detect on the fly */
#define PHG_IsDetectOnTheFly()			PhgRunTimeParams.PhgIsDetectOnTheFly

	/* Are we modeling coherent scatter in both object and tomo */
#define PHG_IsModelCoherent()			PhgRunTimeParams.PhgIsModelCoherent

	/* Are we modeling coherent scatter in the tomograph */
#define PHG_IsModelCoherentInTomo()		PhgRunTimeParams.PhgIsModelCoherentInTomo

	/* Are we modeling coherent scatter in the object*/
#define PHG_IsModelCoherentInObj()		PhgRunTimeParams.PhgIsModelCoherentInObj

	/* Are we modeling polarization */
#define PHG_IsModelPolarization() 		PhgRunTimeParams.PhgIsModelPolarization


/* PROGRAM TYPES */

/* When changing the following list also change phgRunTimeParamLabels in PhgParams.c.
The two lists must have the same order */ 
typedef enum {
	PhgEn_simulate_stratification,
	PhgEn_simulate_forced_detection,
	PhgEn_forced_non_absorbtion,	/* see also PhgEn_forced_non_absorption below */
	PhgEn_acceptance_angle,
	PhgEn_num_to_simulate,
	PhgEn_simulate_SPECT,
	PhgEn_simulate_PETCoincidencesOnly,
	PhgEn_simulate_PETCoincPlusSingles,
	PhgEn_simulate_MultiEmission,
	PhgEn_adjust_for_positron_range,
	PhgEn_adjust_for_collinearity,
	PhgEn_model_coherent,
	PhgEn_model_coherent_in_tomo,
	PhgEn_model_coherent_in_obj,
	PhgEn_minimum_energy,
	PhgEn_photon_energy,
	PhgEn_weight_window_ratio,
	PhgEn_point_source_voxels,
	PhgEn_line_source_voxels,
	PhgEn_random_seed,
	PhgEn_isotope,
	PhgEn_isotope_data_file,
	PhgEn_length_of_scan,
	PhgEn_model_polarization,
	SubObjEn_object,
	SubObjEn_num_slices,
	SubObjEn_slice,
	SubObjEn_slice_number,
	SubObjEn_zMin,
	SubObjEn_zMax,
	SubObjEn_xMin,
	SubObjEn_xMax,
	SubObjEn_yMin,
	SubObjEn_yMax,
	SubObjEn_num_X_bins,
	SubObjEn_num_Y_bins,
	SubObjEn_num_act_X_bins,
	SubObjEn_num_act_Y_bins,
	SubObjEn_num_att_X_bins,
	SubObjEn_num_att_Y_bins,
	CylPosEn_target_cylinder,
	CylPosEn_zMin,
	CylPosEn_zMax,
	CylPosEn_centerX,
	CylPosEn_centerY,
	CylPosEn_radius,
	SubObjEn_activity_indexes,
	SubObjEn_activity_table,
	SubObjEn_activity_index_trans,
	SubObjEn_activity_image,
	SubObjEn_attenuation_indexes,
	SubObjEn_attenuation_table,
	SubObjEn_attenuation_index_trans,
	SubObjEn_attenuation_image,
	SubObjEn_coherent_scatter_table,
	ProdTblEn_productivity_input_table,
	ProdTblEn_productivity_output_table,
	PhoTrkEn_forced_detection_table,
	PhoHStatEn_statistics_file,
	PhgEn_bin_params_file,
	PhgEn_collimator_params_file,
	PhgEn_detector_params_file,
	PhgEn_tomograph_params_file,
	PhoHFileEn_history_file,
	PhoHFileEn_history_params_file,
	PhgEn_forced_non_absorption,	/* correct spelling!, replacing PhgEn_forced_non_absorbtion */
	PhgEn_NULL					/* NULL must always be left last when adding to list,
								it is used to end loops */
}PhgEn_RunTimeParamsTy;


/* Define binning enumerates */
/* When changing the following list also change phgBinParamLabels in PhgParams.c.
The two lists must have the same order */
typedef enum {
	PhgBinEn_num_z_bins,
	PhgBinEn_num_pa_bins,
	PhgBinEn_num_td_bins,
	PhgBinEn_num_aa_bins,
	PhgBinEn_num_e_bins,
	PhgBinEn_num_e1_bins,
	PhgBinEn_num_e2_bins,
	PhgBinEn_scatter_param,
	PhgBinEn_scatter_random_param,
	PhgBinEn_accept_randoms,
	PhgBinEn_num_theta_bins,
	PhgBinEn_num_phi_bins,
	PhgBinEn_num_xr_bins,
	PhgBinEn_num_yr_bins,
	PhgBinEn_num_tof_bins,
	PhgBinEn_min_z,
	PhgBinEn_max_z,
	PhgBinEn_min_td,
	PhgBinEn_max_td,
	PhgBinEn_min_pa,
	PhgBinEn_max_pa,
	PhgBinEn_min_e,
	PhgBinEn_max_e,
	PhgBinEn_min_s,
	PhgBinEn_max_s,
	PhgBinEn_max_theta,
	PhgBinEn_min_xr,
	PhgBinEn_max_xr,
	PhgBinEn_min_yr,
	PhgBinEn_max_yr,
	PhgBinEn_min_tof,
	PhgBinEn_max_tof,
	PhgBinEn_bin_by_crystal,
	PhgBinEn_sum_according_to_type,
	PhgBinEn_weight_image_type,
	PhgBinEn_count_image_type,
	PhgBinEn_add_to_existing_img,
	PhgBinEn_weight_image_path,
	PhgBinEn_weight_squared_image_path,
	PhgBinEn_count_image_path,
	PhgBinEn_do_ssrb,
	PhgBinEn_do_msrb,
	PhgBinEn_detector_radius,
	PhgBinEn_history_file,
	PhgBinEn_history_params_file,
	PhgBinEn_binPETasSPECT,
	PhgBinEn_NULL				/* NULL must always be left last when adding to list,
								it is used to end loops */
}PhgEn_BinParamsTy;

/* Binning parameter types */
/* The following data structure is used for ordering the bin parameters. */
/* Update PHGBIN_NUM_DIMENSIONS if this list is updated */
typedef enum {	PhgBinEn_Null,
				PhgBinEn_Energy1,
				PhgBinEn_Energy2,
				PhgBinEn_Scatter1,
				PhgBinEn_Scatter2,
				PhgBinEn_Z1,
				PhgBinEn_Z2,
				PhgBinEn_TD,
				PhgBinEn_AA,
				PhgBinEn_TOF,
				PhgBinEn_THETA,
				PhgBinEn_PHI,
				PhgBinEn_XR,
				PhgBinEn_YR,
				PhgBinEn_Crystal1,
				PhgBinEn_Crystal2
} PHG_BinEnDimensionsTy;
#define PHGBIN_NUM_DIMENSIONS 16

typedef struct {
	LbFourByte		numDimensions;		/* Number of dimensions user has specified */
	LbUsFourByte	numZBins;			/* Number of z bins */
	LbUsFourByte	numPABins;			/* Number of polar axis bins */
	LbUsFourByte	numTDBins;			/* Number of transaxial distance bins */
	LbUsFourByte	numAABins;			/* Number of azimuthal angle bins */
	LbUsFourByte	numTOFBins;			/* Number of time-of-flight bins */
	LbUsFourByte	numE1Bins;			/* Number of photon 1 energy bins */
	LbUsFourByte	numE2Bins;			/* Number of photon 2 energy bins */
	LbUsFourByte	numS1Bins;			/* Number of photon 1 scatter bins */
	LbUsFourByte	numS2Bins;			/* Number of photon 2 scatter bins */
	LbUsFourByte	numPHIBins;			/* Number of azimuthal angles */
	LbUsFourByte	numThetaBins;			/* Number of elevation angles */
	LbUsFourByte	numXRBins;				/* Number of Xr bins */
	LbUsFourByte	numYRBins;			/* Number of Yr bins */
	LbUsFourByte	numCrystalBins;		/* Number of crystal bins (in one dimension - square for PET) */
	LbUsFourByte	numImageBins;		/* Number of bins in the image */
	LbUsFourByte	scatterRandomParam;	/* Scatter/Random binning parameter indicates how
											scatter and random coincidences will be binned
										*/
										
	Boolean			acceptRandoms;		/* specifies whether random events will be
										accepted or rejected */
	
	double			minZ;				/* Minimum Z axis value */
	double			maxZ;				/* Maximum Z axis value */
	double			minPA;				/* Minimum Polar angle value */
	double			maxPA;				/* Maximum Polar angle value */
	double			minTD;				/* Minimum Transaxial distance value */
	double			maxTD;				/* Maximum Transaxial distance value */
	double			minAA;				/* Minimum Azimuthal angle value */
	double			maxAA;				/* Maximum Azimuthal angle value */
	double			minTOF;				/* Minimum time-of-flight value, nanoseconds */
	double			maxTOF;				/* Maximum time-of-flight value, nanoseconds */
	double			minE;				/* Minimum Energy value */
	double			maxE;				/* Maximum Energy value */
	LbUsFourByte	minS;				/* Minimum scatters */
	LbUsFourByte	maxS;				/* Maximum scatters */
	double			minTheta;				/* Minimum Elevation angle value */
	double			maxTheta;				/* Maximum Elevation angle value */
	double			minPhi;				/* Minimum phi angle value */
	double			maxPhi;				/* Maximum phi angle value */
	double			minXR;					/* Minimum Xr */
	double			maxXR;					/* Maximum Xr */
	double			minYR;					/* Minimum Yr */
	double			maxYR;					/* Maximum Yr */
	
	Boolean			addToExistingImg;	/* Add the results of this run to pre-existing image */
	Boolean			doCounts;			/* Create image of counts */
	Boolean			doWeights;			/* Create image of weights */
	Boolean			doWeightsSquared;	/* Create image of weights squared */
	Boolean			sumAccordingToType;	/* Sum data according to output type? */
	
	/* These fields are computed from fields above, and are here for convenience */
	double			zRange;				/* Range of z values */
	double			eRange;				/* Range of e values */
	double			sRange;				/* Range of s values */
	double			tdRange;			/* Range of transaxial distance values */
	double			tofRange;			/* Range of time-of-flight values */
	double			thetaRange;			/* Range of theta values in 3DRP binning */
	double			phiRange;			/* Range of phi values in 3DRP binning */
	double			xrRange;			/* Range of X r values */
	double			yrRange;			/* Range of Y r values */
	
	LbUsFourByte	scatter2CIsize;		/* Size of bin element for count image */
	LbUsFourByte	scatter2WIsize;		/* Size of bin element for weight image */
	LbUsFourByte	scatter2WISsize;	/* Size of bin element for weight squared image */
	
	LbUsFourByte	scatter1CIsize;		/* Size of bin element for count image */
	LbUsFourByte	scatter1WIsize;		/* Size of bin element for weight image */
	LbUsFourByte	scatter1WISsize;	/* Size of bin element for weight squared image */
	
	LbUsFourByte	energy2CIsize;		/* Size of bin element for count image */
	LbUsFourByte	energy2WIsize;		/* Size of bin element for weight image */
	LbUsFourByte	energy2WISsize;		/* Size of bin element for weight squared image */
	
	LbUsFourByte	energy1CIsize;		/* Size of bin element for count image */
	LbUsFourByte	energy1WIsize;		/* Size of bin element for weight image */
	LbUsFourByte	energy1WISsize;		/* Size of bin element for weight squared image */
	
	LbUsFourByte	aaCIsize;			/* Size of bin element for count image */
	LbUsFourByte	aaWIsize;			/* Size of bin element for weight image */
	LbUsFourByte	aaWISsize;			/* Size of bin element for weight squared image */
	
	LbUsFourByte	tdCIsize;			/* Size of bin element for count image */
	LbUsFourByte	tdWIsize;			/* Size of bin element for weight image */
	LbUsFourByte	tdWISsize;			/* Size of bin element for weight squared image */
	
	LbUsFourByte	tofCIsize;			/* Size of bin element for count image */
	LbUsFourByte	tofWIsize;			/* Size of bin element for weight image */
	LbUsFourByte	tofWISsize;			/* Size of bin element for weight squared image */
	
	LbUsFourByte	paCIsize;			/* Size of bin element for count image */
	LbUsFourByte	paWIsize;			/* Size of bin element for weight image */
	LbUsFourByte	paWISsize;			/* Size of bin element for weight squared image */
	
	LbUsFourByte	z2CIsize;			/* Size of bin element for count image */
	LbUsFourByte	z2WIsize;			/* Size of bin element for weight image */
	LbUsFourByte	z2WISsize;			/* Size of bin element for weight squared image */
	
	LbUsFourByte	z1CIsize;			/* Size of bin element for count image */
	LbUsFourByte	z1WIsize;			/* Size of bin element for weight image */
	LbUsFourByte	z1WISsize;			/* Size of bin element for weight squared image */
	
	LbUsFourByte	phiCIsize;			/* Size of bin element for count image */
	LbUsFourByte	phiWIsize;			/* Size of bin element for weight image */
	LbUsFourByte	phiWISsize;			/* Size of bin element for weight squared image */
	
	LbUsFourByte	thetaCIsize;		/* Size of bin element for count image */
	LbUsFourByte	thetaWIsize;		/* Size of bin element for weight image */
	LbUsFourByte	thetaWISsize;		/* Size of bin element for weight squared image */
	
	LbUsFourByte	xrCIsize;			/* Size of bin element for count image */
	LbUsFourByte	xrWIsize;			/* Size of bin element for weight image */
	LbUsFourByte	xrWISsize;			/* Size of bin element for weight squared image */
	
	LbUsFourByte	yrCIsize;			/* Size of bin element for count image */
	LbUsFourByte	yrWIsize;			/* Size of bin element for weight image */
	LbUsFourByte	yrWISsize;			/* Size of bin element for weight squared image */
	
	LbUsFourByte	crystal2CIsize;		/* Size of bin element for count image */
	LbUsFourByte	crystal2WIsize;		/* Size of bin element for weight image */
	LbUsFourByte	crystal2WISsize;	/* Size of bin element for weight squared image */
	
	LbUsFourByte	crystal1CIsize;		/* Size of bin element for count image */
	LbUsFourByte	crystal1WIsize;		/* Size of bin element for weight image */
	LbUsFourByte	crystal1WISsize;	/* Size of bin element for weight squared image */
	
	LbUsEightByte	countImageSize;		/* Size of the count image buffer */
	LbUsEightByte	weightImageSize;	/* Size of the weight image buffer */
	LbUsEightByte	weightSquImageSize;	/* Size of the weight squared image */
	
	LbUsFourByte	weight_image_type;	/* Type of data for weight image */
	LbUsFourByte	count_image_type;	/* Type of data for count image */
	Boolean			doCrystal;			/* Should do binning by crystal */
	Boolean			doSSRB;				/* Should do ssrb binning */
	Boolean			doMSRB;				/* Should do msrb binning */
	double			image_radius;		/* Image radius for reconstruction image */
	double			detector_radius;	/* Detector radius for reconstruction image */
	
	PHG_BinEnDimensionsTy	PhgBinDimensions[PHGBIN_NUM_DIMENSIONS]; /* Ordering of the dimensions */
	char			weightImgFilePath[PATH_LENGTH];				/* Weight image file  path */
	char			weightSquImgFilePath[PATH_LENGTH];			/* Weight squared image file  path */
	char			countImgFilePath[PATH_LENGTH];				/* Count image file  path */
	char			history_path[PATH_LENGTH];					/* Name of history file */
	char			history_params_path[PATH_LENGTH];			/* Name of history parameters */
	Boolean			isHistoryFile;								/* Flag for doing a history file */
	Boolean			isHistoryParamsFile;						/* Flag for doing a custom history file */
	
	Boolean			isBinPETasSPECT;		/* flag for binning a PET simulation for single photon */

} PHG_BinParamsTy;

/* Runtime params */
typedef struct {

/* PROGRAM GLOBALS */
LbUsEightByte	Phg_EventsToSimulate;			/* Number of events to simulate */
LbUsFourByte	Phg_EventsToSimulateOld;		/* Old 4-byte number of events for reading in old headers */
Boolean			PhgIsCalcEventsToSimulate;		/* Are we calculating the number of decays to simulate from time and activity? */
LbUsFourByte	PhgRandomSeed;					/* Our random seed */
float			Phg_LengthOfScan;				/* The length of our scan */
double			Phg_AcceptanceAngle;			/* The user requested acceptance angle */
double			Phg_SineOfAcceptanceAngle;		/* The sine of the acceptance angle */
double			PhgMinimumEnergy;				/* Minimum energy value for photons */
double			PhgMinWWRatio;					/* Minimum ratio for weight windowing */
double			PhgMaxWWRatio;					/* Maximum ratio for weight windowing  */
PHG_Nuclide		PhgNuclide;						/* The simulation nuclide */

Boolean			PhgIsForcedDetection;			/* Are we doing forced detection? */
Boolean			PhgIsStratification;			/* Are we doing stratification? */
Boolean			PhgIsForcedNonAbsorbtion;		/* Are we preventing absorbtion? */
Boolean			PhgIsSPECT;						/* Are we simulating SPECT? */
Boolean			PhgIsPETCoincidencesOnly;		/* Are we simulating PET coincidences only? */
Boolean			PhgIsPETCoincPlusSingles;		/* Are we simulating PET coincidences and singles? */
Boolean			PhgIsMultiEmission;				/* Are we simulating a multi-emission isotope? */
Boolean			PhgIsPET;						/* Are we simulating PET? (any PET option) */
Boolean			PhgIsPETCoincidences;			/* Are we simulating PET coincidences? (any PET coincidence option) */
Boolean			PhgIsPETSingles;				/* Are we simulating PET singles? (any PET singles option) */
Boolean			PhgIsHistoryFile;				/* Are we creating a history file? */
Boolean			PhgIsHistoryParamsFile;			/* Do we have a history parameters file? */
Boolean			PhgIsAdjForPosRange;			/* Are we adjusting for positron range? */
Boolean			PhgIsAdjForCollinearity;		/* Are we adjusting for collinearity? */
Boolean			PhgIsComputedProductivityTbl;	/* Are we using a computed productivity table? */
Boolean			PhgIsVoxelPointSource;			/* Should we force all decay locations to the center of the voxel? */
Boolean			PhgIsVoxelLineSource;			/* Should we force all decay locations to a line in the center of the voxel? */
Boolean			PhgIsBinOnTheFly;				/* Do we bin on the fly? */
Boolean			PhgIsCollimateOnTheFly;			/* Do we collimate on the fly? */
Boolean			PhgIsDetectOnTheFly;			/* Do we detect on the fly? */
Boolean			PhgIsModelCoherent;				/* Do we model coherent scatter in tomo and obj? */
Boolean			PhgIsModelCoherentInTomo;		/* Do we model coherent scatter in tomo only? */
Boolean			PhgIsModelCoherentInObj;		/* Do we model coherent scatter in obj only? */
Boolean			PhgIsModelPolarization;			/* Do we model polarization? */

char			PhgParamFilePath[PATH_LENGTH];						/* Our param file path */

char			PhgSubObjActivityIndexFilePath[PATH_LENGTH];		/* Activity index file */
char			PhgSubObjActivityTableFilePath[PATH_LENGTH];		/* Activity table file */
char			PhgSubObjActIndexTransFilePath[PATH_LENGTH];		/* Activity index to table translation file */
char			PhgSubObjActImgFilePath[PATH_LENGTH];				/* Attenuation image file  path */
char			PhgSubObjAttenIndexFilePath[PATH_LENGTH];			/* Attenuation index file */
char			PhgSubObjAttenTableFilePath[PATH_LENGTH];			/* Attenuation table file */
char			PhgSubObjAttIndexTransFilePath[PATH_LENGTH];		/* Attenuation index to table translation file */
char			PhgSubObjAttImgFilePath[PATH_LENGTH];				/* Actvity image file  path */
char			PhgSubObjCohScatTableFilePath[PATH_LENGTH];			/* Coherent scatter angles file path */
char			PhgProdTblInputTableFilePath[PATH_LENGTH];			/* Productivity table file */
char			PhgProdTblOutputTableFilePath[PATH_LENGTH];			/* Productivity table file */
char			PhgPhoTrkForcedDetectionFilePath[PATH_LENGTH];		/* Forced Detection table file */
char			PhgPhoHStatFilePath[PATH_LENGTH];					/* Statistics file  path */
char			PhgPhoHFileHistoryFilePath[PATH_LENGTH];			/* History file  path */
char			PhgPhoHParamsFilePath[PATH_LENGTH];					/* History file parameters path */
char			PhgBinParamsFilePath[PHG_MAX_PARAM_FILES][PATH_LENGTH];	/* Binning parameters file  path */
char			PhgCollimatorParamsFilePath[PHG_MAX_PARAM_FILES][PATH_LENGTH];			/* Collimator parameters file  path */
char			PhgDetectorParamsFilePath[PHG_MAX_PARAM_FILES][PATH_LENGTH];				/* Detector parameters file  path */
char			PhgTomoParamsFilePath[PHG_MAX_PARAM_FILES][PATH_LENGTH];	/* List of tomograph files */
} PhgRunTimeParamsTy;

LOCALE	PhgRunTimeParamsTy	PhgRunTimeParams;
LOCALE	LbUsFourByte		PhgNumTomoFiles;

/* MACROS */

/*********************************************************************************
*
*			Name:		PHGGetSineOfAccAngle
*
*			Summary:	Return sine of the acceptance angle.
*			Arguments:
*				
*			Function return: double, the sine of the acceptance angle.
*
*********************************************************************************/
#define PHGGetSineOfAccAngle() PhgRunTimeParams.Phg_SineOfAcceptanceAngle

/*********************************************************************************
*
*			Name:		PHGGetAccAngle
*
*			Summary:	Return the acceptance angle.
*			Arguments:
*				
*			Function return: double, the  acceptance angle.
*
*********************************************************************************/
#define PHGGetAccAngle() PhgRunTimeParams.Phg_AcceptanceAngle

/*********************************************************************************
*
*			Name:		PHGGetLengthOfScan
*
*			Summary:	Return length of the simulation.
*			Arguments:
*				
*			Function return: LbUsFourByte, the length of the simulation.
*
*********************************************************************************/
#define PHGGetLengthOfScan() PhgRunTimeParams.Phg_LengthOfScan

/* PROTOTYPES */
PhgEn_RunTimeParamsTy	PhgLookupRunTimeParamLabel(char *label);
Boolean					PhgGetRunTimeParams(void);
PhgEn_BinParamsTy  		PhgLookupBinParamLabel(char *label);
Boolean					PhgGetBinParams(char *ParamsName, PHG_BinParamsTy *binParams);
#undef LOCALE	
#endif /* PHG_PARAMS_HDR */

