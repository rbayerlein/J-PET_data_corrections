/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1994-2001 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			ColParams.h
*     Revision Number:		1.0
*     Date last revised:	25 December 2000
*     Programmer:			Steven Vannoy
*     Date Originated:		11, July, 1994
*
*     Module Overview:	This is the global include file for the ColParams.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*	  Global variables defined:   
*				ColRunTimeParams
*
*	  Global macros defined:
*********************************************************************************/
#ifndef COL_PARAMS_HDR
#define COL_PARAMS_HDR

#ifdef COL_PARAMS
	#define	LOCALE
#else
	#define LOCALE	extern
#endif


/* PROGRAM CONSTANTS */
#define PATH_LENGTH				256				/* Length of file name paths */

/* OPTION MACROS */
#define	COL_IsExtendTarget()			(ColRunTimeParams[ColCurParams].CollimatorDepth != 0)		/* Do we extend the target cylinder? */

/* The collimator module should only create a history file if it is being called from the PHG
	or the BIN.PHG utility. If not, it is probably being called by the DETECTOR utility which
	may use the history file created by the collimator.
*/

#define	COL_IsDoHistory()				(ColRunTimeParams[ColCurParams].DoHistory)		/* Do we create a history file? */
#define	COL_IsDoCustomHistory()			(ColRunTimeParams[ColCurParams].DoCustomHistory)	/* Do we create a custom history file? */

/* PROGRAM TYPES */

/* The following enumeration defines all of the collimator paramter
	labels used in the GetParams functions
*/
typedef enum {
	ColEn_collimator_type,
	ColEn_collimator_depth,
	ColEn_layers_list,
	ColEn_segment_list,
	ColEn_segment,
	ColEn_seg_type,
	ColEn_material,
	ColEn_inner_radius,
	ColEn_outer_radius,
	ColEn_inner_min_z,
	ColEn_outer_min_z,
	ColEn_inner_max_z,
	ColEn_outer_max_z,
	ColEn_history_file,
	ColEn_history_params_file,
	ColEn_hole_geometry,		/* UNC Param */
	ColEn_focal_length,			/* UNC Param */
	ColEn_radius_of_rotation,	/* UNC Param */
	ColEn_thickness,			/* UNC Param */
	ColEn_hole_radius,			/* UNC Param */
	ColEn_septal_thickness,		/* UNC Param */
	ColEn_min_z,				/* UNC Param */
	ColEn_max_z,				/* UNC Param */
	ColEn_start_angle,			/* UNC Param */
	ColEn_stop_angle,			/* UNC Param */
	ColEn_sum_all_views,		/* UNC Param */
	ColEn_num_views,			/* UNC Param */
	ColEn_depth,
	ColEn_num_slats,			/* Slat Param */
	ColEn_start,				/* Slat Param */
	ColEn_end,					/* Slat Param */
	ColEn_num_layers,			/* Slat Param */
	ColEn_NULL			
}ColEn_RunTimeParamsTy;


/* PROGRAM GLOBAL TYPES */
/* The following struct contains the run-time parameters for the collimator.
	Note that the first field is a flag indicating the struct has been
	initialized. When writing a history file, the module PhoHFile checks
	to see if collimation is being done. If it is not, then the module
	clears out the run time parameters for collimation so that the header
	will reflect the fact that the collimator parameter data has not been used.
*/
typedef struct {
Boolean					Initialized;							/* Indicates structure is initialized */
Boolean					DoHistory;								/* Do we create a history file */
Boolean					DoCustomHistory;						/* Do we create a custom history file */
char					ColHistoryFilePath[PATH_LENGTH];		/* Collimated History file  path */
char					ColHistoryParamsFilePath[PATH_LENGTH];	/* Collimated History parameters file  path */
ColEn_CollimatorTypeTy	ColType;								/* What type of collimator */
Col_SimplePETTy			SimplePETCol;							/* A "simple" PET collimator  */
Col_Monte_Carlo_PET_Ty	*MCPETCol;								/* A Monte Carlo PET collimtaor */
Col_UNC_SPECT_Ty		UNCSPECTCol;							/* A UNC style SPECT collimator */
Col_Slat_Ty				SlatCol;								/* A slat collimator */
} ColRunTimeParamsTy;

/* PROGRAM GLOBALS */
LOCALE	ColRunTimeParamsTy	ColRunTimeParams[PHG_MAX_PARAM_FILES];
LOCALE	LbFourByte			ColNumParams;
LOCALE	LbFourByte			ColCurParams;

/* MACROS */

/* PROTOTYPES */
ColEn_RunTimeParamsTy	ColLookupRunTimeParamLabel(char *label);
Boolean					ColGetRunTimeParams(void);
#undef LOCALE
#endif
