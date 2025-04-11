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
*     Module Name: 			DetParams.h
*     Revision Number:		1.1
*     Date last revised:	August 2006
*     Programmer:			Steven Vannoy
*     Date Originated:		11, July, 1994
*
*     Module Overview:	This is the global include file for the DetParams.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*	  Global variables defined:   
*				DetRunTimeParams
*
*	  Global macros defined:
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		March 2010
*
*			Revision description:	Added block-position-in-tomo and
*						active layer info to DetBlockTomoBlockTy and
*						DetBlockTomoLayerTy.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		July 2007
*
*			Revision description:	Changed input of block positioning 
*						within rings to cylindrical coords.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		August 2006
*
*			Revision description:	Added block detector parameters.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:	Added time-of-flight parameters.
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
*********************************************************************************/
#ifndef DET_PARAMS_HDR
#define DET_PARAMS_HDR

#ifdef DET_PARAMS
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* PROGRAM CONSTANTS */	
#define PATH_LENGTH				256		/* Length of file name paths */
#define DET_MAX_TRANSAXIAL_CUTS	20
#define DET_MAX_AXIAL_CUTS		20
#define DET_MAX_LAYERS			3
#define	DET_MAX_BLOCKS			100		/* Per ring */
#define DET_MAX_RINGS			5
		
/* OPTION MACROS */
#define	DET_IsGaussBlurr()			(DetRunTimeParams[DetCurParams].IsGaussBlurr != 0)		/* Do we extend the target cylinder? */
#define	DET_IsDoHistory()			(DetRunTimeParams[DetCurParams].DoHistory)				/* Do we create a history file? */
#define	DET_IsDoCustomHistory()		(DetRunTimeParams[DetCurParams].DoCustomHistory)				/* Do we create a history file? */

/* PROGRAM TYPES */

/*	This list of enumerates specifies the type of detector being used.
	This allows one monster data structure to house all dector definitions.
	The values on this list must correspond, in order, to detRunTimeParamLabels
	in DetParams.c.
*/
typedef enum {	/* Our label table */
	DetEn_detector_type,
	DetEn_photon_tof_fwhm,
	DetEn_energy_resolution_percentage,
	DetEn_reference_energy_keV,
	DetEn_do_forced_interaction,
	DetEn_history_file,
	DetEn_history_params_file,
	DetEn_poly_axial_cut_pos_list,
	DetEn_poly_transaxial_cut_pos_list,
	DetEn_poly_axial_cut_start,
	DetEn_poly_axial_cut_end,
	DetEn_poly_transaxial_cut_start,
	DetEn_poly_transaxial_cut_end,
	DetEn_poly_scintillator,
	DetEn_poly_radial_depth,
	DetEn_poly_num_detectors_transaxially,
	DetEn_poly_num_detectors_axially,
	DetEn_poly_layer_info_list,
	DetEn_poly_num_layers,
	DetEn_poly_cut_material,
	DetEn_poly_axial_side_cover_thickness,
	DetEn_poly_axial_side_cover_material,
	DetEn_poly_trnsax_side_cvr_thickness,
	DetEn_poly_trnsax_side_cvr_material,
	DetEn_poly_front_cover_thickness,
	DetEn_poly_front_cover_material,
	DetEn_poly_transaxial_position,
	DetEn_poly_axial_min,
	DetEn_poly_axial_max,
	DetEn_poly_block_info_list,
	DetEn_poly_transaxial_length,
	DetEn_poly_axial_length,
	DetEn_poly_num_blocks,
	DetEn_poly_inner_radius,
	DetEn_poly_num_rings,
	DetEn_poly_ring_info_list,
	DetEn_plnr_num_layers,
	DetEn_plnr_layer_material,
	DetEn_plnr_layer_depth,
	DetEn_plnr_layer_is_active,
	DetEn_plnr_layer_info_list,
	DetEn_plnr_axial_length,
	DetEn_plnr_transaxial_length,
	DetEn_plnr_inner_radius,
	DetEn_plnr_num_views,
	DetEn_plnr_min_angle,
	DetEn_plnr_max_angle,
	DetEn_cyln_num_layers,
	DetEn_cyln_layer_material,
	DetEn_cyln_layer_inner_radius,
	DetEn_cyln_layer_outer_radius,
	DetEn_cyln_layer_is_active,
	DetEn_cyln_layer_info_list,
	DetEn_cyln_gap_material,
	DetEn_cyln_transaxial_gap,
	DetEn_cyln_axial_gap,
	DetEn_cyln_num_transaxial,
	DetEn_cyln_num_axial,
	DetEn_cyln_min_z,
	DetEn_cyln_max_z,
	DetEn_cyln_ring_info_list,
	DetEn_cyln_num_rings,
	DetEn_blocktomo_position_algorithm,
	DetEn_blocktomo_num_rings,
	DetEn_blocktomo_ring_description_list,
	DetEn_blocktomo_ring_parameter_file,
	DetEn_blocktomo_ring_axial_shift,
	DetEn_blocktomo_ring_transaxial_rotation,
	DetEn_ring_x_inner_radius,
	DetEn_ring_x_outer_radius,
	DetEn_ring_y_inner_radius,
	DetEn_ring_y_outer_radius,
	DetEn_ring_z_minimum,
	DetEn_ring_z_maximum,
	DetEn_ring_num_blocks_in_ring,
	DetEn_ring_block_description_list,
	DetEn_ring_block_parameter_file,
	DetEn_ring_block_radial_position,
	DetEn_ring_block_angular_position,
	DetEn_ring_block_z_position,
	DetEn_ring_block_transaxial_orientation,
	DetEn_block_reference_x,
	DetEn_block_reference_y,
	DetEn_block_reference_z,
	DetEn_block_x_minimum,
	DetEn_block_x_maximum,
	DetEn_block_y_minimum,
	DetEn_block_y_maximum,
	DetEn_block_z_minimum,
	DetEn_block_z_maximum,
	DetEn_block_num_layers,
	DetEn_block_layer_info_list,
	DetEn_block_layer_inner_x,
	DetEn_block_layer_outer_x,
	DetEn_block_layer_num_y_changes,
	DetEn_block_layer_y_change,
	DetEn_block_layer_num_z_changes,
	DetEn_block_layer_z_change,
	DetEn_block_layer_material_elements,
	DetEn_block_material_info_list,
	DetEn_block_material_index,
	DetEn_block_material_is_active,
	DetEn_randoms_history_file,
	DetEn_coincidence_timing_window_in_ns,
	DetEn_triples_processing_method,
	DetEn_NULL	/* must always end with null--loops terminate using this */
}DetEn_RunTimeParamsTy;

/* Detected position algorithm for block detectors */
typedef enum {
	DetEn_BlockSnapCentroidToXtalCenter,	/* use center of crystal nearest to energy-weighted centroid */
	DetEn_BlockUseCentroid,			/* use energy-weighted centroid */
	DetEn_BlockPosAlgoNULL			/* must always end with null--loops terminate using this */
} DetEn_BlockDetectedPositionAlgoTy;	

/* Triples processing method type */
typedef enum {
	DetEn_DeleteTriples	/* delete photons when three or more appear in single window */
} DetEn_TriplesMethodTy;	


/* This struct defines parameters for detectors within a block */
typedef struct {
	double	CutStart;	/* Beginning of detector cut, relative to block coordinates */
	double	CutEnd;		/* End of detector cut, relative to block coordinates */
} DetPositionTy;

/* This struct defines the material and active status of a detector crystal */
typedef struct {
	LbFourByte	Material;	/* The material index */
	Boolean		IsActive;	/* Active status */
} DetCrystalTy;

/* This struct defines the layer info for a polygonal detector system */
typedef struct {
	LbUsFourByte	NumMatAxially;
	LbUsFourByte	NumMatTransaxially;
	double			RadialDepth;
	DetPositionTy	AxialCutPos[DET_MAX_AXIAL_CUTS];
	DetPositionTy	TransaxialCutPos[DET_MAX_TRANSAXIAL_CUTS];
	DetCrystalTy	Crystals[DET_MAX_AXIAL_CUTS][DET_MAX_TRANSAXIAL_CUTS];
} DetPolyLayerInfoTy;

/*	This struct defines parameters for each block within the polygonal
	detector system
*/
typedef struct {
	double				TransaxialPos;
	LbUsFourByte		NumLayers;
	double				TransaxialLen;
	double				AxialMin;
	double				AxialMax;
	LbUsFourByte		FrontCvrMaterial;
	double				FrontCvrThickness;
	LbUsFourByte		TransSideCvrMaterial;
	double				TransSideCvrThickness;
	LbUsFourByte		AxialSideCvrMaterial;
	double				AxialSideCvrThickness;
	LbUsFourByte		CutMaterial;
	DetPolyLayerInfoTy	*LayerInfo;
} DetPolyBlockTy;

/*	This structure defines parameters for each ring within the polygonal
	detector system.
*/
typedef struct {
	LbUsFourByte	NumBlocks;
	double			MinAxialPos;
	double			MaxAxialPos;
	DetPolyBlockTy	*BlockInfo;
} DetPolRingTy;

/* This structure defines a polygonal detector system */
typedef struct {
	double			InnerRadius;
	LbUsFourByte	NumRings;
	DetPolRingTy	*RingInfo;
}DetPolyDetTy;

/* This struct defines parameters for each layer in a planar detector system */
typedef struct {
	double			LayerDepth;
	LbUsFourByte	LayerMaterial;
	Boolean			IsActive;
} DetPlnrLayerInfoTy;

/* This struct defines parameters for the planer detector system */
typedef struct {
	LbUsFourByte		NumLayers;
	double 				InnerRadius;
	double				TransaxialLength;
	double				AxialLength;
	LbFourByte			NumViews;
	double				MinAngle;
	double				MaxAngle;
	DetPlnrLayerInfoTy	*LayerInfo;
} DetPlnrDetTy;


/* This struct defines parameters for cylindrical detector system layers */
typedef struct {
	LbUsFourByte		LayerMaterial;
	double				InnerRadius;
	double				OuterRadius;
	Boolean				IsActive;
} DetCylnLayerInfoTy;

/* This struct defines parameters for cylindrical detector systems */
typedef struct {
	LbUsFourByte		NumLayers;
	DetCylnLayerInfoTy	*LayerInfo;
	double				MinZ;
	double				MaxZ;
	LbUsFourByte		NumAxial;
	LbUsFourByte		NumTransaxial;
	double				AxialGap;
	double				TransaxialGap;
	LbUsFourByte		GapMaterial;
} DetCylnRingTy;

typedef struct {
	LbUsFourByte		NumRings;
	DetCylnRingTy		*RingInfo;
} DetCylnDetTy;

/* The following struct define parameters for block tomograph detector systems */

/* The following structure defines the base element of a block detector tomograph.
 It is a right-rectangular box of either active or inactive homogeneous material.
 The dimensions are defined in the parent Layer and Block structures. */
typedef struct {
	LbFourByte			MaterialIndex;	/* index for attenuation */
	Boolean				IsActive;		/* Active (scintillating) status */
	LbFourByte			crystalNumInBlock;	/* for active elements, the SimSET-assigned
								crystal number within this block. */
	LbFourByte			crystalNumInTomo;	/* for active elements, the SimSET-assigned
								crystal number within tomograph. */
								
	/* the following three fields duplicate those in DetRunTimeParamsTy (and will
	 initially be copied from them), but will be  used to allow users to more
	 easily simulate resolution characteristics that vary across the block */
	double				PhotonTimeFWHM;	/* Photon time resolution (nanoseconds) */
	double				EnergyResolutionPercentage;		/* Energy resolution in percentage */
	double				ReferenceEnergy;	/* Energy resolution in percentage */
} DetBlockTomoElementTy;

/* The following struct defines one layer of a block. All coordinates are
 in the block's own coordinate system. */
typedef struct {
	double				InnerX;		/* the coord of face nearest tomo center */ 
	double				OuterX;		/* the coord of face furthest from tomo center */
	LbFourByte			NumYChanges;/* number of material changes in transaxial direction */
	double				*YChanges;	/* coords of material changes in transaxial direction */
	LbFourByte			NumZChanges;/* number of material changes in axial direction */
	double				*ZChanges;	/* coords of material changes in axial direction */
	LbUsFourByte		NumElements;/* number of material elements in layer */
	LbUsFourByte		NumActiveElementsInLayer;	/* number of active elements in the layer */
	DetBlockTomoElementTy	*ElementInfo; /* definition of each element */
} DetBlockTomoLayerTy;

/* The following struct defines one block of a tomograph All coordinates are
 in the block's own coordinate system. */
typedef struct {
	char				BlockParamFilePath[PATH_LENGTH];	/* Parameters for block file  path */
	double				RadialPosition;	/* cylindrical coord radial position of block within ring */
	double				AngularPosition;/* cylindrical coord angular position of block within ring (in radians) */
	double				AngularPositionDegrees;	/* angular position in degrees (as input) */
	double				ZPosition;	/* z position of block within ring */
	double				TransaxialOrientation;	/* transaxial direction, relative to vector
							to tomograph center line, the block points  (in degrees) */
	double				AngularPositionTomo;	/* angular position of block within tomograph in radians
							(AngularPosition and ring TransaxialRotation) */
	double				XPositionTomo;	/* x position of block within tomograph (computed from cyl. coord
							and ring TransaxialRotation) */
	double				YPositionTomo;	/* y position of block within tomograph (computed from cyl. coord
							and ring TransaxialRotation) */
	double				ZPositionTomo;	/* z position of block within tomograph (computed from cyl. coord
							and ring AxialShift) */
	double				BlockFaceAngle;	/* the angle a direction vector perpendicular to the 
							block face makes with the x-axis, in radians (calculated from AngularPositionTomo
							and TransaxialOrientation) */
	double				CosBlockFace;	/* cosine of the block face angle */
	double				SinBlockFace;	/* sine of the block face angle */
	double				XRef;	/* x coord of point within block that above positioning refers to */
	double				YRef;	/* y coord of point within block that above positioning refers to */
	double				ZRef;	/* z coord of point within block that above positioning refers to */
	double				XMin;	/* minimum x extent of block */
	double				XMax;	/* maximum x extent of block */
	double				YMin;	/* minimum y extent of block */
	double				YMax;	/* maximum y extent of block */
	double				ZMin;	/* minimum z extent of block */
	double				ZMax;	/* maximum z extent of block */
	LbUsFourByte		NumLayers;	/* number of layers in block in x-direction. */
	LbUsFourByte		NumActiveElementsInBlock;	/* number of active material elements in block */
	LbFourByte			FirstActiveElementInBlock;	/* the crystalNumInTomo of the first
										active element in block */
	LbUsFourByte		NumActiveLayers;			/* number of layers with > 0 active elements */
	Boolean				*ActiveLayers;				/* true for layers with > 0 active elements */
	DetBlockTomoLayerTy	*LayerInfo;	/* definition of each layer */
} DetBlockTomoBlockTy;

/* The following struct defines one axial ring of a tomograph */
typedef struct {
	char				RingParamFilePath[PATH_LENGTH];	/* Parameters for ring file  path */
	double				AxialShift;	/* z-shift to apply to defined ring */
	double				TransaxialRotation; /* transaxial rotation to apply to defined ring */
	/* The following x- and y- radii are used to define elliptical cylinders.  All
	 blocks must lay, in their entirety, outside the inner cylinder and inside
	 the outer cylinder. Similarly they must lay above the min Z and below the max Z */
	double				XInnerRadius;
	double				XOuterRadius;
	double				YInnerRadius;
	double				YOuterRadius;
	double				MinZ;
	double				MaxZ;
	LbUsFourByte		NumBlocks;	/* number of blocks making up this ring */
	LbUsFourByte		NumActiveElementsInRing;	/* number of active material elements in ring */
	LbFourByte			FirstActiveElementInRing;	/* the crystalNumInTomo of the first
										active element in ring */
	DetBlockTomoBlockTy	*BlockInfo;	/* definition of each block */
} DetBlockTomoRingTy;

/* Using the above structures, the following structure defines a block tomograph */
typedef struct {
	LbUsFourByte		NumRings;	/* number of rings in tomograph */
	LbUsFourByte		NumActiveElementsInTomo;	/* number of active material elements in tomograph */
	DetBlockTomoRingTy	*RingInfo;	/* definition of each ring */
	DetEn_BlockDetectedPositionAlgoTy	BlockDetectedPositionAlgo;	/* Detected position algorithm for block detectors */
} DetBlockTomoDetTy;

	
/* Runtime params */
/* The following struct contains the run-time parameters for the detector.
	Note that the first field is a flag indicating the struct has been
	initialized. When writing a history file, the module PhoHFile checks
	to see if detection is being done. If it is not, then the module
	clears out the run time parameters for detection so that the header
	will reflect the fact that the detector parameter data has not been used.
*/
typedef struct {
Boolean					Initialized;							/* Initialization flag */
Boolean					DoForcedInteraction;					/* force at least one interaction to occur in detector? */
Boolean					DoHistory;								/* Are we creating a history file */
Boolean					DoCustomHistory;						/* Are we creating a history file */
Boolean					DoRandomsProcessing;					/* Create randoms processed list mode? */
double					PhotonTimeFWHM;							/* Photon time resolution (nanoseconds) */
double					EnergyResolutionPercentage;				/* Energy resolution in percentage */
double					ReferenceEnergy;						/* Energy resolution in percentage */
double					CoincidenceTimingWindowNS;				/* coincidence timing window in nanoseconds */
DetEn_TriplesMethodTy	TriplesMethod;							/* method for triples processing */
char					DetHistoryFilePath[PATH_LENGTH];		/* Detected History file  path */
char					DetHistoryParamsFilePath[PATH_LENGTH];	/* Detected History parameters file  path */
char					DetRandomsHistoryFilePath[PATH_LENGTH];	/* Randoms-added history parameters file  path */
DetEn_DetectorTypeTy	DetectorType;							/* The type of detector */
DetPolyDetTy			PolyDetector;
DetPlnrDetTy			PlanarDetector;
DetCylnDetTy			CylindricalDetector;
DetBlockTomoDetTy		BlockTomoDetector;
} DetRunTimeParamsTy;

/* PROGRAM GLOBALS */
LOCALE	DetRunTimeParamsTy	DetRunTimeParams[PHG_MAX_PARAM_FILES];
LOCALE	LbFourByte			DetNumParams;
LOCALE	LbFourByte			DetCurParams;

/* MACROS */

/* PROTOTYPES */
DetEn_RunTimeParamsTy	DetLookupRunTimeParamLabel(char *label);
Boolean					DetGetRunTimeParams(void);
#undef LOCALE
#endif /* DET_PARAMS_HDR */
