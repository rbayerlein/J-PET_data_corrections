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
*			Module Name:		PhgHdr.h
*			Revision Number:	1.1
*			Date last revised:	2007
*			Programmer:			Steven Vannoy
*			Date Originated:	24 August 1992
*
*			Module Overview:	Definitions for PhgHdr.c.
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
*			Revision description:
*						- support for randoms and eight-byte number of decays
*						- binning by random-state and crystal number.
*						- binning PET data as SPECT
*
*********************************************************************************/
#ifndef PHG_HEADER_HDR
#define PHG_HEADER_HDR

#ifdef PHG_HEADER
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/*	CONSTANTS	*/
#define PHG_HDR_HEADER_SIZE		32768		/* The size of our headers */
#define PHG_HDR_HEADER_VERSION	2.0			/* The current header version */

#define HDR_PHG_EVENTS_TO_SIM_ID			10000 /* replaces 10004, tag for 4 byte int */
#define	HDR_PHG_HEADER_SIZE_ID				10001
#define HDR_PHG_HEADER_KIND_ID				10002
#define HDR_PHG_HEADER_VERS_ID				10003
#define HDR_PHG_EVENTS_TO_SIM_OLD4BYTE_ID	10004 /* replaced by 10000, tag for 8 byte int */
#define HDR_PHG_RAND_SEED_ID				10005
#define HDR_PHG_LEN_OF_SCAN_ID				10006
#define HDR_PHG_ACC_ANGLE_ID				10007
#define HDR_PHG_ACC_ANGLE_SINE_ID			10008
#define HDR_PHG_MIN_ENERGY_ID				10009
#define HDR_PHG_MAX_WW_RATIO_ID				10010
#define HDR_PHG_ISOTOPE_ID					10011
#define HDR_PHG_PHOTON_ENERGY_KEV_ID		10012
#define HDR_PHG_POSITRON_ENERGY_ID			10013
#define HDR_PHG_IS_FORCED_DETECTION_ID		10014
#define HDR_PHG_IS_STRATIFICATION_ID		10015
#define HDR_PHG_IS_NON_ABSORPTION_ID		10016
#define HDR_PHG_IS_SPECT_ID					10017	/* see also 10049-10053 below */
#define HDR_PHG_IS_PET_COINC_ONLY_ID		10018	/* see also 10049-10053 below */
#define HDR_PHG_IS_HISTORY_FILE_ID			10019
#define HDR_PHG_IS_ADJ_FOR_POS_RANGE_ID		10020
#define HDR_PHG_IS_ADJ_FOR_COLIN_ID			10021
#define HDR_PHG_IS_COMPUTED_PROD_TBL_ID		10022
#define HDR_PHG_IS_VOXEL_PTSRC_ID			10023
#define HDR_PHG_IS_BIN_ON_THE_FLY_ID		10024
#define HDR_PHG_IS_COL_ON_THE_FLY_ID		10025
#define HDR_PHG_IS_DET_ON_THE_FLY_ID		10026

#define HDR_PHG_PHGPARAM__ID				10027
#define HDR_PHG_SUBOBJACTIVITYINDEX_ID		10028
#define HDR_PHG_SUBOBJACTIVITYTABLE_ID		10029
#define HDR_PHG_SUBOBJACTINDEXTRANS_ID		10030
#define HDR_PHG_SUBOBJACTIMG_ID				10031
#define HDR_PHG_SUBOBJATTENINDEX_ID			10032
#define HDR_PHG_SUBOBJATTENTABLE_ID			10033
#define HDR_PHG_SUBOBJATTINDEXTRANS_ID		10034
#define HDR_PHG_SUBOBJATTIMG_ID				10035
#define HDR_PHG_PRODTBLINPUTTABLE_ID		10036
#define HDR_PHG_PRODTBLOUTPUTTABLE_ID		10037
#define HDR_PHG_PHOHSTAT_ID					10038
#define HDR_PHG_PHOHFILEHISTORY_ID			10039
#define HDR_PHG_BINPARAMS_ID				10040
#define HDR_PHG_COLLIMATORPARAMS_ID			10041
#define HDR_PHG_DETECTORPARAMS_ID			10042
#define HDR_PHG_BINWEIGHTIMG_ID				10043
#define HDR_PHG_BINWEIGHTSQUIMG_ID			10044
#define HDR_PHG_BINCOUNTIMG_ID				10045
#define HDR_PHG_MIN_WW_RATIO_ID				10046

#define HDR_PHG_IS_VOXEL_LNSRC_ID			10047
#define HDR_PHG_IS_POLARIZATION_ID			10048

#define HDR_PHG_IS_PET_COINC_PLUS_SINGLES_ID		10049	/* see also 10017-10018 above */
#define HDR_PHG_IS_MULTIEMISSION_ID			10050	/* see also 10017-10018 above */
#define HDR_PHG_IS_PET_ID				10051	/* see also 10017-10018 above */
#define HDR_PHG_IS_PET_COINC_ID				10052	/* see also 10017-10018 above */
#define HDR_PHG_IS_PET_SINGLES_ID			10053	/* see also 10017-10018 above */
#define HDR_PHG_IS_CALC_EVENTS_TO_SIM_ID	10054

#define	HDR_COL_IS_INITIALIZED_ID			20001
#define HDR_COL_IS_DO_HISTORY_ID			20002
#define HDR_COL_HISTORY_ID					20003
#define HDR_COL_TYPE_ID						20004
#define HDR_COL_DEPTH_ID					20005
#define HDR_COL_UNC_HOLE_GEOM_ID			20006
#define HDR_COL_UNC_RAD_OF_ROTATION_ID		20007
#define HDR_COL_UNC_THICKNESS_ID			20008
#define HDR_COL_UNC_HOLE_RADIUS_ID			20009
#define HDR_COL_UNC_SEPTAL_THICKNESS_ID		20010
#define HDR_COL_UNC_MIN_Z_ID				20011
#define HDR_COL_UNC_MAX_Z_ID				20012
#define HDR_COL_UNC_START_ANGLE_ID			20013
#define HDR_COL_UNC_STOP_ANGLE_ID			20014
#define HDR_COL_UNC_SUM_ALL_VIEWS_ID		20015
#define HDR_COL_UNC_NUM_VIEWS_ID			20016

#define HDR_DET_IS_INITIALIZED_ID			30001
#define HDR_DET_IS_DO_HISTORY_ID			30002
#define HDR_DET_HISTORY_ID					30003
#define	HDR_DET_ENERGY_RESOLUTION_PER_ID	30004
#define HDR_DET_REFERENCE_ENERGY_ID			30005
#define HDR_DET_TYPE_ID						30006
#define HDR_DET_IS_FORCED_INTERACTION_ID	30007
#define HDR_DET_PHOTON_TIME_FWHM_ID			30008
#define HDR_DET_COINC_TIMING_WINDOW_ID		30009
#define HDR_DET_TRIPLES_METHOD_ID			30010
#define HDR_DET_RANDOMS_HISTORY_FILE_ID		30011

#define HDR_CYL_TARGET_RADIUS_ID			40001
#define HDR_CYL_TARGET_MIN_Z_ID				40002
#define HDR_CYL_TARGET_MAX_Z_ID				40003
#define HDR_CYL_TARGET_CENTER_X_ID			40004
#define HDR_CYL_TARGET_CENTER_Y_ID			40005

#define HDR_CYL_CRITZONE_RADIUS_ID			40006
#define HDR_CYL_CRITZONE_MIN_Z_ID			40007
#define HDR_CYL_CRITZONE_MAX_Z_ID			40008
#define HDR_CYL_CRITZONE_CENTER_X_ID		40009
#define HDR_CYL_CRITZONE_CENTER_Y_ID		40010

#define HDR_CYL_OBJCYL_RADIUS_ID			40011
#define HDR_CYL_OBJCYL_MIN_Z_ID				40012
#define HDR_CYL_OBJCYL_MAX_Z_ID				40013
#define HDR_CYL_OBJCYL_CENTER_X_ID			40014
#define HDR_CYL_OBJCYL_CENTER_Y_ID			40015

#define HDR_CYL_LIMITCYL_RADIUS_ID			40016
#define HDR_CYL_LIMITCYL_MIN_Z_ID			40017
#define HDR_CYL_LIMITCYL_MAX_Z_ID			40018
#define HDR_CYL_LIMITCYL_CENTER_X_ID		40019
#define HDR_CYL_LIMITCYL_CENTER_Y_ID		40020

#define	HDR_BIN1_NUM_SIMULATIONS_ID			50001
#define HDR_BIN1_EVENTS_TO_SIMULATION_ID	50002
#define HDR_BIN1_NUM_PHOTONS_ID				50003
#define HDR_BIN1_NUM_DECAYS_ID				50004

#define HDR_BIN_NUM_Z_BINS_ID				60001
#define HDR_BIN_NUM_PA_BINS_ID				60002
#define HDR_BIN_NUM_TD_BINS_ID				60003
#define HDR_BIN_NUM_AA_BINS_ID				60004
#define HDR_BIN_NUM_E1_BINS_ID				60005
#define HDR_BIN_NUM_E2_BINS_ID				60006
#define HDR_BIN_NUM_S1_BINS_ID				60007
#define HDR_BIN_NUM_S2_BINS_ID				60008
#define HDR_BIN_NUM_IMAGE_BINS_ID			60009
#define HDR_BIN_SCATTER_PARAM_ID			60010
#define HDR_BIN_MIN_Z_ID					60011
#define HDR_BIN_MAX_Z_ID					60012
#define HDR_BIN_MIN_PA_ID					60013
#define HDR_BIN_MAX_PA_ID					60014
#define HDR_BIN_MIN_TD_ID					60015
#define HDR_BIN_MAX_TD_ID					60016
#define HDR_BIN_MIN_AA_ID					60017
#define HDR_BIN_MAX_AA_ID					60018
#define HDR_BIN_MIN_ENERGY_ID				60019
#define HDR_BIN_MAX_ENERGY_ID				60020
#define HDR_BIN_MIN_SCAT_ID					60021
#define HDR_BIN_MAX_SCAT_ID					60022
#define HDR_BIN_ADD_TO_EXISTING_ID			60023
#define HDR_BIN_DO_COUNTS_ID				60024
#define HDR_BIN_DO_WEIGHTS_ID				60025
#define HDR_BIN_DO_WEIGHTS_SQU_ID			60026
#define HDR_BIN_Z_RANGE_ID					60027
#define HDR_BIN_E_RANGE_ID					60028
#define HDR_BIN_S_RANGE_ID					60029
#define HDR_BIN_TD_RANGE_ID					60030
#define HDR_BIN_SCATTER2CI_SIZE_ID			60031
#define HDR_BIN_SCATTER2WI_SIZE_ID			60032
#define HDR_BIN_SCATTER2WIS_SIZE_ID			60033
#define HDR_BIN_SCATTER1CI_SIZE_ID			60034
#define HDR_BIN_SCATTER1WI_SIZE_ID			60035
#define HDR_BIN_SCATTER1WIS_SIZE_ID			60036
#define HDR_BIN_ENERGY2CI_SIZE_ID			60037
#define HDR_BIN_ENERGY2WI_SIZE_ID			60038
#define HDR_BIN_ENERGY2WIS_SIZE_ID			60039

#define HDR_BIN_ENERGY1CI_SIZE_ID			60040
#define HDR_BIN_ENERGY1WI_SIZE_ID			60041
#define HDR_BIN_ENERGY1WIS_SIZE_ID			60042
#define HDR_BIN_AACI_SIZE_ID				60043
#define HDR_BIN_AAWI_SIZE_ID				60044
#define HDR_BIN_AAWIS_SIZE_ID				60045
#define HDR_BIN_TDCI_SIZE_ID				60046
#define HDR_BIN_TDWI_SIZE_ID				60047
#define HDR_BIN_TDWIS_SIZE_ID				60048
#define HDR_BIN_PACI_SIZE_ID				60049
#define HDR_BIN_PAWI_SIZE_ID				60050
#define HDR_BIN_PAWIS_SIZE_ID				60051

#define HDR_BIN_Z2CI_SIZE_ID				60052
#define HDR_BIN_Z2WI_SIZE_ID				60053
#define HDR_BIN_Z2WIS_SIZE_ID				60054

#define HDR_BIN_Z1CI_SIZE_ID				60055
#define HDR_BIN_Z1WI_SIZE_ID				60056
#define HDR_BIN_Z1WIS_SIZE_ID				60057

#define HDR_BIN_COUNT_IMG_SIZE_ID			60058
#define HDR_BIN_WEIGHT_IMG_SIZE_ID			60059
#define HDR_BIN_WEIGHT_SQU_IMG_SIZE_ID		60060

#define HDR_BIN_WEIGHT_IMG_TYPE_ID			60061
#define HDR_BIN_COUNT_IMG_TYPE_ID			60062
#define HDR_BIN_DIMENSION_ORDER_ID			60063	

#define HDR_BIN_PHICI_SIZE_ID				60064
#define HDR_BIN_PHIWI_SIZE_ID				60065
#define HDR_BIN_PHIWIS_SIZE_ID				60066

#define HDR_BIN_THETACI_SIZE_ID				60067
#define HDR_BIN_THETAWI_SIZE_ID				60068
#define HDR_BIN_THETAWIS_SIZE_ID			60069

#define HDR_BIN_XRCI_SIZE_ID				60070
#define HDR_BIN_XRWI_SIZE_ID				60071
#define HDR_BIN_XRWIS_SIZE_ID				60072

#define HDR_BIN_YRCI_SIZE_ID				60073
#define HDR_BIN_YRWI_SIZE_ID				60074
#define HDR_BIN_YRWIS_SIZE_ID				60075

#define HDR_BIN_NUM_PHI_BINS_ID				60076
#define HDR_BIN_NUM_THETA_BINS_ID			60077
#define HDR_BIN_NUM_XR_BINS_ID				60078
#define HDR_BIN_NUM_YR_BINS_ID				60079
#define HDR_BIN_MIN_THETA_ID				60080
#define HDR_BIN_MAX_THETA_ID				60081
#define HDR_BIN_MIN_PHI_ID					60082
#define HDR_BIN_MAX_PHI_ID					60083
#define HDR_BIN_MIN_XR_ID					60084
#define HDR_BIN_MAX_XR_ID					60085
#define HDR_BIN_MIN_YR_ID					60086
#define HDR_BIN_MAX_YR_ID					60087
#define HDR_BIN_DOSSRB_ID					60088
#define HDR_BIN_DOMSRB_ID					60089
#define HDR_BIN_IMG_RADIUS_ID				60090
#define HDR_BIN_DET_RADIUS_ID				60091
#define HDR_BIN_WTIMG_PATH_ID				60092
#define HDR_BIN_WTSIMG_PATH_ID				60093
#define HDR_BIN_CTIMG_PATH_ID				60094
#define HDR_BIN_THETA_RANGE_ID				60095
#define HDR_BIN_PHI_RANGE_ID				60096
#define HDR_BIN_XR_RANGE_ID					60097
#define HDR_BIN_YR_RANGE_ID					60098

#define HDR_BIN_SUM_WEIGHTS_ID				60099
#define HDR_BIN_SUM_WEIGHTS_SQ_ID			60100

#define HDR_BIN_PET_AS_SPECT_ID             60101

#define HDR_BIN_NUM_TOF_BINS_ID				60102
#define HDR_BIN_MIN_TOF_ID					60103
#define HDR_BIN_MAX_TOF_ID					60104
#define HDR_BIN_TOFCI_SIZE_ID				60105
#define HDR_BIN_TOFWI_SIZE_ID				60106
#define HDR_BIN_TOFWIS_SIZE_ID				60107
#define HDR_BIN_DO_CRYSTAL_BINS_ID			60108
#define HDR_BIN_NUM_CRYSTAL_BINS_ID			60109
#define HDR_BIN_CRYSTAL1CI_SIZE_ID			60110
#define HDR_BIN_CRYSTAL1WI_SIZE_ID			60111
#define HDR_BIN_CRYSTAL1WIS_SIZE_ID			60112
#define HDR_BIN_CRYSTAL2CI_SIZE_ID			60113
#define HDR_BIN_CRYSTAL2WI_SIZE_ID			60114
#define HDR_BIN_CRYSTAL2WIS_SIZE_ID			60115

#define	HDR_ATT_HAS_BEEN_CORRECTED_ID		70001

#define HDR_HISTORY_FILE_IS_SORTED_ID		70101
#define HDR_PHG_TIMESORT_PARAMS_ID			70102

#define HDR_HISTORY_FILE_IS_RANDOMS_ADDED_ID	70201
#define HDR_PHG_ADDRAND_PARAMS_ID			70202

/* GLOBAL TYPES */
/* PROTOTYPES */

Boolean		PhgHdrCreateRunTimeHeader(PhoHFileHdrKindTy hdrKind,
				PhoHFileHdrTy *hdrTyPtr, PHG_BinParamsTy *binParams);
Boolean		PhgHdrGtParams(FILE *headerFile, PhoHFileHdrTy *paramsPtr, LbHdrHkTy *headerHk);
void		PhgHdrFrHeader(LbHdrHkTy *headerHk);
Boolean		PhgHdrMkHeader(FILE *headerFile, PhoHFileHdrTy *paramsPtr,  LbHdrHkTy *headerHk);
Boolean		PhgHdrStFields(PhoHFileHdrTy *paramsPtr, LbHdrHkTy *headerHk);
Boolean		PhgHdrUpHeader(FILE *headerFile, PhoHFileHdrTy *paramsPtr, LbHdrHkTy *headerHk);
Boolean		PhgHdrGtAttenuationCorrected(LbHdrHkTy *headerHk, Boolean *attCorrected);
Boolean		PhgHdrStAttenuationCorrected(LbHdrHkTy *headerHk);
Boolean		PhgHdrStFile(LbHdrHkTy *headerHk, FILE *headerFl);
Boolean 	PhgHdrGtHeaderHk(FILE *headerFl, LbHdrHkTy *headerHk);
Boolean		PhgHdrGtField(LbHdrHkTy *headerHk, LbUsFourByte fieldID, void *fieldData);
LbFourByte	PhgHdrGtFieldSize(LbUsFourByte fieldID);


#undef LOCALE
#endif /* PHG_HEADER_HDR */
