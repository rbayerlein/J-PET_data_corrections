/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992 - 1997 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			phg.h
*     Revision Number:		1.0
*     Date last revised:	August 19, 1997
*     Programmer:			Steven Vannoy
*     Date Originated:		August 28, 1992
*
*     Module Overview:	This is the global include file for the phg.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*	  Global variables defined:   
*
*	  Global macros defined:
*				PHGGetSineOfAccAngle
*
*********************************************************************************/
#ifndef PHG_MAIN_HDR
#define PHG_MAIN_HDR

#include <signal.h>

#ifdef DGUX
	#include <siginfo.h>
	#include <ieeefp.h>
	#include <ucontext.h>
#endif

#ifdef PHG_MAIN
	#define	LOCALE
#else
	#define LOCALE	extern
#endif


/* PROGRAM CONSTANTS */

/* OPTION MACROS */
#define PHG_IsDebugOptions()			LbFgIsSet(PhgOptions, LBFlag0)	/* Did user supply debug options? */

#define	PHG_NumFlags	1												/* Number of flags defined */

/* DEBUGGING OPTIONS */
#define PHGDEBUG_FixedDirection()			LbFgIsSet(PhgDebugOptions, LBFlag0) /* Should we pick a fixed direction */
#define PHGDEBUG_BinFirstDetPosition()		LbFgIsSet(PhgDebugOptions, LBFlag1) /* Should we bin up first interaction points */
#define PHGDEBUG_ReadFromRandSeedFile()		LbFgIsSet(PhgDebugOptions, LBFlag2) /* Should we read from the random seed file? */
#define PHGDEBUG_WriteToRandSeedFile()		LbFgIsSet(PhgDebugOptions, LBFlag3) /* Should we write to the random seed file? */
#define PHGDEBUG_ResetRandSeedFile()		LbFgIsSet(PhgDebugOptions, LBFlag4) /* Should we reset the random seed file? */
#define PHGDEBUG_BinCentroidPosition()		LbFgIsSet(PhgDebugOptions, LBFlag5) /* Should we bin the centroid position within the detector? */
#define PHGDEBUG_DetAbsorbOnly()			LbFgIsSet(PhgDebugOptions, LBFlag6) /* Should we detector make all interactions absorptions? */		
#define PHGDEBUG_ObjCohOnly()				LbFgIsSet(PhgDebugOptions, LBFlag7) /* Should object interactions be coherent only */

/* PROGRAM TYPES */
/* The following is simply the local globals that used to be in PhgBin.c
that needed to be consolidated in order to allow multiple binning parameter
settings, 'fields' is hence the monicker the all fall under.
*/

typedef struct {
	void			*countImage;		/* The count image */
	void			*weightImage;		/* The weight image */
	void			*weightSquImage;	/* The weight squared image */
} PHG_BinDataTy;

typedef struct {
Boolean				IsInitialized;			/* Our initialization flag */
Boolean				ParamsIsInitialized;			/* Our initialization flag */
LbUsEightByte		TotBluePhotons;					/* Counter for output report */
LbUsEightByte		TotPinkPhotons;					/* Counter for output report */
LbUsEightByte		AccBluePhotons;					/* Counter for output report */
LbUsEightByte		AccPinkPhotons;					/* Counter for output report */
LbUsEightByte		NumCoincidences;				/* Counter for output report */
LbUsEightByte		NumAcceptedCoincidences;		/* Counter for output report */
double				AccCoincidenceWeight;			/* Counter for output report */
double				AccCoincidenceSquWeight;		/* Counter for output report */
double				StartAccCoincidenceWeight;		/* Counter for starting value output report */
double				StartAccCoincidenceSquWeight;	/* Counter for starting value  output report */
PhoHFileHdrTy		CountImgHdr;					/* Header for count image file */
PhoHFileHdrTy		WeightImgHdr;					/* Header for weight image file */
PhoHFileHdrTy		WeightSquImgHdr;				/* Header for weight squared image file */
LbHdrHkTy			CountImgHdrHk;					/* Hook for count image file header */
LbHdrHkTy			WeightImgHdrHk;					/* Hook for weight image file header */
LbHdrHkTy			WeightSquImgHdrHk;				/* Hook for weight squared image fileheader  */
FILE				*CountFile;						/* File pointer to count image */
FILE				*WeightFile;					/* File pointer to weight image */
FILE				*WeightSquFile;					/* File pointer to weight squared image */
double				WeightRatio;					/* Scale factor for adding to existing images feature */
PhoHFileHkTy		historyFileHk;					/* Hook to history file */
} PHG_BinFieldsTy;

/* PROGRAM GLOBALS */
LOCALE	PHG_BinParamsTy			PhgBinParams[PHG_MAX_PARAM_FILES];		/* Our parameters */
LOCALE	PHG_BinDataTy			PhgBinData[PHG_MAX_PARAM_FILES];		/* Our histogram data */
LOCALE	PHG_BinFieldsTy			PhgBinFields[PHG_MAX_PARAM_FILES];		/* Various items used for binning */
LOCALE	LbUsFourByte			PhgNumBinParams;							/* The number of binning parameters */

#ifdef PHG_DEBUG
LOCALE		FILE		*PhgDebugDumpFile;				/* For for dumping data structers */
#endif
LOCALE	LbUsFourByte	PhgOptions;						/* User specified options */
LOCALE	LbUsFourByte	PhgDebugOptions;				/* Options affecting debug output */

LOCALE	char			PhgExecutionDateStr[33];		/* String for execution date conversion */

/* PROTOTYPES */
Boolean	PhgRun(int argc, char *argv[]);
void	PhgAbort(char *abortStr, Boolean dumpTrkPhoton);
void	PhgDumpPhoton(PHG_TrackingPhoton *trkPhotonPtr);
void	PhgDumpDecay(PHG_Decay *decayPtr);
void	PhgPrintParams(int argc, char *argv[], Boolean randFromClock);


#ifdef DGUX

void					PhgSigHandler(int sigNum, siginfo_t *sigInfoPtr,
							ucontext_t *sigContextPtr);
#endif

#undef LOCALE
#endif /* PHG_MAIN_HDR */
