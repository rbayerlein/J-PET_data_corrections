/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1993-2004 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		PhgBin.h
*			Revision Number:	1.1
*			Date last revised:	2007
*			Programmer:			Steven Vannoy
*			Date Originated:	28 July 1993
*
*			Module Overview:	Definitions for PhgBin.c.
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
*						- Moved binning option type to PhgParams.h
*						- new prototypes
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:	Added support for TOF binning.
*
*********************************************************************************/
#ifndef PHG_BIN_HDR
#define PHG_BIN_HDR

#ifdef PHG_BIN
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */
#define 		PHG_BIN_WEIGHT_TYPE_R4	2	/* Four byte real */
#define 		PHG_BIN_WEIGHT_TYPE_R8	3	/* Eight byte real */

#define 		PHG_BIN_COUNT_TYPE_I1	0	/* One byte integer */
#define 		PHG_BIN_COUNT_TYPE_I2	1	/* Two byte integer */
#define 		PHG_BIN_COUNT_TYPE_I4	2	/* Four byte integer */
#define			PATH_LENGTH		256			/* Length of file name paths */

/* TYPES */


/* GLOBALS */

/* PROTOTYPES */
Boolean		PhgBinInitParams(char *ParamsName, PHG_BinParamsTy *binParams, PHG_BinDataTy *binData, PHG_BinFieldsTy *binFields);
Boolean		PhgBinInitialize(char *ParamsName, PHG_BinParamsTy *binParams, PHG_BinDataTy *binData, PHG_BinFieldsTy *binFields);
void		PhgBinPrintParams(char *ParamsName, PHG_BinParamsTy *binParams, PHG_BinFieldsTy *binFields);
void		PhgBinPrintReport(PHG_BinParamsTy *binParams, PHG_BinFieldsTy *binFields);
void 		PhgBinCompPetDA(PHG_TrackingPhoton *bluePhotonPtr, PHG_TrackingPhoton *pinkPhotonPtr, 
				double *anglePtr, double *distancePtr);
void 		PhgBinCompSpectDA(PHG_TrackingPhoton *photonPtr, double *anglePtr,
				double *distancePtr);
void		PhgBinPETPhotons(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData, PHG_BinFieldsTy *binFields,
					PHG_Decay *decay,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons);
void		PhgBinSPECTPhotons(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData, PHG_BinFieldsTy *binFields,
					PHG_Decay *decay,
					PHG_TrackingPhoton *photons, LbUsFourByte numPhotons);
void		PhgBinTerminate(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData, PHG_BinFieldsTy *binFields);
Boolean		PhgBinOpenImage(PHG_BinParamsTy *binParams, PHG_BinFieldsTy *binFields,
				PhoHFileHdrKindTy hdrKind, char *imageName,
				FILE **imageFile, PhoHFileHdrTy *headerPtr,
				LbHdrHkTy *headerHkPtr, LbUsFourByte dataSize, void *dataPtr);

#undef LOCALE
#endif /* PHG_BIN_HDR */
