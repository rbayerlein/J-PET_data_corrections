/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1993-2012 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		PhgUsrBin.h
*			Revision Number:	2.0
*			Date last revised:	1 February 2012
*			Programmer:			Steven Vannoy
*			Date Originated:	28 July 1993
*
*			Module Overview:	Definitions for PhgUsrBin.c.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:		none
*
*			Global variables defined:		
*					BinUsrInitializeFPtr;
*					BinUsrModPETPhotonsFPtr;
*					BinUsrModPETPhotons2FPtr;
*					BinUsrModSPECTPhotonsFPtr;
*					BinUsrModSPECTPhotons2FPtr;
*					BinUsrTerminateFPtr;
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
*			Revision date:		1 February 2012
*
*			Revision description:	Changed form of user functions to pointers
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*			Revision description:
*						- binning by random-state and crystal number.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:	Added time-of-flight support.
*
*********************************************************************************/

#ifndef PHG_USR_BIN_HDR
#define PHG_USR_BIN_HDR

#ifdef PHG_USR_BIN
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */

/* DEFINITIONS */
typedef 	void BinUsrParamsFType(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData);
typedef 	Boolean BinUsrPETTrackingFType(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData,
											PHG_Decay			*decay,
											PHG_TrackingPhoton	*bluePhoton,
											PHG_TrackingPhoton	*pinkPhoton);
typedef 	Boolean BinUsrPETTrackingFType2(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData,
											PHG_Decay			*decay,
											PHG_TrackingPhoton	*bluePhoton,
											PHG_TrackingPhoton	*pinkPhoton,
											LbUsFourByte *angleIndex,
											LbUsFourByte *distIndex,
											LbUsFourByte *tofIndex,
											LbUsFourByte *scatter1Index,
											LbUsFourByte *scatter2Index,
											LbUsFourByte *crystal1Index,
											LbUsFourByte *crystal2Index,
											LbUsFourByte *energy1Index,
											LbUsFourByte *energy2Index,
											LbUsFourByte *zDownIndex,
											LbUsFourByte *zUpIndex,
											LbUsFourByte *thetaIndex,
											LbUsFourByte *phiIndex,
											LbUsFourByte *xrIndex,
											LbUsFourByte *yrIndex,
											LbUsEightByte *imageIndex,
											double		 *coincidenceWt,
											double		 *coincidenceSqWt);
typedef 	Boolean BinUsrSPECTTrackingFType(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData,
											PHG_Decay			*decay,
											PHG_TrackingPhoton	*photon);
typedef 	Boolean BinUsrSPECTTrackingFType2(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData,
											PHG_Decay			*decay,
											PHG_TrackingPhoton	*photon,
											LbUsFourByte *angleIndex,
											LbUsFourByte *distIndex,
											LbUsFourByte *scatterIndex,
											LbUsFourByte *energyIndex,
											LbUsFourByte *crystalIndex,
											LbUsFourByte *zIndex,
											LbUsEightByte *imageIndex);

/* GLOBALS */
extern BinUsrParamsFType			*BinUsrInitializeFPtr;
extern BinUsrPETTrackingFType		*BinUsrModPETPhotonsFPtr;
extern BinUsrPETTrackingFType2		*BinUsrModPETPhotonsF2Ptr;
extern BinUsrSPECTTrackingFType		*BinUsrModSPECTPhotonsFPtr;
extern BinUsrSPECTTrackingFType2	*BinUsrModSPECTPhotonsF2Ptr;
extern BinUsrParamsFType			*BinUsrTerminateFPtr;

/* PROTOTYPES */
			
#undef LOCALE
#endif /* PHG_USR_BIN_HDRS */
