/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1992-1997 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		PhoHStat.h
*			Revision Number:	1.0
*			Date last revised:	July 10, 1997
*			Programmer:			Steven Vannoy
*			Date Originated:	19 August 1992
*
*			Module Overview:	Definitions for PhoHStat.c.
*
*			References:			'Photon History Stat Processes' PHG design.
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
*********************************************************************************/
#ifndef PHO_HIST_STAT_HDR
#define PHO_HIST_STAT_HDR
#ifdef PHOTON_HISTORY_STAT
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* TYPES */

/* PROTOTYPES */
void			PhoHStatIncDetectedPhotonStats(PHG_TrackingPhoton *trackingPhotonPtr);
void			PhoHStatIncTotEscPhotons(void);
void			PhoHStatIncTotInvalPhoLocations(void);
void			PhoHStatIncTotLowEnPhotons(void);
void			PhoHStatIncTotAbsorbedPhotons(PHG_TrackingPhoton *trackingPhotonPtr);
void			PhoHStatIncTotPrimOnlyScatter(void);
void			PhoHStatInitPhgStatistics(void);
void			PhoHStatUpdateForDetAttemptWeight(double hitPhoton_Weight);
void			PhoHStatUpdateForDetHitWeight(double hitPhoton_Weight);
void			PhoHStatUpdateRoulettedPhotons(LbUsTwoByte	numRouletted);	
void			PhoHStatUpdateSplitPhotons(LbUsTwoByte numSplits);
void			PhoHStatUpdateStartedPhotons(PHG_TrackingPhoton *trackingPhotonPtr);
Boolean			PhoHStatWrite(void);
LbUsFourByte	PhoHStat_GetTotInvalPhoLocations(void);
#undef LOCALE
#endif /* PHO_HIST_STAT_HDR */
