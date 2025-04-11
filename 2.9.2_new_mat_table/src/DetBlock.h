/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 2002-2009 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/

/*********************************************************************************
*
*			Module Name:		DetBlock.h
*			Revision Number:	2.14
*			Date last revised:	13 May 2009
*			Programmer:			Steven Gillispie
*			Date Originated:	9 August 2002
*
*			Module Overview:	Definitions for DetBlock.c
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:	
*
*			Global variables defined:	None
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

#ifndef DETECTOR_BLOCK_HDR
#define DETECTOR_BLOCK_HDR

#ifdef DETECTOR_BLOCK
	#define	LOCALE
#else
	#define LOCALE	extern
#endif


/* CONSTANTS */

/* TYPES */

/* GLOBALS */

/* PROTOTYPES */
Boolean 		DetBlocValidateBlocks(void);
LbUsFourByte	DetBlocDivideZones(LbUsFourByte		maxZBlocks);

double			DetBlocGtTruncatedFreePaths(PHG_TrackingPhoton *photonPtr, 
					LbUsFourByte ringNum, LbUsFourByte detNum, double *weight);
Boolean 		DetBlocProjectToDetector(PHG_TrackingPhoton *photonPtr, 
					LbFourByte *ringNum, LbFourByte *detNum);
void			DetBlocFindNextInteraction(PHG_TrackingPhoton *photonPtr, 
					LbFourByte *curRingNum, LbFourByte *curDetNum, 
					detEn_ActionTy *actionType, double *fpToGo, 
					LbUsFourByte *detMaterial, Boolean *isActive);
void			DetBlocFindDetPosition(
					PHG_TrackingPhoton *photonPtr, LbUsFourByte interactionIndex);
void			DetBlocFreeData(void);

#undef LOCALE
#endif /* DETECTOR_BLOCK_HDR */
