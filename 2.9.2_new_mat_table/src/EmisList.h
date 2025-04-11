/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992 Department of Radiology						*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			EmisList.h
*     Revision Number:		1.0
*     Date last revised:	Mondayb November 29, 1993
*     Programmer:			Steven Vannoy
*     Date Originated:		Friday September 18, 1992
*
*     Module Overview:	Definitions for Emission LIst Processes.
*
*     References:       'Emission List Gen Processes' PHG design.  
*
**********************************************************************************
*
*     Global functions defined:
*
*     Global variables defined:   
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
*						- support for eight-byte number of decays
*
*********************************************************************************/
#ifdef EMIS_LIST
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

#include <stdio.h>
#include <math.h>

#include "LbTypes.h"
#include "LbMacros.h"

#include "Photon.h"


/* CONSTANTS */
/* TYPES */

/* GLOBALS */
LOCALE	PHG_Decay			EmisListNewDecay;	   		/* The new decay */
LOCALE	PHG_TrackingPhoton	EmisListBluePhoton;			/* The "current" blue photon */
LOCALE	PHG_TrackingPhoton	EmisListPinkPhoton;			/* The "current" pink photon */
LOCALE	PHG_TrackingPhoton	EmisListFDPhoton;			/* Photon used for forced detection */
LOCALE	LbEightByte			EmisListNumCreated;			/* Number of photons created for tracking */
LOCALE	char				EmisListIsotopeDataFilePath[256];	/* Path to isotope data file */

#ifdef PHG_DEBUG
LOCALE	    LbEightByte		EmisListBreakNumber;		/* Used for breakpoints */
LOCALE		LbUsFourByte	EmisListRestartSlice;		/* Slice number necessary for restarting */
LOCALE		LbUsFourByte	EmisListRestartAngle;		/* Angle number necessary for restarting */
LOCALE		LbUsFourByte	EmisListRestartVoxel;		/* Voxel number necessary for restarting */
	
#endif

/* PROTOTYPES */
void 	EmisListCreateDetectedPhoton(PHG_TrackingPhoton	*trackingPhotonPtr);
Boolean EmisListCreatePhotonList(void);	
void	EmisListDoDetection(PHG_TrackingPhoton	*trackingPhotonPtr);
Boolean EmisListInitialize(void);
void	EmisListTerminate(void);
void	EmisLisTrackPhoton(PHG_TrackingPhoton *trackingPhotonPtr);
void	EmisListDoComptonInteraction(PHG_TrackingPhoton *trackingPhotonPtr);
void 	EmisListDoCoherent(PHG_TrackingPhoton	*trackingPhotonPtr, LbUsFourByte materialIndex);
#undef LOCALE
