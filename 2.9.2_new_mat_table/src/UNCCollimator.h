/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1994-2011 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/

/*********************************************************************************
*
*			Module Name:		UNCCollimator.h
*			Revision Number:	1.1
*			Date last revised:	8 December 2011
*			Programmer:			Steven Vannoy
*			Date Originated:	11 July 1994
*
*			Module Overview:	Definitions for UNCCollimator.c.
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
*********************************************************************************/
#ifndef UNC_COL_HDR
#define UNC_COL_HDR

#ifdef UNC_COLLIMATOR
	#define	LOCALE
#else
	#define LOCALE	extern
#endif


/* CONSTANTS */

/* TYPES */
/* GLOBALS */

/* MACROS */		

/*********************************************************************************
*
*			Name:		UNCCOLIsConeBeam
*
*			Summary:	Return true if collimator hole type is cone beam.
*			Arguments:
*				
*			Function return: true or false.
*
*********************************************************************************/
#define UNCCOLIsConeBeam() (ColRunTimeParams[ColCurParams].UNCSPECTCol.HoleGeometry == CONE)

/*********************************************************************************
*
*			Name:		UNCCOLGetFocalLength
*
*			Summary:	Return the collimator focal length.
*			Arguments:
*				
*			Function return: true or false.
*
*********************************************************************************/
#define UNCCOLGetFocalLength() (ColRunTimeParams[ColCurParams].UNCSPECTCol.FocalLength)

/*********************************************************************************
*			Name:		UNCCOLGetCollZMin
*
*			Summary:	Return the collimator focal length.
*			Arguments:
*				
*			Function return: true or false.
*
*********************************************************************************/
#define UNCCOLGetCollZMin() (ColRunTimeParams[ColCurParams].UNCSPECTCol.MinZ)

/*********************************************************************************
*			Name:		UNCCOLGetCollZMax
*
*			Summary:	Return the collimator focal length.
*			Arguments:
*				
*			Function return: true or false.
*
*********************************************************************************/
#define UNCCOLGetCollZMax() (ColRunTimeParams[ColCurParams].UNCSPECTCol.MaxZ)

/* PROTOTYPES */
Boolean			UNCColInitialize(void);
void			UNCColPrintParams(void);
void			UNCColPrintReport(void);
void			UNCCollimate(PHG_Decay *decayPtr, PHG_TrackingPhoton *photons,
					LbUsFourByte numPhotons, CollimatedPhotonsTy *colPhotonsPtr);
void			UNCColTerminate(void);
#undef LOCALE
#endif /* UNC_COL_HDR */
