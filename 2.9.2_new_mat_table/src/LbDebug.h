/*********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992-2012 Department of Radiology	       		*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/*********************************************************************************
*
*     Module Name:        LbDebug.h
*     Revision Number:    1.1
*     Date last revised:  1 October 2012
*     Programmer:         Steven Vannoy
*     Date Originated:    Tuesday, July 21 81992
*
*     Module Overview:    This module declars all constants, types, variables
*						and prototypes for the Debug library.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*     Global variables defined:   
*
*********************************************************************************/

#include "LbTypes.h"


/* PROTOTYPES */
void	LbDbEnter(Boolean *continuePtr);
Boolean	LbDbInit(char *progName, LbUsFourByte initFlags, FILE *outputFile);
void	LbDbOffer(Boolean *continuePtr);
void	LbDbTerminate(void);
