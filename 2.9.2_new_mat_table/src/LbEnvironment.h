/*********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992-2012 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/*********************************************************************************
*
*     Module Name:        LbError.h
*     Revision Number:    1.1
*     Date last revised:  1 October 2012
*     Programmer:         Steven Vannoy
*     Date Originated:    Tuesday, July 21, 1992
*
*     Module Overview:	This file defines types, variables, macros, and functions
*						for the environment library.
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
#ifndef LIB_ENVIRONMENT
#define LIB_ENVIRONMENT

#include "LbTypes.h"


/* CONSTANTS */
#define		LBEnMxArgLen	32
#define		LBEnMxSwitchLen	32
/* PROTOTYPES */
Boolean	LbEnGetOptions(int argc, char **argv, char **optStr,
				LbUsFourByte *optFlags, char optArgs[][LBEnMxArgLen],
				LbUsFourByte optArgFlags,
				LbUsFourByte *firstArg);
				
#endif /* LIB_ENVIRONMENT */
