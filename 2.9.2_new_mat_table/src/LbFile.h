/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 2003-2012 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/

/*********************************************************************************
*
*     Module Name:        LbFile.h
*     Revision Number:    1.5
*     Date last revised:  15 October 2012
*     Programmer:         Steven Gillispie
*     Date Originated:    7 May 2003
*
*     Module Overview:    This module declares all constants, types, variables
*							and prototypes for the File library.
*
*     References:         None
*
**********************************************************************************
*
*     Global functions defined:		
*			LbFlSetDir
*			LbFlFileOpen
*			LbFlFGetS
*
*     Global variables defined:   	None
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

#ifndef LIB_FILE
#define LIB_FILE

#include "LbTypes.h"


/* CONSTANTS */

/* FLAGS */

/* TYPES */

/* PROTOTYPES */
#ifdef MACGUIOS
	LbTwoByte LbFlSetDir(void);
#endif
FILE* LbFlFileOpen(char *path,  char *mode);
char* LbFlFGetS(char *str, int size, FILE *stream);

#endif /* LIB_FILE */
