/*********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992 - 2012 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/*********************************************************************************
*
*     Module Name:        LbParamFile.h
*     Revision Number:    1.2
*     Date last revised:  1 October 2012
*     Programmer:         Steven Vannoy
*     Date Originated:    Monday, March 1, 1993
*
*     Module Overview:    This module declars all constants, types, variables
*						and prototypes for the Parameter File library.
*
*     References:         
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
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		18 June 2012
*
*			Revision description:	Added LbPfGetNextToken
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*			Revision description:
*						- support for randoms and eight-byte number of decays
*
*********************************************************************************/
#ifndef	LBPF_PARAM
#define LBPF_PARAM

#include "LbTypes.h"


/* CONSTANTS */
#define LBPF_PARAM_LEN		256
#define LBPF_LABEL_LEN		64

/* TYPES */
typedef	struct {
	
	LbUsFourByte	curLine;		/* The current line number */
	FILE 			*fileRef;		/* The file hook */
	char			name[255];		/* Name of the parameter file */
} LbPfHkTy;		/* A parameter file hook */

/* NOTE: These types are position dependant, do not change them */
typedef enum {
	LbPfEnBoolean,
	LbPfEnComment,
	LbPfEnChar,
	LbPfEnString,
	LbPfEnInteger,
	LbPfEnReal,
	LbPfEnList,
	LbPfEnEnum,
	LbPfEnLongLong }
	LbPfEnPfTy;
	

/* PROTOTYPES */
Boolean	LbPfOpen(char *filePathPtr, LbUsFourByte flags,
			LbPfHkTy *pfFileHkPtr);
Boolean	LbPfGetParam(LbPfHkTy *pfFileHk, void *paramBufferPtr,
			LbPfEnPfTy *paramTypePtr, LbUsFourByte *paramSizePtr,
			char *labelPtr, Boolean *isEOF);
void	LbPfClose(LbPfHkTy *pfFileHk);
Boolean	LbPfGetNextToken( 	char 			*paramBufferPtr,
							char 			*tokenPtr, 
							LbUsFourByte 	*currentCharIndexPtr,
							LbUsFourByte 	*tokenLenPtr,
							Boolean 		*isLastTokenPtr);
#endif /* LBPF_PARAM */
