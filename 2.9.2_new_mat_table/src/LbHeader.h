/*********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1995-2012 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/*********************************************************************************
*
*     Module Name:        LbHeader.h
*     Revision Number:    1.1
*     Date last revised:  1 October 2012
*     Programmer:         Steven Vannoy
*     Date Originated:     Friday, April 14, 1995
*
*     Module Overview:	This module provides routines for managing generic
*						"headers".
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*			LbHdrNew
*			LbHdrClear
*			LbHdrGtElem
*			LbHdrStElem
*			LbHdrRead
*			LbHdrWrite
*			LbHdrOpen
*
*     Global variables defined:   
*
*********************************************************************************/
#ifndef	LBHEADER
#define LBHEADER

#include "LbTypes.h"


/* CONSTANTS */

/* TYPES */
typedef	struct {
	FILE 			*fileRef;		/* The file hook */
	LbFourByte		headerSize;		/* The size of the header */
	void			*headerData;	/* The actual header data */	
} LbHdrHkTy;		/* A header hook */


/* PROTOTYPES */
void	LbHdrStNull(LbHdrHkTy *headerHk);
void	LbHdrFree(LbHdrHkTy *headerHk);
Boolean	LbHdrStFile(LbHdrHkTy *headerHk, FILE *headerFl);
Boolean	LbHdrNew(FILE *headerFile, LbFourByte size, LbHdrHkTy *headerHk);
Boolean	LbHdrClear(LbHdrHkTy *headerHk);
Boolean	LbHdrGtElem(LbHdrHkTy *headerHk, LbFourByte elemID, LbFourByte elemSize,
			void *elemData);
Boolean	LbHdrStElem(LbHdrHkTy *headerHk, LbFourByte elemID, LbFourByte elemSize,
			void *elemData);
Boolean	LbHdrRead(LbHdrHkTy *headerHk);
Boolean	LbHdrWrite(LbHdrHkTy *headerHk);
Boolean	LbHdrOpen(FILE *headerFile, LbFourByte size, LbHdrHkTy *headerHk);
#endif /* LBHEADER */
