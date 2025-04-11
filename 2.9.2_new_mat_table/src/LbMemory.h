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
*     Module Name:        LbMemory.h
*     Revision Number:    1.1
*     Date last revised:  1 October 2012
*     Programmer:         Steven Vannoy
*     Date Originated:    12 November 1992
*
*     Module Overview:    This module declars all constants, types, variables
*						and prototypes for the Memory library.
*
*     References:         
*
**********************************************
*
*		Modified:		24 Jan 2019
*		Programmer:		R Harrison
*		Purpose:		Added a function to allocate memory blocks with size
*						greater than can be allocated with a 4-byte size
*						(LbMmAlloc8).
*		 
**********************************************************************************
*
*     Global functions defined:
*
*     Global variables defined:   
*
*********************************************************************************/

#ifdef LB_MEMORY
	#define LOCALE
#else	
	#define LOCALE extern
#endif

#include "LbTypes.h"


#define LBMMFg_Accounting	LBFlag0			/* Perform memory accounting */


#ifdef LB_DEBUG
struct  lbMmMemBin {
	void				*thePtr;			/* The address of memory allocated */
	LbUsFourByte		byteCount;			/* Size, in bytes, of memory allocated */
	LbUsFourByte		allocNumber;		/* i'th call to LbMmAlloc */
	struct lbMmMemBin	*nextBinPtr;		/* Next bin in the list */
};

typedef struct lbMmMemBin		lbMmMemBinTy;
typedef lbMmMemBinTy	*lbMmMemBinPtr;

LOCALE lbMmMemBinPtr		lbMmTrackBin;			/* Used to track down memory problems */
LOCALE LbUsFourByte		lbMmTrackNum;			/* Used to track down memory problems */

#endif

/* PROTOTYPES */
void	LbMmCheckTrackNum(void);
void	LbMmCheckList(void);
void	*LbMmAlloc(LbUsFourByte bytesToAlloc);
void	*LbMmAlloc8(LbUsEightByte bytesToAlloc);
void	LbMmFree(void **memPtr);
Boolean	LbMmInit(LbUsFourByte flags);
void	LbMmTerminate(void);
#undef LOCALE
