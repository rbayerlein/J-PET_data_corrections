/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1995-2001 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name: 			calc.attenuation.h
*     Revision Number:		1.0
*     Date last revised:	May 25, 1995
*     Programmer:			Steven Vannoy
*     Date Originated:		May 25, 1995
*
*     Module Overview:	This is the global include file for the calc.attenuation
*						utility.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*	  Global variables defined:   
*
*	  Global macros defined:
*
*********************************************************************************/
#ifndef CALC_ATTN_HDR
#define CALC_ATTN_HDR


#ifdef PHG_MAIN
	#define	LOCALE
#else
	#define LOCALE	extern
#endif


/* PROGRAM CONSTANTS */

/* OPTION MACROS */
/* DEBUGGING OPTIONS */
/* PROGRAM TYPES */

/* PROGRAM GLOBALS */

/* PROTOTYPES */
Boolean ClcAtnInitialize(int argc, char *argv[]);
Boolean	ClcAtnCalcAttenuation(void **atnAry);

#undef LOCALE
#endif /* CALC_ATTN_HDR */
