/*********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1994-2012 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/*********************************************************************************
*
*     Module Name:        LbConvert.h
*     Revision Number:    1.1
*     Date last revised:  1 October 2012
*     Programmer:         Steven Vannoy
*     Date Originated:    Tuesday, October 18, 1994
*
*     Module Overview:    This module declars all constants, types, variables
*						and prototypes for the Conversion library.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*     	LbCvCloseFile
*     	LbCvFilterLbCvFile
*     	LbCvGetHeader
*     	LbCvGetClump
*     	LbCvOpenFile
*     	LbCvProcessLbCvFile
*     	LbCvSetClump
*     	LbCvSetHeader
*
*     Global variables defined:   
*
*********************************************************************************/

#include "LbTypes.h"


/* FLAGS */
/* CONSTANTS */
/* TYPES */
typedef	enum {	LbCvEn_char,			/* 0  */
				LbCvEn_uchar,			/* 1  */
				LbCvEn_short,			/* 2  */
				LbCvEn_ushort,			/* 3  */
				LbCvEn_int,				/* 4  */
				LbCvEn_uint,			/* 5  */
				LbCvEn_long,			/* 6  */
				LbCvEn_ulong,			/* 7  */
				LbCvEn_float,			/* 8  */
				LbCvEn_double,			/* 9  */
				LbCvEn_LbOneByte,		/* 10 */
				LbCvEn_LbUsOneByte,		/* 11 */
				LbCvEn_LbTwoByte,		/* 12 */
				LbCvEn_LbUsTwoByte,		/* 13 */
				LbCvEn_LbFourByte,		/* 14 */
				LbCvEn_LbUsFourByte,	/* 15 */
				LbCvEn_ASCII			/* 16 */
			} LbCvEnDataType;
			
			
/* MACROS */

/* PROTOTYPES */
Boolean	LbCvConvert(FILE *inFile, FILE *outFile, LbCvEnDataType inType,
			LbCvEnDataType outType, LbUsFourByte numPerLine,
			Boolean isScientific);
			
void	LbCvSwap(char *srcBuf, char *dstBuf, LbUsFourByte numBytes);
