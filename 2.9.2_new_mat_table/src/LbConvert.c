/*********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1994-2013 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/*********************************************************************************
*
*     Module Name:        LbConvert.c
*     Revision Number:    1.4
*     Date last revised:  23 July 2013
*     Programmer:         Steven Vannoy
*     Date Originated:    Tuesday, October 18, 1994
*
*     Module Overview:    This module provides routines for converting binary
*							files from one type to another.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*     	LbCvConvert
*
*     Global variables defined:   
*
*********************************************************************************/

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbConvert.h"

/*	CONSTANTS */

/*	LOCAL TYPES */				
/*	LOCAL PROTOTYPES */
/*  LOCAL GLOBALS */
/*	LOCAL MACROS */

/*****************************************
*		LbCvConvert
*
*	Arguments:
*		FILE			*inFile			- The input file
*		FILE			*outFile 		- The output file
*		LbCvEnDataType	inType			- The "type" of the input data
*		LbCvEnDataType	outType			- The "type" of the output data
*		LbUsFourByte	numPerLine		- Numbers per line for ASCII conversion
*		Boolean			isScientific	- Scientific notation for ASCII conversion
*
*	Purpose:	Convert the data in the input file to the type requested
*				for the output file.
*
*				This function does not position the file pointer. Hence you
*				can call it with the input and output file's at any position.
*				The conversion takes place from the current point of inputFile
*				to the end-of-file for inputFile
*
*	Returns:	True unless an error occurs.
*
******************************************/
Boolean	LbCvConvert(FILE *inFile, FILE *outFile, LbCvEnDataType inType,
			LbCvEnDataType outType, LbUsFourByte numPerLine, Boolean isScientific)
{
	double			inBuff[8];			/* Input buffer, holds largest supported type double for allignment */
	char			outBuff[20];		/* Output buffer, holds long ASCII conversion */
	Boolean			okay = false;		/* Process Flag */
	Boolean			endOfFile =false;	/* End of file flag */
	char			errStr[1024];		/* Storage for an error string */
	LbUsFourByte	bytesRead = 0;		/* Number of bytes read */
	LbUsFourByte	inSize = 0;			/* Size of input data type */
	LbUsFourByte	outSize = 0;		/* Size of output data type */
	LbUsFourByte	curConv = 0;		/* Current conversion */

	do {
	
		/* To avoid user errors, clear numPerLine if outType != ASCII */
		if (outType != LbCvEn_ASCII)
			numPerLine = 0;
			
		
		/* Determine size of input value */
		switch (inType) {
			
			case LbCvEn_char:
				inSize = sizeof(char);
				break;
			
			case LbCvEn_uchar:
				inSize = sizeof(unsigned char);
				break;
			
			case LbCvEn_short:
				inSize = sizeof(short);
				break;
			
			case LbCvEn_ushort:
				inSize = sizeof(unsigned short);
				break;
			
			case LbCvEn_int:
				inSize = sizeof(int);
				break;
			
			case LbCvEn_uint:
				inSize = sizeof(unsigned int);
				break;
			
			case LbCvEn_long:
				inSize = sizeof(long);
				break;
			
			case LbCvEn_ulong:
				inSize = sizeof(unsigned long);
				break;
			
			case LbCvEn_float:
				inSize = sizeof(float);
				break;
			
			case LbCvEn_double:
				inSize = sizeof(double);
				break;
			
			case LbCvEn_LbOneByte:
				inSize = sizeof(LbOneByte);
				break;
			
			case LbCvEn_LbUsOneByte:
				inSize = sizeof(LbUsOneByte);
				break;
			
			case LbCvEn_LbTwoByte:
				inSize = sizeof(LbTwoByte);
				break;
			
			case LbCvEn_LbUsTwoByte:
				inSize = sizeof(LbUsTwoByte);
				break;
			
			case LbCvEn_LbFourByte:
				inSize = sizeof(LbFourByte);
				break;
			
			case LbCvEn_LbUsFourByte:
				inSize = sizeof(LbUsFourByte);
				break;
			
			case LbCvEn_ASCII:
				inSize = sizeof(char);
				break;
				
			default:
				ErStGeneric("Invalid value passed as input data type! (LbCvConvert)");
				goto FAIL;
		}
		
		/* Determine size of output value NOTE: this value is undetermined if it is 
			ASCII until the value is read and converted to a string
		 */
		switch (outType) {
			
			case LbCvEn_char:
				outSize = sizeof(char);
				break;
			
			case LbCvEn_uchar:
				outSize = sizeof(unsigned char);
				break;
			
			case LbCvEn_short:
				outSize = sizeof(short);
				break;
			
			case LbCvEn_ushort:
				outSize = sizeof(unsigned short);
				break;
			
			case LbCvEn_int:
				outSize = sizeof(int);
				break;
			
			case LbCvEn_uint:
				outSize = sizeof(unsigned int);
				break;
			
			case LbCvEn_long:
				outSize = sizeof(long);
				break;
			
			case LbCvEn_ulong:
				outSize = sizeof(unsigned long);
				break;
			
			case LbCvEn_float:
				outSize = sizeof(float);
				break;
			
			case LbCvEn_double:
				outSize = sizeof(double);
				break;
			
			case LbCvEn_LbOneByte:
				outSize = sizeof(LbOneByte);
				break;
			
			case LbCvEn_LbUsOneByte:
				outSize = sizeof(LbUsOneByte);
				break;
			
			case LbCvEn_LbTwoByte:
				outSize = sizeof(LbTwoByte);
				break;
			
			case LbCvEn_LbUsTwoByte:
				outSize = sizeof(LbUsTwoByte);
				break;
			
			case LbCvEn_LbFourByte:
				outSize = sizeof(LbFourByte);
				break;
			
			case LbCvEn_LbUsFourByte:
				outSize = sizeof(LbUsFourByte);
				break;
			
			case LbCvEn_ASCII:
				/* This size gets determined later */;
				break;
				
			default:
				ErStGeneric("Invalid value passed as output data type! (LbCvConvert)");
				goto FAIL;
		}
		
		/* Loop through the files */
		endOfFile = false;
		do {
			
			/* Read the input value and check for an error */
			if (fread(inBuff, inSize, 1, inFile) != 1) {
				
				/* See if we are at the end of file, if so we are done */
				if (feof(inFile) != 0) {
					endOfFile = true;
					break;
				}
				else {
					sprintf(errStr, "\nError on read, bytes read = %ld, attempt to read %ld.\n",
						(unsigned long)bytesRead, (unsigned long)inSize);
					ErStFileError(errStr);
					goto FAIL;
				}
			}
			
			/* Update number of bytes read */
			bytesRead += inSize;
			
			/* Perform the conversion */
			switch (outType) {
				
				case LbCvEn_char:
					switch (inType) {
						
						case LbCvEn_char:
							*((char *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((char *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((char *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((char *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((char *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((char *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((char *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((char *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((char *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((char *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((char *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((char *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((char *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((char *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((char *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((char *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_uchar:
					switch (inType) {
						
						case LbCvEn_char:
							*((unsigned char *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((unsigned char *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((unsigned char *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((unsigned char *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((unsigned char *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((unsigned char *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((unsigned char *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((unsigned char *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((unsigned char *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((unsigned char *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((unsigned char *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((unsigned char *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((unsigned char *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((unsigned char *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((unsigned char *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((unsigned char *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_short:
					switch (inType) {
						
						case LbCvEn_char:
							*((short *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((short *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((short *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((short *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((short *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((short *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((short *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((short *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((short *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((short *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((short *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((short *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((short *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((short *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((short *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((short *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_ushort:
					switch (inType) {
						
						case LbCvEn_char:
							*((unsigned short *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((unsigned short *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((unsigned short *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((unsigned short *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((unsigned short *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((unsigned short *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((unsigned short *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((unsigned short *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((unsigned short *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((unsigned short *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((unsigned short *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((unsigned short *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((unsigned short *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((unsigned short *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((unsigned short *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((unsigned short *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_int:
					switch (inType) {
						
						case LbCvEn_char:
							*((int *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((int *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((int *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((int *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((int *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((int *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((int *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((int *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((int *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((int *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((int *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((int *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((int *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((int *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((int *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((int *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_uint:
					switch (inType) {
						
						case LbCvEn_char:
							*((unsigned int *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((unsigned int *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((unsigned int *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((unsigned int *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((unsigned int *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((unsigned int *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((unsigned int *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((unsigned int *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((unsigned int *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((unsigned int *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((unsigned int *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((unsigned int *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((unsigned int *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((unsigned int *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((unsigned int *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((unsigned int *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_long:
					switch (inType) {
						
						case LbCvEn_char:
							*((long *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((long *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((long *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((long *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((long *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((long *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((long *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((long *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((long *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((long *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((long *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((long *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((long *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((long *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((long *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((long *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_ulong:
					switch (inType) {
						
						case LbCvEn_char:
							*((unsigned long *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((unsigned long *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((unsigned long *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((unsigned long *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((unsigned long *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((unsigned long *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((unsigned long *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((unsigned long *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((unsigned long *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((unsigned long *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((unsigned long *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((unsigned long *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((unsigned long *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((unsigned long *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((unsigned long *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((unsigned long *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_float:
					switch (inType) {
						
						case LbCvEn_char:
							*((float *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((float *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((float *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((float *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((float *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((float *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((float *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((float *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((float *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((float *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((float *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((float *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((float *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((float *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((float *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((float *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_double:
					switch (inType) {
						
						case LbCvEn_char:
							*((double *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((double *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((double *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((double *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((double *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((double *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((double *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((double *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((double *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((double *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((double *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((double *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((double *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((double *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((double *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((double *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_LbOneByte:
					switch (inType) {
						
						case LbCvEn_char:
							*((LbOneByte *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((LbOneByte *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((LbOneByte *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((LbOneByte *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((LbOneByte *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((LbOneByte *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((LbOneByte *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((LbOneByte *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((LbOneByte *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((LbOneByte *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((LbOneByte *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((LbOneByte *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((LbOneByte *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((LbOneByte *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((LbOneByte *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((LbOneByte *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				

				case LbCvEn_LbUsOneByte:
					switch (inType) {
						
						case LbCvEn_char:
							*((LbUsOneByte *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((LbUsOneByte *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((LbUsOneByte *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((LbUsOneByte *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((LbUsOneByte *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((LbUsOneByte *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((LbUsOneByte *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((LbUsOneByte *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((LbUsOneByte *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((LbUsOneByte *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((LbUsOneByte *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((LbUsOneByte *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((LbUsOneByte *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((LbUsOneByte *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((LbUsOneByte *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((LbUsOneByte *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_LbTwoByte:
					switch (inType) {
						
						case LbCvEn_char:
							*((LbTwoByte *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((LbTwoByte *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((LbTwoByte *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((LbTwoByte *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((LbTwoByte *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((LbTwoByte *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((LbTwoByte *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((LbTwoByte *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((LbTwoByte *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((LbTwoByte *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((LbTwoByte *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((LbTwoByte *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((LbTwoByte *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((LbTwoByte *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((LbTwoByte *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((LbTwoByte *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_LbUsTwoByte:
					switch (inType) {
						
						case LbCvEn_char:
							*((LbUsTwoByte *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((LbUsTwoByte *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((LbUsTwoByte *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((LbUsTwoByte *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((LbUsTwoByte *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((LbUsTwoByte *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((LbUsTwoByte *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((LbUsTwoByte *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((LbUsTwoByte *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((LbUsTwoByte *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((LbUsTwoByte *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((LbUsTwoByte *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((LbUsTwoByte *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((LbUsTwoByte *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((LbUsTwoByte *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((LbUsTwoByte *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_LbFourByte:
					switch (inType) {
						
						case LbCvEn_char:
							*((LbFourByte *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((LbFourByte *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((LbFourByte *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((LbFourByte *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((LbFourByte *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((LbFourByte *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((LbFourByte *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((LbFourByte *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((LbFourByte *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((LbFourByte *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((LbFourByte *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((LbFourByte *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((LbFourByte *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((LbFourByte *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((LbFourByte *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((LbFourByte *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
				
				case LbCvEn_LbUsFourByte:
					switch (inType) {
						
						case LbCvEn_char:
							*((LbUsFourByte *)outBuff) =  *((char *)inBuff);
							break;
						
						case LbCvEn_uchar:
							*((LbUsFourByte *)outBuff) =  *((unsigned char *)inBuff);
							break;
						
						case LbCvEn_short:
							*((LbUsFourByte *)outBuff) =   *((short *)inBuff);
							break;
						
						case LbCvEn_ushort:
							*((LbUsFourByte *)outBuff) =   *((unsigned short *)inBuff);
							break;
						
						case LbCvEn_int:
							*((LbUsFourByte *)outBuff) =   *((int *)inBuff);
							break;
						
						case LbCvEn_uint:
							*((LbUsFourByte *)outBuff) =   *((unsigned int *)inBuff);
							break;
						
						case LbCvEn_long:
							*((LbUsFourByte *)outBuff) =   *((long int *)inBuff);
							break;
						
						case LbCvEn_ulong:
							*((LbUsFourByte *)outBuff) =   *((unsigned long int *)inBuff);
							break;
						
						case LbCvEn_float:
							*((LbUsFourByte *)outBuff) =   *((float *)inBuff);
							break;
						
						case LbCvEn_double:
							*((LbUsFourByte *)outBuff) =   *((double *)inBuff);
							break;
						
						case LbCvEn_LbOneByte:
							*((LbUsFourByte *)outBuff) =  *((LbOneByte *)inBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							*((LbUsFourByte *)outBuff) =   *((LbUsOneByte *)inBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							*((LbUsFourByte *)outBuff) =   *((LbTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							*((LbUsFourByte *)outBuff) =   *((LbUsTwoByte *)inBuff);
							break;
						
						case LbCvEn_LbFourByte:
							*((LbUsFourByte *)outBuff) =  *((LbFourByte *)inBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							*((LbUsFourByte *)outBuff) =   *((LbUsFourByte *)inBuff);
							break;
						
						case LbCvEn_ASCII:
							/* Do nothing; it is treated as a string later as LbCvEn_ASCII */
							break;
					}
					break;
					
				case LbCvEn_ASCII:
					switch (inType) {
						
						case LbCvEn_ASCII:
							/* Convert new lines to whatever they aren't */
							if (*((char *)inBuff) == 13)
								*((char *)outBuff) = 10;
							else if (*((char *)inBuff) == 10)
								*((char *)outBuff) = 13;
							else
								*((char *)outBuff) = *((char *)inBuff);

							outSize = 1;
							break;
							
						case LbCvEn_char:
							sprintf(outBuff, "%d\t", *((char *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_uchar:
							sprintf(outBuff, "%d\t", *((unsigned char *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_short:
							sprintf(outBuff, "%d\t", *((short *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_ushort:
							sprintf(outBuff, "%hd\t", *((unsigned short *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_int:
							sprintf(outBuff, "%d\t", *((int *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_uint:
							sprintf(outBuff, "%d\t", *((unsigned int *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_long:
							sprintf(outBuff, "%ld\t", *((long int *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_ulong:
							sprintf(outBuff, "%ld\t", *((unsigned long int *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_float:
							if (isScientific) {
								sprintf(outBuff, "%3.3e\t", *((float *)inBuff));
							}
							else {
								sprintf(outBuff, "%f\t", *((float *)inBuff));
							}
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_double:
							if (isScientific) {
								sprintf(outBuff, "%3.3e\t", *((double *)inBuff));
							}
							else {
								sprintf(outBuff, "%f\t", *((double *)inBuff));
							}
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_LbOneByte:
							sprintf(outBuff, "%d\t", *((LbOneByte *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_LbUsOneByte:
							sprintf(outBuff, "%d\t", *((LbUsOneByte *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_LbTwoByte:
							sprintf(outBuff, "%d\t", *((LbTwoByte *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_LbUsTwoByte:
							sprintf(outBuff, "%hd\t", *((LbUsTwoByte *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_LbFourByte:
							sprintf(outBuff, "%ld\t", (long)*((LbFourByte *)inBuff));
							outSize = strlen(outBuff);
							break;
						
						case LbCvEn_LbUsFourByte:
							sprintf(outBuff, "%ld\t", (unsigned long)*((LbUsFourByte *)inBuff));
							outSize = strlen(outBuff);
							break;
					}
					break;
			}
			
			/* Increment our conversion count */
			curConv++;
			
			/* Write the output value */
			if (fwrite(outBuff, outSize, 1, outFile) != 1) {
				ErStFileError("Unable to write to output file (LbCvConvert)");
				goto FAIL;
			}
			
			/* Write a new line character if appropriate */
			if (numPerLine != 0) {
				if ((curConv % numPerLine) == 0) {
					if (fwrite("\n", 1, 1, outFile) != 1) {
						ErStFileError("Unable to write new line (LbCvConvert)");
						goto FAIL;
					}
				}
			}
			
		} while (endOfFile == false);
		
		okay = true;
		FAIL:;
	} while (false);
	
	return (okay);
}

/*****************************************
*		LbCvSwap
*
*	Arguments:
*		char			*srcBuf		- The source data
*		char			*dstBuf		- The destination data
*		LbUsFourByte	numBytes	- Number of bytes to swap
*
*	Purpose:	Swap the order of the bytes from srcBuf into dstBuf.
*
*				This function does not position the file pointer. Hence you
*
*	Returns:	None.
*
******************************************/
void	LbCvSwap(char *srcBuf, char *dstBuf, LbUsFourByte numBytes)
{
	LbUsFourByte byteIndex;
		
	/* Swap the bytes */
	for (byteIndex = 0; byteIndex < numBytes; byteIndex++) {
		dstBuf[(numBytes-1)-byteIndex] = srcBuf[byteIndex];
	}
}
