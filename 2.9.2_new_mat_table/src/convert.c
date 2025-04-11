/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1994-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		convert.c
*			Revision Number:	1.2
*			Date last revised:	3 June 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	October 18, 1994
*
*			Module Overview:	Performs a type conversion from one file to the
*								next. For example, converts a file of four byte
*								reals to two byte integers.
*
*			References:			None.
*
**********************************************************************************
*
*			Global functions defined:
*
*			Global macros defined:
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
*			Revision description:
*
*********************************************************************************/

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbMemory.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbError.h"
#include "LbInterface.h"
#include "LbConvert.h"

/* LOCAL CONSTANTS */

/* LOCAL TYPES */

/* LOCAL GLOBALS */

/* PROTOTYPES */
Boolean			Convert(int argc, char *argv[]);

/* FUNCTIONS */

/* MACROS */
#define Is_InType()			LbFgIsSet(runTimeOptions, LBFlag0)	/* Did user supply runtime option of input type? */
#define Is_OutType()			LbFgIsSet(runTimeOptions, LBFlag1)	/* Did user supply runtime option of output type? */
#define Is_HdrSize()			LbFgIsSet(runTimeOptions, LBFlag2)	/* Did user supply runtime option of header size? */
#define Is_CopyHdr()			LbFgIsSet(runTimeOptions, LBFlag3)	/* Did user supply runtime option of copying header? */
#define Is_NumPerLine()			LbFgIsSet(runTimeOptions, LBFlag4)	/* Did user supply runtime option of numbers per line? */
#define Is_Scientific()			LbFgIsSet(runTimeOptions, LBFlag5)	/* Did user supply runtime option of scientific notation? */


/**********************
*	Convert
*
*	Purpose:	Execute the program.
*
*	Result:	True unless an error occurs.
***********************/
Boolean Convert(int argc, char *argv[])
{
	Boolean				okay = false;		/* Process Loop */
	LbUsFourByte		localArgc;			/* Unsigned version of argc */
	char				inputPath[256];		/* Path name for the input file */
	char				outputPath[256];	/* Path to the output file */
	char				*hdrPtr=0;			/* Pointer to header buffer */
	FILE				*inFile=0;			/* Our data file */
	FILE				*outFile=0;			/* Converted file */
	LbCvEnDataType		inType;				/* Type of input data */
	LbCvEnDataType		outType;			/* Type of output data */
	LbUsFourByte		numPerLine;			/* Num per line for ascii conversion */
	LbUsFourByte		hdrSize;			/* Size of header to skip */
	LbUsFourByte		runTimeOptions=0;	/* The runtime options specified */
	
	/* The following variables are for getting run time options from
		the command line 
	*/
	#define	NUM_FLAGS	6
	
	char	*knownOptions[] = {"i:o:s:hn:e"};
	
	char				optArgs[NUM_FLAGS][LBEnMxArgLen];
	LbUsFourByte		optArgFlags = LBFlag0 + LBFlag1 + LBFlag2 + LBFlag3 + LBFlag4 + LBFlag5;
	LbUsFourByte		argIndex;


		
		
	do { /* Process Loop */
			
		/* Get our runtime options */
		if (!LbEnGetOptions(argc, argv, knownOptions,
				&runTimeOptions, optArgs, optArgFlags, &argIndex)) {
	
			break;
		}
		
		/* Get an unsigned version of argc */
		if (argc < 0) {
			localArgc = 0;
		}
		else {
			localArgc = argc;
		}
		
		/* Run through the options, verifying they are all there (that are necessary) */
		{
			/* Get the input data type */
			if (Is_InType()) {
				
				inType = (LbCvEnDataType) atoi(optArgs[0]);
			}
			else {
				ErStGeneric("You must specify an input type with the -i option");
				break;
			}
		
			/* Get the output data type */
			if (Is_OutType()) {
				
				outType = (LbCvEnDataType) atoi(optArgs[1]);
			}
			else {
				ErStGeneric("You must specify an output type with the -o option");
				break;
			}
		
			/* Get the hdr size to skip */
			if (Is_HdrSize()) {
				
				hdrSize = atoi(optArgs[2]);
			}
			else {
				hdrSize = 0;
			}
		
			/* Get the numbers per line for ASCII conversion */
			if (Is_NumPerLine()) {
				
				numPerLine = atoi(optArgs[4]);
			}
			else {
				numPerLine = 0;
			}
		}
	
		/* Get the name of the input file */
		if ((argIndex != 0) && (argIndex < localArgc) && (argv[argIndex] != 0)) {
			strcpy(inputPath, argv[argIndex]);
	
			/* Go to next argument */
			argIndex++;
		}
		else {
			ErStGeneric("You must supply an input file name");
			break;
		}
	
		/* Get the name of the output file */
		if ((argv[argIndex] != 0)&& (argIndex < localArgc) && (argv[argIndex] != 0)) {
			strcpy(outputPath, argv[argIndex]);
	
			/* Go to next argument */
			argIndex++;
		}
		else {
			ErAbort("You must supply an output file name");
		}
		
		/* Open data file file */
		if ((inFile = LbFlFileOpen(inputPath, "rb")) == 0) {
			ErStFileError("Unable to open input file.");
			break;
		}
		
		/* Turn off buffering */
		setbuf(inFile, 0);
		

		/* Open the output file */
		if ((outFile = LbFlFileOpen(outputPath, "wb")) == 0) {
			ErStFileError("Unable to create/open output file");
			break;
		}		
		
		
		/* Seek past header if requested */
		if ((hdrSize != 0) && (Is_CopyHdr() == false)) {

			if (fseek(inFile, hdrSize, SEEK_SET) != 0) {
				ErStFileError("\nUnable to seek past header.");
				break;
			}
		}
		else if ((hdrSize != 0) && (Is_CopyHdr() == true)) {
		
			/* Allocate a buffer to copy the header into */
			if ((hdrPtr = (char *) LbMmAlloc(hdrSize)) == 0) {
				break;
			}
			
			/* Read the input header */
			if (fread(hdrPtr, hdrSize, 1, inFile) != 1) {
				ErStFileError("Unable to read header from input file");
				break;
			}
			
			/* Write the header */
			if (fwrite(hdrPtr, hdrSize, 1, outFile) != 1) {
				ErStFileError("Unable to write header to output file");
				break;
			}
			
			
			/* Free our allocated memory */
			LbMmFree((void **)&hdrPtr);
		}

		/* Do the conversion */
		if (LbCvConvert(inFile, outFile, inType, outType, numPerLine, Is_Scientific())
				== false) {
			break;
		}

		/* Close the i/o files */
		fclose(inFile);
		fclose(outFile);
		inFile = outFile = 0;
		
		okay = true;
	} while (false);
	
	
	/* Do house cleaning that might have gotten skiped */
	{
		if (hdrPtr != 0)
			LbMmFree((void **)&hdrPtr);
			
		if (inFile != 0)
			fclose(inFile);
			
		if (outFile != 0)
			fclose(outFile);
	}
		

	/* Quit the program */
	return (okay);
}

