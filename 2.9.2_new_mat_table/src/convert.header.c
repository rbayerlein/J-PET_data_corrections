/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1994-2001 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		convert.header.c
*			Revision Number:	1.0
*			Date last revised:	21 April, 1995
*			Programmer:			Steven Vannoy
*			Date Originated:	21 April, 1995
*
*			Module Overview:	Converts a file with 
*								an "old" style header to a new style one.
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
#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbMemory.h"
#include "LbInterface.h"
#include "LbParamFile.h"
#include "LbConvert.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhoHFile.h"
#include "PhgHdr.h"
#include "phg.h"
#include "PhgBin.h"

/* LOCAL CONSTANTS */

/* This dictates the size of the units we read from disk.
It will be applied to a "generic" buffer which then gets
interpreted based on the type of the binned data. Hence,
it must remain a power of 2 >= 8.
*/
#define BUFF_SIZE 2048

/* LOCAL TYPES */

/* LOCAL GLOBALS */
static	Boolean	canceled;	/* Global cancelation flag */
static	PhoHFileHdrTy		oldHdr;					/* Our old header */
static	LbHdrHkTy			newHdr;					/* The new header */			

/* PROTOTYPES */
Boolean			ConvertHeader(int argc, char *argv[]);

/* FUNCTIONS */

/**********************
*	ConvertHeader
*
*	Purpose:	Conver this file to have a new
*				style header.
*
*	Result:	True unless an error occurs.
***********************/
Boolean ConvertHeader(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */
	char				errStr[1024];			/* Error string buffer */
	FILE				*oldFile;				/* The output file */
	FILE				*newFile;				/* The output file */
	LbUsFourByte		numRead;				/* Number of bytes read from final file */
	void				*oldData;				/* Buffer for data portion of file */

	do { /* Process Loop */
			
		/* Verify command line contains two arguments */
		if (argc < 3) {
			/* Ask for the file name */
			ErAbort("\nThis program requires two arguments, an input file and an output file.\n");
		}
		
		/* Open the first file, the one that has the old style header */
		if ((oldFile = LbFlFileOpen(argv[1], "rb")) == 0) {
			sprintf(errStr, "Unable to open input file\n'%s'.", argv[1]);
			ErStFileError(errStr);
			goto FAIL;
		}
		
		/* Open the second file, the one that will have the old style header */
		if ((newFile = LbFlFileOpen(argv[2], "w+b")) == 0) {
			sprintf(errStr, "Unable to open output file\n'%s'.", argv[2]);
			ErStFileError(errStr);
			goto FAIL;
		}
		
		/* Turn buffering off, it shouldn't be necessary but bugs have been found in the
			DG i/o library that make it so.
		*/
		setbuf(oldFile, 0);
		setbuf(newFile, 0);
		
		/* Read the header from the final file */
		if (fread(&oldHdr, sizeof(PhoHFileHdrTy), 1, oldFile) != 1) {
			sprintf(errStr, "\nUnable to read header size from input file '%s'.\n",
				argv[1]);
			ErStFileError(errStr);
			goto FAIL;
		}
		
		/* Verify header is from the old style  version */
		if (oldHdr.H.HdrVersion != PHOHFILE_OLD_HEADER_VERSION) {
			sprintf(errStr, "File '%s' has a header version of %3.2f\n"
				"the current version is %3.2f.\n"
				"You will need a different version of convert.header to perform this operation.\n",
				argv[1], oldHdr.H.HdrVersion, PHOHFILE_OLD_HEADER_VERSION);
			ErStGeneric(errStr);
			goto FAIL;
		}
		
		/* Change the header version */
		oldHdr.H.HdrVersion = PHG_HDR_HEADER_VERSION;
		oldHdr.H.HdrSize = PHG_HDR_HEADER_SIZE;
		
		/* Create the new header and write to the file */
		if (!PhgHdrMkHeader(newFile, &oldHdr, &newHdr)) {
			goto FAIL;
		}
		
		/* Allocate memory buffers */
		if ((oldData = LbMmAlloc(BUFF_SIZE)) == 0) {
			goto FAIL;
		}
			
			
		/* Loop by reading a buffer from the final file and processing until there
			isn't a full buffer left (the partial buffer is handled next)
		 */
		while ((numRead = fread(oldData, 1, BUFF_SIZE, oldFile)) == BUFF_SIZE) {
						
			
			/* Write the data to the output file */
			if (fwrite(oldData, BUFF_SIZE, 1, newFile) != 1) {
				ErStFileError("\nUnable to write updated data to final file.");
				goto FAIL;
			}
		}
		
		/* Process partial buffer if there is  one */
		if (numRead != 0) {
			
			/* Write the data to the output file */
			if (fwrite(oldData, numRead, 1, newFile) != 1) {
				ErStFileError("\nUnable to write updated data to final file.");
				goto FAIL;
			}
		}
		/* If we are here and haven't read to EOF there is an error */
		else if (feof(oldFile) == 0) {
			sprintf(errStr, "Unable to read from file '%s'.\n", argv[1]);
			ErStFileError(errStr);
			goto FAIL;
		}
		
	
		/* Close the files */
		fclose(oldFile);
		fclose(newFile);
		
		okay = true;
		FAIL:;
	} while (false);
	
	/* If error set due to cancellation, handle it, otherwise pass it on */
	if (!okay) {
		if (canceled) {
			ErHandle("User canceled convert.header.", false);
			okay = true;
		}
	}
	
	/* Return the status */
	return (okay);
}
