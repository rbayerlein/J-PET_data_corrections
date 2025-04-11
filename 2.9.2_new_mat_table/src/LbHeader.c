/*********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1995-2013 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/*********************************************************************************
*
*     Module Name:        LbHeader.c
*     Revision Number:    1.2
*     Date last revised:  4 June 2013
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
*
*     Global variables defined:   
*
*********************************************************************************/

#include	<stdio.h>
#include	<string.h>
#include	<ctype.h>

#include "SystemDependent.h"

#include	"LbTypes.h"
#include	"LbError.h"
#include	"LbMemory.h"
#include	"LbHeader.h"

#ifdef _POSIX_SOURCE
	#include <unistd.h>
#endif

/*	LOCAL CONSTANTS */

/*  LOCAL GLOBALS */
char			lbHdErrStr[1024];	/* For creating detailed error messages */
				
/*	LOCAL MACROS */
/*	LOCAL FUNCTIONS */
			
/*	LOCAL MACROS */

/*	FUNCTIONS */
/*****************************************
*		LbHdrStNull
*
*	Purpose:	Sets a header hook to null,
*				equivilent to initializing
*				a variable to zero.
*	Arguments:
*		LbHdrHkTy	*headerHk	- The header hook.
*
*	Returns:	None.
*
******************************************/
void	LbHdrStNull(LbHdrHkTy *headerHk)
{
		/* Initialize the header hook */
		headerHk->fileRef = 0;
		headerHk->headerSize = -1;
		headerHk->headerData = 0;
}

/*****************************************
*		LbHdrFree
*
*	Purpose:	Free's the memory associated with a
*				header hook.
*	Arguments:
*		LbHdrHkTy	*headerHk	- The header hook.
*
*	Returns:	None.
*
******************************************/
void	LbHdrFree(LbHdrHkTy *headerHk)
{
		/* Initialize the header hook */
		headerHk->fileRef = 0;
		headerHk->headerSize = -1;
		if (headerHk->headerData != 0)
			LbMmFree((void **)&headerHk->headerData );
}

/*****************************************
*		LbHdrStFile
*
*	Purpose:	Associate the header with the
*				given file.
*	Arguments:
*		LbHdrHkTy	*headerHk	- The header hook.
*		FILE		*headerFl	- The new file
*
*	Returns:	TRUE unless an error occurs.
*
******************************************/
Boolean	LbHdrStFile(LbHdrHkTy *headerHk, FILE *headerFl)
{
	Boolean okay = false; /* Process flag */
	
	do {
		/* Initialize the header hook */
		headerHk->fileRef = headerFl;
		
		/* Write the header to the file */
		if (LbHdrWrite(headerHk) == false)
			break;
			
		okay = true;
	} while (false);
	
	return(okay);

}

/*****************************************
*		LbHdrNew
*
*	Purpose:	Create a new header associated
*				with the given file.
*	Arguments:
*		FILE		*headerFile	- The file for the header
*		LbFourByte	size		- The size of the header
*		LbHdrHkTy	*headerHk	- The new header hook.
*
*	Returns:	True unless an error occurrs.
*
******************************************/
Boolean	LbHdrNew(FILE *headerFile, LbFourByte size,
			LbHdrHkTy *headerHk)
{
	Boolean		okay = false;	/* Process Flag */
	signed char flag = -1;		/* "Empty" field flag */
	size_t		numWritten;		/* Number of bytes written */

	do { /* Process Loop */
	
		/* Initialize the header hook */
		headerHk->fileRef = headerFile;
		headerHk->headerSize = size;
		
		/* Allocate memory for empty header */
		if ((headerHk->headerData = LbMmAlloc(size)) == 0) {
			break;
		}
		
		/* Initialize the header to "empty" field value */
		memset(headerHk->headerData, flag, size);

		/* Set the current position within the file to the beginning */
		if (fseek(headerFile, 0, SEEK_SET) != 0) {
			ErStFileError("Unable to seek to beginning of header file.");
			break;
		}
		
		/* Write the empty field flag */
		if ((numWritten = fwrite(headerHk->headerData, size, 1,
				headerFile)) != 1){

			ErStFileError("Unable to write empty header to file.");
			break;
		}
		
		okay = true;
	} while (false);
	
	return (okay);
}

/*****************************************
*		LbHdrOpen
*
*	Purpose:	Open a new header associated
*				with the given file.
*	Arguments:
*		FILE		*headerFile	- The file for the header
*		LbFourByte	size		- The size of the header
*		LbHdrHkTy	*headerHk	- The new header hook.
*
*	Returns:	None.
*
******************************************/
Boolean	LbHdrOpen(FILE *headerFile, LbFourByte size,
			LbHdrHkTy *headerHk)
{
	Boolean		okay = false;	/* Process Flag */
	size_t		numRead;		/* Number of bytes read */

	do { /* Process Loop */
	
		/* Initialize the header hook */
		headerHk->fileRef = headerFile;
		headerHk->headerSize = size;
		
		/* Allocate memory for empty header */
		if ((headerHk->headerData = LbMmAlloc(size)) == 0) {
			break;
		}
		
		/* Initialize the header to "empty" field value */
		memset(headerHk->headerData, -1, size);

		/* Set the current position within the file to the beginning */
		if (fseek(headerFile, 0, SEEK_SET) != 0) {
			ErStFileError("Unable to seek to beginning of header file.");
			break;
		}
		
		/* Read the header */
		if ((numRead = fread(headerHk->headerData, size, 1,
				headerFile)) != 1){

			ErStFileError("Unable to read empty header from file.");
			break;
		}
		
		okay = true;
	} while (false);
	
	return (okay);
}

/*****************************************
*		LbHdrWrite
*
*	Purpose:	Write header associated
*				with the given hook.
*	Arguments:
*		LbHdrHkTy	*headerHk	- The header hook.
*
*	Returns:	None.
*
******************************************/
Boolean	LbHdrWrite(LbHdrHkTy *headerHk)
{
	Boolean		okay = false;	/* Process Flag */
	size_t		numWritten;		/* Number of bytes written */
	
	do { /* Process Loop */
		
		/* Set the current position within the file to the beginning */
		if (fseek(headerHk->fileRef, 0, SEEK_SET) != 0) {
			ErStFileError("Unable to seek to beginning of header file.");
			break;
		}
		
		/* Write the header */
		if ((numWritten = fwrite(headerHk->headerData, headerHk->headerSize, 1,
				headerHk->fileRef)) != 1){

			ErStFileError("Unable to write header to file.");
			break;
		}
		
		okay = true;
	} while (false);
	
	return (okay);
}

/*****************************************
*		LbHdrRead
*
*	Purpose:	Read header associated
*				with the given hook.
*	Arguments:
*		LbHdrHkTy	*headerHk	- The header hook.
*
*	Returns:	None.
*
******************************************/
Boolean	LbHdrRead(LbHdrHkTy *headerHk)
			
{
	Boolean		okay = false;	/* Process Flag */
	size_t		numRead;		/* Number of bytes read */
	long		curPos;			/* Current file position */
	
	do { /* Process Loop */
	
		/* Get the current file position */
		if ((curPos = ftell(headerHk->fileRef)) == -1L){
			ErStFileError("Unable to get current file position.");
			break;
		}
		
		/* Set the current position within the file to the beginning */
		if (fseek(headerHk->fileRef, 0, SEEK_SET) != 0) {
			ErStFileError("Unable to seek to beginning of header file.");
			break;
		}
		
		/* Read the header */
		if ((numRead = fread(headerHk->headerData, headerHk->headerSize, 1,
				headerHk->fileRef)) != 1){

			ErStFileError("Unable to read header from file.");
			break;
		}
		
		/* Seek back to original position */
		if (fseek(headerHk->fileRef, curPos, SEEK_SET) != 0) {
			ErStFileError("Unable to seek to original position in file (LbHdrRead).");
			break;
		}
		
		okay = true;
	} while (false);
	
	return (okay);
}

/*****************************************
*		LbHdrGtElem
*
*	Purpose:	Return the data associated
*				with the specified element
*				ID.
*	Arguments:
*		LbHdrHkTy	*headerHk	- The header hook.
*		LbFourByte	elemID		- The element's ID.
*		LbFourByte	elemSize	- The element's size.
*		void		*elemData	- Buffer for element's data.
*
*	Returns:	None.
*
******************************************/
Boolean	LbHdrGtElem(LbHdrHkTy *headerHk, LbFourByte elemID,
			LbFourByte elemSize, void *elemData)
			
{
	Boolean		okay = false;	/* Process Flag */
	LbFourByte	curID;			/* Current element ID */
	LbFourByte	curOffset;		/* Current element offset */
	LbFourByte	curSize;		/* Current element size */
	unsigned char *bytePtr;		/* Buffer pointer for search */
	do { /* Process Loop */
	
		/* Initialize items for search */
		bytePtr = (unsigned char *) headerHk->headerData;
		curOffset = 0;
		
		/* Search for element ID */
		do {
			/* Get the ID, this method prevents conflicts with byte alignment
				and pointer issues since some machines would require four
				byte values to be alligned on even boundaries, our method
				does not require this.
			*/
			curID = (bytePtr[curOffset] << 24);
			curID |= (bytePtr[curOffset+1] << 16);
			curID |= (bytePtr[curOffset+2] << 8);
			curID |= bytePtr[curOffset+3];
	
			/* If we have an element ID, get it's size */
			if (curID != -1) {
				/* Increment to element size */
				curOffset += sizeof(LbFourByte);
				
				/* Get the size of the element */
				curSize = (bytePtr[curOffset] << 24) | (bytePtr[curOffset+1] << 16)
					| (bytePtr[curOffset+2] << 8) | bytePtr[curOffset+3];


				/* Increment the offset */
				curOffset += sizeof(LbFourByte) + curSize;
			}
		} while ((curID != -1) && (curID != elemID));
		
		/* If curID == -1, element not found */
		if (curID == -1) {
			ErStPrimError(ERMgCdHeader, ERErCdHdElemNotFound,
				"Requested element ID not found.");
			break;
		}

		/* Verify that the size is write */
		if (curSize != elemSize) {
			ErStGeneric("Requested element size does not match element size in header.");
			break;
		}
		
		/* Since we found it, we have already indexed past the data so we
			must backup
		 */
		curOffset -= curSize;
		
		/* We made it here so copy the data */
		bytePtr = (unsigned char *)headerHk->headerData + curOffset;
		
		memcpy(elemData, bytePtr, elemSize);
		
		okay = true;
	} while (false);
	
	return (okay);
}

/*****************************************
*		LbHdrStElem
*
*	Purpose:	Set the data associated
*				with the specified element
*				ID.
*	Arguments:
*		LbHdrHkTy	*headerHk	- The header hook.
*		LbFourByte	elemID		- The element's ID.
*		LbFourByte	elemSize	- The element's size.
*		void		*elemData	- Buffer for element's data.
*
*	Returns:	None.
*
******************************************/
Boolean	LbHdrStElem(LbHdrHkTy *headerHk, LbFourByte elemID,
			LbFourByte elemSize, void *elemData)
			
{
	Boolean		okay = false;	/* Process Flag */
	LbFourByte	curID;			/* Current element ID */
	LbFourByte	curOffset;		/* Current element offset */
	LbFourByte	curSize;		/* Current element size */
	unsigned char *bytePtr;		/* Buffer pointer for search */
	LbFourByte	neededSize;		/* Size needed to add field to header */
	
	do { /* Process Loop */
	
		/* Initialize items for search */
		bytePtr = (unsigned char *) headerHk->headerData;
		
		#ifdef LB_DEBUG
			if (bytePtr == 0){
				ErAbort("Called LbHdrStElem with null header data");
			}
		#endif
		curOffset = 0;
		
		/* Search for element ID */
		do {
			/* Get the ID, this method prevents conflicts with byte alignment
				and pointer issues since some machines would require four
				byte values to be alligned on even boundaries, our method
				does not require this.
			*/
			curID = (bytePtr[curOffset] << 24);
			curID |= (bytePtr[curOffset+1] << 16);
			curID |= (bytePtr[curOffset+2] << 8);
			curID |= bytePtr[curOffset+3];
	
			/* If we have an element ID, get it's size */
			if (curID != -1) {
				/* Increment to element size */
				curOffset += sizeof(LbFourByte);
				
				/* Get the size of the element */
				curSize = (bytePtr[curOffset] << 24) | (bytePtr[curOffset+1] << 16)
					| (bytePtr[curOffset+2] << 8) | bytePtr[curOffset+3];

				/* Increment the offset */
				curOffset += sizeof(LbFourByte) + curSize;
			}
		} while ((curID != -1) && (curID != elemID));
		
		/* If curID == -1, element not found, verify there is room and
			set ID/size fields
		 */
		if (curID == -1) {
		
			/* If there isn't room flag it and bail */
			neededSize = curOffset+(sizeof(LbFourByte)*2)+elemSize;
			if (neededSize >= headerHk->headerSize){
				ErStGeneric("Not enough room in header to add new field.");
				break;
			}
			
			/* Set the ID */
			bytePtr[curOffset] = ((elemID & 0xFF000000) >> 24);
			bytePtr[curOffset+1] = ((elemID & 0x00FF0000) >> 16);
			bytePtr[curOffset+2] = ((elemID & 0x0000FF00) >> 8);
			bytePtr[curOffset+3] = (elemID & 0x000000FF);
			
			/* Increment to size field */
			curOffset += sizeof(LbFourByte);
			
			/* Set the element size */
			bytePtr[curOffset] = ((elemSize & 0xFF000000) >> 24);
			bytePtr[curOffset+1] = ((elemSize & 0x00FF0000) >> 16);
			bytePtr[curOffset+2] = ((elemSize & 0x0000FF00) >> 8);
			bytePtr[curOffset+3] = (elemSize & 0x000000FF);
			
			/* Increment to data field */
			curOffset += sizeof(LbFourByte);
			bytePtr = (unsigned char *)headerHk->headerData + curOffset;
		}
		else {
			
			/* Verify sizes match */
			if (curSize != elemSize) {
				ErStGeneric("Element size does not match ID size in header.");
				break;
			}
			
			/* Since we found it, we have already indexed past the data so we
				must backup
			 */
			curOffset -= curSize;
			bytePtr = (unsigned char *)headerHk->headerData + curOffset;
		}
		
		/* We made it here so copy the data */
		memcpy(bytePtr, elemData, elemSize);
		
		okay = true;
	} while (false);
	
	return (okay);
}

