/*********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992-2013 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/*********************************************************************************
*
*     Module Name:        LbMemory.c
*     Revision Number:    1.5
*     Date last revised:  23 July 2013
*     Programmer:         Steven Vannoy
*     Date Originated:    12 November, 1992
*
*     Module Overview:	This module provides routines for memory management.						
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
*     	LbMmAlloc
*		LbMmFree
*		LbMmInit
*		LbMmTerminate
*
*     Global variables defined:   
*
*********************************************************************************/
#define LB_MEMORY

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbDebug.h"
#include "LbMemory.h"
#include "LbInterface.h"

/*	CONSTANTS */

/*  LOCAL GLOBALS */
char	lbMmErrStr[1024];			/* Error string storage */

#undef LB_DEBUG
#ifdef LB_DEBUG

#define	MAX_CONC_ALLOCATION		3000		/* Maximum supported concurrent allocations */

#define	LBMMCHECKSUM_SIZE			8

	
	lbMmMemBinPtr			lbMmMemList = 0;		/* A list of the memory blocks allocated */
	lbMmMemBinPtr			lbMmNewMemBin = 0;		/* Used to allocate new memory bins */
	Boolean					lbMmIsInited = false;	/* Our initialization flag */
	LbUsFourByte			lbMmInitFlags;			/* Flags we were initialized with */
	LbUsFourByte			lbMmBytesAllocated = 0;		/* Count of bytes allocated */
	LbUsFourByte			lbMmCurBytesAllocated = 0;	/* Current count of bytes allocated */
	LbUsFourByte			lbMmTotBytesAllocated = 0;	/* Total count of bytes allocated */
	LbFourByte				lbMmCurNumAllocs = 0;		/* Current number of allocated blocks */
	LbUsFourByte			lbMmTotNumAllocs = 0;		/*	Count of calls to allocate memory
														Not decremented on frees, it is a
														statistic on memory usage
													*/
	LbUsFourByte			lbMmMaxResident = 0;	/* Maximum amount of memory allocated
														at any given time.
													*/
LbUsTwoByte				lbMmAllocNumList[MAX_CONC_ALLOCATION];	/* This list contains
																the allocation numbers
																of the current memory list
																it is used to track down
																corruption of the memory list.
															
																Note this list is in the opposite
																order the allcation numbers will
																be found in the linked-list memory
																list for performance sake.
															*/
															
lbMmMemBinTy	dummyBin;

	char				lbMmCheckSum[] = {'c', 'h', 'e', 'c', 'k', 's', 'u', 'm'};
	
#endif

/*	LOCAL MACROS */

#ifdef LB_DEBUG
/*********************************************************************************
*		lbMmIsAccounting
*
*	Arguments:
*			LbUsFourByte	flags	- Were we initialized with accounting on.
*
*	Purpose:	Check to see if we are performing memory accounting.
*
*	Returns:	True if accounting is turned on.
*
*********************************************************************************/
#define lbMmIsAccounting(flags) LbFgIsSet(flags, LBMMFg_Accounting)
#endif

/*	LOCAL FUNCTIONS */
/*	LOCAL MACROS */

/*	FUNCTIONS */
#ifdef LB_DEBUG
/*********************************************************************************
*		LbMmCheckTrackNum
*
*	Arguments:
*
*	Purpose:	Verifies that the track bin and track num match.
*
*
*********************************************************************************/
void	LbMmCheckTrackNum()
{
	LbTwoByte		dummy = 0;
	char			*checkSumPtr;
	
	/* Verify the check sum is valid */
	checkSumPtr = (char *)lbMmTrackBin->thePtr;
	if (memcmp(&checkSumPtr[lbMmTrackBin->byteCount], lbMmCheckSum, LBMMCHECKSUM_SIZE) != 0) {
		LbInPrintf("\nFound memory checksum corrupted\n");
		dummy /= dummy;
	}

	if (lbMmTrackBin->allocNumber != lbMmTrackNum) {
		LbInPrintf("Corruption of memory list occured, lbMmTrackBin->allocNumber  = %d, lbMmTrackNum = %ld, line = %d\n",
			lbMmTrackBin->allocNumber, (unsigned long)lbMmTrackNum, __LINE__);
			
		dummy /= dummy;
	}
}
#endif
/*********************************************************************************
*		LbMmCheckList
*
*	Arguments:
*
*	Purpose:	Runs through the memory list looking for curruption of the list.
*
*
*********************************************************************************/
void	LbMmCheckList()
{
	#ifndef LB_DEBUG
	return;
	#endif
	
	#ifdef LB_DEBUG
	LbTwoByte		index, dummy = 0;
	lbMmMemBinPtr	curMemBin = lbMmMemList;		/* Used to search for node  */
	char			*checkSumPtr;
	
	/* Run through list comparing allocation numbers */
	for (index = lbMmCurNumAllocs-1; index >= 0; index--) {

		/* Verify the check sum is valid */
		checkSumPtr = (char *)curMemBin->thePtr;
		if (memcmp(&checkSumPtr[curMemBin->byteCount], lbMmCheckSum, LBMMCHECKSUM_SIZE) != 0) {
			LbInPrintf("\nFound memory checksum corrupted, index = %d\n", index);
			dummy /= dummy;
		}
			
		if (lbMmAllocNumList[index] != curMemBin->allocNumber) {
			LbInPrintf("\nFound memory list discrepency, index = %d\n", index);
			dummy /= dummy;
		}
		curMemBin = curMemBin->nextBinPtr;
	}
	
	#endif
}

/*********************************************************************************
*		LbMmAlloc
*
*	Arguments:
*			LbUsFourByte	bytesToAlloc	- Size of memory to allocated.
*
*	Purpose:	Allocate memory.
*
*	Returns:	Ptr to memory block, 0 if failure.
*
*********************************************************************************/
void	* LbMmAlloc(LbUsFourByte bytesToAlloc)
{
	void *newMemPtr;
	
	do { /* Process Loop */
	
		/* If we are not debuging just allocate and leave */
		#ifndef LB_DEBUG
		
			/* Attempt to allocate the memory */
			if ((newMemPtr = calloc(1, bytesToAlloc)) == 0){
				sprintf(lbMmErrStr, "Unable to allocate request for %ld bytes."
					"\nSystem error number (errno) is %d.\n",
					(unsigned long)bytesToAlloc, errno);
				ErStGeneric(lbMmErrStr);
			}
			else {
				/* Clear the memory */
				memset(newMemPtr, '\0', bytesToAlloc);
			}
			return (newMemPtr);
		#endif
			
		#ifdef LB_DEBUG
		/* Verify we've been initialized */
		if (!lbMmIsInited) {
			LbInPrintf("\nYou must initialize the memory manager to call LbMmAlloc.");
			break;
		}
		
		/* Attempt to allocate the memory */
		if ((newMemPtr = calloc(1, bytesToAlloc+LBMMCHECKSUM_SIZE)) == 0){
		
			/* Create error message, more detail is available if in debug mode */				
			if (lbMmIsAccounting(lbMmInitFlags)) {
				sprintf(lbMmErrStr, "Unable to allocate request for %ld bytes.\n"
					"System error number (errno) is %d.\n"
					"Current allocation = %ld.\n",
					(unsigned long)(bytesToAlloc+LBMMCHECKSUM_SIZE), 
					errno, (unsigned long)lbMmCurBytesAllocated);
			}
			else {
				sprintf(lbMmErrStr, "Unable to allocate request for %ld bytes."
					"\nSystem error number (errno) is %d.\n",
					(unsigned long)(bytesToAlloc+LBMMCHECKSUM_SIZE), errno);
			}
		}
		ErStGeneric(lbMmErrStr);
		break;

		/* Clear the memory */
		memset(newMemPtr, '\0', bytesToAlloc);
		
		/* Set our checksum */
		tempPtr = (char *)newMemPtr;
		memcpy((void *)&tempPtr[bytesToAlloc], (void *)lbMmCheckSum, LBMMCHECKSUM_SIZE);

		/* Doing accounting then update structures */
		if (lbMmIsAccounting(lbMmInitFlags)) {
				
			/* Allocate a new instance of the memory list */
			if ((lbMmNewMemBin = (lbMmMemBinPtr) malloc(sizeof(lbMmMemBinTy))) == 0){
				sprintf(lbMmErrStr, "Unable to allocate memory list node.");
				ErStGeneric(lbMmErrStr);
				break;
			}
			
			/* Store the allocation info */
			lbMmNewMemBin->thePtr = newMemPtr;
			lbMmNewMemBin->byteCount = bytesToAlloc;
			lbMmNewMemBin->allocNumber = lbMmTotNumAllocs + 1;
			
			/* If we have set up the tracking bins number then we'll save
				a reference to this bin in a separate ptr.
			*/
			if (lbMmTrackNum == lbMmNewMemBin->allocNumber)
				lbMmTrackBin = lbMmNewMemBin;
			
			/* Add the node to the list */
			if (lbMmMemList == 0) {
				lbMmNewMemBin->nextBinPtr = 0;
				lbMmMemList = lbMmNewMemBin;
			}
			else {
				lbMmNewMemBin->nextBinPtr = lbMmMemList;
				lbMmMemList = lbMmNewMemBin;
			}
			
			/* Save the allocation number in our secondary list */
			lbMmAllocNumList[lbMmCurNumAllocs] = lbMmNewMemBin->allocNumber;
			
			/* Increment our counters */
			lbMmCurBytesAllocated += bytesToAlloc;
			lbMmTotBytesAllocated += bytesToAlloc;
			lbMmCurNumAllocs++;
			lbMmTotNumAllocs++;
			if (lbMmCurBytesAllocated > lbMmMaxResident) {
				lbMmMaxResident = lbMmCurBytesAllocated;
			}
		}
	#endif
		
	} while (false);
	
	return (newMemPtr);
	
	/* The following call can be jumped to in debugger if necessary */
	#ifdef LB_DEBUG
		if (1) {
			LbMmCheckList();
			return (newMemPtr);
		}
	#endif
}

/*********************************************************************************
*		LbMmAlloc8
*
*	Arguments:
*			LbUsEightByte	bytesToAlloc	- Size of memory to allocated.
*
*	Purpose:	Allocate memory - with an eight byte size, for larger memory
*				blocks.
*
*	Returns:	Ptr to memory block, 0 if failure.
*
*********************************************************************************/
void	* LbMmAlloc8(LbUsEightByte bytesToAlloc)
{
	void *newMemPtr;
	
	do { /* Process Loop */
	
		/* If we are not debuging just allocate and leave */
		#ifndef LB_DEBUG
		
			/* Attempt to allocate the memory */
			if ((newMemPtr = calloc(1, bytesToAlloc)) == 0){
				sprintf(lbMmErrStr, "Unable to allocate request for %ld bytes."
					"\nSystem error number (errno) is %d.\n",
					(unsigned long)bytesToAlloc, errno);
				ErStGeneric(lbMmErrStr);
			}
			else {
				/* Clear the memory */
				memset(newMemPtr, '\0', bytesToAlloc);
			}
			return (newMemPtr);
		#endif
			
		#ifdef LB_DEBUG
		/* Verify we've been initialized */
		if (!lbMmIsInited) {
			LbInPrintf("\nYou must initialize the memory manager to call LbMmAlloc.");
			break;
		}
		
		/* Attempt to allocate the memory */
		if ((newMemPtr = calloc(1, bytesToAlloc+LBMMCHECKSUM_SIZE)) == 0){
		
			/* Create error message, more detail is available if in debug mode */				
			if (lbMmIsAccounting(lbMmInitFlags)) {
				sprintf(lbMmErrStr, "Unable to allocate request for %ld bytes.\n"
					"System error number (errno) is %d.\n"
					"Current allocation = %ld.\n",
					(unsigned long)(bytesToAlloc+LBMMCHECKSUM_SIZE), 
					errno, (unsigned long)lbMmCurBytesAllocated);
			}
			else {
				sprintf(lbMmErrStr, "Unable to allocate request for %ld bytes."
					"\nSystem error number (errno) is %d.\n",
					(unsigned long)(bytesToAlloc+LBMMCHECKSUM_SIZE), errno);
			}
		}
		ErStGeneric(lbMmErrStr);
		break;

		/* Clear the memory */
		memset(newMemPtr, '\0', bytesToAlloc);
		
		/* Set our checksum */
		tempPtr = (char *)newMemPtr;
		memcpy((void *)&tempPtr[bytesToAlloc], (void *)lbMmCheckSum, LBMMCHECKSUM_SIZE);

		/* Doing accounting then update structures */
		if (lbMmIsAccounting(lbMmInitFlags)) {
				
			/* Allocate a new instance of the memory list */
			if ((lbMmNewMemBin = (lbMmMemBinPtr) malloc(sizeof(lbMmMemBinTy))) == 0){
				sprintf(lbMmErrStr, "Unable to allocate memory list node.");
				ErStGeneric(lbMmErrStr);
				break;
			}
			
			/* Store the allocation info */
			lbMmNewMemBin->thePtr = newMemPtr;
			lbMmNewMemBin->byteCount = bytesToAlloc;
			lbMmNewMemBin->allocNumber = lbMmTotNumAllocs + 1;
			
			/* If we have set up the tracking bins number then we'll save
				a reference to this bin in a separate ptr.
			*/
			if (lbMmTrackNum == lbMmNewMemBin->allocNumber)
				lbMmTrackBin = lbMmNewMemBin;
			
			/* Add the node to the list */
			if (lbMmMemList == 0) {
				lbMmNewMemBin->nextBinPtr = 0;
				lbMmMemList = lbMmNewMemBin;
			}
			else {
				lbMmNewMemBin->nextBinPtr = lbMmMemList;
				lbMmMemList = lbMmNewMemBin;
			}
			
			/* Save the allocation number in our secondary list */
			lbMmAllocNumList[lbMmCurNumAllocs] = lbMmNewMemBin->allocNumber;
			
			/* Increment our counters */
			lbMmCurBytesAllocated += bytesToAlloc;
			lbMmTotBytesAllocated += bytesToAlloc;
			lbMmCurNumAllocs++;
			lbMmTotNumAllocs++;
			if (lbMmCurBytesAllocated > lbMmMaxResident) {
				lbMmMaxResident = lbMmCurBytesAllocated;
			}
		}
	#endif
		
	} while (false);
	
	return (newMemPtr);
	
	/* The following call can be jumped to in debugger if necessary */
	#ifdef LB_DEBUG
		if (1) {
			LbMmCheckList();
			return (newMemPtr);
		}
	#endif
}

/*********************************************************************************
*		LbMmFree
*
*	Arguments:
*		void	*memPtr	- Pointer to memory to be freed.
*	Purpose:	Free memory.
*
*	Returns:	None.
*
*********************************************************************************/
/* Turn compiler optimization OFF */
#ifdef WINNT
	#pragma optimize( "", off )
#endif

void	LbMmFree(void **memPtr)
{
	#ifdef LB_DEBUG
		Boolean			found;				/* LCV For allocation number table */
		LbUsFourByte	allocIndex;			/* LCV For allocation number table */
		LbUsFourByte	index;
		lbMmMemBinPtr	nextMemBin = 0;		/* Used to search for node being released */
		lbMmMemBinPtr	curMemBin = 0;		/* Used to search for node being released */
		lbMmMemBinPtr	binToDelete = 0;	/* Temp var for node deletion */
	#endif
	
	do { /* Process Loop */
	
		#ifdef LB_DEBUG
			/* Verify we've been initialized */
			if (!lbMmIsInited) {
				LbInPrintf("\nYou must initialize the memory manager to call LbMmFree.");
				break;
			}

			if (lbMmIsAccounting(lbMmInitFlags)) {

				/* First check the list for integrety */
				LbMmCheckList();

				/* Decrement count */
				lbMmCurNumAllocs--;
				
				/* Verify they don't try to delete from an empty list */
				if (lbMmMemList == 0) {
					ErAbort("Attempt to remove memory not in mem list");
				}
				
				/* If we are removing the first node, it is easy */
				if (lbMmMemList->thePtr == *memPtr) {
					/* Find the block in our check list */
					found = false;
					for (allocIndex = 0; allocIndex <= lbMmCurNumAllocs; allocIndex++) {
						if (!found && 
								(lbMmAllocNumList[allocIndex] == lbMmMemList->allocNumber))
							found = true;
							
						if (found)
							lbMmAllocNumList[allocIndex] = lbMmAllocNumList[allocIndex+1];
					}
					binToDelete = lbMmMemList;
					lbMmMemList = lbMmMemList->nextBinPtr;
					/* Decrement our current memory count */
					lbMmCurBytesAllocated -= binToDelete->byteCount;
					free(binToDelete);
				}
				else {
				
					/* Setup local variables to search through the list */
					nextMemBin = lbMmMemList->nextBinPtr;
					curMemBin = lbMmMemList;
					index = 1;
					
					if (nextMemBin == 0) {
						ErAbort("Attempt to remove memory not in mem list");
					}
					
					/* Find this pointer and make the memory bin available */
					do {
						if (nextMemBin->thePtr == *memPtr) {
							/* Save the bin to be deleted */
							binToDelete = nextMemBin;
							
							/* Remove/compact the allocation number table */
							found = false;
							for (allocIndex = 0; allocIndex <= lbMmCurNumAllocs; allocIndex++) {
								if (!found && 
										(lbMmAllocNumList[allocIndex] == binToDelete->allocNumber))
									found = true;
									
								if (found)
									lbMmAllocNumList[allocIndex] = lbMmAllocNumList[allocIndex+1];
							}

							/* Update the current bin's nextBinPtr */
							curMemBin->nextBinPtr = nextMemBin->nextBinPtr;

							/* Decrement our current memory count */
							lbMmCurBytesAllocated -= binToDelete->byteCount;
							
							/* Free the deleted bin */
							free(binToDelete);
							
							/* Finish the loop */
							break;
						}
						else {
						
							/* Go to the next bin */
							curMemBin = nextMemBin;
							nextMemBin = nextMemBin->nextBinPtr;
							index++;
							
							/* Verify we didn't get to the end of the list without finding the bin */
							if ((nextMemBin == 0) || (index > lbMmCurNumAllocs)) {
								ErAbort("Attempt to remove memory not in mem list");
							}
						}
					} while (true);
				}
		}
			
		#endif
	
		/* Free the memory */
		if (*memPtr != 0)
			free(*memPtr);
		
		/* Clear the pointer */
		*memPtr = (void *) 0;
		
	} while (false);
}

/*  Turn compiler optimization back ON */
#ifdef WINNT
	#pragma optimize( "", on )
#endif


/*********************************************************************************
*		LbMmInit
*
*	Purpose:	Initialize the manager.
*	Arguments:
*		LbUsFourByte	flags	- Behavior modification flags
*
*	Returns:	TRUE if initialization succeeds.
*
*********************************************************************************/
Boolean	LbMmInit(LbUsFourByte flags)
{
	/* Remove compiler warning */
	#ifndef LB_DEBUG
		if (flags) {}
	#endif
	
	do	/* Abort Loop */
	{
		#ifdef LB_DEBUG

			/* Clear intense debugging vars */
			dummyBin.allocNumber = 1;
			lbMmTrackBin = &dummyBin;
			lbMmTrackNum = 1;
		
			/* Indicate we've been initialized */
			lbMmIsInited = true;
			lbMmInitFlags = flags;

			if (lbMmIsAccounting(lbMmInitFlags)) {
				
				/* Clear our counters */
				lbMmCurBytesAllocated = 0;
				lbMmTotBytesAllocated = 0;
				lbMmCurNumAllocs = 0;
				lbMmTotNumAllocs = 0;
				lbMmMaxResident = 0;
			}
		
			/* Indicate we've been initialized */
			lbMmIsInited = true;
			lbMmInitFlags = flags;
	
		#endif
			
	} while (false);
	
	#ifdef LB_DEBUG
		return (lbMmIsInited);
	#else
		return (true);
	#endif
}

/*********************************************************************************
*		LbMmTerminate
*
*	Arguments:
*	Purpose:	Terminate the memory manager.
*
*	Returns:	None.
*
*********************************************************************************/
void	LbMmTerminate()
{
	
	#ifdef LB_DEBUG
	{
		lbMmMemBinPtr	nextMemBin = 0;		/* Used to search for node being released */
		lbMmMemBinPtr	binToDelete = 0;	/* Temp var for node deletion */
		LbFourByte	numFrees = 0;			/* Temp to help with debugging */
		
		if (lbMmIsInited)
		{
			/* Verify no memory was left allocated */
			if (lbMmIsAccounting(lbMmInitFlags)) {
			
				/* See if we freed all of the memory allocated */
				if (lbMmCurBytesAllocated != 0) {
					LbInPrintf(
					"\nMemory allocations != memory frees.\nPtrs left allocated = %ld\nBytes left allocated = %ld\n."
					"Memory being freed\n",
					(unsigned long)lbMmCurNumAllocs, (unsigned long)lbMmCurBytesAllocated);
					
					/* First check the list for integrety */
					LbMmCheckList();

					/* Now free any residual memory */
					binToDelete = lbMmMemList;
					nextMemBin = binToDelete->nextBinPtr;
					while (binToDelete != 0) {
						free((void *) binToDelete->thePtr);
						free((void *) binToDelete);
						numFrees++;
						binToDelete = nextMemBin;
						if (binToDelete != 0)
							nextMemBin = binToDelete->nextBinPtr;
					}
				}
				
			}
			
			/* The following is just a kludge to get the MW debugger
				to display these local globals when a break point is
				set in this routine
			*/
			lbMmTotNumAllocs = 0;		
			lbMmMaxResident = 0;
			numFrees = -1;
		}
	
		/* Terminate the manager */
		lbMmIsInited = false;
	}
	#endif
	
}
#undef LB_MEMORY
