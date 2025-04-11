/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 2005-2012 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/

/*********************************************************************************
*
*			Module Name:		LbSort.h
*			Revision Number:	1.3
*			Date last revised:	27 August 2012
*			Programmer:			Steven Gillispie 
*			Date Originated:	22 February 2005
*
*			Module Overview:	Declarations for sorted list utilities
*
*			References:			None
*
**********************************************************************************
*
*			Global functions defined:	None
*
*			Global variables defined:	None
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
*			Revision date:		27 August 2012
*
*			Revision description:	Corrected sort indicators to be signed
*
*********************************************************************************/

#ifndef LB_SORT
#define LB_SORT

#include "SystemDependent.h"
#include "LbTypes.h"


/* CONSTANTS */
#define		LBSORT_KEYSIZE	8		/* Maximum byte length of a sorted list item sort key */


/* TYPES */
typedef union	{			/* Sorted list item sort key type */
	LbUsOneByte			as1Byte;
	LbUsOneByte			asByte[LBSORT_KEYSIZE];
} LbSortKeyType;
typedef LbSortKeyType		*LbSortKeyPtr;			/* Pointer to a sort key */

typedef LbFourByte (LbSortKeyComparator)(LbSortKeyPtr keyPtr1, LbSortKeyPtr keyPtr2);
	/* A function of this type is needed by the LbSort functions since they 
		cannot know how to compare the specific type of keys being used. */

typedef LbUsFourByte	LbSortItem;			/* Sorted list item type (actual type is private) */
typedef LbSortItem		*LbSortItemPtr;		/* Pointer to sorted list item */

typedef LbUsFourByte	LbSortList;			/* Sorted list (actual type is private) */
typedef LbSortList		*LbSortListPtr;		/* Pointer to sorted list */


/* PROTOTYPES */
void			LbSortGetItemSortKey(LbSortItemPtr theListItemPtr, 
					LbSortKeyPtr itsSortKeyPtr, 
					LbUsFourByte sortKeyLength);
LbUsFourByte	LbSortGetItemIndex(LbSortItemPtr theListItemPtr);
void			LbSortSetItemIndex(LbSortItemPtr theListItemPtr, LbUsFourByte newValue);
void*			LbSortGetItemDataPtr(LbSortItemPtr theListItemPtr);
void			LbSortSetItemDataPtr(LbSortItemPtr theListItemPtr, void *newPointer);

LbSortListPtr	LbSortNewList(LbUsFourByte numItems, 
					LbSortKeyPtr defaultKeyPtr, 
					LbUsFourByte sortKeyLength, 
					LbSortKeyComparator *keyComparatorPtr);
Boolean			LbSortInsert(LbSortListPtr theListPtr, 
					LbSortKeyPtr dataSortKeyPtr, 
					LbUsFourByte indexData, void *dataPtr);
Boolean			LbSortDelete(LbSortListPtr theListPtr, LbSortItemPtr theItemPtr);
Boolean			LbSortFirst(LbSortListPtr theListPtr, 
					LbSortItemPtr *firstItemPtr, 
					LbUsFourByte *index, void **dataPtr);
Boolean			LbSortLast(LbSortListPtr theListPtr, 
					LbSortItemPtr *lastItemPtr, 
					LbUsFourByte *index, void **dataPtr);
Boolean			LbSortNext(LbSortListPtr theListPtr, 
					LbSortItemPtr startItemPtr, LbSortItemPtr *nextItemPtr, 
					LbUsFourByte *index, void **dataPtr);
Boolean			LbSortPrev(LbSortListPtr theListPtr, 
					LbSortItemPtr startItemPtr, LbSortItemPtr *prevItemPtr, 
					LbUsFourByte *index, void **dataPtr);
Boolean			LbSortFind(LbSortListPtr theListPtr, 
					LbSortKeyPtr minSortKeyPtr, LbSortItemPtr *sortItemPtr);
void			LbSortDisposeList(LbSortListPtr *theListPtrPtr);


#endif /* LB_SORT */
