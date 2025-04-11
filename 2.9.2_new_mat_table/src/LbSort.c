/*********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                 (C) Copyright 2005 Department of Radiology	   		    	*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/*********************************************************************************
*
*		Module Name:		LbSort.c
*		Revision Number:	1.4
*		Date last revised:	12 August 2005
*		Programmer:			Steven Gillispie 
*		Date Originated:	22 February 2005
*
*		Module Overview:	Functions for sorted list utilities
*
*		References:			Cormen, Leiserson, & Rivest (1990), chaps 13 & 14:
*								sort trees and red-black trees
*
**********************************************************************************
*
*		Global functions defined:	
*			LbSortGetItemSortKey
*			LbSortGetItemIndex
*			LbSortSetItemIndex
*			LbSortGetItemDataPtr
*			LbSortSetItemDataPtr
*			LbSortNewList
*			LbSortInsert
*			LbSortDelete
*			LbSortFirst
*			LbSortLast
*			LbSortNext
*			LbSortPrev
*			LbSortFind
*			LbSortDisposeList
*
*		Global variables defined:	none
*
**********************************************************************************
*
*		Revision Section (Also update version number, if relevant)
*
*		Programmer(s):		
*
*		Revision date:		
*
*		Revision description:
*
*********************************************************************************/

#include "LbSort.h"


#include <string.h>

#include "LbMemory.h"


/*	LOCAL CONSTANTS */
#define		LBSORT_Black	0				/* RBTree black node color */
#define		LBSORT_Red		1				/* RBTree red node color */


/*  LOCAL TYPES */
struct RBTreeNode  {		/* Red-Black sort tree node */
	LbSortKeyType		sortKey;			/* Node sort key */
	LbUsFourByte		index;				/* Integer data for the node */
	void	 			*dataPtr;			/* Pointer data for the node */
	LbUsFourByte		nodeColor;			/* Node color (red or black) */
	struct RBTreeNode	*parent;			/* Pointer to the parent node */
	struct RBTreeNode	*leftBranch;		/* Pointer to the left descending branch */
	struct RBTreeNode	*rightBranch;		/* Pointer to the right descending branch */
};
typedef struct RBTreeNode	RBTreeNode;
typedef RBTreeNode			*RBTreeNodePtr;

typedef struct  {			/* Sorted list (sort tree) data structure */
	RBTreeNode			rbTreeNilNode;		/* Sort tree nil node (avoids NULL dereferences) */
	RBTreeNodePtr		rbTreeNil;			/* Pointer to rbTreeNilNode */
	RBTreeNodePtr		rbTreeStartPtr;		/* Pointer to start of sort tree memory */
	RBTreeNodePtr		rbTreePtr;			/* Red-Black sort tree root node */
	RBTreeNodePtr		rbTreeAvailNodes;	/* List of available nodes within tree */
	LbSortKeyType		defaultSortKey;		/* Default (initial) sort key value */
	LbUsFourByte		sortKeyLength;		/* Byte length of used sort keys */
	LbSortKeyComparator	*keyComparatorPtr;	/* Pointer to sort key comparator */
} LbSortedList;
typedef LbSortedList		*LbSortedListPtr;


/*  LOCAL GLOBALS */

/*	LOCAL MACROS */

/*	LOCAL FUNCTIONS */
Boolean			lbSortNewBlock(LbSortedListPtr theListPtr, LbUsFourByte blockItems, 
					RBTreeNodePtr *blockStartPtr, RBTreeNodePtr *lastNodePtr);
RBTreeNodePtr	lbSortRBTreeNewNode(LbSortedListPtr theListPtr);
void			lbSortRBTreeFreeNode(LbSortedListPtr theListPtr, RBTreeNodePtr returnedNode);
void			lbSortRBTreeLeftRotate(LbSortedListPtr theListPtr, RBTreeNodePtr theNodePtr);
void			lbSortRBTreeRightRotate(LbSortedListPtr theListPtr, RBTreeNodePtr theNodePtr);


/*	FUNCTIONS */

/*********************************************************************************
*
*	Name:			LbSortGetItemSortKey
*
*	Summary:		Return the sorted list item's sort key.
*
*	Arguments:		
*		LbSortItemPtr		theListItemPtr	- The sorted list item pointer.
*		LbSortKeyPtr		itsSortKeyPtr	- Pointer to returned sort key.
*		LbUsFourByte		sortKeyLength	- Bytes in the sort key.
*
*	Function return: None.
*
*********************************************************************************/
void LbSortGetItemSortKey(LbSortItemPtr theListItemPtr, 
							LbSortKeyPtr itsSortKeyPtr, 
							LbUsFourByte sortKeyLength)
{
	/* Copy theListItemPtr's sort key to the supplied memory location */
	memcpy(itsSortKeyPtr, 
			&(((RBTreeNodePtr)theListItemPtr)->sortKey), 
			sortKeyLength);
}

/*********************************************************************************
*
*	Name:			LbSortGetItemIndex
*
*	Summary:		Return the sorted list item's index.
*
*	Arguments:		
*		LbSortItemPtr		theListItemPtr	- The sorted list item pointer.
*
*	Function return: The LbUsFourByte index value.
*
*********************************************************************************/
LbUsFourByte LbSortGetItemIndex(LbSortItemPtr theListItemPtr)
{
	return (((RBTreeNodePtr)theListItemPtr)->index);
}

/*********************************************************************************
*
*	Name:			LbSortSetItemIndex
*
*	Summary:		Set the sorted list item's index value.
*
*	Arguments:		
*		LbSortItemPtr		theListItemPtr	- The sorted list item pointer.
*		LbUsFourByte		newValue		- The new index value.
*
*	Function return: None.
*
*********************************************************************************/
void LbSortSetItemIndex(LbSortItemPtr theListItemPtr, LbUsFourByte newValue)
{
	((RBTreeNodePtr)theListItemPtr)->index = newValue;
}

/*********************************************************************************
*
*	Name:			LbSortGetItemDataPtr
*
*	Summary:		Return the sorted list item's data pointer.
*
*	Arguments:		
*		LbSortItemPtr		theListItemPtr	- The sorted list item pointer.
*
*	Function return: The void* data pointer.
*
*********************************************************************************/
void* LbSortGetItemDataPtr(LbSortItemPtr theListItemPtr)
{
	return (((RBTreeNodePtr)theListItemPtr)->dataPtr);
}

/*********************************************************************************
*
*	Name:			LbSortSetItemDataPtr
*
*	Summary:		Set the sorted list item's data pointer.
*
*	Arguments:		
*		LbSortItemPtr		theListItemPtr	- The sorted list item pointer.
*		void				*newPointer		- The new data pointer.
*
*	Function return: None.
*
*********************************************************************************/
void LbSortSetItemDataPtr(LbSortItemPtr theListItemPtr, void *newPointer)
{
	((RBTreeNodePtr)theListItemPtr)->dataPtr = newPointer;
}


/*********************************************************************************
*
*	Name:			LbSortNewList
*
*	Summary:		Return a new, initialized sorted list.
*
*	Arguments:		
*		LbUsFourByte		numItems			- Items in a list block.
*		LbSortKeyPtr		defaultKeyPtr		- Initial item sort key pointer.
*		LbUsFourByte		sortKeyLength		- Bytes in a list item sort key.
*		LbSortKeyComparator	*keyComparatorPtr	- Compare function for the sort keys.
*
*	Function return: LbSortListPtr.
*
*********************************************************************************/
LbSortListPtr LbSortNewList(LbUsFourByte numItems, 
							LbSortKeyPtr defaultKeyPtr, 
							LbUsFourByte sortKeyLength, 
							LbSortKeyComparator *keyComparatorPtr)
{
	LbSortedListPtr		newListPtr = NULL;		/* Pointer to new returned list */
	RBTreeNodePtr		blockStartPtr;			/* Pointer to sort tree block */
	RBTreeNodePtr		lastNodePtr;			/* Ignored */
	
	
	/* Do initial check for valid parameters */
	if (sortKeyLength > LBSORT_KEYSIZE) {
		/* Invalid parameters; do nothing */
	}
	else {
		/* Allocate memory for the sorted list data structure */
		newListPtr = (LbSortedListPtr) LbMmAlloc(sizeof(LbSortedList));
		
		if (newListPtr) {
			/* Initialize the new sorted list */
			
			/* Set sort key size, default sort key, and key comparator */
			newListPtr->sortKeyLength = sortKeyLength;
			memcpy(&(newListPtr->defaultSortKey), defaultKeyPtr, sortKeyLength);
			newListPtr->keyComparatorPtr = keyComparatorPtr;
			
			/* Set rbTreeNil to point to rbTreeNilNode */
			newListPtr->rbTreeNil = &(newListPtr->rbTreeNilNode);
			
			/* Initialize the sentinel nil node to null values */
			memcpy(&(newListPtr->rbTreeNilNode.sortKey), defaultKeyPtr, sortKeyLength);
			newListPtr->rbTreeNilNode.index = -1;
			newListPtr->rbTreeNilNode.dataPtr = NULL;
			newListPtr->rbTreeNilNode.nodeColor = LBSORT_Black;
			newListPtr->rbTreeNilNode.parent = newListPtr->rbTreeNil;
			newListPtr->rbTreeNilNode.leftBranch = newListPtr->rbTreeNil;
			newListPtr->rbTreeNilNode.rightBranch = newListPtr->rbTreeNil;
			
			/* Set the tree to an empty tree */
			newListPtr->rbTreeStartPtr = newListPtr->rbTreeNil;
			newListPtr->rbTreePtr = newListPtr->rbTreeNil;
			newListPtr->rbTreeAvailNodes = newListPtr->rbTreeNil;
			
			/* Create the first block of the sort tree */
			if (lbSortNewBlock(newListPtr, numItems, 
									&blockStartPtr, &lastNodePtr)) {
				/* Finish initializing the sort tree */
				newListPtr->rbTreeStartPtr = blockStartPtr;
				newListPtr->rbTreeAvailNodes = blockStartPtr;
			}
			else {
				/* Couldn't allocate memory, so free the sorted list */
				LbMmFree((void **) &newListPtr);
			}
		}
	}
	
	/* Return the new list */
	return ((LbSortListPtr) newListPtr);
}

/*********************************************************************************
*
*	Name:			lbSortNewBlock
*
*	Summary:		Return a new, initialized sorted list block.
*
*	Arguments:		
*		LbSortedListPtr		theListPtr		- The relevant sorted list.
*		LbUsFourByte		blockItems		- Number of items in a list block.
*		RBTreeNodePtr		*blockStartPtr	- Returned pointer to new block.
*		RBTreeNodePtr		*lastNodePtr	- Returned pointer to last block node.
*
*	Function return: True if successful.
*
*********************************************************************************/
Boolean lbSortNewBlock(LbSortedListPtr theListPtr, LbUsFourByte blockItems, 
						RBTreeNodePtr *blockStartPtr, RBTreeNodePtr *lastNodePtr)
{
	Boolean				result;					/* Result of function */
	RBTreeNodePtr		listBlockPtr;			/* Pointer to sort tree block */
	RBTreeNodePtr		blockNodePtr;			/* Node pointer into listBlockPtr */
	RBTreeNodePtr		blockNodeParPtr;		/* Pointer to blockNodePtr parent */
	LbUsFourByte		i;						/* Generic index variable */
	
	
	result = true;
	
	/* Allocate memory for the new sort tree block */
	if ((listBlockPtr = (RBTreeNodePtr) 
			LbMmAlloc(blockItems*sizeof(RBTreeNode))) == NULL) {
		result = false;
	}
	
	if (result) {
		/* Record the start and end nodes */
		*blockStartPtr = listBlockPtr;
		*lastNodePtr = listBlockPtr + blockItems - 1;
		
		/* Link the sort tree block nodes into a list (start at end) */
		blockNodePtr = *lastNodePtr;
		blockNodeParPtr = blockNodePtr - 1;
		blockNodePtr->rightBranch = theListPtr->rbTreeNil;
		for (i=1; i<blockItems; i++) {
			blockNodePtr->parent = blockNodeParPtr;
			blockNodeParPtr->rightBranch = blockNodePtr;
			blockNodePtr = blockNodeParPtr;
			blockNodeParPtr--;
		}
	}
	
	return (result);
}

/*********************************************************************************
*
*	Name:			lbSortRBTreeNewNode
*
*	Summary:		Obtain a new node for the Red-Black sort tree.
*
*	Arguments:		
*		LbSortedListPtr		theListPtr		- The relevant sorted list.
*
*	Function return: RBTreeNodePtr; NULL if couldn't be created.
*
*********************************************************************************/
RBTreeNodePtr lbSortRBTreeNewNode(LbSortedListPtr theListPtr)
{
	RBTreeNodePtr		theNodePtr = NULL;		/* Returned node pointer */
	
	
	if (theListPtr->rbTreeAvailNodes != theListPtr->rbTreeNil) {
		/* Get one from the used pile */
		theNodePtr = theListPtr->rbTreeAvailNodes;
		theListPtr->rbTreeAvailNodes = theListPtr->rbTreeAvailNodes->rightBranch;
	}
	else {
		/* No more node space--shouldn't happen if list is big enough */
		theNodePtr = NULL;
		
		/* NOTE:  If dynamic list resizing were implemented, at this point 
			lbSortNewBlock would be called and the new block linked into 
			a list of node blocks and the list of available nodes.  This would 
			also require that the list memory be implemented as a linked list 
			in LbSortNewList and disposed of in LbSortDisposeList. */
	}
	
	if (theNodePtr != NULL) {
		/* Initialize the node */
		memcpy(&(theNodePtr->sortKey), &(theListPtr->defaultSortKey), 
					theListPtr->sortKeyLength);
		theNodePtr->index = 0;
		theNodePtr->dataPtr = NULL;
		theNodePtr->nodeColor = LBSORT_Red;
		theNodePtr->parent = theListPtr->rbTreeNil;
		theNodePtr->leftBranch = theListPtr->rbTreeNil;
		theNodePtr->rightBranch = theListPtr->rbTreeNil;
	}
	
	return (theNodePtr);
}

/*********************************************************************************
*
*	Name:			lbSortRBTreeFreeNode
*
*	Summary:		Return a deleted node back to the Red-Black sort tree free list.
*					NOTE:  The node is assumed to already be removed from the tree.
*
*	Arguments:		
*		LbSortedListPtr		theListPtr		- The relevant sorted list.
*		RBTreeNodePtr		returnedNode	- Returned tree node.
*
*	Function return: None.
*
*********************************************************************************/
void lbSortRBTreeFreeNode(LbSortedListPtr theListPtr, RBTreeNodePtr returnedNode)
{
	/* Add the node to the list of available nodes */
	theListPtr->rbTreeAvailNodes->parent = returnedNode;
	returnedNode->rightBranch = theListPtr->rbTreeAvailNodes;
	theListPtr->rbTreeAvailNodes = returnedNode;
}

/*********************************************************************************
*
*	Name:			lbSortRBTreeLeftRotate
*
*	Summary:		Rotate nodes leftward in the Red-Black sort tree.
*
*	Arguments:		
*		LbSortedListPtr		theListPtr		- The relevant sorted list.
*		RBTreeNodePtr		theNodePtr		- The node to be rotated leftward.
*
*	Function return: None.
*
*********************************************************************************/
void lbSortRBTreeLeftRotate(LbSortedListPtr theListPtr, RBTreeNodePtr theNodePtr)
{
	RBTreeNodePtr	parPtr;			/* Pointer to parent of theNodePtr */
	RBTreeNodePtr	rightPtr;		/* Pointer to theNodePtr's right subtree */
	
	
	parPtr = theNodePtr->parent;
	rightPtr = theNodePtr->rightBranch;
	
	theNodePtr->rightBranch = rightPtr->leftBranch;
	if (rightPtr->leftBranch != theListPtr->rbTreeNil) {
		rightPtr->leftBranch->parent = theNodePtr;
	}
	rightPtr->parent = parPtr;
	if (parPtr == theListPtr->rbTreeNil) {
		/* New root */
		theListPtr->rbTreePtr = rightPtr;
	}
	else {
		if (theNodePtr == parPtr->leftBranch) {
			parPtr->leftBranch = rightPtr;
		}
		else {
			parPtr->rightBranch = rightPtr;
		}
	}
	rightPtr->leftBranch = theNodePtr;
	theNodePtr->parent = rightPtr;
}

/*********************************************************************************
*
*	Name:			lbSortRBTreeRightRotate
*
*	Summary:		Rotate nodes rightward in the Red-Black sort tree.
*
*	Arguments:		
*		LbSortedListPtr		theListPtr		- The relevant sorted list.
*		RBTreeNodePtr		theNodePtr		- The node to be rotated rightward.
*
*	Function return: None.
*
*********************************************************************************/
void lbSortRBTreeRightRotate(LbSortedListPtr theListPtr, RBTreeNodePtr theNodePtr)
{
	RBTreeNodePtr	parPtr;			/* Pointer to parent of theNodePtr */
	RBTreeNodePtr	leftPtr;		/* Pointer to theNodePtr's left subtree */
	
	
	parPtr = theNodePtr->parent;
	leftPtr = theNodePtr->leftBranch;
	
	theNodePtr->leftBranch = leftPtr->rightBranch;
	if (leftPtr->rightBranch != theListPtr->rbTreeNil) {
		leftPtr->rightBranch->parent = theNodePtr;
	}
	leftPtr->parent = parPtr;
	if (parPtr == theListPtr->rbTreeNil) {
		/* New root */
		theListPtr->rbTreePtr = leftPtr;
	}
	else {
		if (theNodePtr == parPtr->rightBranch) {
			parPtr->rightBranch = leftPtr;
		}
		else {
			parPtr->leftBranch = leftPtr;
		}
	}
	leftPtr->rightBranch = theNodePtr;
	theNodePtr->parent = leftPtr;
}

/*********************************************************************************
*
*	Name:			LbSortInsert
*
*	Summary:		Insert the given data into the sorted list.
*
*	Arguments:
*		LbSortListPtr		theListPtr		- The relevant sorted list.
*		LbSortKeyPtr		dataSortKeyPtr	- The data sort key pointer.
*		LbUsFourByte		indexData		- Integer data to insert.
*		void				*dataPtr		- Pointer data to insert.
*
*	Function return: True if successful.
*
*********************************************************************************/
Boolean LbSortInsert(LbSortListPtr theListPtr, 
						LbSortKeyPtr dataSortKeyPtr, 
						LbUsFourByte indexData, void *dataPtr)
{
	LbSortedListPtr	localListPtr;	/* theListPtr as LbSortedListPtr */
	RBTreeNodePtr	newNodePtr;		/* Pointer to newly added node */
	LbSortKeyComparator
					*CompareKeys;	/* Sort key comparator */
	RBTreeNodePtr	curPtr;			/* Pointer to current sort tree node */
	RBTreeNodePtr	parPtr;			/* Pointer to parent of curPtr */
	RBTreeNodePtr	grParPtr;		/* Pointer to parent of parPtr */
	RBTreeNodePtr	auntPtr;		/* Pointer to other child of grParPtr */
	
	
	localListPtr = (LbSortedListPtr) theListPtr;
	
	/* Get a new node, if possible */
	newNodePtr = lbSortRBTreeNewNode(localListPtr);
	
	if (newNodePtr != NULL) {
		/* Fill in the new node's data */
		memcpy(&(newNodePtr->sortKey), dataSortKeyPtr, localListPtr->sortKeyLength);
		newNodePtr->index = indexData;
		newNodePtr->dataPtr = dataPtr;
		
		if (localListPtr->rbTreePtr == localListPtr->rbTreeNil) {
			/* First node into the tree */
			localListPtr->rbTreePtr = newNodePtr;
			newNodePtr->nodeColor = LBSORT_Black;
		}
		else {
			/* Search the sort tree for the proper insertion point */
			CompareKeys = localListPtr->keyComparatorPtr;
			curPtr = localListPtr->rbTreePtr;
			do {
				parPtr = curPtr;
				if ((*CompareKeys)(dataSortKeyPtr, &(curPtr->sortKey)) == -1) {
					/* Left branch */
					curPtr = curPtr->leftBranch;
				}
				else {
					/* Right branch; note: includes == case */
					curPtr = curPtr->rightBranch;
				}
			} while (curPtr != localListPtr->rbTreeNil);
			
			/* Set the inserted node's parent */
			newNodePtr->parent = parPtr;
			
			/* Set the parent's child branch */
			if ((*CompareKeys)(dataSortKeyPtr, &(parPtr->sortKey)) == -1) {
				/* Left branch */
				parPtr->leftBranch = newNodePtr;
			}
			else {
				/* Right branch; note: includes == case */
				parPtr->rightBranch = newNodePtr;
			}
			
			/* Adjust the node colors to accomodate the new node */
			newNodePtr->nodeColor = LBSORT_Red;
			curPtr = newNodePtr;
			grParPtr = parPtr->parent;
			while ((curPtr != localListPtr->rbTreePtr) && 
					(parPtr->nodeColor == LBSORT_Red)) {
				if (parPtr == grParPtr->leftBranch) {
					auntPtr = grParPtr->rightBranch;
					if (auntPtr->nodeColor == LBSORT_Red) {
						parPtr->nodeColor = LBSORT_Black;
						auntPtr->nodeColor = LBSORT_Black;
						grParPtr->nodeColor = LBSORT_Red;
						curPtr = grParPtr;
						parPtr = curPtr->parent;
						grParPtr = parPtr->parent;
					}
					else {
						if (curPtr == parPtr->rightBranch) {
							curPtr = parPtr;
							lbSortRBTreeLeftRotate(localListPtr, curPtr);
							parPtr = curPtr->parent;
							grParPtr = parPtr->parent;
						}
						parPtr->nodeColor = LBSORT_Black;
						grParPtr->nodeColor = LBSORT_Red;
						lbSortRBTreeRightRotate(localListPtr, grParPtr);
						parPtr = curPtr->parent;
						grParPtr = parPtr->parent;
					}
				}
				else {
					auntPtr = grParPtr->leftBranch;
					if (auntPtr->nodeColor == LBSORT_Red) {
						parPtr->nodeColor = LBSORT_Black;
						auntPtr->nodeColor = LBSORT_Black;
						grParPtr->nodeColor = LBSORT_Red;
						curPtr = grParPtr;
						parPtr = curPtr->parent;
						grParPtr = parPtr->parent;
					}
					else {
						if (curPtr == parPtr->leftBranch) {
							curPtr = parPtr;
							lbSortRBTreeRightRotate(localListPtr, curPtr);
							parPtr = curPtr->parent;
							grParPtr = parPtr->parent;
						}
						parPtr->nodeColor = LBSORT_Black;
						grParPtr->nodeColor = LBSORT_Red;
						lbSortRBTreeLeftRotate(localListPtr, grParPtr);
						parPtr = curPtr->parent;
						grParPtr = parPtr->parent;
					}
				}
			}
			localListPtr->rbTreePtr->nodeColor = LBSORT_Black;
		}
	}
	
	return (newNodePtr != NULL);
}

/*********************************************************************************
*
*	Name:			LbSortDelete
*
*	Summary:		Remove the given item from the sorted list.
*
*	Arguments:
*		LbSortListPtr		theListPtr		- The relevant sorted list.
*		LbSortItemPtr		theItemPtr		- Pointer to the sorted list item.
*
*	Function return: True if successful.
*
*********************************************************************************/
Boolean LbSortDelete(LbSortListPtr theListPtr, LbSortItemPtr theItemPtr)
{
	LbSortedListPtr	localListPtr;	/* theListPtr as LbSortedListPtr */
	RBTreeNodePtr	theNodePtr;		/* RBTreeNodePtr version of theItemPtr */
	RBTreeNodePtr	delPtr;			/* Pointer to the node to remove */
	RBTreeNodePtr	curPtr;			/* Pointer to current sort tree node */
	RBTreeNodePtr	parPtr;			/* Pointer to parent of curPtr */
	RBTreeNodePtr	childPtr;		/* Pointer to child of delPtr */
	RBTreeNodePtr	rightPtr;		/* Pointer to right sibling of curPtr */
	RBTreeNodePtr	leftPtr;		/* Pointer to left sibling of curPtr */
	
	
	localListPtr = (LbSortedListPtr) theListPtr;
	theNodePtr = (RBTreeNodePtr) theItemPtr;
	
	/* Determine the memory node to delete */
	if ((theNodePtr->leftBranch == localListPtr->rbTreeNil) || 
			(theNodePtr->rightBranch == localListPtr->rbTreeNil)) {
		delPtr = theNodePtr;
	}
	else {
		/* Must find next largest sort key node */
		curPtr = theNodePtr->rightBranch;
		while (curPtr->leftBranch != localListPtr->rbTreeNil) {
			curPtr = curPtr->leftBranch;
		}
		delPtr = curPtr;
	}
	
	/* Remove the delPtr node from the tree */
	if (delPtr->leftBranch != localListPtr->rbTreeNil) {
		childPtr = delPtr->leftBranch;
	}
	else {
		childPtr = delPtr->rightBranch;
	}
	childPtr->parent = delPtr->parent;	/* childPtr can be rbTreeNil */
	if (delPtr->parent == localListPtr->rbTreeNil) {
		/* New root */
		localListPtr->rbTreePtr = childPtr;
	}
	else {
		if (delPtr == delPtr->parent->leftBranch) {
			delPtr->parent->leftBranch = childPtr;
		}
		else {
			delPtr->parent->rightBranch = childPtr;
		}
	}
	if (delPtr != theNodePtr) {
		/* Transfer the substituted data to its new position */
		memcpy(&(theNodePtr->sortKey), &(delPtr->sortKey), 
				localListPtr->sortKeyLength);
		theNodePtr->index = delPtr->index;
		theNodePtr->dataPtr = delPtr->dataPtr;
	}
	
	if (delPtr->nodeColor == LBSORT_Black) {
		/* Fix up the node colors by pushing the black color up the tree */
		curPtr = childPtr;
		parPtr = curPtr->parent;
		while ((curPtr != localListPtr->rbTreePtr) && 
					(curPtr->nodeColor == LBSORT_Black)) {
			if (curPtr == parPtr->leftBranch) {
				rightPtr = parPtr->rightBranch;
				if (rightPtr->nodeColor == LBSORT_Red) {
					rightPtr->nodeColor = LBSORT_Black;
					parPtr->nodeColor = LBSORT_Red;
					lbSortRBTreeLeftRotate(localListPtr, parPtr);
					rightPtr = parPtr->rightBranch;
				}
				if ((rightPtr->leftBranch->nodeColor == LBSORT_Black) && 
						(rightPtr->rightBranch->nodeColor == LBSORT_Black)) {
					rightPtr->nodeColor = LBSORT_Red;
					curPtr = parPtr;
					parPtr = parPtr->parent;
				}
				else {
					if (rightPtr->rightBranch->nodeColor == LBSORT_Black) {
						rightPtr->leftBranch->nodeColor = LBSORT_Black;
						rightPtr->nodeColor = LBSORT_Red;
						lbSortRBTreeRightRotate(localListPtr, rightPtr);
						rightPtr = parPtr->rightBranch;
					}
					rightPtr->nodeColor = parPtr->nodeColor;
					parPtr->nodeColor = LBSORT_Black;
					rightPtr->rightBranch->nodeColor = LBSORT_Black;
					lbSortRBTreeLeftRotate(localListPtr, parPtr);
					curPtr = localListPtr->rbTreePtr;
				}
			}
			else {
				leftPtr = parPtr->leftBranch;
				if (leftPtr->nodeColor == LBSORT_Red) {
					leftPtr->nodeColor = LBSORT_Black;
					parPtr->nodeColor = LBSORT_Red;
					lbSortRBTreeRightRotate(localListPtr, parPtr);
					leftPtr = parPtr->leftBranch;
				}
				if ((leftPtr->rightBranch->nodeColor == LBSORT_Black) && 
						(leftPtr->leftBranch->nodeColor == LBSORT_Black)) {
					leftPtr->nodeColor = LBSORT_Red;
					curPtr = parPtr;
					parPtr = parPtr->parent;
				}
				else {
					if (leftPtr->leftBranch->nodeColor == LBSORT_Black) {
						leftPtr->rightBranch->nodeColor = LBSORT_Black;
						leftPtr->nodeColor = LBSORT_Red;
						lbSortRBTreeLeftRotate(localListPtr, leftPtr);
						leftPtr = parPtr->leftBranch;
					}
					leftPtr->nodeColor = parPtr->nodeColor;
					parPtr->nodeColor = LBSORT_Black;
					leftPtr->leftBranch->nodeColor = LBSORT_Black;
					lbSortRBTreeRightRotate(localListPtr, parPtr);
					curPtr = localListPtr->rbTreePtr;
				}
			}
		}
		curPtr->nodeColor = LBSORT_Black;
	}
	
	/* Make sure rbTreeNil's parent is set to rbTreeNil */
	localListPtr->rbTreeNil->parent = localListPtr->rbTreeNil;
	
	/* Return the deleted node to the available list */
	lbSortRBTreeFreeNode(localListPtr, delPtr);
	
	/* At this time, no error seems possible */
	return (true);
}

/*********************************************************************************
*
*	Name:			LbSortFirst
*
*	Summary:		Return the first sorted list item, index data, and data pointer. 
*
*	Arguments:
*		LbSortListPtr		theListPtr		- The relevant sorted list.
*		LbSortItemPtr		*firstItemPtr	- Pointer to the first sorted list item.
*		LbUsFourByte		*index			- First index data.
*		void				**dataPtr		- First data pointer.
*
*	Function return: True if successful.
*
*********************************************************************************/
Boolean LbSortFirst(LbSortListPtr theListPtr, 
						LbSortItemPtr *firstItemPtr, 
						LbUsFourByte *index, void **dataPtr)
{
	Boolean			found;			/* Item was found; function return value */
	LbSortedListPtr	localListPtr;	/* theListPtr as LbSortedListPtr */
	RBTreeNodePtr	curPtr;			/* Pointer to current sort node */
	
	
	localListPtr = (LbSortedListPtr) theListPtr;
	
	if (localListPtr->rbTreePtr == localListPtr->rbTreeNil) {
		/* Tree is empty */
		*firstItemPtr = (LbSortItemPtr)(localListPtr->rbTreeNil);
		*index = 0;
		*dataPtr = NULL;
		found = false;
	}
	else {
		curPtr = localListPtr->rbTreePtr;
		
		/* Proceed down the left branches */
		while (curPtr->leftBranch != localListPtr->rbTreeNil) {
			curPtr = curPtr->leftBranch;
		}
		
		/* Return the found item */
		*firstItemPtr = (LbSortItemPtr) curPtr;
		*index = curPtr->index;
		*dataPtr = curPtr->dataPtr;
		found = true;
	}
	
	return (found);
}

/*********************************************************************************
*
*	Name:			LbSortLast
*
*	Summary:		Return the last sorted list item, index data, and data pointer. 
*
*	Arguments:
*		LbSortListPtr		theListPtr		- The relevant sorted list.
*		LbSortItemPtr		*lastItemPtr	- Pointer to the last sorted list item.
*		LbUsFourByte		*index			- Last index data.
*		void				**dataPtr		- Last data pointer.
*
*	Function return: True if successful.
*
*********************************************************************************/
Boolean LbSortLast(LbSortListPtr theListPtr, 
					LbSortItemPtr *lastItemPtr, 
					LbUsFourByte *index, void **dataPtr)
{
	Boolean			found;			/* Item was found; function return value */
	LbSortedListPtr	localListPtr;	/* theListPtr as LbSortedListPtr */
	RBTreeNodePtr	curPtr;			/* Pointer to current sort node */
	
	
	localListPtr = (LbSortedListPtr) theListPtr;
	
	if (localListPtr->rbTreePtr == localListPtr->rbTreeNil) {
		/* Tree is empty */
		*lastItemPtr = (LbSortItemPtr)(localListPtr->rbTreeNil);
		*index = 0;
		*dataPtr = NULL;
		found = false;
	}
	else {
		curPtr = localListPtr->rbTreePtr;
		
		/* Proceed down the right branches */
		while (curPtr->rightBranch != localListPtr->rbTreeNil) {
			curPtr = curPtr->rightBranch;
		}
		
		/* Return the found item */
		*lastItemPtr = (LbSortItemPtr) curPtr;
		*index = curPtr->index;
		*dataPtr = curPtr->dataPtr;
		found = true;
	}
	
	return (found);
}

/*********************************************************************************
*
*	Name:			LbSortNext
*
*	Summary:		Find the successor item to startItemPtr in the sorted list.
*					Return the sorted list item, index data, and data pointer. 
*
*	Arguments:
*		LbSortListPtr		theListPtr		- The relevant sorted list.
*		LbSortItemPtr		startItemPtr	- Pointer to the start sorted list item.
*		LbSortItemPtr		*nextItemPtr	- Pointer to the next sorted list item.
*		LbUsFourByte		*index			- Found index data.
*		void				**dataPtr		- Found data pointer.
*
*	Function return: True if successful.
*
*********************************************************************************/
Boolean LbSortNext(LbSortListPtr theListPtr, 
					LbSortItemPtr startItemPtr, LbSortItemPtr *nextItemPtr, 
					LbUsFourByte *index, void **dataPtr)
{
	Boolean			found = false;	/* Item was found; function return value */
	LbSortedListPtr	localListPtr;	/* theListPtr as LbSortedListPtr */
	RBTreeNodePtr	startNodePtr;	/* RBTreeNodePtr version of startItemPtr */
	RBTreeNodePtr	curPtr;			/* Pointer to current sort node */
	RBTreeNodePtr	parPtr;			/* Pointer to parent of curPtr sort node */
	
	
	localListPtr = (LbSortedListPtr) theListPtr;
	startNodePtr = (RBTreeNodePtr) startItemPtr;
	if (startNodePtr->rightBranch != localListPtr->rbTreeNil) {
		/* Next item is minimum of right subtree */
		curPtr = startNodePtr->rightBranch;
		while (curPtr->leftBranch != localListPtr->rbTreeNil) {
			curPtr = curPtr->leftBranch;
		}
		
		found = true;
		*nextItemPtr = (LbSortItemPtr) curPtr;
		*index = curPtr->index;
		*dataPtr = curPtr->dataPtr;
	}
	else {
		/* Next item is first ancestor whose left branch contains minNode */
		curPtr = startNodePtr;
		parPtr = startNodePtr->parent;
		while ((parPtr != localListPtr->rbTreeNil) && 
				(curPtr == parPtr->rightBranch)) {
			curPtr = parPtr;
			parPtr = parPtr->parent;
		}
		
		if (parPtr == localListPtr->rbTreeNil) {
			/* There is no successor node to minNode in the tree */
			found = false;
			*nextItemPtr = (LbSortItemPtr)(localListPtr->rbTreeNil);
			*index = 0;
			*dataPtr = NULL;
		}
		else {
			found = true;
			*nextItemPtr = (LbSortItemPtr) parPtr;
			*index = parPtr->index;
			*dataPtr = parPtr->dataPtr;
		}
	}
	
	return (found);
}

/*********************************************************************************
*
*	Name:			LbSortPrev
*
*	Summary:		Find the predecessor item to startItemPtr in the sorted list.
*					Return the sorted list item, index data, and data pointer. 
*
*	Arguments:
*		LbSortListPtr		theListPtr		- The relevant sorted list.
*		LbSortItemPtr		startItemPtr	- Pointer to the start sorted list item.
*		LbSortItemPtr		*prevItemPtr	- Pointer to the previous sorted list item.
*		LbUsFourByte		*index			- Found index data.
*		void				**dataPtr		- Found data pointer.
*
*	Function return: True if successful.
*
*********************************************************************************/
Boolean LbSortPrev(LbSortListPtr theListPtr, 
					LbSortItemPtr startItemPtr, LbSortItemPtr *prevItemPtr, 
					LbUsFourByte *index, void **dataPtr)
{
	Boolean			found = false;	/* Item was found; function return value */
	LbSortedListPtr	localListPtr;	/* theListPtr as LbSortedListPtr */
	RBTreeNodePtr	startNodePtr;	/* RBTreeNodePtr version of startItemPtr */
	RBTreeNodePtr	curPtr;			/* Pointer to current sort node */
	RBTreeNodePtr	parPtr;			/* Pointer to parent of curPtr sort node */
	
	
	localListPtr = (LbSortedListPtr) theListPtr;
	startNodePtr = (RBTreeNodePtr) startItemPtr;
	if (startNodePtr->leftBranch != localListPtr->rbTreeNil) {
		/* Prev item is maximum of left subtree */
		curPtr = startNodePtr->leftBranch;
		while (curPtr->rightBranch != localListPtr->rbTreeNil) {
			curPtr = curPtr->rightBranch;
		}
		
		found = true;
		*prevItemPtr = (LbSortItemPtr) curPtr;
		*index = curPtr->index;
		*dataPtr = curPtr->dataPtr;
	}
	else {
		/* Prev item is first ancestor whose right branch contains minNode */
		curPtr = startNodePtr;
		parPtr = startNodePtr->parent;
		while ((parPtr != localListPtr->rbTreeNil) && 
				(curPtr == parPtr->leftBranch)) {
			curPtr = parPtr;
			parPtr = parPtr->parent;
		}
		
		if (parPtr == localListPtr->rbTreeNil) {
			/* There is no predecessor node to minNode in the tree */
			found = false;
			*prevItemPtr = (LbSortItemPtr)(localListPtr->rbTreeNil);
			*index = 0;
			*dataPtr = NULL;
		}
		else {
			found = true;
			*prevItemPtr = (LbSortItemPtr) parPtr;
			*index = parPtr->index;
			*dataPtr = parPtr->dataPtr;
		}
	}
	
	return (found);
}

/*********************************************************************************
*
*	Name:			LbSortFind
*
*	Summary:		Find the last list item in the sorted list <= minSortKey.
*					Return the sorted list item.
*
*	Arguments:
*		LbSortListPtr		theListPtr		- The relevant sorted list.
*		LbSortKeyPtr		minSortKeyPtr	- The minimum sort key pointer.
*		LbSortItemPtr		*sortItemPtr	- Pointer to the sorted list item.
*
*	Function return: True if successful.
*
*********************************************************************************/
Boolean LbSortFind(LbSortListPtr theListPtr, 
					LbSortKeyPtr minSortKeyPtr, LbSortItemPtr *sortItemPtr)
{
	Boolean			found = false;	/* Data was found; function return value */
	LbSortedListPtr	localListPtr;	/* theListPtr as LbSortedListPtr */
	LbSortKeyComparator
					*CompareKeys;	/* Sort key comparator */
	RBTreeNodePtr	curPtr;			/* Pointer to current sort node */
	RBTreeNodePtr	parPtr;			/* Pointer to parent of curPtr sort node */
	RBTreeNodePtr	nextPtr;		/* Pointer to successor of curPtr sort node */
	LbUsFourByte	index;			/* Ignored */
	void			*dataPtr;		/* Ignored */
	
	
	localListPtr = (LbSortedListPtr) theListPtr;
	CompareKeys = localListPtr->keyComparatorPtr;
	
	/* First, find any largest node <= minSortKey */
	parPtr = localListPtr->rbTreeNil;
	curPtr = localListPtr->rbTreePtr;
	while ((curPtr != localListPtr->rbTreeNil) && 
			((*CompareKeys)(minSortKeyPtr, &(curPtr->sortKey)) != 0)) {
		parPtr = curPtr;
		if ((*CompareKeys)(minSortKeyPtr, &(curPtr->sortKey)) == -1) {
			/* Left branch */
			curPtr = curPtr->leftBranch;
		}
		else {
			/* Right branch */
			curPtr = curPtr->rightBranch;
		}
	}
	if (curPtr == localListPtr->rbTreeNil) {
		/* Ran out of nodes (minimum sort key isn't in sort tree), 
			so parent is desired node */
		curPtr = parPtr;
	}
	
	if (curPtr == localListPtr->rbTreeNil) {
		/* Sort tree was empty */
		*sortItemPtr = (LbSortItemPtr)(localListPtr->rbTreeNil);
		found = false;
	}
	else {
		/* As long as the node key equals minSortKey, get its successor */
		nextPtr = curPtr;
		while ((*CompareKeys)(&(nextPtr->sortKey), minSortKeyPtr) == 0) {
			curPtr = nextPtr;
			if (! LbSortNext(theListPtr, 
						(LbSortItemPtr)curPtr, (LbSortItemPtr*)(&nextPtr), 
						&index, &dataPtr)) {
				/* There is no successor sort node to this one in the tree */
				break;
			}
		}
		*sortItemPtr = (LbSortItemPtr) curPtr;
		found = true;
	}
	
	return (found);
}

/*********************************************************************************
*
*	Name:			LbSortDisposeList
*
*	Summary:		Dispose of an existing sorted list.
*
*	Arguments:		
*		LbSortListPtr		theListPtrPtr	- Address of the sorted list pointer.
*
*	Function return: None.
*
*********************************************************************************/
void LbSortDisposeList(LbSortListPtr *theListPtrPtr)
{
	LbSortedListPtr		theListPtr;				/* The sorted list pointer */
	
	
	if (theListPtrPtr != NULL) {
		theListPtr = (LbSortedListPtr)(*theListPtrPtr);
		if (theListPtr != NULL) {
			if (theListPtr->rbTreeStartPtr != NULL) {
				/* Free the list memory */
				LbMmFree((void **) &(theListPtr->rbTreeStartPtr));
			}
			
			/* Free the list itself */
			LbMmFree((void **) theListPtrPtr);
		}
	}
}
