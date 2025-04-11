/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992 Department of Radiology						*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name:        LbMacros.c
*     Revision Number:    1.0
*     Date last revised:  Wednesday, July 22, 1992
*     Programmer:         Steven Vannoy
*     Date Originated:    Tuesday, July 21, 1992
*
*     Module Overview:	This module provides flag manipulation macros.					
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
#ifndef LB_MACROS_HDR
#define LB_MACROS_HDR

/**********************
*	LbFgSet
*
*	Arguments:	flag		- Flag variable to set.
*				bitField	- bitField to set in variable.
*
*	Purpose:	Set the specified bit in the flag variable.
*
*	Result:	None.
***********************/
#define	LbFgSet(flag, bitField)	((flag) |= (bitField))

/**********************
*	LbFgClear
*
*	Arguments:	flag		- Flag variable to set.
*				bitField	- Bit field to clear in variable.
*
*	Purpose:	Clear the specified bit field in the flag variable.
*
*	Result:	None.
***********************/
#define	LbFgClear(flag, bitField)	((flag) &= ~((bitField)))

/**********************
*	LbFgIsSet
*
*	Arguments:	flag		- Flag variable to set.
*				bitField	- Bit field to check.
*
*	Purpose:	Check the flag variable for a set bit.
*
*	Result:	true if any flag in the bitField is set.
***********************/
#define	LbFgIsSet(flag, bitField)	((((flag) & (bitField)) != 0) ? true : false)

#endif /* LB_MACROS_HDR */
