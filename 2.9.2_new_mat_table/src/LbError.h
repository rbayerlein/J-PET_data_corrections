/*********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992-2012 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/*********************************************************************************
*
*     Module Name:        LbError.h
*     Revision Number:    1.1
*     Date last revised:  1 October 2012
*     Programmer:         Steven Vannoy
*     Date Originated:    Tuesday, July 21, 1992
*
*     Module Overview:	This module provides definitions for the error library.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*     	LbEnGetOptions
*
*     Global variables defined:   
*
*********************************************************************************/
#ifndef LIB_ERROR
#define LIB_ERROR

#include "LbTypes.h"


/* FLAGS */
#define	ERFgEcho		LBFlag0

/* CONSTANTS */
#define	ERMgCdNull		0
#define ERMgCdFile		1
#define	ERMgCdGeneric	2
#define ERMgCdHeader	3

#define	ERMxErStrLen	255

/* TYPES */
typedef	LbUsFourByte	ErMgCdTy;

typedef enum{
	ERErCdNull,
	ERErCdGeneric,
	ERErCdOutOfMem,
	ERErCdUserCancel,
	
	/* File Manager */
	ERErCdFlNotFound,
	ERErCdFlNtOpBlkWrite,
	ERErCdFlAtGrPastLim,
	ERErCdFlBadBuff,
	ERErCdFlAccessDenied,
	ERErCdFlPermanent,
	ERErCdFlEndOfFile,
	
	/* Environment Manager */
	ERErCdEnSwNotUnique,
	ERErCdEnSwNotFound,
	
	/* List Mode Manager */
	ERErCdLmBadClumpNum,
	ERErCdLmBadFileType,
	
	/* Header Manager */
	ERErCdHdElemNotFound
} ErErCdTy;

/*	MACROS */
/**********************
*	ErStGeneric
*
*	Arguments:	
*				char	*errorString	- The error string.
*
*	Purpose:	Set primitive error with generic codes.
*
*	Result:	None.
***********************/
#define ErStGeneric(errorString)	ErStPrimError(ERMgCdGeneric, ERErCdGeneric, (errorString))

/**********************
*	ErStCancel
*
*	Arguments:	
*				char	*errorString	- The error string.
*
*	Purpose:	Set primitive error with generic manager and cancel code.
*
*	Result:	None.
***********************/
#define ErStCancel(errorString)	ErStPrimError(ERMgCdGeneric, ERErCdUserCancel, (errorString))


/* PROTOTYPES */
void	ErAbort(char *errorString);
void	ErBreak(char *errorString);
void	ErAlert(char *errorString, Boolean getResponse);
void	ErIgnoreBegin(void);
void	ErIgnoreEnd(void);
Boolean	ErInit(LbUsFourByte initFlags, FILE *outputFile);
void	ErClear(void);
void	ErClearIf(ErMgCdTy managerCode, ErErCdTy errorCode,
			Boolean *cleared);
void	ErClearIfManager(ErMgCdTy managerCode, Boolean *cleared);
void	ErHandle(char *errorString, Boolean getResponse);
void	ErNewMgCode(ErMgCdTy *managerCode);
void	ErWhatError(Boolean *isError, ErMgCdTy  *managerCode,
			ErErCdTy *errorCode, char *errorStr);
Boolean	ErIsInError(void);
void	ErStFileError(char *errorString);
void	ErStPrimError(ErMgCdTy manager, ErErCdTy errorCode, char *errorString);
void	ErTerminate(void);
#endif /* LIB_ERROR */
