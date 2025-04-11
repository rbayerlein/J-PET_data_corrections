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
*     Module Name:        LbInterface.h
*     Revision Number:    1.3
*     Date last revised:  1 October 2012
*     Programmer:         Steven Vannoy
*     Date Originated:    Tuesday, July 21 81992
*
*     Module Overview:    This module declars all constants, types, variables
*						and prototypes for the Debug library.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*     Global variables defined:   
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
*			Revision date:		12 December 2011
*
*			Revision description:	Changed functions to support 4 choices.
*
*********************************************************************************/

#ifndef LIB_INTERFACE
#define LIB_INTERFACE

#include <stdio.h>

#include "LbTypes.h"


/* CONSTANTS */
#define	LBINYes				1	
#define LBINNo				2	
#define LBINCancel			3	
#define LBINQuit			4	
#define LBINUserChoice		5
#define LBINMaxErrorPrompt	60

/* FLAGS */
#define	LBINFg_Numeric		LBFlag0		/* Input should be numeric */
#define LBINFg_Unique		LBFlag1		/* Input should not be sub-string matched */

/* TYPES */
typedef Boolean (*LbInFilterFnTy)(char *inputString, Boolean *isValid,
			char *errorPrompt, LbUserData userData);

/* PROTOTYPES */
LbTwoByte	LbInAskNew(LbUsFourByte flags, char *prompt, LbTwoByte defaultChoice,
				Boolean restrictToChoices, Boolean *canceled, char *choice1,
				char *choice2, char *choice3, char *choice4,
				char *response);
LbTwoByte	LbInAsk(char *prompt, LbTwoByte defaultChoice, Boolean restrictToChoices,
				Boolean *canceled, char *choice1, char *choice2, char *choice3, char *choice4,
				char *response);
double		LbInAskDouble(char *prompt, LbTwoByte defaultChoice, Boolean isSigned,
				LbTwoByte decimalPrec,
				Boolean restrictToChoices,
				LbTwoByte numChoices, Boolean *canceled, 
				double choice1, double choice2, double choice3, double choice4);
double		LbInAskDoubleInRange(char *prompt, LbTwoByte decimalPrec,
				Boolean isSigned, Boolean *canceled, double min, double max);
LbTwoByte	LbInAskFiltered(LbUsFourByte flags, char *prompt, LbTwoByte defaultChoice,
				Boolean *canceled, char *choice1, char *choice2, char *choice3, char *choice4,
				char *response, LbInFilterFnTy filterFunction, Boolean *userResult,
				LbUserData userData);
double		LbInAskFloat(char *prompt, LbTwoByte defaultChoice, Boolean isSigned,
				LbTwoByte decimalPrec,
				Boolean restrictToChoices,
				LbTwoByte numChoices, Boolean *canceled, 
				double choice1, double choice2, double choice3, double choice4);
LbFourByte	LbInAskFourByte(char *prompt, LbTwoByte defaultChoice, Boolean isSigned,
				Boolean restrictToChoices, LbTwoByte numChoices, Boolean *canceled, 
				LbFourByte choice1, LbFourByte choice2, LbFourByte choice3, LbFourByte choice4);
LbFourByte	LbInAskFourByteInRange(char *prompt,
				Boolean isSigned, Boolean *canceled, LbFourByte min, LbFourByte max);
LbTwoByte	LbInAskTwoByte(char *prompt, LbTwoByte defaultChoice, Boolean isSigned,
				Boolean restrictToChoices, LbTwoByte numChoices, Boolean *canceled, 
				LbTwoByte choice1, LbTwoByte choice2, LbTwoByte choice3, LbTwoByte choice4);
LbTwoByte	LbInAskYesNo(char *prompt, LbTwoByte defaultChoice);
LbTwoByte	LbInAskYesNoCancel(char *prompt, LbTwoByte defaultChoice);
LbTwoByte	LbInAskYesNoQuit(char *prompt, LbTwoByte defaultChoice);
char		*LbInGetInput(char *inputBuffer, int bufferLength);
Boolean		LbInMatchStr(char *strToCheck, char *strToMatch);
LbFourByte	LbInPrintf(char *fmt, ...);
#endif /* LIB_INTERFACE */
