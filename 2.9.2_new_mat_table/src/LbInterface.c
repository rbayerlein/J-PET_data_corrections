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
*     Module Name:        LbInterface.c
*     Revision Number:    1.7
*     Date last revised:  23 July 2013
*     Programmer:         Steven Vannoy
*     Date Originated:    Tuesday, July 21, 1992
*
*     Module Overview:    This module provides routines for performing user interface
*						tasks.
*						MACHINE DEPENDANT						
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*		LbInAsk
*		LbInAskDouble
*		LbInAskDoubleInRange
*		LbInAskFiltered
*		LbInAskFloat
*		LbInAskFourByte
*		LbInAskFourByteInRange
*		LbInAskTwoByte
*		LbInAskYesNo
*		LbInAskYesNoCancel
*		LbInAskYesNoQuit
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
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:	removed calls to gets, replaced with
*							calls to lbInGetS.
*
*********************************************************************************/

#include	<stdio.h>
#include	<string.h>
#include	<ctype.h>
#include	<stdarg.h>

#include "SystemDependent.h"

/*	CONSTANTS */
#define LBIN_MAXLINE	65


#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbInterface.h"

/*	LOCAL TYPES */
typedef struct {
	Boolean	RestrictToChoices;
} lbInAskFilterTy;

typedef struct {
	Boolean	RestrictToChoices;
	Boolean	IsSigned;
	double	Min;
	double	Max;
} lbInFloatFilterTy;

typedef struct {
	Boolean	RestrictToChoices;
	Boolean	IsSigned;
	LbFourByte	Min;
	LbFourByte	Max;
} lbInFourFilterTy;

/*  LOCAL GLOBALS */
char		*lbInYesStr = "Yes";			/* Yes string */
char		*lbInNoStr = "No";				/* No string */
char		*lbInCancelStr = "Cancel";		/* Cancelation string */
char		*lbInQuitStr = "Quit";			/* Quit string */
char		lbInPrintfString[2048];			/* Buffer for LbInPrintf emulation */

/*	LOCAL MACROS */
#define		lbInTwoFilter	lbInFourFilter
#define		lbInTwoFilterTy	lbInFourFilterTy
/**********************
*	LBINFg_IsNumeric
*
*	Arguments:	flags	- Our modification flag.
*
*	Purpose:	See if we are dealing with a numeric input.
*
*	Result:	TRUE if numeric flag set.
***********************/
#define LBINFg_IsNumeric(flags)	LbFgIsSet((flags), LBINFg_Numeric)

/**********************
*	LBINFg_IsUnique
*
*	Arguments:	flags	- Our modification flag.
*
*	Purpose:	See if we should not do sub-string matching.
*
*	Result:	TRUE if numeric flag set.
***********************/
#define LBINFg_IsUnique(flags)	LbFgIsSet((flags), LBINFg_Unique)

/*	LOCAL FUNCTIONS */
void	lbInAcknowledgeError(char *errorPrompt);
Boolean	lbInAskFilter(char *inputString, Boolean *isValid, char *errorPrompt,
			LbUserData userData);
Boolean	lbInDoubleFilter(char *inputString, Boolean *isValid, char *errorPrompt,
			LbUserData userData);
Boolean	lbInFloatFilter(char *inputString, Boolean *isValid, char *errorPrompt,
			LbUserData userData);
Boolean	lbInDblRngFilter(char *inputString, Boolean *isValid, char *errorPrompt,
			LbUserData userData);
Boolean	lbInFourFilter(char *inputString, Boolean *isValid, char *errorPrompt,
			LbUserData userData);
Boolean	lbInFourRngFilter(char *inputString, Boolean *isValid,
			char *errorPrompt, LbUserData userData);
void	lbInBuildChoiceStr(char *choiceStr, LbTwoByte defaultChoice,
			char *choice1, char *choice2, char *choice3, char *choice4);
void	lbInPrintQuestion(char *question, LbTwoByte defaultChoice,
			char *choice1, char *choice2, char *choice3, char *choice4);
char	*lbInGetS( char *str, int size );

#ifdef COCOA
/* Utility function for C defined in GUI code */
/* NOTE:  Must match actual definition */
extern void WriteStrCocoa( char* theString );
#endif


/*	LOCAL MACROS */

/*	FUNCTIONS */
/*****************************************
*		lbInAcknowledgeError
*
*	Arguments:
*				char		*errorPrompt	- Prompt to display.
*
*	Purpose:	Present the user with a prompt indicating they have entered a bogus value.
*				Wait for them to acknowledge the error.
*
*	Returns:	Nothing.
*
******************************************/
void	lbInAcknowledgeError(char *errorPrompt)
{
	char responseStr[32];	/* Response to New Line, 32 is for safety */
	
	/* Print the error string */
	LbInPrintf("%s [New Line] to try again:", errorPrompt);
	
	/* Get the response */
	lbInGetS(responseStr, 32);
}

/*****************************************
*		LbInAsk
*
*	Arguments:
*				char		*prompt				- Prompt to display.
*				LbTwoByte	defaultChoice		- Number of choice that is the defaultChoice.
*				Boolean		restrictToChoices	- Flag to restrict input to given choices.
*				Boolean		*canceled			- Set true if user cancels input.
*				char		*choice1			- First choice given.
*				char		*choice2			- Second choice given.
*				char		*choice3			- Third choice given.
*				char		*choice4			- Fourth choice given.
*				char		*response			- User's input if not given choice.
*
*	Purpose:	Present the user with a prompt and get their response.
*
*	Returns:	The number of the choice selected or LBINUserChoice if
*				user enters value that is not a given choice.
*
******************************************/
LbTwoByte	LbInAsk(char *prompt, LbTwoByte defaultChoice, Boolean restrictToChoices,
				Boolean *canceled, char *choice1, char *choice2, char *choice3, char *choice4,
				char *response)
{
	/* Call new Ask routine with zero for flags */
	return(LbInAskNew(0, prompt, defaultChoice, restrictToChoices, canceled, choice1,
		choice2, choice3, choice4, response));
}
/*****************************************
*		LbInAskNew
*
*	Arguments:
*				LbUsFourByte	flags			- Behavior modification flags.
*				char			*prompt				- Prompt to display.
*				LbTwoByte		defaultChoice		- Number of choice that is the defaultChoice.
*				Boolean			restrictToChoices	- Flag to restrict input to given choices.
*				Boolean			*canceled			- Set true if user cancels input.
*				char			*choice1			- First choice given.
*				char			*choice2			- Second choice given.
*				char			*choice3			- Third choice given.
*				char			*choice4			- Fourth choice given.
*				char			*response			- User's input if not given choice.
*
*	Purpose:	Present the user with a prompt and get their response.
*
*	Returns:	The number of the choice selected or LBINUserChoice if
*				user enters value that is not a given choice.
*
******************************************/
LbTwoByte	LbInAskNew(LbUsFourByte flags, char *prompt, LbTwoByte defaultChoice,
				Boolean restrictToChoices, Boolean *canceled, char *choice1,
				char *choice2, char *choice3, char *choice4,
				char *response)
{
	static char		dummyResponseBuff[100];	/* Dummy response buffer */
	lbInAskFilterTy	filterInfo;				/* Our filter's info */
	LbTwoByte		result = -1;			/* Calculated result */
	
	#ifdef LB_DEBUG
		/* Verify they gave choices if they restrict to choices */
		if ((restrictToChoices == true) && (choice2 == 0)
				&& (choice3 == 0) && (choice4 == 0)) {
			ErAbort("Why restrict to choices without supplying them to LbInAskFiltered");
		}
			
	#endif
	
	/* Initialize filter info */
	filterInfo.RestrictToChoices = restrictToChoices;
	
	/* If restricting to choices, pass dummy response buffer */
	if ((restrictToChoices == true) && (response == 0))
		response = dummyResponseBuff;
	
	/* Use generalized filter ask */
	result = LbInAskFiltered(flags, prompt, defaultChoice, canceled, 
		choice1, choice2, choice3, choice4,
		response, lbInAskFilter, 0, (LbUserData) &filterInfo);
		
	return (result);
}

/*****************************************
*		LbInAskDouble
*
*	Arguments:
*				char		*prompt				- Prompt to display.
*				LbTwoByte	defaultChoice		- Number of choice that is the defaultChoice.
*				Boolean		isSigned			- Is input isSigned value?
*				LbTwoByte	decimalPrec			- How many digits of decimal precision?
*				Boolean		restrictToChoices	- Flag to restrict input to given choices.
*				LbTwoByte	numChoices			- The number of choices being supplied.
*				Boolean		*canceled			- Set true if user cancels input.
*				double		choice1				- First choice given.
*				double		choice2				- Second choice given.
*				double		choice3				- Third choice given.
*				double		choice4				- Fourth choice given.
*
*	Purpose:	Present the user with a prompt and get their response.
*
*	Returns:	The number entered if user did not cancel.
*
******************************************/
double	LbInAskDouble(char *prompt, LbTwoByte defaultChoice, Boolean isSigned,
				LbTwoByte decimalPrec,
				Boolean restrictToChoices,
				LbTwoByte numChoices, Boolean *canceled, 
				double choice1, double choice2, double choice3, double choice4)
{
	LbTwoByte			userResult;			/* The result that comes back from LbInAsk */
	char				formatStr[8];		/* Format string for conversion */
	char				choice1Buff[128];	/* String buffer of choice1 */
	char				choice2Buff[128];	/* String buffer of choice2 */
	char				choice3Buff[128];	/* String buffer of choice3 */
	char				choice4Buff[128];	/* String buffer of choice4 */
	char				*choice1Ptr;		/* String conversion of choice1 */
	char				*choice2Ptr;		/* String conversion of choice2 */
	char				*choice3Ptr;		/* String conversion of choice3 */
	char				*choice4Ptr;		/* String conversion of choice4 */
	char				userStr[128];		/* String result for user choice */
	double 				result;				/* Float conversion of string input */
	lbInFloatFilterTy	filterInfo;			/* Our filter info */
	
	#ifdef LB_DEBUG
		/* Verify they gave choices if they restrict to choices */
		if ((restrictToChoices == true) && (numChoices <= 1)) {
			ErAbort("Why restrict to choices without supplying them to LbInAskFiltered");
		}
	#endif

	/* Setup format string */
	sprintf(formatStr, "%%.%dg", decimalPrec);

	/* Convert choices to strings */
	{
		/* Setup choice strings */
		if (numChoices >= 1) {
			choice1Ptr = choice1Buff;
			sprintf(choice1Ptr, formatStr, choice1);
		}
		else
			choice1Ptr = 0;
		
		if (numChoices >= 2) {
			choice2Ptr = choice2Buff;
			sprintf(choice2Ptr, formatStr, choice2);
		}
		else
			choice2Ptr = 0;
			
		if (numChoices >= 3) {
			choice3Ptr = choice3Buff;
			sprintf(choice3Ptr, formatStr, choice3);
		}
		else
			choice3Ptr = 0;
			
		if (numChoices >= 4) {
			choice4Ptr = choice4Buff;
			sprintf(choice4Ptr, formatStr, choice4);
		}
		else
			choice4Ptr = 0;
	}
	
	/* Initialize filter information */
	filterInfo.RestrictToChoices = restrictToChoices;
	filterInfo.IsSigned = isSigned;
	filterInfo.Min = filterInfo.Max = 0;
	
	/* Use generalized filter ask */
	userResult = LbInAskFiltered(LBINFg_Numeric, prompt, defaultChoice, canceled, 
		choice1Ptr, choice2Ptr, choice3Ptr, choice4Ptr,
		userStr, lbInDoubleFilter, 0, (LbUserData) &filterInfo);

	/* Check results */
	if ((canceled == 0) || *canceled == false)
	{
		if (userResult == 1)
			result = choice1;
		else if (userResult == 2)
			result = choice2;
		else if (userResult == 3)
			result = choice3;
		else if (userResult == 4)
			result = choice4;
		else if (userResult == LBINUserChoice)
			result =  atof(userStr);
		#ifdef LB_DEBUG
			else
				ErAbort("Unexepected result in LbInAskDouble");
		#endif
	}
	
	return (result);
}

/*****************************************
*		LbInAskDoubleInRange
*
*	Arguments:
*				char		*prompt				- Prompt to display.
*				Boolean		isSigned			- Is value signed.
*				LbTwoByte	decimalPrec			- The significance.
*				Boolean		*canceled			- Set true if user cancels input.
*				double		min					- Minimum range.
*				double		max					- Maximum range.
*				double		defaultChoice		- Default choice.
*
*	Purpose:	Present the user with a prompt and get their response.
*
*	Returns:	The value entered.
*
******************************************/
double	LbInAskDoubleInRange(char *prompt, LbTwoByte decimalPrec,
				Boolean isSigned, Boolean *canceled, double min, double max)
{
	char				formatStr[8];		/* Format for decimal precision */
	char				minStrBuff[32];		/* String buffer of min choice */
	char				maxStrBuff[32];		/* String buffer of max choice */
	char				userStr[32];		/* String result for user choice */
	double				result;				/* Calculated result */
	lbInFloatFilterTy	filterInfo;			/* Our filter's info */
	LbTwoByte			userResult;			/* The result that comes back from LbInAskFiltered */
	
	
	/* Setup format string */
	sprintf(formatStr, "%%.%df", decimalPrec);

	/* Convert choices to strings */
	{
		/* Setup choice strings */
		sprintf(minStrBuff, formatStr, min);
		sprintf(maxStrBuff, formatStr, max);
	}
	
	/* Initialize filter info */
	filterInfo.Min = min;
	filterInfo.Max = max;
	filterInfo.IsSigned = isSigned;
	filterInfo.RestrictToChoices = false;
	
	/* Use generalized filter ask */
	userResult = LbInAskFiltered(LBINFg_Numeric, prompt, 0, canceled, minStrBuff, maxStrBuff,
		0, 0, userStr, lbInDblRngFilter, 0, (LbUserData) &filterInfo);

	/* Check results */
	if ((canceled == 0) || *canceled == false)
	{
		if (userResult == 1)
			result = min;
		else if (userResult == 2)
			result = max;
		else if (userResult == LBINUserChoice)
			result = (double) atof(userStr);
		#ifdef LB_DEBUG
			else
				ErAbort("Unexepected result in LbInAskDoubleInRange");
		#endif
	}
	
	return (result);
}

/*****************************************
*		lbInDblRngFilter
*
*	Arguments:
*				char			*inputString		- The input string.
*				Boolean			*isValid			- Result of filter.
*				char			*errorPrompt		- String to chastise user.
*				LbUserData		userData			- Our user data.
*
*	Purpose:	Filter input for double, in given range.
*
*	Returns:	True unless an error occurs.
*
******************************************/
Boolean	lbInDblRngFilter(char *inputString, Boolean *isValid, char *errorPrompt,
			LbUserData userData)
{
	Boolean				okay = false;
	double				testValue;
	lbInFloatFilterTy	*filterInfo = (lbInFloatFilterTy *) userData;
	
	do	{ /* Process Loop */
	
		/* First make sure it is a valid foat */
		if (!lbInFloatFilter(inputString, isValid, errorPrompt, userData))
			break;
			
		/* If still valid here, check the range */
		if (*isValid == true) {
		
			/* Convert string to value */
			testValue = atof(inputString);
			
			if ((testValue < filterInfo->Min) || (testValue > filterInfo->Max)) {
			
				/* Value is invalid */
				strcpy(errorPrompt, "Input is out of range");
				*isValid = false;
			}
		}
		
		okay = true;
	} while (false);

	return (okay);
}

/*****************************************
*		lbInDoubleFilter
*
*	Arguments:
*				char			*inputString		- The input string.
*				Boolean			*isValid			- Result of filter.
*				char			*errorPrompt		- String to chastise user.
*				LbUserData		userData	- Our user data.
*
*	Purpose:	Filter input for doubles.
*
*	Returns:	True unless an error occurs.
*
******************************************/
Boolean	lbInDoubleFilter(char *inputString, Boolean *isValid, char *errorPrompt,
			LbUserData userData)
{
	Boolean				okay = false;
	lbInFloatFilterTy	*filterInfo = (lbInFloatFilterTy *) userData;
	LbTwoByte			parseIndex;
	
	do	{ /* Process Loop */
	
		/* If not restricting to choices verify number is valid */
		if (filterInfo->RestrictToChoices == false) {
			
			/* Assume value is good */
			*isValid = true;
			
			/* If not signed, and number contains sing value, its bogus */
			if ((filterInfo->IsSigned == false) && (inputString[0] == '-')){
				strcpy(errorPrompt, "Negative values not allowed");
				*isValid = false;
			}
			else {
				/* Make sure it is not alphabetic */
				parseIndex = 0;
				while (inputString[parseIndex] != 0) {
					if ((isdigit(inputString[parseIndex]) == false)
							&& (inputString[parseIndex] != '.')
							&& (inputString[parseIndex] != '-')) {
							
						/* Value is invalid */
						strcpy(errorPrompt, "Input must be numeric");
						*isValid = false;
						break;
					}
					
					parseIndex++;
				}
			}
		}
		else {
			*isValid = false;
			strcpy(errorPrompt, "You must answer");
		}
				
		okay = true;
	} while (false);

	return (okay);
}

/*****************************************
*		lbInAskFilter
*
*	Arguments:
*				char			*inputString		- The input string.
*				Boolean			*isValid			- Result of filter.
*				char			*errorPrompt		- String to chastise user.
*				LbUserData		userData	- Our user data.
*
*	Purpose:	Filter input for general ask function.
*
*	Returns:	True unless an error occurs.
*
******************************************/
Boolean	lbInAskFilter(char *inputString, Boolean *isValid, char *errorPrompt,
			LbUserData userData)
{
	Boolean	okay = false;
	lbInAskFilterTy	*filterInfo = (lbInAskFilterTy *) userData;
	
	do	{ /* Process Loop */
	
		/* If not restricting to choices its okay, otherewise it is not */
		if ((filterInfo->RestrictToChoices == false) ||
			((filterInfo->RestrictToChoices == true) && (strlen(inputString) == 0))) {
			
			*isValid = true;
		}
		else {
			*isValid = false;
			strcpy(errorPrompt, "You must answer ");
		}
		
		okay = true;
	} while (false);

	return (okay);
}

/*****************************************
*		LbInAskFiltered
*
*	Arguments:
*				LbUsFourByte	flags				- Modification flags.
*				char			*prompt				- Prompt to display.
*				LbTwoByte		defaultChoice		- Number of choice that is the defaultChoice.
*				Boolean			*canceled			- Set true if user cancels input.
*				char			*choice1			- First choice given.
*				char			*choice2			- Second choice given.
*				char			*choice3			- Third choice given.
*				char			*choice4			- Fourth choice given.
*				char			*response			- User's input if not given choice.
*				LbInFilterFnTy	filterFunction		- The filter function.
*				Boolean			*userResult			- Result of filter function.
*				LbUserData		userData			- The user's data.
*
*	Purpose:	Present the user with a prompt and get their response.
*				Allow caller to filter input values.
*
*	Returns:	The number of the choice selected or LBINUserChoice if
*				user enters value that is not a given choice.
*
******************************************/
LbTwoByte	LbInAskFiltered(LbUsFourByte flags,
				char *prompt, LbTwoByte defaultChoice,
				Boolean *canceled, char *choice1, char *choice2, char *choice3, char *choice4,
				char *response, LbInFilterFnTy filterFunction, Boolean *userResult,
				LbUserData userData)
{
	Boolean	  	finished = false;					/* Finished flag */
	Boolean		goofed = false;						/* Goofed flag */
	Boolean		isValid;							/* Filter evaluation */
	Boolean		filterResult;						/* Result of their filter */
	char		lbInInputBuffer[255];				/* Input buffer */
	char		choiceStr[80];						/* Concatenated error string */
	char		errorBuffer[LBINMaxErrorPrompt];	/* Buffer for error string */
	LbTwoByte	inputLen;							/* Length of input string */
	LbTwoByte	result = -1;						/* Calculated result */

	#ifdef LB_DEBUG
		/* Verify the gave us a prompt */
		if (prompt == 0)
			ErAbort("You must supply a prompt to LbInAskFiltered");
					
		/* Verify they gave use a destination */
		if (response == 0)
			ErAbort("You must supply space for result in LbInAskFiltered");
	#endif
	
	/* Assume no cancelation if allowed */
	if (canceled != 0)
		*canceled = false;
		
	/* Print the question */
	lbInPrintQuestion(prompt, defaultChoice, choice1, choice2, choice3, choice4);
		
	do	{
		
		/* Get the answer */
		(void) LbInGetInput(lbInInputBuffer, 255);
		
		/* Get length of input */
		inputLen = strlen(lbInInputBuffer);

		/* See what the answered */
		{
			/* See if the took the default */
			if (inputLen == 0) {
				/* If they supplied a default choice then this is it */
				if (defaultChoice != 0) {
					result = defaultChoice;
					goofed = false;
					
					/* Copy into their response buffer if they supplied one */
					if (response != 0) {
						if (defaultChoice == 1)
							strcpy(response, choice1);
						else if (defaultChoice == 2)
							strcpy(response, choice2);
						else if (defaultChoice == 3)
							strcpy(response, choice3);
						else
							strcpy(response, choice4);
					}
					finished = true;
				}
				else {
					strcpy(errorBuffer, "You must answer or cancel");
					goofed = true;
				}
			}
			else if ((canceled != 0) && LbInMatchStr(lbInInputBuffer, lbInCancelStr))
			{
				/* They supplied cancel flag and the user canceled */
				*canceled = true;
				finished = true;
			}
			else if (!LBINFg_IsUnique(flags) && !LBINFg_IsNumeric(flags)) {
			
				if ((choice1 != 0) && LbInMatchStr(lbInInputBuffer, choice1)) {
					/* They took first choice */
					result = 1;
					strcpy(response, choice1);
					finished = true;
				}
				else if ((choice2 != 0) && LbInMatchStr(lbInInputBuffer, choice2)) {
					/* They took second choice */
					result = 2;
					strcpy(response, choice2);
					finished = true;
				}
				else if ((choice3 != 0) && LbInMatchStr(lbInInputBuffer, choice3)) {
					/* They took third choice */
					result = 3;
					strcpy(response, choice3);
					finished = true;
				}
				else if ((choice4 != 0) && LbInMatchStr(lbInInputBuffer, choice4)) {
					/* They took fourth choice */
					result = 4;
					strcpy(response, choice4);
					finished = true;
				}
			}
			else if (!LBINFg_IsNumeric(flags)) {

				/* This clause is the same as the one above, only the matching must
					be exact, not a substring.
				*/
				if ((choice1 != 0) && (strcmp(lbInInputBuffer, choice1) == 0)) {
					/* They took first choice */
					result = 1;
					strcpy(response, choice1);
					finished = true;
				}
				else if ((choice2 != 0) &&  (strcmp(lbInInputBuffer, choice2) == 0)) {
					/* They took second choice */
					result = 2;
					strcpy(response, choice2);
					finished = true;
				}
				else if ((choice3 != 0) &&  (strcmp(lbInInputBuffer, choice3) == 0)) {
					/* They took third choice */
					result = 3;
					strcpy(response, choice3);
					finished = true;
				}
				else if ((choice4 != 0) &&  (strcmp(lbInInputBuffer, choice4) == 0)) {
					/* They took fourth choice */
					result = 4;
					strcpy(response, choice4);
					finished = true;
				}
			}	
			else if (LBINFg_IsNumeric(flags)) {

				/* This clause is the same as the one above, only the matching must
					be exact, not a substring.
				*/
				if ((choice1 != 0) && (strcmp(lbInInputBuffer, choice1) == 0)) {
					/* They took first choice */
					result = 1;
					strcpy(response, choice1);
					finished = true;
				}
				else if ((choice2 != 0) &&  (strcmp(lbInInputBuffer, choice2) == 0)) {
					/* They took second choice */
					result = 2;
					strcpy(response, choice2);
					finished = true;
				}
				else if ((choice3 != 0) &&  (strcmp(lbInInputBuffer, choice3) == 0)) {
					/* They took third choice */
					result = 3;
					strcpy(response, choice3);
					finished = true;
				}
				else if ((choice4 != 0) &&  (strcmp(lbInInputBuffer, choice4) == 0)) {
					/* They took fourth choice */
					result = 4;
					strcpy(response, choice4);
					finished = true;
				}
			}	
			
			/* If were still aren't finished then do filtering here */
			if (finished == false){
				/* Call the filter function */
				filterResult = (*filterFunction)(lbInInputBuffer, &isValid, errorBuffer,
					userData);
					
				/* See if filter failed */
				if (filterResult == false) {
					#ifdef LB_DEBUG
						/* Verify if filter could fail, they gave us storage */
						if (userResult == 0)
							ErAbort("If your filter can fail, you must give me storage");
					#endif
					*userResult = false;
					finished = true;
					goofed = false;
				}
				
				/* If valid answer, we are done */
				if (isValid) {
					strcpy(response, lbInInputBuffer);
					finished = true;
					goofed = false;
					result = LBINUserChoice;
				}
				else
					goofed = true;
			
			}

			/* See if they goofed */
			if ((finished == false) && (goofed == true)) {
				lbInBuildChoiceStr(choiceStr, defaultChoice, 
					choice1, choice2, choice3, choice4);
					
				/* Print answers again */
				LbInPrintf("%s %s:", errorBuffer, choiceStr);
			}
		}

	} while (!finished);
	
	return (result);
}

/*****************************************
*		LbInAskFloat
*
*	Arguments:
*				char		*prompt				- Prompt to display.
*				LbTwoByte	defaultChoice		- Number of choice that is the defaultChoice.
*				Boolean		isSigned			- Is input isSigned value?
*				LbTwoByte	decimalPrec			- How many digits of decimal precision?
*				Boolean		restrictToChoices	- Flag to restrict input to given choices.
*				LbTwoByte	numChoices			- The number of choices being supplied.
*				Boolean		*canceled			- Set true if user cancels input.
*				double		choice1				- First choice given.
*				double		choice2				- Second choice given.
*				double		choice3				- Third choice given.
*				double		choice4				- Fourth choice given.
*
*	Purpose:	Present the user with a prompt and get their response.
*
*	Returns:	The number entered if user did not cancel.
*
******************************************/
double	LbInAskFloat(char *prompt, LbTwoByte defaultChoice, Boolean isSigned,
				LbTwoByte decimalPrec,
				Boolean restrictToChoices,
				LbTwoByte numChoices, Boolean *canceled, 
				double choice1, double choice2, double choice3, double choice4)
{
	LbTwoByte			userResult;			/* The result that comes back from LbInAsk */
	char				formatStr[8];		/* Format string for conversion */
	char				choice1Buff[128];	/* String buffer of choice1 */
	char				choice2Buff[128];	/* String buffer of choice2 */
	char				choice3Buff[128];	/* String buffer of choice3 */
	char				choice4Buff[128];	/* String buffer of choice4 */
	char				*choice1Ptr;		/* String conversion of choice1 */
	char				*choice2Ptr;		/* String conversion of choice2 */
	char				*choice3Ptr;		/* String conversion of choice3 */
	char				*choice4Ptr;		/* String conversion of choice4 */
	char				userStr[128];		/* String result for user choice */
	double 				result;				/* Float conversion of string input */
	lbInFloatFilterTy	filterInfo;			/* Our filter info */
	
	#ifdef LB_DEBUG
		/* Verify they gave choices if they restrict to choices */
		if ((restrictToChoices == true) && (numChoices <= 1)) {
			ErAbort("Why restrict to choices without supplying them to LbInAskFiltered");
		}
	#endif

	/* Setup format string */
	sprintf(formatStr, "%%.%df", decimalPrec);

	/* Convert choices to strings */
	{
		/* Setup choice strings */
		if (numChoices >= 1) {
			choice1Ptr = choice1Buff;
			sprintf(choice1Ptr, formatStr, choice1);
		}
		else
			choice1Ptr = 0;
		
		if (numChoices >= 2) {
			choice2Ptr = choice2Buff;
			sprintf(choice2Ptr, formatStr, choice2);
		}
		else
			choice2Ptr = 0;
			
		if (numChoices >= 3) {
			choice3Ptr = choice3Buff;
			sprintf(choice3Ptr, formatStr, choice3);
		}
		else
			choice3Ptr = 0;
			
		if (numChoices >= 4) {
			choice4Ptr = choice4Buff;
			sprintf(choice4Ptr, formatStr, choice4);
		}
		else
			choice4Ptr = 0;
	}
	
	/* Initialize filter information */
	filterInfo.RestrictToChoices = restrictToChoices;
	filterInfo.IsSigned = isSigned;
	filterInfo.Min = filterInfo.Max = 0;
	
	/* Use generalized filter ask */
	userResult = LbInAskFiltered(LBINFg_Numeric, prompt, defaultChoice, canceled, 
		choice1Ptr, choice2Ptr, choice3Ptr, choice4Ptr,
		userStr, lbInFloatFilter, 0, (LbUserData) &filterInfo);

	/* Check results */
	if ((canceled == 0) || *canceled == false)
	{
		if (userResult == 1)
			result = choice1;
		else if (userResult == 2)
			result = choice2;
		else if (userResult == 3)
			result = choice3;
		else if (userResult == 4)
			result = choice4;
		else if (userResult == LBINUserChoice)
			result = (double) atof(userStr);
		#ifdef LB_DEBUG
			else
				ErAbort("Unexepected result in LbInAskFloat");
		#endif
	}
	
	return (result);
}

/*****************************************
*		LbInAskFourByte
*
*	Arguments:
*				char			*prompt				- Prompt to display.
*				LbTwoByte		defaultChoice		- Number of choice that is the defaultChoice.
*				Boolean			isSigned				- Is input isSigned value?
*				Boolean			restrictToChoices	- Flag to restrict input to given choices.
*				LbTwoByte		numChoices			- The number of choices being supplied.
*				Boolean			*canceled			- Set true if user cancels input.
*				LbFourByte		choice1				- First choice given.
*				LbFourByte		choice2				- Second choice given.
*				LbFourByte		choice3				- Third choice given.
*				LbFourByte		choice4				- Fourth choice given.
*
*	Purpose:	Present the user with a prompt and get their response.
*
*	Returns:	The number entered if user did not cancel.
*
******************************************/
LbFourByte	LbInAskFourByte(char *prompt, LbTwoByte defaultChoice, Boolean isSigned,
				Boolean restrictToChoices, LbTwoByte numChoices, Boolean *canceled, 
				LbFourByte choice1, LbFourByte choice2, LbFourByte choice3, LbFourByte choice4)
{
	char				choice1Buff[32];	/* String buffer of choice1 */
	char				choice2Buff[32];	/* String buffer of choice2 */
	char				choice3Buff[32];	/* String buffer of choice3 */
	char				choice4Buff[32];	/* String buffer of choice4 */
	char				*choice1Ptr;		/* String conversion of choice1 */
	char				*choice2Ptr;		/* String conversion of choice2 */
	char				*choice3Ptr;		/* String conversion of choice3 */
	char				*choice4Ptr;		/* String conversion of choice4 */
	char				userStr[32];		/* String result for user choice */
	lbInFourFilterTy	filterInfo;			/* Our filter info */
	LbFourByte 			result = -1;
	LbTwoByte			userResult;			/* The result that comes back from LbInAsk */
	
	#ifdef LB_DEBUG
		/* Verify they gave choices if they restrict to choices */
		if ((restrictToChoices == true) && (numChoices <= 1)) {
			ErAbort("Why restrict to choices without supplying them to LbInAskFiltered");
		}
	#endif

	/* Convert choices to strings */
	{
		/* Setup choice strings */
		if (numChoices >= 1) {
			choice1Ptr = choice1Buff;
			sprintf(choice1Ptr, "%ld", (long)choice1);
		}
		else
			choice1Ptr = 0;
		
		if (numChoices >= 2) {
			choice2Ptr = choice2Buff;
			sprintf(choice2Ptr, "%ld", (long)choice2);
		}
		else
			choice2Ptr = 0;
			
		if (numChoices >= 3) {
			choice3Ptr = choice3Buff;
			sprintf(choice3Ptr, "%ld", (long)choice3);
		}
		else
			choice3Ptr = 0;
			
		if (numChoices >= 4) {
			choice4Ptr = choice4Buff;
			sprintf(choice4Ptr, "%ld", (long)choice4);
		}
		else
			choice4Ptr = 0;
	}
	
	/* Initialize filter information */
	filterInfo.RestrictToChoices = restrictToChoices;
	filterInfo.IsSigned = isSigned;
	filterInfo.Min = filterInfo.Max = 0;
	
	/* Use generalized filter ask */
	userResult = LbInAskFiltered(LBINFg_Numeric, prompt, defaultChoice, canceled, 
		choice1Ptr, choice2Ptr, choice3Ptr, choice4Ptr,
		userStr, lbInFourFilter, 0, (LbUserData) &filterInfo);

	/* Check results */
	if ((canceled == 0) || *canceled == false)
	{
		if (userResult == 1)
			result = choice1;
		else if (userResult == 2)
			result = choice2;
		else if (userResult == 3)
			result = choice3;
		else if (userResult == 4)
			result = choice4;
		else if (userResult == LBINUserChoice)
			result = (long) atol(userStr);
		#ifdef LB_DEBUG
			else
				ErAbort("Unexepected result in LbInAskFourByte");
		#endif
	}
	
	return (result);
}

/*****************************************
*		LbInAskFourByteInRange
*
*	Arguments:
*				char		*prompt				- Prompt to display.
*				Boolean		isSigned			- Is value signed.
*				Boolean		*canceled			- Set true if user cancels input.
*				LbFourByte	min					- Minimum range.
*				LbFourByte	max					- Maximum range.
*
*	Purpose:	Present the user with a prompt and get their response.
*
*	Returns:	The value entered.
*
******************************************/
LbFourByte	LbInAskFourByteInRange(char *prompt,
				Boolean isSigned, Boolean *canceled, LbFourByte min, LbFourByte max)
{
	char				minStrBuff[32];		/* String buffer of min choice */
	char				maxStrBuff[32];		/* String buffer of max choice */
	char				userStr[32];		/* String result for user choice */
	LbFourByte			result;				/* Calculated result */
	lbInFourFilterTy	filterInfo;			/* Our filter's info */
	LbTwoByte			userResult;			/* The result that comes back from LbInAskFiltered */
	
	
	/* Convert choices to strings */
	{
		/* Setup choice strings */
		sprintf(minStrBuff, "%ld", (long)min);
		sprintf(maxStrBuff, "%ld", (long)max);
	}
	
	/* Initialize filter info */
	filterInfo.Min = min;
	filterInfo.Max = max;
	filterInfo.IsSigned = isSigned;
	filterInfo.RestrictToChoices = false;
	
	/* Use generalized filter ask */
	userResult = LbInAskFiltered(LBINFg_Numeric, prompt, 0, canceled, minStrBuff, maxStrBuff,
		0, 0, userStr, lbInFourRngFilter, 0, (LbUserData) &filterInfo);

	/* Check results */
	if ((canceled == 0) || *canceled == false)
	{
		if (userResult == 1)
			result = min;
		else if (userResult == 2)
			result = max;
		else if (userResult == LBINUserChoice)
			result = (double) atol(userStr);
		#ifdef LB_DEBUG
			else
				ErAbort("Unexepected result in LbInAskFourByteInRange");
		#endif
	}
	
	return (result);
}

/*****************************************
*		LbInAskTwoByte
*
*	Arguments:
*				char			*prompt				- Prompt to display.
*				LbTwoByte		defaultChoice		- Number of choice that is the defaultChoice.
*				Boolean			isSigned				- Is input isSigned value?
*				Boolean			restrictToChoices	- Flag to restrict input to given choices.
*				LbTwoByte		numChoices			- The number of choices being supplied.
*				Boolean			*canceled			- Set true if user cancels input.
*				LbTwoByte		choice1				- First choice given.
*				LbTwoByte		choice2				- Second choice given.
*				LbTwoByte		choice3				- Third choice given.
*				LbTwoByte		choice4				- Fourth choice given.
*
*	Purpose:	Present the user with a prompt and get their response.
*
*	Returns:	The number entered if user did not cancel.
*
******************************************/
LbTwoByte	LbInAskTwoByte(char *prompt, LbTwoByte defaultChoice, Boolean isSigned,
				Boolean restrictToChoices, LbTwoByte numChoices, Boolean *canceled, 
				LbTwoByte choice1, LbTwoByte choice2, LbTwoByte choice3, LbTwoByte choice4)
{
	char				choice1Buff[32];	/* String buffer of choice1 */
	char				choice2Buff[32];	/* String buffer of choice2 */
	char				choice3Buff[32];	/* String buffer of choice3 */
	char				choice4Buff[32];	/* String buffer of choice4 */
	char				*choice1Ptr;		/* String conversion of choice1 */
	char				*choice2Ptr;		/* String conversion of choice2 */
	char				*choice3Ptr;		/* String conversion of choice3 */
	char				*choice4Ptr;		/* String conversion of choice4 */
	char				userStr[32];		/* String result for user choice */
	lbInTwoFilterTy		filterInfo;			/* Our filter info */
	LbTwoByte 			result = -1;
	LbTwoByte			userResult;			/* The result that comes back from LbInAsk */
	
	#ifdef LB_DEBUG
		/* Verify they gave choices if they restrict to choices */
		if ((restrictToChoices == true) && (numChoices <= 1)) {
			ErAbort("Why restrict to choices without supplying them to LbInAskFiltered");
		}
	#endif

	/* Convert choices to strings */
	{
		/* Setup choice strings */
		if (numChoices >= 1) {
			choice1Ptr = choice1Buff;
			sprintf(choice1Ptr, "%d", choice1);
		}
		else
			choice1Ptr = 0;
		
		if (numChoices >= 2) {
			choice2Ptr = choice2Buff;
			sprintf(choice2Ptr, "%d", choice2);
		}
		else
			choice2Ptr = 0;
			
		if (numChoices >= 3) {
			choice3Ptr = choice3Buff;
			sprintf(choice3Ptr, "%d", choice3);
		}
		else
			choice3Ptr = 0;
			
		if (numChoices >= 4) {
			choice4Ptr = choice4Buff;
			sprintf(choice4Ptr, "%d", choice4);
		}
		else
			choice4Ptr = 0;
	}
	
	/* Initialize filter information */
	filterInfo.RestrictToChoices = restrictToChoices;
	filterInfo.IsSigned = isSigned;
	
	/* Use generalized filter ask */
	userResult = LbInAskFiltered(LBINFg_Numeric, prompt, defaultChoice, canceled, 
		choice1Ptr, choice2Ptr, choice3Ptr, choice4Ptr,
		userStr, lbInTwoFilter, 0, (LbUserData) &filterInfo);

	/* Check results */
	if ((canceled == 0) || *canceled == false)
	{
		if (userResult == 1)
			result = choice1;
		else if (userResult == 2)
			result = choice2;
		else if (userResult == 3)
			result = choice3;
		else if (userResult == 4)
			result = choice4;
		else if (userResult == LBINUserChoice)
			result = (short) atoi(userStr);
		#ifdef LB_DEBUG
			else
				ErAbort("Unexepected result in LbInAskTwoByte");
		#endif
	}
	
	return (result);
}

/*****************************************
*		LbInAskYesNo
*
*	Arguments:
*				char			*prompt				- Prompt to display.
*				LbTwoByte		defaultChoice		- Number of choice that is the defaultChoice.
*
*	Purpose:	Present the user with a prompt and get their response.
*				Restrict answer to Yes or No.
*
*	Returns:	LBINYes, LBINNo
*
******************************************/
LbTwoByte	LbInAskYesNo(char *prompt, LbTwoByte defaultChoice)
{
	LbTwoByte	result = -1;		/* Calculated result */

	/* Just use general ask routine */
	result = LbInAsk(prompt, defaultChoice, true, 0, lbInYesStr, lbInNoStr, 0, 0, 0);
	
	/* Check results */
	if (result == 1)
		result = LBINYes;
	else if (result == 2)
		result = LBINNo;
	#ifdef LB_DEBUG
		else
			ErAbort("Invalid result from LbInAsk in LbInAskYesNo");
	#endif
		
	return (result);
}

/*****************************************
*		LbInAskYesNoCancel
*
*	Arguments:
*				char			*prompt				- Prompt to display.
*				LbTwoByte		defaultChoice		- Number of choice that is the defaultChoice.
*
*	Purpose:	Present the user with a prompt and get their response.
*				Restrict answer to Yes, No, or Cancel.
*
*	Returns:	LBINYes, LBINNo, LBINCancel
*
******************************************/
LbTwoByte	LbInAskYesNoCancel(char *prompt, LbTwoByte defaultChoice)
{
	LbTwoByte result = -1;

	/* Just use general ask routine */
	result = LbInAsk(prompt, defaultChoice, true, 0, lbInYesStr, lbInNoStr, lbInCancelStr, 0, 0);
	
	/* Check results */
	if (result == 1)
		result = LBINYes;
	else if (result == 2)
		result = LBINNo;
	else if (result == 3)
		result = LBINCancel;
	#ifdef LB_DEBUG
		else
			ErAbort("Invalid result from LbInAsk in LbInAskYesNoCancel");
	#endif
	
	return (result);
}

/*****************************************
*		LbInAskYesNoQuit
*
*	Arguments:
*				char			*prompt				- Prompt to display.
*				LbTwoByte		defaultChoice		- Number of choice that is the defaultChoice.
*
*	Purpose:	Present the user with a prompt and get their response.
*				Restrict answer to Yes, No, or Quit.
*
*	Returns:	LBINYes, LBINNo, LBINQuit
*
******************************************/
LbTwoByte	LbInAskYesNoQuit(char *prompt, LbTwoByte defaultChoice)
{
	LbTwoByte result = -1;

	/* Just use general ask routine */
	result = LbInAsk(prompt, defaultChoice, true, 0, lbInYesStr, lbInNoStr, lbInQuitStr, 0, 0);
	
	/* Check results */
	if (result == 1)
		result = LBINYes;
	else if (result == 2)
		result = LBINNo;
	else if (result == 3)
		result = LBINQuit;
	#ifdef LB_DEBUG
		else
			ErAbort("Invalid result from LbInAsk in LbInAskYesNoQuit");
	#endif
	
	return (result);
}

/*****************************************
*		lbInFloatFilter
*
*	Arguments:
*				char			*inputString		- The input string.
*				Boolean			*isValid			- Result of filter.
*				char			*errorPrompt		- String to chastise user.
*				LbUserData		userData	- Our user data.
*
*	Purpose:	Filter input for doubles.
*
*	Returns:	True unless an error occurs.
*
******************************************/
Boolean	lbInFloatFilter(char *inputString, Boolean *isValid, char *errorPrompt,
			LbUserData userData)
{
	Boolean				okay = false;
	double				testValue;
	lbInFloatFilterTy	*filterInfo = (lbInFloatFilterTy *) userData;
	LbTwoByte			parseIndex;
	
	do	{ /* Process Loop */
	
		/* If not restricting to choices verify number is valid */
		if (filterInfo->RestrictToChoices == false) {
			
			/* Assume value is good */
			*isValid = true;
			
			/* If not signed, and number contains sing value, its bogus */
			if ((filterInfo->IsSigned == false) && (inputString[0] == '-')){
				strcpy(errorPrompt, "Negative values not allowed");
				*isValid = false;
			}
			else {
				/* Make sure it is not alphabetic */
				parseIndex = 0;
				while (inputString[parseIndex] != 0) {
					if ((isdigit(inputString[parseIndex]) == false)
							&& (inputString[parseIndex] != '.')
							&& (inputString[parseIndex] != '-')) {
							
						/* Value is invalid */
						strcpy(errorPrompt, "Input must be numeric");
						*isValid = false;
						break;
					}
					
					parseIndex++;
				}
			}
		}
		else {
			*isValid = false;
			strcpy(errorPrompt, "You must answer");
		}
		
		/* If still valid, make sure its not going out of range, or loosing precision */
		if (*isValid) {
			testValue = (double) fabs(atof(inputString));
			if (testValue > LBFLOAT_MAX) {
				strcpy(errorPrompt, "Value out of range.");
				*isValid = false;
			}
		}
		
		okay = true;
	} while (false);

	return (okay);
}

/*****************************************
*		lbInFourFilter
*
*	Arguments:
*				char			*inputString		- The input string.
*				Boolean			*isValid			- Result of filter.
*				char			*errorPrompt		- String to chastise user.
*				LbUserData		userData	- Our user data.
*
*	Purpose:	Filter input for four bytes.
*
*	Returns:	True unless an error occurs.
*
******************************************/
Boolean	lbInFourFilter(char *inputString, Boolean *isValid, char *errorPrompt,
			LbUserData userData)
{
	Boolean				okay = false;
	lbInFourFilterTy	*filterInfo = (lbInFourFilterTy *) userData;
	LbTwoByte			parseIndex;
	
	do	{ /* Process Loop */
	
		/* If not restricting to choices verify number is valid */
		if (filterInfo->RestrictToChoices == false) {
			
			/* Assume value is good */
			*isValid = true;
			
			/* If not signed, and number contains sign value, its bogus */
			if ((filterInfo->IsSigned == false) && (inputString[0] == '-')){
				strcpy(errorPrompt, "Negative values not allowed");
				*isValid = false;
			}
			else {
				/* Make sure it is not alphabetic */
				parseIndex = 0;
				while (inputString[parseIndex] != 0) {
					if ((isdigit(inputString[parseIndex]) == false)
							&& (inputString[parseIndex] != '-')) {
							
						/* Value is invalid */
						strcpy(errorPrompt, "Input must be numeric");
						*isValid = false;
						break;
					}
					
					parseIndex++;
				}
			}
		}
		else {
			*isValid = false;
			strcpy(errorPrompt, "You must answer");
		}
		
		okay = true;
	} while (false);

	return (okay);
}

/*****************************************
*		lbInFourRngFilter
*
*	Arguments:
*				char			*inputString		- The input string.
*				Boolean			*isValid			- Result of filter.
*				char			*errorPrompt		- String to chastise user.
*				LbUserData		userData			- Our user data.
*
*	Purpose:	Filter input for four byte, in given range.
*
*	Returns:	True unless an error occurs.
*
******************************************/
Boolean	lbInFourRngFilter(char *inputString, Boolean *isValid,
			char *errorPrompt, LbUserData userData)
{
	Boolean				okay = false;
	lbInFourFilterTy	*filterInfo = (lbInFourFilterTy *) userData;
	
	do	{ /* Process Loop */
	
		/* First make sure it is a valid foat */
		if (!lbInFourFilter(inputString, isValid, errorPrompt, userData))
			break;
			
		/* If still valid here, check the range */
		if (*isValid == true) {
			if ((atol(inputString) < filterInfo->Min) || (atol(inputString) > filterInfo->Max)) {
			
				/* Value is invalid */
				strcpy(errorPrompt, "Input is out of range");
				*isValid = false;
			}
		}
		
		okay = true;
	} while (false);

	return (okay);
}

/*****************************************
*		lbInBuildChoiceStr
*
*	Arguments:
*				char		*choiceStr			- Concatenated choice string.
*				LbTwoByte	defaultChoice		- Number of choice that is the defaultChoice.
*				char		*choice1			- First choice given.
*				char		*choice2			- Second choice given.
*				char		*choice3			- Third choice given.
*				char		*choice4			- Fourth choice given.
*
*	Purpose:	Build the choice string.
*
*	Returns:	Nothing.
*
******************************************/
void	lbInBuildChoiceStr(char *choiceStr, LbTwoByte defaultChoice,
				char *choice1, char *choice2, char *choice3, char *choice4)
{
	char		*defaultString;				/* Holder for default string */
	char		*emptyString = "";				/* Empty string */
	
	
	/* Set default string */
	switch (defaultChoice) {
		case 1:
			defaultString = choice1;
			break;
		
		case 2:
			defaultString = choice2;
			break;
			
		case 3:
			defaultString = choice3;
			break;
			
		case 4:
			defaultString = choice4;
			break;
			
		default:
			defaultString = 0;
			break;
	}
	
	/* Set choices to empty strings if not supplied */
	{
		if (!choice1)
			choice1 = emptyString;
			
		if (!choice2)
			choice2 = emptyString;
			
		if (!choice3)
			choice3 = emptyString;
			
		if (!choice4)
			choice4 = emptyString;
	}
		
	/* See if no choices */
	if (strlen(choice1) == 0) {
		choiceStr[0] = 0;
	}
	else if (strlen(choice2) == 0) {
		sprintf(choiceStr, " [%s]", choice1);
	}
	else if (strlen(choice3) == 0) {
		if (defaultString != 0)
			sprintf(choiceStr, " (%s..%s) [%s]", choice1, choice2, defaultString);
		else
			sprintf(choiceStr, " (%s..%s)", choice1, choice2);		
	}
	else if (strlen(choice4) == 0) {
		if (defaultString != 0)
			sprintf(choiceStr, " (%s,%s,%s) [%s]", choice1, choice2, choice3, defaultString);
		else
			sprintf(choiceStr, " (%s,%s,%s)", choice1, choice2, choice3);		
	}
	else {
		if (defaultString != 0)
			sprintf(choiceStr, " (%s,%s,%s,%s) [%s]", choice1, choice2, choice3, choice4, defaultString);
		else
			sprintf(choiceStr, " (%s,%s,%s,%s)", choice1, choice2, choice3, choice4);
	}
}

/*****************************************
*		LbInGetInput
*
*	Arguments:
*				char	*inputBuffer		- Storage for user input.
*				int		bufferLength		- number of bytes in inputBuffer
*
*	Purpose:	Read user input from standard in (the keyboard).
*
*	Returns:	lbInInputBuffer
*
******************************************/
char	*LbInGetInput(char *inputBuffer, int bufferLength)
{
	lbInGetS( inputBuffer, bufferLength );
	
	/* If we are running in MPW, we have to strip off the line */
	#ifdef MPW
	{
		LbTwoByte	replaceIndex = 0;		/* Index for replacement chars */
		LbTwoByte	length;					/* Length of input */
		
		/* Find the : */
		length = strlen(inputBuffer);
		length--;
		
		while (inputBuffer[length] != ':')
		{
			/* Make sure we didn't get weird input */
			if (length == 0)
				ErAbort("Weird input from LbInGetInput");
			
			/* Go to next character */
			length--;
		}
		
		/* Increment past colon */
		length++;
		
		/* Copy remainder to front of string */
		do
		{
			/* Copy the current character to the beginning of the buffer */
			inputBuffer[replaceIndex++] = inputBuffer[length];
			
		} while (inputBuffer[length++] != 0);
	}
	#endif
	
	return (inputBuffer);
}

/*****************************************
*		LbInMatchStr
*
*	Arguments:
*				char	*strToCheck		- String that will be matched.
*				char	*strToMatch		- String that we are checking against
*
*	Purpose:	Check to strings to see strToCheck is a valid subset of
*			strToMatch.
*
*	Returns:	TRUE if they match.
*
******************************************/
Boolean	LbInMatchStr(char *strToCheck, char *strToMatch)
{
	Boolean		stringsMatch = false;	/* Flag for whether they match or not */
	LbTwoByte	checkLen;				/* Length of string we are checking */
	LbTwoByte	matchLen;				/* Length of match string */
	LbTwoByte	index;					/* Index for checking */

	/* Get string lengths */
	checkLen = strlen(strToCheck);
	matchLen = strlen(strToMatch);
	
	do	/* Process Loop */
	{
		/* If check string is longer than match, they don't */
		if (checkLen > matchLen)
			break;
		
		/* Loop through comparing characters */
		index = 0;
		do
		{
			/* See if current character matches */
			stringsMatch = (toupper(strToCheck[index]) == toupper(strToMatch[index]));
			
			/* Go to next char */
			index++;
			
		} while (stringsMatch && (index < checkLen));
		
	} while (false);
	return (stringsMatch);
}

/*****************************************
*		LbInPrintf
*
*	Arguments:
*				char *fmt		- Format for string conversion.
*				Variable argument list
*
*	Purpose:	Emulate stdio "LbInPrintf" function.
*
*	Returns:	Number of characters printed, or -1 on error.
*
******************************************/
LbFourByte	LbInPrintf(char *fmt, ...)
{
	LbFourByte	numChars = -1;	/* Number of characters printed */
	va_list		args;
	
	/* Convert variable length argument list and process
		using vsprintf
	*/
	va_start(args, fmt);
	vsprintf(lbInPrintfString, fmt, args);
	va_end(args);
	
	/* Compute string length */
	numChars = strlen(lbInPrintfString);
	
	/* Print the string */
	#ifdef COCOA
		/* Cocoa program (GUI) */
		WriteStrCocoa( lbInPrintfString );
	#else
		/* Command-line OS */
		printf("%s", lbInPrintfString);
	#endif
	
	return (numChars);
}

/*****************************************
*		lbInPrintQuestion
*
*	Arguments:
*				char		*question			- Prompt to display.
*				LbTwoByte	defaultChoice		- Number of choice that is the defaultChoice.
*				char		*choice1			- First choice given.
*				char		*choice2			- Second choice given.
*				char		*choice3			- Third choice given.
*				char		*choice4			- Fourth choice given.
*
*	Purpose:	Present the user with a prompt.
*
*	Returns:	Nothing.
*
******************************************/
void	lbInPrintQuestion(char *prompt, LbTwoByte defaultChoice,
				char *choice1, char *choice2, char *choice3, char *choice4)
{
	char		choiceStr[LBIN_MAXLINE];	/* Line built up for choices */	
	
	/* Build up the choice string */
	lbInBuildChoiceStr(choiceStr, defaultChoice, choice1, choice2, choice3, choice4);
	
	/* See if we must go to two lines */
	if ((strlen(prompt) + strlen(choiceStr) + 1) > LBIN_MAXLINE)
		LbInPrintf("%s\n%s: ", prompt, choiceStr);
	else
		LbInPrintf("%s%s: ", prompt, choiceStr);
		
}

/*****************************************
*		lbInGetS
*
*	Arguments:
*				char		*str	- Returned input string.
*				int			size	- Size of string memory.
*
*	Purpose:	Read str from stdin, up to size (equiv to gets).
*
*	Returns:	char* to str
*
******************************************/
char	*lbInGetS( char *str, int size )
{
	char		*result;		/* Function result */
	char		*p;				/* Pointer into str */
	int			i;				/* Index through str */
	char		c;				/* Read-in character */
	
	
	if ( size <= 0 ) {
		/* Can't accept 0 size */
		result = NULL;
		ErAbort("Attempt to get zero length string in lbInGetS.");
	}
	else if ( str == NULL ) {
		/* str must be initialized */
		result = NULL;
		ErAbort("Attempt to fill uninitialized string in lbInGetS.");
	}
	else {
		result = str;	/* Non-null initialization */
		*str = '\0';	/* String initialization */
		p = str;
		for ( i=0; i<size; i++ ) {
			/* Read a character from stdin */
			c = getc( stdin );
			
			if ( (int) c == 10 ) {	/* LF */
				/* End of line */
				*p = '\0';
				
				break;
			}
			else if ( (int) c == 13 ) {	/* CR */
				/* End of line */
				*p = '\0';
				
				break;
			}
			else {
				/* Normal character */
				*p++ = c;
			}
		}
		
	}
	
	return ( result );
}
