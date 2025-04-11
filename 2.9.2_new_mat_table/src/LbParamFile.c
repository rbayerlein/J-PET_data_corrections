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
*     Module Name:        LbParamFile.c
*     Revision Number:    1.4
*     Date last revised:  23 July 2013
*     Programmer:         Steven Vannoy
*     Date Originated:    Monday, March 1, 1993
*
*     Module Overview:	This module provides routines for parsing a parameter
*						file.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*     	LbPfOpen
*		LbPfGetParam
*		LpPfClose
*		LbPfGetNextToken
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
*			Revision date:		19 June 2012
*
*			Revision description:	Added LbPfGetNextToken
*									Modified lbPfGetToken to use LbPfGetNextToken
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*
*			Revision description:
*						- support for eight-byte integers
*
*********************************************************************************/

#include	<stdio.h>
#include	<string.h>
#include	<ctype.h>

#include	"SystemDependent.h"

#include	"LbTypes.h"
#include	"LbError.h"
#include	"LbFile.h"
#include	"LbParamFile.h"

/*	LOCAL CONSTANTS */
#define LBPF_PARAM_LEN_MAX		256
#define LBPF_COMMENT_INDEX		1

/*  LOCAL GLOBALS */
char			*lbPfTokenStrings[] = {"BOOL", "#", "CHAR", "STR", "INT", "REAL", "NUM_ELEMENTS_IN_LIST", "ENUM", "LONGLONG"};
char			*lbPfTokenStringAltList[] = {"BOOL", "#", "CHAR", "STR", "INT", "REAL", "LIST"};
char			lbPfErrStr[1024];
LbUsTwoByte		lbPfNumTokenTypes = 9;
LbUsTwoByte		lbPfNumAltTokenTypes = 7;
LbUsFourByte	lbPfParamFileCurLine;

/*	LOCAL MACROS */
/*	LOCAL FUNCTIONS */
Boolean	lbPfParseParamString(char *paramBufferPtr,
			LbPfEnPfTy *paramTypePtr, LbUsFourByte *paramSizePtr,
			char *labelPtr);

Boolean	lbPfGetToken(char *paramBufferPtr,
			char *tokenPtr, LbUsFourByte *currentCharPtr, 
			LbUsFourByte *tokenLenPtr,
			Boolean *isLastTokenPtr);

Boolean	lbPfLineIsBlank(char *paramBufferPtr);

/*	LOCAL MACROS */

/*	FUNCTIONS */
/*****************************************
*		LbPfClose
*
*	Purpose:	Close the parameter file.
*	Arguments:
*		LbPfHkTy		*pfFileHkPtr 		- The parameter file hook.
*
*	Returns:	None.
*
******************************************/
void	LbPfClose(LbPfHkTy *pfFileHk)
{
	/* Close the file */
	fclose(pfFileHk->fileRef);

	/* Clear their hook */
	pfFileHk->fileRef = 0;	
	pfFileHk->curLine = 0;	
}

/*****************************************
*		LbPfGetParam
*
*	Purpose:	Get the next parameter.
*	Arguments:
*		LbPfHkTy		*pfFileHk 		- The parameter file hook.
*		void			*paramBufferPtr	- The parameter value buffer.
*		LbPfEnPfTy		*paramTypePtr	- The parameter type.
*		LbUsFourByte	*paramSizePtr	- The size of the paramater value data.
*		char			*labelPtr		- The parameters label.
*		Boolean			*isEOF			- End of file flag.
*
*	Returns:	True unless an error or end of file occurs.
*
******************************************/
Boolean	LbPfGetParam(LbPfHkTy *pfFileHk, void *paramBufferPtr,
			LbPfEnPfTy *paramTypePtr, LbUsFourByte *paramSizePtr,
			char *labelPtr, Boolean *isEOF)
{
	Boolean		okay = false;		/* Process flag */
	Boolean		gotString = false;	/* Sub-process flag */
	char		errString[1024];	/* Buffer for error strings */
	
	do { /* While parameter is a comment */

		okay = false;
		
		do { /* Process Loop */
		
			/* Always assume we are not at the end of file */
			*isEOF = false;
			
			/* Get the next non-empty line from the file */
			do {
				
				/* Set local flag */
				gotString = false;
				
				/* Attempt to read a line from the file */
				if ((LbFlFGetS((char *)paramBufferPtr,
					LBPF_PARAM_LEN_MAX, pfFileHk->fileRef)) == NULL) {
				
					/* See if we failed due to end of file */
					if (feof(pfFileHk->fileRef) != 0) {
						*isEOF = true;
					}
					else {
						/* Register an error */
						sprintf(errString, "Error reading parameter from file (line %ld).",
							(unsigned long)pfFileHk->curLine);
						ErStGeneric(errString);
					}
					
					break;
				}
				
				/* Set local flag true */
				gotString = true;
				
				/* Update hook's current line */
				pfFileHk->curLine++;
				
				/* Set local global for current line */
				lbPfParamFileCurLine = pfFileHk->curLine;
				
			} while (lbPfLineIsBlank((char *)paramBufferPtr));
	
			/* If we didn't get a string, bolt */
			if (gotString == false)
				break;
	
			/* Parse the parameter string, note possibility of failure due to invalid parameter */
			if (!lbPfParseParamString((char *) paramBufferPtr, paramTypePtr,
					paramSizePtr, labelPtr)) {
				
				sprintf(errString, "Error processing parameter file '%s'\n", pfFileHk->name);
				ErAlert(errString, false);	
				break;
			}
	  
			okay = true;
		} while (false);

		/* See if we failed */
		if (!okay)
			break;
    } while (*paramTypePtr == LbPfEnComment);

	return (okay);
}

/*****************************************
*		LbPfOpen
*
*	Purpose:	Open the parameter file.
*	Arguments:
*		char			*filePathPtr		- Path to the parameter file.
*		LbUsFourByte	flags				- Modifier flags.
*		LbPfHkTy		*pfFileHkPtr 		- The parameter file hook.
*
*	Returns:	True unless an error or end of file occurs.
*
******************************************/
Boolean	LbPfOpen(char *filePathPtr, LbUsFourByte flags, LbPfHkTy *pfFileHk)
{
	Boolean		okay = false;	/* Process flag */
	
	
	/* Eliminate compiler warning about unused parameter */
	if (flags) {};
	
	do { /* Process Loop */
		
		/* Open their file */
		if ((pfFileHk->fileRef = LbFlFileOpen(filePathPtr, "r")) == NULL) {
			sprintf(lbPfErrStr,"Unable to open parameter file named '%s'.",
				filePathPtr);
			ErStGeneric(lbPfErrStr);
			break;
		}
	
		/* Reset current line */
		pfFileHk->curLine = 0;	

		/* Save the file name for error reporting */
		strncpy(pfFileHk->name, filePathPtr, 255);
		
		okay = true;
		
	} while (false);
	
	return (okay);
}

/*****************************************
*		LbPfGetNextToken
*
*	Purpose:	Get the next token from the supplied parameter buffer.
*
*	Arguments:
*		char			*paramBufferPtr			- The parameter value buffer.
*		char			*tokenPtr				- The token buffer.
*		LbUsFourByte	*currentCharIndexPtr	- The current char in the parameter buffer.
*		LbUsFourByte	*tokenLenPtr			- The length of the token.
*		Boolean			*isLastTokenPtr			- Is this the last token on the line.
*
*	Returns:	True unless an error or end of file occurs.
*
******************************************/
Boolean	LbPfGetNextToken( 	char 			*paramBufferPtr,
							char 			*tokenPtr, 
							LbUsFourByte 	*currentCharIndexPtr,
							LbUsFourByte 	*tokenLenPtr,
							Boolean 		*isLastTokenPtr)

{
	Boolean		okay = false;			/* Process flag */
	char		ch;						/* Character currently being considered */
	char		quote = '\0';			/* Current quote character */
	
	do { /* Process Loop */
		
		/* Clear variable flags */
		*isLastTokenPtr = false;
		*tokenLenPtr = 0;
		tokenPtr[0] = '\0';
		
		/* Eat the white space up to the first token */
		while (((paramBufferPtr[*currentCharIndexPtr] == ' ') ||
				(paramBufferPtr[*currentCharIndexPtr] == '\t') ||
				(paramBufferPtr[*currentCharIndexPtr] == '\n')) &&
				(*currentCharIndexPtr < LBPF_PARAM_LEN_MAX)) {
				
			(*currentCharIndexPtr)++;
		}
		
		/* If we reached the end of the buffer without a token it is an error */
		if (*currentCharIndexPtr == LBPF_PARAM_LEN_MAX) {
			break;
		}
		
		/* Copy characters into the token buffer until the next white space is reached */
		*tokenLenPtr = 0;
		do {
			ch = paramBufferPtr[*currentCharIndexPtr];
			
			/* Handle special situations */
			if (ch == '\\' && paramBufferPtr[(*currentCharIndexPtr)+1]) {
				/* Take the quoted character:  \ch  */
				ch = paramBufferPtr[*(++currentCharIndexPtr)];
			}
			else {
				/* Quotations */
				if ((ch == '"') || (ch == '\'')) {
					if (!quote) {
						/* Begin quote */
						quote = ch;
					}
					else if (ch == quote) {
						/* End quote */
						quote = '\0';
					}
				}
			}
			
			/* Copy the current character */
			tokenPtr[*tokenLenPtr] = ch;
			
			(*tokenLenPtr)++;
			(*currentCharIndexPtr)++;
			
			/* See if we reached the end of the line */
			if (paramBufferPtr[*currentCharIndexPtr] == '\n' ||
				    (paramBufferPtr[*currentCharIndexPtr] == '\0')) {
				*isLastTokenPtr = true;
			}
			
		} while (( quote || ((paramBufferPtr[*currentCharIndexPtr] != ' ') &&
								(paramBufferPtr[*currentCharIndexPtr] != '\t')))  &&
				(*isLastTokenPtr == false) &&
				(*currentCharIndexPtr <= LBPF_PARAM_LEN_MAX));
		
		/* If we went past our line length, we are in an error situation */
		if (*currentCharIndexPtr == LBPF_PARAM_LEN_MAX) {
			break;
		}
		
		/* Terminate our token string */
		tokenPtr[*tokenLenPtr] = '\0';
		
		okay = true;
		
	} while (false);
	
	return (okay);
}

/*****************************************
*		lbPfParseParamString
*
*	Purpose:	Parse the parameter string.
*	Arguments:
*		char			*paramBufferPtr	- The parameter value buffer.
*		LbPfEnPfTy		*paramTypePtr	- The parameter type.
*		LbUsFourByte	*paramSizePtr	- The size of the paramater value data.
*		char			*labelPtr		- The parameters label.
*
*	Returns:	True unless an error or end of file occurs.
*
******************************************/
Boolean	lbPfParseParamString(char *paramBufferPtr,
			LbPfEnPfTy *paramTypePtr, LbUsFourByte *paramSizePtr,
			char *labelPtr)
{
	Boolean			okay = false;					/* Process flag */
	Boolean			isLastToken;					/* Last token flag */
	Boolean			switchOkay;						/* Test for switch okay */
	char			tokenBuff[LBPF_PARAM_LEN];		/* Token buffer */
	char			errString[1024];				/* String for creating error strings */
	double			dblValue;						/* Temp double value */
	LbUsFourByte	currentChar;					/* Character counter */
	LbUsFourByte	tokenLen;						/* Current token counter */
	LbUsFourByte	index;							/* Loop control variable */
	LbUsFourByte	stringIndex;					/* LCV for parsing strings */
	
	do { /* Process Loop */
	
		/* Clear our local vars */
		currentChar = 0;
		isLastToken = false;
		
		/* Get the first token */
		if (!lbPfGetToken(paramBufferPtr, tokenBuff, &currentChar,
				&tokenLen, &isLastToken)) {
				
			break;
		}
		
		/* Determine the token type */
		for (index = 0; index < lbPfNumTokenTypes; index++) {
			
			/* Compare against current string */
			if (tokenBuff[0] == *lbPfTokenStrings[LBPF_COMMENT_INDEX]){
				*paramTypePtr = LbPfEnComment;
				break;
			}
			else if (strcmp(tokenBuff, lbPfTokenStrings[index]) == 0) {
				*paramTypePtr = (LbPfEnPfTy) index;
				break;
			}
		}
		
		/* if we didn't find a match, look through alternate list */
		if (index == lbPfNumTokenTypes) {
			/* Determine the token type */
			for (index = 0; index < lbPfNumTokenTypes; index++) {
				
				/* Compare against current string */
				if (tokenBuff[0] == *lbPfTokenStringAltList[LBPF_COMMENT_INDEX]){
					*paramTypePtr = LbPfEnComment;
					break;
				}
				else if (index < lbPfNumAltTokenTypes) {
						if (strcmp(tokenBuff, lbPfTokenStringAltList[index]) == 0) {
						*paramTypePtr = (LbPfEnPfTy) index;
						break;
					}
				}
			}
		}
		
		/* If we get through this without a match, its an error */
		if (index >= lbPfNumTokenTypes) {
			sprintf(errString, "Invalid token type found at beginning of line number %ld:\n token = %s\n from line = %s.",
				(unsigned long)lbPfParamFileCurLine, tokenBuff, paramBufferPtr);
			ErStGeneric(errString);
			break;
		}
		
		/* If its the last token, its only valid if its a comment */
		if ((isLastToken == true) && (*paramTypePtr != LbPfEnComment)) {
			sprintf(errString, "Incomplete parameter specification on line %ld.",
				(unsigned long)lbPfParamFileCurLine);
			ErStGeneric(errString);
			break;
		}
		
		/* Get the label or the comment */
		if (!lbPfGetToken(paramBufferPtr, tokenBuff, &currentChar,
				&tokenLen, &isLastToken)) {
				
			break;
		}
		
		/* If its the last token and its not a comment, 
		   or the token is empty, its an error */
		if ((tokenLen == 0) ||
				((isLastToken == true) && (*paramTypePtr != LbPfEnComment))) {
			sprintf(errString, "Incomplete parameter specification on line %ld.",
				(unsigned long)lbPfParamFileCurLine);
			ErStGeneric(errString);
			break;
		}
		
		/* If its not a comment, continue to parse the line */
		if (*paramTypePtr != LbPfEnComment) {
		
			/* Copy the label into the label buffer */
			strncpy(labelPtr, tokenBuff, LBPF_LABEL_LEN);			
			
			/* Get the equal sign */
			if (!lbPfGetToken(paramBufferPtr, tokenBuff, &currentChar,
					&tokenLen, &isLastToken)) {
					
				break;
			}
			
			/* If its the last token, or the token is empty, or its not an equal sign, its an error */
			if ((tokenLen == 0) || (isLastToken == true) ||
					(tokenBuff[0] != '=')) {
				sprintf(errString, "Parameter label not followed by an equal sign on line %ld.",
					(unsigned long)lbPfParamFileCurLine);
				ErStGeneric(errString);
				break;
			}
			
			/* Get the parameter value */
			if (!lbPfGetToken(paramBufferPtr, tokenBuff, &currentChar,
					&tokenLen, &isLastToken)) {
					
				break;
			}
			
			/* If the token is empty, its an error */
			if (tokenLen == 0) {
				sprintf(errString, "Null token provided for parameter length on line %ld.",
					(unsigned long)lbPfParamFileCurLine);
				ErStGeneric(errString);
				break;
			}
			
		}
		
		/* Finish up now that you have all of the fields */
		switchOkay = true;
		switch (*paramTypePtr) {
		
			case LbPfEnBoolean:
				
				/* Set true or false value */
				if (toupper((int)tokenBuff[0]) == 'T'){
					*((Boolean *) paramBufferPtr) = true;
				}
				else if (toupper((int)tokenBuff[0]) == 'F'){
					*((Boolean *) paramBufferPtr) = false;
				}
				else {
					/* Flag boolean error */
					switchOkay = false;
					sprintf(errString, "Err - Line %ld: BOOL parameters must be assigned TRUE or FALSE, not [%s].",
						(unsigned long)lbPfParamFileCurLine, tokenBuff);
					ErStGeneric(errString);
					break;
				}

					
				*paramSizePtr = sizeof(Boolean);
				
				break;
				
			case LbPfEnComment:
			
				/* Just set the length */
				*paramSizePtr = tokenLen;
				break;
				
			case LbPfEnChar:
				/* Copy the char */
				paramBufferPtr[0] = tokenBuff[0];
				
				*paramSizePtr = 1;
				
				break;
				
			case LbPfEnString:
				
				/* Verify that strings is surrounded by "" */
				if ((tokenBuff[0] != '"') ||
						(tokenBuff[tokenLen-1] != '"')){

					switchOkay = false;
					sprintf(errString, "Err - line %ld: String parameters must be enclosed in quotes.\n\t[%s]",
						(unsigned long)lbPfParamFileCurLine, tokenBuff);
					ErStGeneric(errString);
					break;
				}

				stringIndex = 0;
				do {
					if (tokenBuff[1+stringIndex] != '"')
						stringIndex++;
						
				} while ((tokenBuff[1+stringIndex] != '"')
					&& (stringIndex+2 <= tokenLen));
				
				if (stringIndex == 0) {
					paramBufferPtr[0] = '\0';
				}
				else {
					/* Copy the string */
					memcpy(paramBufferPtr, tokenBuff+1, stringIndex);
					paramBufferPtr[stringIndex] = '\0';
				}
				
				/* Set the length */
				*paramSizePtr = stringIndex;
				
				break;
				
			case LbPfEnEnum:
				
				/* Copy the string */
				memcpy(paramBufferPtr, tokenBuff, tokenLen);
				paramBufferPtr[tokenLen] = '\0';
				
				/* Set the length */
				*paramSizePtr = tokenLen;
				
				break;
				
			case LbPfEnInteger:
				
				/* Verify that they provided numeric information.
					Note that this is somewhat problematic, we'll
					assume that if the last char is numeric then
					the parameter is numeric.  This is necessary
					because of the potential for using the '-'
					sign for entering negative values.
				 */
				if (isdigit(tokenBuff[tokenLen-1]) == 0) {
					switchOkay = false;
					sprintf(errString, "You specified an 'INTEGER' type for"
				 	" parameter '%s'  but supplied a non-numeric string '%s'\n",
				 	labelPtr, tokenBuff);
					ErStGeneric(errString);
					break;
				}
				
				/* Convert to integer */
				*((LbUsFourByte *) paramBufferPtr) = atoi(tokenBuff);
				
				/* Set the length */
				*paramSizePtr = sizeof(LbUsFourByte);
				
				break;
				
			case LbPfEnReal:
				
				
				/* Convert to real */
				dblValue = atof(tokenBuff);
				
				*((double *) paramBufferPtr) = dblValue;
				/* Set the length */
				*paramSizePtr = sizeof(double);
				
				break;
				
			case LbPfEnList:
				
				/* Store the number of elements in the list */
				*((LbUsFourByte *) paramBufferPtr) = atoi(tokenBuff);
				
				/* Set the length */
				*paramSizePtr = sizeof(LbUsFourByte);
				
				break;
				
			case LbPfEnLongLong:
				
				/* Verify that they provided numeric information.
					Note that this is somewhat problematic, we'll
					assume that if the last char is numeric then
					the parameter is numeric.  This is necessary
					because of the potential for using the '-'
					sign for entering negative values.
				 */
				if (isdigit(tokenBuff[tokenLen-1]) == 0) {
					switchOkay = false;
					sprintf(errString, "You specified an 'INTEGER' type for"
				 	" parameter '%s'  but supplied a non-numeric string '%s'\n",
				 	labelPtr, tokenBuff);
					ErStGeneric(errString);
					break;
				}
				
				/* Convert to integer */
				*((LbUsEightByte *) paramBufferPtr) = strtoll(tokenBuff, (char**)NULL, 10);
				
				/* Set the length */
				*paramSizePtr = sizeof(LbUsEightByte);
				
				break;
				
			}
			/* Break if switch not okay */
			if (!switchOkay)
				break;

		okay = true;
	} while (false);
	
	return (okay);
}

/*****************************************
*		lbPfGetToken
*
*	Purpose:	Get the next token off the line.
*	Arguments:
*		char			*paramBufferPtr			- The parameter value buffer.
*		char			*tokenPtr				- The token buffer.
*		LbUsFourByte	*currentCharIndexPtr	- The current char in the parameter buffer.
*		LbUsFourByte	*tokenLenPtr			- The length of the token.
*		Boolean			*isLastTokenPtr			- Is this the last token on the line.
*
*	Returns:	True unless an error or end of file occurs.
*
******************************************/
Boolean	lbPfGetToken(char *paramBufferPtr,
			char *tokenPtr, LbUsFourByte *currentCharIndexPtr,
			LbUsFourByte *tokenLenPtr,
			Boolean *isLastTokenPtr)
{
	Boolean		okay = false;			/* Process flag */
	char		errString[1024];		/* Buffer for error strings */
	
	do { /* Process Loop */
		
		/* Clear variable flags */
		*isLastTokenPtr = false;
		*tokenLenPtr = 0;
		tokenPtr[0] = '\0';
		
		/* Eat the white space up to the first token */
		while (((paramBufferPtr[*currentCharIndexPtr] == ' ') ||
				(paramBufferPtr[*currentCharIndexPtr] == '\t') ||
				(paramBufferPtr[*currentCharIndexPtr] == '\n')) &&
				(*currentCharIndexPtr < LBPF_PARAM_LEN_MAX)) {
				
			(*currentCharIndexPtr)++;
		}
		
		/* If we reached the end of the buffer without a token it is an error */
		if (*currentCharIndexPtr == LBPF_PARAM_LEN_MAX) {
			ErStGeneric("No token found in parameter file.");
			break;
		}
		
		/* Get the next token */
		okay = LbPfGetNextToken( paramBufferPtr, tokenPtr, 
									currentCharIndexPtr, tokenLenPtr, isLastTokenPtr);
		
		/* If we went past our line length, we are in an error situation */
		if (*currentCharIndexPtr == LBPF_PARAM_LEN_MAX) {
			sprintf(errString, "Err line %ld: Token too long in parameter file.",
				(unsigned long)lbPfParamFileCurLine);
			ErStGeneric(errString);
			break;
		}
		
	} while (false);
	
	return (okay);
}

/*****************************************
*		lbPfLineIsBlank
*
*	Purpose:	Test line for non white-space characters.
*	Arguments:
*		char			*paramBufferPtr			- The string we are checking.
*	Returns:	True if the line contains only white-space.
*
******************************************/
Boolean	lbPfLineIsBlank(char *paramBufferPtr)
{
	Boolean			lineIsBlank = true;		/* Process flag */
	LbUsFourByte	charIndex;				/* Current char */
	
	
	/* Loop through characters until non white-space, or end of line */
	charIndex = 0;
	while((paramBufferPtr[charIndex] != '\0') &&
		(paramBufferPtr[charIndex] != '\n')){ /* String search */
	
		/* See if cur char is any non white-space char */
		if (isspace(paramBufferPtr[charIndex]) == false) {
			lineIsBlank = false;
			break;
		}
		
		charIndex++;
	}
	
	return (lineIsBlank);
}
