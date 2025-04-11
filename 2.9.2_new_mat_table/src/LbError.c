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
*     Module Name:        LbError.c
*     Revision Number:    1.1
*     Date last revised:  23 July 2013
*     Programmer:         Steven Vannoy
*     Date Originated:    Tuesday, July 21, 1992
*
*     Module Overview:    This module provides routines for handling errors.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*     	ErAbort
*		ErAlert
*		ErIgnoreBegin
*		ErIgnoreEnd
*		ErInit
*		ErClear
*		ErClearIf
*		ErClearIfManager
*		ErHandle
*		ErNewMgCode
*		ErIsInError
*		ErStFileError
*		ErStPrimError
*		ErTerminate
*
*     Global variables defined:   
*
*********************************************************************************/
#define LB_ERROR

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbDebug.h"
#include "LbInterface.h"

/*	CONSTANTS */

/*  LOCAL GLOBALS */
char			erCurrentErStr[ERMxErStrLen];
Boolean			erIsInError = false;
Boolean			erIsIgnoring = false;
ErMgCdTy		erCurrentMgCode = ERMgCdNull;
ErErCdTy		erCurrentErCode = ERErCdNull;
ErMgCdTy		erNextErMgCode = ERMgCdGeneric + 1;
FILE			*erOutputFile = 0;
LbUsFourByte	erInitFlags;

#ifdef LB_DEBUG
	Boolean	erIsInited = false;
#endif

/*	LOCAL MACROS */
/**********************
*	ERFgEcho
*
*	Arguments:	erInitFlags	- Our initialization flag.
*
*	Purpose:	See if we are echoing error messages to stdout.
*
*	Result:	TRUE if echo flag set.
***********************/
#define Is_EchoStdOut(erInitFlags)	LbFgIsSet((erInitFlags), ERFgEcho)

/*****************************************
*		ErAbort
*
*	Arguments:
*				char			*errorString	- String to display.
*	Purpose:	Display an error message to the user.
*				ABORT the program!!!
*
*	Returns:	None.
*
******************************************/
void	ErAbort(char *errorString)
{
	Boolean	continueDebug;	/* Should we continue (debugging only) */
	
	do { /* Process Loop */
	
		#ifdef LB_DEBUG
			/* Verify manager is initialized */
			if (!erIsInited) {
	
				LbInPrintf("\nYou must initialize error manager before using!");
				exit(0);
			}
			
		#endif
		
		/* Write error message to output file */
		fprintf(erOutputFile, "\n%s\n", errorString);
		
		/* Write error messages to output file */
		if (erIsInError) {
			fprintf(erOutputFile, "\n\tBecause '%s'\n", erCurrentErStr);
		
			/* Cear all error parameters */
			ErClear();
		}
		
		/* Offer to debug if we are dealing with standard output */
		if (erOutputFile == stdout) {
			continueDebug = false;
			LbDbOffer(&continueDebug);
			
			/* If they want to continue, break out of the loop and avoid the exit */
			if (continueDebug)
				break;
		}
		
		/* Tell them we are leaving now */
		fprintf(erOutputFile, "\nProgram is terminating.");
		
		/* See if echoing to stdout */
		if (Is_EchoStdOut(erInitFlags))
			LbInPrintf("\n%s", errorString);
	
		/* If we are on the mac we "throw" the error */
		#ifdef __USE_TWRP__
			LbWpExit();
		#endif
		
		/* EXIT THE PROGRAM!!!!!!*/
		exit(0);
	} while (false);
}

/*****************************************
*		ErBreak
*
*	Arguments:
*				char			*errorString	- String to display.
*	Purpose:	Display an error message to the user.
*				Allow aborting or continuing.
*
*	Returns:	None.
*
******************************************/
void	ErBreak(char *errorString)
{
	Boolean	continueDebug;	/* Should we continue to debug */
	
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
		
	#endif
	
	/* Write error message to output file */
	fprintf(erOutputFile, "\nBreakpoint:\n\t%s\n", errorString);
	
	/* Offer to debug if we are dealing with standard output */
	if (erOutputFile == stdout)
		LbDbOffer(&continueDebug);
}

/*****************************************
*		ErAlert
*
*	Arguments:
*				char			*errorString	- String to display.
*				Boolean			getResponse		- Flag for getting response.
*	Purpose:	Display an error message to the user.
*
*	Returns:	None.
*
******************************************/
void	ErAlert(char *errorString, Boolean getResponse)
{

	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
	#endif
	
	/* Write error message to output file */
	fprintf(erOutputFile, "\n%s", errorString);
	
	/* See if echoing to stdout */
	if (Is_EchoStdOut(erInitFlags))
		LbInPrintf("\n%s", errorString);
		
	/* See if getting response from user */
	if (getResponse) {
		/* Write prompt */
		LbInPrintf("\nPlease hit return to acknowledge error.");
		
		/* Get the input */
		(void) getc(stdin);
	}
}

/*****************************************
*		ErIgnoreBegin
*
*	Arguments:
*	Purpose:	Begin ignoring error registration.
*
*	Returns:	None.
*
******************************************/
void	ErIgnoreBegin()
{

	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
		
		/* Verify we are in an error situation */
		if (!erIsInError) {

			fprintf(erOutputFile, "\nYou cannot call ErIgnoreBegin when not in error!");
			exit(0);
		}
	#endif
	
	/* Remember we are ignoring */
	erIsIgnoring = true;
}

/*****************************************
*		ErIgnoreEnd
*
*	Arguments:
*	Purpose:	End ignoring error registration.
*
*	Returns:	None.
*
******************************************/
void	ErIgnoreEnd()
{
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
		
		/* Make sure we are ignoring */
		if (!erIsIgnoring) {

			fprintf(erOutputFile, "\nYou must be ignoring to end ignoring!");
			exit(0);
		}
	#endif

	/* Clear our ignoring flag */
	erIsIgnoring = false;
}

/*****************************************
*		ErInit
*
*	Arguments:
*				LbUsFourByte	initFlags 	- Initialization flags.
*				FILE			*outputFile	- User's output file.
*	Purpose:	Initialize the manager.
*
*	Returns:	TRUE if initialization succeeds.
*
******************************************/
Boolean	ErInit(LbUsFourByte initFlags, FILE *outputFile)
{
	Boolean okay = false;
	
	#ifdef LB_DEBUG
		/* Verify we are not doing this twice */
		if (erIsInited) {
			LbInPrintf("\nYou can't call ErInit twice.");
			exit(0);
		}		
		/* If they gave us an output file, try to verify it */
		if (outputFile) {
			/* Just do an ftell to get current position */
			if (ftell(outputFile) == -1) {
				LbInPrintf("\nUnable to use given output file for ErInit.");
				exit(0);
			}
		}
		else {
			/* Make sure they don't say echo to user since that's the default */
			if (Is_EchoStdOut(initFlags)) {
				LbInPrintf("\nYou cannot set echo flag without supplying file.");
				exit(0);
			}
		}
	#endif

	do	/* Abort Loop */
	{
		/* Clear current error string */
		erCurrentErStr[0] = 0;

		/* Set our output file */
		if (outputFile)
			erOutputFile = outputFile;
		else
			erOutputFile = stdout;
		
		/* Save their init flags */
		erInitFlags = initFlags;
		
		/* Indicate we've been initialized */
		#ifdef LB_DEBUG
			erIsInited = true;
			
		#endif
		
		okay = true;
	} while (false);
	
	return (okay);
}

/*****************************************
*		ErClear
*
*	Arguments:
*	Purpose:	Clear the current error.
*
*	Returns:	None.
*
******************************************/
void	ErClear()
{
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
		
		/* Make sure we are in error condition */
		if (!erIsInError) {

			fprintf(erOutputFile, "\nYou must be in error condition to clear!");
			exit(0);
		}
	#endif
	
	/* Clear the error condtion */
	erIsInError = false;
	erCurrentMgCode = ERMgCdNull;
	erCurrentErCode = ERErCdNull;
	erCurrentErStr[0] = 0;
}

/*****************************************
*		ErClearIf
*
*	Arguments:
*				ErMgCdTy	managerCode	- Manager code to match.
*				ErErCdTy	errorCode	- Error code to match.
*				Boolean		*cleared		- TRUE if error matched and was cleared.
*	Purpose:	Clear the current error it it matches the
*				manager/error code.
*
*	Returns:	None.
*
******************************************/
void	ErClearIf(ErMgCdTy managerCode, ErErCdTy errorCode, Boolean *cleared)
{
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
		
		/* Make sure we are in error condition */
		if (!erIsInError) {

			fprintf(erOutputFile, "\nYou must be in error condition to clear!");
			exit(0);
		}
	#endif
	
	/* See if manager and code match */
	if ((managerCode == erCurrentMgCode) && (errorCode == erCurrentErCode)) {

		/* Clear the current error */
		ErClear();
		*cleared = true;
	}
	else
		*cleared = false;
}

/*****************************************
*		ErClearIfManager
*
*	Arguments:
*				ErMgCdTy	managerCode		- Manager code to match.
*				Boolean		*cleared		- TRUE if error matched and was cleared.
*	Purpose:	Clear the current error it it matches the
*				manager code.
*
*	Returns:	None.
*
******************************************/
void	ErClearIfManager(ErMgCdTy managerCode, Boolean *cleared)
{
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
		
		/* Make sure we are in error condition */
		if (!erIsInError) {

			fprintf(erOutputFile, "\nYou must be in error condition to clear!");
			exit(0);
		}
	#endif
	
	/* See if manager and code match */
	if (managerCode == erCurrentMgCode) {

		/* Clear the current error */
		ErClear();
		*cleared = true;
	}
	else
		*cleared = false;
}

/*****************************************
*		ErHandle
*
*	Arguments:
*				char	*errorString	- High level string to display.
*				Boolean	getResponse		- Flag to require user response.
*	Purpose:	Handle the current error by displayint high/low
*				level messages to the user and clearing the
*				error condition.
*
*	Returns:	None.
*
******************************************/
void	ErHandle(char *errorString, Boolean getResponse)
{
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
		
		/* Make sure we are in an error condition */
		if (!erIsInError) {

			LbInPrintf("\nYou must be in an error to handle one!\n");
			exit(0);
		}
		
		/* Make sure we are not ignoring errors */
		if (erIsIgnoring) {

			LbInPrintf("\nYou cannot handle an error while ignoring!\n");
			exit(0);
		}
	#endif
	
	/* Write error messages to output file */
	fprintf(erOutputFile, "\n%s\nBecause '%s'\n", errorString, erCurrentErStr);
	
	/* See if echoing to stdout or requiring response */
	if (Is_EchoStdOut(erInitFlags) || getResponse)
		LbInPrintf("\n%s\n\t%s\n", errorString, erCurrentErStr);
		
	/* See if getting response from user */
	if (getResponse) {
		/* Write prompt */
		LbInPrintf("\nPlease hit return to acknowledge error.");
		
		/* Get the input */
		(void) getc(stdin);
	}
	
	/* Cear all error parameters */
	ErClear();
}

/*****************************************
*		ErNewMgCode
*
*	Arguments:
*				ErMgCdTy	*newMgCode	- Storage for new code.
*	Purpose:	Return a new unique manager code.
*
*	Returns:	None.
*
******************************************/
void	ErNewMgCode(ErMgCdTy *newMgCode)
{
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
	#endif
	
	/* Give them next available manager code */
	*newMgCode = erNextErMgCode;
	
	/* Increment next available manager */
	erNextErMgCode++;
}

/*****************************************
*		ErWhatError
*
*	Arguments:
*				Boolean		*isInError		- Error condition flag.
*				ErMgCdTy	*managerCode	- Manager of current error.
*				ErErCdTy	*errorCode		- Code for current error.
*				char		*errorStr		- Error string.
*	Purpose:	Get current error information.
*
*	Returns:	None.
*
******************************************/
void	ErWhatError(Boolean *isInError, ErMgCdTy *managerCode,
				ErErCdTy *errorCode, char *errorStr)
{
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
	#endif
	
	/* See if we are in an error condition */
	if (erIsInError) {

		/* Copy parameters */
		{
			*isInError = true;
			*managerCode = erCurrentMgCode;
			*errorCode = erCurrentErCode;
			
			/* If they gave us string space, copy it */
			if (errorStr)
				strcpy(errorStr, erCurrentErStr);
		}
			
	}
	else
	{
		/* We are not in an error condition */
		*isInError = false;
		*managerCode = ERMgCdNull;
		*errorCode = ERErCdNull;
	}
}
/*****************************************
*		ErIsInError
*
*	Arguments:
*	Purpose:	Indicate if we are in an error
*				situation.
*
*	Returns:	True if low level error is set.
*
******************************************/
Boolean	ErIsInError()
{
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
	#endif
	return(erIsInError);
}

/*****************************************
*		ErStFileError
*
*	Arguments:
*				char		*errorString	- Primitive error string.
*	Purpose:	Sets a primitive error after a file operation has
*				failed. The error string is concatenated with the
*				error string associated with errno.
*
*	Returns:	None.
*
******************************************/
void	ErStFileError(char *errorString)
{
	#ifdef LB_DEBUG
		Boolean	continueDebug;		/* Should we continue to debug */
	#endif
	
	char	bigBuffer[1028];	/* Big buffer for error string */
	char	*fileErStr;			/* Pointer for storing error string from file system */
	
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
		
		/* Make sure we are not in an error condtion and not ignoring */
		if (erIsInError && !erIsIgnoring) {

			fprintf(erOutputFile, "\nYou can't register error messages twice!");
			exit(0);
		}
		
		/* Make sure they gave us a string */
		if (strlen(errorString) == 0) {

			/* Enter the debugger */
			LbDbEnter(&continueDebug);

			LbInPrintf("\nYou must supply a string for primitive registration!");
			exit(0);
		}
		
	#endif
	
	/* Concatenate the messages if not SUN_OS */
	#ifndef SUN_OS
		fileErStr = strerror(errno);
		sprintf(bigBuffer, "%s, %s", errorString, fileErStr);
	#else
		sprintf(bigBuffer, "%s, error number = %d", errorString, errno);
	#endif
	
	/* Register error as generic */
	ErStGeneric(bigBuffer);
}

/*****************************************
*		ErStPrimError
*
*	Arguments:
*				ErMgCdTy	managerCode		- Manager. of error.
*				ErErCdTy	errorCode		- Code of error.
*				char		*errorString	- Primitive error string.
*	Purpose:	Set the current error condition.
*
*	Returns:	None.
*
******************************************/
void	ErStPrimError(ErMgCdTy managerCode, ErErCdTy errorCode,
			char *errorString)
{
	#ifdef LB_DEBUG
		Boolean	continueDebug;	/* Should we continue to debug */
		char	message[1024];	/* For building programmer message */
	#endif
	
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
		
		/* Make sure we are not in an error condtion and not ignoring */
		if (erIsInError && !erIsIgnoring) {

			sprintf(message, "Dear Programmer,\n\tYou can not register two primitive error"
				" messages at one time.\n\tYou have registered the following message:\n\t\t '%s'\n"
				"\twithout handling the error.\n\tYou are currently trying to register this message:\n\t"
				" '%s'.\n\n", erCurrentErStr, errorString);
			fprintf(erOutputFile, message);
			exit(0);
		}
		
		/* Make sure they gave us a string */
		if (strlen(errorString) == 0) {

			/* Enter the debugger */
			if (errorCode != ERErCdUserCancel)
				LbDbEnter(&continueDebug);

			LbInPrintf("\nYou must supply a string for primitive registration!");
			exit(0);
		}
		
	#endif
	
	/* If not ignoring, register the error */
	if (!erIsIgnoring) {
		erIsInError = true;
		erCurrentMgCode = managerCode;
		erCurrentErCode = errorCode;
		strncpy(erCurrentErStr, errorString, ERMxErStrLen);
	}

}

/*****************************************
*		ErTerminate
*
*	Arguments:
*	Purpose:	Terminate the error manager.
*
*	Returns:	None.
*
******************************************/
void	ErTerminate()
{
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!erIsInited) {

			LbInPrintf("\nYou must initialize error manager before using!");
			exit(0);
		}
		
		/* Make sure we are not terminating with dangling error */
		if (erIsInError) {

			LbInPrintf("\nYou must clear errors before terminating!");
			exit(0);
		}

		/* Terminate the manager */
		erIsInited = false;
	#endif
	
}
#undef LB_ERROR
