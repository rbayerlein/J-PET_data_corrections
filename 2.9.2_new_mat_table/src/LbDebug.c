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
*     Module Name:        LbDebug.c
*     Revision Number:    1.2
*     Date last revised:  6 June 2013
*     Programmer:         Steven Vannoy
*     Date Originated:    Tuesday, July 21, 1992
*
*     Module Overview:	This module provides routines for debuging programs.
*						MACHINE DEPENDANT						
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*     	LbDbEnter
*		LbDbInit
*		LbDbOffer
*		LbDbTerminate
*
*     Global variables defined:   
*
*********************************************************************************/

#include	<stdio.h>
#include	<string.h>

#include "SystemDependent.h"

#include	"LbTypes.h"
#include	"LbDebug.h"
#include	"LbInterface.h"

/*	CONSTANTS */

/*  LOCAL GLOBALS */
Boolean			lbDbIsInited = false;
char			lbDbProgName[255];
FILE			*lbDbOutputFile = 0;
LbUsFourByte	lbDbInitFlags = 0;

/*	LOCAL MACROS */
/*	LOCAL FUNCTIONS */
/*	LOCAL MACROS */

/*	FUNCTIONS */
#ifdef DGUX
/*****************************************
*		LbDbEnter
*
*	Purpose:	Enter the debugger.
*	Arguments:
*		Boolean	*continuePtr - Should we continue to debug.
*
*	Returns:	None.
*
******************************************/
void	LbDbEnter(Boolean *continuePtr)
{
	#ifdef LB_DEBUG
		extern					int errno;
		char					parent_pid_string[32];
		int 					pid, parent_pid;
		long					key;
		struct	dg_process_info	procInfo;
#endif /* LB_DEBUG */
	
	do { /* Process Loop */
	#ifdef LB_DEBUG
	
		/* Assume we are continuing (Programmer can set to false via the debugger to change value) */
		*continuePtr = true;
		
		/* Verify we've been initialized */
		if (!lbDbIsInited) {
			fprintf(lbDbOutputFile, "\nYou must initialize the debug manager to call LbDbEnter.");
			break;
		}
		
		/* Check our parent to see if we are already in the debugger */
		parent_pid = getppid();
		
		key = DG_PROCESS_INFO_INITIAL_KEY;
		
		/* Get info on parent */
		if (errno = dg_process_info(DG_PROCESS_INFO_SELECTOR_PID,
					parent_pid, DG_PROCESS_INFO_CMD_NAME_ONLY,
				&key, &procInfo, DG_PROCESS_INFO_CURRENT_VERSION) != 1) {
				
			fprintf(stderr, "Can't get info on parent so can't attach to debugger.\n");
			break;
		}
		
		/* See if parent is mxdb */
		if (strcmp(DEBUGGER_NAME, procInfo.cmd) == 0){

			#ifdef WHEN_IT_WORKS
				/* At this point we have to kludge to get to mxdb */
				parent_pid /= 0;
			#endif

			/* We are a son of mxdb so send a signal to stop (Note to user, you must tell mxdb to cat the usr1 signal) */
			if (kill(parent_pid, SIGUSR1) != 0) {
				fprintf(stderr, "Can't call kill so can't attach to debugger.\n");
				break;
			}
				
		}
		else {
			/* We are not a son of mxdb so do this (This code copied from Mxdb manual) */
			pid = vfork();
			if (pid == 0) {
				parent_pid = getppid();
				(void) execl(DEBUGGER_PATHNAME, DEBUGGER_NAME, "-p",
					itoa(parent_pid, parent_pid_string),
					lbDbProgName, (char *) 0);
				(void) fprintf(stderr, "Exec of \"%s\" failed.  Errno is %d.\n",
					DEBUGGER_PATHNAME, errno);
					_exit(255);
			} 
			else if (pid == -1) {
				fprintf(stderr, "Can't vfork, so can't attach to debugger.\n");
				break;
			}
			
			/* Sleep for a few seconds to give the debugger time to attach */
			sleep(10);
		}
	#endif
	} while (false);
}
#endif /* DGUX */

#ifdef AOS_VS
/*****************************************
*		LbDbEnter
*
*	Arguments:
*	Purpose:	Enter the debugger.
*
*	Returns:	None.
*
******************************************/
void	LbDbEnter()
{
	int	errNum;
	
	do { /* Process Loop */
	#ifdef LB_DEBUG
		/* Verify we've been initialized */
		if (!lbDbIsInited) {
			fprintf(lbDbOutputFile, "\nYou must initialize the debug manager to call LbDbEnter.");
			break;
		}

		/* Call swat */
		if ((errNum = call_swat("@CONSOLE")) != 0)
			fprintf(lbDbOutputFile, "\nError trying to enter swat = %d\n", errNum);

	#endif
	} while (false);
}
#endif /* AOS_VS */

#ifdef MPW
/*****************************************
*		LbDbEnter
*
*	Arguments:
*		Boolean			*continuePtr	- Flag to continue.
*
*	Purpose:	Enter the debugger.
*
*	Returns:	None.
*
******************************************/
void	LbDbEnter(Boolean *continuePtr)
{
	DebugStr ((ConstStr255Param) "\pEntering Debugger...");
}
#endif /* MPW */

#ifdef __MWERKS__
/*****************************************
*		LbDbEnter
*
*	Arguments:
*		Boolean			*continuePtr	- Flag to continue.
*
*	Purpose:	Enter the debugger.
*
*	Returns:	None.
*
******************************************/
void	LbDbEnter(Boolean *continuePtr)
{
	if (*continuePtr) {};	/* Avoid unused parameter compiler warning */
	
	DebugStr ((ConstStr255Param) "\pEntering Debugger...");
}
#endif /* __MWERKS__ */

#ifndef DGUX
#ifndef AOS_VS
#ifndef MPW
#ifndef __MWERKS__
/*****************************************
*		LbDbEnter
*
*	Arguments:
*		Boolean			*continuePtr	- Flag to continue.
*	Purpose:	Enter the debugger.
*
*	Returns:	None.
*
******************************************/
void	LbDbEnter(Boolean *continuePtr)
{
	LbInPrintf("\nDebugger not available");
	*continuePtr = false;
}
#endif 
#endif 
#endif
#endif

/*****************************************
*		LbDbInit
*
*	Arguments:
*				char			*progName	- Name of program being debuged.
*				LbUsFourByte	initFlags 	- Initialization flags.
*				FILE			*outputFile	- User's output file.
*	Purpose:	Initialize the manager.
*
*	Returns:	TRUE if initialization succeeds.
*
******************************************/
Boolean	LbDbInit(char *progName, LbUsFourByte initFlags, FILE *outputFile)
{
	Boolean okay = false;
	
	/* Do nothing if not debugging or already initialized */
	#ifndef LB_DEBUG
		return(true);
	#endif
	

	do	/* Abort Loop */
	{
		/* If initializing twice, we are done */
		if (lbDbIsInited)
		{
			goto SUCCESS;
		}
		
		/* If they gave us an output file, try to verify it */
		if (outputFile)
		{
			/* Just do an ftell to get current position */
			if (ftell(outputFile) == -1)
			{
				LbInPrintf("\nUnable to use given output file for LbDbInit.");
				goto FAIL;
			}

			lbDbOutputFile = outputFile;
		}
		else {
		
			/* Set our output file */
			lbDbOutputFile = stdout;
		}

		/* Save their init flags */
		lbDbInitFlags = initFlags;
		
		/* Copy over the program name */
		strcpy(lbDbProgName, progName);
		
		/* Indicate we've been initialized */
		lbDbIsInited = true;
			
		SUCCESS:	okay = true;
		FAIL:;
	} while (false);
	
	return (okay);
}

/*****************************************
*		LbDbOffer
*
*	Purpose:	Allow the user to enter
*			 The debugger if they wish.
*	Arguments:
*		Boolean	*continuePtr	- Should we continue to debug.
*	Returns:	None.
*
******************************************/
void	LbDbOffer(Boolean *continuePtr)
{
	#ifdef LB_DEBUG
		Boolean		canceled;	/* User input cancelation flag */
		LbTwoByte	response;	/* User response */
	#endif
	
	if (*continuePtr) {};		/* Avoid unused parameter compiler warning */
	
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!lbDbIsInited)
		{
			LbInPrintf("\nYou must initialize debug manager before using!\n");
			exit(0);
		}
		
		/* Pop the question */	
		response = LbInAsk("Select an action: ", 3, true,
					&canceled, "continue", "debugger", "quit", 0,
					0);
					
		/* Debug if they want */
		if (response == 2)
			LbDbEnter(continuePtr);
		else if (response == 3) {
			LbInPrintf("Program terminating.\n");
			exit (0);
		}
	#endif
}

/*****************************************
*		LbDbTerminate
*
*	Arguments:
*	Purpose:	Terminate the debug manager.
*
*	Returns:	None.
*
******************************************/
void	LbDbTerminate()
{
	#ifdef LB_DEBUG
		/* Verify manager is initialized */
		if (!lbDbIsInited)
		{
			LbInPrintf("\nYou must initialize debug manager before using!");
			exit(0);
		}
		
		/* Terminate the manager */
		lbDbIsInited = false;
	#endif
	
}
