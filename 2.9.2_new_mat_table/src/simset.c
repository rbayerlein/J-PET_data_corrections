/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1995-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		simset.c
*			Revision Number:	1.17
*			Date last revised:	4 June 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	25 January 1995
*
*			Module Overview:	This is the "main" module for the SimSET software
*								package.
*
*			References:			None.
*
**********************************************************************************
*
*			Global functions defined:
*
*			Global macros defined:
*				
*			Global variables defined:		none
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
*			Revision date:		23 October 2012
*
*			Revision description:	Moved contents of main to new runSimSET function.
*									Made main dependent on OS (for GUIOS systems).
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		22 August 2012
*
*			Revision description:	Changed getArgs to use LbPfGetNextToken function.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		16 April 2012
*
*			Revision description:	Checked for empty command string (non-Unix).
*									Added local command string retrieval (non-Unix).
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		September 2005
*
*			Revision description:	Added resampledecaytime utility (resampt)
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		September 2005
*
*			Revision description:	Added addRandoms utility (adrand)
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		11 August 2005
*
*			Revision description:	Added timesort utility (tmsort)
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		25 November 2003
*
*			Revision description:	Added changes for Mac CW 9 compatibility:
*										Modified changeDir for modern Mac OS
*
*********************************************************************************/

#define UTIL_MAIN			

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>

#ifdef __MWERKS__
	#include <Carbon.h>
#endif

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"

#ifdef MPW
	#pragma segment SIMSET_MAIN
#endif


/* LOCAL CONSTANTS */


/* LOCAL TYPES */

/* The following type "func" and the table "SimsetFuncsTy" are used
	to create a jump table of functions.
*/
typedef Boolean (*func)(int argc, char *argv[]);

/* This struct is used to create a labeled function table */
typedef struct {
	char *funcName;		/* Name of the function */
	func funcPtr;		/* Address of the function */
} SimsetFuncElemTy;


/* PROTOTYPES */
#ifndef COCOA
int				main(int argc, char *argv[]);
#endif
int				runSimSET(int argc, char *argv[]);
Boolean			phgbin(int argc, char *argv[]);
Boolean			Convert(int argc, char *argv[]);
Boolean			BuildAtt(int argc, char *argv[]);
Boolean			CombineBin(int argc, char *argv[]);
Boolean			CombineHist(int argc, char *argv[]);
Boolean			ReverseBytes(int argc, char *argv[]);
Boolean			DisplayHeader(int argc, char *argv[]);
Boolean			makeindexfile(int argc, char *argv[]);
Boolean			ttest(int argc, char *argv[]);
Boolean			stripheader(int argc, char *argv[]);
Boolean			migrate(int argc, char *argv[]);
Boolean			PhgSwap(int argc, char *argv[]);
Boolean			PrintHeader(int argc, char *argv[]);
Boolean			ConvertHeader(int argc, char *argv[]);
Boolean			CalcAttenuation(int argc, char *argv[]);
Boolean			AttnCorrect(int argc, char *argv[]);
Boolean			Collapse(int argc, char *argv[]);
Boolean			Extract(int argc, char *argv[]);
Boolean			Collapse3d(int argc, char *argv[]);
Boolean			PhgRun(int argc, char *argv[]);
Boolean			bcomp(int argc, char *argv[]);
Boolean			BreakpointSwap(int argc, char *argv[]);
Boolean			ExtractLines(int argc, char *argv[]);
Boolean			BuildCoh(int argc, char *argv[]);
Boolean			line3d(int argc, char *argv[]);
Boolean			Scale(int argc, char *argv[]);
Boolean			reorder(int argc, char *argv[]);
Boolean			convertCoh(int argc, char *argv[]);
Boolean			changeDir(int argc, char *argv[]);
Boolean			tmsort(int argc, char *argv[]);
Boolean			adrand(int argc, char *argv[]);
Boolean			resampt(int argc, char *argv[]);

#ifdef GUIOS
int getArgs(char ***argvPtr);
#endif
void setup_args(int *argcPtr, char *argv[], char *progName);


/* LOCAL GLOBALS */

char	errString[256];	/* Our standard error string */

/*
	The following table is list of functions and their names.
	It associates a string with each "module" of the simset
	package.  It is no particular order, but if you ad to 
	this table you must increase the constant NUM_FUNCS
	and you must always be sure that PHG_INDEX refers
	to the function PhgRun
*/
#define NUM_FUNCS	31
#define PHG_INDEX	0
SimsetFuncElemTy SimsetFuncTbl[] = {

	{"phg",					PhgRun},	/* Don't put another function here unless you change PHG_INDEX, rather add to bottom of the list */
	{"phgbin",				phgbin},
	{"convert",				Convert},
	{"buildatt",			BuildAtt},
	{"combinebin",			CombineBin},
	{"combinehist",			CombineHist},
	{"reversebytes",		ReverseBytes},
	{"displayheader",		DisplayHeader},
	{"makeindexfile",		makeindexfile},
	{"ttest",				ttest},
	{"stripheader",			stripheader},
	{"migrate",				migrate},
	{"phgswap",				PhgSwap},
	{"printheader",			PrintHeader},
	{"convertheader",		ConvertHeader},
	{"calcattenuation",		CalcAttenuation},
	{"attncorrect",			AttnCorrect},
	{"collapse",			Collapse},
	{"xdata",				Extract},
	{"collapse3d",			Collapse3d},
	{"bcomp",				bcomp},
	{"breakpointswap",		BreakpointSwap},
	{"extractlines",		ExtractLines},
	{"buildcoh",			BuildCoh},
	{"line3d",				line3d},
	{"scale",				Scale},
	{"cd",					changeDir},
	{"timesort",			tmsort},
	{"addrandoms",			adrand},
	{"resampledecaytime",	resampt},
	{"bin",					phgbin}		/* Remember to update NUM_FUNCS when adding function */
};


#ifdef GUIOS
	#ifndef COCOA
		/* Definitions and variables for the user command string dialog */
		#pragma options align=mac68k
		typedef struct {
			Handle		h;
			Rect		box;
			char		kind;
		} ItemListEntry;
		
		#define			kNumDlogItems	4
		static struct {
			short			count;
			ItemListEntry	item[kNumDlogItems];
		} gDlogItemList = {
			kNumDlogItems-1,
			{
				{0, {60,382,80,452},	kButtonDialogItem},
				{0, {60,301,80,371},	kButtonDialogItem},
				{0, {20,20,36,92},		kStaticTextDialogItem | kItemDisableBit},
				{0, {21,103,37,449},	kEditTextDialogItem}
			}
		};
		#pragma options align=reset
	#endif
	
	/* Local argv variables for the command string */
	#define			MAX_ARGS		25
	static char			gArgStr[256];
	static char* 		gArgv[MAX_ARGS + 1];
#endif


/* FUNCTIONS */

/**********************
*	changeDir
*
*	Purpose:	Change the working directory (Mac Only).
*
*	Result:	True unless an error occurs.
***********************/
Boolean changeDir(int argc, char *argv[])

{
	Boolean		okay = true;
	
	
	/* Avoid unused parameter compiler warnings */
	if (argc) {};
	if (argv) {};
	
	#ifdef GUIOS
		#ifdef COCOA
		{
			char*		paramFilePtr;		/* Pointer to parameter file name */
			int			dirEnd;				/* End of directory string */
			Boolean		hasQuotes;			/* Whether file name uses quotes */
			int			i;					/* Index through file name */
			char		dirStr[256];		/* Working directory */
			int			err;				/* Result of system call */
			
			/* Set the current directory to that of the parameter file */
			paramFilePtr = argv[argc-1];
			dirEnd = strlen(paramFilePtr) - 1;
			hasQuotes = ((paramFilePtr[dirEnd] == '\'') || (paramFilePtr[dirEnd] == '"'));
			for (i=dirEnd; i>=0; i--) {
				if (paramFilePtr[i] == '/') {
					break;
				}
				else {
					dirEnd--;
				}
			}
			strncpy(dirStr, "cd ", 3);
			strncpy(&(dirStr[3]), paramFilePtr, dirEnd);
			dirEnd += 3;
			if (hasQuotes) {
				dirStr[dirEnd++] = paramFilePtr[strlen(paramFilePtr)-1];
			}
			dirStr[dirEnd] = '\0';
			err = system(dirStr);
			err = system(dirStr);
		}
		#else
			/* Select a working directory */
			okay = (LbFlSetDir() == noErr);
		#endif
	#else
		LbInPrintf("This function only works on graphical computer interfaces");
	#endif
	
	return(okay);
}


#ifdef MACGUIOS
/**********************
*	getArgs
*
*	Purpose:	Return argc and set argv, as obtained from the user (Non-Unix only).
*
*	Result:		Int -- new value of argc.
***********************/
int getArgs(char ***argvPtr)

{
	int						numArgs = 0;					/* Count of arguments (= argc) */
	#ifndef COCOA
		Rect				dlogBounds = {66,44,166,516};	/* Dialog bounds */
		DialogPtr			theDialogPtr;					/* New dialog ptr */
		Str255				dialogText;						/* Text sent to Dialog Manager */
		DialogItemIndex		theItem;						/* Clicked dialog item */
		Str255				cmdStr255;						/* User command string */
		int					cmdLen;							/* Length of command string */
	#endif
	char					cmdStr[256];					/* Obtained command string */
	
	
	cmdStr[0] = '\0';
	#ifdef COCOA
		/* argvPtr already contains the command line text */
		if (argvPtr) {
			strncpy(&(cmdStr[0]), (char*)(*argvPtr), 256);
		};
	#else
		/* Make a new dialog to display */
		{
			Handle		itemsHdl;		/* gDlogItemList as handle */
			OSErr		err;			/* Result of ptr->hdl conversion */
			
			
			/* Convert the dialog item list */
			err = PtrToHand(&gDlogItemList, &itemsHdl, sizeof(gDlogItemList));
			
			/* Properly center and place the dialog's bounding rectangle */
			{
				#ifdef OSX
					CGRect		cgScreenRect;		/* CG device screen rectangle */
				#else
					GDHandle	theMainDevice;		/* Screen device */
					Rect*		screenRectPtr;		/* Screen device rectangle ptr */
				#endif
				Rect			screenRect;			/* Screen device rectangle */
				short			screenHeight;		/* Screen height */
				short			screenWidth;		/* Screen width */
				short			dlogHeight;			/* Dialog height */
				short			dlogWidth;			/* Dialog width */
				
				/* Get the screen rectangle */
				#ifdef OSX
					cgScreenRect = CGDisplayBounds(CGMainDisplayID());
					screenRect.top = cgScreenRect.origin.y + GetMBarHeight();
					screenRect.left = cgScreenRect.origin.x;
					screenRect.bottom = cgScreenRect.origin.y + cgScreenRect.size.height;
					screenRect.right = screenRect.left + cgScreenRect.size.width;
				#else
					theMainDevice = GetMainDevice();
					screenRectPtr = &((*theMainDevice)->gdRect);
					screenRect.top = screenRectPtr->top + GetMBarHeight();
					screenRect.left = screenRectPtr->left;
					screenRect.bottom = screenRectPtr->bottom;
					screenRect.right = screenRectPtr->right;
				#endif
				
				/* Get heights and widths */
				screenHeight = screenRect.bottom - screenRect.top;
				screenWidth = screenRect.right - screenRect.left;
				dlogHeight = dlogBounds.bottom - dlogBounds.top;
				dlogWidth  = dlogBounds.right  - dlogBounds.left;
				
				/* Set newly centered bounds */
				dlogBounds.top = screenRect.top + screenHeight/5;
				dlogBounds.left = screenRect.left + (screenWidth - dlogWidth)/2;
				dlogBounds.bottom = dlogBounds.top + dlogHeight;
				dlogBounds.right = dlogBounds.left + dlogWidth;
			}
			
			/* Make and show the new dialog */
			strncpy((char *)(&dialogText[1]), "SimSET", 6+1);
			dialogText[0] = 6;
			theDialogPtr = NewDialog( NULL, &dlogBounds, dialogText, 1, kMovableModalWindowClass, 
									(WindowPtr)-1, 0, 0, itemsHdl );
			{
				short 			kind;				/* Dialog item kind */
				Handle			item;				/* Dialog item */
				Rect 			box;				/* Dialog item box */
				Boolean 		isDefault = true;	/* Dialog item default */
				
				GetDialogItem(theDialogPtr, kStdOkItemIndex, &kind, &item, &box);
				strncpy((char *)(&dialogText[1]), "OK", 2+1);
				dialogText[0] = 2;
				#ifdef OSX
					SetControlTitleWithCFString((ControlHandle)item, CFSTR("OK"));
				#else
					SetControlTitle((ControlHandle)item, dialogText);
				#endif
				SetControlData((ControlHandle) item, kControlEntireControl,
					kControlPushButtonDefaultTag, sizeof(isDefault), &isDefault);
				GetDialogItem(theDialogPtr, kStdCancelItemIndex, &kind, &item, &box);
				strncpy((char *)(&dialogText[1]), "Quit", 4+1);
				dialogText[0] = 4;
				#ifdef OSX
					SetControlTitleWithCFString((ControlHandle)item, CFSTR("Quit"));
				#else
					SetControlTitle((ControlHandle)item, dialogText);
				#endif
				GetDialogItem(theDialogPtr, 3, &kind, &item, &box);
				strncpy((char *)(&dialogText[1]), "Command:", 8+1);
				dialogText[0] = 8;
				SetDialogItemText(item, dialogText);
			}
		}
		#ifdef OSX
		{
			/* Attempt to reset the cursor */
			CGDisplayHideCursor(kCGDirectMainDisplay);
			CGDisplayShowCursor(kCGDirectMainDisplay);
		}
		#else
		{
			Cursor			theArrow;				/* Cursor type */
			SetCursor(GetQDGlobalsArrow(&theArrow));
		}
		#endif
		ShowWindow(GetDialogWindow(theDialogPtr));
		
		/* Get the user's command string */
		while (true) {
			/* Run the dialog until released */
			ModalDialog(NULL, &theItem);
			switch(theItem) {
				case kStdOkItemIndex:
				{
					/* User clicked OK */
					
					/* Retrieve the command string */
					{
						short 			kind;				/* Dialog item kind */
						Handle			item;				/* Dialog item */
						Rect 			box;				/* Dialog item box */
						
						GetDialogItem(theDialogPtr, 4, &kind, &item, &box);
						GetDialogItemText(item, (StringPtr)cmdStr255);
					}
					
					/* Convert cmdStr255 from Pascal string to C string */
					cmdLen = cmdStr255[0];
					memmove(cmdStr, cmdStr255+1, cmdLen);
					cmdStr[cmdLen] = '\0';
					if (cmdLen < 255) {
						/* Add extra protection */
						cmdStr[cmdLen+1] = '\0';
					}
					
					DisposeDialog(theDialogPtr);
					
					break;
				}
				
				
				case kStdCancelItemIndex:
				{
					/* User cancelled */
					
					DisposeDialog(theDialogPtr);
					
					exit(0);
					break;
				}
			}
			
			if (theItem == kStdOkItemIndex) {
				/* Dialog is done */
				break;
			}
		}
	#endif
	
	
	/* Parse cmdStr into gArgStr and gArgv */
	{
		LbUsFourByte	currentCharIndex;	/* Index into cmdStr */
		LbUsFourByte	curArgPos;			/* Index into gArgStr */
		Boolean			okay;				/* Got next token */
		LbUsFourByte	tokenLen;			/* Next token length */
		Boolean			isLastToken;		/* Last token ends string */
		
		numArgs = 0;
		currentCharIndex = 0;
		curArgPos = 0;
		do {
			okay = LbPfGetNextToken( (char*)cmdStr, &(gArgStr[curArgPos]), 
										&currentCharIndex, &tokenLen, 
										&isLastToken );
			if (okay && (gArgStr[curArgPos] != '\0')) {
				gArgv[numArgs++] = &(gArgStr[curArgPos]);
				curArgPos += tokenLen + 1;
			}
		} while (okay && (! isLastToken));
		*argvPtr = gArgv;
	}
	
	return (numArgs);
}
#endif


#ifndef COCOA
/**********************
*	main
*
*	Purpose:	Execute the program.
*
*	Result:	Always returns zero indicating no error.
***********************/

int main(int argc, char *argv[])

{
	int			result;					/* Result of running SimSET */
	
	
	/* Run SimSET */
	result = runSimSET(argc, argv);
	
	
	/* Quit the program */
	return (0);
}
#endif


/**********************
*	runSimSET
*
*	Purpose:	Execute the SimSET program.
*
*	Result:	Always returns zero indicating no error.
***********************/
int runSimSET(int argc, char *argv[])
{
	Boolean		okay = false;				/* Process Loop */
	#ifdef GUIOS
	Boolean		firstTime = true;
	#endif
	char		progName[64];
	LbTwoByte	index1;

	/* Turn off buffering to stdio */
	setbuf(stdout, 0);
	
	#ifdef GUIOS
		DO_IT_AGAIN:;
		
		/* For non-Unix systems we need to get the command line string */
		argc = getArgs(&argv);
		
		/* Check for empty command string */
		if (argc <= 0) {
			LbInPrintf("\nYou must enter at least the name of a SimSET program.\n");
			goto FAIL_WD_INIT;
		}
		
		/* Select a working directory */
		if (firstTime) {
			if (changeDir(argc, argv) == false) {
				goto FAIL_WD_INIT;
			}
			firstTime = false;
		}
	#endif
	

	/* We setup arguments before initializing managers in this special case
		because it is a silly text manipulation task that only requires stdlib.
		Doing it first helps in generating better error messages without a lot
		of extra stuff.
	*/
	setup_args(&argc, argv, progName);

	/* Try to init the debug manager */
	if (!LbDbInit("simset (main)", 0, 0)){
		LbInPrintf("\nUnable to initialize debug manager.\n");
		goto FAIL_LB_INIT;
	}

	/* Try to init the error manager */
	if (!ErInit(0, 0)) {
		LbInPrintf("\nUnable to do error library initialization!");
		goto FAIL_LB_INIT;
	}
	
	/* Try to init the memory manager, use LBMMFg_Accounting for debugging info
		otherwise pass zero
	 */
	if (!LbMmInit(0)) {
		LbInPrintf("\nUnable to do memory library initialization!");
		goto FAIL_LB_INIT;
	}

	do {
	

		/* Loop through the function table to find which one the user requested */
		for (index1 = 0; index1 < NUM_FUNCS; index1++) {
			if (strcmp(progName, SimsetFuncTbl[index1].funcName) == 0) {
				okay = SimsetFuncTbl[index1].funcPtr(argc, argv);
				break;
			}
		}
		
		/* See if we did not get a match */
		if (index1 == NUM_FUNCS) {
		
			/* If the first 3 chars are "phg" we'll assume they want
				to run the PHG and give it a shot
			*/
			if (strncmp(progName, "phg", 3) == 0) {
				okay = SimsetFuncTbl[PHG_INDEX].funcPtr(argc, argv);
			}
	 		else if ((progName != 0) && (strlen(progName) != 0)) {
	 		
	 			/* If we are in a GUI system then we may just want to quit */
	 			#ifdef GUIOS
	 				if (progName[0] == 'q') {
	 					okay = true;
	 					break;
	 				}
	 			#endif
	 			
 				LbInPrintf("\nUtility name '%s' is not supported.\n", argv[0]);
 				LbInPrintf("Supported utilities are:\n");
 				for (index1 = 0; index1 < NUM_FUNCS; index1++){
 					LbInPrintf("\t%s\n", SimsetFuncTbl[index1].funcName);
 				}
 				LbInPrintf("\n");
	 			sprintf(errString, "\nUnsupported utility: '%s'.\n", argv[0]);
	 			ErStGeneric(errString);
	 		}
	 		else {
	 			ErStGeneric("\nThis program requires a command line argument.");
	 		}
 		}
		
	} while (false);
	
	/* Handle error situation if one exists */
	if (!okay) {
		sprintf(errString, "An error occured while executing '%s'.",
			progName);
		ErHandle(errString, false);
	}
	
	/* Terminate the error library */
	ErTerminate();
	
	/* Terminate the memory library */
	LbMmTerminate();

	#ifdef GUIOS
		#ifndef COCOA
			if (progName[0] != 'q')
				goto DO_IT_AGAIN;
		#endif
	#endif
	
	/* If we can't initialize the library, we'll jump to here and exit */
	FAIL_LB_INIT:;
	FAIL_WD_INIT:;
	

	/* Quit SimSET */
	return (0);
}


/**********************
*	setup_args
*
*	Purpose:	Sets up the main arguments.
*
*	Result:	void.
***********************/
void setup_args(int *argcPtr, char *argv[], char *progName)
{
	LbTwoByte	index1, index2;
	LbUsTwoByte	length;

	
	#ifdef WINNT
		LbTwoByte argIndex;
		/*
			We don't want the name of the program from WINNT, we
			want the line they typed in. Since in WINNT the first
			argument is always the command to start the program
		*/
		for (argIndex = 0; argIndex < (*argcPtr)-1; argIndex++){		/* was argc - not argcPtr - before Hopkins changes */
			argv[argIndex] = argv[argIndex+1];
		}
		*argcPtr--;		/* was argc - not argcPtr - before Hopkins changes */

		/* 	We need to peel off any part of a WINNT (or DOS) path
			that might be attached to the program name.
		*/
		length = strlen(argv[0]);
		for (index1 = length-1; index1 >= 0; index1--) {
			if (argv[0][index1] == '\\')
				break;
		}
	#else		
		int*		dummyPtr;			/* Removes compiler warning */
		
		/* Avoid unused parameter compiler warnings */
		dummyPtr = argcPtr;
		
		/* We need to peel off any part of a UNIX path that might
			be attached to the program name.
		*/
		length = strlen(argv[0]);
		for (index1 = length-1; index1 >= 0; index1--) {
			if (argv[0][index1] == '/')
				break;
		}
	#endif
	
	/* Our index is either -1, or on the / */
	index1++;
	
	/* Compute length of just the name */
	length = length - index1;
	
	/* Copy just the program name (not the path) into progName */
	for (index2 = 0; index2 < length; index2++){
		progName[index2] = argv[0][index1];
		index1++;
	}
	progName[index2] = '\0';
	
			
	/* Now convert program name to lower case for matching */
	length = strlen(progName);
	for (index1 = 0; index1 < length; index1++) {
		progName[index1] = tolower((int)progName[index1]);
	}
}
