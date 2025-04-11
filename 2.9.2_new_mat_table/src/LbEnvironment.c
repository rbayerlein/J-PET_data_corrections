/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992-2013 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name:        LbEnvironment.c
*     Revision Number:    1.6
*     Date last revised:  24 July 2013
*     Programmer:         Steven Vannoy
*     Date Originated:    Tuesday, July 21, 1992
*
*     Module Overview:	This module provides routines for getting a program's
*						environment variables.
*						MACHINE DEPENDENT						
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
**********************************************************************************
*
*	  Revision Section (Also update version number, if relevant)
*
*	  Programmer(s):		
*
*	  Revision date:		
*
*	  Revision description:
*
**********************************************************************************
*
*	  Revision Section (Also update version number, if relevant)
*
*	  Programmer(s):		Steven Gillispie
*
*	  Revision date:		24 July 2013
*
*	  Revision description:	Combined all into one definition, except for specials
*
**********************************************************************************
*
*	  Revision Section (Also update version number, if relevant)
*
*	  Programmer(s):		Steven Gillispie
*
*	  Revision date:		18 September 2012
*
*	  Revision description:	Changed to use -1 instead of EOF from getopt 
*								according to 1992 POSIX standards change
*
*********************************************************************************/

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"

#include "LbMacros.h"
#include "LbError.h"
#include "LbEnvironment.h"
#include "LbInterface.h"

/* PROTOTYPES */
#ifdef GUIOS
LbFourByte lbEnParseOptions(int argc, char **argv, char *optStr);
#endif
#ifdef	WINNT
LbFourByte lbEnParseOptions(int argc, char **argv, char *optStr);
#endif

/* LOCAL GLOBALS */
#ifdef GUIOS
char			*lbEnOptArg = 0;
LbFourByte		lbEnCurOpt = -1;
LbUsFourByte	lbEnOptInd = 0;
#endif
#ifdef WINNT
char			*lbEnOptArg = 0;
LbFourByte		lbEnCurOpt = -1;
LbUsFourByte	lbEnOptInd = 0;
#endif


/* DEF for all except later special exceptions (AOS_VS, GUIOS, WINNT) */
#ifndef AOS_VS
#ifndef GUIOS
#ifndef WINNT
/*****************************************
*		LbEnGetOptions
*
*	Arguments:
*				int				argc			- Number of arguments.
*				char			**argv			- Command line.
*				char			**optStr		- Valid options.
*				LbUsFourByte	*optFlags		- Flag storage for set options.
*				char			**optArgs		- Option arguments.
*				LbUsFourByte	optArgFlags		- Option argument flags.
*				LbUsFourByte	*firstArg		- Index into argv for first argument.
*
*	Purpose:	Determine which options where set on the
*				command line.
*
*	Returns:	TRUE if no error occurs.
*
******************************************/
Boolean	LbEnGetOptions(int argc, char **argv, char **optStr,
				LbUsFourByte *optFlags, char optArgs[][LBEnMxArgLen],
				LbUsFourByte optArgFlags,
				LbUsFourByte *firstArg)
{
	extern	char *optarg;
	extern	int optind, opterr;
	
	Boolean 		okay = false;
	LbUsFourByte	index, colonCount = 0;
	LbFourByte		curOpt, optBitNum = -1;
	
	
	if (optArgFlags) {};		/* Eliminate unused parameter compiler warning */
	
	do	/* Process Loop */
	{
		/* Disable error message */
		opterr = 0;
		
		/* Loop through arguments, setting flags in appropriate order */
		while ((curOpt = getopt(argc, argv, *optStr)) != -1)
		{
			/* See if the option is in our string */
			index = 0;
			colonCount = 0;
			while ((*optStr)[index])
			{
				/* See if we have a colon */
				if ((*optStr)[index] == ':')
					colonCount++;
					
				/* Check the string */
				if ((*optStr)[index] == curOpt)
				{
					/* Save the index */
					optBitNum = index - colonCount;

					/* See if there is an argument specifier */
					if ((*optStr)[index+1] == ':')
					{
						/* Copy in the argument if they supplied space */
						#ifdef IF_LBDEBUG
							if (!optArgs)
							{
								LbInPrintf("You must supply storage for arguments");
								break;
							}
						#endif

						/* Save the switch value */
						strcpy(optArgs[optBitNum], optarg);
					}
					break;
				}
				index++;
			}
			/* It is an error if we didn't find the option */
			if (curOpt == '?')
			{
				ErStGeneric("Unknown option specified!");
				break;
			}
			
			/* Set our bit */
			if (optBitNum == 0)
				LbFgSet(*optFlags, LBFlag0);
			else
				LbFgSet(*optFlags, (2 << (optBitNum - 1)));
		}
		/* If we broke before -1, it is an error */
		if (curOpt != -1)
			break;
		
		/* Store the index for argv into firstArg */
		*firstArg = optind;
		
		okay = true;
	} while (false);
	return (okay);
}
#endif
#endif
#endif  /* Normal options, not special exceptions (as below) */

#ifdef	AOS_VS
/*****************************************
*		LbEnGetOptions
*
*	Arguments:
*				int				argc			- Number of arguments.
*				char			**argv			- Command line.
*				char			**optStr		- Valid options.
*				UsFourByte		*optFlags		- Flag storage for set options.
*				char			**optArgs		- Option arguments.
*				LbUsFourByte	optArgFlags		- Which options take arguments.
*				LbUsTwoByte		*firstArg		- Index into argv for first argument.
*
*	Purpose:	Determine which options where set on the
*				command line.
*
*	Returns:	TRUE if no error occurs.
*
******************************************/
Boolean	LbEnGetOptions(int argc, char **argv, char **optStr,
				LbUsFourByte *optFlags, char optArgs[][LBEnMxArgLen],
				LbUsFourByte optArgFlags, LbUsFourByte *firstArg)
{
	Boolean 	okay = false, badArg = false;
	SWITCH		sw;
	LbFourByte	curOpt;
	
	do	/* Process Loop */
	{
		/* Clear out our switch space */
		zero((char *)&sw, sizeof(sw));
		
		/* Initialize the switch */
		sw.sw_start = argv[0];
		sw.sw_names = optStr;
		
		/* Process each value */
		while((curOpt = getswitch(&sw)) != SW_EOF)
		{
			/* See if invalid switch */
			if ((curOpt == SW_NOT_FOUND) || (curOpt == SW_NOT_UNIQUE))
			{
				badArg = true;
				break;
			}
			
			/* Set flag for this option */
			{
				/* See if this is the first option */
				if (curOpt == 0)
					LbFgSet(*optFlags, LBFlag0);
				else
					LbFgSet(*optFlags, (2 << (curOpt -1)));
			}
			
			/* Get switch value if it exists */
			{
				#ifdef LB_DEBUG
					if (sw.sw_arg && !optArgs)
					{
						LbInPrintf("You must supply storage for arguments.");
						break;
					}
				#endif
				/* If switch has value, copy it */
				if (sw.sw_arg)
				{
					/* Make sure the caller wants to take an argument */
					if (curOpt == 0)
						badArg = !(LbFgIsSet(optArgFlags, LBFlag0));
					else
						badArg = !(LbFgIsSet(optArgFlags, (2 << (curOpt -1))));
						
					if (badArg)
						break;
						
					/* Save the switch value */
					strncpy(optArgs[curOpt], sw.sw_arg, sw.sw_arglen);
				}
			}
		}
		/* Store the index for argv into firstArg */
		*firstArg = 1;

		/* If we got a bad argument bolt */
		if (badArg)
			break;
			
		okay = true;
	} while (false);
	
	return (okay);
}
#endif /* AOS_VS */

#ifdef	GUIOS
/*****************************************
*		LbEnGetOptions
*
*	Arguments:
*				int				argc			- Number of arguments.
*				char			**argv			- Command line.
*				char			**optStr		- Valid options.
*				UsFourByte		*optFlags		- Flag storage for set options.
*				char			**optArgs		- Option arguments.
*				LbUsFourByte	optArgFlags		- Which options take arguments.
*				LbUsTwoByte		*firstArg		- Index into argv for first argument.
*
*	Purpose:	Determine which options where set on the
*				command line.
*
*	Returns:	TRUE if no error occurs.
*
******************************************/
Boolean	LbEnGetOptions(int argc, char **argv, char **optStr,
				LbUsFourByte *optFlags, char optArgs[][LBEnMxArgLen],
				LbUsFourByte optArgFlags, LbUsFourByte *firstArg)
{
	Boolean 		okay = false;
	LbUsFourByte	index, colonCount = 0;
	LbFourByte		curOpt, optBitNum = -1;
	
	if (optArgFlags) {};		/* Eliminate unused parameter compiler warning */
	
	do	/* Process Loop */
	{
		/* Loop through arguments, setting flags in appropriate order */
		while ((curOpt = lbEnParseOptions(argc, argv, *optStr)) != -1)
		{
			/* See if the option is in our string */
			index = 0;
			colonCount = 0;
			while ((*optStr)[index])
			{
				/* See if we have a colon */
				if ((*optStr)[index] == ':')
					colonCount++;
					
				/* Check the string */
				if ((*optStr)[index] == curOpt)
				{
					/* Save the index */
					optBitNum = index - colonCount;

					/* See if there is an argument specifier */
					if ((*optStr)[index+1] == ':')
					{
						/* Copy in the argument if they supplied space */
						#ifdef IF_LBDEBUG
							if (!optArgs) {
								ErAbort("You must supply storage for arguments");
							}
						#endif

						/* If argument not supplied, its an error */
						if (lbEnOptArg == 0) {
							ErStGeneric("Program option requires argument.");
							goto FAIL;
						}
						/* Save the switch value */
						strcpy(optArgs[optBitNum], lbEnOptArg);
					}
					break;
				}
				index++;
			}
			/* It is an error if we didn't find the option */
			if (curOpt == '?')
			{
				ErStGeneric("Unknown option specified!");
				break;
			}
			
			/* Set our bit */
			if (optBitNum == 0)
				LbFgSet(*optFlags, LBFlag0);
			else
				LbFgSet(*optFlags, (2 << (optBitNum - 1)));
		}
		/* If we broke before -1, it is an error */
		if (curOpt != -1)
			break;
		
		/* Store the index for argv into firstArg */
		*firstArg = lbEnOptInd;
		
		okay = true;
		
		FAIL:;
	} while (false);
	return (okay);
}
/*****************************************
*		lbEnParseOptions
*
*	Arguments:
*				int				argc			- Number of arguments.
*				char			**argv			- Command line.
*				char			*optStr			- Valid options.
*
*	Purpose:	Determine which options where set on the
*				command line.
*
*				NOTE: This function was written one afternoon to provide
*				functionality under MPW.  It was modeled after the GEN_UNIX
*				function getopt and hence is weird, and doesn't do error
*				handling well.  It would be good to rewrite this function
*				and use it for both the MPW and GEN_UNIX versions.
*
*	Returns:	Next option letter in argv that matches a letter in optstring.
*
******************************************/
LbFourByte lbEnParseOptions(int argc, char **argv, char *optStr)
{
	char		*curArgString;
	char		curArg;
	LbTwoByte 	curIndex = 0;
	
	/* if first time called, Set option index to first argument after program name */
	if (lbEnCurOpt == -1)
		lbEnOptInd = 1;
	
	do	{ /* Process Loop */
	
		/* If we have looked at all arguments indicate -1 and bolt */
		if (lbEnOptInd == (LbUsFourByte)argc) {
		
			lbEnCurOpt = -1;
			break;
		}

		/* If this program argument is not an option variable, bolt */
		if (*argv[lbEnOptInd] != '-') {
		
			lbEnCurOpt = -1;
			break;
		}
		
		/* Index to current argument to check */
		curArgString = argv[lbEnOptInd];
		curArg = curArgString[1];
		curIndex = 0;
		
		/* Loop through possible arguments to find match */
		while (optStr[curIndex]) {
			
			/* if string matches, bolt */
			if (optStr[curIndex] == curArg) {
				lbEnCurOpt = curArg;
				break;
			}
				
			/* Check next char */
			curIndex++;
			if (optStr[curIndex] == ':')
				curIndex++;
		}
		
		/* See if we found the thing */
		if (optStr[curIndex] != 0) {
		
			lbEnOptInd++;

			/* Now get arguments if there are any */
			if (optStr[curIndex+1] == ':')
			{
				/* If no more arguments, but expecting one, bolt */
				if (argv[lbEnOptInd] == 0) {
					lbEnOptArg = 0;
					break;
				}
				
				/* Index to next argument on the line */
				curIndex = 0;
				curArgString = argv[lbEnOptInd];
				
				/* Make sure next thing is not another argument */
				if (curArgString[0] == '-') {
					lbEnCurOpt = '?';
					break;
				}
				else {
					lbEnOptArg = argv[lbEnOptInd];
					lbEnOptInd++;
				}
					
			}
		}
		else {
			/* Option not found, bolt */
			lbEnCurOpt = '?';
			break;
		}
	
	} while (false);
	
	return (lbEnCurOpt);
}
#endif /* GUIOS */

#ifdef	WINNT
/*****************************************
*		LbEnGetOptions
*
*	Arguments:
*				int				argc			- Number of arguments.
*				char			**argv			- Command line.
*				char			**optStr		- Valid options.
*				UsFourByte		*optFlags		- Flag storage for set options.
*				char			**optArgs		- Option arguments.
*				LbUsFourByte	optArgFlags		- Which options take arguments.
*				LbUsTwoByte		*firstArg		- Index into argv for first argument.
*
*	Purpose:	Determine which options where set on the
*				command line.
*
*	Returns:	TRUE if no error occurs.
*
******************************************/
Boolean	LbEnGetOptions(int argc, char **argv, char **optStr,
				LbUsFourByte *optFlags, char optArgs[][LBEnMxArgLen],
				LbUsFourByte optArgFlags, LbUsFourByte *firstArg)
{
	Boolean 		okay = false;
	LbUsFourByte	index, colonCount = 0;
	LbFourByte		curOpt, optBitNum = -1;
	
	if (optArgFlags) {};		/* Eliminate unused parameter compiler warning */
	
	do	/* Process Loop */
	{
		/* Loop through arguments, setting flags in appropriate order */
		while ((curOpt = lbEnParseOptions(argc, argv, *optStr)) != -1)
		{
			/* See if the option is in our string */
			index = 0;
			colonCount = 0;
			while ((*optStr)[index])
			{
				/* See if we have a colon */
				if ((*optStr)[index] == ':')
					colonCount++;
					
				/* Check the string */
				if ((*optStr)[index] == curOpt)
				{
					/* Save the index */
					optBitNum = index - colonCount;

					/* See if there is an argument specifier */
					if ((*optStr)[index+1] == ':')
					{
						/* Copy in the argument if they supplied space */
						#ifdef IF_LBDEBUG
							if (!optArgs) {
								ErAbort("You must supply storage for arguments");
							}
						#endif

						/* If argument not supplied, its an error */
						if (lbEnOptArg == 0) {
							ErStGeneric("Program option requires argument.");
							goto FAIL;
						}
						/* Save the switch value */
						strcpy(optArgs[optBitNum], lbEnOptArg);
					}
					break;
				}
				index++;
			}
			/* It is an error if we didn't find the option */
			if (curOpt == '?')
			{
				ErStGeneric("Unknown option specified!");
				break;
			}
			
			/* Set our bit */
			if (optBitNum == 0)
				LbFgSet(*optFlags, LBFlag0);
			else
				LbFgSet(*optFlags, (2 << (optBitNum - 1)));
		}
		/* If we broke before -1, it is an error */
		if (curOpt != -1)
			break;
		
		/* Store the index for argv into firstArg */
		*firstArg = lbEnOptInd;
		
		okay = true;
		
		FAIL:;
	} while (false);
	return (okay);
}
/*****************************************
*		lbEnParseOptions
*
*	Arguments:
*				int				argc			- Number of arguments.
*				char			**argv			- Command line.
*				char			*optStr			- Valid options.
*
*	Purpose:	Determine which options where set on the
*				command line.
*
*				NOTE: This function was written one afternoon to provide
*				functionality under MPW.  It was modeled after the GEN_UNIX
*				function getopt and hence is weird, and doesn't do error
*				handling well.  It would be good to rewrite this function
*				and use it for both the MPW and GEN_UNIX versions.
*
*	Returns:	Next option letter in argv that matches a letter in optstring.
*
******************************************/
LbFourByte lbEnParseOptions(int argc, char **argv, char *optStr)
{
	char		*curArgString;
	char		curArg;
	LbTwoByte 	curIndex = 0;
	
	/* if first time called, Set option index to first argument after program name */
	if (lbEnCurOpt == -1)
		lbEnOptInd = 1;
	
	do	{ /* Process Loop */
	
		/* If we have looked at all arguments indicate -1 and bolt */
		if (lbEnOptInd == argc) {
		
			lbEnCurOpt = -1;
			break;
		}

		/* If this program argument is not an option variable, bolt */
		if (*argv[lbEnOptInd] != '-') {
		
			lbEnCurOpt = -1;
			break;
		}
		
		/* Index to current argument to check */
		curArgString = argv[lbEnOptInd];
		curArg = curArgString[1];
		curIndex = 0;
		
		/* Loop through possible arguments to find match */
		while (optStr[curIndex]) {
			
			/* if string matches, bolt */
			if (optStr[curIndex] == curArg) {
				lbEnCurOpt = curArg;
				break;
			}
				
			/* Check next char */
			curIndex++;
			if (optStr[curIndex] == ':')
				curIndex++;
		}
		
		/* See if we found the thing */
		if (optStr[curIndex] != 0) {
		
			lbEnOptInd++;

			/* Now get arguments if there are any */
			if (optStr[curIndex+1] == ':')
			{
				/* If no more arguments, but expecting one, bolt */
				if (argv[lbEnOptInd] == 0) {
					lbEnOptArg = 0;
					break;
				}
				
				/* Index to next argument on the line */
				curIndex = 0;
				curArgString = argv[lbEnOptInd];
				
				/* Make sure next thing is not another argument */
				if (curArgString[0] == '-') {
					lbEnCurOpt = '?';
					break;
				}
				else {
					lbEnOptArg = argv[lbEnOptInd];
					lbEnOptInd++;
				}
					
			}
		}
		else {
			/* Option not found, bolt */
			lbEnCurOpt = '?';
			break;
		}
	
	} while (false);
	
	return (lbEnCurOpt);
}
#endif /* WINNT */
