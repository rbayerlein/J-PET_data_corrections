/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 2003-2013 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
********************************************************************************/

/********************************************************************************
*
*     Module Name:        LbFile.c
*     Revision Number:    1.16
*     Date last revised:  4 June 2013
*     Programmer:         Steven Gillispie
*     Date Originated:    May 7, 2003
*
*     Module Overview:    This module provides routines for performing file I/O
*							tasks.
*						  MACHINE DEPENDENT.
*
*     References:         
*
*********************************************************************************
*
*     Global functions defined:
*		LbFlSetDir (Mac only)
*		LbFlFileOpen
*		LbFlFGetS
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
*			Revision date:		17 February 2012
*
*			Revision description:	Changed LbFlFileOpen to use @simset file name 
*										symbol
*
********************************************************************************/

#include	<stdio.h>
#include	<string.h>
#include	<ctype.h>
#include	<errno.h>

#include "SystemDependent.h"

#include "LbFile.h"

#include "LbTypes.h"
#include "LbError.h"

#ifdef __MWERKS__
	#include	<Gestalt.h>
	#include	<Navigation.h>
#endif


/*	CONSTANTS */

/*	LOCAL TYPES */

/*  LOCAL GLOBALS */
#define			kMaxPathLength  512				/* Maximum pathname length */
static char		kSimSETTop[] = {"@simset/"};	/* File name symbol for kSimSETPath/ */
#ifdef OSX
static char		gDefaultDir[kMaxPathLength];	/* Pathname of default directory */
#endif

/*	LOCAL MACROS */

/*	LOCAL FUNCTIONS */
#ifdef MACGUIOS
	char* lbFlUnix2MacPath( char *pathnameUnix,  char *pathnameMac );
#endif


/*	FUNCTIONS */
#ifdef MACGUIOS
#ifndef COCOA
/*********************************************************************************
*
*	Name:			LbFlSetDir
*
*	Summary:		Set the default directory.
*					MAC ONLY.
*					A Carbon-compatible OS is required.
*
*	Arguments:		None.
*
*	Function return:		LbTwoByte error code.
*
*********************************************************************************/

LbTwoByte LbFlSetDir( void )

{
	LbTwoByte			err = noErr;	/* Result of function */
	#ifdef OSX
		NavDialogCreationOptions
						theOptions;		/* Nav Svcs dialog options */
		NavDialogRef	theDLog;		/* Created dialog for Nav Svcs ChooseFolder */
		NavUserAction	userResponse;	/* User response from Nav Svcs ChooseFolder */
		NavReplyRecord	reply;			/* User selection info from Nav Svcs ChooseFolder */
	#else
		long			carbonVersion;	/* BCD version of Carbon */
		NavReplyRecord	reply;			/* User response from NavChooseFolder */
	#endif
	
	
	#ifdef OSX
		err = NavGetDefaultDialogCreationOptions( &theOptions );
		if ( err == noErr ) {
			theOptions.optionFlags = kNavDefaultNavDlogOptions & ~kNavAllowMultipleFiles;
			theOptions.windowTitle = CFSTR( "Default directory selection" );
			err = NavCreateChooseFolderDialog( &theOptions, NULL, NULL, NULL, &theDLog );
			if ( err == noErr ) {
				err = NavDialogRun( theDLog );
				
				if ( err == noErr ) {
					userResponse = NavDialogGetUserAction( theDLog );
					if ( (userResponse != kNavUserActionNone) && 
							(userResponse != kNavUserActionCancel) ) {
						err = NavDialogGetReply( theDLog, &reply );
						
						if ( err == noErr ) {
							if ( reply.validRecord ) {
								AEKeyword		theKey;			/* Keyword of folder spec */
								DescType		theType;		/* Type of folder spec */
								FSRef			theFolderRef;	/* User-selected folder ref */
								Size			theSize;		/* Size of folder spec */
								
								/* Get the returned folder */
								err = AEGetNthPtr( &reply.selection, 
												1, typeFSRef, &theKey, &theType, 
												(Ptr)&theFolderRef, sizeof(FSRef), &theSize );
								if ( err == noErr ) {
									/* Set the default folder */
									err = FSRefMakePath( &theFolderRef, 
															(UInt8*)gDefaultDir, 
															kMaxPathLength-1 );
								}
							}
							NavDisposeReply( &reply );
						}
					}
				}
				
				NavDialogDispose( theDLog );
			}
		}
	#else
		if ( ( Gestalt( gestaltCarbonVersion, &carbonVersion ) == noErr )  &&  
				( carbonVersion >= 0x100 ) ) {
			/* Have at least Carbon 1.0, which includes Navigation Services */
			err = NavChooseFolder( NULL, &reply, NULL, NULL, NULL, NULL );
			if ( err == noErr ) {
				if ( reply.validRecord ) {
					AEKeyword		theKey;			/* Keyword of folder spec */
					DescType		theType;		/* Type of folder spec */
					FSSpec			theFolderSpec;	/* User-selected folder spec */
					Size			theSize;		/* Size of folder spec */
					
					/* Get the returned folder */
					err = AEGetNthPtr( &reply.selection, 
									1, typeFSS, &theKey, &theType, 
									(Ptr)&theFolderSpec, sizeof(FSSpec), &theSize );
					if ( err == noErr ) {
						/* Set the default folder */
						err = HSetVol( NULL, 
										theFolderSpec.vRefNum, 
										theFolderSpec.parID );
					}
				}
				NavDisposeReply( &reply );
			}
			else {
				/* Default folder will remain as before; original is app folder */
			}
		}
		else {
			/* Can't run on this version of the Mac OS */
			ErAbort( "This version of the Mac OS is too old for SimSET." );
		}
	#endif
	
	return( err );
}


/*********************************************************************************
*
*	Name:			lbFlUnix2MacPath
*
*	Summary:		Convert pathnameUnix in Unix format to pathnameMac in Mac format.
*					MAC ONLY.
*
*	Arguments:		
*		char*				pathnameUnix	- Unix file path/name.
*		char*				pathnameMac	- Mac file path/name.
*
*	Function return:		char* to pathnameMac.
*
*********************************************************************************/

char* lbFlUnix2MacPath( char *pathnameUnix,  char *pathnameMac )
{
	char		localPath[kMaxPathLength];	/* Local pathname */
	char		*unixPtr;					/* Pointer into pathnameUnix */
	char		*macPtr;					/* Pointer into pathnameMac */
	char		*endPtr;					/* Last Unix path position */
	size_t		nameLength;					/* Length of dir/file name */
	
	
	/* Make local copy */
	strncpy( localPath, pathnameUnix, kMaxPathLength );
	
	/* Check for simple pathname */
	if ( ! strchr( localPath, '/' ) ) {
		/* Simple file name with no directories */
		strcpy( pathnameMac, localPath );
	}
	else {
		unixPtr = localPath;
		macPtr = pathnameMac;
		
		if ( localPath[0] == '/' ) {
			/* Full, absolute pathname */
			/* Mac pathname has no initial character */
			/* NOTE:  Mac cannot have a file at the topmost level */
			unixPtr++;
		}
		else {
			/* Add current directory label to rule out absolute path */
			*macPtr++ = ':';
		}
		
		while ( unixPtr ) {
			/* Check for . or .. */
			if ( unixPtr[0] == '.' ) {
				if ( unixPtr[1] == '/' ) {
					/* Current directory; can be ignored */
					unixPtr++;
					unixPtr++;
				}
				else if ( ( unixPtr[1] == '.' ) && ( unixPtr[2] == '/' ) ) {
					/* Upper directory */
					*macPtr++ = ':';
					unixPtr++;
					unixPtr++;
					unixPtr++;
				}
				else {
					/* Dot file name */
					/* Handle normally */
					endPtr = strchr( unixPtr, '/' );
					if ( endPtr ) {
						/* Copy directory name */
						nameLength = endPtr-unixPtr;
						strncpy( macPtr, unixPtr, nameLength );
						macPtr += nameLength;
						*macPtr++ = ':';
						unixPtr += nameLength + 1;
					}
					else {
						/* Copy file name */
						strcpy( macPtr, unixPtr );
						unixPtr = endPtr;
					}
				}
			}
			else {
				/* Handle normally */
				endPtr = strchr( unixPtr, '/' );
				if ( endPtr ) {
					/* Copy directory name */
					nameLength = endPtr-unixPtr;
					strncpy( macPtr, unixPtr, nameLength );
					macPtr += nameLength;
					*macPtr++ = ':';
					unixPtr += nameLength + 1;
				}
				else {
					/* Copy file name */
					strcpy( macPtr, unixPtr );
					unixPtr = endPtr;
				}
			}
		}
	}
	
	/* Courtesy return */
	return( pathnameMac );
}
#endif
#endif


/*********************************************************************************
*
*	Name:			LbFlFileOpen
*
*	Summary:		Open a file referred to by *path according to *mode.
*					MACHINE DEPENDENT.
*
*	Arguments:		
*		char*				path	- File path/name.
*		char*				mode	- read/write options.
*
*	Function return:		FILE* to opened file.
*
*********************************************************************************/

FILE*	LbFlFileOpen( char *path,  char *mode )

{
	char		pathStart[kMaxPathLength];	/* Holder of start of file name path */
	int			length;						/* Length of kSimSETTop */
	int			i;							/* Index variable */
	char		usedPath[kMaxPathLength];	/* File name path actually used */
	char		localPath[kMaxPathLength];	/* Local pathname */
	FILE*		theFile;					/* Returned value */
	int			err;						/* errno */
	
	
	/* Check for, and if appropriate substitute for, kSimSETTop */
	strncpy( pathStart, path, kMaxPathLength-1 );
	pathStart[kMaxPathLength-1] = '\0';
	length = strlen(kSimSETTop);
	for (i=1; i<length; i++) {
		/* Make case-insensitive (all lower case) */
		pathStart[i] = tolower( pathStart[i] );
	}
	if ( strncmp( pathStart, kSimSETTop, strlen(kSimSETTop) ) == 0 ) {
		/* The file name starts with kSimSETTop */
		strncpy( usedPath, kSimSET_Path, kMaxPathLength-1 );
		usedPath[kMaxPathLength-1] = '\0';
		strncat( usedPath, &(path[strlen(kSimSETTop)-1]), kMaxPathLength-strlen(usedPath) );
	}
	else {
		/* The file name does not start with kSimSETTop */
		strncpy( usedPath, path, kMaxPathLength );
	}
	
	/* Remove any quotes */
	if ( (usedPath[0] == '\'') || (usedPath[0] == '"') ) {
		strncpy( usedPath, &(usedPath[1]), kMaxPathLength );
		usedPath[strlen(usedPath)-1] = '\0';
	}
	
	/* Perform a system-appropriate file open */
	#ifndef MACGUIOS
	{
		strncpy( localPath, usedPath, kMaxPathLength );
	}
	#else
	{
		/* Mac OS file open */
		#ifdef OSX
		{
			if ( usedPath[0] == '/' ) {
				/* Path starts with top level, so is good as is */
				strncpy( localPath, usedPath, kMaxPathLength );
			}
			else {
				/* Combine the set default directory with the supplied path */
				strncpy( localPath, gDefaultDir, kMaxPathLength );
				strncat( localPath, "/", kMaxPathLength-strlen(localPath) );
				strncat( localPath, usedPath, kMaxPathLength-strlen(localPath) );
			}
		}
		#else
			#ifdef COCOA
				strncpy( localPath, usedPath, kMaxPathLength );
			#else
				/* Convert from Unix pathnames to Mac pathnames */
				lbFlUnix2MacPath( usedPath, localPath );
			#endif
		#endif
	}
	#endif
	
	/* Open the file */
	theFile = fopen( localPath, mode );
	err = errno;
	
	/* Return the file */
	return( theFile );
}


/*********************************************************************************
*
*	Name:			LbFlFGetS
*
*	Summary:		Read str from stream, up to size (equiv to fgets).
*
*	Arguments:		
*		char*				str		- Returned input string.
*		int					size	- Size of string memory.
*		FILE*				stream	- Stream to read str from.
*
*	Function return:		char* to str.
*
*********************************************************************************/

char*	LbFlFGetS( char *str, int size, FILE *stream )
{
	char		*result;		/* Function result */
	char		*p;				/* Pointer into str */
	int			i;				/* Index through str */
	char		c;				/* Read-in character */
	
	
	if ( size <= 0 ) {
		/* Can't accept 0 size */
		result = NULL;
	}
	else {
		result = str;	/* Non-null initialization */
		p = str;
		for ( i=0; i<(size-1); i++ ) {
			/* Read a character from stream */
			c = getc( stream );
			
			if ( c == 10 ) {	/* LF */
				/* End of line */
				*p++ = '\n';
				*p = '\0';
				result = str;
				
				break;
			}
			else if ( c == 13 ) {	/* CR */
				/* End of line */
				*p++ = '\n';
				*p = '\0';
				result = str;
				
				/* Check for trailing LF */
				c = getc( stream );
				if ( c == 10 ) {
					/* Discard it */
				}
				else {
					/* Replace the non-LF */
					ungetc( c, stream );
				}
				
				break;
			}
			else if ( c == EOF ) {	/* EOF */
				/* End of file or error */
				if ( ferror( stream ) ) {
					/* Error */
					result = NULL;
				}
				else {
					/* End of file */
					if ( p == str ) {
						/* No characters ever read */
						result = NULL;
					}
					else {
						/* End of string */
						*p = '\0';
						result = str;
					}
				}
				
				break;
			}
			else {
				/* Normal character */
				*p++ = c;
			}
		}
		
		if ( result ) {
			/* May have run out of space before end condition */
			*p = '\0';
			result = str;
		}
	}
	
	return ( result );
}
