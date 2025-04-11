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
*			Module Name:		migrate.c
*			Revision Number:	1.2
*			Date last revised:	3 June 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	25 January 1995
*
*      migrate -- move a text file from a file system that has one end of line
*	  			format to a file system that has a different format
*
*		SYNOPSIS
*   		   migrate [-t u|m] file [file2]
*
*		DESCRIPTION
*    	  "Migrate" really justs translates carriage returns to new lines, or vice
*	 	 versa. It is called migrate, because if you pass the second file argument,
*	 	 the original file is left alone and the "translated" data is written
*	 	 to file2.
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
*********************************************************************************/

#include 	<stdio.h>
#include 	<string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbMemory.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbError.h"
#include "LbInterface.h"

/* Local Constants */
#define BUFF_SIZE		2048
#define MAC_EOL			13
#define UNIX_EOL		10


/* Local Macros */
#define Is_GoToMac()		LbFgIsSet(runTimeOptions, LBFlag0)	/* Did user specify to go to macintosh format? */
#define Is_GoToUnix()		LbFgIsSet(runTimeOptions, LBFlag1)	/* Did user specify to go to unix format? */
#define Is_InPlace()		LbFgIsSet(runTimeOptions, LBFlag2)	/* Did user specify convert in place? */

Boolean	migrate(int argc, char *argv[]);

	
/* FUNCTIONS */
/**********************
*	migrate
*
*	Purpose:	Execute the program.
*
*	Result:	Always returns zero indicating no error.
***********************/

Boolean migrate(int argc, char *argv[]) {

	Boolean			okay = false;		/* Process Flag */
	LbUsFourByte	localArgc;			/* Unsigned version of argc */
	char			inBuff[BUFF_SIZE];	/* Our input buffer */
	char			errStr[1024];		/* For error strings */
	FILE			*inputFile = 0;		/* The source file */
	FILE			*outputFile = 0;	/* The output file */
	LbUsFourByte	index;				/* Current char index */
	LbUsFourByte	numRead;			/* Current bytes read */
	
	/* The following variables are for command line options */
	#define	NUM_FLAGS	3
	char	*knownOptions[] = {"mui"};
	char				optArgs[NUM_FLAGS][LBEnMxArgLen];
	LbUsFourByte		optArgFlags = LBFlag0 + LBFlag1 + LBFlag2;
	LbUsFourByte		argIndex;
	LbUsFourByte		runTimeOptions=0;	/* The runtime options specified */


		
		
	do { /* Process Loop */
	
		/* Get our runtime options */
		if (!LbEnGetOptions(argc, argv, knownOptions,
				&runTimeOptions, optArgs, optArgFlags, &argIndex)) {
	
			ErStGeneric("Unable to get run time options from LbEnGetOptions");
			break;
		}
	
		/* Run through the options, verifying they are all there (that are necessary) */
		{
			
			/* Verify they specified an end of line format */
			if (!Is_GoToMac() && !Is_GoToUnix()){
				ErStGeneric("\nYou must specify a destination end of line format, -m for Macintosh or -u for UNIX\n");
				break;
			}
						
			/* Verify that they did not supply more than two files if translate in place is
				not selected; this indicates they don't understand what they are doing.
			*/
			if (!Is_InPlace() && (((argc - argIndex) > 2) || ((argc - argIndex) == 0))) {
				ErStGeneric("\nWrong number of arguments, when not converting in place you can only have a source file and a destination file.\n");
				break;
			}
						
			/* Verify that there is a destination file if there is no "in-place" request.
			*/
			if (!Is_InPlace() && (((argc - argIndex) != 2))) {
				ErStGeneric("\nWrong number of arguments, you must specify a destination file if not migrating in-place (-i).\n");
				break;
			}
		}
		
		/* Process the files on the command line */
		if (argc < 0) {
			localArgc = 0;
		}
		else {
			localArgc = argc;
		}
		while (argIndex < localArgc) {
			
			/* Open files according to "in place" options */
			if (Is_InPlace()){
			
				/* Open input file with read/write priviledge */
				if ((inputFile = LbFlFileOpen(argv[argIndex], "r+b")) == 0){
					sprintf(errStr,"\nUnable to open file '%s'\n",argv[argIndex]);
					ErStFileError(errStr);
					goto FAIL;
				}
				
				/* Set output file = to input file */
				outputFile = inputFile;
				
				/* Increment argIndex */
				argIndex++;
			}
			else {
			
				/* Open input file with read privilidge */
				if ((inputFile = LbFlFileOpen(argv[argIndex], "rb")) == 0){
					sprintf(errStr,"\nUnable to open file '%s'\n",argv[argIndex]);
					ErStFileError(errStr);
					goto FAIL;
				}
				
				/* Go to next command line argument */
				argIndex++;
				
				/* Create/Open output file */
				if ((outputFile = LbFlFileOpen(argv[argIndex], "wb+")) == 0) {
					sprintf(errStr, "\nUnable to create/open destination file '%s'\n",argv[argIndex]);
					ErStFileError(errStr);
					goto FAIL;
				}
				
				/* Go to next command line argument */
				argIndex++;
			} 
		
			/* Turn buffering off of input file */
			setbuf(inputFile, 0);
			
			/* Process current file */
			do {
			
				/* Read a block of data */
				numRead = fread(inBuff, 1, BUFF_SIZE, inputFile);
				
				/* Process the block if it was read properly */
				if (numRead != 0) {
				
					/* Filter the data */
					for (index = 0; index < numRead; index++){
						if (Is_GoToMac() && (inBuff[index] == UNIX_EOL))
							inBuff[index] = MAC_EOL;
						else if (Is_GoToUnix() && (inBuff[index] == MAC_EOL))
							inBuff[index] = UNIX_EOL;
					}
					
					/* If in place, seek back  */
					if (Is_InPlace()) {
						
						/* Seek back a block */
						if (fseek(inputFile, -numRead, SEEK_CUR) != 0) {
							sprintf(errStr,"\nUnable to seek backwards in file '%s'\n", argv[argIndex-1]);
							ErStFileError(errStr);
							goto FAIL;
						}
					}

					/* Write the filtered block */
					if (fwrite(inBuff, 1, numRead, outputFile) != numRead) {
						sprintf(errStr, "\nUnable to write filtered data to file '%s'\n", argv[argIndex-1]);
						ErStFileError(errStr);
						goto FAIL;
					}

					/* Flush the output, although you might think this isn't necessary
						you would be surprised
				    */
					fflush(outputFile);
				}
				
				/* Verify that we reached end of file, and not an error */
				#ifndef MPW
				if ((numRead == 0) && (feof(inputFile) == 0)){
				
					/* We were unable to read, but are not at end of file */
					ErStFileError("\nUnable to read from the input file\n");
					goto FAIL;
				}
				#endif
				
			#ifndef MPW
			} while (feof(inputFile) == 0);
			#else
			} while (numRead != 0);
			#endif
			
			/* Close input files */
			{
				if (inputFile != 0) {
					fclose(inputFile);
					inputFile = 0;
				}
				
				if (outputFile != 0) {
					fclose(outputFile);
					outputFile = 0;
				}
			}
		
		} /* End of loop processing files */
		
		okay = true;
		FAIL:;
	} while (false);
	
	
	/* Close input files if, they were left open */
	if (inputFile != 0) {
		fclose(inputFile);
		inputFile = 0;
	}
	
	if (outputFile != 0) {
		fclose(outputFile);
		outputFile = 0;
	}
	
	
	return (okay);
}
