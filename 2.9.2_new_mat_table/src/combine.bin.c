/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1994-2012 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		combine.bin.c
*			Revision Number:	1.3
*			Date last revised:	22 October 2012
*			Programmer:			Steven Vannoy
*			Date Originated:	15 December, 1994
*
*			Module Overview:	Combines a series of history files into one.
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
*			Programmer(s):		Karl G. Baum
*				
*			Revision date:		31 March, 2007
*			Revision description: Fixed several bugs that prevented the utility
*							from correctly combining more than 2 histograms.
*							This also addressed the corrupt header bug.
*
*********************************************************************************/

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbMemory.h"
#include "LbInterface.h"
#include "LbParamFile.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhoHFile.h"
#include "PhgHdr.h"
#include "phg.h"
#include "PhgBin.h"

/* LOCAL CONSTANTS */

/* This dictates the size of the units we read from disk.
It will be applied to a "generic" buffer which then gets
interpreted based on the type of the binned data. Hence,
it must remain a power of 2 >= 8.
*/
#define BUFF_SIZE 2048

/* LOCAL TYPES */

/* LOCAL GLOBALS */
static	Boolean	canceled;	/* Global cancelation flag */ 
static	PhoHFileHdrTy		newHeader;				/* Our Header */
static	PhoHFileHdrTy		finalHeader;			/* Header of final file */

/* PROTOTYPES */
Boolean			CombineBin(int argc, char *argv[]);

/* FUNCTIONS */
	

/**********************
*	CombineBin
*
*	Purpose:	Combine the history files.
*
*	Result:	True unless an error occurs.
***********************/
Boolean CombineBin(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */
	Boolean				deleteInputs = false;	/* Will we delete input files */
	char				errStr[1024];			/* Error string buffer */
	double				finalWeightRatio;		/* Ratio for scaling of pre-existing images */
	double				finalWeightSquRatio;	/* Ratio for scaling of pre-existing images */
	double				newWeightRatio;			/* Ratio for scaling of pre-existing images */
	double				newWeightSquRatio;		/* Ratio for scaling of pre-existing images */
	FILE				*inputFile;				/* Current input file */
	FILE				*finalFile;				/* The output file */
	LbUsFourByte		index;					/* Current element in buffer */
	LbUsFourByte		curFileIndex;			/* Current file index */
	LbUsFourByte		numToProcess;			/* Number of files to process */
	LbUsFourByte		elemsPerBuff;			/* Number of elements per buffer */
	LbUsFourByte		elemSize;				/* Size of buffer elements */
	LbUsFourByte		numRead;				/* Number of bytes read from final file */
	LbHdrHkTy			headerHk;				/* The header hook */
	void				*finalData;				/* Buffer for result of sum */
	void				*newData;				/* Buffer for data being added */
	

	do { /* Process Loop */
			
		/* Verify command line contains two arguments */
		if (argc < 3) {
			/* Ask for the file name */
			ErAbort("\nThis program requires at least 2 file names for input.\n"
				"The result of the combining is stored in the first argument to the program.\n"
				"Hence you must save a copy of it if you want your original file afterwards.");
		}
		
		/* See if we will delete input files on the fly */
		if (LbInAskYesNo("Do you want the input files deleted", LBINYes)
					== LBINYes) {
			deleteInputs = true;
		}
		
		/* Open the first file, the one that will have the final content */
		if ((finalFile = LbFlFileOpen(argv[1], "r+b")) == 0) {
			sprintf(errStr, "Unable to open output file\n'%s'.", argv[1]);
			ErStFileError(errStr);
			break;
		}
		
		/* Turn buffering off, it shouldn't be necessary but bugs have been found in the
			DG i/o library that make it so.
		*/
		setbuf(finalFile, 0);
				
		/* With the path to the PHG output file, get a hook to
			the file's header and an initialized run time parameters
			structure
		*/
		if (!PhgHdrGtParams(finalFile, &finalHeader, &headerHk)) {
			break;
		}
		
		/* Verify header is from the current version */
		if (finalHeader.H.HdrVersion != PHG_HDR_HEADER_VERSION) {
			sprintf(errStr, "File '%s' has a header version of %3.2f\n"
				"the current version is %3.2f.\n"
				"You will need a different version of combine.bin to perform this operation.\n",
				argv[1], finalHeader.H.HdrVersion, PHG_HDR_HEADER_VERSION);
			ErStGeneric(errStr);
			break;
		}
		
		/* Compute the number of elements to be processed per buffer */
		/* This requires parsing the type of the data and then the size of the elements */
		switch(finalHeader.H.HdrKind) {
		
			/* Data is a weight image */
			case PhoHFileEn_BIN_WT:
			
				/* Determine what size of real numbers are used */
				switch(finalHeader.H.BinRunTimeParams.weight_image_type) {

					/* Data is four byte reals */
					case PHG_BIN_WEIGHT_TYPE_R4:

						/* Buffer elements are 4 byte reals */
						elemSize = 4;
						elemsPerBuff = BUFF_SIZE/elemSize;
						break;

					/* Data is eight byte reals */
					case PHG_BIN_WEIGHT_TYPE_R8:

						/* Buffer elements are 8 byte reals */
						elemSize = 8;
						elemsPerBuff = BUFF_SIZE/elemSize;
						break;
					
					default:
						sprintf(errStr, "This file has an unknown weight_image_type.\n");
						ErStGeneric(errStr);
						goto FAIL;
				}
				
				break;
			
			/* Data is a weight squared image */	
			case PhoHFileEn_BIN_WTSQ:
			
				/* Determine what size of real numbers are used */
				switch(finalHeader.H.BinRunTimeParams.weight_image_type) {

					/* Data is four byte reals */
					case PHG_BIN_WEIGHT_TYPE_R4:

						/* Buffer elements are 4 byte reals */
						elemSize = 4;
						elemsPerBuff = BUFF_SIZE/elemSize;
						break;

					/* Data is eight byte reals */
					case PHG_BIN_WEIGHT_TYPE_R8:

						/* Buffer elements are 8 byte reals */
						elemSize = 8;
						elemsPerBuff = BUFF_SIZE/elemSize;
						break;
					
					default:
						sprintf(errStr, "This file has an unknown weight_image_type.\n");
						ErStGeneric(errStr);
						goto FAIL;
				}
			
				break;
			
			/* Data is a count image */
			case PhoHFileEn_BIN_CT:
			
				/* Determine what size of integers are used */
				switch(finalHeader.H.BinRunTimeParams.count_image_type) {

					/* Data is one byte integers */
					case PHG_BIN_COUNT_TYPE_I1:

						/* Buffer elements are 1 byte integers */
						elemSize = 1;
						elemsPerBuff = BUFF_SIZE;
						break;

					/* Data is two byte integers */
					case PHG_BIN_COUNT_TYPE_I2:

						/* Buffer elements are 2 byte integers */
						elemSize = 2;
						elemsPerBuff = BUFF_SIZE/elemSize;
						break;

					/* Data is four byte integers */
					case PHG_BIN_COUNT_TYPE_I4:

						/* Buffer elements are 4 byte integers */
						elemSize = 4;
						elemsPerBuff = BUFF_SIZE/elemSize;
						break;

					
					default:
						sprintf(errStr, "This file has an unknown weight_image_type.\n");
						ErStGeneric(errStr);
						goto FAIL;
				}
			
				break;
				
			default:
				sprintf(errStr, "This file has an uknown header type.\n");
				ErStGeneric(errStr);
				goto FAIL;
		}
		
		/* Allocate memory buffers */
		if ((finalData = LbMmAlloc(BUFF_SIZE)) == 0) {
			goto FAIL;
		}
		
		if ((newData = LbMmAlloc(BUFF_SIZE)) == 0) {
			goto FAIL;
		}
		
		/* Compute loop controls */
		curFileIndex = 2;
		numToProcess = argc - 1;
		
		/* Open first input file */
		if ((inputFile = LbFlFileOpen(argv[curFileIndex], "rb")) == 0) {
			sprintf(errStr, "Unable to open input file\n'%s'.", argv[curFileIndex]);
			ErStFileError(errStr);
			break;
		}

		/* Read the header in the first input file */
		{
		
			/* With the path to the PHG output file, get a hook to
				the file's header and an initialized run time parameters
				structure
			*/
			if (!PhgHdrGtParams(finalFile, &newHeader, &headerHk)) {
				break;
			}
			
			/* Verify header is from the current version */
			if (newHeader.H.HdrVersion != PHG_HDR_HEADER_VERSION) {
				sprintf(errStr, "File '%s' has a header version of %3.2f\n"
					"the current version is %3.2f.\n"
					"You will need a different version of combine.bin to perform this operation.\n",
					argv[curFileIndex], newHeader.H.HdrVersion, PHG_HDR_HEADER_VERSION);
				ErStGeneric(errStr);
				break;
			}
		}

			
		/* Process the files */
		for (curFileIndex = 2; curFileIndex <= numToProcess; curFileIndex++){

			/* Compute ratios for existing data, only applied to weight images */
			{
				/* The weight ratio is proportional to the events simulated in the
					final image, and the number of events simulated in the new image.
				*/
				/* Update: 3/31/2007 Karl G. Baum
				    Changed to use event count that was getting updated during each
					iteration of the for loop. */
				finalWeightRatio = (double) finalHeader.H.PhgRunTimeParams.Phg_EventsToSimulate/
					    (finalHeader.H.PhgRunTimeParams.Phg_EventsToSimulate + newHeader.H.PhgRunTimeParams.Phg_EventsToSimulate);
				/* finalWeightRatio = (double) finalHeader.H.SumEventsToSimulate/
						(finalHeader.H.SumEventsToSimulate + newHeader.H.SumEventsToSimulate);*/
				
				/* The weight squared ratio is is simply the weightRatio squared */
				finalWeightSquRatio = finalWeightRatio * finalWeightRatio;
			}
			
			/* Compute ratios for existing data, only applied to weight images */
			{
				/* The weight ratio is proportional to the events simulated in the
					final image, and the number of events simulated in the new image.
				*/
				/* Update: 3/31/2007 Karl G. Baum
				    Changed to use event count that was getting updated during each
					iteration of the for loop. */
				newWeightRatio = (double) newHeader.H.PhgRunTimeParams.Phg_EventsToSimulate/
					    (finalHeader.H.PhgRunTimeParams.Phg_EventsToSimulate + newHeader.H.PhgRunTimeParams.Phg_EventsToSimulate);
				/* newWeightRatio = (double) newHeader.H.SumEventsToSimulate/
						(finalHeader.H.SumEventsToSimulate + newHeader.H.SumEventsToSimulate);*/
				
				/* The weight squared ratio is is simply the weightRatio squared */
				newWeightSquRatio = newWeightRatio * newWeightRatio;
			}
			
			/* Update: 3/31/07 Karl G. Baum
			    Reset to beginning of the histogram data.  This addresses the problem where
				only the first two histograms get merged, because we do not go back to the
				beginning of the output file for the merging of additional histograms.
			*/
   			if (fseek(finalFile, PHG_HDR_HEADER_SIZE, SEEK_SET) != 0) {
   				ErStFileError("\nUnable to reset to beginning final data file.");
   				break;
   			}

			/* Loop by reading a buffer from the final file and processing until there
				isn't a full buffer left (the partial buffer is handled next)
			 */
			while ((numRead = fread(finalData, 1, BUFF_SIZE, finalFile)) == BUFF_SIZE) {
			
				/* Read the current "new" file */
				if (fread(newData, 1, BUFF_SIZE, inputFile) != BUFF_SIZE) {
					sprintf(errStr, "Error reading files, appearantly not the same size!\n");
					ErStFileError(errStr);
					goto FAIL;
				}
				
				/* Loop through, adding elements together */
				for (index = 0; index < elemsPerBuff; index++) {

					/* The math is done based on the type of data */
					switch(finalHeader.H.HdrKind) {
					
						/* Handle a weight image */
						case PhoHFileEn_BIN_WT:
						
							/* Determine what size of real number was used */
							switch(finalHeader.H.BinRunTimeParams.weight_image_type) {
		
								/* Data is four byte real */
								case PHG_BIN_WEIGHT_TYPE_R4:
									
									/* First scale the existing data */
									((float *)finalData)[index] = (finalWeightRatio *
										((float *)finalData)[index])
									 	+ (newWeightRatio * ((float *)newData)[index]);
									
									break;
			
								/* Data is eight byte real */
								case PHG_BIN_WEIGHT_TYPE_R8:
			
									/* First scale the existing data */
									((double *)finalData)[index] = (finalWeightSquRatio *
										((double *)finalData)[index])
										+ (newWeightSquRatio * ((double *)newData)[index]);

									break;
							}
							break;
								
						/* Handle a weight squared image */
						case PhoHFileEn_BIN_WTSQ:
						
							/* Determine what size of real number was used */
							switch(finalHeader.H.BinRunTimeParams.weight_image_type) {
			
								/* Data is four byte real */
								case PHG_BIN_WEIGHT_TYPE_R4:
									
									/* First scale the existing data */
									((float *)finalData)[index] = (finalWeightRatio *
										((float *)finalData)[index])
									 	+ (newWeightRatio * ((float *)newData)[index]);
									
									break;
			
								/* Data is eight byte real */
								case PHG_BIN_WEIGHT_TYPE_R8:
			
									/* First scale the existing data */
									((double *)finalData)[index] = (finalWeightSquRatio *
										((double *)finalData)[index])
										+ (newWeightSquRatio * ((double *)newData)[index]);

									break;
							}
							break;
						
						/* Handle a count image */
						case PhoHFileEn_BIN_CT:
						
							/* Determine what size of integer was used */
							switch(finalHeader.H.BinRunTimeParams.count_image_type) {
			
								/* Data is one byte integer */
								case PHG_BIN_COUNT_TYPE_I1:
									
									/* Update the image */
									((LbUsOneByte *)finalData)[index] +=
										((LbUsOneByte *)newData)[index];
									
									break;
			
								/* Data is two byte integer */
								case PHG_BIN_COUNT_TYPE_I2:
									
									/* Update the image */
									((LbUsTwoByte *)finalData)[index] +=
										((LbUsTwoByte *)newData)[index];

									break;
				
								/* Data is four byte integer */
								case PHG_BIN_COUNT_TYPE_I4:
									
									/* Update the image */
									((LbUsFourByte *)finalData)[index] +=
										((LbUsFourByte *)newData)[index];
									break;
							}
							break;
						
						/* All other cases are an error */
						case PhoHFileEn_PHGOLD:
						case PhoHFileEn_COLOLD:
						case PhoHFileEn_DETOLD:
						case PhoHFileEn_BIN:
						case PhoHFileEn_PHG2625:
						case PhoHFileEn_COL2625:
						case PhoHFileEn_DET2625:
						case PhoHFileEn_PHG:
						case PhoHFileEn_COL:
						case PhoHFileEn_DET:
							
							sprintf(errStr, "combinebin only works with SimSET''s binned output files.\n");
							ErStFileError(errStr);
							
							break;
					}
				}
				
				/* Seek back one buffer in final file */
				if (fseek(finalFile, -BUFF_SIZE, SEEK_CUR) != 0) {
					sprintf(errStr, "Unable to seek back in final file.\n");
					ErStFileError(errStr);
					goto FAIL;
				}
				
				/* Write the data to the output file */
				if (fwrite(finalData, BUFF_SIZE, 1, finalFile) != 1) {
					ErStFileError("\nUnable to write updated data to final file.");
					goto FAIL;
				}
			}
			
			/* Process partial buffer if there is  one */
			if (numRead != 0) {
				
				/* Recalculate the number of elements in this buffer */
				elemsPerBuff = numRead/elemSize;
				
			
				/* Read the current "new" file */
				if (fread(newData, 1, numRead, inputFile) != numRead) {
					sprintf(errStr, "Error reading files, appearantly not the same size!\n");
					ErStFileError(errStr);
					goto FAIL;
				}

				/* Loop through, adding elements together */
				for (index = 0; index < elemsPerBuff; index++) {
				
					/* The math is done based on the type of data */
					switch(finalHeader.H.HdrKind) {
					
						/* Handle a weight image */
						case PhoHFileEn_BIN_WT:
						
							/* Determine what size of real number was used */
							switch(finalHeader.H.BinRunTimeParams.weight_image_type) {
			
								/* Data is four byte real */
								case PHG_BIN_WEIGHT_TYPE_R4:
									
									/* First scale the existing data */
									((float *)finalData)[index] = (finalWeightRatio *
										((float *)finalData)[index])
									 	+ (newWeightRatio * ((float *)newData)[index]);
									
									break;
				
								/* Data is eight byte real */
								case PHG_BIN_WEIGHT_TYPE_R8:
			
									/* First scale the existing data */
									((double *)finalData)[index] = (finalWeightSquRatio *
										((double *)finalData)[index])
										+ (newWeightSquRatio * ((double *)newData)[index]);
										
									break;
							}
							break;
						
						/* Handle a weight squared image */
						case PhoHFileEn_BIN_WTSQ:
						
							/* Determine what size of real number was used */
							switch(finalHeader.H.BinRunTimeParams.weight_image_type) {

								/* Data is four byte real */
								case PHG_BIN_WEIGHT_TYPE_R4:
									
									/* First scale the existing data */
									((float *)finalData)[index] = (finalWeightRatio *
										((float *)finalData)[index])
									 	+ (newWeightRatio * ((float *)newData)[index]);
									
									break;
						
								/* Data is eight byte real */
								case PHG_BIN_WEIGHT_TYPE_R8:
			
									/* First scale the existing data */
									((double *)finalData)[index] = (finalWeightSquRatio *
										((double *)finalData)[index])
										+ (newWeightSquRatio * ((double *)newData)[index]);
									break;
							}
						
							break;
						
						/* Handle a count image */
						case PhoHFileEn_BIN_CT:
						
							/* Determine what size of integer was used */
							switch(finalHeader.H.BinRunTimeParams.count_image_type) {
			
								/* Data is one byte integer */
								case PHG_BIN_COUNT_TYPE_I1:
									
									/* Update the image */
									((LbUsOneByte *)finalData)[index] +=
										((LbUsOneByte *)newData)[index];
									
									break;
			
								/* Data is two byte integer */
								case PHG_BIN_COUNT_TYPE_I2:
									
									/* Update the image */
									((LbUsTwoByte *)finalData)[index] +=
										((LbUsTwoByte *)newData)[index];

									break;
			
								/* Data is four byte integer */
								case PHG_BIN_COUNT_TYPE_I4:
									
									/* Update the image */
									((LbUsFourByte *)finalData)[index] +=
										((LbUsFourByte *)newData)[index];
									break;
							}
						
							break;
						
						/* All other cases are an error */
						case PhoHFileEn_PHGOLD:
						case PhoHFileEn_COLOLD:
						case PhoHFileEn_DETOLD:
						case PhoHFileEn_BIN:
						case PhoHFileEn_PHG2625:
						case PhoHFileEn_COL2625:
						case PhoHFileEn_DET2625:
						case PhoHFileEn_PHG:
						case PhoHFileEn_COL:
						case PhoHFileEn_DET:
							
							sprintf(errStr, "combinebin only works with SimSET''s binned output files.\n");
							ErStFileError(errStr);
							
							break;
					}
				}
				
				/* Seek back one buffer in final file */
				if (fseek(finalFile, -numRead, SEEK_CUR) != 0) {
					sprintf(errStr, "Unable to seek back in final file.\n");
					ErStFileError(errStr);
					goto FAIL;
				}
				
				/* Write the data to the output file */
				if (fwrite(finalData, numRead, 1, finalFile) != 1) {
					ErStFileError("\nUnable to write updated data to final file.");
					goto FAIL;
				}
			}
			/* If we are here and haven't read to EOF there is an error */
			else if (feof(finalFile) == 0) {
				sprintf(errStr, "Unable to read from file '%s'.\n", argv[curFileIndex]);
				ErStFileError(errStr);
				goto FAIL;
			}
			
			/* Close the input file */
			fclose(inputFile);
			
			/* Delete input file if requested */
			if (deleteInputs == true) {
				if (remove(argv[curFileIndex]) != 0) {
					sprintf(errStr, "Unable to delete file '%s'.\n", argv[curFileIndex]);
					ErStFileError(errStr);
					goto FAIL;
				}
			}
			
			/* Update header */
			finalHeader.H.NumPhotons += newHeader.H.NumPhotons;
			finalHeader.H.NumDecays += newHeader.H.NumDecays;
			finalHeader.H.PhgRunTimeParams.Phg_EventsToSimulate +=
				newHeader.H.PhgRunTimeParams.Phg_EventsToSimulate;
			
			/* Open next file */
			if (curFileIndex < numToProcess) {
								
				/* Open input file */
				if ((inputFile = LbFlFileOpen(argv[curFileIndex+1], "rb")) == 0) {
					sprintf(errStr, "Unable to open input file\n'%s'.", argv[curFileIndex+1]);
					ErStFileError(errStr);
					break;
				}
			}

		} /* End FOR each input file */

   		/* Reset to zero and write the header */
   		if (fseek(finalFile, 0, SEEK_SET) != 0) {
   			ErStFileError("\nUnable to reset to beginning final data file.");
   			break;
   		}

		/* Write the header to the output file */
		/* Update: 3/31/07 Karl G. Baum
			 Fix bug causing combined histograms to have corrupt headers.
		*/
		if (PhgHdrUpHeader(finalFile, &finalHeader, &headerHk) == false) {
			ErStFileError("\nUnable to write updated header to final data file.");
			break;
		}
		/* if (fwrite(&finalHeader, sizeof(finalHeader), 1, finalFile) != 1) {
			ErStFileError("\nUnable to write updated header to final data file.");
			break;
		}*/

		/* Close the final file */
		fclose(finalFile);
		
		okay = true;
		FAIL:;
	} while (false);
	
	/* If error set due to cancellation, handle it, otherwise pass it on */
	if (!okay) {
		if (canceled) {
			ErHandle("User canceled combine.bin.", false);
			okay = true;
		}
	}
	
	/* Return the status */
	return (okay);
}
