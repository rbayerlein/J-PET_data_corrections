/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1994-2001 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		phg.swap.c
*			Revision Number:	1.0
*			Date last revised:	2007
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
*			Programmer(s):		
*
*			Revision date:		
*			Revision description:
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
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*			Revision description:
*						- support for randoms and eight-byte number of decays
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
#include "LbConvert.h"
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
static	PhoHFileHdrTy		nwHdr;				/* Our Header */
static	PhoHFileHdrTy		fhdr;			/* Header of final file */

/* PROTOTYPES */
Boolean			PhgSwap(int argc, char *argv[]);

/* FUNCTIONS */

/**********************
*	PhgSwap
*
*	Purpose:	Combine the history files.
*
*	Result:	True unless an error occurs.
***********************/
Boolean PhgSwap(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */
	char				errStr[1024];			/* Error string buffer */
	FILE				*finalFile;				/* The output file */
	LbUsFourByte		index;					/* Current element in buffer */
	LbUsFourByte		fileIndex;				/* Current element in buffer */
	LbUsFourByte		numToProcess;			/* Number of files to process */
	LbUsFourByte		elemSize;				/* Size of buffer elements */
	LbUsFourByte		elemsPerBuff;			/* Number of "elements" per data buffer read */
	LbFourByte			numRead;				/* Number of bytes read from final file */
	LbHdrHkTy			headerHk;				/* The header hook */
	void				*finalData;				/* Buffer for result of sum */
	void				*newData;				/* Buffer for data being added */
	
	
	do { /* Process Loop */
			
		/* Verify command line contains two arguments */
		if (argc < 2) {
			/* Ask for the file name */
			ErAbort("\nThis program requires a file name for input.\n");
		}
		
		/* Loop through all files on the command line */
		numToProcess = argc - 1;
		for (fileIndex = 1; fileIndex <= numToProcess; fileIndex++) {
		
			/* Open the first file, the one that will have the final content */
			if ((finalFile = LbFlFileOpen(argv[fileIndex], "r+b")) == 0) {
				sprintf(errStr, "Unable to open output file\n'%s'.", argv[fileIndex]);
				ErStFileError(errStr);
				goto FAIL;
			}
			
			/* Turn buffering off, it shouldn't be necessary but bugs have been found in the
				DG i/o library that make it so.
			*/
			setbuf(finalFile, 0);
			
			/* With the path to the PHG output file, get a hook to
				the file's header and an initialized run time parameters
				structure
			*/
			if (!PhgHdrGtParams(finalFile, &fhdr, &headerHk)) {
				break;
			}
			
			/* Verify header is from the current version */
			if (fhdr.H.HdrVersion != PHG_HDR_HEADER_VERSION) {
				sprintf(errStr, "File '%s' has a header version of %3.2f\n"
					"the current version is %3.2f.\n"
					"You will need a different version of phg.swap to perform this operation.\n",
					argv[1], fhdr.H.HdrVersion, PHG_HDR_HEADER_VERSION);
				ErStGeneric(errStr);
				goto FAIL;
			}
			
			/* Compute the number of elements to be processed per buffer */
			/* This requires parsing the type of the data and then the size of the elements */
			switch(fhdr.H.HdrKind) {
			
				/* Data is a weight image */
				case PhoHFileEn_BIN_WT:
				
					/* Determine what size of real numbers are used */
					switch(fhdr.H.BinRunTimeParams.weight_image_type) {
	
						/* Data is four byte reals */
						case PHG_BIN_WEIGHT_TYPE_R4:
	
							/* Buffer elements are 4 byte reals */
							elemSize = 4;
							break;
	
						/* Data is eight byte reals */
						case PHG_BIN_WEIGHT_TYPE_R8:
	
							/* Buffer elements are 8 byte reals */
							elemSize = 8;
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
					switch(fhdr.H.BinRunTimeParams.weight_image_type) {
	
						/* Data is four byte reals */
						case PHG_BIN_WEIGHT_TYPE_R4:
	
							/* Buffer elements are 4 byte reals */
							elemSize = 4;
							break;
	
						/* Data is eight byte reals */
						case PHG_BIN_WEIGHT_TYPE_R8:
	
							/* Buffer elements are 8 byte reals */
							elemSize = 8;
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
					switch(fhdr.H.BinRunTimeParams.count_image_type) {
	
						/* Data is one byte integers */
						case PHG_BIN_COUNT_TYPE_I1:
	
							/* Buffer elements are 1 byte integers */
							elemSize = 1;
							break;
	
						/* Data is two byte integers */
						case PHG_BIN_COUNT_TYPE_I2:
	
							/* Buffer elements are 2 byte integers */
							elemSize = 2;
							break;
	
						/* Data is four byte integers */
						case PHG_BIN_COUNT_TYPE_I4:
	
							/* Buffer elements are 4 byte integers */
							elemSize = 4;
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
			
			/* If element size != 1, then swap the file data */
			if (elemSize != 1) {
				/* Allocate memory buffers */
				if ((finalData = LbMmAlloc(BUFF_SIZE)) == 0) {
					goto FAIL;
				}
				
				if ((newData = LbMmAlloc(BUFF_SIZE)) == 0) {
					goto FAIL;
				}
				
				/* Compute elements per buffer */
				elemsPerBuff = BUFF_SIZE/elemSize;
				
				/* Loop by reading a buffer from the final file and processing until there
					isn't a full buffer left (the partial buffer is handled next)
				 */
				while ((numRead = fread(finalData, 1, BUFF_SIZE, finalFile)) == BUFF_SIZE) {
					/* Loop through, swapping element bytes */
					for (index = 0; index < elemsPerBuff; index++) {
					
						LbCvSwap((char *)finalData+(index*elemSize),
							(char *)newData+(index*elemSize), elemSize);
					}
					
					/* Seek back one buffer in final file */
					if (fseek(finalFile, -BUFF_SIZE, SEEK_CUR) != 0) {
						sprintf(errStr, "Unable to seek back in final file.\n");
						ErStFileError(errStr);
						goto FAIL;
					}
					
					/* Write the data to the output file */
					if (fwrite(newData, BUFF_SIZE, 1, finalFile) != 1) {
						ErStFileError("\nUnable to write updated data to final file.");
						goto FAIL;
					}
				}
				
				/* Process partial buffer if there is  one */
				if (numRead != 0) {
					
					/* Recalculate the number of elements in this buffer */
					elemsPerBuff = numRead/elemSize;
	
					/* Loop through, adding elements together */
					for (index = 0; index < elemsPerBuff; index++) {
						
						LbCvSwap((char *)finalData+(index*elemSize),
							(char *)newData+(index*elemSize), elemSize);
					}
					
					/* Seek back one buffer in final file */
					if (fseek(finalFile, -numRead, SEEK_CUR) != 0) {
						sprintf(errStr, "Unable to seek back in final file.\n");
						ErStFileError(errStr);
						goto FAIL;
					}
					
					/* Write the data to the output file */
					if (fwrite(newData, numRead, 1, finalFile) != 1) {
						ErStFileError("\nUnable to write updated data to final file.");
						goto FAIL;
					}
				}
				/* If we are here and haven't read to EOF there is an error */
				else if (feof(finalFile) == 0) {
					sprintf(errStr, "Unable to read from file '%s'.\n", argv[fileIndex]);
					ErStFileError(errStr);
					goto FAIL;
				}
				
			}
			
			/* Copy all of the data from the old header to the new, that way we get the one byte
				fields for free
			*/
			memcpy((void *)&nwHdr, (const void *)&fhdr, sizeof(fhdr));
			
			/* Now swap the bytes in the header, one field at a time */
			{
				LbCvSwap((char *)&fhdr.H.HdrSize, (char *)&nwHdr.H.HdrSize,
					sizeof(nwHdr.H.HdrSize));
				LbCvSwap((char *)&fhdr.H.HdrKind, (char *)&nwHdr.H.HdrKind,
					sizeof(nwHdr.H.HdrKind));
				LbCvSwap((char *)&fhdr.H.HdrVersion, (char *)&nwHdr.H.HdrVersion,
					sizeof(nwHdr.H.HdrVersion));
	
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.Phg_EventsToSimulate, (char *)&nwHdr.H.PhgRunTimeParams.Phg_EventsToSimulate,
					sizeof(nwHdr.H.PhgRunTimeParams.Phg_EventsToSimulate));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgRandomSeed, (char *)&nwHdr.H.PhgRunTimeParams.PhgRandomSeed,
					sizeof(nwHdr.H.PhgRunTimeParams.PhgRandomSeed));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.Phg_LengthOfScan, (char *)&nwHdr.H.PhgRunTimeParams.Phg_LengthOfScan,
					sizeof(nwHdr.H.PhgRunTimeParams.Phg_LengthOfScan));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.Phg_AcceptanceAngle, (char *)&nwHdr.H.PhgRunTimeParams.Phg_AcceptanceAngle,
					sizeof(nwHdr.H.PhgRunTimeParams.Phg_AcceptanceAngle));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.Phg_SineOfAcceptanceAngle, (char *)&nwHdr.H.PhgRunTimeParams.Phg_SineOfAcceptanceAngle,
					sizeof(nwHdr.H.PhgRunTimeParams.Phg_SineOfAcceptanceAngle));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgMinimumEnergy, (char *)&nwHdr.H.PhgRunTimeParams.PhgMinimumEnergy,
					sizeof(nwHdr.H.PhgRunTimeParams.PhgMinimumEnergy));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgMinWWRatio, (char *)&nwHdr.H.PhgRunTimeParams.PhgMinWWRatio,
					sizeof(nwHdr.H.PhgRunTimeParams.PhgMinWWRatio));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgMaxWWRatio, (char *)&nwHdr.H.PhgRunTimeParams.PhgMaxWWRatio,
					sizeof(nwHdr.H.PhgRunTimeParams.PhgMaxWWRatio));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgNuclide.isotope, (char *)&nwHdr.H.PhgRunTimeParams.PhgNuclide.isotope,
					sizeof(nwHdr.H.PhgRunTimeParams.PhgNuclide.isotope));	
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgNuclide.photonEnergy_KEV, (char *)&nwHdr.H.PhgRunTimeParams.PhgNuclide.photonEnergy_KEV,
					sizeof(nwHdr.H.PhgRunTimeParams.PhgNuclide.photonEnergy_KEV));		
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgNuclide.maximumPositronEnergy, (char *)&nwHdr.H.PhgRunTimeParams.PhgNuclide.maximumPositronEnergy,
					sizeof(nwHdr.H.PhgRunTimeParams.PhgNuclide.maximumPositronEnergy));
				
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsForcedDetection, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsForcedDetection,
					sizeof(nwHdr.H.PhgRunTimeParams.PhgIsForcedDetection));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsStratification, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsStratification,
					sizeof(nwHdr.H.PhgRunTimeParams.PhgIsStratification));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsForcedNonAbsorbtion, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsForcedNonAbsorbtion,
					sizeof(nwHdr.H.PhgRunTimeParams.PhgIsForcedNonAbsorbtion));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsSPECT, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsSPECT,
					sizeof(nwHdr.H.PhgRunTimeParams.PhgIsSPECT));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsPETCoincidencesOnly, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsPETCoincidencesOnly,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsPETCoincidencesOnly));		
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsPETCoincPlusSingles, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsPETCoincPlusSingles,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsPETCoincPlusSingles));		
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsMultiEmission, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsMultiEmission,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsMultiEmission));		
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsPET, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsPET,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsPET));		
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsPETCoincidences, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsPETCoincidences,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsPETCoincidences));		
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsPETSingles, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsPETSingles,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsPETSingles));		
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsHistoryFile, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsHistoryFile,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsHistoryFile));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsAdjForPosRange, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsAdjForPosRange,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsAdjForPosRange));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsAdjForCollinearity, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsAdjForCollinearity,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsAdjForCollinearity));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsComputedProductivityTbl, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsComputedProductivityTbl,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsComputedProductivityTbl));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsVoxelPointSource, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsVoxelPointSource,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsVoxelPointSource));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsBinOnTheFly, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsBinOnTheFly,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsBinOnTheFly));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsCollimateOnTheFly, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsCollimateOnTheFly,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsCollimateOnTheFly));
				LbCvSwap((char *)&fhdr.H.PhgRunTimeParams.PhgIsDetectOnTheFly, (char *)&nwHdr.H.PhgRunTimeParams.PhgIsDetectOnTheFly,
					sizeof(fhdr.H.PhgRunTimeParams.PhgIsDetectOnTheFly));
	
				LbCvSwap((char *)&fhdr.H.ColRunTimeParams.Initialized, (char *)&nwHdr.H.ColRunTimeParams.Initialized,
					sizeof(nwHdr.H.ColRunTimeParams.Initialized));
				LbCvSwap((char *)&fhdr.H.ColRunTimeParams.DoHistory, (char *)&nwHdr.H.ColRunTimeParams.DoHistory,
					sizeof(fhdr.H.ColRunTimeParams.DoHistory));
				LbCvSwap((char *)&fhdr.H.ColRunTimeParams.ColType, (char *)&nwHdr.H.ColRunTimeParams.ColType,
					sizeof(fhdr.H.ColRunTimeParams.ColType));
				
				/* Determine which type of collimator was used and convert its fields
				*/
				switch (fhdr.H.ColRunTimeParams.ColType) {
				
					case ColEn_simple_pet:
						LbCvSwap((char *)&fhdr.H.ColRunTimeParams.SimplePETCol.Depth, (char *)&nwHdr.H.ColRunTimeParams.SimplePETCol.Depth,
					sizeof(fhdr.H.ColRunTimeParams.SimplePETCol.Depth));
						break;
						
					case ColEn_monte_carlo_pet:
						
						/* Just note that the actual collimator description is not stored in the header. Hence,
							we don't have anything to do here.
						*/
						nwHdr.H.ColRunTimeParams.MCPETCol = 0;;
						break;
						
					case ColEn_simple_spect:
						break;
						
					case ColEn_unc_spect:
						
						LbCvSwap((char *)&fhdr.H.ColRunTimeParams.UNCSPECTCol.HoleGeometry, (char *)&nwHdr.H.ColRunTimeParams.UNCSPECTCol.HoleGeometry,
							sizeof(nwHdr.H.ColRunTimeParams.UNCSPECTCol.HoleGeometry));
						LbCvSwap((char *)&fhdr.H.ColRunTimeParams.UNCSPECTCol.RadiusOfRotation, (char *)&nwHdr.H.ColRunTimeParams.UNCSPECTCol.RadiusOfRotation,
							sizeof(nwHdr.H.ColRunTimeParams.UNCSPECTCol.RadiusOfRotation));
						LbCvSwap((char *)&fhdr.H.ColRunTimeParams.UNCSPECTCol.Thickness, (char *)&nwHdr.H.ColRunTimeParams.UNCSPECTCol.Thickness,
							sizeof(nwHdr.H.ColRunTimeParams.UNCSPECTCol.Thickness));
						LbCvSwap((char *)&fhdr.H.ColRunTimeParams.UNCSPECTCol.HoleRadius, (char *)&nwHdr.H.ColRunTimeParams.UNCSPECTCol.HoleRadius,
							sizeof(nwHdr.H.ColRunTimeParams.UNCSPECTCol.HoleRadius));
						LbCvSwap((char *)&fhdr.H.ColRunTimeParams.UNCSPECTCol.SeptalThickness, (char *)&nwHdr.H.ColRunTimeParams.UNCSPECTCol.SeptalThickness,
							sizeof(nwHdr.H.ColRunTimeParams.UNCSPECTCol.SeptalThickness));
						LbCvSwap((char *)&fhdr.H.ColRunTimeParams.UNCSPECTCol.MinZ, (char *)&nwHdr.H.ColRunTimeParams.UNCSPECTCol.MinZ,
							sizeof(nwHdr.H.ColRunTimeParams.UNCSPECTCol.MinZ));
						LbCvSwap((char *)&fhdr.H.ColRunTimeParams.UNCSPECTCol.MaxZ, (char *)&nwHdr.H.ColRunTimeParams.UNCSPECTCol.MaxZ,
							sizeof(nwHdr.H.ColRunTimeParams.UNCSPECTCol.MaxZ));
						LbCvSwap((char *)&fhdr.H.ColRunTimeParams.UNCSPECTCol.StartAngle, (char *)&nwHdr.H.ColRunTimeParams.UNCSPECTCol.StartAngle,
							sizeof(nwHdr.H.ColRunTimeParams.UNCSPECTCol.StartAngle));
						LbCvSwap((char *)&fhdr.H.ColRunTimeParams.UNCSPECTCol.StopAngle, (char *)&nwHdr.H.ColRunTimeParams.UNCSPECTCol.StopAngle,
							sizeof(nwHdr.H.ColRunTimeParams.UNCSPECTCol.StopAngle));
						LbCvSwap((char *)&fhdr.H.ColRunTimeParams.UNCSPECTCol.SumAllViews, (char *)&nwHdr.H.ColRunTimeParams.UNCSPECTCol.SumAllViews,
							sizeof(nwHdr.H.ColRunTimeParams.UNCSPECTCol.SumAllViews));
						LbCvSwap((char *)&fhdr.H.ColRunTimeParams.UNCSPECTCol.NumViews, (char *)&nwHdr.H.ColRunTimeParams.UNCSPECTCol.NumViews,
							sizeof(nwHdr.H.ColRunTimeParams.UNCSPECTCol.NumViews));
						break;
					default:
						break;
				}
				
				LbCvSwap((char *)&fhdr.H.DetRunTimeParams.Initialized, (char *)&nwHdr.H.DetRunTimeParams.Initialized,
					sizeof(nwHdr.H.DetRunTimeParams.Initialized));
				LbCvSwap((char *)&fhdr.H.DetRunTimeParams.DoHistory, (char *)&nwHdr.H.DetRunTimeParams.DoHistory,
					sizeof(nwHdr.H.DetRunTimeParams.DoHistory));
				
				/* Determine which type of detector was used and convert its fields
				*/
				switch (fhdr.H.DetRunTimeParams.DetectorType) {
				
					case DetEn_simple_pet:
						LbCvSwap((char *)&fhdr.H.DetRunTimeParams.EnergyResolutionPercentage, (char *)&nwHdr.H.DetRunTimeParams.EnergyResolutionPercentage,
					sizeof(nwHdr.H.DetRunTimeParams.EnergyResolutionPercentage));
						LbCvSwap((char *)&fhdr.H.DetRunTimeParams.ReferenceEnergy, (char *)&nwHdr.H.DetRunTimeParams.ReferenceEnergy,
					sizeof(nwHdr.H.DetRunTimeParams.ReferenceEnergy));
						break;
						
					case DetEn_simple_spect:
						break;
						
					case DetEn_unc_spect:
						break;
						
					default:
						break;
				}


				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.numZBins, (char *)&nwHdr.H.BinRunTimeParams.numPABins, 
					sizeof(nwHdr.H.BinRunTimeParams.numPABins));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.numPABins, (char *)&nwHdr.H.BinRunTimeParams.numPABins, 
					sizeof(nwHdr.H.BinRunTimeParams.numPABins));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.numTDBins, (char *)&nwHdr.H.BinRunTimeParams.numTDBins, 
					sizeof(nwHdr.H.BinRunTimeParams.numTDBins));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.numAABins, (char *)&nwHdr.H.BinRunTimeParams.numAABins, 
					sizeof(nwHdr.H.BinRunTimeParams.numAABins));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.numE1Bins, (char *)&nwHdr.H.BinRunTimeParams.numE1Bins, 
					sizeof(nwHdr.H.BinRunTimeParams.numE1Bins));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.numE2Bins, (char *)&nwHdr.H.BinRunTimeParams.numE2Bins, 
					sizeof(nwHdr.H.BinRunTimeParams.numE2Bins));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.numS1Bins, (char *)&nwHdr.H.BinRunTimeParams.numS1Bins, 
					sizeof(nwHdr.H.BinRunTimeParams.numS1Bins));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.numS2Bins, (char *)&nwHdr.H.BinRunTimeParams.numS2Bins, 
					sizeof(nwHdr.H.BinRunTimeParams.numS2Bins));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.numImageBins, (char *)&nwHdr.H.BinRunTimeParams.numImageBins, 
					sizeof(nwHdr.H.BinRunTimeParams.numImageBins));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.scatterRandomParam, (char *)&nwHdr.H.BinRunTimeParams.scatterRandomParam, 
					sizeof(nwHdr.H.BinRunTimeParams.scatterRandomParam));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.minZ, (char *)&nwHdr.H.BinRunTimeParams.minZ, 
					sizeof(nwHdr.H.BinRunTimeParams.minZ));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.maxZ, (char *)&nwHdr.H.BinRunTimeParams.maxZ, 
					sizeof(nwHdr.H.BinRunTimeParams.maxZ));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.minPA, (char *)&nwHdr.H.BinRunTimeParams.minPA, 
					sizeof(nwHdr.H.BinRunTimeParams.minPA));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.maxPA, (char *)&nwHdr.H.BinRunTimeParams.maxPA, 
					sizeof(nwHdr.H.BinRunTimeParams.maxPA));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.minTD, (char *)&nwHdr.H.BinRunTimeParams.minTD, 
					sizeof(nwHdr.H.BinRunTimeParams.minTD));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.maxTD, (char *)&nwHdr.H.BinRunTimeParams.maxTD, 
					sizeof(nwHdr.H.BinRunTimeParams.maxTD));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.minAA, (char *)&nwHdr.H.BinRunTimeParams.minAA, 
					sizeof(nwHdr.H.BinRunTimeParams.minAA));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.maxAA, (char *)&nwHdr.H.BinRunTimeParams.maxAA, 
					sizeof(nwHdr.H.BinRunTimeParams.maxAA));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.minE, (char *)&nwHdr.H.BinRunTimeParams.minE, 
					sizeof(nwHdr.H.BinRunTimeParams.minE));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.maxE, (char *)&nwHdr.H.BinRunTimeParams.maxE, 
					sizeof(nwHdr.H.BinRunTimeParams.maxE));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.minS, (char *)&nwHdr.H.BinRunTimeParams.minS, 
					sizeof(nwHdr.H.BinRunTimeParams.minS));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.maxS, (char *)&nwHdr.H.BinRunTimeParams.maxS, 
					sizeof(nwHdr.H.BinRunTimeParams.maxS));
				LbCvSwap((char *)& fhdr.H.BinRunTimeParams.addToExistingImg, (char *)&nwHdr.H.BinRunTimeParams.addToExistingImg,
					sizeof(nwHdr.H.BinRunTimeParams.addToExistingImg));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.doCounts, (char *)&nwHdr.H.BinRunTimeParams.doCounts,
					sizeof(nwHdr.H.BinRunTimeParams.doCounts));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.doWeights, (char *)&fhdr.H.BinRunTimeParams.doWeights,
					sizeof(nwHdr.H.BinRunTimeParams.doWeights));
				LbCvSwap((char *)&fhdr.H.BinRunTimeParams.doWeightsSquared, (char *)&fhdr.H.BinRunTimeParams.doWeightsSquared,
		sizeof(nwHdr.H.BinRunTimeParams.doWeightsSquared));
	
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.zRange, (char *)&nwHdr.H.BinRunTimeParams.zRange, 
		sizeof(nwHdr.H.BinRunTimeParams.zRange));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.eRange, (char *)&nwHdr.H.BinRunTimeParams.eRange, 
		sizeof(nwHdr.H.BinRunTimeParams.eRange));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.sRange, (char *)&nwHdr.H.BinRunTimeParams.sRange, 
		sizeof(nwHdr.H.BinRunTimeParams.sRange));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.tdRange, (char *)&nwHdr.H.BinRunTimeParams.tdRange, 
		sizeof(nwHdr.H.BinRunTimeParams.tdRange));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.scatter2CIsize, (char *)&nwHdr.H.BinRunTimeParams.scatter2CIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.scatter2CIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.scatter2WIsize, (char *)&nwHdr.H.BinRunTimeParams.scatter2WIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.scatter2WIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.scatter2WISsize, (char *)&nwHdr.H.BinRunTimeParams.scatter2WISsize, 
		sizeof(nwHdr.H.BinRunTimeParams.scatter2WISsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.scatter1CIsize, (char *)&nwHdr.H.BinRunTimeParams.scatter1CIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.scatter1CIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.scatter1WIsize, (char *)&nwHdr.H.BinRunTimeParams.scatter1WIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.scatter1WIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.scatter1WISsize, (char *)&nwHdr.H.BinRunTimeParams.scatter1WISsize, 
		sizeof(nwHdr.H.BinRunTimeParams.scatter1WISsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.energy2CIsize, (char *)&nwHdr.H.BinRunTimeParams.energy2CIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.energy2CIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.energy2WIsize, (char *)&nwHdr.H.BinRunTimeParams.energy2WIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.energy2WIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.energy2WISsize, (char *)&nwHdr.H.BinRunTimeParams.energy2WISsize, 
		sizeof(nwHdr.H.BinRunTimeParams.energy2WISsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.energy1CIsize, (char *)&nwHdr.H.BinRunTimeParams.energy1CIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.energy1CIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.energy1WIsize, (char *)&nwHdr.H.BinRunTimeParams.energy1WIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.energy1WIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.energy1WISsize, (char *)&nwHdr.H.BinRunTimeParams.energy1WISsize, 
		sizeof(nwHdr.H.BinRunTimeParams.energy1WISsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.aaCIsize, (char *)&nwHdr.H.BinRunTimeParams.aaCIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.aaCIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.aaWIsize, (char *)&nwHdr.H.BinRunTimeParams.aaWIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.aaWIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.aaWISsize, (char *)&nwHdr.H.BinRunTimeParams.aaWISsize, 
		sizeof(nwHdr.H.BinRunTimeParams.aaWISsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.tdCIsize, (char *)&nwHdr.H.BinRunTimeParams.tdCIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.tdCIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.tdWIsize, (char *)&nwHdr.H.BinRunTimeParams.tdWIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.tdWIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.tdWISsize, (char *)&nwHdr.H.BinRunTimeParams.tdWISsize, 
		sizeof(nwHdr.H.BinRunTimeParams.tdWISsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.paCIsize, (char *)&nwHdr.H.BinRunTimeParams.paCIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.paCIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.paWIsize, (char *)&nwHdr.H.BinRunTimeParams.paWIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.paWIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.paWISsize, (char *)&nwHdr.H.BinRunTimeParams.paWISsize, 
		sizeof(nwHdr.H.BinRunTimeParams.paWISsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.z2CIsize, (char *)&nwHdr.H.BinRunTimeParams.z2CIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.z2CIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.z2WIsize, (char *)&nwHdr.H.BinRunTimeParams.z2WIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.z2WIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.z2WISsize, (char *)&nwHdr.H.BinRunTimeParams.z2WISsize, 
		sizeof(nwHdr.H.BinRunTimeParams.z2WISsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.z1CIsize, (char *)&nwHdr.H.BinRunTimeParams.z1CIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.z1CIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.z1WIsize, (char *)&nwHdr.H.BinRunTimeParams.z1WIsize, 
		sizeof(nwHdr.H.BinRunTimeParams.z1WIsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.z1WISsize, (char *)&nwHdr.H.BinRunTimeParams.z1WISsize, 
		sizeof(nwHdr.H.BinRunTimeParams.z1WISsize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.countImageSize, (char *)&nwHdr.H.BinRunTimeParams.countImageSize, 
		sizeof(nwHdr.H.BinRunTimeParams.countImageSize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.weightImageSize, (char *)&nwHdr.H.BinRunTimeParams.weightImageSize, 
		sizeof(nwHdr.H.BinRunTimeParams.weightImageSize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.weightSquImageSize, (char *)&nwHdr.H.BinRunTimeParams.weightSquImageSize, 
		sizeof(nwHdr.H.BinRunTimeParams.weightSquImageSize));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.weight_image_type, (char *)&nwHdr.H.BinRunTimeParams.weight_image_type, 
		sizeof(nwHdr.H.BinRunTimeParams.weight_image_type));
	LbCvSwap((char *)&fhdr.H.BinRunTimeParams.count_image_type, (char *)&nwHdr.H.BinRunTimeParams.count_image_type, 
		sizeof(nwHdr.H.BinRunTimeParams.count_image_type));
				
				LbCvSwap((char *)&fhdr.H.TargetCylinder.radius, (char *)&nwHdr.H.TargetCylinder.radius,
					sizeof(nwHdr.H.TargetCylinder.radius));
				LbCvSwap((char *)&fhdr.H.TargetCylinder.zMin, (char *)&nwHdr.H.TargetCylinder.zMin,
					sizeof(nwHdr.H.TargetCylinder.zMin));
				LbCvSwap((char *)&fhdr.H.TargetCylinder.zMax, (char *)&nwHdr.H.TargetCylinder.zMax,
					sizeof(nwHdr.H.TargetCylinder.zMax));
				LbCvSwap((char *)&fhdr.H.TargetCylinder.centerX, (char *)&nwHdr.H.TargetCylinder.centerX,
					sizeof(nwHdr.H.TargetCylinder.centerX));
				LbCvSwap((char *)&fhdr.H.TargetCylinder.centerY, (char *)&nwHdr.H.TargetCylinder.centerY,
					sizeof(nwHdr.H.TargetCylinder.centerY));
				
				LbCvSwap((char *)&fhdr.H.CriticalZone.radius, (char *)&nwHdr.H.CriticalZone.radius,
					sizeof(nwHdr.H.CriticalZone.radius));
				LbCvSwap((char *)&fhdr.H.CriticalZone.zMin, (char *)&nwHdr.H.CriticalZone.zMin,
					sizeof(nwHdr.H.CriticalZone.zMin));
				LbCvSwap((char *)&fhdr.H.CriticalZone.zMax, (char *)&nwHdr.H.CriticalZone.zMax,
					sizeof(nwHdr.H.CriticalZone.zMax));
				LbCvSwap((char *)&fhdr.H.CriticalZone.centerX, (char *)&nwHdr.H.CriticalZone.centerX,
					sizeof(nwHdr.H.CriticalZone.centerX));
				LbCvSwap((char *)&fhdr.H.CriticalZone.centerY, (char *)&nwHdr.H.CriticalZone.centerY,
					sizeof(nwHdr.H.CriticalZone.centerY));
							
				LbCvSwap((char *)&fhdr.H.ObjectCylinder.radius, (char *)&nwHdr.H.ObjectCylinder.radius,
					sizeof(nwHdr.H.ObjectCylinder.radius));
				LbCvSwap((char *)&fhdr.H.ObjectCylinder.zMin, (char *)&nwHdr.H.ObjectCylinder.zMin,
					sizeof(nwHdr.H.ObjectCylinder.zMin));
				LbCvSwap((char *)&fhdr.H.ObjectCylinder.zMax, (char *)&nwHdr.H.ObjectCylinder.zMax,
					sizeof(nwHdr.H.ObjectCylinder.zMax));
				LbCvSwap((char *)&fhdr.H.ObjectCylinder.centerX, (char *)&nwHdr.H.ObjectCylinder.centerX,
					sizeof(nwHdr.H.ObjectCylinder.centerX));
				LbCvSwap((char *)&fhdr.H.ObjectCylinder.centerY, (char *)&nwHdr.H.ObjectCylinder.centerY,
					sizeof(nwHdr.H.ObjectCylinder.centerY));
			
				LbCvSwap((char *)&fhdr.H.LimitCylinder.radius, (char *)&nwHdr.H.LimitCylinder.radius,
					sizeof(nwHdr.H.LimitCylinder.radius));
				LbCvSwap((char *)&fhdr.H.LimitCylinder.zMin, (char *)&nwHdr.H.LimitCylinder.zMin,
					sizeof(nwHdr.H.LimitCylinder.zMin));
				LbCvSwap((char *)&fhdr.H.LimitCylinder.zMax, (char *)&nwHdr.H.LimitCylinder.zMax,
					sizeof(nwHdr.H.LimitCylinder.zMax));
				LbCvSwap((char *)&fhdr.H.LimitCylinder.centerX, (char *)&nwHdr.H.LimitCylinder.centerX,
					sizeof(nwHdr.H.LimitCylinder.centerX));
				LbCvSwap((char *)&fhdr.H.LimitCylinder.centerY, (char *)&nwHdr.H.LimitCylinder.centerY,
					sizeof(nwHdr.H.LimitCylinder.centerY));
			
				LbCvSwap((char *)&fhdr.H.NumSimulations, (char *)&nwHdr.H.NumSimulations,
					sizeof(nwHdr.H.NumSimulations));
			
				LbCvSwap((char *)&fhdr.H.SumEventsToSimulate, (char *)&nwHdr.H.SumEventsToSimulate,
					sizeof(nwHdr.H.SumEventsToSimulate));
			
				LbCvSwap((char *)&fhdr.H.NumPhotons, (char *)&nwHdr.H.NumPhotons,
					sizeof(nwHdr.H.NumPhotons));
			
				LbCvSwap((char *)&fhdr.H.NumDecays, (char *)&nwHdr.H.NumDecays,
					sizeof(nwHdr.H.NumDecays));
			}

	
			/* Write the new header to the output file */
			if (PhgHdrUpHeader(finalFile, &nwHdr, &headerHk) != true) {
				goto FAIL;
			}
	
			/* Close the final file */
			fclose(finalFile);
		}
		
		okay = true;
		FAIL:;
	} while (false);
	
	/* If error set due to cancellation, handle it, otherwise pass it on */
	if (!okay) {
		if (canceled) {
			ErHandle("User canceled phg.swap.", false);
			okay = true;
		}
	}
	
	/* Return the status */
	return (okay);
}
