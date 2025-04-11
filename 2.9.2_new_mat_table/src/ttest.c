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
*     Module Name: 			ttest.c
*     Revision Number:		1.5
*     Date last revised:	23 July 2013
*     Programmer:			Steven Vannoy
*     Date Originated:		Wednesday: September 9: 1992
*
*     Module Overview:	This is the main module for the ttest utility.
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*
*     Global variables defined:
*
*	  Global macros defined:
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

#include <stdio.h>
#include <signal.h>

#include "SystemDependent.h"

#include "LbTypes.h"

#include "LbError.h"
#include "LbDebug.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"

/* Local Constants */
#define	BUFF_SIZE	512				/* Number of elements read at a time */

/* Local types */
typedef struct {
	LbUsFourByte	countsA;
	LbUsFourByte	countsB;
	float			weightsA;
	float			weightsB;
	float			weightsqA;
	float			weightsqB;
} validBinTy;

	
/* Local Globals */
static	char		errString[1024];	/* For building error messages */
static	validBinTy	*validBins;

/* Prototypes */
Boolean	ttest(int argc, char *argv[]);
Boolean	getTTestFiles(char inputNames[6][1024], FILE *inputFiles[6], Boolean *cancelPtr);
Boolean getNumValidBins(LbUsFourByte hdrSize, LbUsFourByte numDataBins,
			LbUsFourByte elemSize, LbUsFourByte minValidCounts,
			FILE *inputFiles[6], LbUsFourByte *numValidPairsPtr, Boolean *cancelPtr);

Boolean doTTest(LbUsFourByte hdrSize, LbUsFourByte numDataBins,
			LbUsFourByte *numValidPairs, LbUsFourByte elemSize,
			LbUsFourByte minValidCounts, 
			FILE *inputFiles[6], double *tTests, LbUsFourByte *tTestBins, double *minTValue,
			double *maxTValue);



/**********************
*	ttest
*
*	Purpose:	Execute the program.
*
*	Result:	None.
***********************/

Boolean	ttest(int argc, char *argv[])
{
	int				dummy;					/* Removes compiler warning */
	char**			dummyPtr;				/* Removes compiler warning */
	Boolean			okay = false;			/* Process flag */
	Boolean			canceled = false;		/* Cancelation flag */
	char			prompt[1024];			/* Prompt for histogram loop */
	char			inputNames[6][1024];	/* The input file names */
	FILE			*inputFiles[6];			/* The input file descriptors */
	LbFourByte		fileSize = 0;			/* Size of the file */
	LbUsFourByte	elemSize = 4;			/* Size of a single data element */
	LbUsFourByte	hdrSize = 0;			/* Size of the input file header */
	LbUsFourByte	numDataBins = 0;		/* Number of data bins */
	LbUsFourByte	numValidPairs = 0;		/* Number of valid bins */
	LbUsFourByte	minValidCounts = 30;	/* Minimum counts for bin to be valid */
	LbUsFourByte	index;					/* Loop Control */
	LbUsFourByte	numHistBins;			/* Number of bins for histogram */
	LbFourByte		*histoGram = 0;			/* Histogram of T-Test Results */
	LbFourByte		possHistIndex;			/* Potential index for setting histogram value */
	LbUsFourByte	histIndex;				/* Non-negative index for setting histogram value */
	LbUsFourByte	sum;					/* Cumulative counter */
	LbUsFourByte	histLessThanMin=0;		/* Count of less than minimum */
	LbUsFourByte	histGreaterThanMax=0;	/* Count of values greater than maximum */
	LbUsFourByte	*tTestBins=0;			/* Bin numbers of valid bins */
	double			tTestRange;				/* Range of T-Test values */
	double			minTValue=FLT_MAX;		/* Minimum T-Test value */
	double			maxTValue=FLT_MIN;		/* Maximum T-Test value */
	double			*tTests = 0;			/* Storage for the T-Tests */
	double			rangeMin;				/* Minimum value for range of t-test values */
	double			rangeMax;				/* Maximum value for range of t-test values */
	double			midBin;					/* Midpoint of bin for printing out results */
	double			binSize;				/* Size of histogram bins */
	
	
	/* Avoid unused parameter compiler warnings */
	dummy = argc;
	dummyPtr = argv;
	
	do	 { /* Process Loop */
	
		/* Open the input files */
		if (getTTestFiles(inputNames, inputFiles, &canceled) == false) {
			ErHandle("Unable to get input files for T-Test\n", false);
			break;
		}
		
		/* See if they want To skip a header */
		hdrSize = LbInAskFourByte("Enter size of header to skip",
			1, false, false, 1, &canceled, 32768, 0, 0, 0);
		
		if (canceled){
			ErHandle("User canceled T-Test\n", false);
			goto CANCEL;
		}

		/* Get the number of bins, assume four byte data */
		{
			/* Seek to end of file */
			if (fseek(inputFiles[0], 0, SEEK_END) != 0) {
				sprintf(errString, "\nUnable to seek to end of file '%s'.\n",
					inputNames[0]);
				ErStFileError(errString);
				break;
			}
			
			/* Get the file size */
			if ((fileSize = ftell(inputFiles[0])) < 0){
				sprintf(errString, "\nUnable to get file position for file '%s'.\n",
					inputNames[0]);
				ErStFileError(errString);
				break;
			}
			
			/* Compute the number of data bins */
			numDataBins = (fileSize - hdrSize)/elemSize;
			
			/* Do a quick error check */
			if (((fileSize - hdrSize) % elemSize) != 0) {
				sprintf(errString, "\nFile does not contain an integral number of %ld byte elements, file '%s'.\n",
					(unsigned long)elemSize, inputNames[0]);
				ErStFileError(errString);
				break;
				
			}

			/* Seek back to beginning of file */
			if (fseek(inputFiles[0], 0, SEEK_SET) != 0) {
				sprintf(errString, "\nUnable to seek to beginning of file '%s'.\n",
					inputNames[0]);
				ErStFileError(errString);
				break;
			}
			
		}
		
		/* Get preliminary statistics */
		if (getNumValidBins(hdrSize, numDataBins, elemSize, minValidCounts,
				inputFiles, &numValidPairs, &canceled) == false) {
			break;
		}
		
		/* See if user decided to bag it */
		if (canceled){
			ErAlert("User canceled T-Test\n", false);
			goto CANCEL;
		}
		
		/* Allocate memory for T-Test results */
		if ((tTests = (double *)LbMmAlloc(numValidPairs*sizeof(double))) == 0) {
			break;
		}
		
		/* Allocate memory for T-Test bins */
		if ((tTestBins = (LbUsFourByte *)LbMmAlloc(numValidPairs*sizeof(LbUsFourByte))) == 0) {
			break;
		}
		
		/* Allocate memory for valid T-Test bins */
		if ((validBins = (validBinTy *)LbMmAlloc(numValidPairs*sizeof(validBinTy))) == 0) {
			break;
		}
		
		/* Perform the test */
		if (doTTest(hdrSize, numDataBins, &numValidPairs, elemSize, minValidCounts,
				inputFiles, tTests, tTestBins, &minTValue, &maxTValue) == false) {
			break;
		}
		
		/* See if they want to print T-Tests to the screen */
		if (LbInAskYesNo("\nDo you want to print T-Tests to stdout?", LBINYes) == LBINYes) {
		
			#define PRINT_EXTENSIVE
			#ifdef PRINT_EXTENSIVE
			{
			double sum, difference, squareRoot, tTestValue;
			
			LbInPrintf("\nBin #\tT-Test\tCts->A\tCts->B\tWts->A\tWts->B\tWts Sq->A\tWts Sq->B");
			for (index = 0; index < numValidPairs; index++) {
				sum = validBins[index].weightsqA+validBins[index].weightsqB;
				difference = validBins[index].weightsA-validBins[index].weightsB;
				squareRoot = sqrt(sum);
				tTestValue = difference/squareRoot;
				LbInPrintf("\n%5d  %2.3f   %2.3f ", tTestBins[index], tTests[index], tTestValue);
			}
			}
			#else
			LbInPrintf("\nBin #\tT-Test");
			for (index = 0; index < numValidPairs; index++) {
				LbInPrintf("\n%5d  %2.3f", tTestBins[index], tTests[index]);
			}
			#endif
			
			LbInPrintf("\n");
		}
		
		/* See if they want to histogram the results */
		sprintf(prompt,"\nDo you want to histogram the results?");
		
		while (LbInAskYesNo(prompt, LBINYes) == LBINYes) {
		
			/* See how many bins they want */
			numHistBins = LbInAskFourByte("Enter number of bins for histogram",
				0, false, false, 0, &canceled, 0, 0, 0, 0);
			
			if (canceled){
				ErHandle("User canceled T-Test\n", false);
				goto CANCEL;
			}
			
			/* Allocate memory for histogram */
			if ((histoGram = (LbFourByte *)LbMmAlloc(numHistBins*sizeof(LbFourByte))) == 0) {
				break;
			}
		
			/* Get the range limits */
			rangeMin = LbInAskFloat("Enter range minimum",
				1, true, 3, false, 1, &canceled, minTValue, 0.0, 0.0, 0.0);

			rangeMax = LbInAskFloat("Enter range maximum",
				1, true, 3, false, 1, &canceled, maxTValue, 0.0, 0.0, 0.0);
			
			/* Compute range */
			tTestRange = rangeMax - rangeMin;
			
			/* Verify you have a valid range */
			if (tTestRange <= 0.0) {
				ErStGeneric("Invalid range for tests, must be > 0");
				break;
			}
			
			/* Create the histogram */
			for (index = 0; index < numValidPairs; index++) {
			
				/* Compute index */
				possHistIndex = (LbFourByte)
					(((tTests[index] - rangeMin) * numHistBins)/tTestRange);
				
				/* Increment count */
				if (possHistIndex >= 0) {
					histIndex = possHistIndex;
					if (histIndex < numHistBins)
						histoGram[histIndex]++;
					else if (histIndex == numHistBins)
						histoGram[histIndex-1]++;
					else if (histIndex > numHistBins)
						histGreaterThanMax++;
				}
				else {
					histLessThanMin++;
				}
			}
			
			/* Compute the histogram bin size */
			binSize = tTestRange/numHistBins;
			
			/* Compute the mid-point of the first bin */
			midBin = binSize/2;
			
			/* Print the histogram out */
			LbInPrintf("\nHistogram of T-Test values over range [%3.2f, %3.2f]\n", rangeMin, rangeMax);
			LbInPrintf("\n%3.2f\t%d (Less than minimum)",rangeMin, histLessThanMin);
			for (index = 0; index < numHistBins; index++) {
				LbInPrintf("\n%3.2f\t%d",(rangeMin+(index*binSize)+midBin), histoGram[index]);
			}
			LbInPrintf("\n%3.2f\t%d (Greater than maximum)",rangeMax, histGreaterThanMax);
			LbInPrintf("\n");

			LbInPrintf("\nCumulative Histogram\n");
			sum = histoGram[0];
			LbInPrintf("\n%3.2f",rangeMin);
			for (index = 1; index < numHistBins; index++) {
				
				sum += histoGram[index];
				
				LbInPrintf("\n%3.2f\t%d",(rangeMin+(index*binSize)+midBin), sum);
			}
			LbInPrintf("\n%3.2f",rangeMax);
			LbInPrintf("\n");
			
			/* Free memory */
			if (histoGram != 0)
				LbMmFree((void **) &histoGram);
				
			sprintf(prompt,"\nDo you want to histogram the results again?");

		}
		
		CANCEL:;
		okay = true;
	} while (false);

	/* Close the input files */
	{
		if (inputFiles[0] != 0) {
			fclose(inputFiles[0]);
			inputFiles[0] = 0;
		}
	
		if (inputFiles[1] != 0) {
			fclose(inputFiles[1]);
			inputFiles[1] = 0;
		}
	
		if (inputFiles[2] != 0) {
			fclose(inputFiles[2]);
			inputFiles[2] = 0;
		}
	
		if (inputFiles[3] != 0) {
			fclose(inputFiles[3]);
			inputFiles[3] = 0;
		}
	
		if (inputFiles[4] != 0) {
			fclose(inputFiles[4]);
			inputFiles[4] = 0;
		}
	
		if (inputFiles[5] != 0) {
			fclose(inputFiles[5]);
			inputFiles[5] = 0;
		}
	}
	
	/* Free Memory */
	{
		if (tTests != 0)
			LbMmFree((void **) &tTests);
			
		if (tTestBins != 0)
			LbMmFree((void **) &tTestBins);
			
		if (histoGram != 0)
			LbMmFree((void **) &histoGram);
			
		if (validBins != 0)
			LbMmFree((void **) &validBins);
			
	}
	
	return(okay);
}

/**********************
*	getTTestFiles
*
*	Purpose: Get the input file names and open
*	the files.
*
*	Arguments:
*		char	**inputNames	- Storage for the file names.
*		FILE	**inputFiles	- Storage for the files.
*		Boolean	 *cancelPtr		- Cancelation flag.
*
*	Result:	True unless an error occurs.
***********************/

Boolean getTTestFiles(char inputNames[6][1024], FILE *inputFiles[6], Boolean *cancelPtr)
{
	Boolean	okay = false;	/* Process Flag */
	
	do { /* Process Loop */
	
		/* Initialize all names to zero */
		inputNames[0][0] = '\0';
		inputNames[1][0] = '\0';
		inputNames[2][0] = '\0';
		inputNames[3][0] = '\0';
		inputNames[4][0] = '\0';
		inputNames[5][0] = '\0';
		
		/* Initialize all files to null */
		inputFiles[0] = 0;
		inputFiles[1] = 0;
		inputFiles[2] = 0;
		inputFiles[3] = 0;
		inputFiles[4] = 0;
		inputFiles[5] = 0;
		
		/* Clear the cancelation flag */
		*cancelPtr = false;

		/* Get simulation 'A' count file */
		do {

			/* Get name of simulation 'A' count file */
			LbInAsk("Enter name of count file for simulation 'A'", 0, false,
					cancelPtr, 0, 0, 0, 0,
					inputNames[0]);
		
			/* Bolt if we *cancelPtr */
			if (*cancelPtr) {
				ErStCancel("User canceled.");
				goto CANCEL;
			}
		
			/* Open the file */
			if ((inputFiles[0] = LbFlFileOpen(inputNames[0], "rb")) == 0) {
				sprintf(errString, "Unable to open file '%s'\n", inputNames[0]);
				ErAlert(errString, false);
			}
		} while ((*cancelPtr == false) && (inputFiles[0] == 0));
		
		/* Check for cancelation */
		if (*cancelPtr == true)
			goto CANCEL;
			
		/* Get simulation 'A' weight file */
		do {

			/* Get name of simulation 'A' weight file */
			LbInAsk("Enter name weight file for simulation 'A'", 0, false,
					cancelPtr, 0, 0, 0, 0,
					inputNames[1]);
		
			/* Bolt if we *cancelPtr */
			if (*cancelPtr) {
				ErStCancel("User canceled.");
				goto CANCEL;
			}
		
		
			/* Open the file */
			if ((inputFiles[1] = LbFlFileOpen(inputNames[1], "rb")) == 0) {
				sprintf(errString, "Unable to open file '%s'\n", inputNames[1]);
				ErAlert(errString, false);
			}

		} while ((*cancelPtr == false) && (inputFiles[1] == 0));
		
		/* Check for cancelation */
		if (*cancelPtr == true)
			goto CANCEL;

		/* Get simulation 'A' weight squared file */
		do {

			/* Get name of simulation 'A' weight squared file */
			LbInAsk("Enter name of weight squared file for simulation 'A'", 0, false,
					cancelPtr, 0, 0, 0, 0,
					inputNames[2]);
		
			/* Bolt if we *cancelPtr */
			if (*cancelPtr) {
				ErStCancel("User canceled.");
				goto CANCEL;
			}
		
		
			/* Open the file */
			if ((inputFiles[2] = LbFlFileOpen(inputNames[2], "rb")) == 0) {
				sprintf(errString, "Unable to open file '%s'\n", inputNames[2]);
				ErAlert(errString, false);
			}

		} while ((*cancelPtr == false) && (inputFiles[2] == 0));
		
		/* Check for cancelation */
		if (*cancelPtr == true)
			goto CANCEL;

		/* Get simulation 'B' count file */
		do {

			/* Get name of simulation 'B' count file */
			LbInAsk("Enter name of count file for simulation 'B'", 0, false,
					cancelPtr, 0, 0, 0, 0,
					inputNames[3]);
		
			/* Bolt if we *cancelPtr */
			if (*cancelPtr) {
				ErStCancel("User canceled.");
				goto CANCEL;
			}
					
		
		
			/* Open the file */
			if ((inputFiles[3] = LbFlFileOpen(inputNames[3], "rb")) == 0) {
				sprintf(errString, "Unable to open file '%s'\n", inputNames[3]);
				ErAlert(errString, false);
			}

		} while ((*cancelPtr == false) && (inputFiles[3] == 0));
		
		/* Check for cancelation */
		if (*cancelPtr == true)
			goto CANCEL;

		/* Get simulation 'B' weight file */
		do {

			/* Get name of simulation 'B' weight file */
			LbInAsk("Enter name of weight file for simulation 'B'", 0, false,
					cancelPtr, 0, 0, 0, 0,
					inputNames[4]);
		
			/* Bolt if we *cancelPtr */
			if (*cancelPtr) {
				ErStCancel("User canceled.");
				goto CANCEL;
			}

		
			/* Open the file */
			if ((inputFiles[4] = LbFlFileOpen(inputNames[4], "rb")) == 0) {
				sprintf(errString, "Unable to open file '%s'\n", inputNames[4]);
				ErAlert(errString, false);
			}

		} while ((*cancelPtr == false) && (inputFiles[4] == 0));
		
		/* Check for cancelation */
		if (*cancelPtr == true)
			goto CANCEL;
			
		/* Get simulation 'B' weight squared file */
		do {

			/* Get name of simulation 'B' weight squared file */
			LbInAsk("Enter name of weight squared file for simulation 'B'", 0, false,
					cancelPtr, 0, 0, 0, 0,
					inputNames[5]);
		
			/* Bolt if we *cancelPtr */
			if (*cancelPtr) {
				ErStCancel("User canceled.");
				goto CANCEL;
			}
			
			
			/* Open the file */
			if ((inputFiles[5] = LbFlFileOpen(inputNames[5], "rb")) == 0) {
				sprintf(errString, "Unable to open file '%s'\n", inputNames[5]);
				ErAlert(errString, false);
			}

		} while ((*cancelPtr == false) && (inputFiles[5] == 0));
		
		CANCEL:;
		okay = true;
	} while (false);
	
	/* If we failed, close the inputFiles that were opened */
	if (okay == false) {
	
		if (inputFiles[0] != 0) {
			fclose(inputFiles[0]);
			inputFiles[0] = 0;
		}
	
		if (inputFiles[1] != 0) {
			fclose(inputFiles[1]);
			inputFiles[1] = 0;
		}
	
		if (inputFiles[2] != 0) {
			fclose(inputFiles[2]);
			inputFiles[2] = 0;
		}
	
		if (inputFiles[3] != 0) {
			fclose(inputFiles[3]);
			inputFiles[3] = 0;
		}
	
		if (inputFiles[4] != 0) {
			fclose(inputFiles[4]);
			inputFiles[4] = 0;
		}
	
		if (inputFiles[5] != 0) {
			fclose(inputFiles[5]);
			inputFiles[5] = 0;
		}
	}
			
	return (okay);
}

/**********************
*	getNumValidBins
*
*	Purpose: Get number of "count" bins which contain >= minValidCounts
*			 for both files
*	the files.
*
*	Arguments:
*		LbUsFourByte	hdrSize				- Size of the files header.
*		LbUsFourByte	numDataBins			- The number of bins.
*		LbUsFourByte	elemSize			- The sizeo of the elements.
*		LbUsFourByte	minValidCounts		- Minimum threshold for validity.
*		FILE			**inputFiles		- Storage for the files.
*		LbUsFourByte	*numValidPairsPtr	- Number of valid bins.
*		Boolean			*cancelPtr			- Cancelation flag.
*
*	Result:	True unless an error occurs.
***********************/

Boolean getNumValidBins(LbUsFourByte hdrSize, LbUsFourByte numDataBins,
			LbUsFourByte elemSize, LbUsFourByte minValidCounts,
			FILE *inputFiles[6], LbUsFourByte *numValidPairsPtr, Boolean *cancelPtr)
{
	Boolean			okay = false;			/* Process flag */
	LbUsFourByte	bufferA[BUFF_SIZE];		/* Buffer for input data */
	LbUsFourByte	bufferB[BUFF_SIZE];		/* Buffer for input data */
	LbUsFourByte	numBuffs = 0;			/* Number of big buffers */
	LbUsFourByte	elemsInPartialBuff = 0;	/* Number of partial buffer elements */
	LbUsFourByte	buffIndex;				/* Index for looping */
	LbUsFourByte	elemIndex;				/* Index for looping */
	LbUsFourByte	numMixedPairs=0;			/* Number where 1 bin is good and the other is not */
	LbUsFourByte	numBadPairs=0;			/* Number where both bins are too small, but non-zero */
	LbUsFourByte	numZeroPairs=0;			/* Number where both bins are zero */
	LbUsFourByte	numValidBinsA=0;		/* Number of valid bins in file 'a' */
	LbUsFourByte	numValidBinsB=0;		/* Number of valid bins in file 'b' */
	
	
	do { /* Process loop */
	
		/* Clear counter */
		*numValidPairsPtr = 0;
		
		/* Set up buffer counts */
		numBuffs = numDataBins/BUFF_SIZE;
		elemsInPartialBuff = numDataBins % BUFF_SIZE;
		
		/* Seek past headers */
		if (hdrSize != 0) {

			/* Seek past the header */
			if (fseek(inputFiles[0], hdrSize, SEEK_SET) != 0) {
				ErStFileError("\nUnable to seek past header.");
				break;
			}

			/* Seek past the header */
			if (fseek(inputFiles[3], hdrSize, SEEK_SET) != 0) {
				ErStFileError("\nUnable to seek past header.");
				break;
			}
		}
		
		/* Loop through big buffers */
		for (buffIndex = 0; buffIndex < numBuffs; buffIndex++) {
		
			/* Read a buffer */
			if (fread(bufferA, elemSize, BUFF_SIZE, inputFiles[0]) != BUFF_SIZE) {
				ErStFileError("\nUnable to read buffer from count file 'A'\n");
				goto FAIL;
			}
		
			/* Read a buffer */
			if (fread(bufferB, elemSize, BUFF_SIZE, inputFiles[3]) != BUFF_SIZE) {
				ErStFileError("\nUnable to read buffer from count file 'B'\n");
				goto FAIL;
			}
			
			/* Loop through getting statistics */
			for (elemIndex = 0; elemIndex < BUFF_SIZE; elemIndex++){
				
				if ((bufferA[elemIndex] == 0) && (bufferB[elemIndex] == 0)) {
					numZeroPairs++;
				}
				else if ((bufferA[elemIndex] >= minValidCounts) && (bufferB[elemIndex] >= minValidCounts)){
					(*numValidPairsPtr)++;
				}
				else if ((bufferA[elemIndex] >= minValidCounts) && (bufferB[elemIndex] < minValidCounts)) {
					numMixedPairs++;
				}
				else if ((bufferA[elemIndex] < minValidCounts) && (bufferB[elemIndex] >= minValidCounts)) {
					numMixedPairs++;
				}
				else if ((bufferA[elemIndex] < minValidCounts) && (bufferB[elemIndex] < minValidCounts)) {
					numBadPairs++;
				}
				if (bufferA[elemIndex] >= minValidCounts){
					(numValidBinsA)++;
				}
				if (bufferB[elemIndex] >= minValidCounts){
					(numValidBinsB)++;
				}
			}
		}
		
		/* Read partial buffers */
		{
			/* Read a buffer */
			if (fread(bufferA, elemSize, elemsInPartialBuff, inputFiles[0]) != elemsInPartialBuff) {
				ErStFileError("\nUnable to read buffer from count file 'A'\n");
				goto FAIL;
			}
		
			/* Read a buffer */
			if (fread(bufferB, elemSize, elemsInPartialBuff, inputFiles[3]) != elemsInPartialBuff) {
				ErStFileError("\nUnable to read buffer from count file 'B'\n");
				goto FAIL;
			}
		}
			
		/* Loop through getting statistics */
		for (elemIndex = 0; elemIndex < elemsInPartialBuff; elemIndex++){
			
			if ((bufferA[elemIndex] == 0) && (bufferB[elemIndex] == 0)) {
				numZeroPairs++;
			}
			else if ((bufferA[elemIndex] >= minValidCounts) && (bufferB[elemIndex] >= minValidCounts)){
				(*numValidPairsPtr)++;
			}
			else if ((bufferA[elemIndex] >= minValidCounts) && (bufferB[elemIndex] < minValidCounts)) {
				numMixedPairs++;
			}
			else if ((bufferA[elemIndex] < minValidCounts) && (bufferB[elemIndex] >= minValidCounts)) {
				numMixedPairs++;
			}
			else if ((bufferA[elemIndex] < minValidCounts) && (bufferB[elemIndex] < minValidCounts)) {
				numBadPairs++;
			}
			if (bufferA[elemIndex] >= minValidCounts){
				(numValidBinsA)++;
			}
			if (bufferB[elemIndex] >= minValidCounts){
				(numValidBinsB)++;
			}
		}
		
		
		/* Report the statistics */
		LbInPrintf("\nNumber of bins checked = %ld\n", (unsigned long)numDataBins);
		LbInPrintf("Number of valid bins in set 'A' = %ld\n", (unsigned long)numValidBinsA);
		LbInPrintf("Number of valid bins in set 'B' = %ld\n", (unsigned long)numValidBinsB);
		LbInPrintf("Number of valid pairs found = %ld\n", (unsigned long)(*numValidPairsPtr));
		LbInPrintf("Number of pairs with one valid and one too small = %ld\n", (unsigned long)numMixedPairs);
		LbInPrintf("Number of pairs with both > 0 but < min = %ld\n", (unsigned long)numBadPairs);
		LbInPrintf("Number of pairs with both zero = %ld\n\n", (unsigned long)numZeroPairs);
		
		/* See if they want to continue */
		if (LbInAskYesNo("\nDo you want to continue", LBINYes) != LBINYes) {
			*cancelPtr = true;
		}
		else
			*cancelPtr = false;
				
		okay = true;
		FAIL:;
	} while (false);
	
	return (okay);
}


/**********************
*	doTTest
*
*	Purpose: Perform the T-Test
*	the files.
*
*	Arguments:
*		LbUsFourByte	hdrSize				- Size of the files header.
*		LbUsFourByte	numDataBins			- The number of bins.
*		LbUsFourByte	*numValidPairs		- Number of valid bins.
*		LbUsFourByte	elemSize			- The sizeo of the elements.
*		LbUsFourByte	minValidCounts		- Minimum threshold for validity.
*		FILE			**inputFiles		- Storage for the files.
*		double			*tTests				- The T-Test results.
*		LbUsFourByte	*tTestBins			- The bin number of the value
*		double			*minTValue			- Minimum T-Test value.
*		double			*maxTValue			- Maximum T-Test value.
*	Result:	True unless an error occurs.
***********************/

Boolean doTTest(LbUsFourByte hdrSize, LbUsFourByte numDataBins,
			LbUsFourByte *numValidPairs, LbUsFourByte elemSize,
			LbUsFourByte minValidCounts, 
			FILE *inputFiles[6], double *tTests, LbUsFourByte *tTestBins, double *minTValue,
			double *maxTValue)
{
	Boolean			okay = false;			/* Process flag */
	LbUsFourByte	numBuffs = 0;			/* Number of big buffers */
	LbUsFourByte	elemsInPartialBuff = 0;	/* Number of partial buffer elements */
	LbUsFourByte	buffIndex;				/* Index for looping */
	LbUsFourByte	elemIndex;				/* Index for looping */
	LbUsFourByte	tIndex=0;				/* Index for T-Tests */
	LbUsFourByte	num0StdDev=0;			/* Number of T-Tests less than one standard deviation different */
	LbUsFourByte	num1StdDev=0;			/* Number of T-Tests one standard deviation different */
	LbUsFourByte	num2StdDev=0;			/* Number of T-Tests one standard deviation different */
	LbUsFourByte	num3StdDev=0;			/* Number of T-Tests one standard deviation different or more */
	LbUsFourByte	binNumber=0;			/* The bin number for the valid test element */
	double			percentDistribution;	/* For reporting where the values fall */
	double			weightDifference;		/* Difference of two weights */
	double			weightSquareSum;		/* Sum of weight squares */
	LbUsFourByte	countsA[BUFF_SIZE];		/* Buffer for input data */
	LbUsFourByte	countsB[BUFF_SIZE];		/* Buffer for input data */
	float			weightsA[BUFF_SIZE];	/* Buffer for input data */
	float			weightsB[BUFF_SIZE];	/* Buffer for input data */
	float			weightsSquA[BUFF_SIZE];	/* Buffer for input data */
	float			weightsSquB[BUFF_SIZE];	/* Buffer for input data */


	do { /* Process loop */
		
		/* Clear count of valid pairs, although this was computed via previous
			routines based on the count image, it will be updated here
			according to the weight images. This prevents any mistakes
			between what is in the count image and what the weights
			actualy are.
		*/
		*numValidPairs = 0;
		
		/* Set up buffer counts */
		numBuffs = numDataBins/BUFF_SIZE;
		elemsInPartialBuff = numDataBins % BUFF_SIZE;
		
		/* Seek to beginning of data */
		{
			/* Seek past the header */
			if (fseek(inputFiles[0], hdrSize, SEEK_SET) != 0) {
				ErStFileError("\nUnable to seek past header.");
				break;
			}
			/* Seek past the header */
			if (fseek(inputFiles[1], hdrSize, SEEK_SET) != 0) {
				ErStFileError("\nUnable to seek past header.");
				break;
			}
			/* Seek past the header */
			if (fseek(inputFiles[2], hdrSize, SEEK_SET) != 0) {
				ErStFileError("\nUnable to seek past header.");
				break;
			}

			/* Seek past the header */
			if (fseek(inputFiles[3], hdrSize, SEEK_SET) != 0) {
				ErStFileError("\nUnable to seek past header.");
				break;
			}

			/* Seek past the header */
			if (fseek(inputFiles[4], hdrSize, SEEK_SET) != 0) {
				ErStFileError("\nUnable to seek past header.");
				break;
			}
			
			/* Seek past the header */
			if (fseek(inputFiles[5], hdrSize, SEEK_SET) != 0) {
				ErStFileError("\nUnable to seek past header.");
				break;
			}
		}
		
		/* Loop through big buffers */
		for (buffIndex = 0; buffIndex < numBuffs; buffIndex++) {
		
			/* Read a buffer */
			if (fread(countsA, elemSize, BUFF_SIZE, inputFiles[0]) != BUFF_SIZE) {
				ErStFileError("\nUnable to read buffer from count file 'A'\n");
				goto FAIL;
			}
		
			/* Read a buffer */
			if (fread(weightsA, elemSize, BUFF_SIZE, inputFiles[1]) != BUFF_SIZE) {
				ErStFileError("\nUnable to read buffer from weight file 'A'\n");
				goto FAIL;
			}
		
			/* Read a buffer */
			if (fread(weightsSquA, elemSize, BUFF_SIZE, inputFiles[2]) != BUFF_SIZE) {
				ErStFileError("\nUnable to read buffer from weight squared file 'A'\n");
				goto FAIL;
			}
		
			/* Read a buffer */
			if (fread(countsB, elemSize, BUFF_SIZE, inputFiles[3]) != BUFF_SIZE) {
				ErStFileError("\nUnable to read buffer from count file 'B'\n");
				goto FAIL;
			}
		
			/* Read a buffer */
			if (fread(weightsB, elemSize, BUFF_SIZE, inputFiles[4]) != BUFF_SIZE) {
				ErStFileError("\nUnable to read buffer from weight file 'B'\n");
				goto FAIL;
			}
		
			/* Read a buffer */
			if (fread(weightsSquB, elemSize, BUFF_SIZE, inputFiles[5]) != BUFF_SIZE) {
				ErStFileError("\nUnable to read buffer from weight squared file 'B'\n");
				goto FAIL;
			}
			
			/* Loop through performing test */
			for (elemIndex = 0; elemIndex < BUFF_SIZE; elemIndex++){
				
				/* See if this buffer meets criteria */
				if ((countsA[elemIndex] >= minValidCounts) && (countsB[elemIndex] >= minValidCounts)){
				
					/* Sum weight squared here for efficiency */
					weightSquareSum = ((double) weightsSquA[elemIndex]) + ((double) weightsSquB[elemIndex]);
					
					/* Verify weights are valid */
					if ((weightsA[elemIndex] > 0.0) && (weightsB[elemIndex] > 0.0) &&
							(weightSquareSum > 0.0)) {
					
						/* Perform the T-Test */
						weightDifference = ((double) weightsA[elemIndex]) - ((double) weightsB[elemIndex]);
						tTests[tIndex] = weightDifference/sqrt(weightSquareSum);
						tTestBins[tIndex] = binNumber;
						validBins[tIndex].countsA = countsA[elemIndex];
						validBins[tIndex].countsB = countsB[elemIndex];
						validBins[tIndex].weightsA = weightsA[elemIndex];
						validBins[tIndex].weightsB = weightsB[elemIndex];
						validBins[tIndex].weightsqA = weightsSquA[elemIndex];
						validBins[tIndex].weightsqB = weightsSquB[elemIndex];
						
						/* Check for max value */
						if (tTests[tIndex] > *maxTValue)
							*maxTValue = tTests[tIndex];
						
						/* Check for minimum value */
						if (tTests[tIndex] < *minTValue)
							*minTValue = tTests[tIndex];
	
						/* Check for standard deviation ranges */
						if ((tTests[tIndex] < -3.0) || (tTests[tIndex] > 3.0))
							num3StdDev++;
						else if ((tTests[tIndex] < -2.0) || (tTests[tIndex] > 2.0))
							num2StdDev++;
						else if ((tTests[tIndex] < -1.0) || (tTests[tIndex] > 1.0))
							num1StdDev++;
						else
							num0StdDev++;
							
						/* Increment tIndex */
						tIndex++;
						
						/* Increment count of valid bins */
						(*numValidPairs)++;
					}
				}
				binNumber++;
			}
		}
		
		/* Read partial buffers */
		{
			/* Read a buffer */
			if (fread(countsA, elemSize, elemsInPartialBuff, inputFiles[0]) != elemsInPartialBuff) {
				ErStFileError("\nUnable to read buffer from count file 'A'\n");
				goto FAIL;
			}
		
			/* Read a buffer */
			if (fread(weightsA, elemSize, elemsInPartialBuff, inputFiles[1]) != elemsInPartialBuff) {
				ErStFileError("\nUnable to read buffer from weight file 'A'\n");
				goto FAIL;
			}
		
			/* Read a buffer */
			if (fread(weightsSquA, elemSize, elemsInPartialBuff, inputFiles[2]) != elemsInPartialBuff) {
				ErStFileError("\nUnable to read buffer from weight squared file 'A'\n");
				goto FAIL;
			}
		
			/* Read a buffer */
			if (fread(countsB, elemSize, elemsInPartialBuff, inputFiles[3]) != elemsInPartialBuff) {
				ErStFileError("\nUnable to read buffer from count file 'B'\n");
				goto FAIL;
			}
		
			/* Read a buffer */
			if (fread(weightsB, elemSize, elemsInPartialBuff, inputFiles[4]) != elemsInPartialBuff) {
				ErStFileError("\nUnable to read buffer from weight file 'B'\n");
				goto FAIL;
			}
		
			/* Read a buffer */
			if (fread(weightsSquB, elemSize, elemsInPartialBuff, inputFiles[5]) != elemsInPartialBuff) {
				ErStFileError("\nUnable to read buffer from weight squared file 'B'\n");
				goto FAIL;
			}
		}
			
		/* Loop through performing test */
		for (elemIndex = 0; elemIndex < elemsInPartialBuff; elemIndex++){
			
			/* See if this buffer meets initial criteria */
			if ((countsA[elemIndex] >= minValidCounts) && (countsB[elemIndex] >= minValidCounts)){
					
				/* Sum weight squared here for efficiency */
				weightSquareSum = ((double) weightsSquA[elemIndex]) + ((double) weightsSquB[elemIndex]);
				
				/* Verify weights are valid */
				if ((weightsA[elemIndex] > 0.0) && (weightsB[elemIndex] > 0.0) &&
						(weightSquareSum > 0.0)) {
				
					/* Perform the T-Test */
					weightDifference = ((double) weightsA[elemIndex]) - ((double) weightsB[elemIndex]);
					tTests[tIndex] = weightDifference/sqrt(weightSquareSum);
					tTestBins[tIndex] = binNumber;
					validBins[tIndex].countsA = countsA[elemIndex];
					validBins[tIndex].countsB = countsB[elemIndex];
					validBins[tIndex].weightsA = weightsA[elemIndex];
					validBins[tIndex].weightsB = weightsB[elemIndex];
					validBins[tIndex].weightsqA = weightsSquA[elemIndex];
					validBins[tIndex].weightsqB = weightsSquB[elemIndex];
					
				
					/* Check for max value */
					if (tTests[tIndex] > *maxTValue)
						*maxTValue = tTests[tIndex];
					
					/* Check for minimum value */
					if (tTests[tIndex] < *minTValue)
						*minTValue = tTests[tIndex];

					/* Check for standard deviation ranges */
					if ((tTests[tIndex] < -3.0) || (tTests[tIndex] > 3.0))
						num3StdDev++;
					else if ((tTests[tIndex] < -2.0) || (tTests[tIndex] > 2.0))
						num2StdDev++;
					else if ((tTests[tIndex] < -1.0) || (tTests[tIndex] > 1.0))
						num1StdDev++;
					else
						num0StdDev++;
	
					/* Increment tIndex */
					tIndex++;
					
					/* Increment count of valid bins */
					(*numValidPairs)++;
				}
			}
			
			binNumber++;
		}
		
		
		/* Report the statistics */
		LbInPrintf("Num pairs with valid weights = %d\n", *numValidPairs);
		LbInPrintf("Minimum T-Test value = %2.3f\n", *minTValue);
		LbInPrintf("Maximum T-Test value = %2.3f\n\n", *maxTValue);
		
		percentDistribution = (((double)num0StdDev/(double)(*numValidPairs)) * 100.0);
		LbInPrintf("0.0 <= |T-Test| <= 1.0  \t= %d => %%%3.2f\n",
			num0StdDev, percentDistribution);
			
		percentDistribution += ((double)((double)num1StdDev/(double)(*numValidPairs)) * 100.0);
		LbInPrintf("0.0 <= |T-Test| <= 2.0  \t= %d => %%%3.2f\n",
			num0StdDev+num1StdDev, percentDistribution);
			
		percentDistribution += ((double)((double)num2StdDev/(double)(*numValidPairs)) * 100.0);
		LbInPrintf("0.0 <= |T-Test| <= 3.0  \t= %d => %%%3.2f\n",
			num0StdDev+num1StdDev+num2StdDev, percentDistribution);
		
		/* If there are values beyond 3 standard deviations, report how many */
		if ((*minTValue < -3.0) || (*maxTValue > 3.0)) {
			percentDistribution += ((double)((double)num3StdDev/(double)(*numValidPairs)) * 100.0);
			LbInPrintf("0.0 <= |T-Test| <= %3.2f \t= %d, => %%%3.2f\n",
				((fabs(*minTValue) < *maxTValue) ? *maxTValue : fabs(*minTValue)), num0StdDev+num1StdDev+num2StdDev+num3StdDev, percentDistribution);
		}
			
				
		okay = true;
		FAIL:;
	} while (false);
	
	return (okay);
}

