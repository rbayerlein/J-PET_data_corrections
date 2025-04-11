/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1992-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		ProdTbl.c
*			Revision Number:	1.3
*			Date last revised:	23 July 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	9 September 1992
*
*			Module Overview:	Productivity Table routines.
*
*			References:			'Productivity Table Processes' PHG design.
*
**********************************************************************************
*
*			Global functions defined:
*				ProdTblAddDetectedProductivity
*				ProdTblAddStartingProductivity
*				ProdTblCloseTable
*				ProdTblCreateTable
*				ProdTblGetZmaxZmin
*				ProdTblGetAngleIndex
*				ProdTblTerminate
*
*			Global macros defined:
*					PRODTBLGetProdTblAngleEnd
*					PRODTBLGetProdTblAngleSize
*					PRODTBLGetProdTblAngleStart
*					PRODTBLGetPrimCellProductivity
*					PRODTBLGetScatCellProductivity
*					PRODTBLGetMaxCellProductivity
*					PRODTBLGetNumAngleCells
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
*			Programmer(s):		Robert Harrison
*
*			Revision date:		Mar 2010
*			Revision description:	Changed variables of ProdTblGetAngleIndex
*							to match those in calling function.
*
*********************************************************************************/

#define	PRODUCTIVITY_TABLE

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"

#include "LbMacros.h"
#include "LbError.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbFile.h"
#include "LbInterface.h"
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ProdTbl.h"
#include "SubObj.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "PhgMath.h"
#include "ColUsr.h"
#include "CylPos.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "PhoHFile.h"
#include "phg.h"
#include "PhgBin.h"

/* LOCAL CONSTANTS */
#define PRODTBL_PRODTBL_FILEPREFIX	"phgprod."						/* Prefix to starting productivity table */

/* LOCAL TYPES */
typedef	double						prodTblSquareWeightTy;
typedef	prodTblSquareWeightTy		prodTblSquareWeightArrayTy[PRODTBL_TOT_STRAT_CELLS];
typedef prodTblSquareWeightArrayTy	*prodTblSquareWeightTableTy;

typedef	double						prodTblWeightTy;
typedef	prodTblSquareWeightTy		prodTblWeightArrayTy[PRODTBL_TOT_STRAT_CELLS];
typedef prodTblSquareWeightArrayTy	*prodTblWeightTableTy;

/* LOCAL GLOBALS */
Boolean						ProdTblIsInitialized = false;
prodTblSquareWeightTableTy	ProdTblStartPrimPhoSquWeights = 0;		/* Table of sum of started, primary photons's weights */
prodTblSquareWeightTableTy	ProdTblStartScatPhoSquWeights = 0;		/* Table of sum of started, scattered photons's weights */
prodTblSquareWeightTableTy	ProdTblDetPrimPhoSquWeights = 0;		/* Table of sum of detected, primary photons's weights */
prodTblSquareWeightTableTy	ProdTblDetScatPhoSquWeights = 0;		/* Table of sum of detected, scattered photons's weights */

prodTblSquareWeightTableTy	ProdTblStartPrimPhoWeights = 0;			/* Table of sum of started, primary photons's weights squared */
prodTblSquareWeightTableTy	ProdTblStartScatPhoWeights = 0;			/* Table of sum of started, scattered photons's weights squared */
prodTblSquareWeightTableTy	ProdTblDetPrimPhoWeights = 0;			/* Table of sum of detected, primary photons's weights squared */
prodTblSquareWeightTableTy	ProdTblDetScatPhoWeights = 0;			/* Table of sum of detected, scattered photons's weights squared */

ProdTblProdTblTy			ProdTblCalculatedPrimProdTbl = 0;		/* Our calculated productivity table for primary photons */
ProdTblProdTblTy			ProdTblCalculatedScatProdTbl = 0;		/* Our calculated productivity table for scattered photons */
ProdTblProdTblTy			ProdTblDetPerBinTbl = 0;				/* Counts of detections per stratification bin */

/* LOCAL MACROS */
#define PRODTBLFgIsScatter(flag)	LbFgIsSet((flag), PRODTBLFg_Scatter)
#define PRODTBLFgIsPrimary(flag)	LbFgIsSet((flag), PRODTBLFg_Primary)

/* PROTOTYPES */

/* FUNCTIONS */
void	prodTblParseProdLine(ProdTblProdTblTy prodTable, LbTwoByte sliceIndex, char *prodLine);
void	prodTblFreeTables(void);

/*********************************************************************************
*
*			Name:		ProdTblAddStartingProductivity
*
*			Summary:	Record the weight of the tracking_photon in the productivity table 
*						as a photon starting from the cell 
*						containing the tracking_photon's position.
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*				LbUsFourByte		whichOne			- Specifies primary or scatter.
*			Function return: None.
*
*********************************************************************************/
void ProdTblAddStartingProductivity(PHG_TrackingPhoton *trackingPhotonPtr,
		LbUsFourByte whichOne)	
{	
	
	#ifdef PHG_DEBUG
		if (ProdTblIsInitialized == false) {
			PhgAbort("You can't call ProdTblAddStartingProductivity without"
				" initializing the ProdTbl module.", false);
		}
	#endif
	
   	/* Add the photon's contribution (Primary vs Scatter based on "whichOne") */
	if (PRODTBLFgIsPrimary(whichOne)) {
		
		/* Sum the squared weights */
		ProdTblStartPrimPhoSquWeights[trackingPhotonPtr->sliceIndex][trackingPhotonPtr->angleIndex] +=
			PHGMATH_Square(trackingPhotonPtr->decay_weight * trackingPhotonPtr->photon_primary_weight);
			
		/* Sum the weights */
		ProdTblStartPrimPhoWeights[trackingPhotonPtr->sliceIndex][trackingPhotonPtr->angleIndex] +=
			trackingPhotonPtr->decay_weight * trackingPhotonPtr->photon_primary_weight;

	}
	else {
		/* Sum the squared weights */
		ProdTblStartScatPhoSquWeights[trackingPhotonPtr->sliceIndex][trackingPhotonPtr->angleIndex] +=
			PHGMATH_Square(trackingPhotonPtr->decay_weight * trackingPhotonPtr->photon_scatter_weight);

		/* Sum the weights */
		ProdTblStartScatPhoWeights[trackingPhotonPtr->sliceIndex][trackingPhotonPtr->angleIndex] +=
			trackingPhotonPtr->decay_weight * trackingPhotonPtr->photon_scatter_weight;
	}
		
}

/*********************************************************************************
*
*			Name:		ProdTblAddDetectedProductivity
*
*			Summary:	Record the weight of the tracking_photon in the productivity table 
*						as a photon detected from the cell 
*						containing the tracking_photon's position.
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*				LbUsFourByte		whichOne			- Which one flag.
*			Function return: None.
*
*********************************************************************************/
void ProdTblAddDetectedProductivity(PHG_TrackingPhoton *trackingPhotonPtr,
		LbUsFourByte whichOne)	
{
	#ifdef PHG_DEBUG
		if (ProdTblIsInitialized == false) {
			PhgAbort("You can't call ProdTblAddDetectedProductivity without"
				" initializing the ProdTbl module.", false);
		}
	#endif
	
	
	/* Increment counts per bin */
	ProdTblDetPerBinTbl[trackingPhotonPtr->sliceIndex].productivity[trackingPhotonPtr->angleIndex].cellProductivity += 1;

	/* Add the photon's contribution */
	if (PRODTBLFgIsPrimary(whichOne)) {

		/* Sum primary squared weight */
		ProdTblDetPrimPhoSquWeights[trackingPhotonPtr->sliceIndex][trackingPhotonPtr->angleIndex] +=
			PHGMATH_Square(trackingPhotonPtr->decay_weight * trackingPhotonPtr->photon_primary_weight);
		
		/* Sum primary weight */
		ProdTblDetPrimPhoWeights[trackingPhotonPtr->sliceIndex][trackingPhotonPtr->angleIndex] +=
			trackingPhotonPtr->decay_weight * trackingPhotonPtr->photon_primary_weight;
	}
	else {
	
		/* Sum scattered squared weight */
		ProdTblDetScatPhoSquWeights[trackingPhotonPtr->sliceIndex][trackingPhotonPtr->angleIndex] +=
			PHGMATH_Square(trackingPhotonPtr->decay_weight * trackingPhotonPtr->photon_scatter_weight);
	
		/* Sum scattered  weight */
		ProdTblDetScatPhoWeights[trackingPhotonPtr->sliceIndex][trackingPhotonPtr->angleIndex] +=
			trackingPhotonPtr->decay_weight * trackingPhotonPtr->photon_scatter_weight;
	}
		
}

/*********************************************************************************
*
*			Name:		ProdTblCloseTable
*
*			Summary:	Complete the productivity table and write it to an output file.
*			Arguments:
*				ProdTblProdTblInfoTy 	*prodTableInfoPtr		- Information
*
*			Function return: true unless an error occurs.
*
*********************************************************************************/
Boolean ProdTblCloseTable(ProdTblProdTblInfoTy *prodTableInfoPtr)	
{
	Boolean			okay = false;	/* Process Flag */
	double			avgPrimProd;	/* Average productivity for primary photons */
	double			avgScatProd;	/* Average productivity for scatter photons */
	double			minAccPrimProd;	/* Minimum acceptable primary productivity (computed) */
	double			minAccScatProd;	/* Minimum acceptable scatter productivity (computed) */
	char			errString[256];	/* Storage to create an error string */
	FILE			*prodFile=0;	/* Output productivity file */
	LbUsFourByte	sliceIndex;		/* LCV for slices */
	LbUsFourByte	angleIndex;		/* LCV for angles */
	
	do {	/* Process Loop */

		#ifdef PHG_DEBUG
			if (ProdTblIsInitialized == false) {
				PhgAbort("You can't call ProdTblCloseTable without"
					" initializing the ProdTbl module.", false);
			}
		#endif
		
	
		/* Compute productivities if requested */
		if (strlen(prodTableInfoPtr->outputFileName) != 0) {
		
			/* Calculate the actual productivities */
			avgPrimProd = 0;
			avgScatProd = 0;
			
			for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
				
				/* Set slice info */
				ProdTblCalculatedPrimProdTbl[sliceIndex].zMin = ProdTblPrimProdTbl[sliceIndex].zMin;
				ProdTblCalculatedPrimProdTbl[sliceIndex].zMax = ProdTblPrimProdTbl[sliceIndex].zMax;
	
				for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
					
					/* Compute primary productivity */
					if ((ProdTblDetPrimPhoSquWeights[sliceIndex][angleIndex] == 0) ||
						    (ProdTblStartPrimPhoSquWeights[sliceIndex][angleIndex] == 0)){
						ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity = 0;
					}
					else {
						ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity =
							PHGMATH_SquareRoot(ProdTblDetPrimPhoSquWeights[sliceIndex][angleIndex] /
							ProdTblStartPrimPhoSquWeights[sliceIndex][angleIndex]);
							
						/* Update running sum for computation of average */
						avgPrimProd += ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity;
					}
					ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary =
						ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary;
						
					ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary =
						ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary;
	
					/* Compute scatter productivity */
					if ((ProdTblDetScatPhoSquWeights[sliceIndex][angleIndex] == 0) ||
						    (ProdTblStartScatPhoSquWeights[sliceIndex][angleIndex] == 0)){
						ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].cellProductivity = 0;
					}
					else {
						ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].cellProductivity =
							PHGMATH_SquareRoot(ProdTblDetScatPhoSquWeights[sliceIndex][angleIndex] /
							ProdTblStartScatPhoSquWeights[sliceIndex][angleIndex]);

						/* Update running sum for computation of average */
						avgScatProd += ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].cellProductivity;
					}
					ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary =
						ProdTblScatProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary;
						
					ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary =
						ProdTblScatProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary;
	
				}
			}
			
			/* Compute averages */
			avgPrimProd = avgPrimProd/(ProdTblNumSlices * ProdTblNumAngleCells);
			avgScatProd = avgScatProd/(ProdTblNumSlices * ProdTblNumAngleCells);
			
			/* Compute minimums */
			minAccPrimProd = avgPrimProd/10;
			minAccScatProd = avgScatProd/10; 
			
			/* Now assure that no productivities are less than computed min */
			for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
				
				for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
					
					/* Compute primary productivity accounting for minimum threshold */
					if (ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity < minAccPrimProd)
						ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity = minAccPrimProd;

					/* Compute scatter productivity accounting for minimum threshold */
					if (ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].cellProductivity < minAccScatProd)
						ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].cellProductivity = minAccScatProd;
				}
			}
			
			#ifdef PHG_DEBUG
				if (PhgDebugDumpFile != 0) {
					for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
						/* Print the blank line */
						fprintf(PhgDebugDumpFile, "\n\n");
			
						/* Print out the slices productivities */
						fprintf(PhgDebugDumpFile, "Slice %ld: Computed Slice productivities (Primary), [productivity startOfBoundary endOfBoundary]\n", (unsigned long)sliceIndex);
						for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
							fprintf(PhgDebugDumpFile, "[%3.2e %3.2e %3.2e] ",
								ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity,
								ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary,
								ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary);
						}
			
						/* Print the blank line */
						fprintf(PhgDebugDumpFile, "\n\n");
			
						/* Print out the slices productivities */
						fprintf(PhgDebugDumpFile, "Slice %ld: Computed Slice productivities (Scatter), [productivity startOfBoundary endOfBoundary]\n", (unsigned long)sliceIndex);
						for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
							fprintf(PhgDebugDumpFile, "[%3.2e %3.2e %3.2e] ",
								ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].cellProductivity,
								ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary,
								ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary);
						}
					}
					
					/* Print a blank line */
					fprintf(PhgDebugDumpFile, "\n\n");
					
					/* Dump the counts per bin */
					for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
						/* Print the blank line */
						fprintf(PhgDebugDumpFile, "\n");
			
						/* Print out the slices productivities */
						fprintf(PhgDebugDumpFile, "Slice %ld: Counts per bin (Primary), [counts startOfBoundary endOfBoundary]\n", (unsigned long)sliceIndex);
						for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
							fprintf(PhgDebugDumpFile, "%ld ",
								(unsigned long) ProdTblDetPerBinTbl[sliceIndex].productivity[angleIndex].cellProductivity);
						}
			
						/* Print the blank line */
						fprintf(PhgDebugDumpFile, "\n\n");
					}
					/* Print the blank line */
					fprintf(PhgDebugDumpFile, "\n\n");
				}
			#endif
			
			/* Open the file */
			if ((prodFile = LbFlFileOpen(prodTableInfoPtr->outputFileName, "w")) == 0) {
				sprintf(errString, "Unable to open file productivity output file named '%s'. ", 
						prodTableInfoPtr->outputFileName);
				ErStFileError(errString);
				break;
			}
			
			/* Write out the acceptance angle */
			fprintf(prodFile, "Acceptance angle = %f\n", prodTableInfoPtr->acceptanceAngle);
			fprintf(prodFile, "Number of slices = %ld\n", (unsigned long)(prodTableInfoPtr->numSlices));
			fprintf(prodFile, "Primary Productivities\n");
			
			/* Print out the productivities */
			for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
				
				fprintf(prodFile, "zMin = %f zMax = %f\n",
					ProdTblPrimProdTbl[sliceIndex].zMin, ProdTblPrimProdTbl[sliceIndex].zMax);
				
				for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
					fprintf(prodFile, "%1.6f %1.6f %1.6f ",
						ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity,
						ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary,
						ProdTblCalculatedPrimProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary);
				}
			   fprintf(prodFile, "\n");
			}
			
			fprintf(prodFile, "\nScatter Productivities\n");
			
			/* Print out the productivities */
			for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
				
				fprintf(prodFile, "zMin = %f zMax = %f\n",
					ProdTblPrimProdTbl[sliceIndex].zMin, ProdTblPrimProdTbl[sliceIndex].zMax);
				
				for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
					fprintf(prodFile, "%1.6f %1.6f %1.6f ",
						ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].cellProductivity,
						ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary,
						ProdTblCalculatedScatProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary);
				}
			   fprintf(prodFile, "\n");
			}
	
		}
		
		okay = true;
	} while (false);
	
	/* Free the square weight tables */
	prodTblFreeTables();

	/* Close the productivity file if it was opened */
	if (prodFile != 0)
		fclose(prodFile);
		
	return (okay);
}

/*********************************************************************************
*
*			Name:		prodTblFreeTables
*
*			Summary:	Free all tables allocated by this module. Note, each
*						table is checked for allocation so this routine can be
*						called as a part of error handling cleanup.
*			Arguments:
*
*			Function return: None.
*
*********************************************************************************/
void prodTblFreeTables()
{
	if (ProdTblStartPrimPhoSquWeights != 0) {
		LbMmFree((void **) &ProdTblStartPrimPhoSquWeights);
	}
	if (ProdTblStartScatPhoSquWeights != 0) {
		LbMmFree((void **) &ProdTblStartScatPhoSquWeights);
	}
	if (ProdTblDetPrimPhoSquWeights != 0) {
		LbMmFree((void **) &ProdTblDetPrimPhoSquWeights);
	}
	if (ProdTblDetScatPhoSquWeights != 0) {
		LbMmFree((void **) &ProdTblDetScatPhoSquWeights);
	}
	if (ProdTblStartPrimPhoWeights != 0) {
		LbMmFree((void **) &ProdTblStartPrimPhoWeights);
	}
	if (ProdTblStartScatPhoWeights != 0) {
		LbMmFree((void **) &ProdTblStartScatPhoWeights);
	}
	if (ProdTblDetPrimPhoWeights != 0) {
		LbMmFree((void **) &ProdTblDetPrimPhoWeights);
	}
	if (ProdTblDetScatPhoWeights != 0) {
		LbMmFree((void **) &ProdTblDetScatPhoWeights);
	}
	if (ProdTblPrimProdTbl != 0) {
		LbMmFree((void **) &ProdTblPrimProdTbl);
	}
	if (ProdTblScatProdTbl != 0) {
		LbMmFree((void **) &ProdTblScatProdTbl);
	}
	if (ProdTblMaxProdTbl != 0) {
		LbMmFree((void **) &ProdTblMaxProdTbl);
	}
	if (ProdTblCalculatedPrimProdTbl != 0) {
		LbMmFree((void **) &ProdTblCalculatedPrimProdTbl);
	}
	if (ProdTblCalculatedScatProdTbl != 0) {
		LbMmFree((void **) &ProdTblCalculatedScatProdTbl);
	}
	if (ProdTblDetPerBinTbl != 0) {
		LbMmFree((void **) &ProdTblDetPerBinTbl);
	}
}
/*********************************************************************************
*
*			Name:		ProdTblCreateTable
*
*			Summary:	Initialize a new productivity table based on the number of slices 
*						in the object and the axial-angle stratification descriptions.
*			Arguments:
*				ProdTblProdTblInfoTy	*prodTableInfoPtr - Info to create the table.
*
*			Function return: true unless an error occurs.
*
*********************************************************************************/
Boolean ProdTblCreateTable(ProdTblProdTblInfoTy *prodTableInfoPtr)	
{
	Boolean			okay = false;			/* Process Flag */
	char			errString[1028];		/* Space for building error strings */
	char			inputBuffer[5140];		/* Must hold 24 float strings */
	double			primaryProd;			/* Productivity for a primary cell */
	double			scatterProd;			/* Productivity for a scattered cell */
	double			currentBinStart;		/* Current value for start of bins */
	double			productiveRange;		/* Range of productive stratification angles */
	double			unProductiveRange;		/* Range of unproductive stratification angles */
	double			productiveBinSize;		/* Size of stratification bin for productive angle bins */
	double			unProductiveBinSize;	/* Size of stratification bin for unproductive angle bins */
	FILE			*prodFile = 0;			/* Productivity description file */
	double			acceptanceAngle;		/* Acceptance angle stored in productivity file */
	LbUsFourByte	sliceIndex;				/* Index into current slice */
	LbUsFourByte	angleIndex;				/* Current angle */
	LbUsFourByte	numSlices;				/* Test var for number of slices */
	unsigned long	numSlicesUSL;
	
	do { /* Process Loop */

		#ifdef PHG_DEBUG
			if (ProdTblIsInitialized == false) {
				PhgAbort("You can't call ProdTblCreateTable without"
					" initializing the ProdTbl module.", false);
			}
		#endif	
		
		/* Save the number of slices */
		ProdTblNumSlices = prodTableInfoPtr->numSlices;

		/* Allocate the primary productivity table */
		if ((ProdTblPrimProdTbl = (ProdTblProdTblTy) LbMmAlloc(
				sizeof(ProdTblProdTblElemTy) * prodTableInfoPtr->numSlices)) == 0) {
			break;
		}
		
		/* Allocate the scatter productivity table */
		if ((ProdTblScatProdTbl = (ProdTblProdTblTy) LbMmAlloc(
				sizeof(ProdTblProdTblElemTy) * prodTableInfoPtr->numSlices)) == 0) {
			break;
		}
		
		/* Allocate the maximum productivity table */
		if ((ProdTblMaxProdTbl = (ProdTblProdTblTy) LbMmAlloc(
				sizeof(ProdTblProdTblElemTy) * prodTableInfoPtr->numSlices)) == 0) {
			break;
		}

		/* Allocate the calculated primary productivity table */
		if ((ProdTblCalculatedPrimProdTbl = (ProdTblProdTblTy) LbMmAlloc(
				sizeof(ProdTblProdTblElemTy) * prodTableInfoPtr->numSlices)) == 0) {
			break;
		}
		
		/* Allocate the calculated scatter productivity table */
		if ((ProdTblCalculatedScatProdTbl = (ProdTblProdTblTy) LbMmAlloc(
				sizeof(ProdTblProdTblElemTy) * prodTableInfoPtr->numSlices)) == 0) {
			break;
		}
		
		/* Allocate the counts per bin productivity table */
		if ((ProdTblDetPerBinTbl = (ProdTblProdTblTy) LbMmAlloc(
				sizeof(ProdTblProdTblElemTy) * prodTableInfoPtr->numSlices)) == 0) {
			break;
		}
		
		/* If file name supplied then do it from the file */
		if (strlen(prodTableInfoPtr->inputFileName) != 0) {
			/* Open the file */
			if ((prodFile = LbFlFileOpen(prodTableInfoPtr->inputFileName, "r")) == 0) {
				sprintf(errString, "Unable to open productivity input file '%s'. ",
						prodTableInfoPtr->inputFileName);
				ErStFileError(errString);
				goto FAIL;
			}
	
			/* Get the acceptance angle */
			LbFlFGetS(inputBuffer, sizeof(inputBuffer), prodFile);
			if (sscanf(inputBuffer, "Acceptance angle = %lf", &acceptanceAngle) != 1){
					sprintf(errString, "Error from sscanf on acceptance angle in (ProdTblCreateTable).\n"
					"Input string is '%s'\n"
					"Productivity Input file name specified is [%s]", inputBuffer, prodTableInfoPtr->inputFileName);
					ErStGeneric(errString);
					goto FAIL;
			}
	
			/* Verify here that the acceptance angle is the same as user specified value */
			if (!PhgMathRealNumAreEqual(acceptanceAngle, prodTableInfoPtr->acceptanceAngle, -7, 0, 0, 0)) {
 					sprintf(errString, "Conflicting acceptance angles in productivity input (%3.2e) versus user parameters (%3.2e).",
						acceptanceAngle, prodTableInfoPtr->acceptanceAngle);
						ErStGeneric(errString);
					goto FAIL;
			}
			
			/* Get number of slices */
			LbFlFGetS(inputBuffer, sizeof(inputBuffer), prodFile);
			if (sscanf(inputBuffer, "Number of slices = %ld", &numSlicesUSL) != 1){
					sprintf(errString, "Error from sscanf on number of slices in (ProdTblCreateTable).\n"
					"Input string is '%s'\n"
					"Productivity Input file name specified is [%s]", inputBuffer, prodTableInfoPtr->inputFileName);
					ErStGeneric(errString);
					goto FAIL;
			}
			numSlices = (LbUsFourByte)numSlicesUSL;
			
			/* Verify that slice numbers match */
			if (numSlices != prodTableInfoPtr->numSlices) {
				ErStGeneric("Number of slices in object does not match productivity table!");
				goto FAIL;
			}
	
			/* Get the Primary Productivities seperator */
			LbFlFGetS(inputBuffer, sizeof(inputBuffer), prodFile);
			if (strcmp(inputBuffer, "Primary Productivities\n") != 0){
					ErStGeneric("Primary Productivities seperator not found in productivity file.\n");
					goto FAIL;
			}
			
			/* Read in productivity values for each slice */
			for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
	
				/* Get zMin/zMax for slice */
				LbFlFGetS(inputBuffer, sizeof(inputBuffer), prodFile);
				if (sscanf(inputBuffer, "zMin = %lf zMax = %lf",
						&(ProdTblPrimProdTbl[sliceIndex].zMin),
						&(ProdTblPrimProdTbl[sliceIndex].zMax)) != 2){
					
					sprintf(errString, "Error from sscanf on zMin/zMax in (ProdTblCreateTable).\n"
					"Input string is '%s'\n"
					"Productivity Input file name specified is [%s]", inputBuffer, prodTableInfoPtr->inputFileName);
					ErStGeneric(errString);
					goto FAIL;
				}

				/* Get the slices productivities */
				LbFlFGetS(inputBuffer, sizeof(inputBuffer), prodFile);
	
				/* Parse the line */
				prodTblParseProdLine(ProdTblPrimProdTbl, sliceIndex, inputBuffer);
			}
		}
		else if (PHG_IsStratification() == false){
			/* Stratification is off so build table accordingly */
			for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
				ProdTblPrimProdTbl[sliceIndex].zMin = SUBOBJGetSliceMinZ(sliceIndex);
				ProdTblPrimProdTbl[sliceIndex].zMax = SUBOBJGetSliceMaxZ(sliceIndex);
				ProdTblPrimProdTbl[sliceIndex].productivity[0].cellProductivity = 1;
				ProdTblPrimProdTbl[sliceIndex].productivity[0].startOfBoundary = -1;
				ProdTblPrimProdTbl[sliceIndex].productivity[0].endOfBoundary = 1;

				ProdTblScatProdTbl[sliceIndex].zMin = 
					ProdTblPrimProdTbl[sliceIndex].zMin;
				ProdTblScatProdTbl[sliceIndex].zMax = 
					ProdTblPrimProdTbl[sliceIndex].zMax;
				ProdTblScatProdTbl[sliceIndex].productivity[0] = 
					ProdTblPrimProdTbl[sliceIndex].productivity[0];
			
			}
		}
		else {
			
			/* Calculate productive bin size as cosine of 90-acceptance angle (converted to radians) */
			productiveRange = 
				PHGMATH_Cosine((90-PhgRunTimeParams.Phg_AcceptanceAngle) * (PHGMATH_PI/180));
			
			/* Unproductive range is everything else */
			unProductiveRange = 1 - productiveRange;

			/* Bins sizes are a function of range and number of bins (from 0 to 90, hence /2) */
			productiveBinSize = (productiveRange/(ProdTblNumAccAngleCells/2));
			unProductiveBinSize = (unProductiveRange/(ProdTblNumUnAccAngleCells/2));
			
			
			/* Stratification is on, but we are building an even productivity table */
			for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {

				/* Set start of stratification cells to -1 */
				currentBinStart = -1;
			
				for (angleIndex = 0; angleIndex < (ProdTblNumUnAccAngleCells/2); angleIndex++) {
					
					ProdTblPrimProdTbl[sliceIndex].zMin = SUBOBJGetSliceMinZ(sliceIndex);
					ProdTblPrimProdTbl[sliceIndex].zMax = SUBOBJGetSliceMaxZ(sliceIndex);
					
					ProdTblScatProdTbl[sliceIndex].zMin = SUBOBJGetSliceMinZ(sliceIndex);
					ProdTblScatProdTbl[sliceIndex].zMax = SUBOBJGetSliceMaxZ(sliceIndex);
					
					ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity = 1;
					ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary = currentBinStart;
					ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary = 
						currentBinStart + unProductiveBinSize;
						
					
					ProdTblScatProdTbl[sliceIndex].productivity[angleIndex] =
						ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex];
						
					currentBinStart += unProductiveBinSize;
				}

				for (angleIndex = (ProdTblNumUnAccAngleCells/2); 
						angleIndex < ((ProdTblNumUnAccAngleCells/2) + ProdTblNumAccAngleCells); angleIndex++) {

					ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity = 1;
					ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary = currentBinStart;
					ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary = 
						currentBinStart + productiveBinSize;
					
					ProdTblScatProdTbl[sliceIndex].productivity[angleIndex] =
						ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex];
						
					currentBinStart += productiveBinSize;
				}

				for (angleIndex = ((ProdTblNumUnAccAngleCells/2) + ProdTblNumAccAngleCells);
						angleIndex < (ProdTblNumUnAccAngleCells + ProdTblNumAccAngleCells); 
						angleIndex++) {
						
					ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity = 1;
					ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary = currentBinStart;
					ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary = 
						currentBinStart + unProductiveBinSize;
					
					ProdTblScatProdTbl[sliceIndex].productivity[angleIndex] =
						ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex];
						
					currentBinStart += unProductiveBinSize;
				}
				
				
			}

		}


		/* If file name supplied then initialize stratification table for scatter
			photons from file .
			(if not, and stratification is on, the scatter table was created with the primary
			table above)
		*/
		if (strlen(prodTableInfoPtr->inputFileName) != 0) {
				
			/* Read the blank line */
			LbFlFGetS(inputBuffer, sizeof(inputBuffer), prodFile);
			
			/* Get the Scatter Productivities seperator */
			LbFlFGetS(inputBuffer, sizeof(inputBuffer), prodFile);
			if (strcmp(inputBuffer, "Scatter Productivities\n") != 0){
					ErStGeneric("Primary Productivities seperator not found in productivity file.\n");
					goto FAIL;
			}
			
			/* Read in productivity values for each slice */
			for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
	
				/* Get zMin/zMax for slice */
				LbFlFGetS(inputBuffer, sizeof(inputBuffer), prodFile);
				if (sscanf(inputBuffer, "zMin = %lf zMax = %3lf",
						&(ProdTblScatProdTbl[sliceIndex].zMin),
						&(ProdTblScatProdTbl[sliceIndex].zMax)) != 2){
					
					sprintf(errString, "Error from sscanf on zMin/zMax for scatter photons in (ProdTblCreateTable).\n"
					"Input string is '%s'\n"
					"Productivity Input file name specified is [%s]", inputBuffer, prodTableInfoPtr->inputFileName);
					ErStGeneric(errString);
					goto FAIL;
				}
				
				/* Get the slices productivities */
				LbFlFGetS(inputBuffer, sizeof(inputBuffer), prodFile);
	
				/* Parse the line */
				prodTblParseProdLine(ProdTblScatProdTbl, sliceIndex, inputBuffer);
			}
		}
		else  if (PHG_IsStratification() == false){
			/* Stratification is off so build table accordingly */
			for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
				ProdTblScatProdTbl[sliceIndex].zMin = SUBOBJGetSliceMinZ(sliceIndex);
				ProdTblScatProdTbl[sliceIndex].zMax = SUBOBJGetSliceMaxZ(sliceIndex);
				ProdTblScatProdTbl[sliceIndex].productivity[0].cellProductivity = 1;
				ProdTblScatProdTbl[sliceIndex].productivity[0].startOfBoundary = -1;
				ProdTblScatProdTbl[sliceIndex].productivity[0].endOfBoundary = 1;
			}
		}

		/* Now initialize maximum productivity table values */
		if (PHG_IsSPECT()) {
			
			/* Table is the same except for cell productivities */
			for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
				
				/* Set zMin/zMax */
				ProdTblMaxProdTbl[sliceIndex].zMin = 
					ProdTblPrimProdTbl[sliceIndex].zMin;
				ProdTblMaxProdTbl[sliceIndex].zMax = 
					ProdTblPrimProdTbl[sliceIndex].zMax;
					
				/* Now set stratification cells */
				for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
					
					ProdTblMaxProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary = 
						ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary;
					ProdTblMaxProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary = 
						ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary;
						
					/* Calculate max productivity */
					ProdTblMaxProdTbl[sliceIndex].productivity[angleIndex].cellProductivity = 
						PHGMATH_Max(ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity,
						ProdTblScatProdTbl[sliceIndex].productivity[angleIndex].cellProductivity);
				}
			}
		}
		else {
			for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
				
				/* Set zMin/zMax */
				ProdTblMaxProdTbl[sliceIndex].zMin = 
					ProdTblPrimProdTbl[sliceIndex].zMin;
				ProdTblMaxProdTbl[sliceIndex].zMax = 
					ProdTblPrimProdTbl[sliceIndex].zMax;
					
				/* Now set stratification cells */
				for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
					
					ProdTblMaxProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary = 
						ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary;
					ProdTblMaxProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary = 
						ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary;
						
					/* Calculate max productivity */
					primaryProd = ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity
						* ProdTblPrimProdTbl[sliceIndex].productivity[(ProdTblNumAngleCells-1) - angleIndex].cellProductivity;
					
					scatterProd = ProdTblScatProdTbl[sliceIndex].productivity[angleIndex].cellProductivity
						* ProdTblScatProdTbl[sliceIndex].productivity[(ProdTblNumAngleCells-1) - angleIndex].cellProductivity;

					ProdTblMaxProdTbl[sliceIndex].productivity[angleIndex].cellProductivity = 
						PHGMATH_Max(primaryProd, scatterProd);
				}
			}
		}

		/* Now allocate square weight tables */
		if ((ProdTblStartPrimPhoSquWeights = (prodTblSquareWeightTableTy) LbMmAlloc(
				sizeof(prodTblSquareWeightArrayTy) * ProdTblNumSlices)) == 0) {
				
			ErStGeneric("Unable to allocate memory for primary photon square weight table.");
			goto FAIL;
		}

		if ((ProdTblStartScatPhoSquWeights = (prodTblSquareWeightTableTy) LbMmAlloc(
				sizeof(prodTblSquareWeightArrayTy) * ProdTblNumSlices)) == 0) {
				
			ErStGeneric("Unable to allocate memory for scattered photon square weight table.");
			goto FAIL;
		}
		
		if ((ProdTblDetPrimPhoSquWeights = (prodTblSquareWeightTableTy) LbMmAlloc(
				sizeof(prodTblSquareWeightArrayTy) * ProdTblNumSlices)) == 0) {
				
			goto FAIL;
		}
		
		if ((ProdTblDetScatPhoSquWeights = (prodTblSquareWeightTableTy) LbMmAlloc(
				sizeof(prodTblSquareWeightArrayTy) * ProdTblNumSlices)) == 0) {
				
			goto FAIL;
		}

		/* Now allocate  weight tables */
		if ((ProdTblStartPrimPhoWeights = (prodTblWeightTableTy) LbMmAlloc(
				sizeof(prodTblWeightArrayTy) * ProdTblNumSlices)) == 0) {
				
			goto FAIL;
		}

		if ((ProdTblStartScatPhoWeights = (prodTblWeightTableTy) LbMmAlloc(
				sizeof(prodTblWeightArrayTy) * ProdTblNumSlices)) == 0) {
				
			goto FAIL;
		}
		
		if ((ProdTblDetPrimPhoWeights = (prodTblWeightTableTy) LbMmAlloc(
				sizeof(prodTblWeightArrayTy) * ProdTblNumSlices)) == 0) {
				
			goto FAIL;
		}
		
		if ((ProdTblDetScatPhoWeights = (prodTblWeightTableTy) LbMmAlloc(
				sizeof(prodTblWeightArrayTy) * ProdTblNumSlices)) == 0) {
				
			goto FAIL;
		}

		/* Now initialize the values */
		for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {

			for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
				ProdTblStartPrimPhoSquWeights[sliceIndex][angleIndex] = 0.0;
				ProdTblStartScatPhoSquWeights[sliceIndex][angleIndex] = 0.0;
				ProdTblDetScatPhoSquWeights[sliceIndex][angleIndex] = 0.0;
				ProdTblDetPrimPhoSquWeights[sliceIndex][angleIndex] = 0.0;

				ProdTblStartPrimPhoWeights[sliceIndex][angleIndex] = 0.0;
				ProdTblStartScatPhoWeights[sliceIndex][angleIndex] = 0.0;
				ProdTblDetScatPhoWeights[sliceIndex][angleIndex] = 0.0;
				ProdTblDetPrimPhoWeights[sliceIndex][angleIndex] = 0.0;
				
				/* This is a little kludgy, but do counts per bin initialization here */
				ProdTblDetPerBinTbl[sliceIndex].productivity[angleIndex].cellProductivity = 0;
			}
		}

		okay = true;
		FAIL:;
	} while (false);
	
	/* Do error handling */
	if (okay == false) {

		/* Free the primary productivity table */
		if (ProdTblPrimProdTbl !=  0) {
			LbMmFree((void **)&ProdTblPrimProdTbl);
		}
		
		/* Free the scatter productivity table */
		if (ProdTblScatProdTbl != 0) {
			LbMmFree((void **)&ProdTblScatProdTbl);
		}
		
		/* Free the maximum productivity table */
		if (ProdTblMaxProdTbl != 0) {
			LbMmFree((void **)&ProdTblMaxProdTbl);
		}

		/* Free the calculated primary productivity table */
		if (ProdTblCalculatedPrimProdTbl != 0) {
			LbMmFree((void **)&ProdTblCalculatedPrimProdTbl);
		}
		
		/* Free the calculated scatter productivity table */
		if (ProdTblCalculatedScatProdTbl != 0) {
			LbMmFree((void **)&ProdTblCalculatedScatProdTbl);
		}
		
		/* Free the counts per bin productivity table */
		if (ProdTblDetPerBinTbl != 0) {
			LbMmFree((void **)&ProdTblDetPerBinTbl);
		}
	}
	
	/* If the productivity file was opened, close it */
	if (prodFile != 0)
		fclose(prodFile);
		
	return (okay);
}

/*********************************************************************************
*
*			Name:		ProdTblInitialize
*
*			Summary:	Initialize the manager.
*			Arguments:
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean ProdTblInitialize()	
{
	do { /* Process Loop */
	
		/* If the library is initialized and we are being called, abort with an
			error message in debug mode. In non-debug mode just ignore the call
			and go on.
		*/
		if (ProdTblIsInitialized == true) {
			#ifdef PHG_DEBUG
			PhgAbort("You can't call ProdTblInitialize twice!", false);
			#endif
			break;
		}
		
		/* Initialize globals */
		if (PHG_IsStratification()){
			ProdTblNumAngleCells = PRODTBL_TOT_STRAT_CELLS;
			ProdTblNumAccAngleCells = PRODTBL_ACC_STRAT_CELLS;
			ProdTblNumUnAccAngleCells = PRODTBL_NOTACC_STRAT_CELLS;
		}
		else {
			ProdTblNumAngleCells = 1;
			ProdTblNumAccAngleCells = 0;
			ProdTblNumUnAccAngleCells = 0;
		}
		ProdTblIsInitialized = true;
	} while (false);
	
	return (ProdTblIsInitialized);
}

/*********************************************************************************
*
*			Name:		ProdTblGetAngleIndex
*
*			Summary:	Return productivity angle index for given z_cosine.
*			Arguments:
*				LbUsFourByte	sliceIndex		- The current slice index.
*				double			z_cosine;		- The the z cosine.
*			Function return: Angle index of z cosine.
*
*********************************************************************************/
LbFourByte ProdTblGetAngleIndex(LbFourByte sliceIndex, double z_cosine)	
{
	LbFourByte angleIndex;
	
	#ifdef PHG_DEBUG
		if (ProdTblIsInitialized == false) {
			PhgAbort("You can't call ProdTblGetAngleIndex without"
				" initializing the ProdTbl module.", false);
		}
	#endif
	
	/* Find the angle index */
	for (angleIndex = 0; (LbUsFourByte)angleIndex < ProdTblNumAngleCells; angleIndex++) {
		if ((z_cosine >= 
				ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary)
				&&
				(z_cosine <=
				ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary)){
			
			break;
		}
	}
	
	return (angleIndex);
}

/*********************************************************************************
*
*			Name:		ProdTblGetZmaxZmin
*
*			Summary:	Return z max and min for a given slice.
*			Arguments:
*				LbUsFourByte	sliceIndex	- The slice we are looking at.
*				double			*zMax;		- The max z value for this slice.
*				double			*zMin;		- The min z value for this slice.
*			Function return: None.
*
*********************************************************************************/
void ProdTblGetZmaxZmin(LbUsFourByte sliceIndex, double *zMax, double *zMin)	
{

	#ifdef PHG_DEBUG
		if (ProdTblIsInitialized == false) {
			PhgAbort("You can't call ProdTblGetZmaxZmin without"
				" initializing the ProdTbl module.", false);
		}
	#endif	
	
	/* Set the values */
	*zMax = ProdTblPrimProdTbl[sliceIndex].zMax;
	*zMin = ProdTblPrimProdTbl[sliceIndex].zMin;
}

/*********************************************************************************
*
*			Name:			prodTblParseProdLine
*
*			Summary:		Parse a line of productivities for the prod table.
*
*			Arguments:
*				ProdTblProdTblTy	prodTable	- The table we are initializing.
*				LbTwoByte			sliceIndex	- The current slice.
*				char				*prodLine	- The productivity description line.
*
*			Function return: None.
*
*********************************************************************************/
void prodTblParseProdLine(ProdTblProdTblTy prodTable, LbTwoByte sliceIndex, char *prodLine)
{
	char				prodString[24];		/* Individual productivity specification */
	LbTwoByte			angleIndex;			/* Index to current productivity cell */
	LbTwoByte			prodIndex;			/* Index into productivity string */
	LbTwoByte			inputIndex;			/* Index into input buffer */

	/* Parse the input string */
	inputIndex = 0;
	prodIndex = 0;
	angleIndex = 0;
	while (prodLine[inputIndex] != '\n') {
	
		/* Verify we have a valid string */
		if (prodLine[inputIndex] == '\0') {
			PhgAbort("Invalid string read for slice productivities (prodTblParseProdLine).",
				false);
		}
		
		/* Get cell productivity; Copy digits until space is found */
		while ((prodString[prodIndex] = prodLine[inputIndex]) != ' ') {
			
			/* Verify we have a valid string */
			if (prodString[prodIndex] == '\0') {
				PhgAbort("Invalid string for productivity (prodTblParseProdLine).",
					false);
			}
			
			prodIndex++;
			inputIndex++;
		}
		
		/* Index past seperator */
		inputIndex++;
		
		/* Terminate the productivity string */
		prodString[prodIndex] = '\0';

		/* Convert index to numeric value and save it */
		prodTable[sliceIndex].productivity[angleIndex].cellProductivity = atof(prodString);
		
		/* Setup for next index */
		prodIndex = 0;
		
		/* Get starting boundary; Copy digits until space is found */
		while ((prodString[prodIndex] = prodLine[inputIndex]) != ' ') {
			
			/* Verify we have a valid string */
			if (prodString[prodIndex] == '\0') {
				PhgAbort("Invalid string for productivity (prodTblParseProdLine).",
					false);
			}
			
			prodIndex++;
			inputIndex++;
		}
		
		/* Index past seperator */
		inputIndex++;
		
		/* Terminate the productivity string */
		prodString[prodIndex] = '\0';

		/* Convert index to numeric value and save it */
		prodTable[sliceIndex].productivity[angleIndex].startOfBoundary = atof(prodString);
		
		/* Setup for next index */
		prodIndex = 0;
		
		/* Get end of boundary; Copy digits until space is found */
		while ((prodString[prodIndex] = prodLine[inputIndex]) != ' ') {
			
			/* Verify we have a valid string */
			if (prodString[prodIndex] == '\0') {
				PhgAbort("Invalid string for productivity (prodTblParseProdLine).",
					false);
			}
			
			prodIndex++;
			inputIndex++;
		}
		
		/* Index past seperator */
		inputIndex++;
		
		/* Terminate the productivity string */
		prodString[prodIndex] = '\0';

		/* Convert index to numeric value and save it */
		prodTable[sliceIndex].productivity[angleIndex].endOfBoundary = atof(prodString);
		
		/* Setup for next index */
		prodIndex = 0;
		angleIndex++;
	}
	
	/* Set number of angle cells to angle index */
	if (angleIndex != 0)
		ProdTblNumAngleCells = angleIndex;
}

/*********************************************************************************
*
*			Name:		ProdTblTerminate
*
*			Summary:	Terminate the manager.
*			Arguments:
*				ProdTblProdTblInfoTy *prodTableInfoPtr	- Information
*			Function return: None.
*
*********************************************************************************/
void ProdTblTerminate(ProdTblProdTblInfoTy *prodTableInfoPtr)	
{
	
	/* Only do things, if we are initialized */
	if ((ProdTblIsInitialized) && (prodTableInfoPtr != 0)) {
	
		/* If we are not in an error situation, close the tables
			forcing the output files to be created.
		*/
		if (ErIsInError() == false) {
			/* Close the productivity table */
			if (ProdTblCloseTable(prodTableInfoPtr) == false) {
				ErHandle("Unable to close productivity table.", false);
			}
		}
		else {
			prodTblFreeTables();
		}
		
	}
	
	ProdTblIsInitialized = false;
}

#ifdef PHG_DEBUG
/*********************************************************************************
*
*			Name:			ProdTblDumpObjects
*
*			Summary:		Dump the objects for debuging purposes.
*
*			Arguments:
*
*			Function return: None.
*
*********************************************************************************/
void ProdTblDumpObjects()	
{
	Boolean					okay = false;		/* Process Flag */
	LbUsTwoByte				sliceIndex;			/* Current slice */
	LbUsFourByte			angleIndex;			/* Current angle index */
	
	#ifdef PHG_DEBUG
		if (ProdTblIsInitialized == false) {
			PhgAbort("You can't call ProdTblDumpObjects without"
				" initializing the ProdTbl module.", false);
		}
	#endif	

	if (PhgDebugDumpFile == 0)
	return;

	do { /* Process Loop */
	
		/* Print out the slice information */
		for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {

			/* Print out the slices productivities */
			fprintf(PhgDebugDumpFile, "Slice productivities (Primary), [productivity startOfBoundary endOfBoundary]\n");
			for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
				fprintf(PhgDebugDumpFile, "[%1.2f %1.2f %1.2f] ",
					ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].cellProductivity,
					ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary,
					ProdTblPrimProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary);
			}
			
			/* Print the blank line */
			fprintf(PhgDebugDumpFile, "\n\n");

			/* Print out the slices productivities */
			fprintf(PhgDebugDumpFile, "Slice productivities (Scatter), [productivity startOfBoundary endOfBoundary]\n");
			for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
				fprintf(PhgDebugDumpFile, "[%1.2f %1.2f %1.2f] ",
					ProdTblScatProdTbl[sliceIndex].productivity[angleIndex].cellProductivity,
					ProdTblScatProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary,
					ProdTblScatProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary);
			}
			
			/* Print the blank line */
			fprintf(PhgDebugDumpFile, "\n\n");

			/* Print out the slices productivities */
			fprintf(PhgDebugDumpFile, "Slice productivities (Max), [productivity startOfBoundary endOfBoundary]\n");
			for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
				fprintf(PhgDebugDumpFile, "[%1.2f %1.2f %1.2f] ",
					ProdTblMaxProdTbl[sliceIndex].productivity[angleIndex].cellProductivity,
					ProdTblMaxProdTbl[sliceIndex].productivity[angleIndex].startOfBoundary,
					ProdTblMaxProdTbl[sliceIndex].productivity[angleIndex].endOfBoundary);
			}
			
			/* Print the blank line */
			fprintf(PhgDebugDumpFile, "\n\n");
			
		} /* Loop to next slice */
		
		/* Print out the slice information */
		for (sliceIndex = 0; sliceIndex < ProdTblNumSlices; sliceIndex++) {
			
		
			/* Print weights */
			fprintf(PhgDebugDumpFile, "Sum of started Primary-Photon's squared weight by stratification cell.\n");
			for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
				fprintf(PhgDebugDumpFile, "%3.2e ", ProdTblStartPrimPhoSquWeights[sliceIndex][angleIndex]);
			}
			
			/* Print the blank line */
			fprintf(PhgDebugDumpFile, "\n\n");
			
			fprintf(PhgDebugDumpFile, "Sum of detected primary photons weights squared by strat cell.\n");
			for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
				fprintf(PhgDebugDumpFile, "%3.2e  ", ProdTblDetPrimPhoSquWeights[sliceIndex][angleIndex]);
			}

			/* Print the blank line */
			fprintf(PhgDebugDumpFile, "\n\n");
			
			/* Print weights */
			fprintf(PhgDebugDumpFile, "Sum of started scattered photons weights squared by strat cell.\n");
			for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
				fprintf(PhgDebugDumpFile, "%3.2e  ", ProdTblStartScatPhoSquWeights[sliceIndex][angleIndex]);
			}
						
			/* Print the blank line */
			fprintf(PhgDebugDumpFile, "\n\n");
			
			fprintf(PhgDebugDumpFile, "Sum of detected scattered photons weights squared by strat cell.\n");
			for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
				fprintf(PhgDebugDumpFile, "%3.2e ", ProdTblDetScatPhoSquWeights[sliceIndex][angleIndex]);
			}
			
			/* Print the blank line */
			fprintf(PhgDebugDumpFile, "\n\n");


			/* Print weights */
			fprintf(PhgDebugDumpFile, "Sum of started primary photons weights by strat cell.\n");
			for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
				fprintf(PhgDebugDumpFile, "%3.2e ", ProdTblStartPrimPhoWeights[sliceIndex][angleIndex]);
			}
			
			/* Print the blank line */
			fprintf(PhgDebugDumpFile, "\n\n");
			
			/* Print weights */
			fprintf(PhgDebugDumpFile, "Sum of started scattered photons weights by strat cell.\n");
			for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
				fprintf(PhgDebugDumpFile, "%3.2e  ", ProdTblStartScatPhoWeights[sliceIndex][angleIndex]);
			}

		
			/* Print the blank line */
			fprintf(PhgDebugDumpFile, "\n\n");
			
			fprintf(PhgDebugDumpFile, "Sum of detected primary photons weights by strat cell\n");
			for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
				fprintf(PhgDebugDumpFile, "%3.2e  ", ProdTblDetPrimPhoWeights[sliceIndex][angleIndex]);
			}
						
			/* Print the blank line */
			fprintf(PhgDebugDumpFile, "\n\n");
			
			fprintf(PhgDebugDumpFile, "Sum of detected scattered photons weights by strat cell.\n");
			for (angleIndex = 0; angleIndex < ProdTblNumAngleCells; angleIndex++) {
				fprintf(PhgDebugDumpFile, "%3.2e  ", ProdTblDetScatPhoWeights[sliceIndex][angleIndex]);
			}

			/* Print the blank line */
			fprintf(PhgDebugDumpFile, "\n\n");
			
		} /* Loop to next slice */

		okay = true;
	} while (false);
	
	if (!okay) {
		ErHandle("Dump failed.", false);
		PhgAbort("Program terminated (ProdTblDumpObjects).",
			false);
	}
}
#endif /* PHG_DEBUG */
#undef	PRODUCTIVITY_TABLE
