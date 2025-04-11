/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1994-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		display.header
*			Revision Number:	1.5
*			Date last revised:	24 July 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	March 10, 1994
*
*			Module Overview:	Displays the header of a PHG history file.
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
#include "LbHeader.h"

#include "Photon.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhgMath.h"
#include "PhoHFile.h"
#include "ProdTbl.h"
#include "PhoTrk.h"
#include "SubObj.h"
#include "EmisList.h"
#include "Collimator.h"
#include "Detector.h"
#include "phg.h"
#include "PhgBin.h"


/* LOCAL CONSTANTS */

/* LOCAL TYPES */

/* LOCAL GLOBALS */
static	float			disHdrVersions[] = {1.00};
static	PhoHFileHdrTy	header;					/* The history file header */

/* PROTOTYPES */
Boolean			DisplayHeader(int argc, char *argv[]);
void			display(PhoHFileHdrTy *headerPtr);

/* FUNCTIONS */

/**********************
*	DisplayHeader
*
*	Purpose:	Setup and display.
*
*	Result:	True unless an error occurs.
***********************/
Boolean DisplayHeader(int argc, char *argv[])
{
	Boolean				okay = false;			/* Process Loop */
	Boolean				canceled = false;		/* Canceled Loop */
	char				inputPath[1024];		/* Path to input file */
	char				errStr[1024];			/* Error string buffer */
	FILE				*historyFile;			/* Our data file */
	LbUsFourByte		headerSize;				/* Size of header in the file. */
	PhoHFileHdrKindTy	headerKind;				/* The type of this header */
	float				headerVersion;			/* The version of  this header */


	/* Clear program variables */
	historyFile = 0;
	
	do { /* Process Loop */
			
		/* Get path to input file */
		if (argc > 1) {
			strcpy(inputPath, argv[1]);
		}
		else {
			(void) LbInAsk("\nPlease enter name of input file", 0, false, &canceled,
				0, 0, 0, 0, inputPath);
				
			if (canceled)
				break;
		}
					
		/* Open data file file */
		if ((historyFile = LbFlFileOpen(inputPath, "rb")) == 0) {
			/* Set error message */
			sprintf(errStr, "Unable to open history file\n'%s'.",
				inputPath);
			ErStFileError(errStr);
			
			/* Exit process loop */
			break;
		}		
		
		/* First, read the header size; it is in the first four bytes */
		if (fread(&headerSize, sizeof(LbUsFourByte), 1, historyFile) != 1) {
			ErStFileError("\nUnable to read header size from input file.");
			break;
		}
		
		/* Verify we know how to deal with this */
		if (headerSize != sizeof(header)) {
			sprintf(errStr, "History file header size = %ld, does not match current version = %ld.",
				(unsigned long)headerSize, sizeof(header));
			ErStGeneric(errStr);
			break;
		}

		/* Read the header kind */
		if (fread(&headerKind, sizeof(PhoHFileHdrKindTy), 1, historyFile) != 1) {
			ErStFileError("\nUnable to read header type from input file.");
			break;
		}
		
		/* Verify we know how to deal with this */
		if (headerKind > PhoHFileEn_BIN) {
			sprintf(errStr, "History file header kind = %d, is unsupported for known versions.",
				headerKind);
			ErStGeneric(errStr);
			break;
		}
		
		/* Read the header version */
		if (fread(&headerVersion, sizeof(float), 1, historyFile) != 1) {
			ErStFileError("\nUnable to read header version from input file.");
			break;
		}
		
		/* Verify we know how to deal with this */
		if (headerVersion != disHdrVersions[0]) {
			sprintf(errStr, "History file header version = %f, is unsupported.",
				headerVersion);
			ErStGeneric(errStr);
			break;
		}
		
		/* Reset to zero and read in the header */
		if (fseek(historyFile, 0, SEEK_SET) != 0) {
			ErStFileError("\nUnable to reset to beginning of history file.");
			break;
		}

		/* Read in the header */
		if (fread(&header, sizeof(header), 1, historyFile) != 1) {
			ErStFileError("\nUnable to read header from history file.");
			break;
		}
		
		/* Display the header */
		display(&header);
		

		/* Set process flag */
		okay = true;
	} while (false);
	
	/* Try to close history file if opened */
	if (historyFile != 0)
		fclose(historyFile);
		
	if (!okay && canceled) {
		ErHandle("\nUser canceled program\n", false);
	}
	
	/* Return the results */
	return (okay);
}
/**********************
*	display
*
*	Purpose:	Display the header.
*	Arguments:
*		PhoHFileHdrTy	*headerPtr		- The header.
*	Result:	None.
***********************/
void display(PhoHFileHdrTy *headerPtr)
{

	/* Print message indicating type of header */
	switch(headerPtr->H.HdrKind) {
	
		case PhoHFileEn_PHG:
	
			fprintf(stdout, "\nHeader is for PHG History File");
			break;
			
		case PhoHFileEn_COL:
	
			fprintf(stdout, "\nHeader is for Collimator History File");
			break;
			
		case PhoHFileEn_DET:
	
			fprintf(stdout, "\nHeader is for Detector History File");
			break;
			
		case PhoHFileEn_BIN_CT:
	
			fprintf(stdout, "\nHeader is for Binning Module Histogram File of Counts");
			break;
			
		case PhoHFileEn_BIN_WT:
	
			fprintf(stdout, "\nHeader is for Binning Module Histogram File of Weights");
			break;
			
		case PhoHFileEn_BIN_WTSQ:
	
			fprintf(stdout, "\nHeader is for Binning Module Histogram File of Sum of Weights Squared");
			break;
			
		case PhoHFileEn_BIN:
	
			fprintf(stdout, "\nHeader is for Binning Module History File");
			break;
			
		default:

			fprintf(stdout, "\nUnsupported value in 'HdrKind' field of header");
			break;
	}	

	/* Print the fields */
	fprintf(stdout, "\nHeader Version\t\t\t\t:%.2f",headerPtr->H.HdrVersion);
	fprintf(stdout, "\nHeader Size\t\t\t\t:%ld",(unsigned long)headerPtr->H.HdrSize);
			
			
			
	fprintf(stdout, "\nEvents to Simulate\t\t\t:%lld",headerPtr->H.PhgRunTimeParams.Phg_EventsToSimulate);
	fprintf(stdout, "\nRandom Seed\t\t\t\t:%ld",(unsigned long)headerPtr->H.PhgRunTimeParams.PhgRandomSeed);
	fprintf(stdout, "\nLength of Scan (seconds)\t\t:%9.3f",headerPtr->H.PhgRunTimeParams.Phg_LengthOfScan);
	fprintf(stdout, "\nAcceptance Angle (degrees)\t\t:%3.2f",headerPtr->H.PhgRunTimeParams.Phg_AcceptanceAngle);
	fprintf(stdout, "\nSine of Acceptance Angle\t\t:%3.2f",headerPtr->H.PhgRunTimeParams.Phg_SineOfAcceptanceAngle);
	fprintf(stdout, "\nMinimum Energy (keV)\t\t\t:%3.1f",headerPtr->H.PhgRunTimeParams.PhgMinimumEnergy);
	fprintf(stdout, "\nMin Weight Window Ratio\t\t\t:%lf",headerPtr->H.PhgRunTimeParams.PhgMinWWRatio);
	fprintf(stdout, "\nMax Weight Window Ratio\t\t\t:%lf",headerPtr->H.PhgRunTimeParams.PhgMaxWWRatio);
	fprintf(stdout, "\nNuclide Photon Energy\t\t\t:%3.1f",headerPtr->H.PhgRunTimeParams.PhgNuclide.photonEnergy_KEV);

	fprintf(stdout, "\nForced Detection\t\t\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsForcedDetection == true) ? "TRUE" : "FALSE");
	fprintf(stdout, "\nStratification\t\t\t\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsStratification == true) ? "TRUE" : "FALSE");
	fprintf(stdout, "\nForced NonAbsorbtion\t\t\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsForcedNonAbsorbtion == true) ? "TRUE" : "FALSE");
	fprintf(stdout, "\nSPECT\t\t\t\t\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsSPECT == true) ? "TRUE" : "FALSE");
	fprintf(stdout, "\nPET coincidences only\t\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsPETCoincidencesOnly == true) ? "TRUE" : "FALSE");
	fprintf(stdout, "\nPET coincidences plus singles\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsPETCoincPlusSingles == true) ? "TRUE" : "FALSE");
/*	fprintf(stdout, "\nMulti emission isotope\t\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsMultiEmission == true) ? "TRUE" : "FALSE"); */
	fprintf(stdout, "\nCreate History File\t\t\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsHistoryFile == true) ? "TRUE" : "FALSE");
	fprintf(stdout, "\nAdjust For Positron Range\t\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsAdjForPosRange == true) ? "TRUE" : "FALSE");
	fprintf(stdout, "\nAdjust For Collinearity\t\t\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsAdjForCollinearity == true) ? "TRUE" : "FALSE");
	fprintf(stdout, "\nPre-Computed Productivities\t\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsComputedProductivityTbl == true) ? "TRUE" : "FALSE");
	fprintf(stdout, "\nPoint Source Voxels\t\t\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsVoxelPointSource == true) ? "TRUE" : "FALSE");
	fprintf(stdout, "\nBinning on the Fly\t\t\t:%s", (headerPtr->H.PhgRunTimeParams.PhgIsBinOnTheFly == true) ? "TRUE" : "FALSE");

	fprintf(stdout, "\nParameter File\t\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgParamFilePath);
	fprintf(stdout, "\nActivity Index File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgSubObjActivityIndexFilePath);
	fprintf(stdout, "\nActivity Table File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgSubObjActivityTableFilePath);
	fprintf(stdout, "\nActivity Index Translation File\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgSubObjActIndexTransFilePath);
	fprintf(stdout, "\nActivity Image File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgSubObjActImgFilePath);
	fprintf(stdout, "\nAttenuation Index File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgSubObjAttenIndexFilePath);
	fprintf(stdout, "\nAttenuation Table File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgSubObjAttenTableFilePath);
	fprintf(stdout, "\nAttenuation Index File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgSubObjAttIndexTransFilePath);
	fprintf(stdout, "\nAttenuation Image File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgSubObjAttImgFilePath);
	fprintf(stdout, "\nAttenuation Image File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgSubObjActivityTableFilePath);
	fprintf(stdout, "\nProductivity Input File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgProdTblInputTableFilePath);
	fprintf(stdout, "\nProductivity Ouput File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgProdTblOutputTableFilePath);
	fprintf(stdout, "\nForced Detection Table File\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgPhoTrkForcedDetectionFilePath);
	fprintf(stdout, "\nStatistics Ouput File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgPhoHStatFilePath);
	fprintf(stdout, "\nPhoton History File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgPhoHFileHistoryFilePath);
	fprintf(stdout, "\nBinning Parameters File\t\t\t:%s",headerPtr->H.PhgRunTimeParams.PhgBinParamsFilePath[PhgNumBinParams]);

	fprintf(stdout, "\nBinning \"count\" Image File\t\t:%s",headerPtr->H.BinRunTimeParams.countImgFilePath);
	fprintf(stdout, "\nBinning \"weight\" Image File\t\t:%s",headerPtr->H.BinRunTimeParams.weightImgFilePath);
	fprintf(stdout, "\nBinning \"weight squared\" Image File\t:%s",headerPtr->H.BinRunTimeParams.weightSquImgFilePath);
	fprintf(stdout, "\nBinning: number of Z bins %ld",(unsigned long)headerPtr->H.BinRunTimeParams.numZBins);
	fprintf(stdout, "\nBinning: number of PA bins %ld",(unsigned long)headerPtr->H.BinRunTimeParams.numPABins);
	fprintf(stdout, "\nBinning: number of TD bins %ld",(unsigned long)headerPtr->H.BinRunTimeParams.numTDBins);
	fprintf(stdout, "\nBinning: number of AA bins %ld",(unsigned long)headerPtr->H.BinRunTimeParams.numAABins);
	fprintf(stdout, "\nBinning: number of Eenergy(1) bins %ld",(unsigned long)headerPtr->H.BinRunTimeParams.numE1Bins);
	fprintf(stdout, "\nBinning: number of Energy(2) bins %ld",(unsigned long)headerPtr->H.BinRunTimeParams.numE2Bins);
	fprintf(stdout, "\nBinning: number of Scatter(1) bins %ld",(unsigned long)headerPtr->H.BinRunTimeParams.numS1Bins);
	fprintf(stdout, "\nBinning: number of Scatter(2) bins %ld",(unsigned long)headerPtr->H.BinRunTimeParams.numS2Bins);
	fprintf(stdout, "\nBinning: number of image bins %ld",(unsigned long)headerPtr->H.BinRunTimeParams.numImageBins);
	fprintf(stdout, "\nBinning: scatter/random parameter %ld",(unsigned long)headerPtr->H.BinRunTimeParams.scatterRandomParam);
	fprintf(stdout, "\nBinning: min Z position %3.3f",headerPtr->H.BinRunTimeParams.minZ);
	fprintf(stdout, "\nBinning: max Z position %3.3f",headerPtr->H.BinRunTimeParams.maxZ);
	fprintf(stdout, "\nBinning: min Polar Aangle %3.3f",headerPtr->H.BinRunTimeParams.minPA);
	fprintf(stdout, "\nBinning: max Polar Angle %3.3f",headerPtr->H.BinRunTimeParams.maxPA);
	fprintf(stdout, "\nBinning: min Transaxial Distance %3.3f",headerPtr->H.BinRunTimeParams.minTD);
	fprintf(stdout, "\nBinning: max Transaxial Distance %3.3f",headerPtr->H.BinRunTimeParams.maxTD);
	fprintf(stdout, "\nBinning: min Azimuthal Angle %3.3f",headerPtr->H.BinRunTimeParams.minAA);
	fprintf(stdout, "\nBinning: max Azimuthal Angle %3.3f",headerPtr->H.BinRunTimeParams.maxAA);
	fprintf(stdout, "\nBinning: min Energy %3.3f",headerPtr->H.BinRunTimeParams.minE);
	fprintf(stdout, "\nBinning: max Energy %3.3f",headerPtr->H.BinRunTimeParams.maxE);
	fprintf(stdout, "\nBinning: min number of scatters %ld",(unsigned long)headerPtr->H.BinRunTimeParams.minS);
	fprintf(stdout, "\nBinning: max number of scatters %ld",(unsigned long)headerPtr->H.BinRunTimeParams.maxS);
	fprintf(stdout, "\nBinning: %s to existing images",
		((headerPtr->H.BinRunTimeParams.addToExistingImg == true) ? "adding" : "not adding"));
	fprintf(stdout, "\nBinning: %s count image",
		((headerPtr->H.BinRunTimeParams.doCounts == true) ? "creating" : "not creating"));
	fprintf(stdout, "\nBinning: %s weight image",
		((headerPtr->H.BinRunTimeParams.doWeights == true) ? "creating" : "not creating"));
	fprintf(stdout, "\nBinning: %s weight squared image",
		((headerPtr->H.BinRunTimeParams.doWeightsSquared == true) ? "creating" : "not creating"));
	fprintf(stdout, "\nBinning: range on Z values %3.3f",headerPtr->H.BinRunTimeParams.zRange);
	fprintf(stdout, "\nBinning: range on energy values %3.3f",headerPtr->H.BinRunTimeParams.eRange);
	fprintf(stdout, "\nBinning: range on scatter values %3.3f",headerPtr->H.BinRunTimeParams.sRange);
	fprintf(stdout, "\nBinning: range on transaxial distance values %3.3f",headerPtr->H.BinRunTimeParams.tdRange);
	fprintf(stdout, "\nBinning: size of scatter(2) space in count image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.scatter2CIsize);
	fprintf(stdout, "\nBinning: size of scatter(2) space  in weight image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.scatter2WIsize);
	fprintf(stdout, "\nBinning: size of scatter(2) space  in weight squared image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.scatter2WISsize);
	fprintf(stdout, "\nBinning: size of scatter(1) space in count image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.scatter1CIsize);
	fprintf(stdout, "\nBinning: size of scatter(1) space  in weight image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.scatter1WIsize);
	fprintf(stdout, "\nBinning: size of scatter(1) space  in weight squared image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.scatter1WISsize);
	fprintf(stdout, "\nBinning: size of energy(2) space in count image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.energy2CIsize);
	fprintf(stdout, "\nBinning: size of energy(2) space  in weight image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.energy2WIsize);
	fprintf(stdout, "\nBinning: size of energy(2) space  in weight squared image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.energy2WISsize);
	fprintf(stdout, "\nBinning: size of energy(1) space in count image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.energy1CIsize);
	fprintf(stdout, "\nBinning: size of energy(1) space  in weight image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.energy1WIsize);
	fprintf(stdout, "\nBinning: size of energy(1) space  in weight squared image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.energy1WISsize);
	fprintf(stdout, "\nBinning: size of azimuthal angle space in count image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.aaCIsize);
	fprintf(stdout, "\nBinning: size of azimuthal angle space  in weight image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.aaWIsize);
	fprintf(stdout, "\nBinning: size of azimuthal angle space  in weight squared image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.aaWISsize);
	fprintf(stdout, "\nBinning: size of transaxial distance space in count image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.tdCIsize);
	fprintf(stdout, "\nBinning: size of transaxial distance space  in weight image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.tdWIsize);
	fprintf(stdout, "\nBinning: size of transaxial distance space  in weight squared image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.tdWISsize);
	fprintf(stdout, "\nBinning: size of polar angle space in count image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.paCIsize);
	fprintf(stdout, "\nBinning: size of polar angle space  in weight image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.paWIsize);
	fprintf(stdout, "\nBinning: size of polar angle space  in weight squared image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.paWISsize);
	fprintf(stdout, "\nBinning: size of z-axis position(2) space in count image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.z2CIsize);
	fprintf(stdout, "\nBinning: size of z-axis position(2) space  in weight image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.z2WIsize);
	fprintf(stdout, "\nBinning: size of z-axis position(2) space  in weight squared image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.z2WISsize);
	fprintf(stdout, "\nBinning: size of z-axis position(1) space in count image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.z1CIsize);
	fprintf(stdout, "\nBinning: size of z-axis position(1) space  in weight image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.z1WIsize);
	fprintf(stdout, "\nBinning: size of z-axis position(1) space  in weight squared image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.z1WISsize);
	fprintf(stdout, "\nBinning: count image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.countImageSize);
	fprintf(stdout, "\nBinning: size of weight image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.weightImageSize);
	fprintf(stdout, "\nBinning: size of weight squared image %ld",(unsigned long)headerPtr->H.BinRunTimeParams.weightSquImageSize);
	fprintf(stdout, "\nBinning: weight image type %ld",(unsigned long)headerPtr->H.BinRunTimeParams.weight_image_type);
	fprintf(stdout, "\nBinning: count image type %ld",(unsigned long)headerPtr->H.BinRunTimeParams.count_image_type);



		/* Print out target cylinder */
	fprintf(stdout, "\n\nTarget Cylinder\n\tradius\t= % 3.2f\n\tzMin\t= % 3.2f\n\tzMax\t= % 3.2f\n",
		headerPtr->H.TargetCylinder.radius, headerPtr->H.TargetCylinder.zMin, headerPtr->H.TargetCylinder.zMax);

	fprintf(stdout, "\tcenterX\t= % 3.2f\tcenterY\t= % 3.2f\n",
		headerPtr->H.TargetCylinder.centerX, headerPtr->H.TargetCylinder.centerY);

	/* Print out object cylinder */
	fprintf(stdout, "\nObject Cylinder\n\tradius\t= % 3.2f\n\tzMin\t\t= % 3.2f\tzMax\t= % 3.2f\n",
		headerPtr->H.ObjectCylinder.radius, headerPtr->H.ObjectCylinder.zMin, headerPtr->H.ObjectCylinder.zMax);

	fprintf(stdout, "\tcenterX\t= % 3.2f\tcenterY\t= % 3.2f\n",
		headerPtr->H.ObjectCylinder.centerX, headerPtr->H.ObjectCylinder.centerY);

	/* Print out limit cylinder */
	fprintf(stdout, "\nLimit Cylinder\n\tradius\t= % 3.2f\n\tzMin\t\t= % 3.2f\tzMax\t= % 3.2f\n",
		headerPtr->H.LimitCylinder.radius, headerPtr->H.LimitCylinder.zMin, headerPtr->H.LimitCylinder.zMax);

	fprintf(stdout, "\tcenterX\t= % 3.2f\tcenterY\t= % 3.2f\n",
		headerPtr->H.LimitCylinder.centerX, headerPtr->H.LimitCylinder.centerY);

	/* Print out critical zone */
	fprintf(stdout, "\nCritical Zone\n\tradius\t= % 3.2f\n\tzMin\t= % 3.2f\tzMax\t= % 3.2f\n",
		headerPtr->H.CriticalZone.radius, headerPtr->H.CriticalZone.zMin, headerPtr->H.CriticalZone.zMax);

	fprintf(stdout, "\tcenterX\t= % 3.2f\tcenterY\t= % 3.2f\n",
		headerPtr->H.CriticalZone.centerX, headerPtr->H.CriticalZone.centerY);


	fprintf(stdout, "\nNumber of Photons in History\t\t:%lld",headerPtr->H.NumPhotons);
	fprintf(stdout, "\nNumber of Decays in History\t\t:%lld",headerPtr->H.NumDecays);
	fprintf(stdout, "\nNumber of Simulations in History\t:%ld",(unsigned long)headerPtr->H.NumSimulations);
	fprintf(stdout, "\nSum of Events to Simulate in History\t:%lld",headerPtr->H.SumEventsToSimulate);

	fprintf(stdout, "\n");

}
