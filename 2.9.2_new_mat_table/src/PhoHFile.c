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
*			Module Name:		PhoHFile.c
*			Revision Number:	1.6
*			Date last revised:	10 June 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	24 August, 1992
*
*			Module Overview:	Photon History File management routines.
*
*			References:			'Photon History File Processes' PHG design.
*
**********************************************************************************
*
*			Global functions defined:	
*				PhoHFileClose
*				PhoHFileCreate
*				PhoHFilePrintParams
*				PhoHFilePrintReport
*				PhoHFileWriteDetections
*				PhoHFileReadFields
*				PhoHFilePrintFields
*				PhoHFileGetRunTimeParams
*				PhoHFileLookupRunTimeParamLabel
*				PhoHFileReadEvent
*				PhoHFileOldReadEvent
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
*			Programmer(s):      Robert Harrison
*
*			Revision date:      1 Jan 2013
*
*			Revision description:   Write out the detector crystal number to
*                           the history file.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		24 February 2012
*
*			Revision description:	Modified PhoHFileClose to first check if 
*										file already closed
*									Corrected old photon type 2 reading in 
*										PhoHFileOldReadEvent
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		2007
*
*			Revision description:	
*				- New custom history file options:  decay position, decay time,
*				detector crystal, detector angle.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		16 July 2004
*
*			Revision description:	Added the number of scatters in the collimator
*					to the number of scatters written and read in standard list
*					mode.  Added PhoHFileOldReadEvent to read old list mode
*					files.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		5 January 2004
*
*			Revision description:	Added PhoHFileReadEvent
*
*********************************************************************************/

#define	PHOTON_HIST_FILE

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbParamFile.h"
#include "LbHeader.h"
#include "LbFile.h"
#include "LbInterface.h"

#include "Photon.h"
#include "PhgParams.h"
#include "PhgMath.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhoHFile.h"
#include "PhgHdr.h"
#include "phg.h"
#include "PhgBin.h"

#define MAX_SCATTERS			20		/* Maximum number of scatters accounted for */

typedef char labelTy[LBPF_LABEL_LEN];

/* Local globals */
static char			*phoHFileOpenMode =	"wb";	/* File access mode */
static char			phoHFileErrString[1024];	/* Error string storage */
static LbUsOneByte	PhoHFileDecayFlag;			/* Identification flags */
static LbUsOneByte	PhoHFilePhotonFlag;			/* Identification flags */

/* see comment above PhoHFileEn_RunTimeParamsTy (in PhoHFile.h) when altering
 this variable.  They must be altered in conjunction */
static labelTy		PhoHFileRunTimeParamLabels[] = {	/* Our label table */
					"x_position",
					  "x_pos_min",
					  "x_pos_max",
					"y_position",
					  "y_pos_min",
					  "y_pos_max",
					"z_position",
					  "z_pos_min",
					  "z_pos_max",
					"x_cosine",
					  "x_cos_min",
					  "x_cos_max",
					"y_cosine",
					  "y_cos_min",
					  "y_cos_max",
					"z_cosine",
					  "z_cos_min",
					  "z_cos_max",
					"scatters_in_object",
					  "obj_scatters_min",
					  "obj_scatters_max",
					"scatters_in_collimator",
					  "col_scatters_min",
					  "col_scatters_max",
					"decay_weight",
					"weight",
					"energy",
					  "energy_min",
					  "energy_max",
					"travel_distance",
					  "travel_distance_min",
					  "travel_distance_max",
					"transaxial_distance",
					  "transaxial_distance_min",
					  "transaxial_distance_max",
					"azimuthal_angle_index",
					  "aa_index_min",
					  "aa_index_max",
					"axial_position",
					  "axial_position_min",
					  "axial_position_max",
					"detector_x_position",
					  "detector_x_position_min",
					  "detector_x_position_max",
					"detector_y_position",
					  "detector_y_position_min",
					  "detector_y_position_max",
					"detector_z_position",
				  	  "detector_z_position_min",
					  "detector_z_position_max",
					"num_detector_interactions",
					  "num_detector_interactions_min",
					  "num_detector_interactions_max",
					"det_interaction_positions",
					"detector_angle",
					"detector_angle_min",
					"detector_angle_max",
					"decay_x_position",
					  "decay_x_pos_min",
					  "decay_x_pos_max",
					"decay_y_position",
					  "decay_y_pos_min",
					  "decay_y_pos_max",
					"decay_z_position",
					  "decay_z_pos_min",
					  "decay_z_pos_max",
					"decay_time",
					"decay_type",
					"detector_crystal",
					""
					};

/* Prototypes */
Boolean	phoHFileWriteFields(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay,
			PHG_TrackingPhoton *photon, LbUsOneByte *numAccepted);

void phoHFileCreateDetectedPhoton(PHG_TrackingPhoton *trackingphoton,
		PHG_DetectedPhoton *detectedphoton);
			
/*********************************************************************************
*
*			Name:			PhoHFileClose
*
*			Summary:		Close the History file.
*
*			Arguments:
*
*			Function return: TRUE unless an error occurs.
*
*********************************************************************************/
Boolean PhoHFileClose(PhoHFileHkTy *hdrHkTyPtr)	
{
	Boolean	okay = false;	/* Success flag */
	
	if (! hdrHkTyPtr->histFile) {
		/* File isn't open, so already closed */
		okay = true;
	}
	else {
		do { /* Process Loop */

			/* Update the image header, if its file exists */
			if (hdrHkTyPtr->headerHk.fileRef) {
				if (PhgHdrUpHeader(hdrHkTyPtr->histFile, &hdrHkTyPtr->header,
						&hdrHkTyPtr->headerHk) == false) {
					break;
				}
			}
			PhgHdrFrHeader(&hdrHkTyPtr->headerHk);
			
			/* Close the binary file */
			if (fclose(hdrHkTyPtr->histFile) != 0) {
				ErStGeneric("Error closing history binary file.");
				break;
			}
			
			okay = true;
		} while (false);
	}
	
	/* Clear the file hook to prevent abuse */
	hdrHkTyPtr->histFile = 0;
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			PhoHFileCreate
*
*			Summary:		Create a new log file.
*
*			Arguments:
*				char				*histFilePath		- The path name of the file to create
*				char				*histParamsFilePath	- The path name of the custom parameters
*				PhoHFileHdrKindTy	HdrKind				- The type of the header 
*				PhoHFileHkTy 		*hdrHkTyPtr			- The header hook created
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean PhoHFileCreate(char *histFilePath, char *histParamsFilePath,
			PhoHFileHdrKindTy HdrKind, PhoHFileHkTy *hdrHkTyPtr)	
{
	Boolean	okay = false;		/* Process flags */
	
	do { /* Process Loop */
		
		/* Set our type flags, used later on for marking what type of event is being written */
		PHG_SetIsADecay(PhoHFileDecayFlag);
		PHG_SetIsAPhoton(PhoHFilePhotonFlag);
		
		/* Create the binary file */
		if ((hdrHkTyPtr->histFile = LbFlFileOpen(histFilePath, phoHFileOpenMode)) == 0) {
			sprintf(phoHFileErrString, "Unable to open history file named '%s'",
				histFilePath);
			ErStFileError(phoHFileErrString);
			break;
		}
				
		/* File is newly created so generate an empty header now */
		if (PhgHdrCreateRunTimeHeader(HdrKind, &(hdrHkTyPtr->header), (PHG_BinParamsTy *) 0) == false) {
			break;
		}
		
		/* Create a new header hook for the file */
		if (PhgHdrMkHeader(hdrHkTyPtr->histFile, &(hdrHkTyPtr->header),
				&(hdrHkTyPtr->headerHk)) == false){
			break;
		}
		
		/* Do custom parameter initialization if they ask for it */
		if (histParamsFilePath[0] != '\0') {
			
			if (PhoHFileGetRunTimeParams(histParamsFilePath, &(hdrHkTyPtr->customParams)) == false) {
				sprintf(phoHFileErrString,"Unable to get custom parameters for history file named '%s'",
					histParamsFilePath);
				ErAlert(phoHFileErrString, false);
				break;
			}
			
			hdrHkTyPtr->doCustom = true;
			strcpy(hdrHkTyPtr->customParamsName, histParamsFilePath);
		}
		else {
			hdrHkTyPtr->doCustom = false;
			hdrHkTyPtr->customParamsName[0] = '\0';
		}
		
		hdrHkTyPtr->bluesReceived = 0;
		hdrHkTyPtr->bluesAccepted = 0;
		hdrHkTyPtr->pinksReceived = 0;
		hdrHkTyPtr->pinksAccepted = 0;
		
		okay = true;
	} while (false);
	return (okay);
}


/*********************************************************************************
*
*			Name:			PhoHFilePrintParams
*
*			Summary:		Print parameters for a given history file.
*
*			Arguments:
*				PhoHFileHkTy 		*hdrHkTyPtr			- The header hook
*
*			Function return: none.
*
*********************************************************************************/
void PhoHFilePrintParams(PhoHFileHkTy *hdrHkTyPtr)	
{
	PhoHFileRunTimeParamsTy	*params = &hdrHkTyPtr->customParams;	/* Just a short cut for accessing the parameters */

	LbInPrintf("\n\n\t\tParameters for history file '%s'\n", hdrHkTyPtr->customParamsName);

	if (hdrHkTyPtr->customParams.doXposition == true) {
		LbInPrintf("\n\t\tWriting X-position");
	}
	if (hdrHkTyPtr->customParams.doXrange == true) {
		LbInPrintf("\n\t\t\tX-position range = [%3.2f, %3.2f] cm", params->xPosMin, params->xPosMax);
	}

	if (hdrHkTyPtr->customParams.doYposition == true) {
		LbInPrintf("\n\t\tWriting Y-position");
	}

	if (hdrHkTyPtr->customParams.doYrange == true) {
		LbInPrintf("\n\t\t\tY-position range = [%3.2f, %3.2f] cm", params->yPosMin, params->yPosMax);
	}

	if (hdrHkTyPtr->customParams.doZposition == true) {
		LbInPrintf("\n\t\tWriting Z-position");
	}

	if (hdrHkTyPtr->customParams.doZrange == true) {
		LbInPrintf("\n\t\t\tZ-position range = [%3.2f, %3.2f] cm", params->zPosMin, params->zPosMax);
					
	}

	if (hdrHkTyPtr->customParams.doXcosine == true) {
		LbInPrintf("\n\t\tWriting cosine-X");
	}

	if (hdrHkTyPtr->customParams.doXCosRange == true) {
		LbInPrintf("\n\t\t\tX-cosine range = [%3.2f, %3.2f] rd", params->xCosMin, params->xCosMax);
	}

	if (hdrHkTyPtr->customParams.doYcosine == true) {
		LbInPrintf("\n\t\tWriting cosine-Y");
	}

	if (hdrHkTyPtr->customParams.doYCosRange == true) {
		LbInPrintf("\n\t\t\tY-cosine range = [%3.2f, %3.2f] rd", params->yCosMin, params->yCosMax);
	}

	if (hdrHkTyPtr->customParams.doZcosine == true) {
		LbInPrintf("\n\t\tWriting cosine-Z");
	}

	if (hdrHkTyPtr->customParams.doZCosRange == true) {
		LbInPrintf("\n\t\t\tZ-cosine range = [%3.2f, %3.2f] rd", params->zCosMin, params->zCosMax);
	}

	if (hdrHkTyPtr->customParams.doScattersInObject == true) {
		LbInPrintf("\n\t\tWriting number of scatters in object");
	}

	if (hdrHkTyPtr->customParams.doObjScatRange == true) {
		LbInPrintf("\n\t\t\tNumber of scatters in object range = [%d, %d] ", params->objScattersMin, params->objScattersMax);
	}

	if (hdrHkTyPtr->customParams.doScattersInCollimator == true) {
		LbInPrintf("\n\t\tWriting number of scatters in collimator");
	}

	if (hdrHkTyPtr->customParams.doColScatRange == true) {
		LbInPrintf("\n\t\t\tNumber of scatters in collimator range = [%d, %d] ", params->colScattersMin, params->colScattersMax);
	}

	if (hdrHkTyPtr->customParams.doDecayWeight == true) {
		LbInPrintf("\n\t\tWriting decay weight");
	}

	if (hdrHkTyPtr->customParams.doWeight == true) {
		LbInPrintf("\n\t\tWriting photon weight");
	}

	if (hdrHkTyPtr->customParams.doEnergy == true) {
		LbInPrintf("\n\t\tWriting energy");
	}

	if (hdrHkTyPtr->customParams.doEnergyRange == true) {
		LbInPrintf("\n\t\t\tEnergy range = [%3.2f, %3.2f] keV", params->energyMin, params->energyMax);
	}

	if (hdrHkTyPtr->customParams.doTravelDistance == true) {
		LbInPrintf("\n\t\tWriting travel distance");
	}

	if (hdrHkTyPtr->customParams.doTravDistRange == true) {
		LbInPrintf("\n\t\t\tTravel distance range = [%3.2f, %3.2f] cm", params->travelDistanceMin, params->travelDistanceMax);
	}

	if (hdrHkTyPtr->customParams.doDecayXposition == true) {
		LbInPrintf("\n\t\tWriting decay's X-position");
	}

	if (hdrHkTyPtr->customParams.doDecayXrange == true) {
		LbInPrintf("\n\t\t\tDecay's X-position range = [%3.2f, %3.2f] cm", params->decayXPosMin, params->decayXPosMax);
	}

	if (hdrHkTyPtr->customParams.doDecayYposition == true) {
		LbInPrintf("\n\t\tWriting decay's Y-position");
	}

	if (hdrHkTyPtr->customParams.doDecayYrange == true) {
		LbInPrintf("\n\t\t\tDecay's Y-position range = [%3.2f, %3.2f] cm", params->decayYPosMin, params->decayYPosMax);
	}

	if (hdrHkTyPtr->customParams.doDecayZposition == true) {
		LbInPrintf("\n\t\tWriting decay's Z-position");
	}

	if (hdrHkTyPtr->customParams.doDecayZrange == true) {
		LbInPrintf("\n\t\t\tDecay's Z-position range = [%3.2f, %3.2f] cm", params->decayZPosMin, params->decayZPosMax);
	}

	if (hdrHkTyPtr->customParams.doDecayTime == true) {
		LbInPrintf("\n\t\tWriting decay time");
	}

	if (hdrHkTyPtr->customParams.doDecayType == true) {
		LbInPrintf("\n\t\tWriting decay type");
	}

	if (hdrHkTyPtr->customParams.doTransaxialDistance == true) {
		LbInPrintf("\n\t\tWriting transaxial distance");
	}

	if (hdrHkTyPtr->customParams.doTransDistRange == true) {
		LbInPrintf("\n\t\t\tTransaxial distance range = [%3.2f, %3.2f] cm", params->transaxialDistanceMin, params->transaxialDistanceMax);
	}

	if (hdrHkTyPtr->customParams.doAzimuthalAngleIndex == true) {
		LbInPrintf("\n\t\tWriting azimuthal angle index");
	}

	if (hdrHkTyPtr->customParams.doAziAngRange == true) {
		LbInPrintf("\n\t\t\tAzimuthal angle index range = [%d, %d] ", params->aaIndexMin, params->aaIndexMax);
	}

	if (hdrHkTyPtr->customParams.doAxialPosition == true) {
		LbInPrintf("\n\t\tWriting axial position");
	}

	if (hdrHkTyPtr->customParams.doAxiPosRange == true) {
		LbInPrintf("\n\t\t\tAxial poition range = [%3.2f, %3.2f] cm", params->axialPositionMin, params->axialPositionMax);
	}

	if (hdrHkTyPtr->customParams.doDetectorXposition == true) {
		LbInPrintf("\n\t\tWriting detector X-position");
	}

	if (hdrHkTyPtr->customParams.doDetectorXposRange == true) {
		LbInPrintf("\n\t\t\tDetector X-position range = [%3.2f, %3.2f] cm", params->detectorXpositionMin, params->detectorXpositionMax);
	}

	if (hdrHkTyPtr->customParams.doDetectorYposition == true) {
		LbInPrintf("\n\t\tWriting detector Y-position");
	}

	if (hdrHkTyPtr->customParams.doDetectorYposRange == true) {
		LbInPrintf("\n\t\t\tDetector Y-position range = [%3.2f, %3.2f] cm", params->detectorYpositionMin, params->detectorYpositionMax);
	}

	if (hdrHkTyPtr->customParams.doDetectorZposition == true) {
		LbInPrintf("\n\t\tWriting detector Z-position");
	}

	if (hdrHkTyPtr->customParams.doDetectorZposRange == true) {
		LbInPrintf("\n\t\t\tDetector Z-position range = [%3.2f, %3.2f] cm", params->detectorZpositionMin, params->detectorZpositionMax);
	}

	if (hdrHkTyPtr->customParams.doDetectorAngle == true) {
		LbInPrintf("\n\t\tWriting detector angle");
	}

	if (hdrHkTyPtr->customParams.doDetectorAngleRange == true) {
		LbInPrintf("\n\t\t\tDetector angle range = [%3.2f, %3.2f] radians", params->detectorAngleMin, params->detectorAngleMax);
	}

	if (hdrHkTyPtr->customParams.doNumDetectorInteractions == true) {
		LbInPrintf("\n\t\tWriting number of detector interactions");
	}

	if (hdrHkTyPtr->customParams.doDetectorCrystal == true) {
		LbInPrintf("\n\t\tWriting detector crystal for detection");
	}

}

/*********************************************************************************
*
*			Name:			PhoHFilePrintReport
*
*			Summary:		Print report for a given history file.
*
*			Arguments:
*				PhoHFileHkTy 		*hdrHkTyPtr			- The header hook
*
*			Function return: none.
*
*********************************************************************************/
void PhoHFilePrintReport(PhoHFileHkTy *hdrHkTyPtr)	
{
	LbInPrintf("\n\n\tNumber of blue photons reaching history file = %lld", hdrHkTyPtr->bluesReceived);
	LbInPrintf("\n\tNumber of blue photons written to history file = %lld", hdrHkTyPtr->bluesAccepted);
	LbInPrintf("\n\tNumber of pink photons reaching history file = %lld", hdrHkTyPtr->pinksReceived);
	LbInPrintf("\n\tNumber of pink photons written to history file = %lld", hdrHkTyPtr->pinksAccepted);
	LbInPrintf("\n");
	
}

/*********************************************************************************
*
*			Name:			PhoHFileWriteDetections
*
*			Summary:		Write the detected decay and photons.
*
*			Arguments:
*				PhoHFileHkTy 		*hdrHkTyPtr		- The history hook
*				PHG_Decay			decay			- The decay that started the process.
*				PHG_TrackingPhoton *bluePhotons		- The blue photons detected.
*				LbUsFourByte 		numBluePhotons	- The number of blue photons.
*				PHG_TrackingPhoton *pinkPhotons		- The  pink photons detected.
*				LbUsFourByte		numPinkPhotons	- The number of pink photons.
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean PhoHFileWriteDetections(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay,
			PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
			PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons)
{
	Boolean				okay = false;
	LbUsFourByte		blueIndex;
	LbUsFourByte		pinkIndex;
	LbUsOneByte			numBlueAccepted = 0;
	LbUsOneByte			numPinkAccepted = 0;
	LbEightByte			firstPos;
	LbEightByte			secondPos;
	PHG_DetectedPhoton	detectedPhoton;
	
	do { /* Process Loop */

		/* Write out the binary flag and record */
		if ((numBluePhotons != 0) || (numPinkPhotons != 0)) {
			if (hdrHkTyPtr->doCustom == false) {

				if (fwrite(&PhoHFileDecayFlag, sizeof(LbUsOneByte), 1, hdrHkTyPtr->histFile) != 1) {
		
					ErStGeneric("Unable to write decay flag to history binary file.");
					goto FAILURE;
				}

				if (fwrite(decay, sizeof(PHG_Decay), 1, hdrHkTyPtr->histFile) != 1) {
		
					ErStGeneric("Unable to write decay to history binary file.");
					goto FAILURE;
				}							
						
				/* Increment count of outputed decays */
				hdrHkTyPtr->header.H.NumDecays += 1;
			}	
		}
	
		if (hdrHkTyPtr->doCustom == false) {

			/* Process coincidences */
			for (blueIndex = 0; blueIndex < numBluePhotons; blueIndex++) {
				
				/* Incremement Count */
				hdrHkTyPtr->bluesReceived++;
			
				/* Write out the binary flag and record */
				if (fwrite(&PhoHFilePhotonFlag, sizeof(LbUsOneByte), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write photon flag to history binary file.");
					goto FAILURE;
				}
				
				/* Convert to detected photon */
				phoHFileCreateDetectedPhoton(&(bluePhotons[blueIndex]),
					&detectedPhoton);

				if (fwrite(&detectedPhoton, sizeof(PHG_DetectedPhoton), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write blue photon to history binary file.");
					goto FAILURE;
				}
				
				hdrHkTyPtr->bluesAccepted++;
			}
		}	/* Do custom header if there is at least one photon to write */
		else if ((numBluePhotons != 0) || (numPinkPhotons != 0)) {
			
			/* Clear counter */
			numBlueAccepted = 0;

			/* Save current location */
			if ((firstPos = ftell(hdrHkTyPtr->histFile)) == -1) {
				ErStFileError("Unable to get current file position when writing photons fields (PhoHFileWriteDetections)");
				goto FAILURE;
			}
			
			/* Seek past count position */
			if (fseek(hdrHkTyPtr->histFile, sizeof(LbUsOneByte), SEEK_CUR) != 0) {
				ErStFileError("Unable to seek while writing photons (PhoHFileWriteDetections)");
				goto FAILURE;
			}
			
			/* Process coincidences */
			for (blueIndex = 0; blueIndex < numBluePhotons; blueIndex++) {
				
				/* Incremement Count */
				hdrHkTyPtr->bluesReceived++;
				
				if (phoHFileWriteFields(hdrHkTyPtr, decay, &bluePhotons[blueIndex],
						&numBlueAccepted) == false) {
						
					goto FAILURE;
				}
				hdrHkTyPtr->bluesAccepted += numBlueAccepted;
			}
			
			/* Save current location */
			if ((secondPos = ftell(hdrHkTyPtr->histFile)) == -1) {
				ErStFileError("Unable to get current file position when writing photons fields (PhoHFileWriteDetections)");
				goto FAILURE;
			}
			/* Seek back to count position */
			if (fseek(hdrHkTyPtr->histFile, firstPos, SEEK_SET) != 0) {
				ErStFileError("Unable to seek to count position while writing photons (PhoHFileWriteDetections)");
				goto FAILURE;
			}
			
			/* Write the number of blue photons */
			if (fwrite(&numBlueAccepted, sizeof(LbUsOneByte), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to write number of blues accepted to history binary file (PhoHFileWriteDetections).");
				goto FAILURE;
			}
			
			/* Seek past blue photons */
			if (fseek(hdrHkTyPtr->histFile, secondPos, SEEK_SET) != 0) {
				ErStFileError("Unable to seek while writing photons (PhoHFileWriteDetections)");
				goto FAILURE;
			}			
		}
		
		/* Now process pink photons */
		if (hdrHkTyPtr->doCustom == false) {
					
			for (pinkIndex = 0; pinkIndex < numPinkPhotons; pinkIndex++) {

				/* Incremement Count */
				hdrHkTyPtr->pinksReceived++;
				
				/* Write out the binary flag and record */
				if (fwrite(&PhoHFilePhotonFlag, sizeof(LbUsOneByte), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write photon flag to history binary file.");
					goto FAILURE;
				}

				/* Convert to detected photon */
				phoHFileCreateDetectedPhoton(&(pinkPhotons[pinkIndex]),
					&detectedPhoton);	
				
				if (fwrite(&detectedPhoton, sizeof(PHG_DetectedPhoton), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write pink photon to history binary file.");
					goto FAILURE;
				}
				
				hdrHkTyPtr->pinksAccepted++;
			}
		}
		/* Custom history file bug fix */
		else if (PHG_IsPET() & ((numBluePhotons != 0) || (numPinkPhotons != 0))) {

		       /* old code replaced by above line */
		       /* else if (PHG_IsPET()) */
			
			/* Clear counter */	
			numPinkAccepted = 0;

			if ((firstPos = ftell(hdrHkTyPtr->histFile)) == -1) {
				ErStFileError("Unable to get current file position when writing photons fields (PhoHFileWriteDetections)");
				goto FAILURE;
			}
			
			/* Seek past count position */
			if (fseek(hdrHkTyPtr->histFile, sizeof(LbUsOneByte), SEEK_CUR) != 0) {
				ErStFileError("Unable to seek while writing photons (PhoHFileWriteDetections)");
				goto FAILURE;
			}
			
			for (pinkIndex = 0; pinkIndex < numPinkPhotons; pinkIndex++) {
			
				/* Incremement Count */
				hdrHkTyPtr->pinksReceived++;

				if (phoHFileWriteFields(hdrHkTyPtr, decay, &pinkPhotons[pinkIndex],
						&numPinkAccepted) == false) {

					goto FAILURE;
				}
				hdrHkTyPtr->pinksAccepted += numPinkAccepted;
			}
			
			/* Save current location */
			if ((secondPos = ftell(hdrHkTyPtr->histFile)) == -1) {
				ErStFileError("Unable to get current file position when writing photons fields (PhoHFileWriteDetections)");
				goto FAILURE;
			}
			
			/* Seek back to count position */
			if (fseek(hdrHkTyPtr->histFile, firstPos, SEEK_SET) != 0) {
				ErStFileError("Unable to seek to count position while writing photons (PhoHFileWriteDetections)");
				goto FAILURE;
			}
			
			/* Write the number of pink photons */
			if (fwrite(&numPinkAccepted, sizeof(LbUsOneByte), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to write number of pinks accepted to history binary file (PhoHFileWriteDetections).");
				goto FAILURE;
			}
			
			/* Seek past pink photons */
			if (fseek(hdrHkTyPtr->histFile, secondPos, SEEK_SET) != 0) {
				ErStFileError("Unable to seek while writing photons (PhoHFileWriteDetections)");
				goto FAILURE;
			}			
		}
		
		/* Increment count of outputed photons */
		if (hdrHkTyPtr->doCustom == false) {
			hdrHkTyPtr->header.H.NumPhotons += (numBluePhotons + numPinkPhotons);
		}
		else {
			hdrHkTyPtr->header.H.NumPhotons += (numBlueAccepted + numPinkAccepted);
		}

		okay = true;
		FAILURE:;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*	Name:			PhoHFileRewriteDetections
*
*	Summary:		Write the decay and photons.
*					Match the blue/pink scheme.
*					Similar to PhoHFileWriteDetections except has 
*					PHG_DetectedPhoton for arguments rather than
*					PHG_TrackingPhoton.  This allows history decays
*					to be copied more easily from file to file.
*
*	Arguments:
*		PhoHFileHkTy 		*hdrHkTyPtr		- The history hook
*		PHG_Decay			decay			- The decay that started the process.
*		PHG_DetectedPhoton *bluePhotons		- The blue photons detected.
*		LbUsFourByte 		numBluePhotons	- The number of blue photons.
*		PHG_DetectedPhoton *pinkPhotons		- The  pink photons detected.
*		LbUsFourByte		numPinkPhotons	- The number of pink photons.
*
*	Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean PhoHFileRewriteDetections(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay,
			PHG_DetectedPhoton *bluePhotons, LbUsFourByte numBluePhotons,
			PHG_DetectedPhoton *pinkPhotons, LbUsFourByte numPinkPhotons)
{
	Boolean				okay = false;
	LbUsFourByte		blueIndex;
	LbUsFourByte		pinkIndex;
	LbUsOneByte			numBlueAccepted = 0;
	LbUsOneByte			numPinkAccepted = 0;
	
	
	do { /* Process Loop */
		
		/* Write out the binary flag and record */
		if ((numBluePhotons != 0) || (numPinkPhotons != 0)) {
			if (fwrite(&PhoHFileDecayFlag, sizeof(LbUsOneByte), 1, hdrHkTyPtr->histFile) != 1) {
	
				ErStGeneric("Unable to write decay flag to history binary file.");
				goto FAILURE;
			}

			if (fwrite(decay, sizeof(PHG_Decay), 1, hdrHkTyPtr->histFile) != 1) {
	
				ErStGeneric("Unable to write decay to history binary file.");
				goto FAILURE;
			}							
					
			/* Increment count of outputed decays */
			hdrHkTyPtr->header.H.NumDecays += 1;
		}
		
		/* Process coincidences */
		for (blueIndex = 0; blueIndex < numBluePhotons; blueIndex++) {
			
			/* Incremement Count */
			hdrHkTyPtr->bluesReceived++;
		
			/* Write out the binary flag and record */
			if (fwrite(&PhoHFilePhotonFlag, sizeof(LbUsOneByte), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to write photon flag to history binary file.");
				goto FAILURE;
			}
			
			/* Write out the detected photon */
			if (fwrite(&(bluePhotons[blueIndex]), sizeof(PHG_DetectedPhoton), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to write blue photon to history binary file.");
				goto FAILURE;
			}
			
			hdrHkTyPtr->bluesAccepted++;
		}
		
		/* Now process pink photons */
		for (pinkIndex = 0; pinkIndex < numPinkPhotons; pinkIndex++) {

			/* Incremement Count */
			hdrHkTyPtr->pinksReceived++;
			
			/* Write out the binary flag and record */
			if (fwrite(&PhoHFilePhotonFlag, sizeof(LbUsOneByte), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to write photon flag to history binary file.");
				goto FAILURE;
			}

			/* Write out the detected photon */
			if (fwrite(&(pinkPhotons[pinkIndex]), sizeof(PHG_DetectedPhoton), 1, hdrHkTyPtr->histFile) != 1) {
				
				ErStGeneric("Unable to write pink photon to history binary file.");
				goto FAILURE;
			}
			
			hdrHkTyPtr->pinksAccepted++;
		}
		
		/* Increment count of outputed photons */
		if (hdrHkTyPtr->doCustom == false) {
			hdrHkTyPtr->header.H.NumPhotons += (numBluePhotons + numPinkPhotons);
		}
		else {
			hdrHkTyPtr->header.H.NumPhotons += (numBlueAccepted + numPinkAccepted);
		}

		okay = true;
		FAILURE:;
	} while (false);
	
	return (okay);
}


/*********************************************************************************
*
*			Name:			phoHFileWriteFields
*
*			Summary:		Write the detected decay and photons.
*
*			Arguments:
*				PhoHFileHkTy 		*hdrHkTyPtr		- The history hook
*				PHG_Decay			*decay			- The decay that started the process.
*				PHG_TrackingPhoton 	*photon			- The blue photons detected.
*				LbUsOneByte			*numAccepted	- Count of accepted (written)
*			Function return: None.
*
*********************************************************************************/
Boolean phoHFileWriteFields(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay,
			PHG_TrackingPhoton *photon, LbUsOneByte *numAccepted)
{
	Boolean			okay = false;
	LbUsFourByte	i, size;
	
	do { /* Process Loop */
	
		#ifdef PHG_DEBUG_INTENSELY
			PhoHFilePrintFields(hdrHkTyPtr, decay, photon);
		#endif
		
		/* First, filter acceptable photons */
		
		if (hdrHkTyPtr->customParams.doXrange == true) {
			if ((photon->location.x_position < hdrHkTyPtr->customParams.xPosMin) ||
					(photon->location.x_position > hdrHkTyPtr->customParams.xPosMax)) {

				goto FILTERED_OUT;

			}						
		}

		if (hdrHkTyPtr->customParams.doYrange == true) {
			if ((photon->location.y_position < hdrHkTyPtr->customParams.yPosMin) ||
					(photon->location.y_position > hdrHkTyPtr->customParams.yPosMax)) {

				goto FILTERED_OUT;

			}							
		}

		if (hdrHkTyPtr->customParams.doZrange == true) {
			if ((photon->location.z_position < hdrHkTyPtr->customParams.zPosMin) ||
					(photon->location.z_position > hdrHkTyPtr->customParams.zPosMax)) {

				goto FILTERED_OUT;

			}							
		}

		if (hdrHkTyPtr->customParams.doXCosRange == true) {
			if ((photon->angle.cosine_x < hdrHkTyPtr->customParams.xCosMin) ||
					(photon->angle.cosine_x > hdrHkTyPtr->customParams.xCosMax)) {

				goto FILTERED_OUT;

			}						
		}

		if (hdrHkTyPtr->customParams.doYCosRange == true) {
			if ((photon->angle.cosine_y < hdrHkTyPtr->customParams.yCosMin) ||
					(photon->angle.cosine_y > hdrHkTyPtr->customParams.yCosMax)) {

				goto FILTERED_OUT;

			}							
		}

		if (hdrHkTyPtr->customParams.doZCosRange == true) {
			if ((photon->angle.cosine_z < hdrHkTyPtr->customParams.zCosMin) ||
					(photon->angle.cosine_z > hdrHkTyPtr->customParams.zCosMax)) {

				goto FILTERED_OUT;

			}							
		}

		if (hdrHkTyPtr->customParams.doObjScatRange == true) {
			if (((LbFourByte)photon->num_of_scatters < hdrHkTyPtr->customParams.objScattersMin) ||
					((LbFourByte)photon->num_of_scatters > hdrHkTyPtr->customParams.objScattersMax)) {

				goto FILTERED_OUT;

			}						
		}

		if (hdrHkTyPtr->customParams.doColScatRange == true) {
			if (((LbFourByte)photon->scatters_in_col < hdrHkTyPtr->customParams.colScattersMin) ||
					((LbFourByte)photon->scatters_in_col > hdrHkTyPtr->customParams.colScattersMax)) {

				goto FILTERED_OUT;

			}						
		}

		if (hdrHkTyPtr->customParams.doEnergyRange == true) {
			if ((photon->energy < hdrHkTyPtr->customParams.energyMin) ||
					(photon->energy > hdrHkTyPtr->customParams.energyMax)) {

				goto FILTERED_OUT;

			}						
		}

		if (hdrHkTyPtr->customParams.doTravDistRange == true) {
			if ((photon->travel_distance < hdrHkTyPtr->customParams.travelDistanceMin) ||
					(photon->travel_distance > hdrHkTyPtr->customParams.travelDistanceMax)) {

				goto FILTERED_OUT;

			}
		}

		if (hdrHkTyPtr->customParams.doDecayXrange == true) {
			if ((decay->location.x_position < hdrHkTyPtr->customParams.decayXPosMin) ||
					(decay->location.x_position > hdrHkTyPtr->customParams.decayXPosMax)) {

				goto FILTERED_OUT;

			}
		}

		if (hdrHkTyPtr->customParams.doDecayYrange == true) {
			if ((decay->location.y_position < hdrHkTyPtr->customParams.decayYPosMin) ||
					(decay->location.y_position > hdrHkTyPtr->customParams.decayYPosMax)) {

				goto FILTERED_OUT;

			}
		}

		if (hdrHkTyPtr->customParams.doDecayZrange == true) {
			if ((decay->location.z_position < hdrHkTyPtr->customParams.decayZPosMin) ||
					(decay->location.z_position > hdrHkTyPtr->customParams.decayZPosMax)) {

				goto FILTERED_OUT;

			}
		}

		if (hdrHkTyPtr->customParams.doTransDistRange == true) {
			if ((photon->transaxialPosition < hdrHkTyPtr->customParams.transaxialDistanceMin) ||
					(photon->transaxialPosition > hdrHkTyPtr->customParams.transaxialDistanceMax)) {

				goto FILTERED_OUT;

			}
		}

		if (hdrHkTyPtr->customParams.doAziAngRange == true) {
			if ((photon->azimuthalAngleIndex < hdrHkTyPtr->customParams.aaIndexMin) ||
					(photon->azimuthalAngleIndex > hdrHkTyPtr->customParams.aaIndexMax)) {

				goto FILTERED_OUT;

			}
		}

		if (hdrHkTyPtr->customParams.doAxiPosRange == true) {
			if ((photon->axialPosition < hdrHkTyPtr->customParams.axialPositionMin) ||
					(photon->axialPosition > hdrHkTyPtr->customParams.axialPositionMax)) {

				goto FILTERED_OUT;

			}
		}

		if (hdrHkTyPtr->customParams.doDetectorXposRange == true) {
			if ((photon->detLocation.x_position < hdrHkTyPtr->customParams.detectorXpositionMin) ||
					(photon->detLocation.x_position > hdrHkTyPtr->customParams.detectorXpositionMax)) {

				goto FILTERED_OUT;

			}
		}

		if (hdrHkTyPtr->customParams.doDetectorYposRange == true) {
			if ((photon->detLocation.y_position < hdrHkTyPtr->customParams.detectorYpositionMin) ||
					(photon->detLocation.y_position > hdrHkTyPtr->customParams.detectorYpositionMax)) {

				goto FILTERED_OUT;

			}
		}

		if (hdrHkTyPtr->customParams.doDetectorZposRange == true) {
			if ((photon->detLocation.z_position < hdrHkTyPtr->customParams.detectorZpositionMin) ||
					(photon->detLocation.z_position > hdrHkTyPtr->customParams.detectorZpositionMax)) {

				goto FILTERED_OUT;

			}
		}

		if (hdrHkTyPtr->customParams.doDetectorAngleRange == true) {
			if ((photon->detectorAngle < hdrHkTyPtr->customParams.detectorAngleMin) ||
					(photon->detectorAngle > hdrHkTyPtr->customParams.detectorAngleMax)) {

				goto FILTERED_OUT;

			}
		}

		/* Now do the actual writing of fields */

		if (hdrHkTyPtr->customParams.doXposition == true) {
				if (fwrite(&photon->location.x_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'x_position' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doYposition == true) {
				if (fwrite(&photon->location.y_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'y_position' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doZposition == true) {
				if (fwrite(&photon->location.z_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'z_position' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doXcosine == true) {
				if (fwrite(&photon->angle.cosine_x, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'cosine_x' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doYcosine == true) {
				if (fwrite(&photon->angle.cosine_y, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'cosine_y' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doZcosine == true) {
				if (fwrite(&photon->angle.cosine_z, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'cosine_z' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doScattersInObject == true) {
				if (fwrite(&photon->num_of_scatters, sizeof(LbUsFourByte), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'scatters in object' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doScattersInCollimator == true) {
				if (fwrite(&photon->scatters_in_col, sizeof(LbUsFourByte), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'scatters in collimator' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doDecayWeight == true) {
			if (fwrite(&decay->startWeight, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to write 'decay weight' to history file.");
				goto FAILURE;
			}
		}

		if (hdrHkTyPtr->customParams.doWeight == true) {
			if (fwrite(&photon->photon_current_weight, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to write 'photon weight' to history file.");
				goto FAILURE;
			}
		}

		if (hdrHkTyPtr->customParams.doEnergy == true) {
				if (fwrite(&photon->energy, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'energy' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doTravelDistance == true) {
				if (fwrite(&photon->travel_distance, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'travel distance' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doDecayXposition == true) {
				if (fwrite(&decay->location.x_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'decay x position' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doDecayYposition == true) {
				if (fwrite(&decay->location.y_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'decay y position' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doDecayZposition == true) {
				if (fwrite(&decay->location.z_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'decay z position' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doDecayTime == true) {
				if (fwrite(&decay->decayTime, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'decay time' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doDecayType == true) {
				if (fwrite(&decay->decayType, sizeof(LbUsOneByte), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'decay type' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doTransaxialDistance == true) {
				if (fwrite(&photon->transaxialPosition, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'transaxial position' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doAzimuthalAngleIndex == true) {
				if (fwrite(&photon->azimuthalAngleIndex, sizeof(LbUsFourByte), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'azimuthal angle index' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doAxialPosition == true) {
				if (fwrite(&photon->axialPosition, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'axial position' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doDetectorXposition == true) {
				if (fwrite(&photon->detLocation.x_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'detector X-axis position' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doDetectorYposition == true) {
				if (fwrite(&photon->detLocation.y_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'detector Y-axis position' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doDetectorZposition == true) {
				if (fwrite(&photon->detLocation.z_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'detector Z-axis position' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doDetectorAngle == true) {
				if (fwrite(&photon->detectorAngle, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'detector angle' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doDetectorCrystal == true) {
				if (fwrite(&photon->detCrystal, sizeof(LbFourByte), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'detector crystal' to history file.");
					goto FAILURE;
				}
		}

		if (hdrHkTyPtr->customParams.doNumDetectorInteractions == true) {
			if (fwrite(&photon->num_det_interactions, sizeof(LbUsFourByte), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to write 'number of detector interactions' to history file.");
				goto FAILURE;
			}
		}

		if (hdrHkTyPtr->customParams.doDetInteractionPos == true) {
			
			if (fwrite(&photon->num_det_interactions, sizeof(LbUsFourByte), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to write 'number of detector interactions' to history file for list of interaction positions.");
				goto FAILURE;
			}
			for (i = 0; i < photon->num_det_interactions; i++) {

				if (fwrite(&photon->det_interactions[i].pos.x_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'x-position of interaction' to history file for list of interaction positions.");
					goto FAILURE;
				}

				if (fwrite(&photon->det_interactions[i].pos.y_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'y-position of interaction' to history file for list of interaction positions.");
					goto FAILURE;
				}

				if (fwrite(&photon->det_interactions[i].pos.z_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'z-position of interaction' to history file for list of interaction positions.");
					goto FAILURE;
				}

				if (fwrite(&photon->det_interactions[i].energy_deposited, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'energy of interaction' to history file for list of interaction positions.");
					goto FAILURE;
				}
				size = sizeof(Boolean);
				
				if (fwrite(&photon->det_interactions[i].isActive, size, 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to write 'is-Active of interaction' to history file for list of interaction positions.");
					goto FAILURE;
				}
			}
		}
		
		(*numAccepted)++;
		
		FILTERED_OUT:;
		okay = true;
		FAILURE:;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			PhoHFileReadFields
*
*			Summary:		Read photon fields from a custom history file.
*
*			Arguments:
*				PhoHFileHkTy 		*hdrHkTyPtr		- The history hook.
*				PHG_Decay			*decay			- The decay that started the process.
*				PHG_TrackingPhoton 	*photon			- The blue photons detected.
*				Boolean				*photonAccepted	- Flag indicating the photon was accepted.
*			Function return: None.
*
*********************************************************************************/
Boolean PhoHFileReadFields(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay,
			PHG_TrackingPhoton *photon, Boolean *photonAccepted)
{
	Boolean			okay = false;
	LbUsFourByte	i;
	
	do { /* Process Loop */


		/* Clear the acceptance flag */
		*photonAccepted = false;
		
		/* Read the fields */

		if (hdrHkTyPtr->customParams.doXposition == true) {

			if (fread(&photon->location.x_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to read 'x_position' from history file.");
				goto FAILURE;
			}
			
			/* Check range */
			if (hdrHkTyPtr->customParams.doXrange == true) {
				if ((photon->location.x_position < hdrHkTyPtr->customParams.xPosMin) ||
						(photon->location.x_position > hdrHkTyPtr->customParams.xPosMax)) {

					goto FILTERED_OUT;

				}						
			}
		}

		if (hdrHkTyPtr->customParams.doYposition == true) {
				if (fread(&photon->location.y_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'y_position' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doYrange == true) {
				if ((photon->location.y_position < hdrHkTyPtr->customParams.yPosMin) ||
						(photon->location.y_position > hdrHkTyPtr->customParams.yPosMax)) {

					goto FILTERED_OUT;

				}							
			}
		}

		if (hdrHkTyPtr->customParams.doZposition == true) {
				if (fread(&photon->location.z_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'z_position' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doZrange == true) {
				if ((photon->location.z_position < hdrHkTyPtr->customParams.zPosMin) ||
						(photon->location.z_position > hdrHkTyPtr->customParams.zPosMax)) {

					goto FILTERED_OUT;

				}							
			}
		}

		if (hdrHkTyPtr->customParams.doXcosine == true) {
				if (fread(&photon->angle.cosine_x, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'cosine_x' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doXCosRange == true) {
				if ((photon->angle.cosine_x < hdrHkTyPtr->customParams.xCosMin) ||
						(photon->angle.cosine_x > hdrHkTyPtr->customParams.xCosMax)) {

					goto FILTERED_OUT;

				}						
			}
		}

		if (hdrHkTyPtr->customParams.doYcosine == true) {
				if (fread(&photon->angle.cosine_y, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'cosine_y' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doYCosRange == true) {
				if ((photon->angle.cosine_y < hdrHkTyPtr->customParams.yCosMin) ||
						(photon->angle.cosine_y> hdrHkTyPtr->customParams.yCosMax)) {

					goto FILTERED_OUT;

				}							
			}
		}

		if (hdrHkTyPtr->customParams.doZcosine == true) {
				if (fread(&photon->angle.cosine_z, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'cosine_z' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doZCosRange == true) {
				if ((photon->angle.cosine_z < hdrHkTyPtr->customParams.zCosMin) ||
						(photon->angle.cosine_z > hdrHkTyPtr->customParams.zCosMax)) {

					goto FILTERED_OUT;

				}							
			}
		}

		if (hdrHkTyPtr->customParams.doScattersInObject == true) {
				if (fread(&photon->num_of_scatters, sizeof(LbUsFourByte), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'scatters in object' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doObjScatRange == true) {
				if (((LbFourByte)photon->num_of_scatters < hdrHkTyPtr->customParams.objScattersMin) ||
						((LbFourByte)photon->num_of_scatters > hdrHkTyPtr->customParams.objScattersMax)) {

					goto FILTERED_OUT;

				}						
			}
		}

		if (hdrHkTyPtr->customParams.doScattersInCollimator == true) {
				if (fread(&photon->scatters_in_col, sizeof(LbUsFourByte), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'scatters in collimator' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doColScatRange == true) {
				if (((LbFourByte)photon->scatters_in_col < hdrHkTyPtr->customParams.colScattersMin) ||
						((LbFourByte)photon->scatters_in_col > hdrHkTyPtr->customParams.colScattersMax)) {

					goto FILTERED_OUT;

				}						
			}
		}

		if (hdrHkTyPtr->customParams.doDecayWeight == true) {
			if (fread(&decay->startWeight, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to read 'decay weight' from history file.");
				goto FAILURE;
			}
		}

		if (hdrHkTyPtr->customParams.doWeight == true) {
			/* read in the current weight field */
			if (fread(&photon->photon_current_weight, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to read 'photon weight' from history file.");
				goto FAILURE;
			}
			/* copy to the scatter or primary weight field */
			if ((photon->num_of_scatters > 0)) {
				photon->photon_scatter_weight = photon->photon_current_weight;
			}
			else {
				photon->photon_primary_weight = photon->photon_current_weight;
			}						
		}

		if (hdrHkTyPtr->customParams.doEnergy == true) {
				if (fread(&photon->energy, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'energy' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doEnergyRange == true) {
				if ((photon->energy < hdrHkTyPtr->customParams.energyMin) &&
						(photon->energy > hdrHkTyPtr->customParams.energyMax)) {

					goto FILTERED_OUT;

				}						
			}
		}

		if (hdrHkTyPtr->customParams.doTravelDistance == true) {
				if (fread(&photon->travel_distance, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'travel distance' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doTravDistRange == true) {
				if ((photon->travel_distance < hdrHkTyPtr->customParams.travelDistanceMin) ||
						(photon->travel_distance > hdrHkTyPtr->customParams.travelDistanceMax)) {

					goto FILTERED_OUT;

				}
			}
		}

		if (hdrHkTyPtr->customParams.doDecayXposition == true) {
				if (fread(&decay->location.x_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'decay x position' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doDecayXrange == true) {
				if ((decay->location.x_position < hdrHkTyPtr->customParams.decayXPosMin) ||
						(decay->location.x_position > hdrHkTyPtr->customParams.decayXPosMax)) {

					goto FILTERED_OUT;

				}
			}
		}

		if (hdrHkTyPtr->customParams.doDecayYposition == true) {
				if (fread(&decay->location.y_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'decay y position' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doDecayYrange == true) {
				if ((decay->location.y_position < hdrHkTyPtr->customParams.decayYPosMin) ||
						(decay->location.y_position > hdrHkTyPtr->customParams.decayYPosMax)) {

					goto FILTERED_OUT;

				}
			}
		}

		if (hdrHkTyPtr->customParams.doDecayZposition == true) {
				if (fread(&decay->location.z_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'decay z position' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doDecayZrange == true) {
				if ((decay->location.z_position < hdrHkTyPtr->customParams.decayZPosMin) ||
						(decay->location.z_position > hdrHkTyPtr->customParams.decayZPosMax)) {

					goto FILTERED_OUT;

				}
			}
		}

		if (hdrHkTyPtr->customParams.doDecayTime == true) {
			/* read in the current weight field */
			if (fread(&decay->decayTime, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to read 'decay time' from history file.");
				goto FAILURE;
			}
		}

		if (hdrHkTyPtr->customParams.doDecayType == true) {
			/* read in the current weight field */
			if (fread(&decay->decayType, sizeof(LbUsOneByte), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to read 'decay type' from history file.");
				goto FAILURE;
			}
		}

		if (hdrHkTyPtr->customParams.doTransaxialDistance == true) {
				if (fread(&photon->transaxialPosition, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'transaxial position' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doTransDistRange == true) {
				if ((photon->transaxialPosition < hdrHkTyPtr->customParams.transaxialDistanceMin) ||
						(photon->transaxialPosition > hdrHkTyPtr->customParams.transaxialDistanceMax)) {

					goto FILTERED_OUT;

				}
			}
		}

		if (hdrHkTyPtr->customParams.doAzimuthalAngleIndex == true) {
				if (fread(&photon->azimuthalAngleIndex, sizeof(LbUsFourByte), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'azimuthal angle index' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doAziAngRange == true) {
				if ((photon->azimuthalAngleIndex < hdrHkTyPtr->customParams.aaIndexMin) ||
						(photon->azimuthalAngleIndex > hdrHkTyPtr->customParams.aaIndexMax)) {

					goto FILTERED_OUT;

				}
			}
		}

		if (hdrHkTyPtr->customParams.doAxialPosition == true) {

			if (fread(&photon->axialPosition, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to read 'axial position' from history file.");
				goto FAILURE;
			}

			if (hdrHkTyPtr->customParams.doAxiPosRange == true) {
				if ((photon->axialPosition < hdrHkTyPtr->customParams.axialPositionMin) ||
						(photon->axialPosition > hdrHkTyPtr->customParams.axialPositionMax)) {

					goto FILTERED_OUT;

				}
			}
		}

		if (hdrHkTyPtr->customParams.doDetectorXposition == true) {
			if (fread(&photon->detLocation.x_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to read 'detector X-axis position' from history file.");
				goto FAILURE;
			}

			if (hdrHkTyPtr->customParams.doDetectorXposRange == true) {
				if ((photon->detLocation.x_position < hdrHkTyPtr->customParams.detectorXpositionMin) ||
						(photon->detLocation.x_position > hdrHkTyPtr->customParams.detectorXpositionMax)) {

					goto FILTERED_OUT;

				}
			}
		}

		if (hdrHkTyPtr->customParams.doDetectorYposition == true) {
			if (fread(&photon->detLocation.y_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to read 'detector Y-axis position' from history file.");
				goto FAILURE;
			}

			if (hdrHkTyPtr->customParams.doDetectorYposRange == true) {
				if ((photon->detLocation.y_position < hdrHkTyPtr->customParams.detectorYpositionMin) ||
						(photon->detLocation.y_position > hdrHkTyPtr->customParams.detectorYpositionMax)) {

					goto FILTERED_OUT;

				}
			}
		}

		if (hdrHkTyPtr->customParams.doDetectorZposition == true) {
				if (fread(&photon->detLocation.z_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'detector Z-axis position' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doDetectorZposRange == true) {
				if ((photon->detLocation.z_position < hdrHkTyPtr->customParams.detectorZpositionMin) ||
						(photon->detLocation.z_position > hdrHkTyPtr->customParams.detectorZpositionMax)) {

					goto FILTERED_OUT;

				}
			}
		}

		if (hdrHkTyPtr->customParams.doDetectorAngle == true) {
				if (fread(&photon->detectorAngle, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'detector angle' from history file.");
					goto FAILURE;
				}

			if (hdrHkTyPtr->customParams.doDetectorAngleRange == true) {
				if ((photon->detectorAngle < hdrHkTyPtr->customParams.detectorAngleMin) ||
						(photon->detectorAngle > hdrHkTyPtr->customParams.detectorAngleMax)) {

					goto FILTERED_OUT;

				}
			}
		}

		if (hdrHkTyPtr->customParams.doDetectorCrystal == true) {
				if (fread(&photon->detCrystal, sizeof(LbFourByte), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'detector crystal' from history file.");
					goto FAILURE;
				}

		}

		if (hdrHkTyPtr->customParams.doNumDetectorInteractions == true) {
			if (fread(&photon->num_det_interactions, sizeof(LbUsFourByte), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to read 'number of detector interactions' from history file.");
				goto FAILURE;
			}

			if ((photon->detLocation.z_position < hdrHkTyPtr->customParams.detectorZpositionMin) ||
					(photon->detLocation.z_position > hdrHkTyPtr->customParams.detectorZpositionMax)) {

				goto FILTERED_OUT;

			}
		}

		if (hdrHkTyPtr->customParams.doDetInteractionPos == true) {
			if (fread(&photon->num_det_interactions, sizeof(LbUsFourByte), 1, hdrHkTyPtr->histFile) != 1) {

				ErStGeneric("Unable to read 'number of detector interactions' from history file for list of interaction positions.");
				goto FAILURE;
			}
			for (i = 0; i < photon->num_det_interactions; i++) {

				if (fread(&photon->det_interactions[i].pos.x_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'x-position of interaction' from history file for list of interaction positions.");
					goto FAILURE;
				}

				if (fread(&photon->det_interactions[i].pos.y_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'y-position of interaction' from history file for list of interaction positions.");
					goto FAILURE;
				}

				if (fread(&photon->det_interactions[i].pos.z_position, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'z-position of interaction' from history file for list of interaction positions.");
					goto FAILURE;
				}

				if (fread(&photon->det_interactions[i].energy_deposited, sizeof(double), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'energy of interaction' from history file for list of interaction positions.");
					goto FAILURE;
				}

				if (fread(&photon->det_interactions[i].isActive, sizeof(Boolean), 1, hdrHkTyPtr->histFile) != 1) {

					ErStGeneric("Unable to read 'is-Active of interaction' from history file for list of interaction positions.");
					goto FAILURE;
				}
			}
		}
		
		
		*photonAccepted = true;
		FILTERED_OUT:;
		okay = true;
		FAILURE:;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			PhoHFilePrintFields
*
*			Summary:		Print the fields of photons.
*
*			Arguments:
*				PhoHFileHkTy 		*hdrHkTyPtr		- Header information 
*				PHG_Decay			*decay			- The decay that started the process.
*				PHG_TrackingPhoton 	*photon			- The blue photons detected.
*			Function return: None.
*
*********************************************************************************/
void PhoHFilePrintFields(PhoHFileHkTy *hdrHkTyPtr, PHG_Decay *decay, PHG_TrackingPhoton *photon)
{
	

		if (hdrHkTyPtr->customParams.doXposition == true) {
				LbInPrintf("\nx position = %3.3e", photon->location.x_position);
		}

		if (hdrHkTyPtr->customParams.doYposition == true) {
				LbInPrintf("\ny position = %3.3e", photon->location.y_position);
		}

		if (hdrHkTyPtr->customParams.doZposition == true) {
				LbInPrintf("\nz position = %3.3e", photon->location.z_position);
		}

		if (hdrHkTyPtr->customParams.doXcosine == true) {
				LbInPrintf("\ncosine x = %3.3e", photon->angle.cosine_x);
		}

		if (hdrHkTyPtr->customParams.doYcosine == true) {
				LbInPrintf("\ncosine y = %3.3e", photon->angle.cosine_y);
		}

		if (hdrHkTyPtr->customParams.doZcosine == true) {
				LbInPrintf("\ncosine z = %3.3e", photon->angle.cosine_z);
		}

		if (hdrHkTyPtr->customParams.doScattersInObject == true) {
				LbInPrintf("\nnumber of scatters = %d", photon->num_of_scatters);
		}

		if (hdrHkTyPtr->customParams.doScattersInCollimator == true) {
				LbInPrintf("\nnumber of scatters in collimator = %d", photon->scatters_in_col);
		}

		if (hdrHkTyPtr->customParams.doDecayWeight == true) {
			LbInPrintf("\ndecay weight = %3.3e", decay->startWeight);
		}

		if (hdrHkTyPtr->customParams.doWeight == true) {
			LbInPrintf("\nphoton weight = %3.3e", photon->photon_current_weight);
		}

		if (hdrHkTyPtr->customParams.doEnergy == true) {
				LbInPrintf("\nenergy = %3.1f", photon->energy);
		}

		if (hdrHkTyPtr->customParams.doTravelDistance == true) {
				LbInPrintf("\ntravel distance = %3.3e", photon->travel_distance);
		}

		if (hdrHkTyPtr->customParams.doDecayXposition == true) {
				LbInPrintf("\n decay x position = %3.3e", decay->location.x_position);
		}

		if (hdrHkTyPtr->customParams.doDecayYposition == true) {
				LbInPrintf("\n decay y position = %3.3e", decay->location.y_position);
		}

		if (hdrHkTyPtr->customParams.doDecayZposition == true) {
				LbInPrintf("\n decay z position = %3.3e", decay->location.z_position);
		}

		if (hdrHkTyPtr->customParams.doDecayTime == true) {
				LbInPrintf("\n decay time = %4.4e", decay->decayTime);
		}

		if (hdrHkTyPtr->customParams.doDecayType == true) {
			if (decay->decayType == PhgEn_SinglePhoton)
				LbInPrintf("\n decay type = single photon");
			else if (decay->decayType == PhgEn_Positron)
				LbInPrintf("\n decay type = positron");
			else if (decay->decayType == PhgEn_PETRandom)
				LbInPrintf("\n decay type = random");
			else if (decay->decayType == PhgEn_Complex)
				LbInPrintf("\n decay type = complex");
			else 
				LbInPrintf("\n WARNING: decay type UNKNOWN.");
		}

		if (hdrHkTyPtr->customParams.doTransaxialDistance == true) {
				LbInPrintf("\ntransaxial distance = %3.3e", photon->transaxialPosition);
		}

		if (hdrHkTyPtr->customParams.doAzimuthalAngleIndex == true) {
				LbInPrintf("\nazimuthal angle = %d", photon->azimuthalAngleIndex);
		}

		if (hdrHkTyPtr->customParams.doAxialPosition == true) {
				LbInPrintf("\naxial position = %3.3e", photon->axialPosition);
		}

		if (hdrHkTyPtr->customParams.doDetectorXposition == true) {
				LbInPrintf("\ndetector x position = %3.3e", photon->detLocation.x_position);
		}

		if (hdrHkTyPtr->customParams.doDetectorYposition == true) {
				LbInPrintf("\ndetector y position = %3.3e", photon->detLocation.y_position);
		}

		if (hdrHkTyPtr->customParams.doDetectorZposition == true) {
				LbInPrintf("\ndetector z position = %3.3e", photon->detLocation.z_position);
		}

		if (hdrHkTyPtr->customParams.doDetectorAngle == true) {
				LbInPrintf("\ndetector angle = %3.3e", photon->detectorAngle);
		}

		if (hdrHkTyPtr->customParams.doDetectorCrystal == true) {
				LbInPrintf("\ndetector crystal = %d", photon->detCrystal);
		}

		if (hdrHkTyPtr->customParams.doNumDetectorInteractions == true) {
			LbInPrintf("\nnumber of detector interactions = %d", photon->num_det_interactions);
		}
		
		LbInPrintf("\n");
}

/*********************************************************************************
*
*			Name:		phoHFileCreateDetectedPhoton
*
*			Summary:	Make a new detected photon (To write to the history file).
*			Arguments:
*				PHG_TrackingPhoton	*trackingphoton	- The tracking photon.
*				PHG_DetectedPhoton	*detectedphoton	- The detected photon.
*				
*			Function return: None.
*
*********************************************************************************/
void phoHFileCreateDetectedPhoton(PHG_TrackingPhoton *trackingphoton,
		PHG_DetectedPhoton *detectedphoton)	
{
	LbUsFourByte		numScatters;	/* sum of the scatters in object and collimator */
	
	/* Convert tracking photon travel distance to "time since creation" */
	detectedphoton->time_since_creation = trackingphoton->travel_distance/PHGMATH_SPEED_OF_LIGHT;

	/* Convert the tracking photon to an detected photon */
	detectedphoton->location.x_position = (float) trackingphoton->location.x_position;
	detectedphoton->location.y_position = (float) trackingphoton->location.y_position;
	detectedphoton->location.z_position = (float) trackingphoton->location.z_position;

	detectedphoton->angle.cosine_x = (float) trackingphoton->angle.cosine_x;
	detectedphoton->angle.cosine_y = (float) trackingphoton->angle.cosine_y;
	detectedphoton->angle.cosine_z = (float) trackingphoton->angle.cosine_z;

	detectedphoton->energy = (float) trackingphoton->energy;
	detectedphoton->transaxialPosition = trackingphoton->transaxialPosition;
	detectedphoton->azimuthalAngleIndex = trackingphoton->azimuthalAngleIndex;
	detectedphoton->detectorAngle = trackingphoton->detectorAngle;
	detectedphoton->detCrystal = trackingphoton->detCrystal;
	
	/* Set the weight */
	detectedphoton->photon_weight = trackingphoton->photon_current_weight;
		
	/* Construct the flags */
	{
		/* Start by setting number of scatters flags */

		/* sum the scatters in the object and collimator--if the user wants to preserve
		the difference between the two, they need to use the custom history file */
		numScatters = trackingphoton->num_of_scatters + trackingphoton->scatters_in_col;
		
		if (numScatters >= MAX_SCATTERS)
			detectedphoton->flags = (MAX_SCATTERS << 2);
		else
			detectedphoton->flags = (numScatters << 2);
		
		/* Now set "blue" flag */
		if (PHG_IsBlue(trackingphoton))
			LbFgSet(detectedphoton->flags, PHGFg_PhotonBlue);
	}
}
					
/*********************************************************************************
*
*			Name:			PhoHFileGetRunTimeParams
*
*			Summary:	Read in the runtime parameters.
*			Arguments:
*
*			Function return: true unless an error occurs.
*
*********************************************************************************/
Boolean PhoHFileGetRunTimeParams(char *paramPath, PhoHFileRunTimeParamsTy *customParams)	
{
	double					paramBuffer[LBPF_PARAM_LEN];
	Boolean					okay = false;					/* Process flag */
	Boolean					isEOF;							/* End of file flag */
	Boolean					switchOkay;						/* Switch flag */
	LbPfEnPfTy				paramType;
	LbUsFourByte			paramSize;
	char					paramLabel[LBPF_LABEL_LEN];
	PhoHFileEn_RunTimeParamsTy	whichParam;
	LbPfHkTy				paramFileHk;
	
	do { /* Process Loop */


		
		/* Attempt to open parameter file */
		if (!LbPfOpen(paramPath, 0, &paramFileHk)) {
			break;
		}
		
		/* Set default boolean values */
		customParams->doXposition = false;
		customParams->doXrange = false;
		customParams->doYposition = false;
		customParams->doYrange = false;
		customParams->doZposition = false;
		customParams->doZrange = false;
		customParams->doXcosine = false;
		customParams->doXCosRange = false;
		customParams->doYcosine = false;
		customParams->doYCosRange = false;
		customParams->doZcosine = false;
		customParams->doZCosRange = false;
		customParams->doScattersInObject = false;
		customParams->doObjScatRange = false;
		customParams->doScattersInCollimator = false;
		customParams->doColScatRange = false;
		customParams->doDecayWeight = false;
		customParams->doEnergy = false;
		customParams->doEnergyRange = false;
		customParams->doTravelDistance = false;
		customParams->doTravDistRange = false;
		customParams->doDecayXposition = false;
		customParams->doDecayXrange = false;
		customParams->doDecayYposition = false;
		customParams->doDecayYrange = false;
		customParams->doDecayZposition = false;
		customParams->doDecayZrange = false;
		customParams->doDecayTime = false;
		customParams->doDecayType = false;
		customParams->doTransaxialDistance = false;
		customParams->doTransDistRange = false;
		customParams->doAzimuthalAngleIndex = false;
		customParams->doAziAngRange = false;
		customParams->doAxialPosition = false;
		customParams->doAxiPosRange = false;
		customParams->doDetectorXposition = false;
		customParams->doDetectorXposRange = false;
		customParams->doDetectorYposition = false;
		customParams->doDetectorYposRange = false;
		customParams->doDetectorZposition = false;
		customParams->doDetectorZposRange = false;
		customParams->doNumDetectorInteractions = false;
		customParams->doNumDetInteractionsRange = false;
		customParams->doDetInteractionPos = false;
		customParams->doDetectorAngle = false;
		customParams->doDetectorAngleRange = false;
		customParams->doDetectorCrystal = false;
		
		customParams->xPosMin = LBDOUBLE_MAX;
		customParams->xPosMax = LBDOUBLE_MIN;
		customParams->yPosMin = LBDOUBLE_MAX;
		customParams->yPosMax = LBDOUBLE_MIN;
		customParams->zPosMin = LBDOUBLE_MAX;
		customParams->zPosMax = LBDOUBLE_MIN;
		customParams->xCosMin = LBDOUBLE_MAX;
		customParams->xCosMax = LBDOUBLE_MIN;
		customParams->yCosMin = LBDOUBLE_MAX;
		customParams->yCosMax = LBDOUBLE_MIN;
		customParams->zCosMin = LBDOUBLE_MAX;
		customParams->zCosMax = LBDOUBLE_MIN;
		customParams->objScattersMin = LBFOURBYTE_MAX;
		customParams->objScattersMax = LBFOURBYTE_MIN;
		customParams->colScattersMin = LBFOURBYTE_MAX;
		customParams->colScattersMax  = LBFOURBYTE_MIN;
		customParams->energyMin = LBDOUBLE_MAX;
		customParams->energyMax = LBDOUBLE_MIN;
		customParams->travelDistanceMin = LBDOUBLE_MAX;
		customParams->travelDistanceMax = LBDOUBLE_MIN;
		customParams->decayXPosMin = LBDOUBLE_MAX;
		customParams->decayXPosMax = LBDOUBLE_MIN;
		customParams->decayYPosMin = LBDOUBLE_MAX;
		customParams->decayYPosMax = LBDOUBLE_MIN;
		customParams->decayZPosMin = LBDOUBLE_MAX;
		customParams->decayZPosMax = LBDOUBLE_MIN;
		customParams->transaxialDistanceMin = LBDOUBLE_MAX;
		customParams->transaxialDistanceMax = LBDOUBLE_MIN;
		customParams->aaIndexMin = LBFOURBYTE_MAX;
		customParams->aaIndexMax = LBFOURBYTE_MIN;
		customParams->axialPositionMin = LBDOUBLE_MAX;
		customParams->axialPositionMax = LBDOUBLE_MIN;
		customParams->detectorXpositionMin = LBDOUBLE_MAX;
		customParams->detectorXpositionMax = LBDOUBLE_MIN;
		customParams->detectorYpositionMin = LBDOUBLE_MAX;
		customParams->detectorYpositionMax = LBDOUBLE_MIN;
		customParams->detectorZpositionMin = LBDOUBLE_MAX;
		customParams->detectorZpositionMax = LBDOUBLE_MIN;
		customParams->detectorAngleMin = LBDOUBLE_MAX;
		customParams->detectorAngleMax = LBDOUBLE_MIN;
		customParams->numDetInteractionsMin = LBFOURBYTE_MAX;
		customParams->numDetInteractionsMax = LBFOURBYTE_MIN;
		customParams->detectorAngleMin = LBDOUBLE_MAX;
		customParams->detectorAngleMax = LBDOUBLE_MIN;
			
		/* Loop through the parameters */
		while(LbPfGetParam(&paramFileHk, (void *)paramBuffer,
				&paramType, &paramSize, paramLabel, &isEOF)) {
				
				/* Find the runtime parameter */
				whichParam = PhoHFileLookupRunTimeParamLabel(paramLabel);
				
				
				switchOkay = true;
				switch (whichParam) {

					case PhoHFileEn_x_position:
						customParams->doXposition =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_x_pos_min:
						customParams->xPosMin =
							*((double *) paramBuffer);
						customParams->doXrange = true;
						break;

					case PhoHFileEn_x_pos_max:
						customParams->xPosMax =
							*((double *) paramBuffer);
						break;

					case PhoHFileEn_y_position:
						customParams->doYposition =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_y_pos_min:
						customParams->yPosMin =
							*((double *) paramBuffer);
						customParams->doYrange = true;
						break;

					case PhoHFileEn_y_pos_max:
						customParams->yPosMax =
							*((double *) paramBuffer);
						break;

					case PhoHFileEn_z_position:
						customParams->doZposition =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_z_pos_min:
						customParams->zPosMin =
							*((double *) paramBuffer);
						customParams->doZrange = true;
						break;

					case PhoHFileEn_z_pos_max:
						customParams->zPosMax =
							*((double *) paramBuffer);
						break;

					case PhoHFileEn_x_cosine:
						customParams->doXcosine =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_x_cos_min:
						customParams->xCosMin =
							*((double *) paramBuffer);
						customParams->doXCosRange = true;
						break;

					case PhoHFileEn_x_cos_max:
						customParams->xCosMax =
							*((double *) paramBuffer);
						break;

					case PhoHFileEn_y_cosine:
						customParams->doYcosine =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_y_cos_min:
						customParams->yCosMin =
							*((double *) paramBuffer);
						customParams->doYCosRange = true;
						break;

					case PhoHFileEn_y_cos_max:
						customParams->yCosMax =
							*((double *) paramBuffer);
						break;

					case PhoHFileEn_z_cosine:
						customParams->doZcosine =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_z_cos_min:
						customParams->zCosMin =
							*((double *) paramBuffer);
						customParams->doZCosRange = true;
						break;

					case PhoHFileEn_z_cos_max:
						customParams->zCosMax =
							*((double *) paramBuffer);
						break;

					case PhoHFileEn_scatters_in_object:
						customParams->doScattersInObject =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_obj_scatters_min:
						customParams->objScattersMin =
							*((LbUsFourByte *) paramBuffer);
						customParams->doObjScatRange = true;
						break;

					case PhoHFileEn_obj_scatters_max:
						customParams->objScattersMax =
							*((LbUsFourByte *) paramBuffer);
						break;

					case PhoHFileEn_scatters_in_collimator:
						customParams->doScattersInCollimator =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_col_scatters_min:
						customParams->colScattersMin =
							*((LbUsFourByte *) paramBuffer);
						customParams->doColScatRange = true;
						break;

					case PhoHFileEn_col_scatters_max:
						customParams->colScattersMax =
							*((LbUsFourByte *) paramBuffer);
						break;

					case PhoHFileEn_decay_weight:
						customParams->doDecayWeight =
							*((Boolean *) paramBuffer);

					case PhoHFileEn_weight:
						customParams->doWeight =
							*((Boolean *) paramBuffer);
						break;
						
					case PhoHFileEn_energy:
						customParams->doEnergy =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_energy_min:
						customParams->energyMin =
							*((double *) paramBuffer);
						customParams->doEnergyRange = true;
						break;

					case PhoHFileEn_energy_max:
						customParams->energyMax =
							*((double *) paramBuffer);
						break;
						
					case PhoHFileEn_travel_distance:
						customParams->doTravelDistance =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_travel_distance_min:
						customParams->travelDistanceMin =
							*((double *) paramBuffer);
						customParams->doTravDistRange = true;
						break;

					case PhoHFileEn_travel_distance_max:
						customParams->travelDistanceMax =
							*((double *) paramBuffer);
						break;
						
					case PhoHFileEn_decay_x_position:
						customParams->doDecayXposition =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_decay_x_pos_min:
						customParams->decayXPosMin =
							*((double *) paramBuffer);
						customParams->doDecayXrange = true;
						break;

					case PhoHFileEn_decay_x_pos_max:
						customParams->decayXPosMax =
							*((double *) paramBuffer);
						break;
						
					case PhoHFileEn_decay_y_position:
						customParams->doDecayYposition =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_decay_y_pos_min:
						customParams->decayYPosMin =
							*((double *) paramBuffer);
						customParams->doDecayYrange = true;
						break;

					case PhoHFileEn_decay_y_pos_max:
						customParams->decayYPosMax =
							*((double *) paramBuffer);
						break;
						
					case PhoHFileEn_decay_z_position:
						customParams->doDecayZposition =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_decay_z_pos_min:
						customParams->decayZPosMin =
							*((double *) paramBuffer);
						customParams->doDecayZrange = true;
						break;

					case PhoHFileEn_decay_z_pos_max:
						customParams->decayZPosMax =
							*((double *) paramBuffer);
						break;
						
					case PhoHFileEn_decay_time:
						customParams->doDecayTime =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_decay_type:
						customParams->doDecayType =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_transaxial_distance:
						customParams->doTransaxialDistance =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_transaxial_distance_min:
						customParams->transaxialDistanceMin =
							*((double *) paramBuffer);
						customParams->doTransDistRange = true;
						break;

					case PhoHFileEn_transaxial_distance_max:
						customParams->transaxialDistanceMax =
							*((double *) paramBuffer);
						break;
						
					case PhoHFileEn_azimuthal_angle_index:
						customParams->doAzimuthalAngleIndex =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_aa_index_min:
						customParams->aaIndexMin =
							*((LbUsFourByte *) paramBuffer);
						customParams->doAziAngRange = true;
						break;

					case PhoHFileEn_aa_index_max:
						customParams->aaIndexMax =
							*((LbUsFourByte *) paramBuffer);
						break;
						
					case PhoHFileEn_axial_position:
						customParams->doAxialPosition =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_axial_position_min:
						customParams->axialPositionMin =
							*((double *) paramBuffer);
						customParams->doAxiPosRange = true;
						break;

					case PhoHFileEn_axial_position_max:
						customParams->axialPositionMax =
							*((double *) paramBuffer);
						break;
						
					case 	PhoHFileEn_detector_x_position:
						customParams->doDetectorXposition =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_detector_x_position_min:
						customParams->detectorXpositionMin =
							*((double *) paramBuffer);
						customParams->doDetectorXposRange = true;
						break;

					case PhoHFileEn_detector_x_position_max:
						customParams->detectorXpositionMax =
							*((double *) paramBuffer);
						break;
						
					case 	PhoHFileEn_detector_y_position:
						customParams->doDetectorYposition =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_detector_y_position_min:
						customParams->detectorYpositionMin =
							*((double *) paramBuffer);
						customParams->doDetectorYposRange = true;
						break;

					case PhoHFileEn_detector_y_position_max:
						customParams->detectorYpositionMax =
							*((double *) paramBuffer);
						break;
						
					case 	PhoHFileEn_detector_z_position:
						customParams->doDetectorZposition =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_detector_z_position_min:
						customParams->detectorZpositionMin =
							*((double *) paramBuffer);
						customParams->doDetectorZposRange = true;
						break;

					case PhoHFileEn_detector_z_position_max:
						customParams->detectorZpositionMax =
							*((double *) paramBuffer);
						break;
						
					case 	PhoHFileEn_num_detector_interactions:
						customParams->doNumDetectorInteractions =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_num_detector_interactions_min:
						customParams->numDetInteractionsMin =
							*((LbFourByte *) paramBuffer);
						break;

					case PhoHFileEn_num_detector_interactions_max:
						customParams->numDetInteractionsMax =
							*((LbFourByte *) paramBuffer);
						customParams->doNumDetInteractionsRange = true;
						break;
						
					case 	PhoHFileEn_det_interaction_positions:
						customParams->doDetInteractionPos =
							*((Boolean *) paramBuffer);
						break;
						
					case 	PhoHFileEn_detector_angle:
						customParams->doDetectorAngle =
							*((Boolean *) paramBuffer);
						break;

					case PhoHFileEn_detector_angle_min:
						customParams->detectorAngleMin =
							*((double *) paramBuffer);
						customParams->doDetectorAngleRange = true;
						break;

					case PhoHFileEn_detector_angle_max:
						customParams->detectorAngleMax =
							*((double *) paramBuffer);
						break;
						
					case 	PhoHFileEn_detector_crystal:
						customParams->doDetectorCrystal =
							*((Boolean *) paramBuffer);
						break;

					default:
						sprintf(phoHFileErrString, "(PhoHFileGetRunTimeParams) Unknown (hence unused) parameter (%s).\n",
                            paramLabel);
						ErAlert(phoHFileErrString, false);
						break;
				}
				if (!switchOkay)
					break;
		}
		/* Close the parameter file */
		LbPfClose(&paramFileHk);

		
		/* See if we quit due to error */
		if (!isEOF)
			break;
		
		/* Verify that if we specified a maximum range value we also specified a minimum */
		if ((customParams->xPosMax != LBDOUBLE_MAX) && (customParams->xPosMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum x-position range in the history parameters, but not a minimum x-position\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
			break;
		}
		if ((customParams->yPosMax != LBDOUBLE_MAX) && (customParams->yPosMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum y-position range in the history parameters, but not a minimum y-position\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}	
		if ((customParams->zPosMax != LBDOUBLE_MAX) && (customParams->zPosMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum z-position range in the history parameters, but not a minimum z-position\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}	
		if ((customParams->xCosMax != LBDOUBLE_MAX) && (customParams->xCosMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum x-cosine range in the history parameters, but not a minimum x-cosine\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}	
		if ((customParams->yCosMax != LBDOUBLE_MAX) && (customParams->yCosMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum y-cosine range in the history parameters, but not a minimum y-cosine\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}	
		if ((customParams->zCosMax != LBDOUBLE_MAX) && (customParams->zCosMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum z-cosine range in the history parameters, but not a minimum z-cosine\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->objScattersMax != LBFOURBYTE_MAX) && (customParams->objScattersMin == LBFOURBYTE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum number of scatters in the object in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->colScattersMax != LBFOURBYTE_MAX) && (customParams->colScattersMin == LBFOURBYTE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum number of scatters in the collimator in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->energyMax != LBDOUBLE_MAX) && (customParams->energyMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum energy in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->travelDistanceMax != LBDOUBLE_MAX) && (customParams->travelDistanceMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum travel distance in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->decayXPosMax != LBDOUBLE_MAX) && (customParams->decayXPosMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum travel distance in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->decayYPosMax != LBDOUBLE_MAX) && (customParams->decayYPosMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum travel distance in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->decayZPosMax != LBDOUBLE_MAX) && (customParams->decayZPosMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum travel distance in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->transaxialDistanceMax != LBDOUBLE_MAX) && (customParams->transaxialDistanceMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum transaxial distance in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->aaIndexMax != LBFOURBYTE_MAX) && (customParams->aaIndexMin == LBFOURBYTE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum azimuthal angle index in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->axialPositionMax != LBDOUBLE_MAX) && (customParams->axialPositionMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum axial position in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->detectorXpositionMax != LBDOUBLE_MAX) && (customParams->detectorXpositionMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum x-position within the detector in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->detectorYpositionMax != LBDOUBLE_MAX) && (customParams->detectorYpositionMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum y-position within the detector in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->detectorZpositionMax != LBDOUBLE_MAX) && (customParams->detectorZpositionMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum z-position within the detector in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		if ((customParams->numDetInteractionsMax != LBFOURBYTE_MAX) && (customParams->numDetInteractionsMin == LBFOURBYTE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum number of detector interactions in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}		
		if ((customParams->detectorAngleMax != LBDOUBLE_MAX) && (customParams->detectorAngleMin == LBDOUBLE_MIN)) {
			sprintf(phoHFileErrString,"You specified a maximum detector angle in the history parameters, but not a minimum\n"
				"Your history parameter range file is named %s, please correct this to either no maximum specification\n"
				"or a specification for both maximum and minimum", paramPath);
			ErStGeneric(phoHFileErrString);
		}
		okay = true;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:			PhoHFileLookupRunTimeParamLabel
*
*			Summary:	Find the given label in the label table and return
*						its enumeration value.
*			Arguments:
*				char		*label	- The label to find.
*			Function return: Enumerated type corresponding to label.
*
*********************************************************************************/
PhoHFileEn_RunTimeParamsTy  PhoHFileLookupRunTimeParamLabel(char *label)	
{
	PhoHFileEn_RunTimeParamsTy	whichParam;			/* Parameter we find */
	LbUsFourByte			labelIndex;			/* Label we are checking */

	/* Set param to null */
	whichParam = PhoHFileEn_NULL;
	
	/* Search table for given label */
	for (labelIndex = 0; labelIndex < (PhoHFileEn_NULL); labelIndex++){
	
		/* See if it matches */
		if (strcmp(label, PhoHFileRunTimeParamLabels[labelIndex]) == 0) {
			whichParam = (PhoHFileEn_RunTimeParamsTy) labelIndex;
			break;
		}
	}
	
	return (whichParam);
}

/*********************************************************************************
*
*			Name:		PhoHFileReadEvent
*
*			Summary:	Read the next event from the history file.
*
*			Arguments:
*				FILE				*historyFile	- The history file.
*				PHG_Decay 			*decayPtr		- Storage for decay.
*				PHG_DetectedPhoton	*photonPtr		- Storage for photon.
*				(Note:  decayPtr can equal photonPtr, if desired.)
*
*			Function return: PhoHFileEventType corresponding to event type.
*
*********************************************************************************/
PhoHFileEventType PhoHFileReadEvent(FILE *historyFile, 
			PHG_Decay *decayPtr, PHG_DetectedPhoton *photonPtr)
{
	PhoHFileEventType	eventType = PhoHFileNullEvent;	/* The event we read */
	LbUsOneByte			flag;				/* Storage for type of event to read */
	LbUsFourByte		decaySize;			/* Size of a decay */
	LbUsFourByte		photonSize;			/* Size of a photon */
	
	do { /* Process Loop */

		/* See what type of event we have */
		if ((fread(&flag, sizeof(LbUsOneByte), 1, historyFile)) != 1) {
		
			/* See if we are not at end of file */
			if (feof(historyFile) == 0) {
				ErAbort("Unable to read event type.");
			}
			else {
				/* We are end of file so just break */
				break;
			}
		}
	
		/* See if we have a decay or a photon */
		if (PHG_IsADecay((LbUsFourByte) flag)) {
			
			/* Set our event type */
			eventType = PhoHFileDecayEvent;
			
			/* Read the decay */
			decaySize = sizeof(PHG_Decay);
			
			if ((fread(decayPtr, decaySize, 1, historyFile)) != 1) {
				ErAbort("Unable to read decay.");
			}
			
		}
		else if (PHG_IsAPhoton((LbUsFourByte) flag)) {
			/* Set our event type */
			eventType = PhoHFilePhotonEvent;
			
			/* Read the photon */
			photonSize = sizeof(PHG_DetectedPhoton);
			
			if ((fread(photonPtr, photonSize, 1, historyFile)) != 1) {
				ErAbort("Unable to read photon.");
			}
		}
		else
			ErAbort("Unexpected type from flag in PhoHFileReadEvent.");
		
	} while (false);
	
	return (eventType);
}

/*********************************************************************************
*
*			Name:		PhoHFileOldReadEvent
*
*			Summary:	Read the next event from an old history file.
*				Old list mode files did not have all the
*				info needed for processing SPECT/DHCI collimator
*				and detector list mode data.  However, this routine
*				can be used for reading other old list mode files
*				(e.g. PHG list mode files) which can be
*				processed.
*
*			Arguments:
*				FILE				*historyFile	- The history file.
*				PHG_Decay 			*decayPtr		- Storage for decay.
*				PHG_DetectedPhoton	*photonPtr		- Storage for photon.
*				(Note:  decayPtr can equal photonPtr, if desired.)
*
*			Function return: PhoHFileEventType corresponding to event type.
*
*********************************************************************************/
PhoHFileEventType PhoHFileOldReadEvent(
			FILE *historyFile, 
			PHG_OldDecay *oldDecayPtr,
			PHG_DetectedPhoton *photonPtr,
			Boolean isOldPhotons1,
			Boolean isOldPhotons2 )
{
	PhoHFileEventType	eventType = PhoHFileNullEvent;	/* The event we read */
	LbUsOneByte			flag;				/* Storage for type of event to read */
	LbUsFourByte		decaySize;			/* Size of a decay */
	LbUsFourByte		photonSize;			/* Size of a photon */
	PHG_OldDetectedPhoton1 oldPhoton1;	/* oldest photon type for reading 2.6.2.4 and earlier list mode */
	PHG_OldDetectedPhoton2 oldPhoton2;	/* old photon type for reading 2.6.2.5 & 2.6.2.6 list mode */
	
	do { /* Process Loop */

		/* See what type of event we have */
		if ((fread(&flag, sizeof(LbUsOneByte), 1, historyFile)) != 1) {
		
			/* See if we are not at end of file */
			if (feof(historyFile) == 0) {
				ErAbort("Unable to read event type.");
			}
			else {
				/* We are end of file so just break */
				break;
			}
		}
	
		/* See if we have a decay or a photon */
		if (PHG_IsADecay((LbUsFourByte) flag)) {
			
			/* Set our event type */
			eventType = PhoHFileDecayEvent;
			
			/* Read the decay */
			decaySize = sizeof(PHG_OldDecay);
			
			if ((fread(oldDecayPtr, decaySize, 1, historyFile)) != 1) {
				ErAbort("Unable to read decay.");
			}
			
		}
		else if (PHG_IsAPhoton((LbUsFourByte) flag)) {
			/* Set our event type */
			eventType = PhoHFilePhotonEvent;
			
			/* Read in a new photon or old photon record */
			if (isOldPhotons1) {
				
				/* Read the photon */
				photonSize = sizeof(PHG_OldDetectedPhoton1);
				
				if ((fread(&oldPhoton1, photonSize, 1, historyFile)) != 1) {
					ErAbort("Unable to read photon.");
				}
				
				/* Assign the new photon record from the old photon record */
				photonPtr->location = oldPhoton1.location;
				photonPtr->angle = oldPhoton1.angle;
				photonPtr->flags = oldPhoton1.flags;
				photonPtr->photon_weight = oldPhoton1.photon_weight;
				photonPtr->energy = oldPhoton1.energy;
				photonPtr->time_since_creation = oldPhoton1.time_since_creation;
				photonPtr->transaxialPosition = oldPhoton1.transaxialPosition;
				photonPtr->azimuthalAngleIndex = oldPhoton1.azimuthalAngleIndex;
				photonPtr->detectorAngle = -1e9; /* a nonsense value - it should never be used */
				photonPtr->detCrystal = -1; /* a nonsense value - it should never be used */
				
			} else if (isOldPhotons2) {
				
				/* Read the photon */
				photonSize = sizeof(PHG_OldDetectedPhoton2);
				
				if ((fread(&oldPhoton2, photonSize, 1, historyFile)) != 1) {
					ErAbort("Unable to read photon.");
				}
				
				/* Assign the new photon record from the old photon record */
				photonPtr->location = oldPhoton2.location;
				photonPtr->angle = oldPhoton2.angle;
				photonPtr->flags = oldPhoton2.flags;
				photonPtr->photon_weight = oldPhoton2.photon_weight;
				photonPtr->energy = oldPhoton2.energy;
				photonPtr->time_since_creation = oldPhoton2.time_since_creation;
				photonPtr->transaxialPosition = oldPhoton2.transaxialPosition;
				photonPtr->azimuthalAngleIndex = oldPhoton2.azimuthalAngleIndex;
				photonPtr->detectorAngle = oldPhoton2.detectorAngle;
				photonPtr->detCrystal = -1; /* a nonsense value - it should never be used */
				
			} else {
			
				/* Read the photon */
				photonSize = sizeof(PHG_DetectedPhoton);
				
				if ((fread(photonPtr, photonSize, 1, historyFile)) != 1) {
					ErAbort("Unable to read photon.");
				}
			
			}
		}
		else
			ErAbort("Unexpected type from flag in PhoHFileOldReadEvent.");
		
	} while (false);
	
	return (eventType);
}

#undef	PHOTON_HIST_FILE
