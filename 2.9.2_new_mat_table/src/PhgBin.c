/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1993-2013 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		PhgBin.c
*			Revision Number:	2.6
*			Date last revised:	23 July 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	28 July 1993
*
*			Module Overview:	Definitions for PhgBin.c.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:
*				PhgBinInitialize
*				PhgBinPhotons
*				PhgBinTerminate
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
*			Revision date:		1 February 2012
*
*			Revision description:	Changed form of user functions to pointers
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
*						- binning by random-state and crystal number.
*						- binning PET data as SPECT
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		February 2004
*
*			Revision description:	Changed PhgBinCompSpectDA to fix bug that
*					caused SimSET to exit with an error when simulating
*					simple SPECT or cylindrical detectors.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:	Added support for time-of-flight binning.
*
*********************************************************************************/

#define PHG_BIN


#include <stdio.h>
#include <memory.h>
#include <limits.h>

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
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "PhgMath.h"
#include "CylPos.h"
#include "PhoHFile.h"
#include "Collimator.h"
#include "Detector.h"
#include "PhgHdr.h"
#include "phg.h"
#include "PhgUsrBin.h"
#include "PhgBin.h"


/* Local Prototypes */
			
void	phgBinIncrementPETImage(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData,
			LbUsEightByte imageIndex, PHG_Decay *decay,
			PHG_TrackingPhoton *bluePhoton, PHG_TrackingPhoton *pinkPhoton,
			double coincidenceWeight, double coincidenceSquWeight);

/* Global variables */
static	char	phgBinErrStr[1024];			/* Storage for creating error strings */
static	double	phgBinDetDiameter;			/* Detector diameter for MSRB */
static	double	phgBinObjDiameter;			/* Image/Object diameter for MSRB */


/*********************************************************************************
*
*			Name:			PhgBinInitParams
*
*			Summary:		Initializes the binning parameters, without allocating
*							any memory, opening files, etc.
*
*			NOTE: This routine is called by PhgBinInitalize. It is available to the
*				  outside world to allow the use of the binning parameters without
*				  doing the binning. Consequently, it should never be called
*				  separately if PhgBinInitialize has been or will be called.
*
*			Arguments:
*				char			*paramsName	 - The name of the parameter file.
*				PHG_BinParamsTy	*binParams	 - User defined binning parameters.
*				PHG_BiDataTy	*binData	 - Storage for binned data.
*				PHG_BinFieldsTy *binFields	 - Various binning information.
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean PhgBinInitParams(char *paramsName, PHG_BinParamsTy *binParams, PHG_BinDataTy *binData, PHG_BinFieldsTy *binFields)
{
	LbUsFourByte	dimensionIndex;		/* Index for processing dimension calculations */
	LbUsFourByte	prevSize;			/* Size of current dimension for intialization */
	LbUsFourByte	prevBins;			/* Number of bins in previous dimension for initialization */
		
	do { 	/* Process Loop */
			
		/* Set certain defaults - note that some are set in PhgParams as well */
		binParams->detector_radius = -1.0;

		/* Initialize size variables in case user doesn't specify all params in parameter file */
		{
			binParams->scatter2CIsize = 0;
			binParams->scatter2WIsize = 0;
			binParams->scatter2WISsize = 0;
			
			binParams->scatter1CIsize = 0;
			binParams->scatter1WIsize = 0;
			binParams->scatter1WISsize = 0;
			
			binParams->energy2CIsize = 0;
			binParams->energy2WIsize = 0;
			binParams->energy2WISsize = 0;
			
			binParams->energy1CIsize = 0;
			binParams->energy1WIsize = 0;
			binParams->energy1WISsize = 0;
			
			binParams->crystal2CIsize = 0;
			binParams->crystal2WIsize = 0;
			binParams->crystal2WISsize = 0;
			
			binParams->crystal1CIsize = 0;
			binParams->crystal1WIsize = 0;
			binParams->crystal1WISsize = 0;
			
			binParams->aaCIsize = 0;
			binParams->aaWIsize = 0;
			binParams->aaWISsize = 0;
			
			binParams->tdCIsize = 0;
			binParams->tdWIsize = 0;
			binParams->tdWISsize = 0;
			
			binParams->tofCIsize = 0;
			binParams->tofWIsize = 0;
			binParams->tofWISsize = 0;
			
			binParams->paCIsize = 0;
			binParams->paWIsize = 0;
			binParams->paWISsize = 0;
			
			binParams->z2CIsize = 0;
			binParams->z2WIsize = 0;
			binParams->z2WISsize = 0;
			
			binParams->z1CIsize = 0;
			binParams->z1WIsize = 0;
			binParams->z1WISsize = 0;
			
			binParams->phiCIsize = 0;
			binParams->phiWIsize = 0;
			binParams->phiWISsize = 0;
			
			binParams->thetaCIsize = 0;
			binParams->thetaWIsize = 0;
			binParams->thetaWISsize = 0;
			
			binParams->xrCIsize = 0;
			binParams->xrWIsize = 0;
			binParams->xrWISsize = 0;
			
			binParams->yrCIsize = 0;
			binParams->yrWIsize = 0;
			binParams->yrWISsize = 0;
		}	
		/* Clear out dimensions */
		for (dimensionIndex = 0; dimensionIndex < PHGBIN_NUM_DIMENSIONS; dimensionIndex++)
			binParams->PhgBinDimensions[dimensionIndex] = PhgBinEn_Null;
			
		/* Get the binning parameters */
		if (PhgGetBinParams(paramsName, binParams) == false) {
			break;
		}
		
		/* Convert elevation angle ranges to radians from degrees */
		binParams->minTheta = PHGMATH_RadiansFromDegrees(binParams->minTheta);
		binParams->maxTheta = PHGMATH_RadiansFromDegrees(binParams->maxTheta);
		
		/* if numE1Bins was set, but numE2Bins wasn't, set numE2Bins = numE1Bins
		if ( (binParams->numE1Bins >= 1) && (binParams->numE2Bins == 0) ) {
			binParams->numE2Bins = binParams->numE1Bins;
		} */
		
		/* Set min and max 'phi' for 3D Binning */
		if(binParams->numPHIBins > 0)
		{
			binParams->minPhi = -(PHGMATH_PI/(2*binParams->numPHIBins));
			binParams->maxPhi = PHGMATH_PI - (PHGMATH_PI/(2*binParams->numPHIBins));
		}
		else {
			binParams->minPhi = 0;
			binParams->maxPhi = 0;
		}
		
		/* If user is doing MSRB binning and they didn't specify a radius, set default */
		if (( binParams->doMSRB == true) && (binParams->detector_radius == -1.0) && PHG_IsDetectOnTheFly()) {
			binParams->detector_radius = DetGtInsideRadius();
		}
		else if (( binParams->doMSRB == true) && (binParams->detector_radius == -1.0)) {
			binParams->detector_radius = CylPosGetTargetRadius();
		}
		else if (( binParams->doMSRB == true) && (binParams->detector_radius == 0.0)) {
			sprintf(phgBinErrStr, "You have specified MSRB binning but your detector radius is %3.2f, check your binning parameters file (PhgBinInitialize)",
				binParams->detector_radius);
			ErStGeneric(phgBinErrStr);
			break;
		}
		
		/* If user is doing MSRB binning Set image radius */
		if ( binParams->doMSRB == true) {
			
			binParams->image_radius = PHGMATH_Max(fabs(binParams->minTD), fabs(binParams->maxTD));
			
			if (binParams->image_radius <= 0.0)
				binParams->image_radius = CylPosGetObjectRadius();
		}

		/* Setup optimization variables for MSRB if doing */
		if (binParams->doMSRB == true) {
			phgBinObjDiameter = binParams->image_radius * 2;
			phgBinDetDiameter = binParams->detector_radius * 2;
			
			/* Just to be REALLY careful */
			#ifdef PHG_DEBUG
				if ((phgBinObjDiameter <= 0.0) || (phgBinDetDiameter <= 0.0)) {
						sprintf(phgBinErrStr, "Computed invalid object diameter and detector diameter of %3.2f and %3.2f respectively (PhgBinInitialize)",
					phgBinObjDiameter, phgBinDetDiameter);
				ErStGeneric(phgBinErrStr);
				break;
			}
			#endif
		}
		
		/* Perform some validation on the input parameters */
		{
			/* Check min/max Z */
			if ((binParams->numZBins > 0) && (binParams->minZ > binParams->maxZ)){
				sprintf(phgBinErrStr, "Invalid values for min/max Z in binning parameters, (min = %3.2f) (max = %3.2f) (PhgBinInitialize)",
					binParams->minZ, binParams->maxZ);
				ErStGeneric(phgBinErrStr);
				break;
			}

			/* Check min/max TD */
			if ((binParams->numTDBins > 0) && (binParams->minTD > binParams->maxTD)){
				sprintf(phgBinErrStr, "Invalid values for min/max TD in binning parameters, (min = %3.2f) (max = %3.2f) (PhgBinInitialize)",
					binParams->minTD, binParams->maxTD);
				ErStGeneric(phgBinErrStr);
				break;
			}

			/* Check min/max Energy*/
			if ((binParams->numE1Bins >= 1) && (binParams->minE > binParams->maxE)){
				sprintf(phgBinErrStr, "Invalid values for min/max Energy in binning parameters, (min = %3.2f) (max = %3.2f) (PhgBinInitialize)",
					binParams->minE, binParams->maxE);
				ErStGeneric(phgBinErrStr);
				break;
			}

			/* Check min/max scatter */
			if (binParams->minS > binParams->maxS){
				sprintf(phgBinErrStr, "Invalid values for min/max Scatters in binning parameters, (min = %ld) (max = %ld) (PhgBinInitialize)",
					(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
				ErStGeneric(phgBinErrStr);
				break;
			}
			
			/* If we are doing SPECT and collimating on the fly then we need
				to verify that NumViews == NumAngleBins, or that NumAngleBins
				== 1
				   8/3/01 - to support singles check is changed to only occur when doing SPECT/UNC: SDV
			*/
			if (binParams->numAABins != 1) {
				if (PHG_IsCollimateOnTheFly() && PHG_IsSPECT() && (ColIsUNC() == true))
					if ((ColGtNumViews() != binParams->numAABins) && (binParams->numAABins > 0)){
						sprintf(phgBinErrStr, "Num views in collimator parameters must match"
							" the number of angle bins in the binning parameters file,"
							" (num views = %ld) (num angle bins = %ld) (PhgBinInitialize)",
							(unsigned long)ColGtNumViews(), (unsigned long)binParams->numAABins);
						ErStGeneric(phgBinErrStr);
						break;					
				}
			}
			
			/* If we are doing planar detection then just verify that the user didn't specify
				conflicting parameters
			*/
			if (binParams->numAABins != 1) {
				if (PHG_IsDetectOnTheFly() && PHG_IsSPECT())
					if ((DetGtNumViews() !=  0) 
						&& (DetGtNumViews() != (LbFourByte)(binParams->numAABins))
						&& (binParams->numAABins != 0)) {
						sprintf(phgBinErrStr, "Num views in detector parameters must match"
							" the number of angle bins in the binning parameters file,"
							" (num views = %ld) (num angle bins = %ld) (PhgBinInitialize)",
							(unsigned long)DetGtNumViews(), (unsigned long)binParams->numAABins);
						ErStGeneric(phgBinErrStr);
						break;					
				}
			}
			
			/* binParams->scatterRandomParam from 4-10 are PET-only binning options */
			if ( (binParams->scatterRandomParam > 3) && (PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
				sprintf(phgBinErrStr, "This scatter-randoms binning option not compatible with"
					" SPECT binning,\n"
					" (scatter_random_param or scatter_param = %ld) (PhgBinInitialize)",
					(unsigned long)binParams->scatterRandomParam);
				ErStGeneric(phgBinErrStr);
				break;					
			}
			
			/* binParams->scatterRandomParam from 6-10 don't make sense with binParams->acceptRandoms == false */
			if ( (binParams->scatterRandomParam > 5) && (!binParams->acceptRandoms) ) {
				sprintf(phgBinErrStr, "This scatter-randoms binning option requires random"
					" coincidence binning,\n"
					" (scatter_random_param or scatter_param = %ld and accept_randoms = false) (PhgBinInitialize)",
					(unsigned long)binParams->scatterRandomParam);
				ErStGeneric(phgBinErrStr);
				break;					
			}
			
			/* Compute scatter range here, numS[1|2]Bins may depend on it */
			binParams->sRange = binParams->maxS - binParams->minS + 1;
				
			/* Compute number of scatter/random bins from scatter/random param */
			if (binParams->scatterRandomParam == 0) {
				binParams->numS1Bins = 1;
				binParams->numS2Bins = 1;
			}
			else if (binParams->scatterRandomParam == 1) {
				binParams->numS1Bins = 2;
				binParams->numS2Bins = 1;
			}
			else if ( (binParams->scatterRandomParam == 2) || (binParams->scatterRandomParam == 3) ) {
				if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
					binParams->numS1Bins = binParams->sRange;
					binParams->numS2Bins = binParams->sRange;
				}
				else {
					binParams->numS1Bins = binParams->sRange;
					binParams->numS2Bins = 1;
				}
			} 
			/* the rest of the options are PET-only */
			else if ( (binParams->scatterRandomParam == 4) || (binParams->scatterRandomParam == 5) ) {
				binParams->numS1Bins = binParams->sRange;
				binParams->numS2Bins = 1;
			}
			/* the rest of these options mirror 1-5, but have an extra bin for randoms */
			else if (binParams->scatterRandomParam == 6) {
				binParams->numS1Bins = 3;
				binParams->numS2Bins = 1;
			}
			/* to add the extra bin for randoms this has to be made 1D */
			else if ( (binParams->scatterRandomParam == 7) || (binParams->scatterRandomParam == 8) ) {
				binParams->numS1Bins = binParams->sRange * binParams->sRange + 1;
				binParams->numS2Bins = 1;
			} 
			else /* binParams->scatterRandomParam == 9,10 */  {
				binParams->numS1Bins = binParams->sRange + 1;
				binParams->numS2Bins = 1;
			}
			
			
			/* Time-of-flight binning cannnot be done with SPECT simulation */
			if ( (binParams->numTOFBins != 0) && (PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
				sprintf(phgBinErrStr, "Time-of-flight binning not compatible with"
					" SPECT binning,"
					" (num tof bins = %ld) (PhgBinInitialize)",
					(unsigned long)binParams->numTOFBins);
				ErStGeneric(phgBinErrStr);
				break;					
			}
			
			/* Check min/max TOF */
			if ((binParams->numTOFBins > 0) && (binParams->minTOF > binParams->maxTOF)){
				sprintf(phgBinErrStr, "Invalid values for min/max TOF in binning parameters, (min = %3.2f) (max = %3.2f) (PhgBinInitialize)",
					binParams->minTOF, binParams->maxTOF);
				ErStGeneric(phgBinErrStr);
				break;
			}

			
			/* Verify elevation angle within range */
			if (binParams->maxTheta > (PHGMATH_PI_DIV2)) {
				sprintf(phgBinErrStr, "The elevlation angle maximum specified in the binning parameters can not be > 90.0 . (PhgBinInitialize)");
				ErStGeneric(phgBinErrStr);
				break;
			}
			
			/* Verify elevation angle within range */
			if (binParams->minTheta < -(PHGMATH_PI_DIV2)) {
				sprintf(phgBinErrStr, "The elevlation angle minimum specified in the binning parameters can not be < -90.0 . (PhgBinInitialize)");
				ErStGeneric(phgBinErrStr);
				break;
			}
			
			/* Verify we don't have conflicting switches for SSRB & MSRB */
			if ((binParams->doSSRB) && (binParams->doMSRB)) {
				ErStGeneric("You have specified both SSRB and MSRB in the binning parameters, this won't work. (PhgBinInitialize)");
				break;
			}
			
			/* Verify we don't have conflicting switches for SSRB & MSRB */
			if ((binParams->doSSRB || binParams->doMSRB) && (binParams->numZBins <= 1)) {
				ErStGeneric("You have specified either SSRB or MSRB without any Z bins, this won't do anything. (PhgBinInitialize)");
				break;
			}
			
		}
		
		/* Initialize memory buffers to zero for error handling */
		binData->countImage = 0;
		binData->weightImage = 0;
		binData->weightSquImage = 0;
				
		/* Tie polar angle bins to one for now - this is currently unused */
		binParams->numPABins = 1;
		
		/* Set azimuthal angle range (depends on simulation type) */
		if (PHG_IsSPECT() || binParams->isBinPETasSPECT) {
			binParams->minAA = 0;
			binParams->maxAA = PHGMATH_2PI;
		}
		else {
			binParams->minAA = 0;
			binParams->maxAA = PHGMATH_PI;
		}
		
		/* Compute remaining ranges */
		binParams->zRange = binParams->maxZ - binParams->minZ;
		binParams->eRange = binParams->maxE - binParams->minE;
		binParams->tdRange = binParams->maxTD - binParams->minTD;
		binParams->tofRange = binParams->maxTOF - binParams->minTOF;
		binParams->phiRange = binParams->maxPhi - binParams->minPhi;
		binParams->thetaRange = binParams->maxTheta - binParams->minTheta;
		binParams->xrRange = binParams->maxXR - binParams->minXR;
		binParams->yrRange = binParams->maxYR - binParams->minYR;
		
		/* Compute dimension sizes based on ordering of parameter file */
		prevSize = 0;
		prevBins = 0;
		
		for (dimensionIndex = 0; dimensionIndex < PHGBIN_NUM_DIMENSIONS; dimensionIndex++){
		
			switch(binParams->PhgBinDimensions[dimensionIndex]) {
			
				case PhgBinEn_TD:
					if (prevSize != 0) {
						binParams->tdCIsize = prevSize * prevBins;
						binParams->tdWIsize = prevSize * prevBins;
						binParams->tdWISsize = prevSize * prevBins;
					}
					else {
						binParams->tdCIsize = 1;
						binParams->tdWIsize = 1;
						binParams->tdWISsize = 1;
					}
					
					prevSize = binParams->tdCIsize;
					prevBins = (binParams->numTDBins > 0) ? binParams->numTDBins : 1;
					break;
			
				case PhgBinEn_THETA:
					if (prevSize != 0) {
						binParams->thetaCIsize = prevSize * prevBins;
						binParams->thetaWIsize = prevSize * prevBins;
						binParams->thetaWISsize = prevSize * prevBins;
					}
					else {
						binParams->thetaCIsize = 1;
						binParams->thetaWIsize = 1;
						binParams->thetaWISsize = 1;
					}
					
					prevSize = binParams->thetaCIsize;
					prevBins = (binParams->numThetaBins > 0) ? binParams->numThetaBins : 1;
					break;
			
			
				case PhgBinEn_PHI:
					if (prevSize != 0) {
						binParams->phiCIsize = prevSize * prevBins;
						binParams->phiWIsize = prevSize * prevBins;
						binParams->phiWISsize = prevSize * prevBins;
					}
					else {
						binParams->phiCIsize = 1;
						binParams->phiWIsize = 1;
						binParams->phiWISsize = 1;
					}
					
					prevSize = binParams->phiCIsize;
					prevBins = (binParams->numPHIBins > 0) ? binParams->numPHIBins : 1;
					break;
			
			
				case PhgBinEn_XR:
					if (prevSize != 0) {
						binParams->xrCIsize = prevSize * prevBins;
						binParams->xrWIsize = prevSize * prevBins;
						binParams->xrWISsize = prevSize * prevBins;
					}
					else {
						binParams->xrCIsize = 1;
						binParams->xrWIsize = 1;
						binParams->xrWISsize = 1;
					}
					
					prevSize = binParams->xrCIsize;
					prevBins = (binParams->numXRBins > 0) ? binParams->numXRBins : 1;
					break;
			
			
				case PhgBinEn_YR:
					if (prevSize != 0) {
						binParams->yrCIsize = prevSize * prevBins;
						binParams->yrWIsize = prevSize * prevBins;
						binParams->yrWISsize = prevSize * prevBins;
					}
					else {
						binParams->yrCIsize = 1;
						binParams->yrWIsize = 1;
						binParams->yrWISsize = 1;
					}
					
					prevSize = binParams->yrCIsize;
					prevBins = (binParams->numYRBins > 0) ? binParams->numYRBins : 1;
					break;
			
				case PhgBinEn_AA:
					if (prevSize != 0) {
						binParams->aaCIsize = prevSize * prevBins;
						binParams->aaWIsize = prevSize * prevBins;
						binParams->aaWISsize = prevSize * prevBins;
					}
					else {
						binParams->aaCIsize = 1;
						binParams->aaWIsize = 1;
						binParams->aaWISsize = 1;
					}
					
					prevSize = binParams->aaCIsize;
					prevBins = (binParams->numAABins > 0) ? binParams->numAABins : 1;
					break;
		
		
				case PhgBinEn_Crystal1:
					if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
						if (prevSize != 0) {
							binParams->crystal2CIsize = prevSize * prevBins;
							binParams->crystal2WIsize = prevSize * prevBins;
							binParams->crystal2WISsize = prevSize * prevBins;
							
							/* The size of the next fastest varying parameter is based on the previous (faster) parameter */
							binParams->crystal1CIsize = ((binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1) * binParams->crystal2CIsize;
							binParams->crystal1WIsize = ((binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1) * binParams->crystal2WIsize;
							binParams->crystal1WISsize = ((binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1) * binParams->crystal2WISsize;
						}
						else {
							binParams->crystal2CIsize = 1;
							binParams->crystal2WIsize = 1;
							binParams->crystal2WISsize = 1;
							
							/* The size of the next fastest varying parameter is based on the previous (faster) parameter */
							binParams->crystal1CIsize = ((binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1) * binParams->crystal2CIsize;
							binParams->crystal1WIsize = ((binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1) * binParams->crystal2WIsize;
							binParams->crystal1WISsize = ((binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1) * binParams->crystal2WISsize;
						}
					}
					else {
						if (prevSize != 0) {
							binParams->crystal2CIsize = 0;
							binParams->crystal2WIsize = 0;
							binParams->crystal2WISsize = 0;
							
							/* The size of the next fastest varying parameter is based on the previous (faster) parameter */
							binParams->crystal1CIsize = prevSize * prevBins;
							binParams->crystal1WIsize = prevSize * prevBins;
							binParams->crystal1WISsize = prevSize * prevBins;
						}
						else {
						
							/* For SPECT we clear E2 just to be complete */
							binParams->crystal2CIsize = 0;
							binParams->crystal2WIsize = 0;
							binParams->crystal2WISsize = 0;
							
							/* The size of the next fastest varying parameter is based on the previous (faster) parameter */
							binParams->crystal1CIsize = 1;
							binParams->crystal1WIsize = 1;
							binParams->crystal1WISsize = 1;
						}
						
					}
					prevSize = binParams->crystal1CIsize;
					prevBins = (binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1;
					break;
			
				case PhgBinEn_Crystal2:
				
					/* We took care of crystal binning on crystal1 above */
					break;
		
						
				case PhgBinEn_TOF:
					if (prevSize != 0) {
						binParams->tofCIsize = prevSize * prevBins;
						binParams->tofWIsize = prevSize * prevBins;
						binParams->tofWISsize = prevSize * prevBins;
					}
					else {
						binParams->tofCIsize = 1;
						binParams->tofWIsize = 1;
						binParams->tofWISsize = 1;
					}
					
					prevSize = binParams->tofCIsize;
					prevBins = (binParams->numTOFBins > 0) ? binParams->numTOFBins : 1;
					break;
		
		
				case PhgBinEn_Energy1:
					if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
						if (prevSize != 0) {
							binParams->energy2CIsize = prevSize * prevBins;
							binParams->energy2WIsize = prevSize * prevBins;
							binParams->energy2WISsize = prevSize * prevBins;
							
							/* The size of the next fastest varying parameter is based on the previous (faster) parameter */
							binParams->energy1CIsize = ((binParams->numE2Bins > 0) ? binParams->numE2Bins : 1) * binParams->energy2CIsize;
							binParams->energy1WIsize = ((binParams->numE2Bins > 0) ? binParams->numE2Bins : 1) * binParams->energy2WIsize;
							binParams->energy1WISsize = ((binParams->numE2Bins > 0) ? binParams->numE2Bins : 1) * binParams->energy2WISsize;
						}
						else {
							binParams->energy2CIsize = 1;
							binParams->energy2WIsize = 1;
							binParams->energy2WISsize = 1;
							
							/* The size of the next fastest varying parameter is based on the previous (faster) parameter */
							binParams->energy1CIsize = ((binParams->numE2Bins > 0) ? binParams->numE2Bins : 1) * binParams->energy2CIsize;
							binParams->energy1WIsize = ((binParams->numE2Bins > 0) ? binParams->numE2Bins : 1) * binParams->energy2WIsize;
							binParams->energy1WISsize = ((binParams->numE2Bins > 0) ? binParams->numE2Bins : 1) * binParams->energy2WISsize;
						}
					}
					else {
						if (prevSize != 0) {
							binParams->energy2CIsize = 0;
							binParams->energy2WIsize = 0;
							binParams->energy2WISsize = 0;
							
							/* The size of the next fastest varying parameter is based on the previous (faster) parameter */
							binParams->energy1CIsize = prevSize * prevBins;
							binParams->energy1WIsize = prevSize * prevBins;
							binParams->energy1WISsize = prevSize * prevBins;
						}
						else {
						
							/* For SPECT we clear E2 just to be complete */
							binParams->energy2CIsize = 0;
							binParams->energy2WIsize = 0;
							binParams->energy2WISsize = 0;
							
							/* The size of the next fastest varying parameter is based on the previous (faster) parameter */
							binParams->energy1CIsize = 1;
							binParams->energy1WIsize = 1;
							binParams->energy1WISsize = 1;
						}
						
					}
					prevSize = binParams->energy1CIsize;
					prevBins = (binParams->numE1Bins > 0) ? binParams->numE1Bins : 1;
					break;
			
				case PhgBinEn_Energy2:
				
					/* We took care of energy on energy1 */
					break;
		
						
				case PhgBinEn_Scatter1:
					if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
						if (prevSize != 0) {
							binParams->scatter2CIsize = prevSize * prevBins;
							binParams->scatter2WIsize = prevSize * prevBins;
							binParams->scatter2WISsize = prevSize * prevBins;
							
							/* The size of the next fastest varying parameter is based on the previous (faster) parameter */
							binParams->scatter1CIsize = ((binParams->numS2Bins > 0) ? binParams->numS2Bins : 1) * binParams->scatter2CIsize;
							binParams->scatter1WIsize = ((binParams->numS2Bins > 0) ? binParams->numS2Bins : 1) * binParams->scatter2WIsize;
							binParams->scatter1WISsize = ((binParams->numS2Bins > 0) ? binParams->numS2Bins : 1) * binParams->scatter2WISsize;
						}
						else {
							binParams->scatter2CIsize = 1;
							binParams->scatter2WIsize = 1;
							binParams->scatter2WISsize = 1;
						
							/* The size of the next fastest varying parameter is based on the previous (faster) parameter */
							binParams->scatter1CIsize = ((binParams->numS2Bins > 0) ? binParams->numS2Bins : 1) * binParams->scatter2CIsize;
							binParams->scatter1WIsize = ((binParams->numS2Bins > 0) ? binParams->numS2Bins : 1) * binParams->scatter2WIsize;
							binParams->scatter1WISsize = ((binParams->numS2Bins > 0) ? binParams->numS2Bins : 1) * binParams->scatter2WISsize;
						}
					}
					else {
						if (prevSize != 0) {
							binParams->scatter2CIsize = 0;
							binParams->scatter2WIsize = 0;
							binParams->scatter2WISsize = 0;
							
							/* The size of the next fastest varying parameter is based on the previous (faster) parameter */
							binParams->scatter1CIsize = prevSize * prevBins;
							binParams->scatter1WIsize = prevSize * prevBins;
							binParams->scatter1WISsize = prevSize * prevBins;
						}
						else {
							binParams->scatter2CIsize = 0;
							binParams->scatter2WIsize = 0;
							binParams->scatter2WISsize = 0;
						
							/* The size of the next fastest varying parameter is based on the previous (faster) parameter */
							binParams->scatter1CIsize = 1;
							binParams->scatter1WIsize = 1;
							binParams->scatter1WISsize = 1;
						}
					}
										
					prevSize = binParams->scatter1CIsize;
					prevBins = ((binParams->numS1Bins > 0) ? binParams->numS1Bins : 1);
					break;
		
				case PhgBinEn_Scatter2:
					/* We took care of scatters on scatter 1 */
					break;
			
				case PhgBinEn_Z1:
					if (prevSize != 0) {
						binParams->z1CIsize = prevSize * prevBins;
						binParams->z1WIsize = prevSize * prevBins;
						binParams->z1WISsize = prevSize * prevBins;
				
						if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) && (binParams->doSSRB == false) && (binParams->doMSRB == false)) {
							binParams->z2CIsize = ((binParams->numZBins > 0) ? binParams->numZBins : 1) * binParams->z1CIsize;
							binParams->z2WIsize = ((binParams->numZBins > 0) ? binParams->numZBins : 1) * binParams->z1WIsize;
							binParams->z2WISsize = ((binParams->numZBins > 0) ? binParams->numZBins : 1) * binParams->z1WISsize;
						}
						else {
							binParams->z2CIsize = binParams->z1CIsize;
							binParams->z2WIsize = binParams->z1WIsize;
							binParams->z2WISsize = binParams->z1WISsize;
						}
					}
					else {
						binParams->z1CIsize = 1;
						binParams->z1WIsize = 1;
						binParams->z1WISsize = 1;
				
						if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) && (binParams->doSSRB == false) && (binParams->doMSRB == false)) {
							binParams->z2CIsize = ((binParams->numZBins > 0) ? binParams->numZBins : 1) * binParams->z1CIsize;
							binParams->z2WIsize = ((binParams->numZBins > 0) ? binParams->numZBins : 1) * binParams->z1WIsize;
							binParams->z2WISsize = ((binParams->numZBins > 0) ? binParams->numZBins : 1) * binParams->z1WISsize;
						}
						else {
							binParams->z2CIsize = binParams->z1CIsize;
							binParams->z2WIsize = binParams->z1WIsize;
							binParams->z2WISsize = binParams->z1WISsize;
						}
					}
					prevSize = binParams->z2CIsize;
					prevBins = ((binParams->numZBins > 0) ? binParams->numZBins : 1);
					break;
					
				case PhgBinEn_Z2:
				
					/* We took care of this on Z1 */
					break;
					
				case PhgBinEn_Null:
					break;
					
				#ifdef PHG_DEBUG
					default:
						ErStGeneric("Invalid binning dimension in array (PhgBinInitialize)");
						goto FAIL;
				#endif
			}
		}
			
		/* We can now compute the sizes of the images, based on last dimension  */
		switch(binParams->PhgBinDimensions[PHGBIN_NUM_DIMENSIONS-1]) {
		
			case PhgBinEn_TD:
				binParams->countImageSize = ((binParams->numTDBins > 0) ? binParams->numTDBins : 1) * binParams->tdCIsize;
				binParams->weightImageSize = ((binParams->numTDBins > 0) ? binParams->numTDBins : 1) * binParams->tdWIsize;
				binParams->weightSquImageSize = ((binParams->numTDBins > 0) ? binParams->numTDBins : 1) * binParams->tdWISsize;
				break;
		
			case PhgBinEn_AA:
				binParams->countImageSize = ((binParams->numAABins > 0) ? binParams->numAABins : 1) * binParams->aaCIsize;
				binParams->weightImageSize = ((binParams->numAABins > 0) ? binParams->numAABins : 1) * binParams->aaWIsize;
				binParams->weightSquImageSize = ((binParams->numAABins > 0) ? binParams->numAABins : 1) * binParams->aaWISsize;
				break;
	
			case PhgBinEn_TOF:
				binParams->countImageSize = ((binParams->numTOFBins > 0) ? binParams->numTOFBins : 1) * binParams->tofCIsize;
				binParams->weightImageSize = ((binParams->numTOFBins > 0) ? binParams->numTOFBins : 1) * binParams->tofWIsize;
				binParams->weightSquImageSize = ((binParams->numTOFBins > 0) ? binParams->numTOFBins : 1) * binParams->tofWISsize;
				break;
	
			case PhgBinEn_Energy1:
				binParams->countImageSize = ((binParams->numE1Bins > 0) ? binParams->numE1Bins : 1) * binParams->energy1CIsize;
				binParams->weightImageSize = ((binParams->numE1Bins > 0) ? binParams->numE1Bins : 1) * binParams->energy1WIsize;
				binParams->weightSquImageSize = ((binParams->numE1Bins > 0) ? binParams->numE1Bins : 1) * binParams->energy1WISsize;
				break;
		
			case PhgBinEn_Energy2:
				binParams->countImageSize = ((binParams->numE2Bins > 0) ? binParams->numE2Bins : 1) * binParams->energy2CIsize;
				binParams->weightImageSize = ((binParams->numE2Bins > 0) ? binParams->numE2Bins : 1) * binParams->energy2WIsize;
				binParams->weightSquImageSize = ((binParams->numE2Bins > 0) ? binParams->numE2Bins : 1) * binParams->energy2WISsize;
				break;
	
			case PhgBinEn_Crystal1:
				binParams->countImageSize = ((binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1) * binParams->crystal1CIsize;
				binParams->weightImageSize = ((binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1) * binParams->crystal1WIsize;
				binParams->weightSquImageSize = ((binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1) * binParams->crystal1WISsize;
				break;
		
			case PhgBinEn_Crystal2:
				binParams->countImageSize = ((binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1) * binParams->crystal2CIsize;
				binParams->weightImageSize = ((binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1) * binParams->crystal2WIsize;
				binParams->weightSquImageSize = ((binParams->numCrystalBins > 0) ? binParams->numCrystalBins : 1) * binParams->crystal2WISsize;
				break;
	
					
			case PhgBinEn_Scatter1:
				binParams->countImageSize = ((binParams->numS1Bins > 0) ? binParams->numS1Bins : 1) * binParams->scatter1CIsize;
				binParams->weightImageSize = ((binParams->numS1Bins > 0) ? binParams->numS1Bins : 1) * binParams->scatter1WIsize;
				binParams->weightSquImageSize = ((binParams->numS1Bins > 0) ? binParams->numS1Bins : 1) * binParams->scatter1WISsize;
				break;
	
			case PhgBinEn_Scatter2:
				binParams->countImageSize = ((binParams->numS2Bins > 0) ? binParams->numS2Bins : 1) * binParams->scatter2CIsize;
				binParams->weightImageSize = ((binParams->numS2Bins > 0) ? binParams->numS2Bins : 1)  * binParams->scatter2WIsize;
				binParams->weightSquImageSize = ((binParams->numS2Bins > 0) ? binParams->numS2Bins : 1)  * binParams->scatter2WISsize;
				break;
		
			case PhgBinEn_Z1:
				binParams->countImageSize = ((binParams->numZBins > 0) ? binParams->numZBins : 1) * binParams->z1CIsize;
				binParams->weightImageSize = ((binParams->numZBins > 0) ? binParams->numZBins : 1)  * binParams->z1WIsize;
				binParams->weightSquImageSize = ((binParams->numZBins > 0) ? binParams->numZBins : 1)  * binParams->z1WISsize;
				break;

			case PhgBinEn_Z2:
				binParams->countImageSize = ((binParams->numZBins > 0) ? binParams->numZBins : 1)  * binParams->z2CIsize;
				binParams->weightImageSize = ((binParams->numZBins > 0) ? binParams->numZBins : 1)  * binParams->z2WIsize;
				binParams->weightSquImageSize = ((binParams->numZBins > 0) ? binParams->numZBins : 1)  * binParams->z2WISsize;
				break;
				

			case PhgBinEn_THETA:
				binParams->countImageSize = ((binParams->numThetaBins > 0) ? binParams->numThetaBins : 1) * binParams->thetaCIsize;
				binParams->weightImageSize = ((binParams->numThetaBins > 0) ? binParams->numThetaBins : 1)  * binParams->thetaWIsize;
				binParams->weightSquImageSize = ((binParams->numThetaBins > 0) ? binParams->numThetaBins : 1)  * binParams->thetaWISsize;
				break;

			case PhgBinEn_PHI:
				binParams->countImageSize = ((binParams->numPHIBins > 0) ? binParams->numPHIBins : 1) * binParams->phiCIsize;
				binParams->weightImageSize = ((binParams->numPHIBins > 0) ? binParams->numPHIBins : 1) * binParams->phiWIsize;
				binParams->weightSquImageSize = ((binParams->numPHIBins > 0) ? binParams->numPHIBins : 1) * binParams->phiWISsize;
				break;
				

			case PhgBinEn_XR:
				binParams->countImageSize = ((binParams->numXRBins > 0) ? binParams->numXRBins : 1) * binParams->xrCIsize;
				binParams->weightImageSize = ((binParams->numXRBins > 0) ? binParams->numXRBins : 1) * binParams->xrWIsize;
				binParams->weightSquImageSize = ((binParams->numXRBins > 0) ? binParams->numXRBins : 1) * binParams->xrWISsize;
				break;
				

			case PhgBinEn_YR:
				binParams->countImageSize = ((binParams->numYRBins > 0) ? binParams->numYRBins : 1) * binParams->yrCIsize;
				binParams->weightImageSize = ((binParams->numYRBins > 0) ? binParams->numYRBins : 1) * binParams->yrWIsize;
				binParams->weightSquImageSize = ((binParams->numYRBins > 0) ? binParams->numYRBins : 1) * binParams->yrWISsize;
				break;
				
			
			default:
				ErStGeneric("Invalid binning dimension in array (PhgBinInitialize)");
				goto FAIL;
		}
		
		/* Finally, the image size is adjusted to the size of the data structure that is
			being used to store it.
		*/
		if ((binParams->sumAccordingToType == true) && (binParams->weight_image_type == PHG_BIN_WEIGHT_TYPE_R4)){
			binParams->weightImageSize *= sizeof(float);
			binParams->weightSquImageSize *= sizeof(float);
		}
		else {
			binParams->weightImageSize *= sizeof(double);
			binParams->weightSquImageSize *= sizeof(double);
		}
		
		/* Currently the countImageSize contains the number of bins, regardless of if we
			are doing count images or not. So we'll set our numImageBins according to this.
		*/
		binParams->numImageBins = binParams->countImageSize;

		/* Now adjust size of count image based on the size of the bins */
		if (binParams->count_image_type == PHG_BIN_COUNT_TYPE_I1){
			binParams->countImageSize *= sizeof(LbUsOneByte);
		}
		else if (binParams->count_image_type == PHG_BIN_COUNT_TYPE_I2){
			binParams->countImageSize *= sizeof(LbUsTwoByte);
		}
		else if (binParams->count_image_type == PHG_BIN_COUNT_TYPE_I4){
			binParams->countImageSize *= sizeof(LbUsFourByte);
		}
			
		binFields->ParamsIsInitialized = true;
		FAIL:;
	} while (false);
	
	return (binFields->ParamsIsInitialized);
}

/*********************************************************************************
*
*			Name:			PhgBinInitialize
*
*			Summary:		Initialize the binning module.
*				char			*paramsName	 - The name of the parameter file.
*				PHG_BinParamsTy	*binParams	 - User defined binning parameters.
*				PHG_BiDataTy	*binData	 - Storage for binned data.
*				PHG_BinFieldsTy *binFields	 - Various binning information.
*			Arguments:
*			Function return: TRUE unless an error occurs.
*
*********************************************************************************/
Boolean PhgBinInitialize(char *paramsName, PHG_BinParamsTy *binParams, PHG_BinDataTy *binData, PHG_BinFieldsTy *binFields)
{
	double			weightRatio;		/* Ratio for scaling of pre-existing images */
	double			weightSquRatio;		/* Ratio for scaling of pre-existing images */
	LbUsEightByte	imageIndex;			/* Index for processing existing files */
	Boolean			cleared;			/* Used for checking error status */
	do {
	
		/* Clear processing  variables */
		binFields->NumCoincidences = 0;
		binFields->NumAcceptedCoincidences = 0;
		binFields->TotBluePhotons = 0;
		binFields->TotPinkPhotons = 0;
		binFields->AccBluePhotons = 0;
		binFields->AccPinkPhotons = 0;
		binFields->AccCoincidenceWeight = 0;
		binFields->AccCoincidenceSquWeight = 0;
		binFields->StartAccCoincidenceWeight = 0;
		binFields->StartAccCoincidenceSquWeight = 0;

		LbHdrStNull(&binFields->CountImgHdrHk);
		binFields->CountFile = 0;
		binData->countImage = 0;
		LbHdrStNull(&binFields->WeightSquImgHdrHk);
		binFields->WeightSquFile = 0;
		binData->weightSquImage = 0;
		LbHdrStNull(&binFields->WeightImgHdrHk);
		binFields->WeightFile = 0;
		binData->weightImage = 0;
		
		/* Initialize parameters */
		if (PhgBinInitParams(paramsName, binParams, binData, binFields) == false)
			break;

		/* Create history file if were are supposed to */
		if (binParams->isHistoryFile) {
			if (PhoHFileCreate(binParams->history_path,
					binParams->history_params_path, PhoHFileEn_BIN, &binFields->historyFileHk) == false) {
					
				sprintf(detErrStr,"Failed to create history file specified in binning parameters file named:\n"
					"'%s'\n"
					"The custom parameters file name is: '%s'\n"
					" (PhgBinInitialize)",
					binParams->history_path, binParams->history_params_path);
				PhgAbort(phgBinErrStr, false);
			}
		}

		/* Allocate image buffers */
		if (binParams->doCounts == true) {
		
			/* Allocate the image buffer */
			if ((binData->countImage = LbMmAlloc8(binParams->countImageSize))
					== 0) {
				break;
			}			
		}
		if (binParams->doWeights == true) {
		
			/* Allocate the image buffer */
			if ((binData->weightImage = LbMmAlloc8(binParams->weightImageSize))
					== 0) {
				break;
			}
		}
		if (binParams->doWeightsSquared == true) {
		
			/* Allocate the image buffer */
			if ((binData->weightSquImage = LbMmAlloc8(binParams->weightSquImageSize))
					== 0) {
				break;
			}
		}
		
		/* We'll clear the ratio's now, they'll be changed later on if necessary */	
		binFields->WeightRatio = 1.0;
				

		/* See if they are binning counts */
		if (binParams->doCounts == true) {
		
			/* Open the count image file */
			if (PhgBinOpenImage(binParams, binFields, PhoHFileEn_BIN_CT,
					binParams->countImgFilePath,
					&binFields->CountFile, &binFields->CountImgHdr,
					&binFields->CountImgHdrHk,
					binParams->countImageSize,
					(void *)binData->countImage) == false) {
				break;
			}
			
			/* Update our count of accepted photons from total in histogram */
			if (binFields->CountImgHdr.H.NumSimulations >= 1) {

				/* Update header fields */
				binFields->CountImgHdr.H.SumEventsToSimulate += PhgRunTimeParams.Phg_EventsToSimulate;
				binFields->CountImgHdr.H.NumSimulations++;

				/* Compute our ratios, this should be consistent each time this routine is called.
					Some good error checking should be performed when doing this step to insure the
					user is not adding to existing images that come from different simulations and
					so on; but its not happening yet.
				*/
				binFields->WeightRatio = (double) PhgRunTimeParams.Phg_EventsToSimulate/
									(double)binFields->CountImgHdr.H.SumEventsToSimulate;

				if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
					binFields->NumAcceptedCoincidences = binFields->CountImgHdr.H.NumPhotons;
				}
				else {
					binFields->AccBluePhotons = binFields->CountImgHdr.H.NumPhotons;
				}
			}

		}

		/* See if they are binning weights */
		if (binParams->doWeights == true) {
						
			/* Open the weight image file */
			if (PhgBinOpenImage(binParams, binFields, PhoHFileEn_BIN_WT,
					binParams->weightImgFilePath,
					&binFields->WeightFile, &binFields->WeightImgHdr,
					&binFields->WeightImgHdrHk,
					(((binParams->sumAccordingToType == false) && (binParams->weight_image_type == PHG_BIN_WEIGHT_TYPE_R4))
						? binParams->weightImageSize/2 : binParams->weightImageSize),
					(void *)binData->weightImage) == false) {
				break;
			}


			/* If this is not the first simulation (addToExistingImages is true)
				we need to scale the existing weights by their proportional contribution.
			*/
			if (binFields->WeightImgHdr.H.NumSimulations >= 1) {

				/* Update header fields */
				binFields->WeightImgHdr.H.SumEventsToSimulate += PhgRunTimeParams.Phg_EventsToSimulate;
				binFields->WeightImgHdr.H.NumSimulations++;

				/* Compute our ratios, this should be consistent each time this routine is called.
					Some good error checking should be performed when doing this step to insure the
					user is not adding to existing images that come from different simulations and
					so on; but its not happening yet.
				*/
				binFields->WeightRatio = (double) PhgRunTimeParams.Phg_EventsToSimulate/
									(double)binFields->WeightImgHdr.H.SumEventsToSimulate;
			
				/* Compute ratios for existing data */
				weightRatio = (double) (binFields->WeightImgHdr.H.SumEventsToSimulate -
						PhgRunTimeParams.Phg_EventsToSimulate)/
						(double)binFields->WeightImgHdr.H.SumEventsToSimulate;
				
				/* Scale the weights */
				for (imageIndex = 0; imageIndex < binParams->numImageBins; imageIndex++) {
					
					if ((binParams->sumAccordingToType == true) && (binParams->weight_image_type == PHG_BIN_WEIGHT_TYPE_R4)) {

							/* Update the image */
							((float *)binData->weightImage)[imageIndex] *= weightRatio;
		
							/* Update our starting point for cumulative coincidence weight */
							binFields->AccCoincidenceWeight += ((float *)binData->weightImage)[imageIndex];
					}
					else {

							/* Update the image */
							((double *)binData->weightImage)[imageIndex] *= weightRatio;
		
							/* Update our starting point for cumulative coincidence weight */
							binFields->AccCoincidenceWeight += ((double *)binData->weightImage)[imageIndex];
					}

				}
				/* Save the starting weights */
				binFields->StartAccCoincidenceWeight = binFields->AccCoincidenceWeight;
				
				/* Get squared weights sum if there is no file to get it from (value may be zero but that is okay) */
				if (binParams->doWeightsSquared == false) {

					if (PhgHdrGtField(&(binFields->WeightImgHdrHk),  HDR_BIN_SUM_WEIGHTS_SQ_ID, (void *) &(binFields->AccCoincidenceSquWeight)) == false) {
						
						/* If the error is merely due to the field not existing, that is okay it is an older header
							and we will just set things to zero and warn the user
						*/
						ErClearIf(ERMgCdHeader, ERErCdHdElemNotFound, &cleared);
						if (!cleared) {
							PhgAbort("Unable to read weight squared value from header.", false);
						}
						else {
							ErAlert("There is no previous weight squared sum, so be aware that this value is starting out zero.", false);
						}
						binFields->AccCoincidenceSquWeight = 0.0;
					}
					/* Save the starting weights squared */
					binFields->StartAccCoincidenceSquWeight = binFields->AccCoincidenceSquWeight;

				}
				
			}
		}

		/* See if they are binning weights squared */
		if (binParams->doWeightsSquared == true) {
			
			/* Open the weight squared image file */
			if (PhgBinOpenImage(binParams, binFields, PhoHFileEn_BIN_WTSQ,
					binParams->weightSquImgFilePath,
					&binFields->WeightSquFile, &binFields->WeightSquImgHdr,
					&binFields->WeightSquImgHdrHk,
							(((binParams->sumAccordingToType == false) && (binParams->weight_image_type == PHG_BIN_WEIGHT_TYPE_R4))
						? binParams->weightSquImageSize/2 : binParams->weightSquImageSize),
					(void *)binData->weightSquImage) == false) {
				break;
			}
			/* If this is not the first simulation (addToExistingImages is true)
				we need to scale the existing weights by their proportional contribution.
			*/
			if (binFields->WeightSquImgHdr.H.NumSimulations >= 1) {

				/* Update header fields */
				binFields->WeightSquImgHdr.H.SumEventsToSimulate += PhgRunTimeParams.Phg_EventsToSimulate;
				binFields->WeightSquImgHdr.H.NumSimulations++;

				/* Compute our ratios, this should be consistent each time this routine is called.
					Some good error checking should be performed when doing this step to insure the
					user is not adding to existing images that come from different simulations and
					so on; but its not happening yet.
				*/
				binFields->WeightRatio = (double) PhgRunTimeParams.Phg_EventsToSimulate/
									(double)binFields->WeightSquImgHdr.H.SumEventsToSimulate;
			
				/* Compute ratios for existing data */
				weightRatio = (double) (binFields->WeightSquImgHdr.H.SumEventsToSimulate -
						PhgRunTimeParams.Phg_EventsToSimulate)/
						(double)binFields->WeightImgHdr.H.SumEventsToSimulate;
				
				weightSquRatio = weightRatio   * weightRatio;
				
				/* Scale the weights squared */
				for (imageIndex = 0; imageIndex < binParams->numImageBins; imageIndex++) {
					
					if ((binParams->sumAccordingToType == true) && (binParams->weight_image_type == PHG_BIN_WEIGHT_TYPE_R4)) {

						/* Update the image */
						((float *)binData->weightSquImage)[imageIndex] *= 
							weightSquRatio;
	
						/* Update our starting point for cumulative coincidence weight */
						binFields->AccCoincidenceSquWeight += ((float *)binData->weightSquImage)[imageIndex];
	
					}
					else {

						/* Update the image */
						((double *)binData->weightSquImage)[imageIndex] *= 
							weightSquRatio;
	
						/* Update our starting point for cumulative coincidence weight */
						binFields->AccCoincidenceSquWeight += ((double *)binData->weightSquImage)[imageIndex];
					}
				}
				/* Save the starting weights squared */
				binFields->StartAccCoincidenceSquWeight = binFields->AccCoincidenceSquWeight;

			}

		}
		
		/* Call the user binning routine */
		if (BinUsrInitializeFPtr) {
			(*BinUsrInitializeFPtr)(binParams, binData);
		}
		
		binFields->IsInitialized = true;
		FAIL:;
	} while (false);
	
	/* If we failed, attempt to free memory that might have been allocated */
	if (!binFields->IsInitialized) {
		if (binData->countImage != 0) {
			LbMmFree(&(binData->countImage));
		}
		if (binData->weightImage != 0) {
			LbMmFree(&(binData->weightImage));
		}
		if (binData->weightSquImage != 0) {
			LbMmFree(&(binData->weightSquImage));
		}
	}
	
	return (binFields->IsInitialized);
}

/*********************************************************************************
*
*			Name:			PhgBinPrintParams
*
*			Summary:		Print the binning parameters.
*				char			*paramsName	- Name of parameters.
*				PHG_BinParamsTy	*binParams	- User defined binning parameters.
*				PHG_BinFieldsTy	*binFields	- Fields for binning
*
*			Arguments:
*			Function return: None.
*
*********************************************************************************/
void PhgBinPrintParams(char *paramsName, PHG_BinParamsTy *binParams,
		PHG_BinFieldsTy	*binFields)
{
	LbFourByte	dimensionIndex;		/* Index for processing dimension calculations */


	LbInPrintf("\n\n Binning Parameters for '%s'\n", paramsName);
	
	if (binParams->doMSRB) {
		LbInPrintf("\n\t Multi-Slice Rebinning is being performed\n");
		LbInPrintf("\t\tObject radius = %3.2f\tDetector radius = %3.2f\n", binParams->image_radius,
			binParams->detector_radius);
	}
	if (binParams->doSSRB) {
		LbInPrintf("\n\t Single-Slice Rebinning is being performed\n");
	}
	
	LbInPrintf("\n\t Binning parameters are ordered by:\n\t slowest varying to fastest (top of this list to bottom).\n");
	
	/* Loop through each dimension */
	for (dimensionIndex = (PHGBIN_NUM_DIMENSIONS-1); dimensionIndex >= (PHGBIN_NUM_DIMENSIONS-binParams->numDimensions); dimensionIndex--){

		/* We can now compute the sizes of the images, based on last dimension  */
		switch(binParams->PhgBinDimensions[dimensionIndex]) {
		
			case PhgBinEn_TD:
				if (binParams->numTDBins >= 1) {
					LbInPrintf("\n\t  Number of transaxial dist. bins = %ld, range = [%3.2f, %3.2f] (cm)",
						(unsigned long)binParams->numTDBins, binParams->minTD, binParams->maxTD);
				}
				break;
		
			case PhgBinEn_TOF:
				if (binParams->numTOFBins >= 1) {
					LbInPrintf("\n\t  Number of time-of-flight bins = %ld, range = [%3.2f, %3.2f] (nanosecs)",
						(unsigned long)binParams->numTOFBins, binParams->minTOF, binParams->maxTOF);
				}
				break;
		
			case PhgBinEn_PHI:
				if (binParams->numPHIBins >= 1) {
					LbInPrintf("\n\t  Number of 'PHI' bins = %ld, range = [%3.2f, %3.2f] (radians)",
						(unsigned long)binParams->numPHIBins,
						binParams->minPhi, binParams->maxPhi);
				}
				break;
		
			case PhgBinEn_THETA:
				if (binParams->numThetaBins >= 1) {
					LbInPrintf("\n\t  Number of 'Theta' bins = %ld, range = [%3.2f, %3.2f] (radians)",
						(unsigned long)binParams->numThetaBins,
						binParams->minTheta, binParams->maxTheta);
				}
				break;
		
			case PhgBinEn_XR:
				if (binParams->numXRBins >= 1) {
					LbInPrintf("\n\t  Number of Xr bins = %ld, range = [%3.2f, %3.2f] (cm)",
						(unsigned long)binParams->numXRBins, binParams->minXR, binParams->maxXR);
				}
				break;
		
			case PhgBinEn_YR:
				if (binParams->numYRBins >= 1) {
					LbInPrintf("\n\t  Number of Yr bins = %ld, range = [%3.2f, %3.2f] (cm)",
						(unsigned long)binParams->numYRBins, binParams->minYR, binParams->maxYR);
				}
				break;
		
			case PhgBinEn_AA:
				if (binParams->numAABins >= 1) {
					LbInPrintf("\n\t  Number of azimuthal angle bins = %ld, range = [%3.2f, %3.2f] (radians)",
						(unsigned long)binParams->numAABins, binParams->minAA, binParams->maxAA);
				}
				break;
	
			case PhgBinEn_Crystal1:
				
				if (binParams->numCrystalBins >= 1) {
	
					/* Print crystal parameters based on type of simulation */
					if (PHG_IsSPECT() || binParams->isBinPETasSPECT){
						LbInPrintf("\n\t  Number of crystal bins = %ld",
							(unsigned long)binParams->numCrystalBins);
					}
					else {
						LbInPrintf("\n\t  Crystal vs. crystal binning (binned in upper triangle):");
						LbInPrintf("\n\t  Number of crystal*crystal bins = %ld * %ld",
							(unsigned long)binParams->numCrystalBins, 
							(unsigned long)binParams->numCrystalBins);
					}
				}
				break;

			case PhgBinEn_Crystal2:
				break;
				
			case PhgBinEn_Energy1:
				
				if ((binParams->numE1Bins >= 1) || (binParams->numE2Bins >= 1)) {

					/* Print energy parameters based on type of simulation */
					if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ){
						LbInPrintf("\n\t  Number of energy (1) bins = %ld, range = [%3.2f, %3.2f] (keV)",
							(unsigned long)binParams->numE1Bins, binParams->minE, binParams->maxE);
		
						LbInPrintf("\n\t  Number of energy (2) bins = %ld, range = [%3.2f, %3.2f] (keV)",
							(unsigned long)binParams->numE2Bins, binParams->minE, binParams->maxE);
					}
					else if (binParams->numE1Bins >= 1) {
						LbInPrintf("\n\t  Number of energy bins = %ld, range = [%3.2f, %3.2f] (keV)",
							(unsigned long)binParams->numE1Bins, binParams->minE, binParams->maxE);
					}
				}
				break;
				
			case PhgBinEn_Energy2:
				break;
	
					
			case PhgBinEn_Scatter1:
				
				if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
					
					if ( binParams->acceptRandoms ) {
						LbInPrintf("\n\t  Random events are being accepted.\n");
					} else {
						LbInPrintf("\n\t  Random events are being rejected.\n");
					}
					
				}
				
				/* scatter/random parameter = 0 */
				if (binParams->scatterRandomParam == 0) {
					if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
						
						if ( binParams->acceptRandoms ) {
							LbInPrintf("\n\t  All scatter/random states (primary, scatter, random) stored together.");
							LbInPrintf("\n\t  Number of scatter/random bins = %ld, range = [%ld, %ld] (# of scatters)",
								(unsigned long)binParams->numS1Bins, 
								(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
						} else {
							LbInPrintf("\n\t  All scatter states (primary, scatter) stored together.");
							LbInPrintf("\n\t  Number of scatter bins = %ld, range = [%ld, %ld] (# of scatters)",
								(unsigned long)binParams->numS1Bins, 
								(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
						}

					}
					else {
						LbInPrintf("\n\t  All scatter states (primary, scatter) stored together.");
						LbInPrintf("\n\t  Number of scatter bins = %ld, range = [%ld, %ld] (# of scatters)",
							(unsigned long)binParams->numS1Bins, 
							(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
					}
				}

				/* scatter/random parameter = 1 */
				else if (binParams->scatterRandomParam == 1){
					LbInPrintf("\n\t  There are two scatter bins,\n\t    the first contains unscattered events,\n"
						"\t    the second contains events whose scatter count\n\t    falls in the range [%ld, %ld]",
						(unsigned long)PHGMATH_Max(binParams->minS, 1), 
						(unsigned long)binParams->maxS);
				}

				/* scatter/random parameter = 2 */
				else if (binParams->scatterRandomParam == 2){
					if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
						LbInPrintf("\n\t  Scatter state is stored as two-dimensions, ");
						LbInPrintf("\n\t  one for the number of scatters of each photon. ");
						LbInPrintf("\n\t  Number of scatter bins (photon 1) = %ld, range = [%ld, %ld] (# of scatters)",
							(unsigned long)binParams->numS1Bins, 
							(unsigned long)binParams->minS, (unsigned long)binParams->maxS);

						LbInPrintf("\n\t  Number of scatter bins (photon 2) = %ld, range = [%ld, %ld] (# of scatters)",
							(unsigned long)binParams->numS2Bins, 
							(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
					}
					else {
						LbInPrintf("\n\t  Number of scatter bins = %ld, range = [%ld, %ld] (# of scatters)",
							(unsigned long)binParams->numS1Bins, 
							(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
					}
				}

				/* scatter/random parameter = 3 */
				else if (binParams->scatterRandomParam == 3){
					if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
						LbInPrintf("\n\t  Scatter state is stored as two-dimensions, ");
						LbInPrintf("\n\t  one for the number of scatters of each photon. ");
						LbInPrintf("\n\t  Photons with number of scatters > maximum");
						LbInPrintf("\n\t  are stored with the maximum number of scatters.");
						LbInPrintf("\n\t  Number of scatter bins (photon 1) = %ld, range = [%ld, >= %ld] (# of scatters)",
							(unsigned long)binParams->numS1Bins, 
							(unsigned long)binParams->minS, (unsigned long)binParams->maxS);

						LbInPrintf("\n\t  Number of scatter bins (photon 2) = %ld, range = [%ld, >= %ld] (# of scatters)",
							(unsigned long)binParams->numS2Bins, 
							(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
					}
					else {
						LbInPrintf("\n\t  Number of scatter bins = %ld, range = [%ld, >= %ld] (# of scatters)",
							(unsigned long)binParams->numS1Bins, 
							(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
					}
				}

				/* scatter/random parameter = 4 */
				else if (binParams->scatterRandomParam == 4){
					LbInPrintf("\n\t  Scatter state is stored as one-dimension, as the ");
					LbInPrintf("\n\t  sum of the number of scatters for the two photons. ");
					LbInPrintf("\n\t  Number of scatter bins = %ld, range = [%ld, %ld] (# of scatters)",
						(unsigned long)binParams->numS1Bins, 
						(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
				}

				/* scatter/random parameter = 5 */
				else if (binParams->scatterRandomParam == 5){
					LbInPrintf("\n\t  Scatter state is stored as one-dimension, as the ");
					LbInPrintf("\n\t  sum of the number of scatters for the two photons. ");
					LbInPrintf("\n\t  Events with sum of scatters > maximum");
					LbInPrintf("\n\t  are stored with the maximum number of scatters.");
					LbInPrintf("\n\t  Number of scatter bins = %ld, range = [%ld, >= %ld] (# of scatters)",
						(unsigned long)binParams->numS1Bins, 
						(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
				}

				/* scatter/random parameter = 6 */
				else if (binParams->scatterRandomParam == 6){
					LbInPrintf("\n\t  There are three scatter/random bins,\n\t"
						"    the first contains primary (true) events,\n"
						"\t    the second contains non-random scatter coincidences "
						"\n\t   whose scatter count falls in the range [%ld, %ld]",
						"\n\t   the third contains all random coincidences.",
						(unsigned long)PHGMATH_Max(binParams->minS, 1), 
						(unsigned long)binParams->maxS);
				}

				/* scatter/random parameter = 7 */
				else if (binParams->scatterRandomParam == 7){
					LbInPrintf("\n\t  Number of scatter/random state bins = %ld",
						(unsigned long)binParams->numS1Bins);
					LbInPrintf("\n\t  Non-random coincidences with number of scatters ");
					LbInPrintf("\n\t  in range = [%ld, %ld] for each photon are binned in ",
						(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
					LbInPrintf("\n\t  the first %ld bins, with bin number ",
						(unsigned long)(binParams->numS1Bins - 1) );
					LbInPrintf("\n\t  (number scatters photon 1) * (number of scatters photon 2).");
					LbInPrintf("\n\t  All random coincidences are stored in the last bin. ");
				}

				/* scatter/random parameter = 8 */
				else if (binParams->scatterRandomParam == 8){
					LbInPrintf("\n\t  Number of scatter/random state bins = %ld",
						(unsigned long)binParams->numS1Bins);
					LbInPrintf("\n\t  Non-random coincidences with number of scatters ");
					LbInPrintf("\n\t  in range = [%ld, %ld] for each photon are binned in ",
						(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
					LbInPrintf("\n\t  the first %ld bins, with bin number = ",
						(unsigned long)(binParams->numS1Bins - 1) );
					LbInPrintf("\n\t  (number scatters photon 1) * (number of scatters photon 2).");
					LbInPrintf("\n\t  Photons with number of scatters > %ld",
						(unsigned long)binParams->maxS);
					LbInPrintf("\n\t  are treated as though they had scattered %ld times.",
						(unsigned long)binParams->maxS);
					LbInPrintf("\n\t  All random coincidences are stored in the last bin. ");
				}

				/* scatter/random parameter = 9 */
				else if (binParams->scatterRandomParam == 9){
					LbInPrintf("\n\t  Scatter state is stored as one-dimension, as the ");
					LbInPrintf("\n\t  sum of the number of scatters for the two photons. ");
					LbInPrintf("\n\t  Number of scatter bins = %ld, range = [%ld, %ld] (# of scatters)",
						(unsigned long)binParams->numS1Bins, 
						(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
				}

				/* scatter/random parameter = 10 */
				else if (binParams->scatterRandomParam == 10){
					LbInPrintf("\n\t  Scatter state is stored as one-dimension, as the ");
					LbInPrintf("\n\t  sum of the number of scatters for the two photons. ");
					LbInPrintf("\n\t  Events with sum of scatters > maximum");
					LbInPrintf("\n\t  are stored with the maximum number of scatters.");
					LbInPrintf("\n\t  Number of scatter bins = %ld, range = [%ld, >= %ld] (# of scatters)",
						(unsigned long)binParams->numS1Bins, 
						(unsigned long)binParams->minS, (unsigned long)binParams->maxS);
				}

				break;
	
			case PhgBinEn_Scatter2:
				break;
		
			case PhgBinEn_Z1:
				
				if (binParams->numZBins >= 1) {
	
					/* Print z parameters based on type of simulation */
					if (PHG_IsPET() && !(binParams->doSSRB) && !(binParams->doMSRB)){
						LbInPrintf("\n\t  Number of Z axis (1) bins = %ld, range [%3.2f, %3.2f] (cm)",
							(unsigned long)binParams->numZBins, binParams->minZ, binParams->maxZ);
		
						LbInPrintf("\n\t  Number of Z axis (2) bins = %ld, range [%3.2f, %3.2f] (cm)",
							(unsigned long)binParams->numZBins, binParams->minZ, binParams->maxZ);
					}
					else {
						LbInPrintf("\n\t  Number of Z axis bins = %ld, range [%3.2f, %3.2f] (cm)",
							(unsigned long)binParams->numZBins, binParams->minZ, binParams->maxZ);
					}
				}
				break;

			case PhgBinEn_Z2:
				break;
				
				default:
					break;
		}
	} /* End of for-each-dimension */
	
	/* Indicate if files are being added to or created from scratch */
	LbInPrintf("\n\n\t  Add to existing images is %s\n", ((binParams->addToExistingImg == true) ? "true." :
		"false."));
	
	/* Indicate if files are being summed as doubles or according to output type */
	if (binParams->sumAccordingToType == false) {
		LbInPrintf("\n\t  Weight data is being summed in double precision");
	}
	else {
		LbInPrintf("\n\t  Weight data is being summed in precision of weight image format");
	}
	
	/* Print out image file paths and formats */
	if (binParams->doWeights == true) {
		switch(binParams->weight_image_type){
				
			case PHG_BIN_WEIGHT_TYPE_R4:
	
				LbInPrintf("\n\t  Weight image format is four byte reals.");
				break;
	
			case PHG_BIN_WEIGHT_TYPE_R8:
	
				LbInPrintf("\n\t  Weight image format is eight byte reals.");
				break;
		}
	
		LbInPrintf("\n\t  Weight image file name is '%s'\n",
			binParams->weightImgFilePath);
	}
	
	if (binParams->doWeightsSquared == true) {
		switch(binParams->weight_image_type){
				
			case PHG_BIN_WEIGHT_TYPE_R4:
	
				LbInPrintf("\n\t  Weight squared image format is four byte reals.");
				break;
	
			case PHG_BIN_WEIGHT_TYPE_R8:
	
				LbInPrintf("\n\t  Weight squared image format is eight byte reals.");
				break;
		}
	
		LbInPrintf("\n\t  Weight squared image file name is '%s'\n",
			binParams->weightSquImgFilePath);
	}

	/* Print out image file paths and formats */
	if (binParams->doCounts == true) {
		switch(binParams->count_image_type){
				
			case PHG_BIN_COUNT_TYPE_I1:
	
				LbInPrintf("\n\t  Count image format is one byte integers.");
				break;
				
			case PHG_BIN_COUNT_TYPE_I2:
	
				LbInPrintf("\n\t  Count image format is two byte integers.");
				break;
				
			case PHG_BIN_COUNT_TYPE_I4:
	
				LbInPrintf("\n\t  Count image format is four byte integers.");
				break;
		}
	
		LbInPrintf("\n\t  Count image file name is '%s'",
			binParams->countImgFilePath);
	}
			
	/* Print out history parameters */
	if (binFields->historyFileHk.doCustom) {
		LbInPrintf("\nHistory file parameters for binning module");
		PhoHFilePrintParams(&binFields->historyFileHk);
	}
}

/*********************************************************************************
*
*			Name:			PhgBinOpenImage
*
*			Summary:		Open the image.
*
*			Arguments:
*				PHG_BinParamsTy		*binParams		- Specific binning parameters
*				PHG_BinFieldsTy 	*binFields		- Various binning information.
*				PhoHFileHdrKindTy	hdrKind			- What type of header
*				FILE				**imageFile		- Handle to opened file
*				PhoHFileHdrTy		*headerPtr		- Header info for image file
*				LbHdrHkTy			*headerHkPtr	- Hook for header
*				LbUsFourByte		dataSize		- The size of the image data
*				void				*dataPtr		- Ptr to the data
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean PhgBinOpenImage(PHG_BinParamsTy *binParams, PHG_BinFieldsTy *binFields,
			PhoHFileHdrKindTy hdrKind, char *imageName,
			FILE **imageFile, PhoHFileHdrTy *headerPtr,
			LbHdrHkTy *headerHkPtr, LbUsFourByte dataSize, void *dataPtr)
{
	PHG_BinFieldsTy*	dummyPtr;			/* Removes compiler warning */
	Boolean			okay = false;			/* Process flag */
	Boolean			imageExists = false;	/* Does the image exist already */
	LbUsFourByte	dataStart;				/* Beginning of data in file */
	LbUsFourByte	dataEnd;				/* End of file position */
	LbUsFourByte	i;						/* Loop control */
	float			singlePrec;				/* For conversion */
	
	/* Avoid unused parameter compiler warning */
	dummyPtr = binFields;
	
	do {	/* Process Loop */
	
		/* See if the image exists already */
		if ((*imageFile = LbFlFileOpen(imageName, "rb")) == 0) {
				imageExists = false;
		}
		else {
			imageExists = true;
			fclose(*imageFile);
		}
		
		/* If the file exists and we are adding to existing images, we'll read in the header */
		if ((imageExists == true) && (binParams->addToExistingImg == true)) {

			/* Open the file for updating */
			if ((*imageFile = LbFlFileOpen(imageName, "r+b")) == 0) {
				ErStFileError("\nUnable to open image file");
				break;
			}
			
			/* Read in the header */
			if (PhgHdrGtParams(*imageFile, headerPtr, headerHkPtr) == false){
				break;
			}
			
			/* Verify the length of the scan is the same as the current one */
			if (headerPtr->H.PhgRunTimeParams.Phg_LengthOfScan != PHGGetLengthOfScan()) {
				ErStGeneric("The length of a pre-existing scan and the current on are different.  These can not be added together (PhgBinOpenImage).");
				break;
			}
			
			/* Verify that the data stored is the same size as what we
				are expecting to process
			*/
			{
				/* Get the current file position */
				dataStart = ftell(*imageFile);
				
				/* Seek to end of file */
				if (fseek(*imageFile, 0, SEEK_END) != 0){
					ErStFileError("Unable to seek to end of existing image (PhgBinOpenImage).");
					break;
				}
				
				/* Get current file position */
				dataEnd = ftell(*imageFile);
				
				/* Check size */
				if ((dataEnd-dataStart) != dataSize) {
					sprintf(phgBinErrStr,
						"Existing image data size (%ld) is not the same as the current image size specification (%ld)!",
							(unsigned long)(dataEnd - dataStart), (unsigned long)dataSize);
					ErStGeneric(phgBinErrStr);
					break;
				}
				
				/* Seek back to beginning of data */
				if (fseek(*imageFile, dataStart, SEEK_SET) != 0){
					ErStFileError("Unable to seek to beginning of image data: (PhgBinOpenImage).");
					break;
				}
				
			}
			/* Read the existing image data 
				NOTE: the data must be converted if it was written as single precision and we are summing in double
			*/
			if ((binParams->sumAccordingToType == false) && (headerPtr->H.BinRunTimeParams.weight_image_type == PHG_BIN_WEIGHT_TYPE_R4)
					&& (hdrKind != PhoHFileEn_BIN_CT)) {
				/* Tell the user this may take a while */
				LbInPrintf("\nReading existing image and converting from single to double precision, this may take a while\n");
				
				for (i = 0; i < dataSize/sizeof(float); i++) {
					if (fread(&singlePrec, sizeof(float), 1, *imageFile) != 1) {
						ErStFileError("Unable to read existing image data: (PhgBinOpenImage).");
						break;
					}
					((double *)dataPtr)[i] = singlePrec;
				}
			}
			else {
				if (fread(dataPtr, dataSize, 1, *imageFile) != 1) {
					ErStFileError("Unable to read existing image data: (PhgBinOpenImage).");
					break;
				}
			}
							
		}
		else  {
		
			/* Just create/open the file */
			if ((*imageFile = LbFlFileOpen(imageName, "wb")) == 0) {
				ErStFileError("Unable to create/open image file");
				break;
			}
			
			/* Reserve file space by seeking to the "will-be" end of file (-1), and writing a byte */
			if (fseek(*imageFile, (dataSize+PHG_HDR_HEADER_SIZE)-1, SEEK_SET) != 0){
					ErStFileError("Unable to seek to end of file for reserving space: (PhgBinOpenImage).");
					break;
			}
				
			/* Write a null value to allocate the image space */
			{
				char emptyValue = 0;
				
				if (fwrite(&emptyValue, 1, 1, *imageFile) != 1) {
					ErStFileError("Unable to write end of file marker: (PhgBinOpenImage).");
					break;
				}		
			}
			
			/* File is newly created so generate an empty header now */
			if (PhgHdrCreateRunTimeHeader(hdrKind, headerPtr, binParams) == false) {
				break;
			}
			
			/* Create a new header hook for the file */
			if (PhgHdrMkHeader(*imageFile, headerPtr, headerHkPtr) == false){
				break;
			}
		}


		okay = true;
	} while (false);
	

	return(okay);
}


/*********************************************************************************
*
*			Name:			PhgBinPETPhotons
*
*			Summary:		Update the binning images with the current batch of
*							detected photons.
*
*			Arguments:
*				PHG_BinParamsTy		*binParams		- User defined binning parameters.
*				PHG_BiDataTy		*binData	 	- Storage for binned data.
*				PHG_BinFieldsTy *binFields			 - Various binning information.
*				PHG_Decay			decay			- The decay that started the process.
*				PHG_TrackingPhoton *bluePhotons		- The blue photons detected.
*				LbUsFourByte 		numBluePhotons	- The number of blue photons.
*				PHG_TrackingPhoton *pinkPhotons		- The  pink photons detected.
*				LbUsFourByte		numPinkPhotons	- The number of pink photons.
*			Function return: None.
*
*********************************************************************************/
void PhgBinPETPhotons(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData, PHG_BinFieldsTy *binFields,
		PHG_Decay *decay,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons)
{
	double				distance = 0;				/* Computed transaxial distance */
	double				angle = 0;					/* Computed azimuthal angle */
	double				tofDifference;			/* difference in time-of-flight for two photons */
	double				coincidenceWeight;		/* Computed coincidence weight */
	double				coincidenceSquWeight;	/* Computed coincidence weight squared */
	double				blueWeight;				/* Weight of blue photon */
	double				pinkWeight;				/* Weight of pink photon */
	double				phi;					/* Angle computed for 3DRP binning */
	double				d;						/* D computed in 3DRP binning */
	double				theta;					/* Theta computed for 3DRP binning */
	double				avgZ;					/* Average Z for SSRB and MSRB */
	double				lowZ;					/* Average Z for MSRB */
	double				highZ;					/* Average Z for MSRB */
	double				xr;						/* xr for MSRB */
	double				yr;						/* yr for MSRB */
	LbUsFourByte		lowZIndex = 0;				/* Index of low Z for MSRB */
	LbUsFourByte		highZIndex = 0;				/* Index of high Z for MSRB */
	LbUsFourByte		zIndex = 0;					/* LCV for MSRB */
	LbUsFourByte		blueIndex;				/* Current blue photon */
	LbUsFourByte		pinkIndex;				/* Current pink index */
	LbUsFourByte		blueScatters;			/* Number of scatters for current blue photon */
	LbUsFourByte		pinkScatters;			/* Number of scatters for current pink photon */
	Boolean				ignoreMaxScatters;		/* set true if the scatter/random param is 3,5,8 or 10 */
	Boolean				ignoreMinScatters;		/* set true if the scatter/random param is 4,5,9 or 10 */
	LbUsFourByte		angleIndex = 0;				/* Index for angle bin */
	LbUsFourByte		distIndex = 0;				/* Index for dist bin */
	LbUsFourByte		tofIndex = 0;				/* Index for time-of-flight bin */
	LbUsFourByte		scatter1Index = 0;			/* Index for scatter 1 bin */
	LbUsFourByte		scatter2Index = 0;			/* Index for scatter 2 bin */
	LbUsFourByte		crystal1Index = 0;			/* Index for crystal 1 bin */
	LbUsFourByte		crystal2Index = 0;			/* Index for crystal 2 bin */
	LbUsFourByte		energy1Index = 0;			/* Index for energy 1 bin */
	LbUsFourByte		energy2Index = 0;			/* Index for energy 2 bin */
	LbUsFourByte		zDownIndex = 0;				/* Z index for photon having min(blue.y,pink.y) */
	LbUsFourByte		zUpIndex = 0;				/* Z index for photon having max(blue.y,pink.y) */
	LbUsFourByte		*z1IndexPtr = 0;			/* Swapping ptr for doing up/down stuff */
	LbUsFourByte		*z2IndexPtr = 0;			/* Swapping ptr for doing up/down stuff */
	LbUsFourByte		thetaIndex = 0;				/* Index for theta in 3DRP */
	LbUsFourByte		phiIndex = 0;				/* Index for PHI in 3DRP */
	LbUsFourByte		xrIndex = 0;				/* Xr index in 3DRP */
	LbUsFourByte		yrIndex = 0;				/* Yr index in 3DRP */
	LbUsEightByte		imageIndex = 0;				/* Index for image */
	double				flip;						/* For 3DRP */
	PHG_TrackingPhoton	bluePhoton;				/* Current blue photon */	
	PHG_TrackingPhoton	pinkPhoton;				/* Current pink photon */	

	/* Compute statistics */
	binFields->NumCoincidences += (numBluePhotons * numPinkPhotons);
	binFields->TotBluePhotons += numBluePhotons;
	binFields->TotPinkPhotons += numPinkPhotons;

	/* the maximum scatter range is ignored when 
	scatterRandomParam == 3,5,8 or 10 */
	ignoreMaxScatters = ( (binParams->scatterRandomParam == 3) ||
				(binParams->scatterRandomParam == 5) ||
				(binParams->scatterRandomParam == 8) ||
				(binParams->scatterRandomParam == 10) );
	
	/* for scatter-random parameters 4,5,9,10 the
	range min is for the sum of blue and pink scatters,
	and thus cannot be applied to individual photons */
	ignoreMinScatters = ( (binParams->scatterRandomParam == 4) ||
				(binParams->scatterRandomParam == 5) ||
				(binParams->scatterRandomParam == 9) ||
				(binParams->scatterRandomParam == 10) );
			
	/* Process coincidences */
	for (blueIndex = 0; blueIndex < numBluePhotons; blueIndex++) {
			

		for (pinkIndex = 0; pinkIndex < numPinkPhotons; pinkIndex++) {
			
			/* Init local blue/pink photon copies */
			bluePhoton = bluePhotons[blueIndex];
			pinkPhoton = pinkPhotons[pinkIndex];
			
			/* Let user modify and/or reject photons */
			if (BinUsrModPETPhotonsFPtr && 
					(*BinUsrModPETPhotonsFPtr)(binParams, binData,
						decay, &bluePhoton, &pinkPhoton) == false) {
				
				/* They rejected it so go to next pink */
				continue;
			}
			
			/* Reject the coincidence if it is random and 
			  we are rejecting randoms */
			if ( (decay->decayType == PhgEn_PETRandom) &&
					(!binParams->acceptRandoms) )
				continue;
			
			/* Get scatter count */
			blueScatters = bluePhoton.num_of_scatters +
				bluePhoton.scatters_in_col;
			
			/* See if this photon fails the acceptance criteria */
			/* NOTE: It seems like failing a blue should go to the
				next blue, not the next pink. However, because we let
				the user get their hands on the photon pair we will
				always call the user routine with every blue and
				every pink
			*/
			{
				if ( (blueScatters < binParams->minS) &&
						(!ignoreMinScatters) )
					continue;

				/* Max scatters is ignored for some scatterRandomParam, this means
					take all photons from max and above, and put them into
					the top bin
				*/					
				if	((blueScatters > binParams->maxS) &&
						(!ignoreMaxScatters))
					continue;
					
				if ((binParams->numE1Bins > 0) && (bluePhoton.energy < binParams->minE))
						continue;
						
				if	((binParams->numE1Bins > 0) && (bluePhoton.energy > binParams->maxE))
					continue;
					
				if ((binParams->numZBins > 0) && (bluePhoton.location.z_position < binParams->minZ))
					continue;
					
				if ((binParams->numZBins > 0) && (bluePhoton.location.z_position > binParams->maxZ))
					continue;
			}

			/* Increment accepted blue photon count */
			binFields->AccBluePhotons++;

			/* Get scatter count */
			pinkScatters = pinkPhoton.num_of_scatters +
				pinkPhoton.scatters_in_col;
			
			/* See if this photon fails the acceptance criteria */
			{
				if ( (pinkScatters < binParams->minS) &&
						(!ignoreMinScatters) )
					continue;
					
				/* Max scatters is ignored for some scatterRandomParam, this means
					take all photons from max and above, and put them into
					the top bin
				*/					
				if	( (pinkScatters > binParams->maxS) &&
						(!ignoreMaxScatters) )
					continue;
					
				if ((binParams->numE2Bins > 0) && (pinkPhoton.energy < binParams->minE))
						continue;
						
				if	((binParams->numE2Bins > 0) && (pinkPhoton.energy > binParams->maxE))
					continue;
					
				if ((binParams->numZBins > 0) && (pinkPhoton.location.z_position < binParams->minZ))
					continue;
					
				if ((binParams->numZBins > 0) && (pinkPhoton.location.z_position > binParams->maxZ))
					continue;
			}
			
			/* Increment accepted pink photon count */
			binFields->AccPinkPhotons++;

			/* Compute distance and angle if they have specified one or more
				of these bins.  If there is just one bin, then this is done
				simply to check the range.
				
				Otherwise distance and angle are not checked against the range
				specified in the binning parameters file.
			*/
			if ((binParams->numAABins >= 1) || (binParams->numTDBins >= 1)) {
				
				/* Compute distance and angle */
				PhgBinCompPetDA(&(bluePhoton), &(pinkPhoton),
					&angle, &distance);

				/* See if distance is outside acceptance range */
				if ((binParams->numTDBins > 0) && ((distance < binParams->minTD) ||
						(distance > binParams->maxTD))) {
					
					/* Go to next photon */
					continue;
				}
				
				/* Compute the angle index for the image */
				if (binParams->numAABins > 1) {
				
					/* Compute the index and then check for range */
					angleIndex = (LbUsFourByte) 
						floor((
						 (angle+(PHGMATH_PI_DIV2)) * binParams->numAABins)/
						 (PHGMATH_PI));
		
					/* Check for boundary */
					if (angleIndex >= binParams->numAABins)
						angleIndex = binParams->numAABins - 1;
				}
				else {
					angleIndex = 0;
				}
				
				
				/* Compute the distance index for the image */
				if (binParams->numTDBins > 1) {
					distIndex = (LbUsFourByte) 
						floor(((distance - binParams->minTD) * binParams->numTDBins)/
							(binParams->tdRange));
					
					/* Check for boundary */
					if (distIndex >= binParams->numTDBins)
						distIndex = binParams->numTDBins - 1;
				}
				else {
					distIndex = 0;
				}
			}
			
			
			/* Compute the time-of-flight index for the image */
			if (binParams->numTOFBins > 1) {
				
				/* tofDifference is in nanoseconds, hence '1E9*' */
				tofDifference = 1.0E9 * (bluePhoton.travel_distance - pinkPhoton.travel_distance) /
								PHGMATH_SPEED_OF_LIGHT;
								
				if ( (tofDifference >= binParams->maxTOF) ||
					(tofDifference <= binParams->minTOF) ) {
					
					/* outside coincidence window, go to next photon */
					continue;
				
				}
				
				/* make sure that positive TOF is always oriented in the same direction,
				+TOF equating to +x. */
				if (bluePhoton.location.x_position < pinkPhoton.location.x_position) {
					tofDifference = -tofDifference;
				} else if ( (bluePhoton.location.x_position == pinkPhoton.location.x_position)
							&&
							(bluePhoton.location.y_position < pinkPhoton.location.y_position) ) {
					tofDifference = -tofDifference;
				}
								
				tofIndex = (LbUsFourByte) 
					floor(((tofDifference - binParams->minTOF) * binParams->numTOFBins)/
						(binParams->tofRange));
				
				/* Check for boundary */
				if (tofIndex >= binParams->numTOFBins)
					tofIndex = binParams->numTOFBins - 1;
			}
			else {
				tofIndex = 0;
			}

			/* Compute scatter indexes based on the value of scatterRandomParam.
				scatterRandomParam == 0, no binning by scatter/random state -
							all events treated the same
				scatterRandomParam == 1, all non-scattered coincidences in one bin;
							one or more scatters (either photon) in the other
				scatterRandomParam == 2, index for each photon computed from 
							(num_scatters - min) and up to (max - min)
				scatterRandomParam == 3, like 2, but NOTE this is handled by setting
							the scatterIndex to the max value if the computed value
							is too high.  Photons with num scatters >= 
							binParams->maxS were passed through the range filter 
							above for scatterRandomParam 3,5,8, and 10.
				scatterRandomParam == 4, coincidences are binned by using 
							(blueScatters + pinkScatters).  The coincidence is
							rejected if the sum is < binParams->minS or >
							binParams->maxS.
				scatterRandomParam == 5, same as 4, except that when the sum 
							> binParams->maxS the scatterIndex is set to the
							max rather than the coincidence being rejected.
				scatterRandomParam == 6, same as 1 except there is a third bin
							for randoms, i.e., the first bin is for true
							unscattered coincidences, the second for non-random
							scattered coincidences, the third for random
							coincidences.
				scatterRandomParam == 7, same as 2 except there is an extra bin
							for randoms, i.e., non-random coincidences are
							binned as in 2, but all random coincidences are
							binned in an extra bin added at the end of the
							scatter state array.  The scatter state array is
							computed using a single linearized index rather than
							using both scatter1Index and scatter2Index to
							allow for the extra bin for randoms.
				scatterRandomParam == 8, same as 3 except there is an extra bin
							for randoms, i.e., non-random coincidences are
							binned as in 3, but all random coincidences are
							binned in an extra bin added at the end of the
							scatter state array.  The scatter state array is
							computed using a single linearized index rather than
							using both scatter1Index and scatter2Index to
							allow for the extra bin for randoms.
				scatterRandomParam == 9, same as 4 except there is an extra bin
							for randoms, i.e., non-random coincidences are
							binned as in 4, but all random coincidences are
							binned in an extra bin added at the end of the
							scatter state array.
				scatterRandomParam == 10, same as 5 except there is an extra bin
							for randoms, i.e., non-random coincidences are
							binned as in 5, but all random coincidences are
							binned in an extra bin added at the end of the
							scatter state array.
			*/
			{
				 
				/* scatterRandomParam == 0: no unscattered/scattered/random binning */
				if (binParams->scatterRandomParam == 0) {
					scatter1Index = 0;
					scatter2Index = 0;
				}
				
				/* scatterRandomParam == 1: unscattered vs scattered binning */
				else if (binParams->scatterRandomParam == 1) {
				
					scatter2Index = 0;	/* only one scatter index used */
					
					if ((blueScatters == 0) && (pinkScatters == 0)) {
						scatter1Index = 0;	/* unscattered */
					}
					else {
						scatter1Index = 1;	/* scattered */
					}
					
				}
				
				/* scatterRandomParam == 2&3: 2D binning by photons' scatter states */
				else if ( (binParams->scatterRandomParam == 2) || 
						(binParams->scatterRandomParam == 3) ) {
					
					/* Compute scatter1 index based on number of scatters */
					if (binParams->sRange != 0) {
						scatter1Index = (blueScatters - binParams->minS);
								
						/* Check boundary */
						if (scatter1Index >= binParams->numS1Bins)
							scatter1Index = binParams->numS1Bins -1;
							
					}
					else {
						scatter1Index = 0;
					}
					
					/* Compute scatter2 index based on number of scatters */
					if (binParams->sRange != 0) {
						scatter2Index = (pinkScatters - binParams->minS);
								
						/* Check boundary */
						if (scatter2Index >= binParams->numS2Bins)
							scatter2Index = binParams->numS2Bins -1;
							
					}
					else {
						scatter2Index = 0;
					}
					
				}
				
				/* scatterRandomParam == 4&5: 1D binning by sum of photons' scatter states */
				else if ( (binParams->scatterRandomParam == 4) || 
						(binParams->scatterRandomParam == 5) ) {
					
					scatter2Index = 0;	/* only one scatter index used */
					
					if (binParams->sRange != 0) {
						scatter1Index = (blueScatters + pinkScatters - binParams->minS);
								
						/* Check boundaries */
						if (blueScatters + pinkScatters < binParams->minS)
							continue;	/* reject; sum of scatters is less than min */
						
						if (scatter1Index >= binParams->numS1Bins) {
							/* sum of scatters greater than max */
							if (binParams->scatterRandomParam == 4) {
								/* reject */
								continue;
							} else {
								/* for (scatterRandomParam == 5) we accept to max bin */
								scatter1Index = binParams->numS1Bins -1;
							}
						}
						
					}
					else {
						scatter1Index = 0;
					}
					
				}
				
				/* scatterRandomParam == 6: true vs scattered vs random binning */
				else if (binParams->scatterRandomParam == 6) {
				
					scatter2Index = 0;	/* only one scatter index used */
					
					if (decay->decayType == PhgEn_PETRandom) {
						scatter1Index = 2;
					}
					else if ((blueScatters == 0) && (pinkScatters == 0)) {
						scatter1Index = 0;
					}
					else {
						scatter1Index = 1;
					}
					
				}
				
				/* scatterRandomParam == 7&8: 2D binning by photons' scatter states
				 plus extra bin for randoms - 2D binning computed here and done as
				 single index because of extra random bin */
				else if ( (binParams->scatterRandomParam == 7) || 
						(binParams->scatterRandomParam == 8) ) {
					
					if (decay->decayType == PhgEn_PETRandom) {
						scatter1Index = binParams->numS1Bins -1;	/* random */
					}
					
					else if (binParams->sRange != 0) {
						
						/* scatter1Index first used for  blue scatter index value */
						scatter1Index = (blueScatters - binParams->minS);
						/* Check boundary */
						if (scatter1Index >= binParams->sRange )
							scatter1Index = binParams->sRange - 1;
						
						/* scatter2Index will be used temporarily to get the pink
						scatters index value, but cleared after final computation
						of scatter1Index */
						scatter2Index = (pinkScatters - binParams->minS);
						/* Check boundary */
						if (scatter2Index >= binParams->sRange )
							scatter2Index = binParams->sRange - 1;
						
						/* now compute a linearized index for 2D array of blue
						 scatter vs pink scatter */
						scatter1Index = scatter2Index + (scatter1Index * binParams->sRange);
						
								
					}
					else {
						scatter1Index = 0;
					}
					
					scatter2Index = 0;	/* only one scatter index used */

				}
				
				/* scatterRandomParam == 9&10: 1D binning by sum of photons' scatter states
				 plus extra bin for randoms */
				else if ( (binParams->scatterRandomParam == 9) || 
						(binParams->scatterRandomParam == 10) ) {
					
					scatter2Index = 0;	/* only one scatter index used */

					if (decay->decayType == PhgEn_PETRandom) {
						scatter1Index = binParams->numS1Bins -1;	/* random */
					}
					
					if (binParams->sRange != 0) {
						scatter1Index = (blueScatters + pinkScatters - binParams->minS);
								
						/* Check boundaries */
						if (blueScatters + pinkScatters < binParams->minS)
							continue;	/* reject; sum of scatters is less than min */
						
						if (scatter1Index >= binParams->numS1Bins - 1) {
							/* sum of scatters greater than max */
							if (binParams->scatterRandomParam == 9) {
								/* reject */
								continue;
							} else {
								/* for (scatterRandomParam == 10) we accept to max bin */
								scatter1Index = binParams->numS1Bins - 2;
							}
						}
					}
					else {
						scatter1Index = 0;
					}
					
				}
				
			}

			/* Compute energy 1 index */
			if ((binParams->numE1Bins != 0) && (binParams->eRange != 0)) {
				energy1Index = (LbUsFourByte) floor((
						(bluePhoton.energy - binParams->minE) *
						binParams->numE1Bins) / binParams->eRange);
						
				/* Check boundary */
				if (energy1Index >= binParams->numE1Bins)
					energy1Index = binParams->numE1Bins -1;
			}
			else {
				energy1Index = 0;
			}
			
			/* Compute energy 2 index */
			if ((binParams->numE2Bins != 0) && (binParams->eRange != 0)) {
				energy2Index = (LbUsFourByte) floor((
						(pinkPhoton.energy - binParams->minE) *
						 binParams->numE2Bins) / binParams->eRange);
						
				/* Check boundary */
				if (energy2Index >= binParams->numE2Bins)
					energy2Index = binParams->numE2Bins -1;
			}
			else {
				energy2Index = 0;
			}
			
			/* Compute crystal indexes */
			if (binParams->numCrystalBins != 0) {
				
				/* if debug is set, check to make sure that the crystal numbers are in range */
				#ifdef PHG_DEBUG
					if ((pinkPhoton.detCrystal < 0) || (bluePhoton.detCrystal < 0)) {
						PhgAbort("Invalid computation of crystal index (1) (PhgBinPETPhotons)", true);
					}
				#endif
						
				#ifdef PHG_DEBUG
					if ((pinkPhoton.detCrystal >= (LbFourByte)(binParams->numCrystalBins)) || 
							(bluePhoton.detCrystal >= (LbFourByte)(binParams->numCrystalBins))) {
						PhgAbort("Invalid computation of crystal index (2) (PhgBinPETPhotons)", true);
					}
				#endif
						
				/* Assign indices */
				if (pinkPhoton.detCrystal > bluePhoton.detCrystal) {
					crystal1Index = bluePhoton.detCrystal;
					crystal2Index = pinkPhoton.detCrystal;
				} else {
					crystal2Index = bluePhoton.detCrystal;
					crystal1Index = pinkPhoton.detCrystal;
				}

			}
			
			/* Do Z axis binning */
			if ((binParams->numZBins > 1) && (binParams->doSSRB == false) && (binParams->doMSRB == false)) {
				/* Pick "up/down" Z indexes */
				if (bluePhoton.location.y_position < pinkPhoton.location.y_position) {
					z1IndexPtr = &zDownIndex;
					z2IndexPtr = &zUpIndex;
				}
				else if (bluePhoton.location.y_position > pinkPhoton.location.y_position) {
					z1IndexPtr = &zUpIndex;
					z2IndexPtr = &zDownIndex;
				}
				else if (bluePhoton.location.x_position < pinkPhoton.location.x_position) {
					z1IndexPtr = &zDownIndex;
					z2IndexPtr = &zUpIndex;
				}
				else {
					z1IndexPtr = &zUpIndex;
					z2IndexPtr = &zDownIndex;
				}
			
				/* Compute z 1 index */
				if (binParams->zRange != 0) {
					*z1IndexPtr = (LbUsFourByte) floor((								
							(bluePhoton.location.z_position - binParams->minZ) *	
							 binParams->numZBins) / binParams->zRange);
							
					/* Check boundary */
					if (*z1IndexPtr >= binParams->numZBins)
						*z1IndexPtr = binParams->numZBins -1;
				}
				else 
					*z1IndexPtr = 0;
				
				/* Compute z 2 index */
				if (binParams->zRange != 0) {
					*z2IndexPtr = (LbUsFourByte) floor((
							(pinkPhoton.location.z_position - binParams->minZ) *
							 binParams->numZBins) / binParams->zRange);
							
					/* Check boundary */
					if (*z2IndexPtr >= binParams->numZBins)
						*z2IndexPtr = binParams->numZBins -1;
				}
				else 
					*z2IndexPtr = 0;
					
				/* Set z low/high for loop */
				lowZIndex = 0;
				highZIndex = 0;
				
			}
			else if (binParams->doSSRB == true) {
			
				/* Note we use the zUpIndex because it exists */
				avgZ = (bluePhoton.location.z_position + pinkPhoton.location.z_position)/2;
				
				/* Check boundaries */
				if ((avgZ < binParams->minZ) || (avgZ > binParams->maxZ))
					continue;
				
				/* Compute index */
				zUpIndex = (LbUsFourByte) floor((
						(avgZ - binParams->minZ) *
						 binParams->numZBins) / binParams->zRange);
				if (zUpIndex == binParams->numZBins)
					zUpIndex--;
					
				#ifdef PHG_DEBUG
					if (zUpIndex >= binParams->numZBins) {
						PhgAbort("Invalid computation of zUpIndex -SSRB (PhgBinPETPhotons)", true);
					}
				#endif
				
				/* Clear zDownIndex since it's not used in SSRB */	 
				zDownIndex = 0;
					
				/* Set z low/high for loop */
				lowZIndex = 0;
				highZIndex = 0;
			}
			else if (binParams->doMSRB == true) {
			
				/* Note we use the zUpIndex because it exists */
				avgZ = (bluePhoton.location.z_position + pinkPhoton.location.z_position)/2;

				
				lowZ = avgZ - ((fabs(pinkPhoton.location.z_position - bluePhoton.location.z_position)
					* phgBinObjDiameter)/(2 * phgBinDetDiameter));
					

				highZ = avgZ + ((fabs(pinkPhoton.location.z_position - bluePhoton.location.z_position)
					* phgBinObjDiameter)/(2 * phgBinDetDiameter));
				
				/* Compute low Z index */
				lowZIndex = (LbUsFourByte) floor((
						(lowZ - binParams->minZ) *
						 binParams->numZBins) / binParams->zRange);
				if (lowZIndex == binParams->numZBins)
					lowZIndex--;
					
				#ifdef PHG_DEBUG
					if (lowZIndex >= binParams->numZBins) {
						PhgAbort("Invalid computation of lowZIndex -MSRB (PhgBinPETPhotons)", true);
					}
				#endif
				
				/* Compute high Z index */
				highZIndex = (LbUsFourByte) floor((
						(highZ - binParams->minZ) *
						 binParams->numZBins) / binParams->zRange);
				if (highZIndex == binParams->numZBins)
					highZIndex--;
					
				#ifdef PHG_DEBUG
					if (highZIndex >= binParams->numZBins) {
						PhgAbort("Invalid computation of highZIndex -MSRB (PhgBinPETPhotons)", true);
					}
				#endif
			}
			else {
				zDownIndex = 0;
				zUpIndex = 0;
			}
			
			/* Do 3DRP parameters */
			{
				if (binParams->numThetaBins || binParams->numPHIBins || binParams->numXRBins || binParams->numYRBins) {

					phi = atan2((pinkPhoton.location.y_position - bluePhoton.location.y_position),
						(pinkPhoton.location.x_position-bluePhoton.location.x_position));
					
					/* Convert to range [-pi/2n, pi - pi/2n] */
					if (binParams->numPHIBins > 0) {
						d = (PHGMATH_PI/(2*binParams->numPHIBins));
					}
					else {
						d = 0;
					}
					
					flip = 1.0;
					
					if (phi < -d) {
						phi += PHGMATH_PI;
						flip = -1.0;
					}
					else if ((PHGMATH_PI - d) < phi) {
						phi -= PHGMATH_PI;
						flip = -1.0;
					}

					theta = flip * atan((pinkPhoton.location.z_position-bluePhoton.location.z_position)/
						PHGMATH_SquareRoot(
							PHGMATH_Square(pinkPhoton.location.x_position-bluePhoton.location.x_position) +
							PHGMATH_Square(pinkPhoton.location.y_position-bluePhoton.location.y_position)));
				}
				else {
					phiIndex = 0;
					thetaIndex = 0;
				}
				
				if (binParams->numPHIBins > 1) {

					/* Check range */
					if ((phi < binParams->minPhi) || (phi > binParams->maxPhi))
						continue;
					
					phiIndex = (LbUsFourByte) floor((
							(phi - binParams->minPhi) *
							 binParams->numPHIBins) / binParams->phiRange);
					if (phiIndex == binParams->numPHIBins)
						phiIndex--;
					
					#ifdef PHG_DEBUG
						if (phiIndex >= binParams->numPHIBins) {
							PhgAbort("Invalid computation of phiIndex (PhgBinPETPhotons)", true);
						}
					#endif
							 
				}
				
				if (binParams->numThetaBins > 1) {

					/* Check range */
					if ((theta < binParams->minTheta) || (theta > binParams->maxTheta))
						continue;
					
					thetaIndex = (LbUsFourByte) floor((
							(theta - binParams->minTheta) *
							 binParams->numThetaBins) / binParams->thetaRange);
					if (thetaIndex == binParams->numThetaBins)
						thetaIndex--;
					
					#ifdef PHG_DEBUG
						if (thetaIndex >= binParams->numThetaBins) {
							PhgAbort("Invalid computation of thetaIndex (PhgBinPETPhotons)", true);
						}
					#endif
							 
				}
					
				if (binParams->numXRBins > 1) {			 
				
					xr = (-(bluePhoton.location.x_position * PHGMATH_Sine(phi))) + 
						(bluePhoton.location.y_position * PHGMATH_Cosine(phi));
							 
					/* Check range */
					if ((xr < binParams->minXR) || (xr > binParams->maxXR))
						continue;
					
					xrIndex = (LbUsFourByte) floor((
							(xr - binParams->minXR) *
							 binParams->numXRBins) / binParams->xrRange);
					if (xrIndex == binParams->numXRBins)
						xrIndex--;
					
					#ifdef PHG_DEBUG
						if (xrIndex >= binParams->numXRBins) {
							PhgAbort("Invalid computation of xrIndex (PhgBinPETPhotons)", true);
						}
					#endif
						
						
				}
				else
					xrIndex = 0;
				
				if (binParams->numYRBins > 1 ) {
					yr = (-(bluePhoton.location.x_position * PHGMATH_Cosine(phi) * PHGMATH_Sine(theta))) - 
						(bluePhoton.location.y_position * PHGMATH_Sine(phi) * PHGMATH_Sine(theta)) +
						(bluePhoton.location.z_position * PHGMATH_Cosine(theta));
							 
					/* Check range */
					if ((yr < binParams->minYR) || (yr > binParams->maxYR))
						continue;
					
					yrIndex = (LbUsFourByte) floor((
							(yr - binParams->minYR) *
							 binParams->numYRBins) / binParams->yrRange);

					if (yrIndex == binParams->numYRBins)
						yrIndex--;
					
					#ifdef PHG_DEBUG
						if (yrIndex >= binParams->numYRBins) {
							PhgAbort("Invalid computation of yrIndex (PhgBinPETPhotons)", true);
						}
					#endif
				}
				else 
					yrIndex = 0;
			}
			
			/* Set the weight variables */
			blueWeight = bluePhoton.photon_current_weight;
			pinkWeight = pinkPhoton.photon_current_weight;

			/* Incremement accepted coincidence count--this is later
			 * decremented if the user rejects the coincidence */
			binFields->NumAcceptedCoincidences++;
			
			/* Although most binning is probably not MSRB, in order to generalize the code
				we have the following loop construct which is only executed once for non-MSRB
				binning
			*/
			
			for (zIndex = lowZIndex; zIndex <= highZIndex; zIndex++) {

				/* Compute the coincidence weight */
				coincidenceWeight = (decay->startWeight * blueWeight *
					pinkWeight) * binFields->WeightRatio;
			
				/* Adjust weight and set indexes if doing MSRB */
				if (binParams->doMSRB == true) {

					/* Adjust weight */
					coincidenceWeight /= ((highZIndex - lowZIndex)+1);
					
					/* Set Indexes */
					zUpIndex = zIndex;
					zDownIndex = 0;
				}
				
				/* Compute the coincidence square weight */
				coincidenceSquWeight = PHGMATH_Square(coincidenceWeight);

				/* Convert index to one dimension */
				imageIndex = (energy2Index * binParams->energy2CIsize) +
					(energy1Index * binParams->energy1CIsize) +
					(distIndex * binParams->tdCIsize) + 
					(tofIndex * binParams->tofCIsize) + 
					(angleIndex * binParams->aaCIsize) +
					(zUpIndex * binParams->z2CIsize) +
					(zDownIndex * binParams->z1CIsize) +
					(scatter2Index * binParams->scatter2CIsize) +
					(scatter1Index * binParams->scatter1CIsize) +
					(crystal2Index * binParams->crystal2CIsize) +
					(crystal1Index * binParams->crystal1CIsize) +
					(phiIndex * binParams->phiCIsize) +
					(thetaIndex * binParams->thetaCIsize) +
					(xrIndex * binParams->xrCIsize) +
					(yrIndex * binParams->yrCIsize);
		
				/* Call the user binning routine and continue with next photon if rejected */
				if (BinUsrModPETPhotonsF2Ptr && 
						(*BinUsrModPETPhotonsF2Ptr)(
							binParams,
							binData,
							decay,
							&bluePhoton,
							&pinkPhoton,
							&angleIndex,
							&distIndex,
							&tofIndex,
							&scatter1Index,
							&scatter2Index,
							&crystal1Index,
							&crystal2Index,
							&energy1Index,
							&energy2Index,
							&zDownIndex,
							&zUpIndex,
							&thetaIndex,
							&phiIndex,
							&xrIndex,
							&yrIndex,
							&imageIndex,
							&coincidenceWeight,
							&coincidenceSquWeight) == false) {
					
					/* the number of accepted coincidences, which was incremented above,
					 * is decremented to delete this coincidence.  However, this is
					 * skipped for MSRB, as the coincidence may be included in another slice */
					if (binParams->doMSRB == false) {
						binFields->NumAcceptedCoincidences--;
					}
					
					continue;
				}

				/* Write detected photons to history file, if requested */
				if (binParams->isHistoryFile) {

					/* Write the photons */
					if (PhoHFileWriteDetections(&binFields->historyFileHk, decay,
							&bluePhoton,
							1,
							&pinkPhoton,
							1) == false) {
						
						/* Abort Program execution */
						PhgAbort("Got failure from PhoHFileWriteDetections (BinPETPhotons).", true);
					}
				}
				
				#ifdef PHG_DEBUG
					if (imageIndex >= (binParams->numImageBins)) {
						sprintf(phgBinErrStr,"Invalid imageIndex computed\n"
							"\timageIndex = %ld\n"
							"\tbinParams->numImageBins = %ld"
							"\tcrystalIndex1 = %ld\n"
							"\tcrystalIndex2 = %ld\n"
							"\tenergyIndex1 = %ld\n"
							"\tenergyIndex2 = %ld\n"
							"\tdistIndex = %ld\n"
							"\ttofIndex = %ld\n"
							"\tangleIndex = %ld\n"
							"\tzIndex1 = %ld\n"
							"\tzIndex2 = %ld\n"
							"\tscatterIndex1 = %ld"
							"\tscatterIndex2 = %ld"
							"\tphiIndex = %ld"
							"\txrIndex = %ld"
							"\tyrIndex = %ld"
							"\n in (PhgBinPETPhotons)\n",
							 (unsigned long)imageIndex, (unsigned long)binParams->numImageBins, 
							 (unsigned long)crystal1Index, (unsigned long)crystal2Index, 
							 (unsigned long)energy1Index, (unsigned long)energy2Index,
							 (unsigned long)distIndex, (unsigned long)tofIndex, 
							 (unsigned long)angleIndex, (unsigned long)zUpIndex,
							 (unsigned long)zDownIndex, (unsigned long)scatter1Index, 
							 (unsigned long)scatter2Index, (unsigned long)phiIndex, 
							 (unsigned long)xrIndex, (unsigned long)yrIndex);
						
						PhgAbort(phgBinErrStr,true);

					}
				#endif

				/* Compute sum statistics */
				binFields->AccCoincidenceWeight += coincidenceWeight;
				binFields->AccCoincidenceSquWeight += coincidenceSquWeight;
				
				
				/* Increment the image */
				phgBinIncrementPETImage(binParams, binData, imageIndex, decay, &bluePhoton, &pinkPhoton, coincidenceWeight,
					coincidenceSquWeight);
			}			
		}
	}
}

/*********************************************************************************
*
*			Name:			phgBinIncrementPETImage
*
*			Summary:		Increment the PET image.
*
*			Arguments:
*				PHG_BinParamsTy		*binParams		- User defined binning parameters.
*				PHG_BiDataTy		*binData	 	- Storage for binned data.
*			Function return: None.
*
*********************************************************************************/
void phgBinIncrementPETImage(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData,
	LbUsEightByte imageIndex, PHG_Decay *decay,
	PHG_TrackingPhoton *bluePhoton, PHG_TrackingPhoton *pinkPhoton, double coincidenceWeight,
	double coincidenceSquWeight)

{
	PHG_Decay*					dummyPtr1;		/* Removes compiler warning */
	PHG_TrackingPhoton*			dummyPtr2;		/* Removes compiler warning */
	
	
	/* Avoid unused parameter compiler warnings */
	dummyPtr1 = decay;
	dummyPtr2 = bluePhoton;
	dummyPtr2 = pinkPhoton;
	
	/* Increment count image if they are creating one */
	if (binParams->doCounts != false) {
	
		/* Increment count image count */
		switch(binParams->count_image_type){
			case PHG_BIN_COUNT_TYPE_I1:
			
				/* Check for overflow */
				if (((LbUsOneByte *)binData->countImage)[imageIndex]
						== LBUSONEBYTE_MAX) {
					PhgAbort("Overflow in count histogram (PhgBinPETPhotons)",
						true);
				}
				
				/* Update the image */
				((LbUsOneByte *)binData->countImage)[imageIndex] += 1;
				break;
				
			case PHG_BIN_COUNT_TYPE_I2:
			
				/* Check for overflow */
				if (((LbUsTwoByte *)binData->countImage)[imageIndex]
						== LBUSTWOBYTE_MAX) {
					PhgAbort("Overflow in count histogram (PhgBinPETPhotons)",
						true);
				}
				
				/* Update the image */
				((LbUsTwoByte *)binData->countImage)[imageIndex] += 1;
				break;
				
			case PHG_BIN_COUNT_TYPE_I4:
			
				/* Check for overflow */
				if (((LbUsFourByte *)binData->countImage)[imageIndex]
						== LBUSFOURBYTE_MAX) {
					PhgAbort("Overflow in count histogram (PhgBinPETPhotons)",
						true);
				}
				
				/* Update the image */
				((LbUsFourByte *)binData->countImage)[imageIndex] += 1;
				break;
		}
	}
	/* Update the weight image if they are computing one */
	if (binParams->doWeights == true) {

		if ((binParams->sumAccordingToType == true) && (binParams->weight_image_type == PHG_BIN_WEIGHT_TYPE_R4)) {

			/* Update the image */
			((float *)binData->weightImage)[imageIndex] += coincidenceWeight;

		}
		else {
			
			/* Update the image */
			((double *)binData->weightImage)[imageIndex] += coincidenceWeight;
		}
	}
	
	/* Update the weight squared image if they are computing one */
	if (binParams->doWeightsSquared == true) {

		/* Increment weight square images  */
		if ((binParams->sumAccordingToType == true) && (binParams->weight_image_type == PHG_BIN_WEIGHT_TYPE_R4)) {
				
				/* We don't check for floating point overflow because
					the hardware will get it
				*/
				((float *)binData->weightSquImage)[imageIndex] += 
					coincidenceSquWeight;
		}
		else {
			
				/* We don't check for floating point overflow because
					the hardware will get it
				*/
				((double *)binData->weightSquImage)[imageIndex] += 
					coincidenceSquWeight;
		}
	}
}

/*********************************************************************************
*
*			Name:			PhgBinSPECTPhotons
*
*			Summary:		Update the binning images with the current batch of
*							detected photons.
*
*			Arguments:
*				PHG_BinParamsTy		*binParams		- User defined binning parameters.
*				PHG_BiDataTy		*binData		- Storage for binned data.
*				PHG_BinFieldsTy 	*binFields		- Various binning information.
*				PHG_Decay			decay			- The decay that started the process.
*				PHG_TrackingPhoton *photons			- The photons detected.
*				LbUsFourByte 		numPhotons		- The number of blue photons.
*			Function return: None.
*
*********************************************************************************/
void PhgBinSPECTPhotons(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData, PHG_BinFieldsTy *binFields,
		PHG_Decay *decay,
		PHG_TrackingPhoton *photons, LbUsFourByte numPhotons)
{
	double			distance;			/* */
	double			angle;				/* */
	double			detectionWeight;	/* */
	double			detectionSquWeight;	/* */
	LbUsFourByte	index;				/* Current blue photon */
	LbUsFourByte	scatters;			/* Number of scatters for current blue photon */
	LbUsFourByte	angleIndex = 0;		/* Index for angle bin */
	LbUsFourByte	distIndex = 0;		/* Index for dist bin */
	LbUsFourByte	scatterIndex = 0;	/* Index for scatter bin */
	LbUsFourByte	energyIndex = 0;	/* Index for energy bin */
	LbUsFourByte	crystalIndex = 0;	/* Index for crystal bin */
	LbUsFourByte	zIndex = 0;			/* Index for z bin */
	LbUsFourByte	imageIndex;			/* Index for image */
	
		
	/* Compute statistics */
	binFields->TotBluePhotons += numPhotons;

	/* Process coincidences */
	for (index = 0; index < numPhotons; index++) {
		
		/* Call the user binning routine and continue with next photon if rejected */
		if (BinUsrModSPECTPhotonsFPtr && 
				(*BinUsrModSPECTPhotonsFPtr)(binParams, binData, decay,
					&(photons[index])) == false) {
			
			continue;
		}
		
		/* Get scatter count */
		scatters = photons[index].num_of_scatters +
			photons[index].scatters_in_col;
		
		/* See if this photon fails criteria */
		{
			if (scatters < binParams->minS)
				continue;
	
			/* Max scatters is ignored if scatterRandomParam == 3, this means
				take all photons from max and above, and put them into
				the top bin
			*/
			if ((scatters > binParams->maxS) &&
					(binParams->scatterRandomParam != 3))
				continue;
				
			if ((binParams->numE1Bins > 0) && (photons[index].energy < binParams->minE))
				continue;
	
			if ((binParams->numE1Bins > 0) && (photons[index].energy > binParams->maxE))
				continue;
	
			if ((binParams->numZBins > 0) && (photons[index].location.z_position <
					binParams->minZ))
				continue;
	
			if ((binParams->numZBins > 0) && (photons[index].location.z_position >
					binParams->maxZ))
				continue;
		}

		/* Compute distance angle */  			
		PhgBinCompSpectDA(&(photons[index]),
			&angle, &distance);
		
		/* See if outside distance range */
		if ((binParams->numTDBins > 0) && ((distance < binParams->minTD) ||
				(distance > binParams->maxTD))) {
			
			/* Go to next photon */
			continue;
		}
		
		/* Compute the angle index for the image 
			NOTE That if collimation is being used in UNC SPECT,
			we have already computed this value.
		*/
		if (PHG_IsCollimateOnTheFly() && (ColIsUNC() == true)) {
			if (binParams->numAABins > 1){
				angleIndex = (LbUsFourByte) angle;
			}
			else 
				angleIndex = 0;
		}
		else if (binParams->numAABins >= 1){
			angleIndex = (LbUsFourByte) 
				floor((angle * binParams->numAABins)/(PHGMATH_2PI));
	
			#ifdef PHG_DEBUG
				if (angleIndex > binParams->numAABins) {
					sprintf(phgBinErrStr,"Got an angle index > number of bins.\n"
						"angleIndex = %ld, angle = %3.4f", (unsigned long)angleIndex, angle);
					PhgAbort(phgBinErrStr,true);
				}
			#endif
					
			/* Check for boundary */
			if (angleIndex == binParams->numAABins)
				angleIndex = binParams->numAABins - 1;
		}
		else {
			angleIndex = 0;
		}
		
		/* Compute the distance index for the image */
		if (binParams->numTDBins > 1){
			distIndex = (LbUsFourByte) 
				floor(((distance - binParams->minTD) * binParams->numTDBins)/
				(binParams->tdRange));
			
			/* Check for boundary */
			if (distIndex == binParams->numTDBins)
				distIndex = binParams->numTDBins - 1;
		}
		else
			distIndex = 0;
			
		/* Compute scatter indexes based on the value of scatterRandomParam.
			scatterRandomParam == 0, no scatter bin
			scatterRandomParam == 1, all primary pairs in one bin; one or more scatters in the other
			scatterRandomParam == 2, index computed from min and max
			scatterRandomParam == 3, NOTE this is handled by setting binParams->maxS and
			letting scatters through the range filter at the top of this routine. What
			happens is the bin index is computed greater than num bins and the index
			gets set to max bins.  Remember, scatterRandomParam == 3 means all photons
			whose num scatters >= minS and num scatters <= maxS go in the computed
			bin, AND all photons whose num scatters is greater than maxS also go
			in the top bin; they are not filtered out above.
		 */
		if (binParams->scatterRandomParam == 0) {
			scatterIndex = 0;
		}
		else if (binParams->scatterRandomParam == 1) {
			if (scatters == 0) {
				scatterIndex = 0;
			}
			else {
				scatterIndex = 1;
			}
		}
		else {
			if (binParams->sRange != 0) {
				scatterIndex = (LbUsFourByte)(scatters - binParams->minS);
						
				/* Check Boundary */
				if (scatterIndex >= binParams->numS1Bins)
					scatterIndex = binParams->numS1Bins -1;
			}
			else {
				scatterIndex = 0;
			}
		}	
		
		/* Compute energy 1 index */
		if ((binParams->eRange != 0) && (binParams->numE1Bins > 1)) {
			energyIndex = (LbUsFourByte) (((photons[index].energy - binParams->minE) *
					 binParams->numE1Bins) / binParams->eRange);
						
				/* Check Boundary */
			if (energyIndex >= binParams->numE1Bins)
				energyIndex = binParams->numE1Bins -1;
		}
		else {
			energyIndex = 0;
		}
					
		/* Compute crystal index */
		if (binParams->numCrystalBins != 0) {
			
			/* if debug is set, check to make sure that the crystal number is in range */
			#ifdef PHG_DEBUG
				if (photons[index].detCrystal < 0) {
					PhgAbort("Invalid computation of crystal index (1) (PhgBinPETPhotons)", true);
				}
			#endif
					
			#ifdef PHG_DEBUG
				if (photons[index].detCrystal >= (LbFourByte)(binParams->numCrystalBins)) {
					PhgAbort("Invalid computation of crystal index (2) (PhgBinPETPhotons)", true);
				}
			#endif
					
			/* Assign index */
			crystalIndex = photons[index].detCrystal;

		}
			
		/* Compute z index */
		if ((binParams->zRange != 0) && (binParams->numZBins > 1)) {
			zIndex = (LbUsFourByte) floor(((photons[index].location.z_position - binParams->minZ) *
					 binParams->numZBins) / binParams->zRange);
						
				/* Check Boundary */
			if (zIndex >= binParams->numZBins)
				zIndex = binParams->numZBins -1;
		}
		else {
			zIndex = 0;
		}
		
		/* Convert index to one dimension */
		imageIndex = (energyIndex * binParams->energy1CIsize) +
			(crystalIndex * binParams->crystal1CIsize) + 
			(distIndex * binParams->tdCIsize) + 
			(angleIndex * binParams->aaCIsize) +
			(zIndex * binParams->z2CIsize) +
			(scatterIndex * binParams->scatter1CIsize);
		
		/* Call the user binning routine and continue with next photon if rejected */
		if (BinUsrModSPECTPhotonsF2Ptr && 
				(*BinUsrModSPECTPhotonsF2Ptr)(binParams, binData, decay,
					&(photons[index]),
					&angleIndex,
					&distIndex,
					&scatterIndex,
					&energyIndex,
					&crystalIndex,
					&zIndex,
					&imageIndex) == false) {
			
			continue;
		}

		/* Write detected photons to history file, if requested */
		if (binParams->isHistoryFile) {

			/* Write the photons */
			if (PhoHFileWriteDetections(&binFields->historyFileHk, decay,
					&photons[index],
					1,
					0,
					0) == false) {
				
				/* Abort Program execution */
				PhgAbort("Got failure from PhoHFileWriteDetections (BinSPECTPhotons).", true);
			}
		}
		
		#ifdef PHG_DEBUG
			if (imageIndex >= binParams->numImageBins) {
				sprintf(phgBinErrStr,"Invalid imageIndex computed\n"
					"\timageIndex = %ld\n"
					"\tbinParams->numImageBins = %ld\n"
					"\tenergyIndex = %ld\n"
					"\tcrystalIndex = %ld\n"
					"\tdistIndex = %ld\n"
					"\tangleIndex = %ld\n"
					"\tzIndex = %ld\n"
					"\tscatterIndex = %ld"
					"\n in (PhgBinSPECTPhotons)\n",
					(unsigned long)imageIndex, (unsigned long)binParams->numImageBins, 
					(unsigned long)energyIndex, (unsigned long)crystalIndex, 
					(unsigned long)distIndex, (unsigned long)angleIndex,
					(unsigned long)zIndex, (unsigned long)scatterIndex);
				PhgAbort(phgBinErrStr,true);
			}
		#endif
		
		/* Increment appropriate images */
		if (binParams->doCounts != false) {

			/* Increment image count */
				/* Increment count image count */
				switch(binParams->count_image_type){
					case PHG_BIN_COUNT_TYPE_I1:
					
						/* Check for overflow */
						if (((LbUsOneByte *)binData->countImage)[imageIndex]
								== LBUSONEBYTE_MAX) {
							PhgAbort("Overflow in count histogram (PhgBinSPECTPhotons)",
								true);
						}
							/* Update the image */
						((LbUsOneByte *)binData->countImage)[imageIndex] += 1;	
						break;
						
					case PHG_BIN_COUNT_TYPE_I2:
					
						/* Check for overflow */
						if (((LbUsTwoByte *)binData->countImage)[imageIndex]
								== LBUSTWOBYTE_MAX) {
							PhgAbort("Overflow in count histogram (PhgBinSPECTPhotons)",
								true);
						}
						
						/* Update the image */
						((LbUsTwoByte *)binData->countImage)[imageIndex] += 1;
						break;
						
					case PHG_BIN_COUNT_TYPE_I4:
					
						/* Check for overflow */
						if (((LbUsFourByte *)binData->countImage)[imageIndex]
								==  LBUSFOURBYTE_MAX) {
							PhgAbort("Overflow in count histogram (PhgBinSPECTPhotons)",
								true);
						}
						
						/* Update the image */
						((LbUsFourByte *)binData->countImage)[imageIndex] += 1;
						break;
				}
		}
		
		/* Compute the weight variables */
		detectionWeight = (decay->startWeight * 
			photons[index].photon_current_weight) * binFields->WeightRatio;
		
		/* Compute the adjusted detection weight squared */			
		detectionSquWeight = PHGMATH_Square(detectionWeight);
		 
		/* Compute sum statistics */
		binFields->AccBluePhotons++;
		binFields->AccCoincidenceWeight += detectionWeight;
		binFields->AccCoincidenceSquWeight += detectionSquWeight;
	
		/* Update the weight image if they are computing one */
		if (binParams->doWeights == true) {
			/* Increment weightimages  */
			if ((binParams->sumAccordingToType == true) && (binParams->weight_image_type == PHG_BIN_WEIGHT_TYPE_R4)) {
				
					/* Update the image */
					((float *)binData->weightImage)[imageIndex] += 
						detectionWeight;
			}
			else {
				
					/* Update the image */
					((double *)binData->weightImage)[imageIndex] += 
						detectionWeight;
			}
		}
		
		/* Update the weight squared image if they are computing one */
		if (binParams->doWeightsSquared == true) {

			/* Increment weightimages  */
			if ((binParams->sumAccordingToType == true) && (binParams->weight_image_type == PHG_BIN_WEIGHT_TYPE_R4)) {
						
					/* Update the image */
					((float *)binData->weightSquImage)[imageIndex] += 
						detectionSquWeight;
			}
			else {
						
					/* Update the image */
					((double *)binData->weightSquImage)[imageIndex] += 
						detectionSquWeight;
			}
		}
	}
}

/*********************************************************************************
*
*			Name:			PhgBinCompSpectDA
*
*			Summary:		Compute distance angle indexes for SPECT event.
*
*			Arguments:
*				PHG_TrackingPhoton *photonPtr		- The detected photon.
*				double 				*anglePtr		- The computed angle.
*				double			 	*distancePtr	- The computed distance.
*			Function return: None.
*
*********************************************************************************/
void PhgBinCompSpectDA(PHG_TrackingPhoton *photonPtr, double *anglePtr,
		double *distancePtr)
{
	/* Compute the transaxial position based on the type of detector/collimator */
	if (DetIsCylindrical() || DetIsBlock()) {
	
		/* The angle is computed from the detected position of the photon
		*/
		*anglePtr = atan2(photonPtr->location.y_position, 
					photonPtr->location.x_position);
		
		if (*anglePtr < 0) {
			
			*anglePtr = *anglePtr + (2 * PHGMATH_PI);
			
		}
		
		/* without a coincident photon or collimation there is no way to determine
		 a distance: we just set it to 0 */
		*distancePtr = 0.0;
		
	}
	else if (PHG_IsCollimateOnTheFly() && (ColIsUNC() == true)) {
	
		/* The transaxial position and angle index were computed in the collimator
			module so we get them here for free
		*/
		*distancePtr = photonPtr->transaxialPosition;
		*anglePtr = photonPtr->azimuthalAngleIndex;
		
	}
	else if (PHG_IsDetectOnTheFly() && (DetIsUNCSpect() || DetIsPlanar())) {
	
		/* The transaxial position and angle index were computed in the detector
			module so we get them here for free
		*/
		*distancePtr = photonPtr->detLocation.y_position;
		*anglePtr = photonPtr->detectorAngle;

	}
	else {
		/* Compute the angle when there is no detector or only the simple
		 spect detector */
		*anglePtr = atan2(photonPtr->angle.cosine_y, photonPtr->angle.cosine_x)
					+PHGMATH_PI;
	
		
		/* Assume the collimator is perpendicular to the direction of travel */
		*distancePtr = ((photonPtr->location.y_position * photonPtr->angle.cosine_x) -
			(photonPtr->location.x_position * photonPtr->angle.cosine_y)) /
			PHGMATH_SquareRoot(1 - PHGMATH_Square(photonPtr->angle.cosine_z));	
	}	
}

void PhgBinCompPetDAOld(PHG_TrackingPhoton *bluePhotonPtr, PHG_TrackingPhoton *pinkPhotonPtr, 
		double *anglePtr, double *distancePtr);

/*********************************************************************************
*
*			Name:			PhgBinCompPetDA
*
*			Summary:		Compute distance angle indexes for PET event.
*
*			Arguments:
*				PHG_TrackingPhoton *bluePhotonPtr		- The detected blue photon.
*				PHG_TrackingPhoton *pinkPhotonPtr		- The detected pink photon.
*				double 				*anglePtr			- The computed angle.
*				double			 	*distancePtr		- The computed distance.
*			Function return: None.
*
*********************************************************************************/
void PhgBinCompPetDA(PHG_TrackingPhoton *bluePhotonPtr, PHG_TrackingPhoton *pinkPhotonPtr, 
		double *anglePtr, double *distancePtr)
{
	double			temp;			/* Temp var for computation */
	PHG_Position	*bLoc, *pLoc;	/* Temps for computation */
	
	/* Set temps for conveience */
	bLoc = &bluePhotonPtr->location;
	pLoc = &pinkPhotonPtr->location;
		
	/* Compute theta1 for the blue photon */
	*anglePtr = atan2(pLoc->y_position - bLoc->y_position,
			  pLoc->x_position - bLoc->x_position);
	
	/* Check for range */
	if (*anglePtr < -(PHGMATH_PI_DIV2)) {
		*anglePtr = *anglePtr + PHGMATH_PI;
	}
	else if (*anglePtr >= (PHGMATH_PI_DIV2)) {
		*anglePtr = *anglePtr - PHGMATH_PI;
	}
	
	
	if (PhgMathRealNumAreEqual(bLoc->x_position,
			pLoc->x_position, -7, 0, 0, 0) == false) {

		/* Compute distance */
		*distancePtr = fabs((bLoc->x_position * PHGMATH_Sine(*anglePtr)) -
			(bLoc->y_position * PHGMATH_Cosine(*anglePtr)));
	
		temp = ((pLoc->y_position - bLoc->y_position)/
			(pLoc->x_position - bLoc->x_position));
		
		if ((bLoc->y_position - (temp * bLoc->x_position)) < 0) {
			*distancePtr = *distancePtr * -1.0;
		}
	}
	else {
		*distancePtr = bLoc->x_position;
	}
	
	/* If debugging do some tests */
	#ifdef PHG_DEBUG_NEW_DA_COMP
	{
		double angle, distance;
		
		/* Compute old way */
		PhgBinCompPetDAOld(bluePhotonPtr, pinkPhotonPtr, &angle, &distance);
		
		/* Compare angles */
		if (PhgMathRealNumAreEqual(*anglePtr,
				angle, -3, 0, 0, 0) == false) {
		
			LbInPrintf("\nNon-matching angles: angle.new = %3.4f, angle.old = %3.4f\n"
				"\t  x1 = %3.4f,  x2 = %3.4f, y1 = %3.4f, y2 = %3.4f",
				*anglePtr, angle, bLoc->x_position, pLoc->x_position,
				bLoc->y_position, pLoc->y_position);
		}
		
		/* Compare distances */
		if (PhgMathRealNumAreEqual(*distancePtr,
				distance, -3, 0, 0, 0) == false) {
		
			LbInPrintf("\nNon-matching distances: distance.new = %3.4f, distance.old = %3.4f\n"
				"\t  x1 = %3.4f,  x2 = %3.4f, y1 = %3.4f, y2 = %3.4f",
				*distancePtr, distance, bLoc->x_position, pLoc->x_position,
				bLoc->y_position, pLoc->y_position);
		}
	}
	#endif
}
void PhgBinCompPetDAOld(PHG_TrackingPhoton *bluePhotonPtr, PHG_TrackingPhoton *pinkPhotonPtr, 
		double *anglePtr, double *distancePtr)
{
	double			theta1;			/* Temp var for computation */
	double			theta2;			/* Temp var for computation */
	double			minAngle;		/* Temp var for computation */
	double			maxAngle;		/* Temp var for computation */
	double			outsideRadius;	/* The outside radius of the cylinder we have reached */
		
	/* Compute theta1 for the blue photon */
	theta1 = atan2(bluePhotonPtr->location.y_position,
			  bluePhotonPtr->location.x_position);

	/* Compute theta2 */
	theta2 = atan2(pinkPhotonPtr->location.y_position,
		  pinkPhotonPtr->location.x_position);
		
	/* Compute the angle */
	*anglePtr = (theta1 + theta2 + PHGMATH_PI)/2;
	
	/* Check for range */
	if (*anglePtr < -(PHGMATH_PI_DIV2)) {
		*anglePtr = *anglePtr + PHGMATH_PI;
	}
	else if (*anglePtr >= (PHGMATH_PI_DIV2)) {
		*anglePtr = *anglePtr - PHGMATH_PI;
	}
	
	/* Order theta1 and theta2 for computation of distance */
	if (theta2 < theta1) {
		minAngle = theta2;
		maxAngle = theta1;
	}
	else {
		minAngle = theta1;
		maxAngle = theta2;
	}
	
	/* Get the outside radius */
	if (PHG_IsDetectOnTheFly()) {
		outsideRadius = DetGtOutsideRadius();
	}
	else if (PHG_IsCollimateOnTheFly()) {
		outsideRadius = ColGtOutsideRadius();
	}
	else {
		outsideRadius = CylPosGetTargetRadius();
	}
	
	/* Compute the distance */
	*distancePtr = fabs(outsideRadius* cos((maxAngle-minAngle)/2));
	
	/* Check for range */
	if ((maxAngle <= 0) ||
			((minAngle < 0) && ((-1 *sin(minAngle)) > sin(maxAngle)))) {
			
		*distancePtr *= -1;
	}
	else if (PhgMathRealNumAreEqual(-maxAngle, minAngle, 
			 -7, 0, 0, 0)) {
		if (cos(maxAngle) < 0)
			*distancePtr *= -1;
	}

}

/*********************************************************************************
*
*			Name:			PhgBinPrintReport
*
*			Summary:		Print a report of the final statistics.
*
*			Arguments:
*				PHG_BinParamsTy		*binParams		- User defined binning parameters.
*				PHG_BinFieldsTy *binFields	- Binning statistics
*			Function return: None.
*
*********************************************************************************/
void PhgBinPrintReport(	PHG_BinParamsTy *binParams,
						PHG_BinFieldsTy *binFields)
						
{
	/* Ignore if in error condition */
	if (ErIsInError())
		return;
		
	/* Tell the user whats happening */
	LbInPrintf("\n\n***************** Processing Binning Data *************\n");
	
	/* Print out our statistics */
	if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ){
		LbInPrintf("\nTotal blue photons reaching the binning module in this simulation = %lld", binFields->TotBluePhotons);
		LbInPrintf("\nTotal pink photons reaching the binning module in this simulation = %lld", binFields->TotPinkPhotons);
		LbInPrintf("\n\nTotal coincidences reaching the binning module in this simulation = %lld", binFields->NumCoincidences);
		LbInPrintf("\nTotal accepted blue photons in this simulation = %lld", binFields->AccBluePhotons);
		LbInPrintf("\nTotal accepted pink photons in this simulation = %lld", binFields->AccPinkPhotons);
		LbInPrintf("\nTotal accepted coincidences in this simulation = %lld",
			binFields->NumAcceptedCoincidences - binFields->CountImgHdr.H.NumPhotons);
		LbInPrintf("\nSum of accepted coincidence weights in this simulation = %3.6e",
			binFields->AccCoincidenceWeight - binFields->StartAccCoincidenceWeight);
		LbInPrintf("\nSum of accepted coincidence squared weights in this simulation = %3.6e",
			binFields->AccCoincidenceSquWeight - binFields->StartAccCoincidenceSquWeight);
			
		/* Print out image totals if we had a non-zero starting value */
		if ((binFields->WeightImgHdr.H.NumPhotons != 0) || (binFields->WeightSquImgHdr.H.NumPhotons != 0) ||(binFields->CountImgHdr.H.NumPhotons != 0)) {
			LbInPrintf("\nTotal accepted coincidences in image = %lld",
				binFields->NumAcceptedCoincidences);
			LbInPrintf("\nSum of accepted coincidence weights in image = %3.6e.", binFields->AccCoincidenceWeight);
			LbInPrintf("\nSum of accepted coincidence squared weights in image = %3.6e.", binFields->AccCoincidenceSquWeight);
		}
		
		/* Report the quality factor */
		if ((binFields->NumAcceptedCoincidences != 0) && (binFields->AccCoincidenceSquWeight != 0)){
			LbInPrintf("\nQuality factor in image data = %3.2e",
				PHGMATH_Square(binFields->AccCoincidenceWeight)/(binFields->NumAcceptedCoincidences * binFields->AccCoincidenceSquWeight));
		}
		LbInPrintf("\n");
	}
	else {
		LbInPrintf("\nTotal photons reaching the binning module in this simulation = %lld", binFields->TotBluePhotons);
		LbInPrintf("\nTotal accepted photons in this simulation = %lld",
			binFields->AccBluePhotons - binFields->CountImgHdr.H.NumPhotons);
		LbInPrintf("\nSum of accepted weights in this simulation = %3.6e",
			binFields->AccCoincidenceWeight - binFields->StartAccCoincidenceWeight);
		LbInPrintf("\nSum of accepted squared weights in this simulation = %3.6e",
			binFields->AccCoincidenceSquWeight - binFields->StartAccCoincidenceSquWeight);

		/* Print out image totals if we had a non-zero starting value */
		if ((binFields->WeightImgHdr.H.NumPhotons != 0) || (binFields->WeightSquImgHdr.H.NumPhotons != 0) ||(binFields->CountImgHdr.H.NumPhotons != 0)) {
			LbInPrintf("\nTotal accepted photons in image = %lld",
				binFields->AccBluePhotons);
			LbInPrintf("\nSum of accepted weights in image = %3.6e",
				binFields->AccCoincidenceWeight);
			LbInPrintf("\nSum of accepted squared weights in image = %3.6e",
				binFields->AccCoincidenceSquWeight);
		}
		/* Report the quality factor */
		if ((binFields->AccBluePhotons != 0) && (binFields->AccCoincidenceSquWeight != 0)){
			LbInPrintf("\nQuality factor in image data = %3.2e",
				PHGMATH_Square(binFields->AccCoincidenceWeight)/(binFields->AccBluePhotons * binFields->AccCoincidenceSquWeight));
		}
		LbInPrintf("\n");
	}
			
	/* Print out history report */
	if (binFields->historyFileHk.doCustom) {
		LbInPrintf("\nHistory file report for binning module");
		PhoHFilePrintReport(&binFields->historyFileHk);
	}
}

/*********************************************************************************
*
*			Name:			PhgBinTerminate
*
*			Summary:		Process the bin buffers and terminate the module.
*
*			Arguments:
*				PHG_BinParamsTy	*binParams	 - User defined binning parameters.
*				PHG_BiDataTy	*binData	 - Storage for binned data.
*				PHG_BinFieldsTy *binFields	 - Various binning information.
*			Function return: None.
*
*********************************************************************************/
void PhgBinTerminate(PHG_BinParamsTy *binParams, PHG_BinDataTy *binData, PHG_BinFieldsTy *binFields)
{
	Boolean			okay = false;	/* Process Flag */
	LbUsFourByte	i;				/* Loop control */
	float			singlePrec;		/* For converting output */
	
	do { /* Process Loop */


		/* If we didn't initialize, skip all of this */
		if (binFields->IsInitialized == false) {
			okay = true;
			break;
		}
		
		/* Call the user binning routine */
		if (BinUsrTerminateFPtr) {
			(*BinUsrTerminateFPtr)(binParams, binData);
		}
		
		/* Process the count image */
		if (binParams->doCounts == true) {
		
			if (!ErIsInError()) {
				
					LbInPrintf("\nWriting count image.");
	
				/* Seek back to beginning */
				if (fseek(binFields->CountFile, 0, SEEK_SET) != 0){
					ErStFileError("Unable to seek to beginning of count image (PhgBinTerminate).");
					break;
				}
	
				/* Update the number of photons */
				if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
					binFields->CountImgHdr.H.NumPhotons = binFields->NumAcceptedCoincidences;
				}
				else {
					binFields->CountImgHdr.H.NumPhotons = binFields->AccBluePhotons;
				}
				
				binFields->CountImgHdr.H.weightSum = binFields->AccCoincidenceWeight;
				binFields->CountImgHdr.H.weightSquSum = binFields->AccCoincidenceSquWeight;
			
				/* Update the number of simulations in this file */
				if (binFields->CountImgHdr.H.NumSimulations == 0)
					binFields->CountImgHdr.H.NumSimulations++;
				
				/* Update the image header */ 
				if (PhgHdrUpHeader(binFields->CountFile, &binFields->CountImgHdr, &binFields->CountImgHdrHk) == false) {
					break;
				}
	
				/* Write the image buffer */ 
				if (fwrite(binData->countImage, binParams->countImageSize, 1, binFields->CountFile) != 1) {
					ErStFileError("Unable to write count image file.");
					break;
				}
			}				
		}

		/* Process the weight image */
		if (binParams->doWeights == true) {
		
		
			if (!ErIsInError()) {
			
				LbInPrintf("\nWriting weight image.");

				/* Update the number of photons */
				if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
					binFields->WeightImgHdr.H.NumPhotons = binFields->NumAcceptedCoincidences;
				}
				else {
					binFields->WeightImgHdr.H.NumPhotons = binFields->AccBluePhotons;
				}
				
				/* Update the number of simulations in this file */
				if (binFields->WeightImgHdr.H.NumSimulations == 0)
					binFields->WeightImgHdr.H.NumSimulations++;
	
				/* Update sum of weights */
				binFields->WeightImgHdr.H.weightSum = binFields->AccCoincidenceWeight;
				binFields->WeightImgHdr.H.weightSquSum = binFields->AccCoincidenceSquWeight;
				
				/* Update the image header */ 
				if (PhgHdrUpHeader(binFields->WeightFile, &binFields->WeightImgHdr, &binFields->WeightImgHdrHk) == false) {
					break;
				}
				
				/* Write the image buffer */
				if ((binParams->sumAccordingToType == false) && (binParams->weight_image_type == PHG_BIN_WEIGHT_TYPE_R4)) {
					
					/* Tell user this may take a while */
					LbInPrintf("\n\tConverting weight file to single precision, this may take a while.");
					for (i = 0; i < binParams->weightImageSize/sizeof(double); i++) {
						singlePrec = ((double *)binData->weightImage)[i];
						if (fwrite((void *) &singlePrec, 1, sizeof(float), binFields->WeightFile) != sizeof(float)) {
							ErStFileError("\nUnable to write weight image file.");
							goto FAIL;
						}
					}
 				}
 				else {
					if (fwrite(binData->weightImage, binParams->weightImageSize, 1, binFields->WeightFile) != 1) {
						ErStFileError("\nUnable to write weight image file.");
						break;
					}
				}
			}
		}
		

		/* Process the weight squared image */
		if (binParams->doWeightsSquared == true) {
		
			if (!ErIsInError()) {
			
				LbInPrintf("\nWriting weight squared image.");

				/* Update the weight-squared of photons */
				if ( !(PHG_IsSPECT() || binParams->isBinPETasSPECT) ) {
					binFields->WeightSquImgHdr.H.NumPhotons = binFields->NumAcceptedCoincidences;
				}
				else {
					binFields->WeightSquImgHdr.H.NumPhotons = binFields->AccBluePhotons;
				}
				
				/* Update the number of simulations in this file */
				if (binFields->WeightSquImgHdr.H.NumSimulations == 0)
					binFields->WeightSquImgHdr.H.NumSimulations++;
	
				/* Update sum of weights squared */
				binFields->WeightSquImgHdr.H.weightSquSum = binFields->AccCoincidenceSquWeight;
				binFields->WeightSquImgHdr.H.weightSum = binFields->AccCoincidenceWeight;
				
				/* Update the image header */ 
				if (PhgHdrUpHeader(binFields->WeightSquFile, &binFields->WeightSquImgHdr, &binFields->WeightSquImgHdrHk) == false) {
					break;
				}
	
				/* Write the image buffer */ 
				if ((binParams->sumAccordingToType == false) && (binParams->weight_image_type == PHG_BIN_WEIGHT_TYPE_R4)) {
					
					/* Tell user this may take a while */
					LbInPrintf("\n\tConverting weight squared file to single precision, this may take a while.");
					for (i = 0; i < binParams->weightImageSize/sizeof(double); i++) {
						singlePrec = ((double *)binData->weightSquImage)[i];
						if (fwrite((void *) &singlePrec, 1, sizeof(float), binFields->WeightSquFile) != sizeof(float)) {
							ErStFileError("\nUnable to write weight square image file.");
							goto FAIL;
						}
					}
 				}
 				else {
					if (fwrite(binData->weightSquImage, binParams->weightSquImageSize, 1, binFields->WeightSquFile) != 1) {
						ErStFileError("\nUnable to write weight squared image file.");
						break;
					}
				}
			}
		}
	
		okay = true;
		FAIL:;
	} while (false);
	
	/* Clear our initialization flags */
	binFields->IsInitialized = false;
	binFields->ParamsIsInitialized = false;

	/* Free the header hook */
	PhgHdrFrHeader(&binFields->CountImgHdrHk);
	
	/* Close the image file */
	if (binFields->CountFile != 0)
		fclose(binFields->CountFile);
	
	/* Free the buffer memory  */
	if (binData->countImage != 0)
		LbMmFree(&(binData->countImage));
	
	/* Free the header hook */
	PhgHdrFrHeader(&binFields->WeightSquImgHdrHk);
	
	/* Close the image file */
	if (binFields->WeightSquFile != 0)
		fclose(binFields->WeightSquFile);
	
	/* Free the buffer memory  */
	if (binData->weightSquImage != 0)
		LbMmFree(&(binData->weightSquImage));
	
		
	/* Free the header hook */
	PhgHdrFrHeader(&binFields->WeightImgHdrHk);

	/* Close the image file */
	if (binFields->WeightFile != 0)
		fclose(binFields->WeightFile);
	
	/* Free the buffer memory  */
	if (binData->weightImage != 0)
		LbMmFree(&(binData->weightImage));
	
	/* Do error handling here */
	if (!okay) {
		ErHandle("Failed to process binning files.", false);
	}
	
	if (!ErIsInError())
		LbInPrintf("\n\n***************** Finished Processing Binning Data *************\n");

}
#undef PHG_BIN
