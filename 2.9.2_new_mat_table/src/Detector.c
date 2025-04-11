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
*			Module Name:		Detector.c
*			Revision Number:	2.4
*			Date last revised:	23 July 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	11 July, 1994
*
*			Module Overview:	Simulates detector functionality.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:
*				DetInitialize
*				DetGaussTimeBlur
*				DetIsSimplePet
*				DetIsSimpleSpect
*				DetIsUNCSpect
*				DetIsPlanar
*				DetIsCylindrical
*				DetIsDualHeaded
*				DetGtOutsideRadius
*				DetPETPhotons
*				DetSPECTPhotons
*				DetTerminate
*
*			Global variables defined:		none
*
*********************************************************************************
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
*			Revision date:		26 January 2012
*
*			Revision description:	Changed form of user functions to pointers
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		R Harrison
*
*			Revision date:		19 Jan 2012
*
*			Revision description:	Added print out of detected photon
*					positioning option for block detectors.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		1 December 2009
*
*			Revision description:	Corrected to allow multiple history files.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		January 2008
*
*			Revision description:	Added block detectors.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		July 2007
*
*			Revision description:	Changed input of block positioning 
*						within rings to cylindrical coords.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		August 2006
*
*			Revision description:	Added block detector parameter display and
*						other block functionality.
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:
*					- Added time-of-flight blurring, including function
*					DetGaussTimeBlur
*					- Zeroed detected event counters to fix bug in history file
*					processing.
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		February 2004
*
*			Revision description:	Added the DetIs... functions.
*
*********************************************************************************/

#define DETECTOR


#include <stdio.h>
#include <memory.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbMemory.h"
#include "LbParamFile.h"
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
#include "ProdTbl.h"
#include "SubObj.h"
#include "EmisList.h"
#include "PhoHStat.h"
#include "PhoHFile.h"
#include "ColUsr.h"
#include "UNCCollimator.h"
#include "Collimator.h"
#include "ColSlat.h"
#include "DetCylinder.h"
#include "Detector.h"
#include "DetUsr.h"
#include "DetPlanar.h"
#include "DetGeometric.h"
#include "DetBlock.h"
#include "phg.h"
#include "PhgBin.h"


/* LOCAL CONSTANTS */
#define GAUSS_MAGIC_NUM	235.4820045		/* Ratio between FWHM and sigma for a Guassian
											distribution * 100. This number is used to
											convert from a user specified percentage
											blur, to a standard deviation for the
											blurring routine.
										*/


									
/* Local Types */


/* Local Prototypes */
void 			detSimplePET(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					DetectedPhotonsTy *detPhotonsPtr);
void 			detPolyPET(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					DetectedPhotonsTy *detPhotonsPtr);
void 			detCylnPET(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
					PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
					DetectedPhotonsTy *detPhotonsPtr);
void			detSimpleSPECT(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
					DetectedPhotonsTy *detPhotonsPtr);
void 			detCylindricalSPECT(PHG_Decay *decayPtr,
					PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
					DetectedPhotonsTy *detPhotonsPtr);

/* Local Global Variables */
/* 	(None for now) */	


/*********************************************************************************
*
*			Name:			DetInitialize
*
*			Summary:		Initialize the binning module.
*
*			Arguments:
*				Boolean		doHistory	- Do we create a history file (overides all)
*			Function return: TRUE unless an error occurs.
*
*********************************************************************************/
Boolean DetInitialize(Boolean doHistory)
{
	Boolean			okay = false;	/* Process Flag */
	LbUsFourByte	ringIndex;		/* Ring index for initialization */
	LbUsFourByte	layerIndex;		/* Layer index for initialization */
	
	do {
		/* Clear processing  variables */
		detData[DetCurParams].detTotBluePhotons = 0;
		detData[DetCurParams].detTotPinkPhotons = 0;
		detData[DetCurParams].detTotAccBluePhotons = 0;
		detData[DetCurParams].detTotAccPinkPhotons = 0;
		detData[DetCurParams].detTotPhotonsDepositingEnergy = 0;
		detData[DetCurParams].detTotPhotonsAbsorbed = 0;
		detData[DetCurParams].detTotWtAbsorbed = 0.0;
		detData[DetCurParams].detTotForcedAbsorptions = 0;
		detData[DetCurParams].detTotWtForcedAbsorbed = 0.0;
		detData[DetCurParams].detTotFirstTimeAbsorptions = 0;
		detData[DetCurParams].detTotWtFirstTimeAbsorbed = 0.0;
		detData[DetCurParams].detTotPhotonsPassingThrough = 0;
		detData[DetCurParams].detTotReachingCrystal = 0;
		detData[DetCurParams].detWeightAdjusted = 0.0;
		detData[DetCurParams].detNumReachedMaxInteractions = 0;

		#ifdef PHG_DEBUG
		detData[DetCurParams].detCylCountCohInteractions = 0;
		detData[DetCurParams].detCylCountInteractions = 0;
		#endif
		
		/* Clear our detector type, so we can tell if it gets initialized */
		DetRunTimeParams[DetCurParams].DetectorType = DetEn_DetType_NULL;
		DetRunTimeParams[DetCurParams].Initialized = false;
		DetRunTimeParams[DetCurParams].DoForcedInteraction = false;
		DetRunTimeParams[DetCurParams].PhotonTimeFWHM = 0.0;
		DetRunTimeParams[DetCurParams].EnergyResolutionPercentage = -1.0;
		DetRunTimeParams[DetCurParams].ReferenceEnergy = 0.0;
		DetRunTimeParams[DetCurParams].DoRandomsProcessing = false;
		DetRunTimeParams[DetCurParams].TriplesMethod = DetEn_DeleteTriples;
		DetRunTimeParams[DetCurParams].CoincidenceTimingWindowNS = 0.0;
		
		/* Get the detector parameters */
		if (DetGetRunTimeParams() == false) {
			break;
		}
		
		/* time-of-flight is only available for PET simulations */
		if ((DetRunTimeParams[DetCurParams].PhotonTimeFWHM != 0.0) && (PHG_IsSPECT())) {
			PhgAbort("You may only specify a time-of-flight resolution for PET.", false);
		}
		
		/* If time-of-flight resolution must be >= 0 */
		if ( DetRunTimeParams[DetCurParams].PhotonTimeFWHM < 0.0 ) {
			PhgAbort("Time-of-flight resolution must be > 0.", false);
		}
		
		/* If they changed the energy resolution, verify they supplied a reference energy */
		if ((DetRunTimeParams[DetCurParams].EnergyResolutionPercentage != -1.0) && (DetRunTimeParams[DetCurParams].ReferenceEnergy == 0)) {
			PhgAbort("If you specify an energy resolution in the detector parameters, you must also specify a reference energy.", false);
		}
		
		/* Make sure that they don't misunderstand the energy resolution via specifying a value too high */
		if (DetRunTimeParams[DetCurParams].EnergyResolutionPercentage > 50.0) {
			PhgAbort("An energy resolution of greater than 50.0% doesn't make sense, it's way too bad.", false);
		}
		
		/* Override do history flag if they want it turned off */
		if (doHistory == false)
			DetRunTimeParams[DetCurParams].DoHistory = doHistory;
		
		/* Set our initialization flag */
		DetRunTimeParams[DetCurParams].Initialized = true;
		
		/* Call the user initialization routine */
		if (DetUsrInitializeFPtr) {
			(*DetUsrInitializeFPtr)(&(DetRunTimeParams[DetCurParams]));
		}
		
		/* Create history file if were are supposed to */
		if (DET_IsDoHistory()) {
			if (PhoHFileCreate(DetRunTimeParams[DetCurParams].DetHistoryFilePath,
					DetRunTimeParams[DetCurParams].DetHistoryParamsFilePath, PhoHFileEn_DET, 
					&(detData[DetCurParams].detHistFileHk)) == false) {
					
				sprintf(detErrStr,"Failed to create history file specified in detector parameters file named:\n"
					"'%s'\n"
					"The custom parameters file name is: '%s'\n"
					" (DetInitialize)",
					DetRunTimeParams[DetCurParams].DetHistoryFilePath, DetRunTimeParams[DetCurParams].DetHistoryParamsFilePath);
				PhgAbort(detErrStr, false);
			}
		}

		/* Do type specific initialization */
		switch(DetRunTimeParams[DetCurParams].DetectorType){
		
			case	DetEn_simple_pet:
			case	DetEn_simple_spect:
				if ((PHG_IsCollimateOnTheFly() == true) && (ColIsSlat() == true)) {
					PhgAbort("You must use a 'real' detector module, ie PLANAR, when using slat colliamators.", false);
				}
				
			case 	DetEn_Polygonal:
				
				break;
			
			case	DetEn_DualHeaded:
			case 	DetEn_Planar:
				
				/* Verify we are not "inside" collimator radius, note that ColGtOutsideRadius works whether collimating or not */
				if (ColGtOutsideRadius() > DetRunTimeParams[DetCurParams].PlanarDetector.InnerRadius) {
					PhgAbort("Your collimator (or target cylinder) radius is greater than the detector inner radius", false);
				}
			
				/* Initialize the inner bounding cylinder */
				detData[DetCurParams].detInBoundCyl.radius = DetRunTimeParams[DetCurParams].PlanarDetector.InnerRadius;
				
				/* If collimating on the fly the z limits are defined by the collimator unless it is slat collimation */
				if (PHG_IsCollimateOnTheFly() && ColIsSlat() == false) {
				
					ColGtZLimits(&detData[DetCurParams].detInBoundCyl.zMin, &detData[DetCurParams].detInBoundCyl.zMax);
					
					/* Verify we match the collimator parameters */ 
					if ((detData[DetCurParams].detInBoundCyl.zMax - detData[DetCurParams].detInBoundCyl.zMin) != DetRunTimeParams[DetCurParams].PlanarDetector.AxialLength) {
						PhgAbort("When doing collimation and detection the detector axial length must match the length of the collimator [specified in the collimator parameter file using minZ and maxZ]", false);
					}
	
				}
				else {
					detData[DetCurParams].detInBoundCyl.zMin = -(DetRunTimeParams[DetCurParams].PlanarDetector.AxialLength/2);
					detData[DetCurParams].detInBoundCyl.zMax = DetRunTimeParams[DetCurParams].PlanarDetector.AxialLength/2;
				}
				
				detData[DetCurParams].detInBoundCyl.centerX = 0.0;
				detData[DetCurParams].detInBoundCyl.centerY = 0.0;
				
				/* Compute the transaxial limit */
				detData[DetCurParams].detPlnTransLimit = DetRunTimeParams[DetCurParams].PlanarDetector.TransaxialLength/2;
								
				/* Compute depth of the detector */
				detData[DetCurParams].detPlnDetectorDepth = 0.0;
				for (layerIndex = 0; layerIndex < DetRunTimeParams[DetCurParams].PlanarDetector.NumLayers; layerIndex++){
					detData[DetCurParams].detPlnDetectorDepth += DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[layerIndex].LayerDepth;
				}
				
				/* Set range to constant value, later it may be parameterized */
				detData[DetCurParams].detViewSize = 180;
				
				
				/* Initiailize "big cylinder */
				detData[DetCurParams].detPlnrBigCylinder = detData[DetCurParams].detInBoundCyl;
				detData[DetCurParams].detPlnrBigCylinder.radius =
					PHGMATH_SquareRoot(PHGMATH_Square(detData[DetCurParams].detInBoundCyl.radius)
					+ PHGMATH_Square(detData[DetCurParams].detPlnTransLimit));

				/* Compute limiting angle for hitting the detector */
				detData[DetCurParams].detPlnrAngularCoverage = atan(DetRunTimeParams[DetCurParams].PlanarDetector.TransaxialLength/
					(2 * DetRunTimeParams[DetCurParams].PlanarDetector.InnerRadius));
				
				/* Shift theta into range of [0,pi] */
				detData[DetCurParams].detPlnrAngularCoverage = fabs(detData[DetCurParams].detPlnrAngularCoverage);
				
				/* Convert to degrees */
				detData[DetCurParams].detPlnrAngularCoverage = PHGMATH_DegreesFromRadians(detData[DetCurParams].detPlnrAngularCoverage);

				/* Check for run-time collimation and the presence of num_views */
				if ((PHG_IsCollimateOnTheFly() == false) && (DetRunTimeParams[DetCurParams].PlanarDetector.NumViews == 0)) {
					ErStGeneric("If you don't use the collimator, you must specify num_views in the "
						"detector parameters (DetInitialize).");
					goto FAILURE;
				}
	
				/* If they specified a number of views and we are collimating on the fly, verify
				   that the number of views is consistent
				   8/3/01 - to support singles check is changed to only occur when doing SPECT/UNC: SDV
				*/
				if (DetRunTimeParams[DetCurParams].PlanarDetector.NumViews != 0) {
					if (PHG_IsCollimateOnTheFly() && PHG_IsSPECT() && (ColIsUNC() == true))
						if ((LbFourByte)(ColGtNumViews()) != DetRunTimeParams[DetCurParams].PlanarDetector.NumViews) {
							sprintf(detErrStr, "Num views in collimator parameters must match"
								" the number of views in the detector parameters file,"
								" (collimator views = %ld) (detector views = %ld) (DetInitialize)",
								(unsigned long)ColGtNumViews(), 
								(long)DetRunTimeParams[DetCurParams].PlanarDetector.NumViews);
							ErStGeneric(detErrStr);
							goto FAILURE;					
					}
				}
			
				/* If they did not specify the number of views and we are collimating on the fly,
					set num views to that of the collimator for convenience
				   8/3/01 - to support singles check is changed to only occur when doing SPECT/UNC: SDV
				*/
				if ((PHG_IsCollimateOnTheFly() == true) && (DetRunTimeParams[DetCurParams].PlanarDetector.NumViews == 0)
						&& (ColIsUNC() == true)) {
					DetRunTimeParams[DetCurParams].PlanarDetector.NumViews = ColGtNumViews();
				}
				/* If they specify continuous rotation, verify min and maximum angle are zero and 360
				*/
				if (DetRunTimeParams[DetCurParams].PlanarDetector.NumViews == -1) {
					if (PHG_IsPET()) {
						if ((PhgMathRealNumAreEqual(DetRunTimeParams[DetCurParams].PlanarDetector.MinAngle,0.0,1,0,0, 0) == false) ||
								(PhgMathRealNumAreEqual(DetRunTimeParams[DetCurParams].PlanarDetector.MaxAngle,PHGMATH_PI,1,0,0,0) == false)) {
							sprintf(detErrStr, "When doing continuous detector angles, you must specify a minimum angle of zero and a maximum angle of 180 degrees.  Your 'num_views' is set to -1 indicating you want continuous angles.");
							ErStGeneric(detErrStr);
							goto FAILURE;
						}
					}
					else {
						if ((PhgMathRealNumAreEqual(DetRunTimeParams[DetCurParams].PlanarDetector.MinAngle,0.0,1,0,0, 0) == false) ||
								(PhgMathRealNumAreEqual(DetRunTimeParams[DetCurParams].PlanarDetector.MaxAngle,PHGMATH_2PI,1,0,0,0) == false)) {
							sprintf(detErrStr, "When doing continuous detector angles, you must specify a minimum angle of zero and a maximum angle of 360 degrees.  Your 'num_views' is set to -1 indicating you want continuous angles.");
							ErStGeneric(detErrStr);
							goto FAILURE;
						}
					}
				}
				
				/* Check min/max angle */
				if (PHG_IsPET() && (DetRunTimeParams[DetCurParams].PlanarDetector.MaxAngle > PHGMATH_PI)) {
					ErStGeneric("For PET, the maximum detector angle can not be greater than 180 degrees");
					goto FAILURE;
				}
				
				/* Compute planar delta, note that if num views == 0, delta is not used*/
				if (DetRunTimeParams[DetCurParams].PlanarDetector.NumViews != 0) {
					detData[DetCurParams].detPlnrDelta = 
						(DetRunTimeParams[DetCurParams].PlanarDetector.MaxAngle - DetRunTimeParams[DetCurParams].PlanarDetector.MinAngle)
						/DetRunTimeParams[DetCurParams].PlanarDetector.NumViews;
				}
								
				/* If we are collimating on the fly with a slat collimator then we need to
					initialize
					
				*/
				if (PHG_IsCollimateOnTheFly() && (ColIsSlat() == true)) {
					ColSlatSetParamsFromDetector(DetRunTimeParams[DetCurParams].PlanarDetector.MinAngle,
						detData[DetCurParams].detInBoundCyl.zMin,
						detData[DetCurParams].detInBoundCyl.zMax);
				}
				break;
				
			case 	DetEn_Cylindrical:
				/*	Compute the bounding cylinders 
					NOTE: I am assuming that the rings have the same inner and outer radius. This
					should be enforced through checking the parameters until some time at which
					we support changes on a per ring basis.
				*/
				
				/* Verify we are not "inside" collimator radius, note that ColGtOutsideRadius works whether collimating or not */
				if (ColGtOutsideRadius() > DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[0].LayerInfo[0].InnerRadius) {
					PhgAbort("Your collimator (or target cylinder) radius is greater than the detector inner radius", false);
				}
				detData[DetCurParams].detInBoundCyl.radius = DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[0].LayerInfo[0].InnerRadius;
				detData[DetCurParams].detInBoundCyl.centerX = 0.0;
				detData[DetCurParams].detInBoundCyl.centerY = 0.0;
				
				detData[DetCurParams].detInBoundCyl.zMin = DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[0].MinZ;
				detData[DetCurParams].detInBoundCyl.zMax = DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[0].MaxZ;

				detData[DetCurParams].detOutBoundCyl = detData[DetCurParams].detInBoundCyl;
				for (ringIndex = 0; ringIndex < DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings;
						ringIndex++){
					for (layerIndex = 0; layerIndex < DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[ringIndex].NumLayers; layerIndex++){
						
						if (detData[DetCurParams].detOutBoundCyl.radius < DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[ringIndex].LayerInfo[layerIndex].OuterRadius) {
							detData[DetCurParams].detOutBoundCyl.radius = DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[ringIndex].LayerInfo[layerIndex].OuterRadius;
						}
					}
				}
				break;
				
			case 	DetEn_Block:
				/*	Compute the bounding cylinders 
					NOTE: These cylinders should probably not be used too much!
					As the cylinders may be different for different rings in the
					tomograph, we compute these cylinders as the minimum (for the
					inner cylinder) or maximum (for the outer cylinder) of the
					individual ring cylinders.
				*/
				
				detData[DetCurParams].detInBoundCyl.radius = MAXFLOAT;
				for (ringIndex = 0; ringIndex < DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings; ringIndex++) {
					
					if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[ringIndex].XInnerRadius <
							detData[DetCurParams].detInBoundCyl.radius) {
						detData[DetCurParams].detInBoundCyl.radius =
								DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[ringIndex].XInnerRadius;
					}
					
					if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[ringIndex].YInnerRadius <
							detData[DetCurParams].detInBoundCyl.radius) {
						detData[DetCurParams].detInBoundCyl.radius =
								DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[ringIndex].YInnerRadius;
					}
					
				}
				
				detData[DetCurParams].detOutBoundCyl.radius = detData[DetCurParams].detInBoundCyl.radius;
				for (ringIndex = 0; ringIndex < DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings; ringIndex++) {
					
					if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[ringIndex].XOuterRadius >
							detData[DetCurParams].detOutBoundCyl.radius) {
						detData[DetCurParams].detOutBoundCyl.radius =
								DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[ringIndex].XOuterRadius;
					}
					
					if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[ringIndex].YOuterRadius <
							detData[DetCurParams].detOutBoundCyl.radius) {
						detData[DetCurParams].detOutBoundCyl.radius =
								DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[ringIndex].YOuterRadius;
					}
					
				}
				
				/* Verify we are not "inside" collimator radius, note that ColGtOutsideRadius works whether collimating or not */
				if (ColGtOutsideRadius() > detData[DetCurParams].detInBoundCyl.radius) {
					PhgAbort("Your collimator (or target cylinder) radius is greater than the detector inner radius", false);
				}
				
				detData[DetCurParams].detInBoundCyl.centerX = 0.0;
				detData[DetCurParams].detInBoundCyl.centerY = 0.0;
				
				detData[DetCurParams].detInBoundCyl.zMin = 
					DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[0].MinZ;
				detData[DetCurParams].detInBoundCyl.zMax = 
					DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings - 1].MaxZ;
				
				/* Do the zone initialization */
				if (DetBlocDivideZones(15) == 0) {	/* Maximum of 15 blocks per zone (per ring) */
					/* Memory failure allocating zone database */
					goto FAILURE;
				}
				
				break;
				
			default:
				break;
		}

		okay = true;
		FAILURE:;
	} while (false);
	
	/* If we failed, free any potential memory allocated */
	if (!okay) {
	
		switch(DetRunTimeParams[DetCurParams].DetectorType){
		
			case 	DetEn_Polygonal:
				{
					LbUsFourByte ringIndex;
					LbUsFourByte blockIndex;
					
					if (DetRunTimeParams[DetCurParams].PolyDetector.RingInfo != 0) {
						for (ringIndex = 0; ringIndex < DetRunTimeParams[DetCurParams].PolyDetector.NumRings; ringIndex++){
						
							if (DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[ringIndex].BlockInfo != 0){
								for (blockIndex = 0; blockIndex < DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[ringIndex].NumBlocks; blockIndex++){
									
									if (DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[ringIndex].BlockInfo[blockIndex].LayerInfo != 0){
										LbMmFree((void **)&DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[ringIndex].BlockInfo[blockIndex].LayerInfo);
									}
								}
								
								LbMmFree((void **) &DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[ringIndex].BlockInfo);
							}
						}
						LbMmFree((void **) &DetRunTimeParams[DetCurParams].PolyDetector.RingInfo);
					}
				}
				
				break;
			
			case	DetEn_DualHeaded:
			case 	DetEn_Planar:
			
				if (DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo != 0) {
					LbMmFree((void **)&DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo);
				}
				break;
				
			case 	DetEn_Cylindrical:
				{
					LbUsFourByte ringIndex;
				
					if (DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo != 0) {
						
						for (ringIndex = 0; ringIndex < DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings; ringIndex++){
								
							if (DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[ringIndex].LayerInfo != 0) {
								LbMmFree((void **)&DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[ringIndex].LayerInfo);
							}
						}
						LbMmFree((void **) &DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo);
					}
				}
				break;
				
			case 	DetEn_Block:
				{
					LbUsFourByte rNum, bNum, lNum;
				
					if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo != 0) {
						
						for (rNum = 0; rNum < DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings; rNum++){
								
							if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo != 0) {
								
								for (bNum = 0; bNum < DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].NumBlocks; bNum++){
									
									if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo != 0) {
										
										for (lNum = 0; lNum < DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].NumLayers; lNum++){
										
											if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo[lNum].YChanges != 0) {
												LbMmFree((void **)&DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo[lNum].YChanges);
											}
										
											if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo[lNum].ZChanges != 0) {
												LbMmFree((void **)&DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo[lNum].ZChanges);
											}
										
											if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo[lNum].ElementInfo != 0) {
												LbMmFree((void **)&DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo[lNum].ElementInfo);
											}
										
										}
										
										LbMmFree((void **)&DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo);
										
									}
									
								}
								
								LbMmFree((void **)&DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo);
								
							}
							
						}
						
						LbMmFree((void **) &DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo);
						
					}
					
					DetBlocFreeData();
				}
				
				break;
				
			default:
				break;
		}

	}
	return (okay);
}

/*********************************************************************************
*
*			Name:			DetPrintParams
*
*			Summary:		Print the simulation parameters.
*
*			Arguments:
*			Function return: None.
*
*********************************************************************************/
void DetPrintParams(void)
{
	LbUsFourByte	layerIndex, r, l;	/* LCV's */
	LbUsFourByte	ringIndex;			/* Loop control variable */
	LbUsFourByte	blockIndex;			/* Loop control variable */
	LbUsFourByte	elementIndex;		/* Loop control variable */
	LbFourByte		yChangeIndex;		/* Loop control variable */
	LbFourByte		zChangeIndex;		/* Loop control variable */
	DetBlockTomoRingTy		*curRingInfoPtr;		/* the current ring info (shortens lines in code!) */
	DetBlockTomoBlockTy		*curBlockInfoPtr;		/* the current block info (shortens lines in code!) */
	DetBlockTomoLayerTy		*curLayerInfoPtr;		/* the current layer info (shortens lines in code!) */
	Boolean         detShowBlocks = false;	/* should internal block parameters be displayed? */
	
	for (DetCurParams = 0; DetCurParams < DetNumParams; DetCurParams++) {
	/* Do type specific reporting */
	LbInPrintf("\n\n Detector Parameters for '%s'\n", PhgRunTimeParams.PhgDetectorParamsFilePath[DetCurParams]);
	
	switch(DetRunTimeParams[DetCurParams].DetectorType){
	
		case 	DetEn_Polygonal:
			LbInPrintf("\n\tDetector type is polygonal\n");
			break;
		
		case	DetEn_DualHeaded:
		case 	DetEn_Planar:
		
			if (DetRunTimeParams[DetCurParams].DetectorType == DetEn_Planar)
				LbInPrintf("\n\tDetector type is planar");
			else
				LbInPrintf("\n\tDetector type is dual-headed planar");
				
			LbInPrintf("\n\tForced interaction is %s",
				((DetRunTimeParams[DetCurParams].DoForcedInteraction == true) ? "on" : "off"));
			LbInPrintf("\n\tNumber of layers in detector is %d",
				DetRunTimeParams[DetCurParams].PlanarDetector.NumLayers);
			for (layerIndex = 0; layerIndex < DetRunTimeParams[DetCurParams].PlanarDetector.NumLayers; layerIndex++){
				LbInPrintf("\n\tLayer %d is '%s', depth = %3.3f, status = %s", layerIndex,
					SubObjGtAttenuationMaterialName(DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[layerIndex].LayerMaterial),
					DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[layerIndex].LayerDepth,
					(DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[layerIndex].IsActive) ? "active" : "not active");
			}
			LbInPrintf("\n\tInner radius of detector is %3.3f",
				DetRunTimeParams[DetCurParams].PlanarDetector.InnerRadius);
			LbInPrintf("\n\tTransaxial length of detector is %3.3f",
				DetRunTimeParams[DetCurParams].PlanarDetector.TransaxialLength);
			LbInPrintf("\n\tAxial length of detector is %3.3f",
				DetRunTimeParams[DetCurParams].PlanarDetector.AxialLength);
			LbInPrintf("\n\tRadial depth of detector is %3.3f",
				detData[DetCurParams].detPlnDetectorDepth);
			if (DetRunTimeParams[DetCurParams].PlanarDetector.NumViews != 0) {
				if (DetRunTimeParams[DetCurParams].PlanarDetector.NumViews > 0) {
					LbInPrintf("\n\tNumber of detector views = %d",
						DetRunTimeParams[DetCurParams].PlanarDetector.NumViews);
				}
			}
			else {
				LbInPrintf("\n\tDetector rotation is continuous");
			}
			
			LbInPrintf("\n\tMinimum detector angle = %3.2f radians",
				DetRunTimeParams[DetCurParams].PlanarDetector.MinAngle);
			LbInPrintf("\n\tMaximum detector angle = %3.2f radians",
				DetRunTimeParams[DetCurParams].PlanarDetector.MaxAngle);
				
			break;
			
		case 	DetEn_Cylindrical:
		
			LbInPrintf("\n\tDetector type is cylindrical");
			LbInPrintf("\n\tForced interaction is %s",
				((DetRunTimeParams[DetCurParams].DoForcedInteraction == true) ? "on" : "off"));
			LbInPrintf("\n\tNumber of rings in detector is %d",
				DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings);
			LbInPrintf("\n\tInner radius of detector is %3.3f",
				detData[DetCurParams].detInBoundCyl.radius);
			LbInPrintf("\n\tOuter radius of detector is %3.3f",
				detData[DetCurParams].detOutBoundCyl.radius);
			LbInPrintf("\n\tAxial length of detector is %3.3f",
				DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings-1].MaxZ-
				DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[0].MinZ);
		
			for (r = 0; r < DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings; r++) {
				LbInPrintf("\n\tRing %d - ", r);
				for (l = 0; l < DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[r].NumLayers; l++) {
					LbInPrintf("\n\t\tlayer %d -", l);
						LbInPrintf("\n\t\t\tdetector material index %d = '%s'"
						"\n\t\t\tinner radius = %3.2f"
						"\n\t\t\touter radius = %3.2f"
						"\n\t\t\tminimum z = %3.2f"
						"\n\t\t\tmaximum z = %3.2f"						
						"\n\t\t\tstatus = %s",
						DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[r].LayerInfo[l].LayerMaterial,
						SubObjGtAttenuationMaterialName(DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[r].LayerInfo[l].LayerMaterial),
						DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[r].LayerInfo[l].InnerRadius,
						DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[r].LayerInfo[l].OuterRadius,
						DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[r].MinZ,
						DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[r].MaxZ,
						(DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[r].LayerInfo[l].IsActive) ? "active" : "not active");					
				}
			}
			break;
				LbInPrintf("\n\tLayer %d is '%s', depth = %3.3f, status = %s", layerIndex,
					SubObjGtAttenuationMaterialName(DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[layerIndex].LayerMaterial),
					DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[layerIndex].LayerDepth,
					(DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo[layerIndex].IsActive) ? "active" : "not active");

		case	DetEn_Block:
			LbInPrintf("\n\tDetector type is block tomograph");
			LbInPrintf("\n\tForced interaction is %s",
				((DetRunTimeParams[DetCurParams].DoForcedInteraction == true) ? "on" : "off"));
			
			/* print out the algorithm for the detected position of photons */
			LbInPrintf("\n\tDetected position algorithm positions photon: %s",
				(( DetRunTimeParams[DetCurParams].BlockTomoDetector.BlockDetectedPositionAlgo == DetEn_BlockSnapCentroidToXtalCenter ) ? 
						"at center of crystal nearest the energy-weighted centroid" : 
						"at energy-weighted centroid"));
					
			LbInPrintf("\n\tNumber of rings in tomograph is %d",
				DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings);
				
			/* iterate through the rings and report the details of the rings, blocks, etc. */
			for (ringIndex = 0; ringIndex < DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings; ringIndex++) {
			
				/*  assign pointer to current ring */
				curRingInfoPtr = &(DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[ringIndex]);
				
				LbInPrintf("\n\tRing %d axial shift is %3.2f",
					ringIndex,
					curRingInfoPtr->AxialShift);
				
				LbInPrintf("\n\tRing %d transaxial rotation is %3.2f",
					ringIndex,
					curRingInfoPtr->TransaxialRotation);
				
				LbInPrintf("\n\tThe geometry of ring %d is defined in the file %s:",
					ringIndex,
					curRingInfoPtr->RingParamFilePath);
				
				if (curRingInfoPtr->XInnerRadius == curRingInfoPtr->YInnerRadius) {
					
					LbInPrintf("\n\t\tRing %d lies outside the circular cylinder of radius %3.2f.",
						ringIndex,
						curRingInfoPtr->XInnerRadius);
					
				} else {
					
					LbInPrintf("\n\t\tRing %d lies outside the elliptical cylinder with x-radius %3.2f",
						ringIndex,
						curRingInfoPtr->XInnerRadius);
					
					LbInPrintf("\n\t\tand y-radius %3.2f.",
						curRingInfoPtr->YInnerRadius);
					
				}
				
				if (curRingInfoPtr->XOuterRadius == curRingInfoPtr->YOuterRadius) {
					
					LbInPrintf("\n\t\tRing %d lies inside the circular cylinder of radius %3.2f.",
						ringIndex,
						curRingInfoPtr->XOuterRadius);
					
				} else {
					
					LbInPrintf("\n\t\tRing %d lies inside the elliptical cylinder with x-radius %3.2f",
						ringIndex,
						curRingInfoPtr->XOuterRadius);
					
					LbInPrintf("\n\t\tand y-radius %3.2f.",
						curRingInfoPtr->YOuterRadius);
					
				}
									
				LbInPrintf("\n\t\tRing %d's minimum axial value is %3.2f.",
					ringIndex,
					curRingInfoPtr->MinZ);
					
				LbInPrintf("\n\t\tRing %d's maximum axial value is %3.2f.",
					ringIndex,
					curRingInfoPtr->MaxZ);
					
				LbInPrintf("\n\t\tRing %d consists of %d blocks.",
					ringIndex,
					curRingInfoPtr->NumBlocks);
				
				if (!detShowBlocks) {
					LbInPrintf("\n\t\tOnly block positions within ring reported below, internal structure of blocks not shown.");
				}
				
				for (blockIndex = 0; blockIndex < curRingInfoPtr->NumBlocks; blockIndex++) {
					
					curBlockInfoPtr = &(curRingInfoPtr->BlockInfo[blockIndex]);
					
					/*
					LbInPrintf("\n\t\tBlock %d from %s at (r,a,z) = (%3.2f, %3.2f, %3.2f), rotated %3.2f degrees.",
							blockIndex,
							curBlockInfoPtr->BlockParamFilePath,
							curBlockInfoPtr->RadialPosition,
							curBlockInfoPtr->AngularPositionDegrees,
							curBlockInfoPtr->ZPosition,
							curBlockInfoPtr->TransaxialOrientation);
					*/
					
					if (detShowBlocks) {
						
						LbInPrintf("\n\t\t\tReference point within block %d for above positioning is at x,y,z coordinates",
								blockIndex);
												
						LbInPrintf("\n\t\t\t\t %3.2f, %3.2f, %3.2f.",
								curBlockInfoPtr->XRef,
								curBlockInfoPtr->YRef,
								curBlockInfoPtr->ZRef);
												
						LbInPrintf("\n\t\t\tThe block x-minimum is %3.2f.", curBlockInfoPtr->XMin);
						LbInPrintf("\n\t\t\tThe block x-maximum is %3.2f.", curBlockInfoPtr->XMax);
												
						LbInPrintf("\n\t\t\tThe block y-minimum is %3.2f.", curBlockInfoPtr->YMin);
						LbInPrintf("\n\t\t\tThe block y-maximum is %3.2f.", curBlockInfoPtr->YMax);
												
						LbInPrintf("\n\t\t\tThe block z-minimum is %3.2f.", curBlockInfoPtr->ZMin);
						LbInPrintf("\n\t\t\tThe block z-maximum is %3.2f.", curBlockInfoPtr->ZMax);
												
						LbInPrintf("\n\t\t\tBlock %d number of layers (in x-direction) is %d.",
								blockIndex,
								curBlockInfoPtr->NumLayers);
						
						for (layerIndex = 0; layerIndex < curBlockInfoPtr->NumLayers; layerIndex++) {
							
							curLayerInfoPtr = &(curBlockInfoPtr->LayerInfo[layerIndex]);
							
							LbInPrintf("\n\t\t\tBlock %d layer %d description:",
									blockIndex,
									layerIndex);
							
							LbInPrintf("\n\t\t\t\tLayer minimum (inner) x boundary is %3.2f.",
									curLayerInfoPtr->InnerX);
							
							LbInPrintf("\n\t\t\t\tLayer maximum (outer) x boundary is %3.2f.",
									curLayerInfoPtr->OuterX);
							
							LbInPrintf("\n\t\t\t\tThere are %d material changes in the y-direction.",
									curLayerInfoPtr->NumYChanges);
							
							LbInPrintf("\n\t\t\t\tY coordinates for changes occur at");
							
							for (yChangeIndex = 0; yChangeIndex < curLayerInfoPtr->NumYChanges; yChangeIndex++) {
								
								LbInPrintf("\n\t\t\t\t\ty = %3.2f",
										curLayerInfoPtr->YChanges[yChangeIndex]);
							
							}
								
							LbInPrintf("\n\t\t\t\tThere are %d material changes in the z-direction.",
									curLayerInfoPtr->NumZChanges);
							
							LbInPrintf("\n\t\t\t\tZ coordinates for changes occur at");
							
							for (zChangeIndex = 0; zChangeIndex < curLayerInfoPtr->NumZChanges; zChangeIndex++) {
								
								LbInPrintf("\n\t\t\t\t\tz = %3.2f",
										curLayerInfoPtr->ZChanges[zChangeIndex]);
								
							}
							
							LbInPrintf("\n\t\t\t\tThe y- and z-material changes create %d material elements:",
									curLayerInfoPtr->NumElements);
							
							for (elementIndex = 0; elementIndex < curLayerInfoPtr->NumElements; elementIndex++) {
							
								LbInPrintf("\n\t\t\t\t\tElement %d consists of %s",
										elementIndex,
										SubObjGtAttenuationMaterialName(curLayerInfoPtr->ElementInfo[elementIndex].MaterialIndex));
								
								if (curLayerInfoPtr->ElementInfo[elementIndex].IsActive) {
									
									LbInPrintf("\n\t\t\t\t\tElement is active (e.g., scintillating).");
									
									LbInPrintf("\n\t\t\t\t\tElement is crystal number %d in block %d, ring %d.",
											curLayerInfoPtr->ElementInfo[elementIndex].crystalNumInBlock,
											blockIndex,
											ringIndex);
											
									LbInPrintf("\n\t\t\t\t\tElement is crystal number %d in tomograph %d, ring %d.",
											curLayerInfoPtr->ElementInfo[elementIndex].crystalNumInTomo,
											blockIndex,
											ringIndex);
											
								} else {
								
									LbInPrintf("\n\t\t\t\t\tElement is not active (e.g., not scintillating).");
									
								}
								
							}
							
							
						}
												
					}
					
				}
					
				
			}

			break;
				
		case	DetEn_simple_pet:
				LbInPrintf("\n\tDetector type is simple PET");
				break;
				
		case	DetEn_simple_spect:
				LbInPrintf("\n\tDetector type is simple SPECT");
				break;
				
		case	DetEn_unc_spect:
				LbInPrintf("\n\tDetector type is UNC SPECT");
				break;
				
			
		default:
			LbInPrintf("\nUnknown detector type in DetPrintParams = %d\n", DetRunTimeParams[DetCurParams].DetectorType);
			break;
	}
	if (DET_DoEnergyBlur() == true) {
		LbInPrintf("\n\tDetector energy is being blurred by a  Gaussian");
		LbInPrintf("\n\tEnergy resolution = %3.1f%%", DetRunTimeParams[DetCurParams].EnergyResolutionPercentage);
		LbInPrintf("\n\tReference energy = %3.1f", DetRunTimeParams[DetCurParams].ReferenceEnergy);
	}
	else {
		LbInPrintf("\n\tDetector is being modelled with perfect energy resolution");
	}
			
	if ( DET_DoTofBlur() == true ) {
		LbInPrintf("\n\tPhoton time-of-detection is being blurred by a Gaussian");
		LbInPrintf("\n\tPhoton time resolution (nanoseconds) = %5.3f",
			DetRunTimeParams[DetCurParams].PhotonTimeFWHM);
		LbInPrintf("\n\tCoincidence time-of-flight resolution (nanoseconds) = %5.3f",
			PHGMATH_SquareRoot( 2.0 * PHGMATH_Square(DetRunTimeParams[DetCurParams].PhotonTimeFWHM) ) );
	}
	else if ( !(PHG_IsSPECT()) ) {
		LbInPrintf("\n\tDetector is being modelled with perfect time-of-flight resolution");
	}
			
	/* Print out history parameters */
	if (detData[DetCurParams].detHistFileHk.doCustom) {
		LbInPrintf("\nHistory file parameters for detector module");
		PhoHFilePrintParams(&(detData[DetCurParams].detHistFileHk));
	}
	}
}

/*********************************************************************************
*
*			Name:			DetPrintReport
*
*			Summary:		Print the final statistics.
*
*			Arguments:
*			Function return: None.
*
*********************************************************************************/
void DetPrintReport(void)
{
	for (DetCurParams = 0; DetCurParams < DetNumParams; DetCurParams++){
	/* Tell the user whats happening */
	LbInPrintf("\n\n***************** Detected Photon Statistics *************\n");
	
	/* Print out our statistics according to type of detector */
	switch(DetRunTimeParams[DetCurParams].DetectorType){
	
		case 	DetEn_Polygonal:				
			break;
		
		
		case 	DetEn_Cylindrical:
			LbInPrintf("\nTotal blue photons entering detector module = %lld.", detData[DetCurParams].detTotBluePhotons);
			LbInPrintf("\nTotal pink photons entering detector module = %lld.", detData[DetCurParams].detTotPinkPhotons);
			LbInPrintf("\n\nTotal blue photons accepted by detector = %lld.", detData[DetCurParams].detTotAccBluePhotons);
			LbInPrintf("\nTotal pink photons accepted by detector = %lld.", detData[DetCurParams].detTotAccPinkPhotons);
			LbInPrintf("\nTotal photons depositing energy in detector = %lld.", detData[DetCurParams].detTotPhotonsDepositingEnergy);
			LbInPrintf("\nTotal photons absorbed in detector = %lld.", detData[DetCurParams].detTotPhotonsAbsorbed);
			LbInPrintf("\nTotal weight absorbed in detector = %3.2e.", detData[DetCurParams].detTotWtAbsorbed);
			LbInPrintf("\nTotal photons forced to absorb in detector = %lld.", detData[DetCurParams].detTotForcedAbsorptions);
			LbInPrintf("\nTotal weight forced to absorb in detector = %3.2e.", detData[DetCurParams].detTotWtForcedAbsorbed);
			LbInPrintf("\nTotal first interaction absorptions = %lld.", detData[DetCurParams].detTotFirstTimeAbsorptions);
			LbInPrintf("\nTotal weight of first interaction absorptions = %3.2e.", detData[DetCurParams].detTotWtFirstTimeAbsorbed);
			LbInPrintf("\nTotal photons passing through crystal without interacting = %lld.", detData[DetCurParams].detTotPhotonsPassingThrough);
		
			#ifdef PHG_DEBUG
			LbInPrintf("\nTotal interactions = %lld", detData[DetCurParams].detCylCountInteractions);
			LbInPrintf("\nTotal coherent interactions = %lld", detData[DetCurParams].detCylCountCohInteractions);
			#endif
			LbInPrintf("\nWeight decremented for FI = %3.2e\n\t", detData[DetCurParams].detWeightAdjusted);
			
			break;
		
		
		case 	DetEn_Block:
			LbInPrintf("\nTotal blue photons entering detector module = %lld.", detData[DetCurParams].detTotBluePhotons);
			LbInPrintf("\nTotal pink photons entering detector module = %lld.", detData[DetCurParams].detTotPinkPhotons);
			LbInPrintf("\nTotal photons depositing energy in detector = %lld.", detData[DetCurParams].detTotPhotonsDepositingEnergy);
			LbInPrintf("\nTotal photons absorbed in detector = %lld.", detData[DetCurParams].detTotPhotonsAbsorbed);
			LbInPrintf("\nTotal weight absorbed in detector = %3.2e.", detData[DetCurParams].detTotWtAbsorbed);
			LbInPrintf("\nTotal photons forced to absorb in detector = %lld.", detData[DetCurParams].detTotForcedAbsorptions);
			LbInPrintf("\nTotal weight forced to absorb in detector = %3.2e.", detData[DetCurParams].detTotWtForcedAbsorbed);
			LbInPrintf("\nTotal first interaction absorptions = %lld.", detData[DetCurParams].detTotFirstTimeAbsorptions);
			LbInPrintf("\nTotal weight of first interaction absorptions = %3.2e.", detData[DetCurParams].detTotWtFirstTimeAbsorbed);
			LbInPrintf("\nTotal photons passing through detectors without interacting = %lld.", detData[DetCurParams].detTotPhotonsPassingThrough);
			LbInPrintf("\nWeight decremented for FI = %3.2e\n\t", detData[DetCurParams].detWeightAdjusted);
			
			break;
		
		
		case	DetEn_DualHeaded:
			LbInPrintf("\nTotal blue photons entering detector module = %lld.", detData[DetCurParams].detTotBluePhotons);
			LbInPrintf("\nTotal pink photons entering detector module = %lld.", detData[DetCurParams].detTotPinkPhotons);
			LbInPrintf("\nTotal photons depositing energy in detector = %lld.", detData[DetCurParams].detTotPhotonsDepositingEnergy);
			LbInPrintf("\nTotal photons absorbed in detector = %lld.", detData[DetCurParams].detTotPhotonsAbsorbed);
			LbInPrintf("\nTotal weight absorbed in detector = %3.2e.", detData[DetCurParams].detTotWtAbsorbed);
			LbInPrintf("\nTotal photons forced to absorb in detector = %lld.", detData[DetCurParams].detTotForcedAbsorptions);
			LbInPrintf("\nTotal weight forced to absorb in detector = %3.2e.", detData[DetCurParams].detTotWtForcedAbsorbed);
			LbInPrintf("\nTotal first interaction absorptions = %lld.", detData[DetCurParams].detTotFirstTimeAbsorptions);
			LbInPrintf("\nTotal weight of first interaction absorptions = %3.2e.", detData[DetCurParams].detTotWtFirstTimeAbsorbed);
			LbInPrintf("\nTotal photons passing through crystal without interacting = %lld.", detData[DetCurParams].detTotPhotonsPassingThrough);
			LbInPrintf("\nWeight decremented for FI = %3.2e\n\t", detData[DetCurParams].detWeightAdjusted);
			
			break;
		
		
		case 	DetEn_Planar:
			LbInPrintf("\nTotal photons entering detector module = %lld.", detData[DetCurParams].detTotBluePhotons);
			LbInPrintf("\nTotal photons reaching the crystal = %lld.", detData[DetCurParams].detTotReachingCrystal);
			LbInPrintf("\nTotal photons depositing energy in detector = %lld.", detData[DetCurParams].detTotPhotonsDepositingEnergy);
			LbInPrintf("\nTotal photons absorbed in detector = %lld.", detData[DetCurParams].detTotPhotonsAbsorbed);
			LbInPrintf("\nTotal weight absorbed in detector = %3.2e", detData[DetCurParams].detTotWtAbsorbed);
			LbInPrintf("\nTotal photons forced to absorb in detector = %lld.", detData[DetCurParams].detTotForcedAbsorptions);
			LbInPrintf("\nTotal weight forced to absorb in detector = %3.2e", detData[DetCurParams].detTotWtForcedAbsorbed);
			LbInPrintf("\nTotal first interaction absorptions = %lld.", detData[DetCurParams].detTotFirstTimeAbsorptions);
			LbInPrintf("\nTotal weight of first interaction absorptions = %3.2e", detData[DetCurParams].detTotWtFirstTimeAbsorbed);
			LbInPrintf("\nTotal photons passing through crystal without interacting = %lld.", detData[DetCurParams].detTotPhotonsPassingThrough);
			LbInPrintf("\nTotal photons reaching maximum number of interactions (%lld) = %lld.", MAX_DET_INTERACTIONS,
				detData[DetCurParams].detNumReachedMaxInteractions);
			LbInPrintf("\nWeight decremented for FI = %3.2e\n\t", detData[DetCurParams].detWeightAdjusted);
			
			#ifdef PHG_DEBUG_YE_IMAGES
				LbInPrintf("\nTotal weight exiting detector through Y limit = %3.3e.", detPlnrExitYWeight);
				LbInPrintf("\nSqrt(w^2) exiting detector through Y limit = %3.3e.", PHGMATH_SquareRoot(detPlnrExitYWeightSq));
				LbInPrintf("\nTotal counts exiting detector through Y limit = %ld.", detPlnrExitYCount);
				LbInPrintf("\nTotal weight exiting detector through Z limit = %3.3e.", detPlnrExitZWeight);
				LbInPrintf("\nSqrt(w^2) exiting detector through Z limit = %3.3e.", PHGMATH_SquareRoot(detPlnrExitZWeightSq));
				LbInPrintf("\nTotal counts exiting detector through Z limit = %ld.", detPlnrExitZCount);
			#endif

			break;
			
			
		default:
			LbInPrintf("\nTotal blue photons reaching detector = %lld.", detData[DetCurParams].detTotBluePhotons);
			LbInPrintf("\nTotal pink photons reaching detector = %lld.", detData[DetCurParams].detTotPinkPhotons);
			LbInPrintf("\n\nTotal blue photons accepted by detector = %lld.", detData[DetCurParams].detTotAccBluePhotons);
			LbInPrintf("\nTotal pink photons accepted by detector = %lld.", detData[DetCurParams].detTotAccPinkPhotons);
			break;
	}
			
	/* Print out history report */
	if (detData[DetCurParams].detHistFileHk.doCustom) {
		LbInPrintf("\nHistory file report for detector module");
		PhoHFilePrintReport(&(detData[DetCurParams].detHistFileHk));
	}
	
	LbInPrintf("\n\n****************************************************************\n");
	}
	fflush(stdout);
}

/*********************************************************************************
*
*			Name:			DetIsSimplePet
*
*			Summary:		Returns true if the detector is simple PET.
*
*			Arguments:
*			Function return: true if detector is simple PET.
*
*********************************************************************************/
Boolean	DetIsSimplePet(void)
{
return((DetRunTimeParams[DetCurParams].DetectorType == DetEn_simple_pet) ? true : false);
}

/*********************************************************************************
*
*			Name:			DetIsSimpleSpect
*
*			Summary:		Returns true if the detector is simple SPECT.
*
*			Arguments:
*			Function return: true if the detector is simple SPECT.
*
*********************************************************************************/
Boolean	DetIsSimpleSpect(void)
{
return((DetRunTimeParams[DetCurParams].DetectorType == DetEn_simple_spect) ? true : false);
}

/*********************************************************************************
*
*			Name:			DetIsUNCSpect
*
*			Summary:		Returns true if the detector is UNC SPECT.
*
*			Arguments:
*			Function return: true if the detector is UNC SPECT.
*
*********************************************************************************/
Boolean	DetIsUNCSpect(void)
{
return((DetRunTimeParams[DetCurParams].DetectorType == DetEn_unc_spect) ? true : false);
}

/*********************************************************************************
*
*			Name:			DetIsPlanar
*
*			Summary:		Returns true if the detector is planar SPECT.
*
*			Arguments:
*			Function return: true if the detector is planar SPECT.
*
*********************************************************************************/
Boolean	DetIsPlanar(void)
{
return((DetRunTimeParams[DetCurParams].DetectorType == DetEn_Planar) ? true : false);
}

/*********************************************************************************
*
*			Name:			DetIsCylindrical
*
*			Summary:		Returns true if the detector is cylindrical.
*
*			Arguments:
*			Function return: true if detector is cylindrical.
*
*********************************************************************************/
Boolean	DetIsCylindrical(void)
{
return((DetRunTimeParams[DetCurParams].DetectorType == DetEn_Cylindrical) ? true : false);
}

/*********************************************************************************
*
*			Name:			DetIsBlock
*
*			Summary:		Returns true if the detector is block.
*
*			Arguments:
*			Function return: true if detector is block.
*
*********************************************************************************/
Boolean	DetIsBlock(void)
{
return((DetRunTimeParams[DetCurParams].DetectorType == DetEn_Block) ? true : false);
}

/*********************************************************************************
*
*			Name:			DetIsDualHeaded
*
*			Summary:		Returns true if the detector is dualed-headed PET.
*
*			Arguments:
*			Function return: true if detector is dual-headed.
*
*********************************************************************************/
Boolean	DetIsDualHeaded(void)
{
return((DetRunTimeParams[DetCurParams].DetectorType == DetEn_DualHeaded) ? true : false);
}

/*********************************************************************************
*
*			Name:			DetGtOutsideRadius
*
*			Summary:		Returns the outside radius of the detector. Originally
*							added for the binning module.
*
*			Arguments:
*			Function return: Outside radius of the detector.
*
*********************************************************************************/
double	DetGtOutsideRadius(void)
{
	double outsideRadius;

	/* If we aren't doing run-time detection just return the collimator radius */
	if (PHG_IsDetectOnTheFly() == false) {
	
			outsideRadius = ColGtOutsideRadius();
	}
	else {
		/* Access the appropriate data structure based on detector type */
		switch(DetRunTimeParams[DetCurParams].DetectorType) {
		
			case DetEn_simple_pet:
			
		
				/* Simple PET detectors don't extend the outside radius
					beyond that of the collimator. So, we'll just return
					the collimator radius.
				*/
				outsideRadius = ColGtOutsideRadius();
				
				break;
		
			case DetEn_simple_spect:
			
		
				/* Simple PET detectors don't extend the outside radius
					beyond that of the collimator. So, we'll just return
					the collimator radius.
				*/
				outsideRadius = ColGtOutsideRadius();
				
				break;
				
			case DetEn_DualHeaded:
				outsideRadius = detData[DetCurParams].detOutBoundCyl.radius;
				break;
				
			case DetEn_Planar:
				outsideRadius = detData[DetCurParams].detOutBoundCyl.radius;
				break;
				
			default:
				PhgAbort("Invalid detector type in DetRunTimeParams (DetGtOutsideRadius)",false);
				break;
				
		}
	}
	return (outsideRadius);
}

/*********************************************************************************
*
*			Name:			DetGtInsideRadius
*
*			Summary:		Returns the inside radius of the detector. Originally
*							added for the binning module.
*
*			Arguments:
*			Function return: Inside radius of the detector.
*
*********************************************************************************/
double	DetGtInsideRadius(void)
{
	double insideRadius;

	/* Verify we have detector information */
	if (PHG_IsDetectOnTheFly() == false) {
		PhgAbort("Your parameters require use of the detector module via detector parameters (DetGtInsideRadius)",false);
	}
	
	/* Access the appropriate data structure based on detector type */
	switch(DetRunTimeParams[DetCurParams].DetectorType) {
	
		case DetEn_simple_pet:
		
	
			/* Simple PET detectors don't extend the outside radius
				beyond that of the collimator. So, we'll just return
				the collimator radius.
			*/
			insideRadius = ColGtOutsideRadius();
			
			break;
	
			/* Simple SPECT detectors don't extend the outside radius
				beyond that of the collimator. So, we'll just return
				the collimator radius.
			*/
		case DetEn_simple_spect:
		
	
			insideRadius = ColGtOutsideRadius();
			
			break;
			
		case DetEn_DualHeaded:
			insideRadius = detData[DetCurParams].detInBoundCyl.radius;
			break;
			
		case DetEn_Planar:
			insideRadius = detData[DetCurParams].detInBoundCyl.radius;
			break;
			
			
		case DetEn_Cylindrical:
			insideRadius = DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[0].LayerInfo[0].InnerRadius;
			break;
			
		case DetEn_Block:
			/* ALERT: this may not be a very useful number for block detectors.
			If you are using it, look at how it is derived in DetInitialize. */
			insideRadius = detData[DetCurParams].detInBoundCyl.radius;
			break;
			
		default:
			PhgAbort("Invalid detector type in DetRunTimeParams (DetGtInsideRadius)",false);
			break;
			
	}
	return (insideRadius);
}

/*********************************************************************************
*
*			Name:			DetPETPhotons
*
*			Summary:		Update the binning images with the current batch of
*							detected photons.
*
*			Arguments:
*				PHG_Decay			decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton *bluePhotons		- The blue photons detected.
*				LbUsFourByte 		numBluePhotons	- The number of blue photons.
*				PHG_TrackingPhoton *pinkPhotons		- The  pink photons detected.
*				LbUsFourByte		numPinkPhotons	- The number of pink photons.
*				DetectedPhotonsTy *detPhotonsPtr	- The accepted detected photons
*			Function return: None.
*
*********************************************************************************/
void DetPETPhotons(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		DetectedPhotonsTy *detPhotonsPtr)
{
	
	/*   Make sure that the number of detected photons gets initialized to 0.  */
	detPhotonsPtr->NumDetectedBluePhotons = 0;
	detPhotonsPtr->NumDetectedPinkPhotons = 0;

	/* For efficiency this routine gets called, but if there are no photons to process,
	then bolt out of here.
	*/
	if ((numBluePhotons == 0) && (numPinkPhotons == 0)) {
		return;
	}
	
	/* Compute statistics */
	detData[DetCurParams].detTotBluePhotons += numBluePhotons;
	detData[DetCurParams].detTotPinkPhotons += numPinkPhotons;

	/* Call appropriate detector routine based on user specified type */
	switch(DetRunTimeParams[DetCurParams].DetectorType) {
	
		case DetEn_simple_pet:
		
			/* Call simple detector module */
			detSimplePET(decayPtr, bluePhotons, numBluePhotons,
				pinkPhotons, numPinkPhotons,
				detPhotonsPtr);
				
			break;
		
		
		case DetEn_DualHeaded:
			/*
			DetDualHeaded(decayPtr, bluePhotons, numBluePhotons,
				pinkPhotons, numPinkPhotons,
				detPhotonsPtr);
			*/
			DetGeometric(DetEn_DualHeaded, 
				decayPtr, bluePhotons, numBluePhotons,
				pinkPhotons, numPinkPhotons,
				detPhotonsPtr);
			break;
		
		
		case DetEn_Polygonal:
			detPolyPET(decayPtr, bluePhotons, numBluePhotons,
				pinkPhotons, numPinkPhotons,
				detPhotonsPtr);
				
			break;
		
		
		case DetEn_Cylindrical:
			/*
			DetCylinder(decayPtr, bluePhotons, numBluePhotons,
				pinkPhotons, numPinkPhotons,
				detPhotonsPtr);
			*/
			DetGeometric(DetEn_Cylindrical, 
				decayPtr, bluePhotons, numBluePhotons,
				pinkPhotons, numPinkPhotons,
				detPhotonsPtr);
			break;
		
		
		case DetEn_Block:
			DetGeometric(DetEn_Block, 
				decayPtr, bluePhotons, numBluePhotons,
				pinkPhotons, numPinkPhotons,
				detPhotonsPtr);
			break;
		
		
		default:
			PhgAbort("Invalid detector type in DetRunTimeParams (DetPETPhotons)", false);
			break;
			
	}
	
	/* Update statistics */
	detData[DetCurParams].detTotAccBluePhotons += detPhotonsPtr->NumDetectedBluePhotons; 
	detData[DetCurParams].detTotAccPinkPhotons += detPhotonsPtr->NumDetectedPinkPhotons; 

	/* Write detected photons to history file, if requested */
	if (DET_IsDoHistory()) {

		/* Write the photons */
		if (PhoHFileWriteDetections(&(detData[DetCurParams].detHistFileHk), decayPtr,
				detData[DetCurParams].detDetectedBluePhotons,
				detData[DetCurParams].detDetectedBluePhotonIndex,
				detData[DetCurParams].detDetectedPinkPhotons,
				detData[DetCurParams].detDetectedPinkPhotonIndex) == false) {
			
			/* Abort Program execution */
			PhgAbort("Got failure from PhoHFileWriteDetections (DetPETPhotons).", true);
		}
	}

}

/*********************************************************************************
*
*			Name:			DetGaussEnergyBlur
*
*			Summary:		Perform a Gaussian blur on the energy.
*
*			Arguments:
*				double		energy - The incoming energy
*
*			Function return: Blurred energy value.
*
*********************************************************************************/
double DetGaussEnergyBlur(double energy)
{
	double				standDev;					/* Standard deviation for gauss blurr */
	double				energyProd;					/* Product of incoming and incident energy */


	/* Compute the standard deviation for the gaussian blurr */
	{
		energyProd = energy *
			DetRunTimeParams[DetCurParams].ReferenceEnergy;
			
		standDev = (DetRunTimeParams[DetCurParams].EnergyResolutionPercentage *
			PHGMATH_SquareRoot(energyProd))/GAUSS_MAGIC_NUM;
	}

	/* Blurr the blue photon */
	return(PhgMathSampleFromGauss(energy, standDev));

}

/*********************************************************************************
*
*			Name:			DetGaussTimeBlur
*
*			Summary:		Perform a Gaussian blur on the photon travel distance
*							(which is used for time-of-flight).
*
*			Arguments:
*				double		travelDistance - The photon travel distance
*
*			Function return: Blurred travel distance value.
*
*********************************************************************************/
double DetGaussTimeBlur(double travelDistance)
{
	double				standDev;					/* Standard deviation for gauss blurr */

	/* Compute the standard deviation for the gaussian blur.
	The FWHM is provided in nanoseconds, which we convert to a standard
	deviation in cm. */
	{
		standDev = ( DetRunTimeParams[DetCurParams].PhotonTimeFWHM *
			1.0E-9 * PHGMATH_SPEED_OF_LIGHT )/(GAUSS_MAGIC_NUM / 100.0);
	}

	/* Blur the blue photon */
	return(PhgMathSampleFromGauss(travelDistance, standDev));

}

/*********************************************************************************
*
*			Name:			detSimplePET
*
*			Summary:		Perform "simple" detection. This routine performs
*							a gaussian blurr on the detected photon's energy value.
*
*			Arguments:
*				PHG_Decay			decayPtr			- The decayPtr that started the process.
*				PHG_TrackingPhoton *bluePhotons		- The blue photons detected.
*				LbUsFourByte 		numBluePhotons	- The number of blue photons.
*				PHG_TrackingPhoton *pinkPhotons		- The  pink photons detected.
*				LbUsFourByte		numPinkPhotons	- The number of pink photons.
*				DetectedPhotonsTy *detPhotonsPtr	- The accepted detected photons
*
*			Function return: None.
*
*********************************************************************************/
void detSimplePET(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		DetectedPhotonsTy *detPhotonsPtr)
{
	LbUsFourByte		photonIndex;				/* Index for current photon */

	/* Clear the counters */
	detPhotonsPtr->NumDetectedBluePhotons = 0;
	detPhotonsPtr->NumDetectedPinkPhotons = 0;
	detData[DetCurParams].detDetectedBluePhotonIndex = 0;
	detData[DetCurParams].detDetectedPinkPhotonIndex = 0;

	do { /* Process Loop */
	
		/* Loop through all blue photons */
		for (photonIndex = 0; photonIndex < numBluePhotons; photonIndex++) {
			
			/* Let user modify and/or reject photons */
			if (DetUsrModPETPhotonsFPtr && 
					(*DetUsrModPETPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
						&(bluePhotons[photonIndex])) == false) {
				
				/* They rejected it so go to next photon */
				continue;
			}
	
			/* Blur energy if it was requested */
			if (DET_DoEnergyBlur() == true) {
			
				bluePhotons[photonIndex].energy =
						DetGaussEnergyBlur(bluePhotons[photonIndex].energy);
								
			}
			
			/* Blur time if it was requested */
			if ( DET_DoTofBlur() == true ) {
			
				bluePhotons[photonIndex].travel_distance =
						DetGaussTimeBlur(bluePhotons[photonIndex].travel_distance);	
									
			}

			/* We made it to here so save the detected photons */
			{
				detPhotonsPtr->DetectedTrkngBluePhotons[detPhotonsPtr->NumDetectedBluePhotons]
					= bluePhotons[photonIndex];

				/* Increment the counters */
				detPhotonsPtr->NumDetectedBluePhotons++;

				/* Update our detected photon block if doing history file */
				if  (DET_IsDoHistory()) {
					DetUpdateDetectedPhotonBlock(&bluePhotons[photonIndex]);
				}
			}
		}
		
		/* Loop through all pink photons */
		for (photonIndex = 0; photonIndex < numPinkPhotons; photonIndex++) {

			/* Let user modify and/or reject photons */
			if (DetUsrModPETPhotonsFPtr && 
					(*DetUsrModPETPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
						&(pinkPhotons[photonIndex])) == false) {
				
				/* They rejected it so go to next pink */
				continue;
			}

			/* Blur energy if it was requested */
			if (DET_DoEnergyBlur() == true) {
			
				pinkPhotons[photonIndex].energy =
						DetGaussEnergyBlur(pinkPhotons[photonIndex].energy);
								
			}
			
			/* Blur time if it was requested */
			if ( DET_DoTofBlur() == true ) {
			
				pinkPhotons[photonIndex].travel_distance =
						DetGaussTimeBlur(pinkPhotons[photonIndex].travel_distance);	
									
			}

			/* We made it to here so save the detected photons */
			{
				detPhotonsPtr->DetectedTrkngPinkPhotons[detPhotonsPtr->NumDetectedPinkPhotons]
					= pinkPhotons[photonIndex];

				/* Increment the counters */
				detPhotonsPtr->NumDetectedPinkPhotons++;

				/* Update our detected photon block if doing history file */
				if  (DET_IsDoHistory()) {
					DetUpdateDetectedPhotonBlock(&pinkPhotons[photonIndex]);
				}
			}
		}
	} while (false);
}



/*********************************************************************************
*
*			Name:			detPolyPET
*
*			Summary:		Perform detection for polygonal detectors
*
*			Arguments:
*				PHG_Decay			decayPtr			- The decayPtr that started the process.
*				PHG_TrackingPhoton *bluePhotons		- The blue photons detected.
*				LbUsFourByte 		numBluePhotons	- The number of blue photons.
*				PHG_TrackingPhoton *pinkPhotons		- The  pink photons detected.
*				LbUsFourByte		numPinkPhotons	- The number of pink photons.
*				DetectedPhotonsTy *detPhotonsPtr	- The accepted detected photons
*
*			Function return: None.
*
*********************************************************************************/
void detPolyPET(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *bluePhotons, LbUsFourByte numBluePhotons,
		PHG_TrackingPhoton *pinkPhotons, LbUsFourByte numPinkPhotons,
		DetectedPhotonsTy *detPhotonsPtr)
{
	LbUsFourByte		photonIndex;				/* Index for current photon */

	/* Clear the counters */
	detPhotonsPtr->NumDetectedBluePhotons = 0;
	detPhotonsPtr->NumDetectedPinkPhotons = 0;
	detData[DetCurParams].detDetectedBluePhotonIndex = 0;
	detData[DetCurParams].detDetectedPinkPhotonIndex = 0;
	
	do { /* Process Loop */
	
		/* Loop through all blue photons */
		for (photonIndex = 0; photonIndex < numBluePhotons; photonIndex++) {
			
			/* Let user modify and/or reject photons */
			if (DetUsrModPETPhotonsFPtr && 
					(*DetUsrModPETPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
						&(bluePhotons[photonIndex])) == false) {
				
				/* They rejected it so go to next photon */
				continue;
			}
	
	
			/* Blur energy if it was requested */
			if (DET_DoEnergyBlur() == true) {
			
				bluePhotons[photonIndex].energy =
						DetGaussEnergyBlur(bluePhotons[photonIndex].energy);
								
			}
			
			/* Blur time if it was requested */
			if ( DET_DoTofBlur() == true ) {
			
				bluePhotons[photonIndex].travel_distance =
						DetGaussTimeBlur(bluePhotons[photonIndex].travel_distance);	
									
			}

			/* We made it to here so save the detected photons */
			{
				detPhotonsPtr->DetectedTrkngBluePhotons[detPhotonsPtr->NumDetectedBluePhotons]
					= bluePhotons[photonIndex];

				/* Increment the counters */
				detPhotonsPtr->NumDetectedBluePhotons++;

				/* Update our detected photon block if doing history file */
				if  (DET_IsDoHistory()) {
					DetUpdateDetectedPhotonBlock(&bluePhotons[photonIndex]);
				}
			}
		}
		
		/* Loop through all pink photons */
		for (photonIndex = 0; photonIndex < numPinkPhotons; photonIndex++) {

			/* Let user modify and/or reject photons */
			if (DetUsrModPETPhotonsFPtr && 
					(*DetUsrModPETPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
						&(pinkPhotons[photonIndex])) == false) {
				
				/* They rejected it so go to next pink */
				continue;
			}

	
			/* Blur energy if it was requested */
			if (DET_DoEnergyBlur() == true) {
			
				pinkPhotons[photonIndex].energy =
						DetGaussEnergyBlur(pinkPhotons[photonIndex].energy);
								
			}
			
			/* Blur time if it was requested */
			if ( DET_DoTofBlur() == true ) {
			
				pinkPhotons[photonIndex].travel_distance =
						DetGaussTimeBlur(pinkPhotons[photonIndex].travel_distance);	
									
			}

			
			/* We made it to here so save the detected photons */
			{
				detPhotonsPtr->DetectedTrkngPinkPhotons[detPhotonsPtr->NumDetectedPinkPhotons]
					= pinkPhotons[photonIndex];

				/* Increment the counters */
				detPhotonsPtr->NumDetectedPinkPhotons++;

				/* Update our detected photon block if doing history file */
				if  (DET_IsDoHistory()) {
					DetUpdateDetectedPhotonBlock(&pinkPhotons[photonIndex]);
				}
			}
		}
	} while (false);
}

/*********************************************************************************
*
*			Name:			DetSPECTPhotons
*
*			Summary:		Update the binning images with the current batch of
*							detected photons.
*
*			Arguments:
*				PHG_Decay			decayPtr			- The decayPtr that started the process.
*				PHG_TrackingPhoton *photons			- The photons detected.
*				LbUsFourByte 		numPhotons		- The number of blue photons.
*				DetectedPhotonsTy *detPhotonsPtr	- The accepted detected photons
*			Function return: None.
*
*********************************************************************************/
void DetSPECTPhotons(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
		DetectedPhotonsTy *detPhotonsPtr)
{

	/*  Make sure that the number of detected photons gets initialized to 0. */
	detPhotonsPtr->NumDetectedBluePhotons = 0;
	
	/* Compute statistics */
	detData[DetCurParams].detTotBluePhotons += numPhotons;

	/* Call appropriate detector routine based on user specified type */
	switch(DetRunTimeParams[DetCurParams].DetectorType) {
		
		case DetEn_simple_spect:
			/* Call simple detector module */
			detSimpleSPECT(decayPtr, photons, numPhotons,
				detPhotonsPtr);
			break;
		
		
		case DetEn_Planar:
			/* Call simple detector module */
			/*
			DetPlanarSPECT(decayPtr, photons, numPhotons,
				detPhotonsPtr);
			*/
			DetGeometric(DetEn_Planar, 
				decayPtr, photons, numPhotons, NULL, 0,
				detPhotonsPtr);
			break;
		
		
		case DetEn_DualHeaded:
			PhgAbort("You must simulate PET, not SPECT, to use dual-headed detectors (DetSPECTPhotons).", false);
		
		
		case DetEn_Cylindrical:
			/* Call cylindrical detector module */
			/*
			detCylindricalSPECT(decayPtr, photons, numPhotons,
				detPhotonsPtr);
			*/
			DetGeometric(DetEn_Cylindrical, 
				decayPtr, photons, numPhotons, NULL, 0,
				detPhotonsPtr);
			break;
		
		
		case DetEn_Block:
			/* Call block detectors module */
			DetGeometric(DetEn_Block, 
				decayPtr, photons, numPhotons, NULL, 0,
				detPhotonsPtr);
			break;
		
		
		default:
			PhgAbort("Invalid detector type in DetRunTimeParams (DetSPECTPhotons)", false);
			break;
		
	}
	
	/* Update statistics */
	detData[DetCurParams].detTotAccBluePhotons += detPhotonsPtr->NumDetectedBluePhotons; 

	/* Write detected photons to history file, if requested */
	if (DET_IsDoHistory()) {

		/* Write the photons */
		if (PhoHFileWriteDetections(&(detData[DetCurParams].detHistFileHk), decayPtr,
				detData[DetCurParams].detDetectedBluePhotons,
				detData[DetCurParams].detDetectedBluePhotonIndex,
				detData[DetCurParams].detDetectedPinkPhotons,
				detData[DetCurParams].detDetectedPinkPhotonIndex) == false) {
			
			/* Abort Program execution */
			PhgAbort("Got failure from PhoHFileWriteDetections trapped in DetSPECTPhotons.",
				true);
		}
	}

}

/*********************************************************************************
*
*			Name:			detSimpleSPECT
*
*			Summary:		Perform "simple" detection. This routine simply
*							performs a gaussian blurr on the detected photon's
*							energy value.
*
*			Arguments:
*				PHG_Decay			decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton *photons			- The blue photons detected.
*				LbUsFourByte 		numPhotons	- The number of blue photons.
*				DetectedPhotonsTy *detPhotonsPtr	- The accepted detected photons
*
*			Function return: None.
*
*********************************************************************************/
void detSimpleSPECT(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
		DetectedPhotonsTy *detPhotonsPtr)
{
	LbUsFourByte		photonIndex;				/* Index for current photon */
	
	
	/* Clear the counters */
	detPhotonsPtr->NumDetectedBluePhotons = 0;
	detPhotonsPtr->NumDetectedPinkPhotons = 0;
	detData[DetCurParams].detDetectedBluePhotonIndex = 0;
	detData[DetCurParams].detDetectedPinkPhotonIndex = 0;

	do { /* Process Loop */
	
		/* Loop through all photons */
		for (photonIndex = 0; photonIndex < numPhotons; photonIndex++) {
			
			/* Let user modify and/or reject photons */
			if (DetUsrModSPECTPhotonsFPtr && 
					(*DetUsrModSPECTPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
						&(photons[photonIndex])) == false) {
				
				/* They rejected it so go to next photon */
				continue;
			}
	
			/* Blurr the blue photon */
			if (DET_DoEnergyBlur()) {
				 photons[photonIndex].energy = 
					DetGaussEnergyBlur( photons[photonIndex].energy );
			}
				
			
			/* We made it to here so save the detected photons */
			{
				detPhotonsPtr->DetectedTrkngBluePhotons[detPhotonsPtr->NumDetectedBluePhotons]
					= photons[photonIndex];

				/* Increment the counters */
				detPhotonsPtr->NumDetectedBluePhotons++;

				/* Update our detected photon block if doing history file */
				if  (DET_IsDoHistory()) {
					DetUpdateDetectedPhotonBlock(&photons[photonIndex]);
				}
			}
		}

	} while (false);
}

/*********************************************************************************
*
*			Name:			detCylindricalSPECT
*
*			Summary:		Perform "cylindrical" detection. This routine exists
*							for the purpose of doing singles type studies
*
*			Arguments:
*				PHG_Decay			decayPtr		- The decayPtr that started the process.
*				PHG_TrackingPhoton *photons			- The blue photons detected.
*				LbUsFourByte 		numPhotons	- The number of blue photons.
*				DetectedPhotonsTy *detPhotonsPtr	- The accepted detected photons
*
*			Function return: None.
*
*********************************************************************************/
void detCylindricalSPECT(PHG_Decay *decayPtr,
		PHG_TrackingPhoton *photons, LbUsFourByte numPhotons,
		DetectedPhotonsTy *detPhotonsPtr)
{
	LbUsFourByte		photonIndex;				/* Index for current photon */

	/* Clear the counters */
	detPhotonsPtr->NumDetectedBluePhotons = 0;
	detPhotonsPtr->NumDetectedPinkPhotons = 0;
	detData[DetCurParams].detDetectedBluePhotonIndex = 0;
	detData[DetCurParams].detDetectedPinkPhotonIndex = 0;
	
	do { /* Process Loop */
	
		/* Loop through all photons */
		for (photonIndex = 0; photonIndex < numPhotons; photonIndex++) {
			
			/* Let user modify and/or reject photons */
			if (DetUsrModSPECTPhotonsFPtr && 
					(*DetUsrModSPECTPhotonsFPtr)(&(DetRunTimeParams[DetCurParams]), decayPtr,
						&(photons[photonIndex])) == false) {
				
				/* They rejected it so go to next photon */
				continue;
			}
	
			/* Track through detector */
			if (DetCylTrack(decayPtr, &(photons[photonIndex])) == false) {
				continue;
			}
			
			/* We made it here so there must be a detection */
			detPhotonsPtr->DetectedTrkngBluePhotons[detPhotonsPtr->NumDetectedBluePhotons]
				= photons[photonIndex];

			/* Increment the counters */
			detPhotonsPtr->NumDetectedBluePhotons++;

			/* Update our detected photon block if doing history file */
			if  (DET_IsDoHistory()) {
				DetUpdateDetectedPhotonBlock(&photons[photonIndex]);
			}
		}

	} while (false);
}

/*********************************************************************************
*
*			Name:			DetGtNumViews
*
*			Summary:		Returns the number of views specified for the planar
*							detector.
*
*			Arguments:
*			Function return: Number of views.
*
*********************************************************************************/
LbFourByte	DetGtNumViews(void)
{
	#ifdef PHG_DEBUG
		/* If we aren't doing run-time detection bolt */
		if (PHG_IsDetectOnTheFly() == false) {
			PhgAbort("You can't call this routine if you aren't modeling detectors! (DetGtNumViews)", false);
		}
	#endif
	

	return (DetRunTimeParams[DetCurParams].PlanarDetector.NumViews);
}

/*********************************************************************************
*
*			Name:			DetTerminate
*
*			Summary:		Process the bin buffers and terminate the module.
*
*			Arguments:
*			Function return: None.
*
*********************************************************************************/

void DetTerminate(void)

{
	Boolean			okay = false;	/* Process Flag */
	
	
	do { /* Process Loop */
		
		/* Call the user termination routine */
		if (DetUsrTerminateFPtr) {
			(*DetUsrTerminateFPtr)(&(DetRunTimeParams[DetCurParams]));
		}
		
		/* Print out some detector statistics if we are not in an error situation */
		if (!ErIsInError()) {
			
			/* Print out our statistics according to type of detector */
			switch(DetRunTimeParams[DetCurParams].DetectorType){
			
				case 	DetEn_Polygonal:				
					break;
					
				case 	DetEn_Cylindrical:
					break;
				
				case 	DetEn_Block:
					break;
				
				case	DetEn_DualHeaded:
					break;
					
				case 	DetEn_Planar:
					break;
					
					
				default:
					break;
			}
	
		}


		/* Close the history file if we created one */
		if (DET_IsDoHistory()) {
			if (!PhoHFileClose(&(detData[DetCurParams].detHistFileHk))) {
				break;
			}
		}

		/* Free memory used by detector */
		switch(DetRunTimeParams[DetCurParams].DetectorType){
		
			case 	DetEn_Polygonal:
				{
					LbUsFourByte ringIndex;
					LbUsFourByte blockIndex;
					
					if (DetRunTimeParams[DetCurParams].PolyDetector.RingInfo != 0) {
						for (ringIndex = 0; ringIndex < DetRunTimeParams[DetCurParams].PolyDetector.NumRings; ringIndex++){
						
							if (DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[ringIndex].BlockInfo != 0){
								for (blockIndex = 0; blockIndex < DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[ringIndex].NumBlocks; blockIndex++){
									
									if (DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[ringIndex].BlockInfo[blockIndex].LayerInfo != 0){
										LbMmFree((void **)&DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[ringIndex].BlockInfo[blockIndex].LayerInfo);
									}
								}
								
								LbMmFree((void **) &DetRunTimeParams[DetCurParams].PolyDetector.RingInfo[ringIndex].BlockInfo);
							}
						}
						LbMmFree((void **) &DetRunTimeParams[DetCurParams].PolyDetector.RingInfo);
					}
				}
				
				break;
			
			case	DetEn_DualHeaded:
			case 	DetEn_Planar:
			
				if (DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo != 0) {
					LbMmFree((void **)&DetRunTimeParams[DetCurParams].PlanarDetector.LayerInfo);
				}
				break;
				
			case 	DetEn_Cylindrical:
				{
					LbUsFourByte ringIndex;
				
					if (DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo != 0) {
						
						for (ringIndex = 0; ringIndex < DetRunTimeParams[DetCurParams].CylindricalDetector.NumRings; ringIndex++){
								
							if (DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[ringIndex].LayerInfo != 0) {
								LbMmFree((void **)&DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo[ringIndex].LayerInfo);
							}
						}
						LbMmFree((void **) &DetRunTimeParams[DetCurParams].CylindricalDetector.RingInfo);
					}
				}
				break;
				
			case 	DetEn_Block:
				{
					LbUsFourByte rNum, bNum, lNum;
				
					if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo != 0) {
						
						for (rNum = 0; rNum < DetRunTimeParams[DetCurParams].BlockTomoDetector.NumRings; rNum++){
								
							if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo != 0) {
								
								for (bNum = 0; bNum < DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].NumBlocks; bNum++){
									
									if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo != 0) {
										
										for (lNum = 0; lNum < DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].NumLayers; lNum++){
										
											if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo[lNum].YChanges != 0) {
												LbMmFree((void **)&DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo[lNum].YChanges);
											}
										
											if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo[lNum].ZChanges != 0) {
												LbMmFree((void **)&DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo[lNum].ZChanges);
											}
										
											if (DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo[lNum].ElementInfo != 0) {
												LbMmFree((void **)&DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo[lNum].ElementInfo);
											}
										
										}
										
										LbMmFree((void **)&DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo[bNum].LayerInfo);
										
									}
									
								}
								
								LbMmFree((void **)&DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo[rNum].BlockInfo);
								
							}
							
						}
						
						LbMmFree((void **) &DetRunTimeParams[DetCurParams].BlockTomoDetector.RingInfo);
						
					}
					
					DetBlocFreeData();
				}
				break;
				
			default:
				break;
		}

			
		okay = true;
	} while (false);
	
	/* Do error handling here */
	if (!okay) {
		ErHandle("Failed to process binning files.", false);
	}
}

/*********************************************************************************
*
*			Name:		DetUpdateDetectedPhotonBlock
*
*			Summary:	Make a new detected photon (To write to the history file).
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*				
*			Function return: None.
*
*********************************************************************************/
void DetUpdateDetectedPhotonBlock(PHG_TrackingPhoton	*trackingPhotonPtr)	
{
	
	/* Store deteced photon in bin */
	if (PHG_IsBlue(trackingPhotonPtr)) {
	
		detData[DetCurParams].detDetectedBluePhotons[detData[DetCurParams].detDetectedBluePhotonIndex] = *trackingPhotonPtr;
		detData[DetCurParams].detDetectedBluePhotonIndex++;
	}
	else {
	
		detData[DetCurParams].detDetectedPinkPhotons[detData[DetCurParams].detDetectedPinkPhotonIndex] = *trackingPhotonPtr;
		detData[DetCurParams].detDetectedPinkPhotonIndex++;
	}

}

#undef DETECTOR
