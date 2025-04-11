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
*			Module Name:		EmisList.c
*			Revision Number:	1.10
*			Date last revised:	23 July 2013
*
*			Programmer:         Steven Vannoy
*			Date Originated:	18 September, 1992
*
*			Module Overview:	This module is responsible for creating decays
*								and controlling the photon tracking.
*
*			References:	        'Emission List Gen Processes' PHG design.
*
**********************************************************************************
*
*			Global functions defined:
*				EmisListCreateDetectedPhoton(PHG_TrackingPhoton	*trackingPhotonPtr);
*				EmisListCreatePhotonList(void);	
*				EmisListDoDetection(PHG_TrackingPhoton	*trackingPhotonPtr);
*				EmisListInitialize(void);
*				EmisListTerminate(void);
*				EmisLisTrackPhoton(PHG_TrackingPhoton *trackingPhotonPtr);
*				EmisListDoComptonInteraction(PHG_TrackingPhoton *trackingPhotonPtr);
*				EmisListDoCoherent(PHG_TrackingPhoton	*trackingPhotonPtr, LbUsFourByte materialIndex);
*
*			Global variables defined:
*				PHG_Decay			EmisListNewDecay	   				The new decay
*				PHG_TrackingPhoton	EmisListBluePhoton					The "current" blue photon
*				PHG_TrackingPhoton	EmisListPinkPhoton					The "current" pink photon
*				PHG_TrackingPhoton	EmisListFDPhoton					Photon used for forced detection
*				LbEightByte			EmisListNumCreated					Number of photons created for tracking
*				char				EmisListIsotopeDataFilePath			Path to isotope data file
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
*			Revision date:		19 December 2012
*
*			Revision description:	Converted timing actions to use LbTiming
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		10 October 2012
*
*			Revision description:	Reduced timing options by system to just two
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
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 July 2004
*
*			Revision description:	
*							- Bug fixes for list mode.
*
*********************************************************************************/

#define EMIS_LIST

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

/* General library include files */
#include "LbTypes.h"
#include "LbError.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbInterface.h"
#include "LbHeader.h"
#include "LbTiming.h"

/* PHG specific include files */
#include "Photon.h"
#include "PhgParams.h"
#include "ColTypes.h"
#include "ColParams.h"
#include "DetTypes.h"
#include "DetParams.h"
#include "CylPos.h"
#include "PhgMath.h"
#include "PhoHFile.h"
#include "PhgHdr.h"
#include "PhoHStat.h"
#include "PhoTrk.h"
#include "ProdTbl.h"
#include "SubObj.h"
#include "Collimator.h"
#include "Detector.h"
#include "EmisList.h"
#include "phg.h"
#include "PhgBin.h"

/* LOCAL CONSTANTS */
#define EMLI_NUM_EMP_RANGE_DISTANCES	1000
/* LOCAL TYPES */

/* LOCAL GLOBALS */
static char					emLiErrStr[1024];					/* For creating error strings */
static Boolean				EmisListIsInitialized = false;			/* Is this module initialized? */
static PHG_TrackingPhoton	*EmisListDetectdTrkngBluePhotons;		/* Blue photons for current decay */
static PHG_TrackingPhoton	*EmisListDetectdTrkngPinkPhotons;		/* Pink photon for current decay */
static PHG_TrackingPhoton	*EmisListTrackedBluePhotons;			/* Blue photons for current decay */
static PHG_TrackingPhoton	*EmisListTrackedPinkPhotons;			/* Pink photon for current decay */
static LbUsFourByte			EmisListCurBluePhotonIndex;				/* Array index for blue photons */
static LbUsFourByte			EmisListCurPinkPhotonIndex;				/* Current pink index */
static CollimatedPhotonsTy	EmisListCollimatedPhotons[PHG_MAX_PARAM_FILES];		/* These are the successfully collimated photons */
static DetectedPhotonsTy	EmisListDetectedPhotons[PHG_MAX_PARAM_FILES];				/* These are the successfully detected photons */
static PhoHFileHkTy			EmisListPHGHistoryFileHk;				/* The PHG History file */
static PHG_Direction		newDecayEmissionAngle;		/* The new decay's emission angle */
#ifdef DO_POS_RANGE_THE_OLD_WAY
static double				emLiEmpiricalRangeDist[EMLI_NUM_EMP_RANGE_DISTANCES];
#endif
static double				emLiIsotopeParams[100];						/* Data for positron range */
static PHG_Direction		emLiBluePolarization;					/* Data for blue polarization */
static PHG_Direction		emLiPinkPolarization;					/* Data for pink polarization */
static double				emLiBluePhiPolar;
static double				emLiPinkPhiPolar;
static Boolean				emLiPolarized;							/* Flag for tracking polarization */

#ifdef PHG_DEBUG
static double				emLiSumSingleCohScatters;			/* Sum of single coherent scatters for validation */
static double				emLiSumSingleCompScatters;			/* Sum of single compton scatters for validation */
static double				emLiSumSingleScatAborptions;		/* Sum of single scatter absorptions for validation */
#ifdef POSRANGE_TEST
static LbUsFourByte			emLiPosRangeDistance[2000];			/* Histogram for positron range distances */
static LbUsFourByte			emLiPosEnergy[2000];				/* Histogram for positron range energies */
#endif
static double				emLisAttWater1000keV;				/* Used for positron range free path calculation */
#define						NUM_MU_DIST_BINS	33
static	LbFourByte			emlsMuDistribution[NUM_MU_DIST_BINS];	/* Histogram of selected mu values for non-collinearity testing */
static	LbFourByte			emLiDBPolarTestBlue[150];
static	LbFourByte			emLiDBPolarTestPink[150];
static	LbFourByte			emLiDBPolarTotalBlue;
static	LbFourByte			emLiDBPolarTotalPink;
static	double				emLiDBPolarPhiMuWt[50][50];
static	double				emLiDBPolarPhiMuWtSq[50][50];
static	LbUsFourByte		emLiDBPolarPhiMuCt[50][50];
#ifdef DO_RH_POS_RANGE
static	double				emLiDbSigma;
#endif
#endif


/* LOCAL MACROS */
/*********************************************************************************
*
*			Name:		EMLIGetMinDeltaT
*
*			Summary:	Compare three delta T values, and return minimum.
*
*			Arguments:
*				dtx		- Delta t with respect to x axis.
*				dty		- Delta t with respect to y axis.
*				dtz		- Delta t with respect to z axis.
*
*			Function return: Minimum of 3 values.
*
*********************************************************************************/
#define	EMLIGetMinDeltaT(dtx, dty, dtz) ((((dtx) <= (dty)) && ((dtx) <= (dtz))) ? (dtx) : \
 (((dty) <= (dtx)) && ((dty) <= (dtz))) ? (dty) : (dtz))


/* PROTOTYPES */
Boolean emLiCreateDecay(void);
Boolean	emLiCreateNewStart(PHG_TrackingPhoton	*trackingPhotonPtr);
void	emLiCreatePhoton(LbUsOneByte flags, PHG_Decay	 *decayPtr,
			PHG_Direction *emissionDirectionPtr,
			LbFourByte sliceIndex, LbFourByte angleIndex,
			LbFourByte xIndex, LbFourByte yIndex);
void	emLiDoEscape(void);
Boolean emLiReadEvent(void);
void 	emLiAdjForNonCollinearity(PHG_TrackingPhoton *trackingPhotonPtr);
void	emLiCreateIsotopeTable(void);	
double	emLiGtIsotopeRange(void);
void	emLiCompPolarization(PHG_TrackingPhoton *bluePhoton, PHG_TrackingPhoton *pinkPhoton);
void	emLiDoPolarizationAdjustCompt(PHG_TrackingPhoton *photonPtr, double mu, double phi,
			PHG_Direction *preDir);
double emLiSamplePositronEnergy(void);
double emLiComputePosRangeWater( double positronEnergy, double *sigmaWater,
			PHG_Direction *positronDirectionPtr );
void emLiPositronTrkRange(	double waterRange,
							PHG_Position startingPos,
							PHG_Direction direction,
							double positronEnergyMeV,
							double sigma,
							PHG_Position *finalPosPtr,
							Boolean *discard,
							LbFourByte *sliceIndex,
							LbFourByte *xIndex,
							LbFourByte *yIndex);

#ifdef __MWERKS__
void EmisListDoMacEvent(void);

/*********************************************************************************
*
*			Name:		EmisListDoMacEvent
*
*			Summary:	Handle Macintosh related events.
*				
*			Function return: None.
*
*********************************************************************************/
void EmisListDoMacEvent()	
{
	(void) SIOUXHandleOneEvent(0);
}
#endif

/*********************************************************************************
*
*			Name:		EmisListDoDetection
*
*			Summary:	Do the things you do when a photon is detected.
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*				
*			Function return: None.
*
*********************************************************************************/
void EmisListDoDetection(PHG_TrackingPhoton	*trackingPhotonPtr)	
{
	LbFourByte		finalAngleIndex;		/* Final angle index */
	LbFourByte		finalSliceIndex;		/* Final sliceIndex */	
	LbFourByte		startIndex;				/* Current start */	
	LbFourByte		startsToProc;			/* Number of starts to process */
	
	/* Save the final indexes */
	finalAngleIndex = trackingPhotonPtr->angleIndex;
	finalSliceIndex = trackingPhotonPtr->sliceIndex;
	
	/* Update detection statistics */
	PhoHStatIncDetectedPhotonStats(trackingPhotonPtr);
	
	/* Get number of starts to process (no more than max starts even though 
		there might have been more than PHG_MAXIMUM_STARTS scatters as reflected
		by trackingPhotonPtr->numStarts)
	*/
	startsToProc = trackingPhotonPtr->numStarts;
	if (startsToProc >= PHG_MAXIMUM_STARTS)
		startsToProc = PHG_MAXIMUM_STARTS -1;
		
	/* Now loop through each start and update its statistics */
	for (startIndex = 0; startIndex < startsToProc; startIndex++) {
		
		/* Set trackingPhoton's slice/angle index */
		trackingPhotonPtr->angleIndex = trackingPhotonPtr->starts_list[startIndex].angleIndex;
		trackingPhotonPtr->sliceIndex = trackingPhotonPtr->starts_list[startIndex].sliceIndex;
		
		/* Update the productivity */
		if (trackingPhotonPtr->num_of_scatters == 0){

			ProdTblAddDetectedProductivity(trackingPhotonPtr,
				PRODTBLFg_Primary);
		}
		else {
			ProdTblAddDetectedProductivity(trackingPhotonPtr,
				PRODTBLFg_Scatter);
		}
	}
	

	/* Restore the final indexes */
	trackingPhotonPtr->angleIndex = finalAngleIndex;
	trackingPhotonPtr->sliceIndex = finalSliceIndex;
	
	EmisListCreateDetectedPhoton(trackingPhotonPtr);
}

/*********************************************************************************
*
*			Name:		EmisListDoCoherent
*
*			Summary:	Do the things you do when a photon interacts in a coherent
*						kind of way. Notably, calculate the new direction.
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*				LbUsFourByte		materialIndex		- The materialthe photon is in.
*				
*			Function return: None.
*
*********************************************************************************/
void EmisListDoCoherent(PHG_TrackingPhoton	*trackingPhotonPtr, LbUsFourByte materialIndex)	
{
	double phi;
	double cosPhi;
	double sinPhi;
	double temp1;
	double mu;
	double magnitude;
	PHG_Direction	origDirection;		/* Incoming direction */
	
	mu = SubObjGetCohTheta2(materialIndex, trackingPhotonPtr->energy);
		
	phi = PhgMathGetRandomNumber() * PHGMATH_2PI;
	
	cosPhi = PHGMATH_Cosine(phi);
	sinPhi = PHGMATH_Sine(phi);
	origDirection = trackingPhotonPtr->angle;

	/* See if we are heading straight out (cosine(z) very close to 1 */
	if (PhgMathRealNumAreEqual(fabs(trackingPhotonPtr->angle.cosine_z), 1.0, -7, 0, 0, 0)) {
	
		/* Calc new direction */
		trackingPhotonPtr->angle.cosine_x = 
			PHGMATH_SquareRoot(1 - PHGMATH_Square(mu))
			* cosPhi;

		trackingPhotonPtr->angle.cosine_y = 
			PHGMATH_SquareRoot(1 - PHGMATH_Square(mu))
			* sinPhi;
		
		/* Note simplification from Kahn's method */
		trackingPhotonPtr->angle.cosine_z = 
			trackingPhotonPtr->angle.cosine_z
			* mu;
	}
	else {
	
		temp1 = (PHGMATH_SquareRoot((1.0-PHGMATH_Square(mu))/(1.0 - PHGMATH_Square(trackingPhotonPtr->angle.cosine_z))));
		
		/* The general case */
		trackingPhotonPtr->angle.cosine_x = (mu * origDirection.cosine_x) +
			(temp1 *
			((origDirection.cosine_x * origDirection.cosine_z * cosPhi) -
			(origDirection.cosine_y * sinPhi)));

		trackingPhotonPtr->angle.cosine_y = (mu *origDirection.cosine_y) +
			(temp1 *
			((origDirection.cosine_y * origDirection.cosine_z * cosPhi) +
			(origDirection.cosine_x * sinPhi)));

		trackingPhotonPtr->angle.cosine_z = (mu *origDirection.cosine_z) -
			(temp1 *
			((1 - PHGMATH_Square(origDirection.cosine_z))) * cosPhi);

	}	
	
	#ifdef HEAVY_PHG_DEBUG
	{
		char	debugStr[1024];
		double	absDiff;
		
		/* Make sure we are close to 1.0 */
		if (PhgMathRealNumAreEqual(PHGMATH_Square(trackingPhotonPtr->angle.cosine_x) +
				PHGMATH_Square(trackingPhotonPtr->angle.cosine_y) +
				PHGMATH_Square(trackingPhotonPtr->angle.cosine_z), 1.0, -7, 0, &absDiff, 0) == false) {
		
			sprintf(debugStr, "\nCalculation of direction cosines not valid. \nSum of squares = %3.2e, 1.0 - Sum of squares = %3.2e",
				PHGMATH_Square(trackingPhotonPtr->angle.cosine_x) +
				PHGMATH_Square(trackingPhotonPtr->angle.cosine_y) +
				PHGMATH_Square(trackingPhotonPtr->angle.cosine_z),
				absDiff);
				
			LbInPrintf("%s", debugStr);
			
		}
	}
	#endif

	/* Re-normalize direction vector to have magnitude 1 */
	magnitude = PHGMATH_SquareRoot(PHGMATH_Square(trackingPhotonPtr->angle.cosine_x) +
		PHGMATH_Square(trackingPhotonPtr->angle.cosine_y) +
		PHGMATH_Square(trackingPhotonPtr->angle.cosine_z));
		
	trackingPhotonPtr->angle.cosine_x = trackingPhotonPtr->angle.cosine_x/magnitude;
	trackingPhotonPtr->angle.cosine_y = trackingPhotonPtr->angle.cosine_y/magnitude;
	trackingPhotonPtr->angle.cosine_z = trackingPhotonPtr->angle.cosine_z/magnitude;


	#ifdef HEAVY_PHG_DEBUG
	{
		char	debugStr[1024];
		double	absDiff;
		
		if (PhgMathRealNumAreEqual((trackingPhotonPtr->angle.cosine_x * origDirection.cosine_x) +
				(trackingPhotonPtr->angle.cosine_y * origDirection.cosine_y) +
				(trackingPhotonPtr->angle.cosine_z * origDirection.cosine_z), mu, -7, 0, &absDiff, 0) == false) {
		
			sprintf(debugStr, "\nInvalid calculation of direction cosines in EmisListDoComptonInteraction. \nx'*x + y'*y + z'*z = %3.2e, theta = %3.2e, dif = %3.2e",
				(trackingPhotonPtr->angle.cosine_x * origDirection.cosine_x) +
				(trackingPhotonPtr->angle.cosine_y * origDirection.cosine_y) +
				(trackingPhotonPtr->angle.cosine_z * origDirection.cosine_z),
				mu, absDiff);
				
			LbInPrintf("%s", debugStr);
		}

	}
	#endif
				    
	/* Update angle index */
	trackingPhotonPtr->angleIndex = ProdTblGetAngleIndex(
		trackingPhotonPtr->sliceIndex,
		trackingPhotonPtr->angle.cosine_z);
	
	#ifdef DEBUGGING_AD
	if (trackingPhotonPtr->num_of_scatters == 0) {
		LbInPrintf("3.3f\t%3.3f\t%3.3f\t%3.3f",
		mu,
		trackingPhotonPtr->angle.cosine_x,trackingPhotonPtr->angle.cosine_y,
		trackingPhotonPtr->angle.cosine_z);
	}
	#endif
	
}

/*********************************************************************************
*
*			Name:		emLiCompPolarization
*
*			Summary:	Compute the polarization vectors
*			Arguments:
*			Function return: None.
*
*********************************************************************************/
void emLiCompPolarization(PHG_TrackingPhoton *bluePhoton, PHG_TrackingPhoton *pinkPhoton)	
{
	double sinPhi, cosPhi, phiPolar, temp;
	PHG_Direction	*blueDir, *pinkDir;
	
	/* Set local shortcuts */
	blueDir = &(bluePhoton->angle);
	pinkDir = &(pinkPhoton->angle);
	
	/* Get a random "phiPolar" */
	phiPolar = PhgMathGetRandomNumber() * PHGMATH_PI;
	
	/* Compute cosine and sine */
	sinPhi = PHGMATH_Sine(phiPolar);
	cosPhi = PHGMATH_Cosine(phiPolar);
	
	/* Check for z cosine close to one and compute special case if so */
	if (PhgMathRealNumAreEqual(fabs(blueDir->cosine_z), 1.0, -7, 0, 0, 0)) {
		emLiBluePolarization.cosine_x = cosPhi;
		emLiBluePolarization.cosine_y = sinPhi;
		emLiBluePolarization.cosine_z = 0.0;
	}
	else {
		/* Compute normal polarization cosine */
		temp = 1/PHGMATH_SquareRoot(1-PHGMATH_Square(blueDir->cosine_z));
		emLiBluePolarization.cosine_x = temp * ((blueDir->cosine_x * blueDir->cosine_z * cosPhi) - (blueDir->cosine_y * sinPhi));
		emLiBluePolarization.cosine_y = temp * ((blueDir->cosine_y * blueDir->cosine_z * cosPhi) + (blueDir->cosine_x * sinPhi));
		emLiBluePolarization.cosine_z = -temp * (1 - PHGMATH_Square(blueDir->cosine_z)) * cosPhi;
	}
	
	/* Compute pink polarization */
	emLiPinkPolarization.cosine_x = (emLiBluePolarization.cosine_y * pinkDir->cosine_z) - (emLiBluePolarization.cosine_z * pinkDir->cosine_y);
	emLiPinkPolarization.cosine_y = (emLiBluePolarization.cosine_z * pinkDir->cosine_x) - (emLiBluePolarization.cosine_x * pinkDir->cosine_z);
	emLiPinkPolarization.cosine_z = (emLiBluePolarization.cosine_x * pinkDir->cosine_y) - (emLiBluePolarization.cosine_y * pinkDir->cosine_x);

	/* Compute pink polar */
	emLiBluePhiPolar = phiPolar;
	emLiPinkPhiPolar = emLiBluePhiPolar + (PHGMATH_PI_DIV2);
	if (emLiPinkPhiPolar >= PHGMATH_PI) {
		emLiPinkPhiPolar -= PHGMATH_PI;
	}
	
	
}

/*********************************************************************************
*
*			Name:		emLiDoPolarizationAdjustCompt
*
*			Summary:	Adjust photon weights due to polarization
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*			Function return: None.
*
*********************************************************************************/
void emLiDoPolarizationAdjustCompt(PHG_TrackingPhoton *photonPtr, double mu, double phi,
		PHG_Direction *preDir)	
{
	double temp;
	double *phiPolar;
	double	cosPhi, sinPhi, eOutNorm, deltaPhi, knTemp, plTemp;
	PHG_Direction *polarDir;
	PHG_Direction *postDir = &(photonPtr->angle);
	
	/* Figure out which polarization we are needing */
	if PHG_IsBlue(photonPtr) {
		phiPolar = &emLiBluePhiPolar;
		polarDir = &emLiBluePolarization;
	}
	else {
		phiPolar = &emLiPinkPhiPolar;
		polarDir = &emLiPinkPolarization;
	}

	/* Compute temporary var */
	temp = PHGMATH_SquareRoot((1-PHGMATH_Square(mu))/(1-PHGMATH_Square(preDir->cosine_z)));
	
	cosPhi = (mu * preDir->cosine_z - postDir->cosine_z)/(temp * (1-PHGMATH_Square(preDir->cosine_z)));
		
	sinPhi = ((postDir->cosine_y - (mu * preDir->cosine_y) - (temp * preDir->cosine_y * preDir->cosine_z * cosPhi)))/
				(temp * preDir->cosine_x);
						
	eOutNorm = photonPtr->energy/511.0;
	
	deltaPhi = phi - (*phiPolar);
	
	knTemp = (eOutNorm + (1/eOutNorm) - (1 - PHGMATH_Square(mu)));
	
	plTemp =  (eOutNorm + (1/eOutNorm) - 2*(1 - PHGMATH_Square(mu))*PHGMATH_Square(PHGMATH_Cosine(deltaPhi)));
	
	photonPtr->photon_scatter_weight *= (plTemp/knTemp);


	#ifdef HEAVY_PHG_DEBUG
	{
		LbUsFourByte muIndex, phiIndex;
		
		muIndex = floor(((mu+1.0)*50)/2);
		if (muIndex > 50) {
			ErAbort("Got bad value for muIndex debugging check EmisListDoComptonInteraction");
		}
		if (muIndex == 50)
			muIndex -= 1;
		
		if (deltaPhi >= PHGMATH_2PI) {
			deltaPhi -= PHGMATH_2PI;
		}
		else if (deltaPhi < 0) {
			deltaPhi += PHGMATH_2PI;
		}
		
		phiIndex = floor(((deltaPhi)*50)/PHGMATH_2PI);
		if (phiIndex > 50) {
			LbInPrintf("\nphiIndex = %d, phi = %3.3f\n", phiIndex, deltaPhi);
			ErAbort("Got bad value for phiIndex debugging check EmisListDoComptonInteraction");
		}
		if (phiIndex == 50)
			phiIndex -= 1;
		
		emLiDBPolarPhiMuWt[muIndex][phiIndex] += photonPtr->photon_scatter_weight;
		emLiDBPolarPhiMuWtSq[muIndex][phiIndex] += PHGMATH_Square(photonPtr->photon_scatter_weight);
		emLiDBPolarPhiMuCt[muIndex][phiIndex] += 1;
	}
	#endif
}

/*********************************************************************************
*
*			Name:		EmisListDoComptonInteraction
*
*			Summary:	Do the things you do when a photon interacts. Notably,
*					Calculate the new direction and energy levels using Kahn's
*					implementation of Klein-Nishina.
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*			Function return: None.
*
*********************************************************************************/
void EmisListDoComptonInteraction(PHG_TrackingPhoton *trackingPhotonPtr)	
{
	Boolean			accepted = false;	/* Loop variable for calculation */
	double			cosPhi, sinPhi;		/* Cosine and sine of phi */
	double			energy;				/* Photon's incoming energy */
	double			magnitude;			/* Magnitude of direction vector */
	double			mu;					/* Cosine of scatter angle */
	double			phi;				/* Scatter angle */
	double			r1, r2, r3;			/* Three random numbers used to calculate new direction and energy */
	double			temp1;				/* Temporary value for calculation */
	double			y;					/* Some intermediate calculation */
	PHG_Direction	origDirection;		/* Incoming direction */

	#ifdef 	HEAVY_PHG_DEBUG
		char 	debugStr[170];
		double	absDiff;		/* Absolute difference */
	#endif
	

	/* If we are doing polarization, compute polarization vectors */
	if (PHG_IsModelPolarization()) {
		if (emLiPolarized == false) {
			emLiCompPolarization(&EmisListBluePhoton, &EmisListPinkPhoton);
			emLiPolarized = true;
		}
	}

	/* Save the original direction */
	origDirection = trackingPhotonPtr->angle;
		
	/* Set local var to energy level */
	energy = (trackingPhotonPtr->energy)/511.0;
	
	/* Begin loop to calculate energy and direction */
	while (!accepted) {
	
		/* Get three random numbers */
		r1 = PhgMathGetRandomNumber();
		r2 = PhgMathGetRandomNumber();
		r3 = PhgMathGetRandomNumber();
				
		/* Make initial test */
		if (r1 <= (((2*energy) + 1)/((2*energy) + 9))) {
			
			/* Do intermediate calculation */
			y = 1 + (2 * energy * r2);
			
			/* Do next test */
			if (r3 <= (4 * ((1/y) - (1/PHGMATH_Square(y))))) {
				
				/* Calculate mu */
				mu = 1 - (2 * r2);
				
				/* Accept these values */
				accepted = true;
			}
		}
		else {
			
			/* Calculate temporary value */
			y = ((2 * energy) + 1)/(1 + (2 * energy * r2));
			
			/* Calculate potential mu */
			mu = 1 - ((1/energy) * (y - 1));
			
			/* Test */
			if (r3 <= (.5 * (PHGMATH_Square(mu) + (1/y)))) {
			
				/* Accept the value */
				accepted = true;
			}
		}
	}

	/* Calculate new energy value */
	trackingPhotonPtr->energy = (trackingPhotonPtr->energy * (1/y));
	
	/* Restrict mu to valid range */
	if (mu < -1.0)
		mu = -1.0;
	else if (mu > 1.0)
		mu = 1.0;
		
	/* Now we must calculate the direction cosines given the cosine of the scatter angle (mu) */
	{
		/* Pick a random angle given the azimuthal direction of scatter */
		phi = PHGMATH_2PI * PhgMathGetRandomNumber();
		cosPhi = PHGMATH_Cosine(phi);
		sinPhi = PHGMATH_Sine(phi);
		
		/* See if we are heading straight out (cosine(z) very close to 1 */
		if (PhgMathRealNumAreEqual(fabs(trackingPhotonPtr->angle.cosine_z), 1.0, -7, 0, 0, 0)) {
		
			/* Calc new direction */
			trackingPhotonPtr->angle.cosine_x = 
				PHGMATH_SquareRoot(1 - PHGMATH_Square(mu))
				* cosPhi;

			trackingPhotonPtr->angle.cosine_y = 
				PHGMATH_SquareRoot(1 - PHGMATH_Square(mu))
				* sinPhi;
			
			/* Note simplification from Kahn's method */
			trackingPhotonPtr->angle.cosine_z = 
				trackingPhotonPtr->angle.cosine_z
				* mu;
		}
		else {
		
			temp1 = (PHGMATH_SquareRoot((1.0-PHGMATH_Square(mu))/(1.0 - PHGMATH_Square(trackingPhotonPtr->angle.cosine_z))));
			
			/* The general case */
			trackingPhotonPtr->angle.cosine_x = (mu * origDirection.cosine_x) +
				(temp1 *
				((origDirection.cosine_x * origDirection.cosine_z * cosPhi) -
				(origDirection.cosine_y * sinPhi)));

			trackingPhotonPtr->angle.cosine_y = (mu *origDirection.cosine_y) +
				(temp1 *
				((origDirection.cosine_y * origDirection.cosine_z * cosPhi) +
				(origDirection.cosine_x * sinPhi)));

			trackingPhotonPtr->angle.cosine_z = (mu *origDirection.cosine_z) -
				(temp1 *
				((1 - PHGMATH_Square(origDirection.cosine_z))) * cosPhi);

		}	
	}
	
	#ifdef HEAVY_PHG_DEBUG
	{
		
		/* Make sure we are close to 1.0 */
		if (PhgMathRealNumAreEqual(PHGMATH_Square(trackingPhotonPtr->angle.cosine_x) +
				PHGMATH_Square(trackingPhotonPtr->angle.cosine_y) +
				PHGMATH_Square(trackingPhotonPtr->angle.cosine_z), 1.0, -7, 0, &absDiff, 0) == false) {
		
			sprintf(debugStr, "\nCalculation of direction cosines not valid. \nSum of squares = %3.2e, 1.0 - Sum of squares = %3.2e",
				PHGMATH_Square(trackingPhotonPtr->angle.cosine_x) +
				PHGMATH_Square(trackingPhotonPtr->angle.cosine_y) +
				PHGMATH_Square(trackingPhotonPtr->angle.cosine_z),
				absDiff);
				
			LbInPrintf("%s", debugStr);
			
		}
	}
	#endif

	/* Re-normalize direction vector to have magnitude 1 */
	magnitude = PHGMATH_SquareRoot(PHGMATH_Square(trackingPhotonPtr->angle.cosine_x) +
		PHGMATH_Square(trackingPhotonPtr->angle.cosine_y) +
		PHGMATH_Square(trackingPhotonPtr->angle.cosine_z));
		
	trackingPhotonPtr->angle.cosine_x = trackingPhotonPtr->angle.cosine_x/magnitude;
	trackingPhotonPtr->angle.cosine_y = trackingPhotonPtr->angle.cosine_y/magnitude;
	trackingPhotonPtr->angle.cosine_z = trackingPhotonPtr->angle.cosine_z/magnitude;
	
	
	#ifdef HEAVY_PHG_DEBUG
	{

		if (PhgMathRealNumAreEqual((trackingPhotonPtr->angle.cosine_x * origDirection.cosine_x) +
				(trackingPhotonPtr->angle.cosine_y * origDirection.cosine_y) +
				(trackingPhotonPtr->angle.cosine_z * origDirection.cosine_z), mu, -7, 0, &absDiff, 0) == false) {
		
			sprintf(debugStr, "\nInvalid calculation of direction cosines in EmisListDoComptonInteraction. \nx'*x + y'*y + z'*z = %3.2e, mu = %3.2e, dif = %3.2e",
				(trackingPhotonPtr->angle.cosine_x * origDirection.cosine_x) +
				(trackingPhotonPtr->angle.cosine_y * origDirection.cosine_y) +
				(trackingPhotonPtr->angle.cosine_z * origDirection.cosine_z),
				mu, absDiff);
				
			LbInPrintf("%s", debugStr);
		}

	}
	#endif

	/* If we are doing polarization, compute polarization vectors */
	if (PHG_IsModelPolarization() && (trackingPhotonPtr->num_of_scatters == 0)) {
			emLiDoPolarizationAdjustCompt(trackingPhotonPtr, mu, phi,&origDirection);
	}
	
	/* Update angle index */
	trackingPhotonPtr->angleIndex = ProdTblGetAngleIndex(
		trackingPhotonPtr->sliceIndex,
		trackingPhotonPtr->angle.cosine_z);

}

/*********************************************************************************
*
*			Name:		EmisListCreateDetectedPhoton
*
*			Summary:	Make a new detected photon (To write to the history file).
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*				
*			Function return: None.
*
*********************************************************************************/
void EmisListCreateDetectedPhoton(PHG_TrackingPhoton	*trackingPhotonPtr)	
{
	
	/* Assign the current photon weight--the weight that will be used hereafter */
	if (trackingPhotonPtr->num_of_scatters > 0) {
		trackingPhotonPtr->photon_current_weight = trackingPhotonPtr->photon_scatter_weight;
	}
	else {
		trackingPhotonPtr->photon_current_weight = trackingPhotonPtr->photon_primary_weight;
	}
		
	/* Store deteced photon in bin */
	if (PHG_IsBlue(trackingPhotonPtr)) {
	
		/* Store the tracked photon */
		EmisListTrackedBluePhotons[EmisListCurBluePhotonIndex] = *trackingPhotonPtr;

		/* Store the photon for use tracking through the tomograph and binning */
		if (PHG_IsCollimateOnTheFly() || PHG_IsDetectOnTheFly() || PHG_IsBinOnTheFly())
			EmisListDetectdTrkngBluePhotons[EmisListCurBluePhotonIndex] = *trackingPhotonPtr;
			
		EmisListCurBluePhotonIndex++;
		
		#ifdef PHG_DEBUG
		{
			char errStr[1024];
			
			if (EmisListCurBluePhotonIndex >= PHG_MAX_DETECTED_PHOTONS) {
				sprintf(errStr, "(EmisListCurBluePhotonIndex = %ld) > (PHG_MAX_DETECTED_PHOTONS = %d)"
					" (EmisListCreateDetectedPhoton)", (unsigned long)EmisListCurBluePhotonIndex, PHG_MAX_DETECTED_PHOTONS);
				PhgAbort(errStr, true);
			}
		}
		#endif
	}
	else {
	
		/* Store the tracked photon */
		EmisListTrackedPinkPhotons[EmisListCurPinkPhotonIndex] = *trackingPhotonPtr;

		/* Store the photon for use tracking through the tomograph and binning */
		if (PHG_IsCollimateOnTheFly() || PHG_IsDetectOnTheFly() || PHG_IsBinOnTheFly())
			EmisListDetectdTrkngPinkPhotons[EmisListCurPinkPhotonIndex] = *trackingPhotonPtr;
				
		EmisListCurPinkPhotonIndex++;
		
		#ifdef PHG_DEBUG
		{
			char errStr[1024];
			
			if (EmisListCurPinkPhotonIndex >= PHG_MAX_DETECTED_PHOTONS) {
				sprintf(errStr, "(EmisListCurPinkPhotonIndex = %ld) > (PHG_MAX_DETECTED_PHOTONS = %d)"
					" (EmisListCreateDetectedPhoton)", (unsigned long)EmisListCurPinkPhotonIndex, PHG_MAX_DETECTED_PHOTONS);
				PhgAbort(errStr, true);
			}
		}
		#endif
	}

}

/*********************************************************************************
*
*			Name:		emLiReadEvent
*
*			Summary:	Reads an event from the history file.
*			Arguments:
*
*				
*			Function return: True unless no photons left.
*
*********************************************************************************/
Boolean emLiReadEvent()	
{
	Boolean				gotOne = false;		/* Flag indicating event was read */
	LbUsOneByte			flag;				/* Storage for type of event to read */
	LbUsFourByte		decaySize;			/* Size of a decay */
	LbUsFourByte		photonSize;			/* Size of a photon */
	PHG_DetectedPhoton	detPhoton;			/* The detected photon */
	PHG_TrackingPhoton	*trPhotonPtr;		/* The tracking photon */
	static FILE 		*historyFile = 0;	/* The history file */
	
	
	do { /* Process Loop */
		
		/* Open the history file if it hasn't been done */
		if (historyFile == 0) {
			if ((historyFile = LbFlFileOpen("the.history", "r+b")) == 0) {
				ErStFileError("Unable to open existing history file.");
				ErAbort("Unable to get events from (emLiReadEvent).");
				break;
			}
			
			/* Seek past the header */
				if (fseek(historyFile, PHG_HDR_HEADER_SIZE, SEEK_SET) != 0){
					ErStFileError("Unable to seek past header (emLiReadEvent).");
					break;
				}
		}
		
		/* Read in the event flag and verify that we start with a decay */
		if ((fread(&flag, sizeof(LbUsOneByte), 1, historyFile)) != 1) {
		
			/* See if we are not at end of file */
			if (feof(historyFile) == 0) {
				PhgAbort("Unable to read event type.", false);
			}
			else {
				/* We are end of file so just break */
				break;
			}
		}
	
		/* In order to be in synch, this must be a decay */
		if (PHG_IsADecay((LbUsFourByte) flag)) {
						
			/* Read the decay */
			decaySize = sizeof(PHG_Decay);
			
			if ((fread(&EmisListNewDecay, decaySize, 1, historyFile)) != 1) {
				ErAbort("Unable to read decay.");
			}
		}
		else {
			PhgAbort("Did not get decay flag when expected from history file (emLiReadEvent)", false);
		}
		
		/* Read in the event flag and verify that we have a photon */
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

		 if (PHG_IsAPhoton((LbUsFourByte) flag)){
			
			/* Read the photon */
			photonSize = sizeof(PHG_DetectedPhoton);
			
			if ((fread(&detPhoton, photonSize, 1, historyFile)) != 1) {
				ErAbort("Unable to read photon.");
			}
		}
		else {
			PhgAbort("Did not get photon flag when expected from history file (emLiReadEvent)", false);
		}
		
		/* Convert to a tracking photon */
		{
			if (LbFgIsSet(detPhoton.flags, PHGFg_PhotonBlue) == true)
				trPhotonPtr = &EmisListBluePhoton;
			else
				trPhotonPtr = &EmisListPinkPhoton;
				
			/* Flags are currently stored in first two bytes, with upper bits
				representing number of scatters
			*/
			trPhotonPtr->flags = (detPhoton.flags & 3);
			trPhotonPtr->num_of_scatters = 
				(detPhoton.flags >> 2);
	
			trPhotonPtr->location.x_position = detPhoton.location.x_position;
			trPhotonPtr->location.y_position = detPhoton.location.y_position;
			trPhotonPtr->location.z_position = detPhoton.location.z_position;
			trPhotonPtr->angle.cosine_x = detPhoton.angle.cosine_x;
			trPhotonPtr->angle.cosine_y = detPhoton.angle.cosine_y;
			trPhotonPtr->angle.cosine_z = detPhoton.angle.cosine_z;
			trPhotonPtr->transaxialPosition = detPhoton.transaxialPosition;
			trPhotonPtr->azimuthalAngleIndex = detPhoton.azimuthalAngleIndex;
			
			SubObjGtPositionIndexes(&trPhotonPtr->location,
			&trPhotonPtr->sliceIndex, &trPhotonPtr->xIndex,
			&trPhotonPtr->yIndex);
			
			trPhotonPtr->angleIndex =  -1;
			trPhotonPtr->origSliceIndex = -1;
			trPhotonPtr->origAngleIndex = -1;
			trPhotonPtr->scatters_in_col = 0; /* We don't know this value */
			
			if (trPhotonPtr->num_of_scatters == 0){
				trPhotonPtr->photon_scatter_weight = 0;
				trPhotonPtr->photon_primary_weight =
					detPhoton.photon_weight;
				trPhotonPtr->photon_current_weight =
					detPhoton.photon_weight;
			}
			else {
				trPhotonPtr->photon_scatter_weight =
					detPhoton.photon_weight;
				trPhotonPtr->photon_current_weight =
					detPhoton.photon_weight;
				trPhotonPtr->photon_primary_weight = 0;
			}
			trPhotonPtr->scatter_target_weight  = 0;
			trPhotonPtr->decay_weight = EmisListNewDecay.startWeight;
			trPhotonPtr->energy = detPhoton.energy;
			trPhotonPtr->travel_distance  =
				detPhoton.time_since_creation * PHGMATH_SPEED_OF_LIGHT;
			trPhotonPtr->numStarts = 0;
			trPhotonPtr->number = (LbUsEightByte) -1;
		}
		
		gotOne = true;
	} while (false);
	
	return(gotOne);
}

/*********************************************************************************
*
*			Name:		emLiCreateDecay
*
*			Summary:	Create a new decay according to voxel decay information.
*			Arguments:
*
*				
*			Function return: True unless no photons left.
*
*********************************************************************************/
Boolean emLiCreateDecay()	
{
	Boolean				gotOne = false;		      	/* Did we get a photon */
	Boolean				noMoreDecays;				/* Did we reach the end of the decays? */
	LbFourByte			sliceIndex;		      		/* Slice decay came from */
	LbFourByte   		angleIndex;	    	   		/* Stratification angle slice came from */
	LbFourByte			xIndex;		      			/* X index into object */
	LbFourByte   		yIndex;	    	   			/* Y index into object */
	PHG_Position		newPosition;
	Boolean				discard;
	double				positronEnergy;			/* the emission energy of the positron */
	double				positronEnergyMeV;		/* positron energy in MeV */
	double				rangeInWater;			/* sampled positron range for water */
	double				sigmaWater;				/* standard dev of range in water */
	PHG_Direction		decayTravelAngle;
	LbFourByte			sliceIndexOr;		      		/* Slice decay came from */
	LbFourByte			xIndexOr;		      			/* X index into object */
	LbFourByte   		yIndexOr;	    	   			/* Y index into object */
	
	
	do { /* Process Loop */
	
		do {		
			/* Loop to get good decay - to eliminate decays where positron range
				takes the decay outside the object */
			noMoreDecays = false;
			#ifdef PHG_DEBUG
				/* Save the current slice/voxel/angle information */
				EmisListRestartSlice = SubObjCurSliceIndex;
				EmisListRestartVoxel = SubObjCurVoxelIndex;
				EmisListRestartAngle = SubObjCurAngleIndex;
				
				/* If requested write the current random seed info to disk to recreate this photon */
				if (PHGDEBUG_WriteToRandSeedFile()) {
					PhgMathWriteSeed(-1);
				}
			#endif

			/* Get next available decay */
			if (SubObjGenVoxAngCellDecay(&EmisListNewDecay, &newDecayEmissionAngle,
			        &sliceIndex, &angleIndex, &xIndex, &yIndex) == false) {

				/* We are out of decays */
				noMoreDecays = true;	/* Later used to break out of the process loop */
				break;
			}
			
			discard = false;	/* Get out of the decay creation loop when positron range is not in use */

			/* See if we need to adjust for positron range */
			if (PHG_IsPET() && PHG_IsRangeAdjust()){
				
				/* Save our starting indexes */
				xIndexOr = xIndex;
				yIndexOr = yIndex;
				sliceIndexOr = sliceIndex;

				/* Compute emission energy of positron  */
				positronEnergy = emLiSamplePositronEnergy();		

/* !RH debug tests */
#ifdef PHG_DEBUG
	#ifdef POSRANGE_TEST_B1or2
	{
		double randE;
		
		randE = PhgMathGetRandomNumber();
		if ( randE < 0.4 ) {
			positronEnergy = 200.0;
		} else {
			positronEnergy = 1000.0;
		}
	}
	#endif
	#ifdef POSRANGE_TEST_B3
	{
			positronEnergy = 1000.0;
	}
	#endif
#endif


#ifdef PHG_DEBUG
	#ifdef POSRANGE_TEST
	{
		/*	!RH debug binning for positron energy testing.
		*	Use the Test A2 binning as the default, as it covers a larger range. */
		LbUsFourByte	index;
		
		#ifdef POSRANGE_TEST_A1
		{
			if (positronEnergy < 1000.0) {
				index = (int) (positronEnergy);
				emLiPosEnergy[index] += 1;
			}
		}
		#else
		{
			if (positronEnergy < 4000.0) {
				index = (int) (positronEnergy/2);
				emLiPosEnergy[index] += 1;
			}
		}
		#endif
	}
	#endif
#endif
				
				/* Convert energy to MeV units */
				positronEnergyMeV = positronEnergy / 1000;
				

				/* Compute range and direction for positron to travel */
				rangeInWater = emLiComputePosRangeWater( positronEnergyMeV, &sigmaWater, &decayTravelAngle);
				
				/* Adjust the position based on range and direction of positron */
				emLiPositronTrkRange(rangeInWater,
							EmisListNewDecay.location,
							decayTravelAngle,
							positronEnergyMeV,
							sigmaWater,
							&newPosition,
							&discard,
							&sliceIndex,
							&xIndex,
							&yIndex);
			
				
#ifdef PHG_DEBUG
	#ifdef POSRANGE_TEST
	{		/* !RH debug binning for positron range test */
	double			distance;
	LbUsFourByte	index;

		if (discard == false) {
			distance = PHGMATH_SquareRoot(PHGMATH_Square(newPosition.x_position)+
				PHGMATH_Square(newPosition.y_position) + PHGMATH_Square(newPosition.z_position));

			if (distance < 4.0) {
				#ifdef POSRANGE_TEST_B3
				{
					if ( fabs(decayTravelAngle.cosine_z) > 0.995 ) {
						index = (distance * 1000)/2.0;
						if (index == 2000)
							index--;
						emLiPosRangeDistance[index] += 1;
					}
				}
				#else
				{
					index = (distance * 1000)/2.0;
					if (index == 2000)
						index--;
					emLiPosRangeDistance[index] += 1;
				}
				#endif
			}
		}
	}
	#endif
#endif

				/* We made it here so update the position */
				EmisListNewDecay.location = newPosition;

				
			}
			
		} while (discard == true);		/* end of decay creation loop */
		
		if (noMoreDecays == true) {
			/* This breaks out of the process loop when SubObjGenVoxAngCellDecay runs out of decays */
			/* We are out of decays */
			break;
		}
		
		
		/* Make a new photon from the generated decay information */
		emLiCreatePhoton((PHGFg_PhotonBlue), &EmisListNewDecay,
			&newDecayEmissionAngle, sliceIndex,
			angleIndex, xIndex, yIndex);


		/* Remember that we created a photon */
		gotOne = true;
		
		/* See if we need to create a coincident photon */
		if (PHG_IsPET()) {
			
			/* Create a second photon in the opposite direction */
			angleIndex = (PRODTBLGetNumAngleCells() - 1) - angleIndex;
			
			emLiCreatePhoton(0, &EmisListNewDecay,
				&newDecayEmissionAngle, sliceIndex,
				angleIndex, xIndex, yIndex);
				
			/* Adjust for non-collinearity if requested	*/
			if (PHG_IsNonCollinearityAdjust()) {
				if (PhgMathGetRandomNumber() < 0.5) {
					emLiAdjForNonCollinearity(&EmisListPinkPhoton);
				}
				else {
					emLiAdjForNonCollinearity(&EmisListBluePhoton);
				}
			}
		}
			   	
	} while (false);
	
	/* Set polarization flag if we are modelling polarization */
	if (gotOne && PHG_IsModelPolarization()) {
		emLiPolarized = false;
	}
	
	return (gotOne);
}

/*********************************************************************************
*
*			Name:		emLiAdjForNonCollinearity
*
*			Summary:	Adjust the tracking_photon's direction to account for 
.						non-collinearity in positron annhilation.
*			Arguments:
*				PHG_TrackingPhoton	*tp	- The tracking photon.
*
*				
*			Function return: None.
*
*********************************************************************************/
void emLiAdjForNonCollinearity(PHG_TrackingPhoton *tp)	
{
	double 			phi;
	double 			cosPhi;
	double			sinPhi;
	double 			mu;
	double 			temp1;
	double			theta;
	PHG_Direction	origDirection;
	
	/* Determine deviation from collinearity  (mean, standard.deviation) */
	theta = PhgMathSampleFromGauss(0.0, 0.0037059);
	mu = PHGMATH_Cosine(theta);

#ifdef PHG_DEBUG
/* Bin the mu value */
{
LbFourByte i =  (LbFourByte) floor((
							(mu + 1.0) *
							 NUM_MU_DIST_BINS) / 2);
emlsMuDistribution[i] += 1;
}
#endif



	/* Save current direction vector */
	origDirection = tp->angle;
	
	/* Now we must calculate the direction cosines given the cosine of the non-collinearity */
	{
		/* Pick a random angle given the azimuthal direction of non-collinearity */
		phi = PHGMATH_2PI * PhgMathGetRandomNumber();
		cosPhi = PHGMATH_Cosine(phi);
		sinPhi = PHGMATH_Sine(phi);
		
		/* See if we are heading straight out cosine(z) very close to 1 */
		if (PhgMathRealNumAreEqual(fabs(tp->angle.cosine_z), 1.0, -7, 0, 0, 0)) {
		
			/* Calc new direction */
			tp->angle.cosine_x = 
				PHGMATH_SquareRoot(1 - PHGMATH_Square(mu))
				* cosPhi;

			tp->angle.cosine_y = 
				PHGMATH_SquareRoot(1 - PHGMATH_Square(mu))
				* sinPhi;
			
			/* Note simplification from Kahn's method */
			tp->angle.cosine_z = 
				tp->angle.cosine_z
					* mu;
		}
		else {
		
			temp1 = (PHGMATH_SquareRoot((1.0-PHGMATH_Square(mu))/(1.0 - PHGMATH_Square(tp->angle.cosine_z))));
			
			/* The general case */
			tp->angle.cosine_x = (mu * origDirection.cosine_x) +
				(temp1 *
				((origDirection.cosine_x * origDirection.cosine_z * cosPhi) -
				(origDirection.cosine_y * sinPhi)));

			tp->angle.cosine_y = (mu *origDirection.cosine_y) +
				(temp1 *
				((origDirection.cosine_y * origDirection.cosine_z * cosPhi) +
				(origDirection.cosine_x * sinPhi)));

			tp->angle.cosine_z = (mu *origDirection.cosine_z) -
				(temp1 *
				((1 - PHGMATH_Square(origDirection.cosine_z))) * cosPhi);

		}
	}
	
	#ifdef PHG_DEBUG
{	
		/* Validate that the computation was valid */
		double dotProd;
		
		dotProd  = (tp->angle.cosine_x * origDirection.cosine_x) +
			(tp->angle.cosine_y * origDirection.cosine_y) +
			(tp->angle.cosine_z * origDirection.cosine_z);
			
		if (PhgMathRealNumAreEqual(mu, dotProd, -4, 0, 0, 0) == false) {
			sprintf(emLiErrStr, "Invalid calculation of 'mu' according to check, \n\t mu = %3.3f and dotProd = %3.3f\n (emLiAdjForNonCollinearity)", mu, dotProd);
			ErAbort(emLiErrStr);
		}
}
	#endif
	 	
}

/*********************************************************************************
*
*			Name:		emLiCreateNewStart
*
*			Summary:	Create a new start for the photon.
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*
*				
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean emLiCreateNewStart(PHG_TrackingPhoton *trackingPhotonPtr)	
{
	Boolean				okay = false;			/* Process Loop */
	
	do { /* Process Loop */
	
		/* Only do this if we haven't reached the maximum number of starts */
		if (trackingPhotonPtr->numStarts < PHG_MAXIMUM_STARTS) {

			/* Initialize the start information */
			trackingPhotonPtr->starts_list[trackingPhotonPtr->numStarts].angleIndex = 
				trackingPhotonPtr->angleIndex;
				
			trackingPhotonPtr->starts_list[trackingPhotonPtr->numStarts].sliceIndex =
				trackingPhotonPtr->sliceIndex;
			
			/*	The following check will never get hit due to the initial check above. It
				used to be here while we were worried about the tracking algorithm, but now
				we are not. It's left here to be reinstated if necessary.
			*/
			/* Verify we haven't gone over the limit */
			if (trackingPhotonPtr->numStarts == PHG_MAXIMUM_STARTS) {
	
				#ifdef PHG_DEBUG
					LbInPrintf("\nPhoton has exceeded maximum number of starts, photon #%lld\n",
						trackingPhotonPtr->number);
				#endif
				break;
			}
		}
		
		/* Increment the number of starts (This continues even though tracking starts does not) */
		trackingPhotonPtr->numStarts++;
			

		okay = true;
	} while (false);
	
	return (okay);
}

/*********************************************************************************
*
*			Name:		emLiCreatePhoton
*
*			Summary:	Make a new tracking photon.
*			Arguments:
*				LbUsOneByte			flags					- Modifying flags.
*				PHG_Decay			*decayPtr				- The decay.
*				PHG_Direction		*emissionDirectionPtr	- Direction of emission.
*				LbFourByte			sliceIndex				- Slice decay came from.
*				LbFourByte			angleIndex				- Angle decay came from.
*				LbFourByte			xIndex					- X index into object.
*				LbFourByte			yIndex					- Y index into object.
*				
*			Function return: None.
*
*********************************************************************************/
void emLiCreatePhoton(LbUsOneByte flags, PHG_Decay	 *decayPtr,
		PHG_Direction *emissionDirectionPtr,
		LbFourByte sliceIndex, LbFourByte angleIndex,
		LbFourByte xIndex, LbFourByte yIndex)	
{
	PHG_TrackingPhoton	*newPhotonPtr;		/* The "new" photon */
	
	EmisListNumCreated++;
	
	/* Set appropriate global photon variable (based on blue/pink flag).
		NOTE, if SPECT then Blue is always set.
	*/
	if (LbFgIsSet(flags, PHGFg_PhotonBlue) == true)
		newPhotonPtr = &EmisListBluePhoton;
	else
		newPhotonPtr = &EmisListPinkPhoton;

	
	/* Initialize fields */
	newPhotonPtr->flags = flags;
	newPhotonPtr->location = decayPtr->location;
	newPhotonPtr->num_of_scatters = 0;
	newPhotonPtr->scatters_in_col = 0;
	newPhotonPtr->photon_scatter_weight = 1.0;
	newPhotonPtr->photon_primary_weight = 1.0;
	newPhotonPtr->photon_current_weight = 0.0;		/* Unused until photon leaves object */
	newPhotonPtr->decay_weight = decayPtr->startWeight;
	newPhotonPtr->scatter_target_weight = 0;
	newPhotonPtr->energy = PhgRunTimeParams.PhgNuclide.photonEnergy_KEV;
	newPhotonPtr->travel_distance = 0;
	newPhotonPtr->angleIndex = angleIndex;
	newPhotonPtr->sliceIndex = sliceIndex;
	newPhotonPtr->origAngleIndex = angleIndex;
	newPhotonPtr->origSliceIndex = sliceIndex;
	newPhotonPtr->xIndex = xIndex;
	newPhotonPtr->yIndex = yIndex;
	newPhotonPtr->angle = *emissionDirectionPtr;
	newPhotonPtr->num_det_interactions = 0;

	/* The following fields are used for SPECT collimation */
	{
		newPhotonPtr->transaxialPosition = -1;
		newPhotonPtr->axialPosition = -1;
		newPhotonPtr->azimuthalAngleIndex = -1;
		newPhotonPtr->detectorAngle = -1;
	}
	
	newPhotonPtr->numStarts = 0;
	#ifdef PHG_DEBUG
		/* Clear the starts_list array so the debugger doesn't show the whole thing */
		memset((void *)newPhotonPtr->starts_list,
			(PHG_MAXIMUM_STARTS * sizeof(PHG_InteractionInfo)), 0);
	#endif
	
	/* Now if not a blue photon, flip the direction */
	if (!PHG_IsBlue(newPhotonPtr)) {

		/* Reverse direction, this photon should go the opposite direction */
		newPhotonPtr->angle.cosine_x *= -1.0;
		newPhotonPtr->angle.cosine_y *= -1.0;
		newPhotonPtr->angle.cosine_z *= -1.0;
	}

	#ifdef PHG_DEBUG
		/* Set tracking photon's number for debugging */
		newPhotonPtr->number = EmisListNumCreated;
	#endif		

	/* Add initial position to the starts_list */
	if (!emLiCreateNewStart(newPhotonPtr)) {
		PhgAbort("Unable to add original start (emLiCreatePhoton).", true);
	}
	
	/* See if scatter productivity is higher */
	if (PRODTBLGetScatCellProductivity(sliceIndex, angleIndex) >
		PRODTBLGetPrimCellProductivity(sliceIndex, angleIndex)) {
		
		/* Mark this photon to be tracked as a scatter photon  */
		PHGSetTrackAsScatter(newPhotonPtr);
		
		/* Set scatter target */
		newPhotonPtr->scatter_target_weight = 
			PRODTBLGetScatCellProductivity(sliceIndex, angleIndex);

		
		/* Record the photon's start in the productivity table
		ProdTblAddStartingProductivity(
			newPhotonPtr, PRODTBLFg_Scatter); */
		
		/* See if we should also track as primary */
		if ((PhgMathGetRandomNumber() * 
			PRODTBLGetScatCellProductivity(sliceIndex, angleIndex)) <=
			PRODTBLGetPrimCellProductivity(sliceIndex, angleIndex)) {
			
			/* Mark this photon to be tracked as a primary photon also */
			PHGSetTrackAsPrimary(newPhotonPtr);
			
			/* Adjust photon weight */
			newPhotonPtr->photon_primary_weight *=
				(PRODTBLGetScatCellProductivity(sliceIndex, angleIndex)/
				PRODTBLGetPrimCellProductivity(sliceIndex, angleIndex));
			
			
			/* Record the photon's start in the productivity table
			ProdTblAddStartingProductivity(newPhotonPtr, 
			    PRODTBLFg_Primary); */
		}
	}
	else {
		/* Mark this photon to be tracked as a primary photon */
		PHGSetTrackAsPrimary(newPhotonPtr);
		
		/* Record the photon's start in the productivity table
		ProdTblAddStartingProductivity(newPhotonPtr, PRODTBLFg_Primary); */
		
		/* See if we should also track as scatter */
		if ((PhgMathGetRandomNumber() * 
			 PRODTBLGetPrimCellProductivity(sliceIndex, angleIndex)) <=
			PRODTBLGetScatCellProductivity(sliceIndex, angleIndex)) {
			
			/* Mark this photon to be tracked as scatter photon also */	
			PHGSetTrackAsScatter(newPhotonPtr);
			
			/* Adjust photon weight */
			newPhotonPtr->photon_scatter_weight *=
			(PRODTBLGetPrimCellProductivity(sliceIndex, angleIndex)/
			PRODTBLGetScatCellProductivity(sliceIndex, angleIndex));
			
			/* Set scatter target */
			newPhotonPtr->scatter_target_weight = 
				PRODTBLGetPrimCellProductivity(sliceIndex, angleIndex);
			
			/* Record the photon's start in the productivity table
			ProdTblAddStartingProductivity(newPhotonPtr, PRODTBLFg_Scatter); */
		}
	}
}


/*********************************************************************************
*
*			Name:		EmisListCreatePhotonList
*
*			Summary:	Create the photon list.
*			Arguments:
*				
*			Function return: true unless an error occurs.
*
*********************************************************************************/
Boolean EmisListCreatePhotonList()	
{
	Boolean 			okay = false;				/* Process flag */
	LbUsFourByte		curBinParams;				/* LCV For bin parameters */
	Boolean 			histFileCreated = false;	/* Did we get history file created */
	LbTmTimingType		calcDecaysStartTime;		/* Time we started calculating decays */
	LbTmTimingType		trackPhotonsStartTime;		/* Time we started tracking photons */
	double				calcDecaysTime;				/* User time to calc decays */
	double				calcDecaysCPUTime;			/* System time to calc decays */
	double				trackPhotonsTime;			/* User time to track photons */
	double				trackPhotonsCPUTime;		/* System time to track photons */
	Boolean				timingValid;				/* Whether timing exists on this system */
	LbUsFourByte		loopV;						/* Loop counter */
	
	
	do { /* Process Loop */

	
		/* Initialize the photon history list if requested */
		if (PHG_IsHist()) {
			if ((histFileCreated = PhoHFileCreate(PhgRunTimeParams.PhgPhoHFileHistoryFilePath,
					PhgRunTimeParams.PhgPhoHParamsFilePath,
					PhoHFileEn_PHG, &EmisListPHGHistoryFileHk)) == false) {
					
				sprintf(emLiErrStr,"Failed to create history file specified in PHG parameters file named:\n"
					"'%s'\n"
					"The custom parameters file name is: '%s'\n"
					" (EmisListCreatePhotonList)",
					PhgRunTimeParams.PhgPhoHFileHistoryFilePath, PhgRunTimeParams.PhgPhoHParamsFilePath);
				PhgAbort(emLiErrStr, false);
			}
					
			/* Print out history parameters */
			if (EmisListPHGHistoryFileHk.doCustom) {
				LbInPrintf("\n\n\tHistory file parameters for PHG history file");
				PhoHFilePrintParams(&EmisListPHGHistoryFileHk);
			}
		}
		
		/* Initialize the phg statistics */
		PhoHStatInitPhgStatistics();
		
		/* Get start timing info */
		LbTmStartTiming(&calcDecaysStartTime);
		
		/* Calculate the number of decays */
		SubObjCalcTimeBinDecays(0, PhgRunTimeParams.Phg_EventsToSimulate);
		
		/* Get timing information */
		timingValid = LbTmStopTiming(&calcDecaysStartTime, &calcDecaysTime, &calcDecaysCPUTime);
		LbTmStartTiming(&trackPhotonsStartTime);
		
		/* Indicate we are about to simulate */
		LbInPrintf("\n\n***************** Simulation Beginning *************\n");

		/* Clear out our counter */
		EmisListNumCreated = 0;
		
		/* Loop through all decays */
		while (emLiCreateDecay()) {
			
			/* Although these get cleared out at the beginning of appropriate routines, it
				is better that they get cleared here.
			*/
			for (ColCurParams = 0; ColCurParams < ColNumParams; ColCurParams++) {
				EmisListCollimatedPhotons[ColCurParams].NumCollimatedBluePhotons = 0;
				EmisListCollimatedPhotons[ColCurParams].NumCollimatedPinkPhotons = 0;
			}
			for (DetCurParams = 0; DetCurParams < DetNumParams; DetCurParams++) {
				EmisListDetectedPhotons[DetCurParams].NumDetectedBluePhotons = 0;			
				EmisListDetectedPhotons[DetCurParams].NumDetectedPinkPhotons = 0;
			}
			EmisListCurBluePhotonIndex	 = 0;	
			EmisListCurPinkPhotonIndex = 0;
					
			
			#ifdef __MWERKS__
				/* Pause for user interaction (on mac) */
				EmisListDoMacEvent();
			#endif
			
			/* Reset current detected photon indexes */
			EmisListCurBluePhotonIndex = 0;
			EmisListCurPinkPhotonIndex = 0;
			
			/* Increment statistics for this start */
			PhoHStatUpdateStartedPhotons(&EmisListBluePhoton);
			
			/* Update primary productivity starts if photon will be tracked as primary */
			if ( PHG_IsTrackAsPrimary(&EmisListBluePhoton) ) {
			
				ProdTblAddStartingProductivity(&EmisListBluePhoton, PRODTBLFg_Primary);
			
			}
			
			/* Update scatter productivity starts if photon will be tracked as scatter */
			if ( PHG_IsTrackAsScatter(&EmisListBluePhoton) ) {
			
				ProdTblAddStartingProductivity(&EmisListBluePhoton, PRODTBLFg_Scatter);
			
			}
				
			/* Track that photon */
			EmisLisTrackPhoton(&EmisListBluePhoton);
			
			/* Track pink photon if we are doing PET and (1) there is a detected blue photon
			 OR (2) we are tracking singles */
			if ( PHG_IsPETCoincPlusSingles() ||
					( PHG_IsPETCoincidencesOnly() && (EmisListCurBluePhotonIndex != 0) )  ) {
				
				/* Increment statistics for this start */
				PhoHStatUpdateStartedPhotons(&EmisListPinkPhoton);
				
				/* Update primary productivity starts if photon will be tracked as primary */
				if ( PHG_IsTrackAsPrimary(&EmisListPinkPhoton) ) {
				
					ProdTblAddStartingProductivity(&EmisListPinkPhoton, PRODTBLFg_Primary);
				
				}
				
				/* Update scatter productivity starts if photon will be tracked as scatter */
				if ( PHG_IsTrackAsScatter(&EmisListPinkPhoton) ) {
				
					ProdTblAddStartingProductivity(&EmisListPinkPhoton, PRODTBLFg_Scatter);
				
				}
				
				/* Track that photon */
				EmisLisTrackPhoton(&EmisListPinkPhoton);
				
			}
			/* Collimate/detect/bin/write photons that made it to the target cylinder */
			if (PHG_IsPETCoincidencesOnly()){
			
				/* For PET check for both blue and pink photons to process detections */
				if ((EmisListCurBluePhotonIndex != 0) &&
						(EmisListCurPinkPhotonIndex != 0)) {
					
					/* Collimate them if requested */
					if (PHG_IsCollimateOnTheFly()) {
						
						/* Loop through possible collimator configurations */
						for (ColCurParams = 0; ColCurParams < ColNumParams; ColCurParams++) {
						
							/* Some collimator models rely on detector data so keep current detector 
							parameters index consistent just in case.
							*/
							DetCurParams = ColCurParams;

							/* Collimate the photons */
							ColPETPhotons(&EmisListNewDecay,
								EmisListDetectdTrkngBluePhotons, EmisListCurBluePhotonIndex,
								EmisListDetectdTrkngPinkPhotons, EmisListCurPinkPhotonIndex,
								&EmisListCollimatedPhotons[ColCurParams]);
									
								/* Restore tracked photons if more than one detector */
								if (ColNumParams > 1) {
									for (loopV = 0; loopV < EmisListCurBluePhotonIndex; loopV++) {
										EmisListDetectdTrkngBluePhotons[loopV] = EmisListTrackedBluePhotons[loopV];
									}
									for (loopV = 0; loopV < EmisListCurPinkPhotonIndex; loopV++) {
										EmisListDetectdTrkngPinkPhotons[loopV] = EmisListTrackedPinkPhotons[loopV];
									}
								}
						}
					}
					
					/* Detect them if requested */
					if (PHG_IsDetectOnTheFly()) {
			
						/* Send collimated photons if it was done */
						if (PHG_IsCollimateOnTheFly()) {
							for (DetCurParams = 0; DetCurParams < DetNumParams; DetCurParams++) {
							
								/* Keep collimator params index current with detector */
								ColCurParams = DetCurParams;
								
								DetPETPhotons(&EmisListNewDecay,
									EmisListCollimatedPhotons[ColCurParams].CollimatedTrkngBluePhotons,
									EmisListCollimatedPhotons[ColCurParams].NumCollimatedBluePhotons,
									EmisListCollimatedPhotons[ColCurParams].CollimatedTrkngPinkPhotons,
									EmisListCollimatedPhotons[ColCurParams].NumCollimatedPinkPhotons,
									&(EmisListDetectedPhotons[DetCurParams]));
							}
						}
						else {
							/* Loop through possible detector configurations */
							for (DetCurParams = 0; DetCurParams < DetNumParams; DetCurParams++) {

								ColCurParams = DetCurParams;
								
								DetPETPhotons(&EmisListNewDecay,
									EmisListDetectdTrkngBluePhotons, EmisListCurBluePhotonIndex,
									EmisListDetectdTrkngPinkPhotons, EmisListCurPinkPhotonIndex,
									&(EmisListDetectedPhotons[DetCurParams]));
									
								/* Restore tracked photons if more than one detector */
								if (DetNumParams > 1) {
									for (loopV = 0; loopV < EmisListCurBluePhotonIndex; loopV++) {
										EmisListDetectdTrkngBluePhotons[loopV] = EmisListTrackedBluePhotons[loopV];
									}
									for (loopV = 0; loopV < EmisListCurPinkPhotonIndex; loopV++) {
										EmisListDetectdTrkngPinkPhotons[loopV] = EmisListTrackedPinkPhotons[loopV];
									}
								}
							}
						}
					}
					
					/* Bin them up if binning is being done */
					if (PHG_IsBinOnTheFly()) {
						for (curBinParams = 0; curBinParams < PhgNumBinParams; curBinParams++) {
							ColCurParams = curBinParams;
							DetCurParams = curBinParams;
							
							/* Pass detected Photons if we created them */
							if (PHG_IsDetectOnTheFly()){
								PhgBinPETPhotons(&PhgBinParams[curBinParams],
									&PhgBinData[curBinParams], &PhgBinFields[curBinParams],
									&EmisListNewDecay,
									EmisListDetectedPhotons[DetCurParams].DetectedTrkngBluePhotons, 
									EmisListDetectedPhotons[DetCurParams].NumDetectedBluePhotons,
									EmisListDetectedPhotons[DetCurParams].DetectedTrkngPinkPhotons,
									EmisListDetectedPhotons[DetCurParams].NumDetectedPinkPhotons);
							}
							else if (PHG_IsCollimateOnTheFly()) {
								
								PhgBinPETPhotons(&PhgBinParams[curBinParams], &PhgBinData[curBinParams], &PhgBinFields[curBinParams],
									&EmisListNewDecay,
									EmisListCollimatedPhotons[ColCurParams].CollimatedTrkngBluePhotons, 
									EmisListCollimatedPhotons[ColCurParams].NumCollimatedBluePhotons,
									EmisListCollimatedPhotons[ColCurParams].CollimatedTrkngPinkPhotons,
									EmisListCollimatedPhotons[ColCurParams].NumCollimatedPinkPhotons);
							}
							else {
								/* Bin up non-collimated photons */
								PhgBinPETPhotons(&PhgBinParams[curBinParams], &PhgBinData[curBinParams], &PhgBinFields[curBinParams],
									&EmisListNewDecay,
									EmisListDetectdTrkngBluePhotons, EmisListCurBluePhotonIndex,
									EmisListDetectdTrkngPinkPhotons, EmisListCurPinkPhotonIndex);
							}
						}
					}
					
					/* Now write them out (Abort on failure) */
					if (PHG_IsHist()) {
						if (PhoHFileWriteDetections(&EmisListPHGHistoryFileHk, &EmisListNewDecay,
								EmisListTrackedBluePhotons, EmisListCurBluePhotonIndex,
								EmisListTrackedPinkPhotons, EmisListCurPinkPhotonIndex) == false) {
							
							/* Abort Program execution */
							PhgAbort("Got failure from PhoHFileWriteDetections (EmisListCreatePhotonList).", true);
						}
					}
				}
			}
			else if (PHG_IsPETCoincPlusSingles()){
			
				/* For PET check for both blue and pink photons to process detections */
				if ((EmisListCurBluePhotonIndex != 0) ||
						(EmisListCurPinkPhotonIndex != 0)) {
					
					/* Collimate them if requested */
					if (PHG_IsCollimateOnTheFly()) {
						
						/* Loop through possible collimator configurations */
						for (ColCurParams = 0; ColCurParams < ColNumParams; ColCurParams++) {
						
							/* Some collimator models rely on detector data so keep current detector 
							parameters index consistent just in case.
							*/
							DetCurParams = ColCurParams;

							/* Collimate the photons */
							ColPETPhotons(&EmisListNewDecay,
								EmisListDetectdTrkngBluePhotons, EmisListCurBluePhotonIndex,
								EmisListDetectdTrkngPinkPhotons, EmisListCurPinkPhotonIndex,
								&EmisListCollimatedPhotons[ColCurParams]);
									
								/* Restore tracked photons if more than one detector */
								if (ColNumParams > 1) {
									for (loopV = 0; loopV < EmisListCurBluePhotonIndex; loopV++) {
										EmisListDetectdTrkngBluePhotons[loopV] = EmisListTrackedBluePhotons[loopV];
									}
									for (loopV = 0; loopV < EmisListCurPinkPhotonIndex; loopV++) {
										EmisListDetectdTrkngPinkPhotons[loopV] = EmisListTrackedPinkPhotons[loopV];
									}
								}
						}
					}
					
					/* Detect them if requested */
					if (PHG_IsDetectOnTheFly()) {
			
						/* Send collimated photons if it was done */
						if (PHG_IsCollimateOnTheFly()) {
							for (DetCurParams = 0; DetCurParams < DetNumParams; DetCurParams++) {
							
								/* Keep collimator params index current with detector */
								ColCurParams = DetCurParams;
								
								DetPETPhotons(&EmisListNewDecay,
									EmisListCollimatedPhotons[ColCurParams].CollimatedTrkngBluePhotons,
									EmisListCollimatedPhotons[ColCurParams].NumCollimatedBluePhotons,
									EmisListCollimatedPhotons[ColCurParams].CollimatedTrkngPinkPhotons,
									EmisListCollimatedPhotons[ColCurParams].NumCollimatedPinkPhotons,
									&(EmisListDetectedPhotons[DetCurParams]));
							}
						}
						else {
							/* Loop through possible detector configurations */
							for (DetCurParams = 0; DetCurParams < DetNumParams; DetCurParams++) {

								ColCurParams = DetCurParams;
								
								DetPETPhotons(&EmisListNewDecay,
									EmisListDetectdTrkngBluePhotons, EmisListCurBluePhotonIndex,
									EmisListDetectdTrkngPinkPhotons, EmisListCurPinkPhotonIndex,
									&(EmisListDetectedPhotons[DetCurParams]));
									
								/* Restore tracked photons if more than one detector */
								if (DetNumParams > 1) {
									for (loopV = 0; loopV < EmisListCurBluePhotonIndex; loopV++) {
										EmisListDetectdTrkngBluePhotons[loopV] = EmisListTrackedBluePhotons[loopV];
									}
									for (loopV = 0; loopV < EmisListCurPinkPhotonIndex; loopV++) {
										EmisListDetectdTrkngPinkPhotons[loopV] = EmisListTrackedPinkPhotons[loopV];
									}
								}
							}
						}
					}
					
					/* Bin them up if binning is being done */
					if (PHG_IsBinOnTheFly()) {
						for (curBinParams = 0; curBinParams < PhgNumBinParams; curBinParams++) {
							ColCurParams = curBinParams;
							DetCurParams = curBinParams;
							
							/* Pass detected Photons if we created them */
							if (PHG_IsDetectOnTheFly()){
							
								if ( PhgBinParams[0].isBinPETasSPECT ) {
									
									PhgBinSPECTPhotons(&PhgBinParams[curBinParams],
										&PhgBinData[curBinParams], &PhgBinFields[curBinParams],
										&EmisListNewDecay,
										EmisListDetectedPhotons[DetCurParams].DetectedTrkngBluePhotons,
										EmisListDetectedPhotons[DetCurParams].NumDetectedBluePhotons);
									
									PhgBinSPECTPhotons(&PhgBinParams[curBinParams],
										&PhgBinData[curBinParams], &PhgBinFields[curBinParams],
										&EmisListNewDecay,
										EmisListDetectedPhotons[DetCurParams].DetectedTrkngPinkPhotons,
										EmisListDetectedPhotons[DetCurParams].NumDetectedPinkPhotons);
									
								} else {
								
									PhgBinPETPhotons(&PhgBinParams[curBinParams],
										&PhgBinData[curBinParams], &PhgBinFields[curBinParams],
										&EmisListNewDecay,
										EmisListDetectedPhotons[DetCurParams].DetectedTrkngBluePhotons, 
										EmisListDetectedPhotons[DetCurParams].NumDetectedBluePhotons,
										EmisListDetectedPhotons[DetCurParams].DetectedTrkngPinkPhotons,
										EmisListDetectedPhotons[DetCurParams].NumDetectedPinkPhotons);
									
								}
								
							}
							else if (PHG_IsCollimateOnTheFly()) {
								
								if ( PhgBinParams[0].isBinPETasSPECT ) {
									
									PhgBinSPECTPhotons(&PhgBinParams[curBinParams],
										&PhgBinData[curBinParams], &PhgBinFields[curBinParams],
										&EmisListNewDecay,
										EmisListCollimatedPhotons[ColCurParams].CollimatedTrkngBluePhotons, 
										EmisListCollimatedPhotons[ColCurParams].NumCollimatedBluePhotons);
									
									PhgBinSPECTPhotons(&PhgBinParams[curBinParams],
										&PhgBinData[curBinParams], &PhgBinFields[curBinParams],
										&EmisListNewDecay,
										EmisListCollimatedPhotons[ColCurParams].CollimatedTrkngPinkPhotons,
										EmisListCollimatedPhotons[ColCurParams].NumCollimatedPinkPhotons);
									
								} else {
								
									PhgBinPETPhotons(&PhgBinParams[curBinParams], &PhgBinData[curBinParams], &PhgBinFields[curBinParams],
										&EmisListNewDecay,
										EmisListCollimatedPhotons[ColCurParams].CollimatedTrkngBluePhotons, 
										EmisListCollimatedPhotons[ColCurParams].NumCollimatedBluePhotons,
										EmisListCollimatedPhotons[ColCurParams].CollimatedTrkngPinkPhotons,
										EmisListCollimatedPhotons[ColCurParams].NumCollimatedPinkPhotons);
									
								}
								
							}
							else {

								/* Bin up non-collimated photons */
								if ( PhgBinParams[0].isBinPETasSPECT ) {
									
									PhgBinSPECTPhotons(&PhgBinParams[curBinParams],
										&PhgBinData[curBinParams], &PhgBinFields[curBinParams],
										&EmisListNewDecay,
										EmisListDetectdTrkngBluePhotons, EmisListCurBluePhotonIndex);
									
									PhgBinSPECTPhotons(&PhgBinParams[curBinParams],
										&PhgBinData[curBinParams], &PhgBinFields[curBinParams],
										&EmisListNewDecay,
										EmisListDetectdTrkngPinkPhotons, EmisListCurPinkPhotonIndex);
									
								} else {
								
									PhgBinPETPhotons(&PhgBinParams[curBinParams], &PhgBinData[curBinParams], &PhgBinFields[curBinParams],
										&EmisListNewDecay,
										EmisListDetectdTrkngBluePhotons, EmisListCurBluePhotonIndex,
										EmisListDetectdTrkngPinkPhotons, EmisListCurPinkPhotonIndex);
									
								}
								
							}
						}
					}
					
					/* Now write them out (Abort on failure) */
					if (PHG_IsHist()) {
						if (PhoHFileWriteDetections(&EmisListPHGHistoryFileHk, &EmisListNewDecay,
								EmisListTrackedBluePhotons, EmisListCurBluePhotonIndex,
								EmisListTrackedPinkPhotons, EmisListCurPinkPhotonIndex) == false) {
							
							/* Abort Program execution */
							PhgAbort("Got failure from PhoHFileWriteDetections (EmisListCreatePhotonList).", true);
						}
					}
				}
			}
			else {
			
				/* For SPECT check for blue photos to process detections */
				if (EmisListCurBluePhotonIndex != 0) {
						
					/* Collimate them if necessary */
					if (PHG_IsCollimateOnTheFly()) {
						for (ColCurParams = 0; ColCurParams < ColNumParams; ColCurParams++) {
							DetCurParams = ColCurParams;
							ColSPECTPhotons(&EmisListNewDecay,
								EmisListDetectdTrkngBluePhotons, EmisListCurBluePhotonIndex,
								&EmisListCollimatedPhotons[ColCurParams]);
									
								/* Restore tracked photons if more than one detector */
								if (ColNumParams > 1) {
									for (loopV = 0; loopV < EmisListCurBluePhotonIndex; loopV++) {
										EmisListDetectdTrkngBluePhotons[loopV] = EmisListTrackedBluePhotons[loopV];
									}
									for (loopV = 0; loopV < EmisListCurPinkPhotonIndex; loopV++) {
										EmisListDetectdTrkngPinkPhotons[loopV] = EmisListTrackedPinkPhotons[loopV];
									}
								}
						}
					}
						
					/* Collimate them if necessary */
					if (PHG_IsDetectOnTheFly()) {
						for (DetCurParams = 0; DetCurParams < DetNumParams; DetCurParams++) {
							ColCurParams = DetCurParams;
							
							if (PHG_IsCollimateOnTheFly()) {
								DetSPECTPhotons(&EmisListNewDecay,
									EmisListCollimatedPhotons[ColCurParams].CollimatedTrkngBluePhotons,
									EmisListCollimatedPhotons[ColCurParams].NumCollimatedBluePhotons,
									&(EmisListDetectedPhotons[DetCurParams]));
							}
							else {
								DetSPECTPhotons(&EmisListNewDecay,
									EmisListDetectdTrkngBluePhotons, EmisListCurBluePhotonIndex,
									&(EmisListDetectedPhotons[DetCurParams]));
									
								/* Restore tracked photons if more than one detector */
								if (DetNumParams > 1) {
									for (loopV = 0; loopV < EmisListCurBluePhotonIndex; loopV++) {
										EmisListDetectdTrkngBluePhotons[loopV] = EmisListTrackedBluePhotons[loopV];
									}
									for (loopV = 0; loopV < EmisListCurPinkPhotonIndex; loopV++) {
										EmisListDetectdTrkngPinkPhotons[loopV] = EmisListTrackedPinkPhotons[loopV];
									}
								}
							}
						}
					}

					/* Bin them up if necessary */
					if (PHG_IsBinOnTheFly()) {

						for (curBinParams = 0; curBinParams < PhgNumBinParams; curBinParams++) {
							ColCurParams = curBinParams;
							DetCurParams = curBinParams;
							
							/* Pass detected photons if we created them */
							if (PHG_IsDetectOnTheFly()) {
								PhgBinSPECTPhotons(&PhgBinParams[curBinParams],
									&PhgBinData[curBinParams], &PhgBinFields[curBinParams],
									&EmisListNewDecay,
									EmisListDetectedPhotons[DetCurParams].DetectedTrkngBluePhotons,
									EmisListDetectedPhotons[DetCurParams].NumDetectedBluePhotons);
							}
							else if (PHG_IsCollimateOnTheFly()) {
								PhgBinSPECTPhotons(&PhgBinParams[curBinParams],
									&PhgBinData[curBinParams], &PhgBinFields[curBinParams],
									&EmisListNewDecay,
									EmisListCollimatedPhotons[ColCurParams].CollimatedTrkngBluePhotons,
									EmisListCollimatedPhotons[ColCurParams].NumCollimatedBluePhotons);
							}
							else {
								/* Bin up non-collimated photons */
								PhgBinSPECTPhotons(&PhgBinParams[curBinParams], &PhgBinData[curBinParams], &PhgBinFields[curBinParams],
								&EmisListNewDecay,
									EmisListDetectdTrkngBluePhotons, EmisListCurBluePhotonIndex);
							}
						}
					}

					/* Now write them out (Abort on failure) */
					if (PHG_IsHist()) {
						if (PhoHFileWriteDetections(&EmisListPHGHistoryFileHk, &EmisListNewDecay,
								EmisListTrackedBluePhotons, EmisListCurBluePhotonIndex,
								EmisListTrackedPinkPhotons, EmisListCurPinkPhotonIndex) == false) {
							
							/* Abort Program execution */
							PhgAbort("Unexpected failure from PhoHFileWriteDetections trapped in EmisListCreatePhotonList.", true);
						}
					}
				}
			}
			
			/* Although these get cleared out at the beginning appropriate routines, they need to get cleared out
				at the end also
			*/
			for (ColCurParams = 0; ColCurParams < ColNumParams; ColCurParams++) {
				EmisListCollimatedPhotons[ColCurParams].NumCollimatedBluePhotons = 0;
				EmisListCollimatedPhotons[ColCurParams].NumCollimatedPinkPhotons = 0;
			}
			for (DetCurParams = 0; DetCurParams < DetNumParams; DetCurParams++) {
				EmisListDetectedPhotons[DetCurParams].NumDetectedBluePhotons = 0;			
				EmisListDetectedPhotons[DetCurParams].NumDetectedPinkPhotons = 0;
			}
			EmisListCurBluePhotonIndex	 = 0;	
			EmisListCurPinkPhotonIndex = 0;
					
		} /* End of loop for tracking photons */
		
		/* Get timing information */
		timingValid = LbTmStopTiming(&trackPhotonsStartTime, &trackPhotonsTime, &trackPhotonsCPUTime);
		
		LbInPrintf("\n***************** Simulation Finished *************\n\n");
		
		/* Tell them how many photons we tracked */
		if (PHG_IsPET()) {
			LbInPrintf("\n\nTracked %lld photons.", SUBOBJGetDecaysProcessed()*2);
		}
		else {
			LbInPrintf("\n\nTracked %lld photons.", SUBOBJGetDecaysProcessed());
		}
		
		/* Print out timing information */
		{
			/* Only do this when we want it */
			if (timingValid) {
				LbInPrintf("\n\nReal time for calculating time bin decays = %3.1f seconds.", 
						calcDecaysTime);
				
				LbInPrintf("\nReal time for tracking photons = %3.1f seconds.", 
						trackPhotonsTime);
				
				LbInPrintf("\nTotal CPU time for calculating time bin decays = %3.1f seconds.",
						calcDecaysCPUTime);
				
				LbInPrintf("\nTotal CPU time for tracking photons  = %3.1f. seconds",
						trackPhotonsCPUTime);
				
				LbInPrintf("\n\nCPU time for simulation =  %3.1f.",
						(calcDecaysCPUTime + trackPhotonsCPUTime));
			}
		}
		
		/* Write the photon statistics */
		if (!PhoHStatWrite()) {
			break;
		}
		
		/* Close the photon history list */
		if (PHG_IsHist()) {
					
			/* Print out history report */
			if (EmisListPHGHistoryFileHk.doCustom) {
				LbInPrintf("\n\n\tHistory file report for PHG history file");
				PhoHFilePrintReport(&EmisListPHGHistoryFileHk);
			}
			
			if (!PhoHFileClose(&EmisListPHGHistoryFileHk)) {
				break;
			}
		}
		
		#ifdef __MWERKS__
			ABORT:;
		#endif
		okay = true;
	} while (false);
	
	/* If an error occured, do best to clean up */
	if (!okay) {
		ErIgnoreBegin();
			if (histFileCreated)
				(void) PhoHFileClose(&EmisListPHGHistoryFileHk);
		ErIgnoreEnd();
	}
	
	return (okay);
}


/*********************************************************************************
*
*			Name:		emLiDoEscape
*
*			Summary:	Do the things you do when a photon escapes.
*			Arguments:
*				
*			Function return: None.
*
*********************************************************************************/
void emLiDoEscape()	
{
	/* Increment escaped photon statistics */
	PhoHStatIncTotEscPhotons();
}

/*********************************************************************************
*
*			Name:		EmisListInitialize
*
*			Summary:	Initialize the emission list manager. This routine is here
*						to fascilitate debug code that needs
*						initialization/termination to work.
*			Arguments:
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean EmisListInitialize()	
{
	
	do { /* Process Loop */
		
		/* Initialize arrays for polarization test */
		#ifdef PHG_DEBUG
		LbFourByte	i, j;
		
		if (PHG_IsModelPolarization()) {
			for (i = 0; i < 150; i++) {
				emLiDBPolarTestBlue[i] = 0;
				emLiDBPolarTestPink[i] = 0;
			}
			emLiDBPolarTotalBlue = 0;
			emLiDBPolarTotalPink = 0;
			
			for (i = 0; i < 50; i++) {
			for (j = 0; j < 50; j++) {
				emLiDBPolarPhiMuWt[i][j] = 0.0;
				emLiDBPolarPhiMuWtSq[i][j] = 0.0;
				emLiDBPolarPhiMuCt[i][j] = 0;
			}
			}
		}
		#endif
#ifdef DO_RH_POS_RANGE
{
FILE	*sigmaFile;
char	sigma[128];

		/* Open the input file */
		if ((sigmaFile = LbFlFileOpen("sigma", "r")) == 0) {
			ErStFileError("Unable to open sigma file should be named 'sigma' (EmisListInitialize).");
			PhgAbort("Terminating", false);
		}
	
		if (LbFlFGetS(sigma, 128, sigmaFile) == NULL) {
			ErStFileError("Unable to read sigma file (EmisListInitialize).");
			PhgAbort("Terminating", false);
		}

		emLiDbSigma = atof(sigma);
}	
#endif
		
		/* You can set EmisListBreakNumber to facilitate stopping on a specific
			photon in EmisLisTrackPhoton. Setting it to -1 avoids the debug
			check.
		*/
		
		/* Allocate memory for detected photons block */
		if ((EmisListDetectdTrkngBluePhotons = (PHG_TrackingPhoton *) 
				LbMmAlloc(sizeof(PHG_TrackingPhoton) * PHG_MAX_DETECTED_PHOTONS)) == 0) {
			break;
		}
		
		/* Allocate memory for detected photons block */
		if ((EmisListDetectdTrkngPinkPhotons = (PHG_TrackingPhoton *) 
				LbMmAlloc(sizeof(PHG_TrackingPhoton) * PHG_MAX_DETECTED_PHOTONS)) == 0) {
			
			LbMmFree((void **)&EmisListDetectdTrkngBluePhotons);
			
			break;
		}
		
		/* Allocate memory for tracked photons block */
		if ((EmisListTrackedBluePhotons = (PHG_TrackingPhoton *) 
				LbMmAlloc(sizeof(PHG_TrackingPhoton) * PHG_MAX_DETECTED_PHOTONS)) == 0) {
			break;
		}
		
		/* Allocate memory for detected photons block */
		if ((EmisListTrackedPinkPhotons = (PHG_TrackingPhoton *) 
				LbMmAlloc(sizeof(PHG_TrackingPhoton) * PHG_MAX_DETECTED_PHOTONS)) == 0) {
			
			LbMmFree((void **)&EmisListDetectdTrkngBluePhotons);
			
			break;
		}
		
		#ifdef PHG_DEBUG
		
			emLiSumSingleCohScatters = 0.0;
			emLiSumSingleCompScatters = 0.0;
			emLiSumSingleScatAborptions = 0.0;
			if (PHG_IsNonCollinearityAdjust())
			{
				int i;
				for (i = 0; i < NUM_MU_DIST_BINS; i++) {
					emlsMuDistribution[i] = 0;
				}
			}
		#endif
		
		
		if (PHG_IsRangeAdjust() == true) {
		{
		char			isotope[32];
		char			rangeParam[32];
		FILE			*rangeFile;
		LbUsFourByte	ki;
		LbUsFourByte	skipI;
		
		#ifdef PHG_DEBUG
		/* Get the attenuation of water at 1000kev */
		if (PHG_IsModelCoherentInObj()) {
			 emLisAttWater1000keV = 0.07061;
		}
		else {
			 emLisAttWater1000keV = 0.070554220889095007;
		}
		#endif
		
		/* Open isotope data file */
		if ( (rangeFile = LbFlFileOpen(EmisListIsotopeDataFilePath, "rb")) == NULL ) {

			PhgAbort("Can't open input file for positron range data file. (EmisListInitialize)", false);
		}
		
		/* Loop until isotope found */
		do {
			/* Read in name of isotope */
			if ( (LbFlFGetS(isotope, 32,rangeFile)) ==  '\0' ) {
					PhgAbort("Didn't read in isotope name from positron range file as expected. (EmisListInitialize)", false);
			}
			
			/* Lop off the new line */
			isotope[strlen(isotope)-1] = '\0';
			
			/* See if this is the isotope we are modelling */
			if (strcmp(phgEn_IsotopeStr[(int) PhgRunTimeParams.PhgNuclide.isotope], isotope) == 0) {

				/* Seek to next isotope */
				for (ki = 0; ki < 100; ki++) {
					/* Read in  isotope */
					if ( (LbFlFGetS(rangeParam, 32, rangeFile)) == '\0' ) {
							PhgAbort("Error reading past unwanted positron range data. (EmisListInitialize)", false);
					}
					emLiIsotopeParams[ki] = atof(rangeParam);
				}
				
				/* Break out of the loop, we are done */
				break;
				
			}
			else {
				/* Seek to next isotope */
				for (skipI = 0; skipI < 100; skipI++) {
					/* Read in name of isotope */
					if ( (LbFlFGetS(rangeParam, 32, rangeFile)) == '\0' ) {
							PhgAbort("Error reading past unwanted positron range data. (EmisListInitialize)", false);
					}
				}
			}
		} while (true);
		fclose(rangeFile);
		}
		}
		EmisListIsInitialized = true;
	} while (false);
	
	return (EmisListIsInitialized);
}

/*********************************************************************************
*
*			Name:		emLiGtIsotopeRange
*
*			Summary:	Returns a "freePaths" range for positron range.
*
*			Arguments:
*			Function return: None.
*
*********************************************************************************/
#ifdef DO_POS_RANGE_THE_OLD_WAY
double emLiGtIsotopeRange()	
{
double			tempr;
LbUsFourByte	rangeIndex;
double			range;

	/* generate a random number from 0 to 1000 and turn it into an index */
	tempr = 1000.0 * PhgMathGetRandomNumber();
	rangeIndex = floor( tempr );

	/* if the index is 0 to 989, get the range by interpolating from the cumulativeD array */
	if ( rangeIndex < 990 ) {
		range = emLiEmpiricalRangeDist[rangeIndex] +
			( (emLiEmpiricalRangeDist[rangeIndex+1] - emLiEmpiricalRangeDist[rangeIndex]) *
			(1- (tempr - rangeIndex)) );
		
	/* if the index >= 990, use the slower decaying part of the bi-exponential density
		to determine the range */
	} else {
			tempr = PhgMathGetRandomNumber();
			range = emLiEmpiricalRangeDist[990] - emLiIsotopeParams.c*log(1-tempr);
	}
	
	/* Adjust for attenuation in water */
	range  *= emLisAttWater1000keV;
	
	return(range);
}
#endif
/*********************************************************************************
*
*			Name:		EmisListTerminate
*
*			Summary:	Terminate the emission list manager. This routine is here
*						to fascilitate debug code that needs
*						initialization/termination to work.
*
*			Arguments:
*			Function return: None.
*
*********************************************************************************/
void EmisListTerminate()	
{
	if (EmisListDetectdTrkngBluePhotons != 0) 
		LbMmFree((void **)&EmisListDetectdTrkngBluePhotons);

	if (EmisListDetectdTrkngPinkPhotons != 0) 
		LbMmFree((void **)&EmisListDetectdTrkngPinkPhotons);

	if (EmisListTrackedBluePhotons != 0) 
		LbMmFree((void **)&EmisListTrackedBluePhotons);

	if (EmisListTrackedPinkPhotons != 0) 
		LbMmFree((void **)&EmisListTrackedPinkPhotons);

	#ifdef PHG_DEBUG_DEBUG_COHERENT
		if (ErIsInError() == false) {
		LbInPrintf("\nWeight of photons leaving object with 1 coherent scatter = %3.4e", emLiSumSingleCohScatters);
		LbInPrintf("\nWeight of photons leaving object with 1 compton scatter = %3.4e", emLiSumSingleCompScatters);
		LbInPrintf("\nWeight of photons absorbed in object with 1 scatter = %3.4e\n", emLiSumSingleScatAborptions);
	}

	#endif
	
	#ifdef PHG_DEBUG
	{	
		#ifdef NONCOLL_TEST
		{
			if (PHG_IsNonCollinearityAdjust())
			{
				int i;
				LbInPrintf("\nDistribution of MU values in non-collinearity computation is\n");
				for (i = 0; i < NUM_MU_DIST_BINS; i++) {
					LbInPrintf("%ld\n", (long)emlsMuDistribution[i]);
				}
			}
			LbInPrintf("\n");
		}
		#endif

		#ifdef POSRANGE_TEST
		{
			/* !RH added following prints... */
			/* Print out positron emission energies and distances */
			if (PHG_IsRangeAdjust())
			{
				LbUsFourByte i;
			
				LbInPrintf("\nDistribution of positron emission energies\n");
				#ifdef POSRANGE_TEST_A1
				        for (i = 0; i < 1000; i++){
				                LbInPrintf("%ld\n", (unsigned long)emLiPosEnergy[i]);
				        }
				#else
				        for (i = 0; i < 2000; i++){
				                LbInPrintf("%ld\n", (unsigned long)emLiPosEnergy[i]);
				        }
				#endif
				
				/* !!RH new POSRANGE test print statements:  */
				LbInPrintf("\n");
				LbInPrintf("\n");
				for (i = 0; i < 2000; i++){
					LbInPrintf("%ld\n", (unsigned long)emLiPosRangeDistance[i]);
				}
				/* !!RH end new  */
				
			}
		}
		#endif
	}
	#endif

	#ifdef PHG_DEBUG	
	{
		LbFourByte i, j;
			
		/* Compute polarization variable if it is being done */
#define DO_TEST_E
#ifdef DO_TEST_A
		if (PHG_IsModelPolarization()) {
		LbInPrintf("\nTotal blue scatters = %d, total pink scatters = %d", emLiDBPolarTotalBlue,emLiDBPolarTotalPink);
		emLiDBPolarTotalPink = 0;

		LbInPrintf("\nBlue followed by Pink Polarization Phi's\n");
			for (i = 0; i < 150; i++) {
				LbInPrintf("%d\t", emLiDBPolarTestBlue[i]);
			}
			LbInPrintf("\n");
			for (i = 0; i < 150; i++) {
				LbInPrintf("%d\t", emLiDBPolarTestPink[i]);
			}
			LbInPrintf("\n");
		}
#endif
#ifdef DO_TEST_E
		if (PHG_IsModelPolarization()) {
		LbInPrintf("\n");
		LbInPrintf("\nmu X phi table weight\n");
			for (i = 0; i < 50; i++) {
			for (j = 0; j < 50; j++) {
				LbInPrintf("%3.4f\t", emLiDBPolarPhiMuWt[i][j]);
			}
			LbInPrintf("\n");
			}
		LbInPrintf("\n");
		LbInPrintf("\nmu X phi table weight squared \n");
			for (i = 0; i < 50; i++) {
			for (j = 0; j < 50; j++) {
				LbInPrintf("%3.4f\t", emLiDBPolarPhiMuWtSq[i][j]);
			}
			LbInPrintf("\n");
			}
		LbInPrintf("\n");
		LbInPrintf("\nmu X phi table counts\n");
			for (i = 0; i < 50; i++) {
			for (j = 0; j < 50; j++) {
				LbInPrintf("%d\t", emLiDBPolarPhiMuCt[i][j]);
			}
			LbInPrintf("\n");
			}
		}
#endif

	}
	#endif
	
	EmisListIsInitialized = false;
	
}

/*********************************************************************************
*
*			Name:		EmisLisTrackPhoton
*
*			Summary:	Track a photon from the photon list.
*			Arguments:
*				PHG_TrackingPhoton	*trackingPhotonPtr	- The tracking photon.
*			Function return: None.
*
*********************************************************************************/
void EmisLisTrackPhoton(PHG_TrackingPhoton *trackingPhotonPtr)	
{
	double				totalFreePathLength;			/* Total free paths */
	double				freePathsUsed;					/* Free paths used on each move */
	double				distanceTraveled;				/* Distance photon travels */
	double				comptonToScatterProbability;	/* Probability of compton to probability of scatter */
	double				interactionProbability;			/* Random number for selecting interaction type */
	double				scatterProbability;				/* Probability of any type of scatter */
	PHG_Position		newPosition;					/* New photon position */
	PhoTrkActionTy		action;							/* Action to be taken based on new position */
	
	#ifdef PHG_DEBUG
	
			/* See if this is the photon we want to break on */
			if ((LbEightByte)(trackingPhotonPtr->number) == EmisListBreakNumber)
				LbInPrintf("Set breakpoint here.\n");
				
			/* Clear the FD photon number for debugging purposes */
			EmisListFDPhoton.number = 0;
	#endif
	
	/* Force detection if requested (FD is ON) */
	if ((PHG_IsForcedDetection()) && 
			(PHG_IsTrackAsPrimary(trackingPhotonPtr))) {
		
		/* Copy the photon */
		{
			EmisListFDPhoton = 	*trackingPhotonPtr;
		}
		
		/* Attempt initial forced detection  */
		PhoTrkAttInitForcedDetection(&EmisListFDPhoton);
		
		#ifdef PHG_DEBUG
			/* Clear the FD photon number for debugging purposes */
			EmisListFDPhoton.number = 0;
		#endif
			
	}


	do { /* Track the photon */

		/* If forced detection is on, and this is not a track as scatter photon we are done */
		if (PHG_IsForcedDetection()  && (!PHG_IsTrackAsScatter(trackingPhotonPtr))) {
			
			/* Terminate the loop */
			break;
		}
	
		/* Force detection if requested (FD is ON) */
		if (PHG_IsForcedDetection()) {
			
			/* Copy the photon */
			{
				EmisListFDPhoton = 	*trackingPhotonPtr;	
			}
			
			/* Attempt forced detection  */
			PhoTrkAttemptForcedDetection(&EmisListFDPhoton);
		}
	
		/* Choose a random number of free path lengths for the photon to travel */
		PhgMathGetTotalFreePaths(&totalFreePathLength);

		/* Determine stopping point of the photon */
		action = PhoTrkCalcNewPosition(trackingPhotonPtr, trackingPhotonPtr->location, trackingPhotonPtr->angle,
			trackingPhotonPtr->energy, totalFreePathLength, &newPosition, &distanceTraveled,
			&freePathsUsed);

		/* Set new position */
		trackingPhotonPtr->location = newPosition;
		
		/* Store the distance traveled */
		trackingPhotonPtr->travel_distance += distanceTraveled;
		
		#ifdef PHG_DEBUG
		/* Compute sum weight of single coherent scatters */
		if ((action != PhoTrkInteract) && (trackingPhotonPtr->num_of_scatters == 1) && (trackingPhotonPtr->energy == PhgRunTimeParams.PhgNuclide.photonEnergy_KEV)) {
			emLiSumSingleCohScatters += (trackingPhotonPtr->photon_scatter_weight * trackingPhotonPtr->decay_weight);
		}
		else if ((action != PhoTrkInteract) && (trackingPhotonPtr->num_of_scatters == 1)) {
			emLiSumSingleCompScatters += (trackingPhotonPtr->photon_scatter_weight * trackingPhotonPtr->decay_weight);
		}
		#endif
		
		/* See if the photon is out of the cylinder and forced detection is on, we are done */
		if ((action != PhoTrkInteract) && (PHG_IsForcedDetection())) {
	
			/* Photon left cylinder with forced detection on, so we are done tracking */
			break;
		}
		else if (action != PhoTrkInteract) {
				
			/* Record escape */
			emLiDoEscape();
		}
		
		/* Deal with action specific tasks */
		switch (action) {
			
			/* If the photon was detected (hit the target) */
			case PhoTrkDetect:
				
				/* Handle the detection */
				if ((trackingPhotonPtr->num_of_scatters > 0) ||
						(PHG_IsTrackAsPrimary(trackingPhotonPtr))) {
					
					EmisListDoDetection(trackingPhotonPtr);
				}

				/* Terminate the switch */
				break;
			
			/* If the photon should be discarded */
			case PhoTrkDiscard:

				/* Terminate the switch */
				break;
				
			/* If the photon interacted */
			case PhoTrkInteract:

				/* If not tracking as scatter, we are done */
				if (!PHG_IsTrackAsScatter(trackingPhotonPtr)) {

					/* Increment statistic */
					PhoHStatIncTotPrimOnlyScatter();

					/* Mark the photon for discarding*/
					action = PhoTrkDiscard;
					
					/* Terminate the switch */
					break;
				}
				
				/* Get interaction type probability */
				comptonToScatterProbability = SubObjGetProbComptToScatter(trackingPhotonPtr);
				scatterProbability = SubObjGetProbScatterInObj(trackingPhotonPtr);
				interactionProbability = PhgMathGetRandomNumber();

				/* See if non-absorption is being forced  */
				if (PHG_IsNoForcedNonAbsorbtion() == false){
				
					/* Adjust the weight */
					trackingPhotonPtr->photon_scatter_weight *= scatterProbability;

					/* See if compton scatter should be done */
					if ((interactionProbability < comptonToScatterProbability)) {
				
						/* Model Compton scatter (Using Klein-Nishina) */
						EmisListDoComptonInteraction(trackingPhotonPtr);
					}
					/* See if we're modelling coherent scatter */
					else if (PHG_IsModelCoherentInObj()) {
						EmisListDoCoherent(trackingPhotonPtr, SUBOBJGetTissueIndex(trackingPhotonPtr));
					}
				}
				else {
					/* See if absorption should be done */
					if (interactionProbability > scatterProbability) {
								
						/* Mark the photon for discarding*/
						action = PhoTrkDiscard;
						
						/* Update the statistics */
						PhoHStatIncTotAbsorbedPhotons(trackingPhotonPtr);

						#ifdef PHG_DEBUG
						if (trackingPhotonPtr->num_of_scatters == 0) {
							emLiSumSingleScatAborptions += (trackingPhotonPtr->photon_primary_weight * trackingPhotonPtr->decay_weight);
						}
						#endif
						
						/* Terminate the switch */
						break;
						
					}	/* See if coherent should be done */
					else if (interactionProbability > (scatterProbability * comptonToScatterProbability	)) {
							EmisListDoCoherent(trackingPhotonPtr, SUBOBJGetTissueIndex(trackingPhotonPtr));
					}
					else {		
						/* Model Compton scatter (Using Klein-Nishina) */
						EmisListDoComptonInteraction(trackingPhotonPtr);
					}
				}
	
				/* See if interaction is cause for discarding */
				if (trackingPhotonPtr->energy < PhgRunTimeParams.PhgMinimumEnergy) {
				
					/* Account for the photon loss */
					PhoHStatIncTotLowEnPhotons();
					
					/* Mark the photon for discard */
					action = PhoTrkDiscard;
					
					/* Terminate the tracking loop */
					break;
				}					
				
				/* Increment the number of scatters */
				trackingPhotonPtr->num_of_scatters++;

				/* Add the new start to the starts_list */
				if (!emLiCreateNewStart(trackingPhotonPtr))
					PhgAbort("Unable to add a new start to the starts list. (EmisListTrackPhoton).", true);

				/* Add starting productivity */
				ProdTblAddStartingProductivity(trackingPhotonPtr,
					PRODTBLFg_Scatter);
					   
				/* Terminate the switch */
				break;
				
			default:
				PhgAbort("Unexpected result returned from PhoTrkCalcNewPosition trapped in EmisListTrackPhoton.",
					true);
				break;
		}
	} while (action == PhoTrkInteract);	
}

/*********************************************************************************
*
*			Name:		emLiSamplePositronEnergy
*
*			Summary:	Returns the emission energy of a positron.
*
*			Arguments:
*			Function return: None.
*
*********************************************************************************/
double emLiSamplePositronEnergy()
{
	double			tempr;
	LbUsFourByte	energyIndex;
	double			energy;
	
	/* generate a random number from 0 to 100 and turn it into an index */
	tempr = 100.0 * PhgMathGetRandomNumber();
	energyIndex = floor( tempr );

	/* get the energy by interpolating from the cumulative Positron Energy array */
	if ( energyIndex > 99 ) {
		energyIndex = 99;
	}
	if ( energyIndex != 0 ) {
		energy = emLiIsotopeParams[energyIndex-1] +
			( (emLiIsotopeParams[energyIndex] - emLiIsotopeParams[energyIndex-1]) *
			(tempr - energyIndex) );
	} else {
			energy = emLiIsotopeParams[0] * tempr;
	}

	return(energy);
}


/*********************************************************************************
*
*			Name:		emLiComputePosRangeWater
*
*			Summary:	Returns a positron range for water sampled from an 
*						energy-dependent Gaussian.  Uses the algorithm given by
*						Palmer and Brownell, IEEE TMI 11:3:373-378.
*			Arguments:
*				double			positronEnergy	The emitted energy of the positron.
*				(the following two arguments are computed and returned)
*				double			*sigmaWaterPtr		The standard deviation of positron range in water for this energy.
*				PHG_Direction	*positronDirectionPtr		The direction to project the positron.
*			Function return: double--the range of the positron in water.
*
*********************************************************************************/
double emLiComputePosRangeWater( double positronEnergy, double *sigmaWaterPtr, PHG_Direction *positronDirectionPtr )
{
	double				rangeInWater;		/* sampled distance for this positron to travel through water */
	double				rangeExtrapolated;	/* the 'maximum' range */
	PHG_Position		sampleRange;		/* sample from 3D Gaussian */
	double				b1Water, b2Water, densityWater;		/* positron range constants for water */

	/* get positron range constants for water */
	SubObjGetWaterPosRangeConstants(&b1Water, &b2Water, &densityWater);

	/* compute the standard deviation of the Gaussian */
	rangeExtrapolated = (0.1 * b1Water * positronEnergy * positronEnergy) / (b2Water + positronEnergy);
	*sigmaWaterPtr = rangeExtrapolated / (2.0 * densityWater);
	
	sampleRange.x_position = PhgMathSampleFromGauss(0.0, *sigmaWaterPtr);
	sampleRange.y_position = PhgMathSampleFromGauss(0.0, *sigmaWaterPtr);
	sampleRange.z_position = PhgMathSampleFromGauss(0.0, *sigmaWaterPtr);
	
	rangeInWater = PHGMATH_SquareRoot( sampleRange.x_position*sampleRange.x_position +
										sampleRange.y_position*sampleRange.y_position +
										sampleRange.z_position*sampleRange.z_position );
	
	positronDirectionPtr->cosine_x = sampleRange.x_position / rangeInWater;
	positronDirectionPtr->cosine_y = sampleRange.y_position / rangeInWater;
	positronDirectionPtr->cosine_z = sampleRange.z_position / rangeInWater;
	
	return(rangeInWater);
}

/*********************************************************************************
*
*			Name:		emLiPositronTrkRange
*
*			Summary:	Tracks a positron, adjusting its range to account for
*						heterogeneous media.
*			Arguments:
*				double				waterRange			- This positron's range in water.
*				PHG_Position		startingPos			- The starting position of the positron.
*				PHG_Direction		direction			- The direction to project the positron.
*				double				positronEnergyMeV	- The photon's energy.
*				double				sigmaWater				- standard deviation of positron range in water.
*				double				*finalPosPtr		- The positron's annihilation point.
*				Boolean				*discard			- Flag for escaping.
*				LbFourByte			*sliceIndex			- Position index.
*				LbFourByte			*xIndex				- Position index.
*				LbFourByte			*yIndex				- Position index.
*
*			Function return: None.
*
*********************************************************************************/
void emLiPositronTrkRange(	double waterRange,
							PHG_Position startingPos,
							PHG_Direction direction,
							double positronEnergyMeV,
							double sigmaWater,
							PHG_Position *finalPosPtr,
							Boolean *discard,
							LbFourByte *sliceIndex,
							LbFourByte *xIndex,
							LbFourByte *yIndex)
							
{
	double			b1, b2;					/* Positron range constants from Palmer and Brownell */
	double			density;				/* Density of material in current cell */
	double			Rex;					/* extrapolated range for this material, energy */
	double			sigma;					/* standard deviation of range for this material, energy */
	double			distanceTrackedR;		/* Total distance tracked */
	double			initialZDistanceR;		/* Initial distance from current z to z boundary */
	double			initialXDistanceR;		/* Initial distance from current x to x boundary */
	double			initialYDistanceR;		/* Initial distance from current y to y boundary */
	double			distToNextXR;			/* Distance to next x crossing */
	double			distToNextYR;			/* Distance to next y crossing */
	double			distToNextZR;			/* Distance to next z crossing */
	double			generalDistToXR;		/* General distance to move relative to x axis */
	double			generalDistToYR;		/* General distance to move relative to y axis */
	double			generalDistToZR;		/* General distance to move relative to z axis */
	double			nextDistR;				/* Next incremental move */
	double			distanceR;				/* Distance about to move */
	double			waterRangeForMoveR;		/* Equivalent range-in-water used to make the next move */
	double			waterRangeUsedR;		/* Equivalent range-in-water used to make the previous moves */
	LbFourByte		newSliceIndexR;			/* New slice index for entering new slice */
	PHG_Position	newPositionR;			/* Used to update slice info */
	double			distToObjectSurfaceR;	/* distance to object surface */
	

	/* Clear counters */
	distanceTrackedR = 0.0;
	distanceR = 0.0;
	waterRangeUsedR = 0.0;
	waterRangeForMoveR = 0.0;
	
	/* Clear distance marker */
	nextDistR = 0.0;
	
	/* Clear discard flag */
	*discard = false;
	
	/* Get initial distance to voxel wall */
	SubObjGetInnerCellDistance(&startingPos, &direction, *sliceIndex,
		 *xIndex,  *yIndex, &initialXDistanceR, &initialYDistanceR, &initialZDistanceR);
	
	/* Avoid div by zero from direction cosines */
	{
		if ((direction.cosine_x <= 0.0000001) && 
				(direction.cosine_x >= -0.0000001)) {				
			direction.cosine_x = 
				(0.0000001 * ((direction.cosine_x < 0) ? -1 : 1));
		}
		if ((direction.cosine_y <= 0.0000001) && 
				(direction.cosine_y >= -0.0000001)) {				
			direction.cosine_y = 
				(0.0000001 * ((direction.cosine_y < 0) ? -1 : 1));
		}
		if ((direction.cosine_z <= 0.0000001) && 
				(direction.cosine_z >= -0.0000001)) {				
			direction.cosine_z = 
				(0.0000001 * ((direction.cosine_z < 0) ? -1 : 1));
		}
	}

	/* Calculate initial distances  (normalized to positive value) */
	{
		distToNextXR = initialXDistanceR/direction.cosine_x;
		distToNextYR = initialYDistanceR/direction.cosine_y;
		distToNextZR = initialZDistanceR/direction.cosine_z;
	}
	
	/* Calculate general distances */
	{
		generalDistToXR = SUBOBJGetSliceAttVoxelWidth(*sliceIndex)/fabs(direction.cosine_x);
		generalDistToYR = SUBOBJGetSliceAttVoxelHeight(*sliceIndex)/fabs(direction.cosine_y);
		generalDistToZR = SUBOBJGetSliceVoxelDepth(*sliceIndex)/fabs(direction.cosine_z);
	}
	
	do { /* Tracking Loop */

		/* Get shortest distance to next cell boundary */
		nextDistR = EMLIGetMinDeltaT(distToNextXR, distToNextYR, distToNextZR);
		
		/* Calculate distance to object's cylindrical surface */
		CylPosCalcDistanceToObjectSurface(&startingPos, &direction, &distToObjectSurfaceR);

		/* See if this will take us out of the object.
			If so, adjust distance to go onto object cylinder
		*/
		if (nextDistR >= distToObjectSurfaceR) {

			/* Set distance to put us onto cylinder */
			nextDistR = distToObjectSurfaceR;		
		}

		/* Get attenutation of current cell, fails if we are out of the object */
		if (!SubObjGetCellPosRangeConstants(*sliceIndex, *xIndex, *yIndex,
							&b1, &b2, &density)) {
			PhgAbort("Attempt to get density and positron range constants for cell not in subobject (PhoTrkCalcRange)",
				true);
		}
		
		/* compute standard deviation of range for this material, energy */
		Rex = (0.1 * b1 * positronEnergyMeV * positronEnergyMeV) / (b2 + positronEnergyMeV);
		sigma = Rex / (2*density);
		
		/* Compute distance to move through this cell */
		distanceR = nextDistR - distanceTrackedR;
		
		/* Store equivalent range in water for distance moved thus far */
		waterRangeUsedR = waterRangeForMoveR;

		/* Compute equivalent range in water used for this move */
		waterRangeForMoveR += (distanceR * (sigmaWater/sigma));
		
		/* See if this move will exhaust range in water available */
		if (waterRangeForMoveR >= waterRange) {

			/* Truncate distance to move */
			distanceR = distanceTrackedR + (waterRange-waterRangeUsedR)/(sigmaWater/sigma);

			/* See if this will take us out of the object.
				If so, bolt
			*/
			if (distanceR >= distToObjectSurfaceR) {

				*discard = true;
				break;	
			}
			
			/* Compute the new position */
			PhoTrkProject(&startingPos, &direction, distanceR, finalPosPtr);
			
			/* Break out of this loop, we are done tracking */
			break;
		}

		/* Increment distance tracked, then see if we need to keep going */
		distanceTrackedR = nextDistR;

		/* See if this will take us out of the object.
			If so, bolt
		*/
		if (distanceTrackedR >= distToObjectSurfaceR) {

			*discard = true;
			break;	
		}
		
		/* See which value we used */
		if (nextDistR == distToNextXR) {
			/* From here on out each x crossing is the same distance */
			distToNextXR += generalDistToXR;
			
			/* Increment our x index */
			*xIndex += ((direction.cosine_x >= 0) ? 1 : -1);
			
			/* See if we leave the object */
			if ((*xIndex < 0) || ((LbUsFourByte)*xIndex >= SubObjObject[*sliceIndex].attNumXBins)) { 
				*discard = true;
				break;
			}
			
		}
		else if (nextDistR == distToNextYR) {
			/* From here on out each y crossing is the same distance */
			distToNextYR += generalDistToYR;
			
			/* Increment our y index */
			/* Notice that Y goes top to bottom so it is reveresed from X */
			 *yIndex += ((direction.cosine_y >= 0) ? -1 : 1);
			
			
			/* See if we leave the object */
			if ((*yIndex < 0) || ((LbUsFourByte)*yIndex >= SubObjObject[*sliceIndex].attNumYBins)) { 
				*discard = true;
				break;
			}
			

		}
		else {
			
			/* We are going to a new slice so update the values */
			newSliceIndexR =  *sliceIndex +
				((direction.cosine_z >= 0) ? 1 : -1);

			/* See if we went out the end of the object */
			if ((newSliceIndexR < 0) || ((LbUsFourByte)newSliceIndexR == SubObjNumSlices)) {
				*discard = true;
				break;
			}

			/* If we didn't break, we entered a new slice of the object; recalculate slice variables */						
			{
				/* Compute the new position */
				PhoTrkProject(&startingPos, &direction, distanceTrackedR, &newPositionR);

				/* Update slice parameters */
				PhoTrkEnterNewSlice(newPositionR, direction,
			   		 *sliceIndex, distanceTrackedR,
					&distToNextXR, &distToNextYR, &distToNextZR,
					&generalDistToXR, &generalDistToYR, &generalDistToZR,
					sliceIndex, 
					xIndex,
					yIndex);
			}
		}
	} while (true);		
}


#undef EMIS_LIST
