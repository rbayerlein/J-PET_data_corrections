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
*			Module Name:		PhgMath.c
*			Revision Number:	1.6
*			Date last revised:	23 July 2013
*			Programmer:			Steven Vannoy
*			Date Originated:	Friday, September 18, 1992
*
*			Module Overview:	PHG Math routines.
*
*			References:			'PHG Physics/Math Functions' PHG design.
*
**********************************************************************************
*
*			Global functions defined:
*				PhgMathInit
*				PhgMathTerminate
*				PhgMathInitRNGFromSeed
*				PhgMathReadSeed
*				PhgMathWriteSeed
*				PhgMathGetRandomNumberOld
*				PhgMathGetRandomNumber
*				PhgMathGetDPRandomNumber
*				PhgMathGetTotalFreePaths
*				PhgMathRealNumAreEqual
*				PhgMathRealNumIsGreater
*				PhgMathSampleFromGauss
*				PhgMathSolveQuadratic
*
*			Global variables defined:		None
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
*			Revision date:		Sept 1, 2006
*
*			Revision description:	Added PhgMathRealNumIsGreater
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		24 May 2005
*
*			Revision description:	Added MT RNG; revised other RNG functions
*
*********************************************************************************/

#include	<stdio.h>
#include	<string.h>
#include	<time.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbError.h"
#include "LbMemory.h"
#include "LbParamFile.h"
#include "LbFile.h"
#include "LbInterface.h"
#include "LbHeader.h"
#include "MT19937.h"

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


/* Compiler switches controlling the RNG (Random Number Generator) in use; 
	be sure ONE and ONLY one is defined:         */
/*	Use local RNGs; first SimSET RNG:  */
/*#define		PHGMATH_USE_ORIG_RNG*/
/*	Use built-in Unix RNG:  */
/*#define		USE_DRAND48*/
/*	Use LC RNG (Linear Congruential RNG); previous SimSET RNG:  */
/* #define		PHGMATH_USE_LC_RNG */
/*	Use MT RNG (Mersenne Twister RNG):  */
#define		PHGMATH_USE_MT_RNG

/* LOCAL CONSTANTS */
#define						PHGMATH_POW_NORMALIZER		7		/*	Used to convert negative power 
																	integers to indexes in the power
																	table.
																*/
#define						PHGMATH_POW_MIN				-7		/*	Minimum power available for using
																	the table.
																*/
#define						PHGMATH_POW_MAX				6		/*	Maximum power available for using
																	the table.
																*/
#ifdef PHGMATH_USE_LC_RNG
#define						PHGMATH_RAND_MODULO		LBFOURBYTE_MAX
#define						PHGMATH_NUM_SEEDS		97
#endif

/* LOCAL TYPES */

/* LOCAL GLOBALS */
Boolean		phgMathIsInited = false;
double		phgMathPowTbl[] = {.0000001, .000001, .00001, .0001, .001, .01, .1, 1, 10, 100, 1000, 10000, 100000, 1000000};
double		phgMathGaussDeviate;				/* Holds pre-computed Gauss deviate */
Boolean		phgMathHaveGaussDeviate = false;	/* Flag for precomputed value */
#ifdef PHGMATH_USE_LC_RNG
LbFourByte	phgMathRandTable[PHGMATH_NUM_SEEDS];
LbFourByte	phgMathNextIntegerRand1;
LbFourByte	phgMathNextIntegerRand2;
#endif

#ifdef PHG_DEBUG
char		phgMathRandSeedPath[] = "phgmath.randseed";
FILE		*phgMathRandSeedFile = NULL;
Boolean		phgMathWriteToRandSeedFile = false;
Boolean		phgMathReadFromRandSeedFile = false;
Boolean		phgMathResetRandSeedFile = false;
#endif

#ifdef __MWERKS__
	double	phgMathRandSeed = 1;
#endif
#ifdef WINNT
	double	phgMathRandSeed = 1;
#endif

/* LOCAL MACROS */

/* PROTOTYPES */
#ifdef PHGMATH_USE_LC_RNG
Boolean phgMathLCReadState(FILE *phgMathRandSeedFile);
Boolean phgMathLCWriteState(FILE *phgMathRandSeedFile);
LbFourByte	phgMathRand0(LbFourByte	*randSeed);
#endif

#ifdef SUN_OS
double drand48(void);
double erand48(unsigned short *seed);
#endif
#ifdef HP_OS
double drand48(void);
double erand48(unsigned short *seed);
#endif
#ifdef DEC_ULTRIX
double drand48(void);
double erand48(unsigned short *seed);
#endif


/* FUNCTIONS */

/*********************************************************************************
*
*			Name:		PhgMathInit
*
*			Summary:	Initialize the math library.  Primarily, the random
*							number generator.
*
*			Arguments:
*				LbFourByte		*randSeed	- Seed for generator
*
*			Function return: True unless an error occurs.
*
*********************************************************************************/
Boolean PhgMathInit(LbFourByte *randSeed)
{
	LbFourByte		rngCount;			/* Count of defined RNGs */
	char			errString[81];		/* Error message */
	time_t			sysTime;			/* System time */
	
	
	do {
		
		/* Check that one and only one RNG is defined */
		rngCount = 0;
		#ifdef PHGMATH_USE_ORIG_RNG
			rngCount++;
		#endif
		#ifdef USE_DRAND48
			rngCount++;
		#endif
		#ifdef PHGMATH_USE_LC_RNG
			rngCount++;
		#endif
		#ifdef PHGMATH_USE_MT_RNG
			rngCount++;
		#endif
		if (rngCount != 1) {
			sprintf(errString, "There must be exactly one RNG, not %ld (PhgMathInit).", (long)rngCount);
			PhgAbort(errString, false);
		}
		
		
		#ifdef PHG_DEBUG
			phgMathRandSeedFile = NULL;
		#endif
		
		
		/* If random seed is not positive, set it from the system clock */
		if (*randSeed <= 0) {
			(void) time(&sysTime);
			*randSeed = ( LbFourByte )(sysTime & 0x0FFFFFFFF);	/* Keep lower 4 bytes */
			if (*randSeed < 0)
				*randSeed *= -1;
		}
		
		/* Initialize the RNG with the seed */
		PhgMathInitRNGFromSeed(*randSeed);
		
		#ifdef PHG_DEBUG
			/* If requested, overwrite above steps by initializing from random seed file */
			if (PHGDEBUG_ReadFromRandSeedFile()) {
				PhgMathReadSeed(-1);
			}
		#endif
		
		
		phgMathIsInited = true;
	} while (false);
	
	return (phgMathIsInited);
}

/*********************************************************************************
*
*			Name:		PhgMathTerminate
*
*			Summary:	Terminate the math library.
*
*			Arguments:	None.
*
*			Function return: None.
*
*********************************************************************************/
void PhgMathTerminate()	
{
	do {
		
		/* Don't do anything if the library is not initialized */
		if (phgMathIsInited == false)
			break;
		
		#ifdef PHG_DEBUG
			if (phgMathRandSeedFile != NULL) {
				fclose(phgMathRandSeedFile);
				phgMathRandSeedFile = NULL;
			}
		#endif
		
	} while (false);
}

/*********************************************************************************
*
*			Name:		PhgMathInitRNGFromSeed
*
*			Summary:	Initialize the RNG using the provided seed value.
*
*			Arguments:
*				LbFourByte		randSeed	- Seed for generator
*
*			Function return: None.
*
*********************************************************************************/
void PhgMathInitRNGFromSeed(LbFourByte randSeed)
{
	do {
		
		/* Original (now considered obsolete) SimSET RNG */
		#ifdef PHGMATH_USE_ORIG_RNG
		{
			/* There doesn't seem to be any universal way to use the seed */
			
			LbFourByte		localRandSeed = randSeed;
				/* To eliminate compiler warning for unused parameter */
		}
		#endif
		
		/* Built-in Unix RNG */
		#ifdef USE_DRAND48
			srand48(randSeed);
		#endif
		
		/* Linear Congruential RNG */
		#ifdef PHGMATH_USE_LC_RNG
		{
			LbFourByte		localSeed;	/* Local copy for recording */
			LbFourByte		i;			/* Index through array */
			
			
			/* Initialize our array of seeds */
			localSeed = randSeed;
			for (i=0; i<PHGMATH_NUM_SEEDS; i++) {
				phgMathRandTable[i] = phgMathRand0(&localSeed);
			}
			
			/* Save the seed value */
			phgMathNextIntegerRand1 = localSeed;
			
			/* Always initialize the secondary seed from the given seed */
			phgMathNextIntegerRand2 = (phgMathNextIntegerRand1+
					(PHGMATH_RAND_MODULO/2));
			if (phgMathNextIntegerRand2 < 0)
				phgMathNextIntegerRand2 *= -1;
			else if (phgMathNextIntegerRand2 == 0)
				phgMathNextIntegerRand2 = 1;
		}
		#endif
		
		/* Mersenne Twister RNG */
		#ifdef PHGMATH_USE_MT_RNG
			init_genrand((unsigned long)randSeed);
		#endif
		
	} while (false);
}

#ifdef PHGMATH_USE_LC_RNG
/*********************************************************************************
*
*		Name:		phgMathLCReadState
*
*		Summary:	Read the random number state vector and indices from the dump file.
*						This routine is used to restart the random sequence
*						from a specific point.
*
*		Arguments:
*				FILE	*phgMathRandSeedFile	- File containing the saved values.
*		
*		Function return: True if successful.
*
*********************************************************************************/
Boolean phgMathLCReadState(FILE *phgMathRandSeedFile)
{
	Boolean		result = true;		/* Function result */
	
	/* This function only has an effect in debug mode */
#ifdef PHG_DEBUG
	LbFourByte numRead;	/* Return value from reading data */
	
	do {
		/* Read in the state values */
		/* NOTE:  File must already be open */
		
		if ((numRead = fread((void *)&phgMathNextIntegerRand1,
				sizeof(phgMathNextIntegerRand1), 1, phgMathRandSeedFile)) != 1) {
			LbInPrintf("Unable to read random seed value 'phgMathNextIntegerRand1'");
			break;
		}
		
		if ((numRead = fread((void *)&phgMathNextIntegerRand2,
				sizeof(phgMathNextIntegerRand2), 1, phgMathRandSeedFile)) != 1) {
			LbInPrintf("Unable to read random seed value 'phgMathNextIntegerRand2'");
			break;
		}
		
		
		if ((numRead = fread((void *)phgMathRandTable,
				sizeof(LbFourByte), PHGMATH_NUM_SEEDS,
				phgMathRandSeedFile)) != PHGMATH_NUM_SEEDS) {
			LbInPrintf("Unable to read random seed value 'phgMathRandTable'");
			break;
		}
	} while (false);
#else
	FILE	*localFile = phgMathRandSeedFile;
		/* To eliminate compiler warning for unused parameter */
#endif
	
	return (result);
}
#endif

#ifdef PHGMATH_USE_LC_RNG
/*********************************************************************************
*
*		Name:		phgMathLCWriteState
*
*		Summary:	Write the random number state vector and indices to the dump file.
*						This routine is used to restart the random sequence
*						from a specific point.
*
*		Arguments:
*				FILE	*phgMathRandSeedFile	- File containing the saved values.
*				
*		Function return: True if successful.
*
*********************************************************************************/
Boolean phgMathLCWriteState(FILE *phgMathRandSeedFile)
{
	Boolean		result = true;		/* Function result */
	
	/* This function only has an effect in debug mode */
#ifdef PHG_DEBUG
	LbFourByte numWritten;	/* Return value from writing data */
	
	do {
		/* Write out the state values */
		/* NOTE:  File must already be open */
		
		if ((numWritten = fwrite((void *)&phgMathNextIntegerRand1,
				sizeof(phgMathNextIntegerRand1), 1, phgMathRandSeedFile)) != 1) {
			LbInPrintf("Unable to write random seed value 'phgMathNextIntegerRand1'");
			break;
		}
		
		if ((numWritten = fwrite((void *)&phgMathNextIntegerRand2,
				sizeof(phgMathNextIntegerRand2), 1, phgMathRandSeedFile)) != 1) {
			LbInPrintf("Unable to write random seed value 'phgMathNextIntegerRand2'");
			break;
		}
		
		
		if ((numWritten = fwrite((void *)phgMathRandTable,
				sizeof(LbFourByte), PHGMATH_NUM_SEEDS,
				phgMathRandSeedFile)) != PHGMATH_NUM_SEEDS) {
			LbInPrintf("Unable to write random seed value 'phgMathRandTable'");
			break;
		}
	} while (false);
#else
	FILE	*localFile = phgMathRandSeedFile;
		/* To eliminate compiler warning for unused parameter */
#endif
	
	return (result);
}
#endif

/*********************************************************************************
*
*			Name:		PhgMathReadSeed
*
*			Summary:	Read the random number seed data from the dump file.
*						This routine is used to restart the random sequence
*						from a specific point.
*			Arguments:
*				LbFourByte	resetCount	- Backup and write over previous values?
*				
*			Function return: None.
*
*********************************************************************************/
void PhgMathReadSeed(LbFourByte	resetCount)	
{
	/* This function only has an effect in debug mode */
#ifdef PHG_DEBUG
	do {
	
		/* If file is not already open, do it */
		if (phgMathRandSeedFile == NULL) {
			
			/* Open the file */
			if ((phgMathRandSeedFile = LbFlFileOpen(phgMathRandSeedPath, "r+b")) == NULL) {
				LbInPrintf("Unable to open random seed file.");
				break;
			}
		}
		else if (resetCount != 0) {

			if (resetCount < 0) {
				if (fseek(phgMathRandSeedFile, 0, SEEK_SET) != 0) {
					LbInPrintf("Unable to seek to beginning of random seed file.");
					break;
				}
			}
			else {
				
				/* Compute offset from current position */
				#ifdef PHGMATH_USE_LC_RNG
					resetCount = 
						-(resetCount * (sizeof(phgMathNextIntegerRand1) +
						sizeof(phgMathNextIntegerRand2) + 
						(sizeof(LbFourByte) * PHGMATH_NUM_SEEDS)));
				#endif
				#ifdef PHGMATH_USE_MT_RNG
					/* Be sure these are correct, since can't see into MT files */
					resetCount = 
						-(resetCount * (sizeof(int) +
						(sizeof(unsigned long) * 624)));
				#endif
				
				if (fseek(phgMathRandSeedFile, resetCount, SEEK_CUR) != 0) {
					LbInPrintf("Unable to seek to beginning of random seed file.");
					break;
				}
			}
			
		}
		
		/* Read in the seed values */
		#ifdef PHGMATH_USE_LC_RNG
			if (! phgMathLCReadState(phgMathRandSeedFile)) {
				/* Error message has already been printed */
				break;
			}
		#endif
		#ifdef PHGMATH_USE_MT_RNG
			if (! MTReadState(phgMathRandSeedFile)) {
				/* Error message has already been printed */
				break;
			}
		#endif
		
	} while (false);
#else
	if (resetCount) {};		/* Eliminate compiler warning about unused parameter */
#endif
}

/*********************************************************************************
*
*			Name:		PhgMathWriteSeed
*
*			Summary:	Write the random number seed data to the dump file.
*						This routine is used to restart the random sequence
*						from a specific point.
*			Arguments:
*				LbFourByte	resetCount	- Backup and write over previous values?
*				
*			Function return: None.
*
*********************************************************************************/
void PhgMathWriteSeed(LbFourByte resetCount)	
{
	/* This function only has an effect in debug mode */
#ifdef PHG_DEBUG
	do {
	
		/* If file is not already open, do it */
		if (phgMathRandSeedFile == NULL) {
			
			/* Open the file */
			if ((phgMathRandSeedFile = LbFlFileOpen(phgMathRandSeedPath, "w+b")) == NULL) {
				LbInPrintf("Unable to open random seed file.");
				break;
			}
		}
		else if (resetCount != 0) {

			if (resetCount < 0) {
				if (fseek(phgMathRandSeedFile, 0, SEEK_SET) != 0) {
					LbInPrintf("Unable to seek to beginning of random seed file.");
					break;
				}
			}
			else {
				
				/* Compute offset from current position */
				#ifdef PHGMATH_USE_LC_RNG
					resetCount = 
						-(resetCount * (sizeof(phgMathNextIntegerRand1) +
						sizeof(phgMathNextIntegerRand2) + 
						(sizeof(LbFourByte) * PHGMATH_NUM_SEEDS)));
				#endif
				#ifdef PHGMATH_USE_MT_RNG
					/* Be sure these are correct, since can't see into MT files */
					resetCount = 
						-(resetCount * (sizeof(int) +
						(sizeof(unsigned long) * 624)));
				#endif
				
				if (fseek(phgMathRandSeedFile, resetCount, SEEK_CUR) != 0) {
					LbInPrintf("Unable to seek to beginning of random seed file.");
					break;
				}
			}
			
		}
		
		/* Write out the seed values */
		#ifdef PHGMATH_USE_LC_RNG
			if (! phgMathLCWriteState(phgMathRandSeedFile)) {
				/* Error message has already been printed */
				break;
			}
		#endif
		#ifdef PHGMATH_USE_MT_RNG
			if (! MTWriteState(phgMathRandSeedFile)) {
				/* Error message has already been printed */
				break;
			}
		#endif
		
	} while (false);
#else
	if (resetCount) {};		/* Eliminate compiler warning about unused parameter */
#endif
}

/*********************************************************************************
*
*			Name:		PhgMathGetRandomNumberOld
*
*			Summary:	Return random number.
*						NOTE:  This function is now obsolete.
*
*			Arguments:	None.
*				
*			Function return: double, a random number between zero and one.
*
*********************************************************************************/
double PhgMathGetRandomNumberOld()	
{
	double	randomNumber;	/* The number we calculate */
	
	#ifdef PHG_DEBUG
		randomNumber = 2;
	#endif
	
	#ifdef __MWERKS__
	{
		double_t randomValue = 1.0;
		
		/* Reduce to between zero and 1.0 */
		randomNumber = (double) (randomx(&randomValue)-1)/(LONG_MAX-1);
	}
	#endif
	
	#ifdef AOS_VS
		randomNumber = dg_rand1();
	#endif
	
	#ifdef GEN_UNIX_OLD_RANDOM
		#ifdef PHG_DEBUG
			randomNumber = erand48(PhgRandSeedArray);
		#else 
			randomNumber = drand48();
		#endif
	#endif
	
	#ifdef PHG_DEBUG
	{
		char	errString[81];
		
		if (randomNumber == 2) {
			PhgAbort("You are running on an unsupported machine (PhgMathGetRandomNumber).",
				false);
		}
		
		if (randomNumber >= 1.0) {
			sprintf(errString, "Random number calculated out of range = %3.2e (PhgMathGetRandomNumber).", randomNumber);
			PhgAbort(errString, true);
		}
	}		
	#endif
	
	return (randomNumber);
}

#ifdef PHGMATH_USE_LC_RNG
/*********************************************************************************
*
*			Name:		phgMathRand0
*
*			Summary:	This is a simple random number generator that is used
*						to create an array of seed values for the "real" generator.
*			Arguments:
*				LbFourByte	randSeed	- Seed for generator
*
*				
*			Function return: A simple random number.
*
*********************************************************************************/
LbFourByte phgMathRand0(LbFourByte *randSeed)	
{
	
	#define MULTIPLIER	69621
	#define FACTOR1		30845
	#define FACTOR2		23902
	
	#ifdef PHG_DEBUG
		/* Verify seed value is valid, this should never happen in this local
			routine, because all calling routines should check it; hence this
			check is in debug mode only
		*/
		if ((*randSeed <= 0.0) || (*randSeed == PHGMATH_RAND_MODULO)){
			PhgAbort("Invalid seed passed to random number generator.", false);
		}
	#endif
	
	/* Compute the random seed */
	*randSeed = MULTIPLIER*(*randSeed % FACTOR1) - FACTOR2*(*randSeed/FACTOR1);
	if (*randSeed < 0)
		*randSeed += PHGMATH_RAND_MODULO;
	
	#undef MULTIPLIER
	#undef FACTOR1
	#undef FACTOR2

	return (*randSeed);
}
#endif

/*********************************************************************************
*
*			Name:		PhgMathGetRandomNumber
*
*			Summary:	This is the SimSET random number generator.
*
*			Arguments:	None.
*				
*			Function return: A random number between 0.0 and 1.0 (exclusive) = (0,1).
*
*********************************************************************************/
double PhgMathGetRandomNumber()
{
	double 		randReal;	/* The floating point random number created */
	
	
	#ifdef PHG_DEBUG
		if (phgMathIsInited == false) {
			PhgAbort("You have to initialize the math library before calling 'PhgMathGetRandomNumber'", false);
		}
		
		/* Allow debugger to reset from file */
		if (phgMathReadFromRandSeedFile == true) {
			PhgMathReadSeed(phgMathResetRandSeedFile);
		}
		
		/* Allow debugger to reset seed file */
		if (phgMathWriteToRandSeedFile == true) {
			PhgMathWriteSeed(phgMathResetRandSeedFile);
		}
	#endif
	
	
	/* Obtain the next random number */
	
	/* Original (now considered obsolete) SimSET RNG */
	#ifdef PHGMATH_USE_ORIG_RNG
		randReal = PhgMathGetRandomNumberOld();
	#endif
	
	/* Built-in Unix RNG */
	#ifdef USE_DRAND48
		randReal = drand48();
	#endif
	
	/* Linear Congruential RNG */
	#ifdef PHGMATH_USE_LC_RNG
	{
		LbFourByte	randIndex;	/* Index into table of random integers */
		
		
		/* Compute an index into the table of random integers */
		randIndex = phgMathNextIntegerRand2 % PHGMATH_NUM_SEEDS;
		
		/* Get a random integer from our precomputed table */
		phgMathNextIntegerRand2 = phgMathRandTable[randIndex];
		
		/* Replace the used random integer with a new one */
		phgMathRandTable[randIndex] = phgMathRand0(&phgMathNextIntegerRand1);
		
		/* Compute a random real number from the random integer */
		randReal = ((double)phgMathNextIntegerRand2) * (1.0 / PHGMATH_RAND_MODULO);
	}
	#endif
	
	/* Mersenne Twister RNG */
	#ifdef PHGMATH_USE_MT_RNG
	{
		LbUsFourByte	u;			/* An unsigned random integer */
		
		do {
			u = (LbUsFourByte) genrand_int32();
		} while (u == 0);
		randReal = ((double)u) * (1.0/4294967296.0);  /* divide by 2^32 */
	}
	#endif
	
	
	#ifdef PHG_DEBUG
		if ((randReal <= 0.0) || (randReal >= 1.0)) {
			PhgAbort("Invalid random number computed, value not between 0.0 & 1.0", false);
		}
	#endif
	
	return (randReal);
}

/*********************************************************************************
*
*			Name:		PhgMathGetDPRandomNumber
*
*			Summary:	Random number generator for double-precision numbers.
*
*			Arguments:	None.
*				
*			Function return: A random number between 0.0 and 1.0 (exclusive) = (0,1).
*
*********************************************************************************/
double PhgMathGetDPRandomNumber()
{
	double 		randReal;	/* The floating point random number created */
	
	
/* Only has different effect for Mersenne Twister RNG */
#ifdef PHGMATH_USE_MT_RNG
	#ifdef PHG_DEBUG
		if (phgMathIsInited == false) {
			PhgAbort("You have to initialize the math library before calling 'PhgMathGetDPRandomNumber'", false);
		}
		
		/* Allow debugger to reset from file */
		if (phgMathReadFromRandSeedFile == true) {
			PhgMathReadSeed(phgMathResetRandSeedFile);
		}
		
		/* Allow debugger to reset seed file */
		if (phgMathWriteToRandSeedFile == true) {
			PhgMathWriteSeed(phgMathResetRandSeedFile);
		}
	#endif
	
	
	/* Obtain the next random number */
	{
		LbUsFourByte	a, b;		/* Intermediate variables for 53-bit random number */
		
		a = ((LbUsFourByte)genrand_int32()) >> 5;
		if (a == 0) {
			/* Must now ensure that b is not also 0 */
			do {
				b = ((LbUsFourByte)genrand_int32()) >> 6;
			} while (b == 0);
		}
		else {
			b = ((LbUsFourByte)genrand_int32()) >> 6;
		}
		randReal = (a*67108864.0+b) * (1.0/9007199254740992.0);
			/*   ( a * 2^26  +  b ) / ( 2^53 )   */
	}
	
	
	#ifdef PHG_DEBUG
		if ((randReal <= 0.0) || (randReal >= 1.0)) {
			PhgAbort("Invalid random number computed, value not between 0.0 & 1.0", false);
		}
	#endif
	
#else
	/* Just use the regular RNG */
	randReal = PhgMathGetRandomNumber();
#endif
	
	
	return (randReal);
}

/*********************************************************************************
*
*			Name:		PhgMathGetTotalFreePaths
*
*			Summary:	Choose a random number of free path lengths for the photon 
*				        to travel, using an exponential distribution.
*			Arguments:
*				double	*totalFreePathLengthPtr.
*
*				
*			Function return: None.
*
*********************************************************************************/
void PhgMathGetTotalFreePaths(double *totalFreePathLengthPtr)	
{
	/* Calculate free path from exponential distribution */
	*totalFreePathLengthPtr = -(PHGMATH_Log(1- PhgMathGetRandomNumber()));
}

/*********************************************************************************
*
*			Name:		PhgMathRealNumAreEqual
*
*			Summary:	Compare to real numbers for equality
*						allowing for user specified tolerance.
*
*			Arguments:
*				double		r1		- First real value.
*					WARNING!: r1 must be > 0 if perMag != 0
*				double		r2		- Second real value.
*				LbOneByte	absMag	- Magnitude of allowable absolute difference.
*				LbOneByte	perMag	- Magnitude of allowable percentage difference.
*				double		*absDif	- Absolute difference.
*				double		*perDif	- Percentage difference.
*				
*
*				
*			Function return: True if values are equal within tolerance,
*				i.e., if (|r1 - r2| <= 10^absMag) AND 
*				(|r1 - r2|/r1 <= 10^perMag), with the second
*				comparison performed only when perMag != 0.
*
*********************************************************************************/
Boolean PhgMathRealNumAreEqual(double r1, double r2, LbOneByte absMag,
			LbOneByte perMag, double *absDifPtr, double *perDifPtr)	
{
	Boolean	realsAreEqual = false;	/* Equality flag */
	double	absDifference;			/* Absolute difference */
	double	absTolerance;			/* Tolerance for absolute difference */
	double	perDifference;			/* Percentage difference */
	double	perTolerance;			/* Tolerance for percentage difference */
	
	do {	/* Process Loop */
		/* Assume no difference */
		absDifference = 0.0;
		perDifference = 0.0;
		
		/* Calculate abslute tolerance */
		#ifdef PHG_DEBUG
			if ((absMag < PHGMATH_POW_MIN) || (absMag > PHGMATH_POW_MAX)) {
				PhgAbort("You have called PhgMathRealNumAreEqual with an absolute magnitude out of range.",
					true);
			}
		#endif
		
		absTolerance = phgMathPowTbl[PHGMATH_POW_NORMALIZER+absMag];

		/* Calculate percentage tolerance */
		#ifdef PHG_DEBUG
			if ((perMag < PHGMATH_POW_MIN) || (perMag > PHGMATH_POW_MAX)) {
				PhgAbort("You have called PhgMathRealNumAreEqual with an percentage magnitude out of range.",
					true);
			}
		#endif
		perTolerance = phgMathPowTbl[PHGMATH_POW_NORMALIZER+perMag];

		/* See if absolute difference is outside tolerance */
		if ((absDifference = fabs(r1 - r2)) > absTolerance) {
			
			/* See if percentage difference is outside tolerance */
			if (perMag != 0) {
				if ((perDifference = fabs(r1 - r2)/r1) > perTolerance)
					break;
			}
			else {
				/* Only one tolerance given, and its out of bounds */
				break;
			}
		}
		
		/* If we made it to here, we are equal */
		realsAreEqual = true;
	} while (false);
	
	/* Save results if they want them */
	if (absDifPtr != 0)
		*absDifPtr = absDifference;
	
	if (perDifPtr != 0)
		*perDifPtr = perDifference;
		
	return (realsAreEqual);
}

/*********************************************************************************
*
*			Name:		PhgMathRealNumIsGreater
*
*			Summary:	Compare two real numbers to see if the first is strictly
*						greater than the second even after allowing for precision
*						errors.  The user specifies a tolerance to account for
*						precision.
*
*			Arguments:
*				double		r1		- First real value.
*				double		r2		- Second real value.
*				LbOneByte	Mag	- Magnitude of allowable  tolerance.
*				double		*Dif	-  difference (computed here).
*				
*
*				
*			Function return: True if r1 is greater than (r2 + tolerance),
*				i.e., if (r1 - r2 > 10^Mag).
*
*********************************************************************************/
Boolean PhgMathRealNumIsGreater(double r1, 			/* First real value */
								double r2, 			/* Second real value */
								LbOneByte Mag,	/* Magnitude of allowable absolute tolerance */
								double *DifPtr	/* Absolute difference (computed here) */
								)	
{
	Boolean	realIsGreater = true;	/* Flag indicating if inequality is satisfied */
	double	Difference;			/* Absolute difference */
	double	Tolerance;			/* Tolerance for absolute difference */
	
	do {	/* Process Loop */

		/* Assume no difference */
		Difference = 0.0;
		
		/* Calculate abslute tolerance */
		#ifdef PHG_DEBUG
			if ((Mag < PHGMATH_POW_MIN) || (Mag > PHGMATH_POW_MAX)) {
				PhgAbort("You have called PhgMathRealNumIsGreater with a magnitude out of range.",
					true);
			}
		#endif
		
		Tolerance = phgMathPowTbl[PHGMATH_POW_NORMALIZER+Mag];

		/* See if absolute difference is > tolerance */
		if ( (Difference = (r1 - r2)) > Tolerance ) {
			
			/* r1 is definitely greater--realIsGreater was preset true */
			break;

		}
		
		/* If we made it to here, we are equal */
		realIsGreater = false;

	} while (false);
	
	/* Save result if they want them */
	if (DifPtr != 0)
		*DifPtr = Difference;
	
	return (realIsGreater);
}

/*********************************************************************************
*
*			Name:		PhgMathSampleFromGauss
*
*			Summary:	Sample from a gaussian distribution with given mean and
*						standard deviation. NOTE: The core of this routine was
*						taken from "Numerical Recipes".
*			Arguments:
*				double		mean		- Mean of the Gauss distribution.
*				double		standDev	- Standard deviation for distribution.
*
*			Function return: Selection from distribution.
*
*********************************************************************************/
double PhgMathSampleFromGauss(double mean,
		double standDev)	
{
	double	sampleValue;	/* Our sampled value */
	double	fac, r, v1, v2;
	
	/* Find converging value if necessary */
	if (phgMathHaveGaussDeviate == false) {
		do {
			/* Pick two uniform numbers in the square extending from -1 to 1 */
			v1 = 2.0 * PhgMathGetRandomNumber() - 1.0;
			v2 = 2.0 * PhgMathGetRandomNumber() - 1.0;
			
			/* Compute radius of circle */
			r = (v1 * v1) + (v2 * v2);
			
		/* Terminate loop if radius is within unit circle */
		} while ((r >= 1.0) || (r == 0.0));
		
		/* Make the Box-Muller transform to get two normal deviates */
		fac = PHGMATH_SquareRoot(-2.0 * PHGMATH_Log(r)/r);
		
		/* Compute first deviate */
		sampleValue = v1 * fac;
		
		/* Store second deviate for next time */
		phgMathGaussDeviate = v2 * fac;
		
		/* Set flag indicating we have a value */
		phgMathHaveGaussDeviate = true;
	}
	else {
		/* Use pre-computed value */
		sampleValue = phgMathGaussDeviate;
		
		/* Clear flag */
		phgMathHaveGaussDeviate = false;
	}
	
	/* Now adapt deviate to our distribution */
	sampleValue = standDev * sampleValue + mean ;
		
	return (sampleValue);
}

/*********************************************************************************
*
*			Name:		PhgMathSolveQuadratic
*
*			Summary:	Compute real roots for quadratic formula.
*			Arguments:
*				double		a			- 
*				double		b
*				double		c
*				double		*minRoot	- Smaller of two roots.
*				double		*maxRoot	- Larger of two roots.
*
*			Function return: Number of roots found.
*
*********************************************************************************/
LbUsFourByte PhgMathSolveQuadratic(double a, double b, double c,
				double *minRoot, double *maxRoot)	
{
	LbUsFourByte numRealRoots = 0;
	double	q;
	double	fourAC;
	double	bSquMinus4AC;
	double	root1, root2;
	
	do /* Process Loop */
	{
		/* Compute 4*a*c */
		fourAC = 4*a*c;
		
		/* Check for single root/no real roots */
		if ((bSquMinus4AC = PHGMATH_Square(b)-fourAC) <= 0.0) {
			
			if ((bSquMinus4AC = PHGMATH_Square(b)-fourAC) < 0.0) {
				/* no real roots */
				break;
			} else if (a == 0) {
				/* degenerative condition, a=b=0, we return no roots */
				break;
			} else {
				/* return single root */
				*minRoot = -b/(2*a);
				numRealRoots = 1;
				break;
			}
		}
				
				
		
		/* Compute temporary 'q' */
		if (b > 0.0) {
			q = -0.5*(b + PHGMATH_SquareRoot(bSquMinus4AC));
		}
		else {
			q = -0.5*(b - PHGMATH_SquareRoot(bSquMinus4AC));
		}
		
		#ifdef PHG_DEBUG
			if (q == 0.0) {
				ErAbort("Got q == 0.0 in PhgMathSolveQuadratic");
			}
		#endif	
		
		if (a == 0.0) {
			*minRoot = c/q;
			numRealRoots = 1;
			break;
		}
		else {
			root1 = q/a;
			root2 = c/q;
			*minRoot = PHGMATH_Min(root1, root2);
			*maxRoot = PHGMATH_Max(root1, root2);
			numRealRoots =2;
			break;
		}
	} while (false);
	
	
	return (numRealRoots);
}
