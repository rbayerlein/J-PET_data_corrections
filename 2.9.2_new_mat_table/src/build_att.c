/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1990-2012 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		build_att.c
*			Revision Number:	1.2
*			Date last revised:	14 December 2012
*			Programmer:			Steven Vannoy
*			Date Originated:	5 January 1993
*
*			Module Overview:	Creates attenuation tables.
*
*			References:			Picard Y, Thompson C.J., Marrett S., "Improving 
*								the Precision and Accuracy of Monte Carlo Simulation 
*								in Positron Emission Tomography", IEEE Trans. on Nuc. 
*								Sci., Vol. 39, NO. 4, 1992
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
*********************************************************************************/

#include <stdio.h>
#include <string.h>

#include "SystemDependent.h"

#include "LbTypes.h"
#include "LbMacros.h"
#include "LbEnvironment.h"
#include "LbFile.h"
#include "LbError.h"
#include "LbInterface.h"
#include "LbMemory.h"

#define	DEBUGGING

/* LOCAL CONSTANTS */
#define			MIN_ENERGY			5.0
#define			PI 	 				3.14159265

#define			AIR_INDEX					0
#define			WATER_INDEX					1
#define			BLOOD_INDEX					2
#define			BONE_INDEX					3
#define			BRAIN_INDEX					4
#define			BRAIN_STEM_INDEX			5
#define			CALCIUM_INDEX				6
#define			CEREBELLUM_INDEX			7
#define			CEREBRUM_INDEX				8
#define			CERBRO_SPIN_FLUID_INDEX		9
#define			FAT_INDEX					10
#define			HEART_INDEX					11
#define			LUNG_INFLATED_INDEX			12
#define			MUSCLE_INDEX				13
#define			SKIN_INDEX					14
#define			LEAD_INDEX			 		15
#define			ALUM_INDEX					16
#define			LUCITE_INDEX				17
#define			TUNGST_INDEX				18
#define			POLYSTY_INDEX				19
#define			SODIUM_ID_INDEX				20

#define			AIR_DE						(3.622 * pow(10.0, 20.0))
#define			WATER_DE					(3.343 * pow(10.0, 23.0))
#define			BLOOD_DE					(3.513 * pow(10.0, 23.0))
#define			BONE_DE						(5.267 * pow(10.0, 23.0))
#define			BRAIN_DE					(3.438 * pow(10.0, 23.0))
#define			BRAIN_STEM_DE				(3.500 * pow(10.0, 23.0))
#define			CALCIUM_DE					(4.628 * pow(10.0, 23.0))
#define			CEREBELLUM_DE				(3.428 * pow(10.0, 23.0))
#define			CEREBRUM_DE					(3.433 * pow(10.0, 23.0))
#define			CERBRO_SPIN_FLUID_DE		(3.420 * pow(10.0, 23.0))
#define			FAT_DE						(3.061 * pow(10.0, 23.0))
#define			HEART_DE					(3.415 * pow(10.0, 23.0))
#define			LUNG_INFLATED_DE			(8.607 * pow(10.0, 22.0))
#define			MUSCLE_DE					(3.445 * pow(10.0, 23.0))
#define			SKIN_DE						(3.639 * pow(10.0, 23.0))
#define			LEAD_DE			 			(2.704 * pow(10.0, 24.0))
#define			ALUM_DE						(7.840 * pow(10.0, 23.0))
#define			LUCITE_DE					(3.833 * pow(10.0, 23.0))
#define			TUNGST_DE					(4.690 * pow(10.0, 24.0))
#define			POLYSTY_DE					(3.380 * pow(10.0, 23.0))
#define			SODIUM_ID_DE				(9.429 * pow(10.0, 23.0))

#define			AIR_A						(7.680)
#define			WATER_A						(6.869 * pow(10.0, 3.0))
#define			BLOOD_A						(7.110 * pow(10.0, 3.0))
#define			BONE_A						(3.989 * pow(10.0, 4.0))
#define			BRAIN_A						(7.031 * pow(10.0, 3.0))
#define			BRAIN_STEM_A				(7.727 * pow(10.0, 3.0))
#define			CALCIUM_A					(1.762 * pow(10.0, 5.0))
#define			CEREBELLUM_A				(7.044 * pow(10.0, 3.0))
#define			CEREBRUM_A					(7.002 * pow(10.0, 3.0))
#define			CERBRO_SPIN_FLUID_A			(7.415 * pow(10.0, 3.0))
#define			FAT_A						(3.725 * pow(10.0, 3.0))
#define			HEART_A						(6.658 * pow(10.0, 3.0))
#define			LUNG_INFLATED_A				(1.779 * pow(10.0, 3.0))
#define			MUSCLE_A					(7.070 * pow(10.0, 3.0))
#define			SKIN_A						(6.618 * pow(10.0, 3.0))
#define			LEAD_A			 			(8.720 * pow(10.0, 6.0))
#define			ALUM_A						(9.211 * pow(10.0, 4.0))
#define			LUCITE_A					(4.972 * pow(10.0, 3.0))
#define			TUNGST_A					(1.554 * pow(10.0, 7.0))
#define			POLYSTY_A					(2.805 * pow(10.0, 3.0))
#define			SODIUM_ID_A					(1.999 * pow(10.0, 6.0))

/* Define low energy "A" values */
#define			LEAD_A_LOW			 		(1.722 * pow(10.0, 6.0))
#define			TUNGST_A_LOW				(4.202 * pow(10.0, 6.0))

/* Define "B" Values */
#define			AIR_B						3.169
#define			WATER_B						3.192
#define			BLOOD_B						3.178
#define			BONE_B						3.075
#define			BRAIN_B						3.171
#define			BRAIN_STEM_B				3.169
#define			CALCIUM_B					3.054
#define			CEREBELLUM_B				3.172
#define			CEREBRUM_B					3.171
#define			CERBRO_SPIN_FLUID_B			3.183
#define			FAT_B						3.202
#define			HEART_B						3.178
#define			LUNG_INFLATED_B				3.173
#define			MUSCLE_B					3.175
#define			SKIN_B						3.182
#define			LEAD_B			 			2.588
#define			ALUM_B						3.131
#define			LUCITE_B					3.202
#define			TUNGST_B					2.648 
#define			POLYSTY_B					3.218  
#define			SODIUM_ID_B					2.791

/* Define low energy "B" values */
#define			LEAD_B_LOW		 			2.549
#define			TUNGST_B_LOW				2.716

/* Define low energy thresholds */
#define			LEAD_LOW_ENERGY				90
#define			TUNGST_LOW_ENERGY			70

#define			NUM_TYPES			21								/* Number of types defined above */
#define			NUM_ATT_VALUES		996								/* From MIN_ENERGY to 1000 kev */


double			densityTbl[NUM_TYPES];								/* Density values */
double			fitParamATbl[NUM_TYPES];							/* Fitting parameter A */
double			fitParamBTbl[NUM_TYPES];							/* Fitting parameter B */
double			*linearAttenuationTbl;
double			*comptonProbability;

char			*tissueNames[] = {"air", "water", "blood", "bone", "brain",
						"brain_stem", "calcium", "cerebellum", "cerebrum",
						"cerebro spinal fluid", "fat", "heart", 
						"lung inflated", "muscle", "skin", "lead",
						"alum", "lucite", "tungsten", "polysterene",
						"sodium_id"};
/* LOCAL GLOBALS */
/* PROTOTYPES */
Boolean	BuildAtt(int argc, char *argv[]);

/* FUNCTIONS */
/**********************
*	BuildAtt
*
*	Purpose:	Execute the program.
*
*	Result:	True unless an error occurs.
***********************/
Boolean BuildAtt(int argc, char *argv[])
{
	Boolean			okay = false;		/* Process Loop */
	double			alpha;
	double			coeff1;
	double			sigma;
	double			tau;
	double			temp1, temp2, temp3, temp4, temp5;
	FILE			*attFile;
	LbUsFourByte	typeIndex;
	LbUsFourByte	valueIndex;
	
	
	/* Remove compiler warnings */
	if (argc) {};
	if (argv) {};
	
	do { /* Process Loop */
		
		/* Create the output file */
		if ((attFile = LbFlFileOpen("phg.attenuation", "w")) == 0) {
			ErAbort("Unable to create/open output file. ");
		}

		/* Allocate memory for linear attenuation values */
		if ((linearAttenuationTbl = (double *) LbMmAlloc(sizeof(double) * NUM_TYPES * NUM_ATT_VALUES)) == 0) {
			ErStGeneric("Unable to allocate memory for linear attenuation table.");
			break;
		}
		
		/* Allocate memory for probability of compton scatter */
		if ((comptonProbability = (double *) LbMmAlloc(sizeof(double) * NUM_TYPES * NUM_ATT_VALUES)) == 0) {
			ErStGeneric("Unable to allocate memory for compton probability table.");
			break;
		}
		
		/* Initialize the density table */
		{
				densityTbl[0] = AIR_DE;
				densityTbl[1] = WATER_DE;
				densityTbl[2] = BLOOD_DE;
				densityTbl[3] = BONE_DE;
				densityTbl[4] = BRAIN_DE;
				densityTbl[5] = BRAIN_STEM_DE;
				densityTbl[6] = CALCIUM_DE;
				densityTbl[7] = CEREBELLUM_DE;
				densityTbl[8] = CEREBRUM_DE;
				densityTbl[9] = CERBRO_SPIN_FLUID_DE;
				densityTbl[10] = FAT_DE;
				densityTbl[11] = HEART_DE;
				densityTbl[12] = LUNG_INFLATED_DE;
				densityTbl[13] = MUSCLE_DE;
				densityTbl[14] = SKIN_DE;
				densityTbl[15] = LEAD_DE;
				densityTbl[16] = ALUM_DE;
				densityTbl[17] = LUCITE_DE;
				densityTbl[18] = TUNGST_DE;
				densityTbl[19] = POLYSTY_DE;				
				densityTbl[20] = SODIUM_ID_DE;

				fitParamATbl[0] = AIR_A;
				fitParamATbl[1] = WATER_A;
				fitParamATbl[2] = BLOOD_A;
				fitParamATbl[3] = BONE_A;
				fitParamATbl[4] = BRAIN_A;
				fitParamATbl[5] = BRAIN_STEM_A;
				fitParamATbl[6] = CALCIUM_A;
				fitParamATbl[7] = CEREBELLUM_A;
				fitParamATbl[8] = CEREBRUM_A;
				fitParamATbl[9] = CERBRO_SPIN_FLUID_A;
				fitParamATbl[10] = FAT_A;
				fitParamATbl[11] = HEART_A;
				fitParamATbl[12] = LUNG_INFLATED_A;
				fitParamATbl[13] = MUSCLE_A;
				fitParamATbl[14] = SKIN_A;
				fitParamATbl[15] = LEAD_A_LOW;
				fitParamATbl[16] = ALUM_A;
				fitParamATbl[17] = LUCITE_A;
				fitParamATbl[18] = TUNGST_A_LOW;
				fitParamATbl[19] = POLYSTY_A;
				fitParamATbl[20] = SODIUM_ID_A;

				fitParamBTbl[0] = AIR_B;
				fitParamBTbl[1] = WATER_B;
				fitParamBTbl[2] = BLOOD_B;
				fitParamBTbl[3] = BONE_B;
				fitParamBTbl[4] = BRAIN_B;
				fitParamBTbl[5] = BRAIN_STEM_B;
				fitParamBTbl[6] = CALCIUM_B;
				fitParamBTbl[7] = CEREBELLUM_B;
				fitParamBTbl[8] = CEREBRUM_B;
				fitParamBTbl[9] = CERBRO_SPIN_FLUID_B;
				fitParamBTbl[10] = FAT_B;
				fitParamBTbl[11] = HEART_B;
				fitParamBTbl[12] = LUNG_INFLATED_B;
				fitParamBTbl[13] = MUSCLE_B;
				fitParamBTbl[14] = SKIN_B;
				fitParamBTbl[15] = LEAD_B_LOW;
				fitParamBTbl[16] = ALUM_B;
				fitParamBTbl[17] = LUCITE_B;
				fitParamBTbl[18] = TUNGST_B_LOW;
				fitParamBTbl[19] = POLYSTY_B;
				fitParamBTbl[20] = SODIUM_ID_B;

		}
		
		
		/* Compute constant coefficient */
		coeff1 = 2 * PI * (pow(2.817938 * pow(10.0, -13.0), 2.0));
		
		/* Compute attenuation and compton scatter probability for all types */
		for (valueIndex = 0; valueIndex < NUM_ATT_VALUES; valueIndex++){
		
			/* Check current energy level, and set new "A" and "B" values if
				necessary
			*/
			{
				/* Check low energy for lead */
				if ((valueIndex + MIN_ENERGY) >= LEAD_LOW_ENERGY) {
					fitParamATbl[LEAD_INDEX] = LEAD_A;
					fitParamBTbl[LEAD_INDEX] = LEAD_B;
			}
				
				/* Check low energy for tungsten */
				if ((valueIndex + MIN_ENERGY) >= TUNGST_LOW_ENERGY) {
					fitParamBTbl[TUNGST_INDEX] = TUNGST_B;
					fitParamATbl[TUNGST_INDEX] = TUNGST_A;
				}
			}
			
			/* Compute alpha, reduced photon energy */
			alpha = (valueIndex+MIN_ENERGY)/511.0;

			/* Compute sigma */
			{
				temp1 = (1+alpha)/(alpha*alpha);
				temp2 = (2*(1+alpha))/(1 + 2*alpha);
				temp3 = log(1+2*alpha)/alpha;
				temp4 = log(1+2*alpha)/(2*alpha);
				temp5 = (1+3*alpha)/((1+2*alpha)*(1+2*alpha));
				
				sigma = coeff1*((temp1*(temp2-temp3)) + (temp4-temp5));
			}	
			for (typeIndex = 0; typeIndex < NUM_TYPES; typeIndex++) {
			
				/* Compute tau for this tissue */
				tau = fitParamATbl[typeIndex] * pow((double)(valueIndex+MIN_ENERGY), -fitParamBTbl[typeIndex]);
				
				/* Compute linear attenuation for this type */
				linearAttenuationTbl[valueIndex + (typeIndex*NUM_ATT_VALUES)] =
					 densityTbl[typeIndex]*sigma+tau;
				
				
				/* Compute compton probability */
				comptonProbability[valueIndex + (typeIndex*NUM_ATT_VALUES)] = (densityTbl[typeIndex]*sigma)/
					linearAttenuationTbl[valueIndex + (typeIndex*NUM_ATT_VALUES)];
			}
		}
		
		/* Put the number of types at the beginning of the file (1 is added for perfect absorber) */
		fprintf(attFile, "%d\n", NUM_TYPES+1);

		/* Print out attenuation and compton scatter probability for all types */
		for (typeIndex = 0; typeIndex < NUM_TYPES; typeIndex++) {
			fprintf(attFile, "%s\n", tissueNames[typeIndex]);
			
			for (valueIndex = 0; valueIndex < NUM_ATT_VALUES; valueIndex++){

/* NOTE: We are currently writing the compton probability twice just create "dummy" tables,
This will eventually need to be replaced with a cumulative compton, coherent probability
*/
				fprintf(attFile, "%4.5lf %3.5lf %3.5lf %3.5lf\n",
					(double) valueIndex+MIN_ENERGY, linearAttenuationTbl[valueIndex + (typeIndex*NUM_ATT_VALUES)],
					comptonProbability[valueIndex + (typeIndex*NUM_ATT_VALUES)],
					comptonProbability[valueIndex + (typeIndex*NUM_ATT_VALUES)]);
			}
		}
		
		/* Now print out the values for our perfect absorber */
		fprintf(attFile, "perfect_absorber\n");
		for (valueIndex = 0; valueIndex < NUM_ATT_VALUES; valueIndex++){

			fprintf(attFile, "%4.5lf %3.5lf %3.5lf %3.5lf\n",
				(double) valueIndex+MIN_ENERGY, 1000000.0, 0.0, 0.0);
		}
		
		/* Close the file */
		fclose(attFile);
		okay = true;
	} while (false);
	
	/* Free the memory */
	if (linearAttenuationTbl != 0)
		LbMmFree((void **)&linearAttenuationTbl);
		
	if (comptonProbability != 0)
		LbMmFree((void **)&comptonProbability);
		
	/* Quit the program */
	return (okay);
}
