/*********************************************************************************

*                                                                                *

*                       Source code developed by the                             *

*           Imaging Research Laboratory - University of Washington               *

*               (C) Copyright 2006-2012 Department of Radiology                  *

*                           University of Washington                             *

*                              All Rights Reserved                               *

*                                                                                *

*********************************************************************************/

 

/*********************************************************************************

*

*			Module Name:		addrandUsrSample.c

*			Revision Number:	1.0

*			Date last revised:	1 July 2012

*			Programmer:			Steven Gillispie

*			Date Originated:	1 February 2012

*

*			Module Overview:	This file is for Example 1 for the user functions.

*								It is a #include file for addrandUsr.c.  The

*								example is very simple - a pair of counters are

*								initialized, incremented, and printed out -

*								but each of the functions also gives some

*								information about where in addrandoms.c the

*								function is called and what other functions

*								could be inserted at that point.

*

*			References:			

*

**********************************************************************************

*

*			Global functions defined:

*

*			Global variables defined:

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

*			Revision date:		1 July 2012

*

*			Revision description:	This file started as a revision of S. 

*								Gillispie's addrandUsr.c version 2.2.

*

*********************************************************************************/





/* Global variables for addrandUsrSample */

/*

##

##

	Local globals are declared within the addrandUsrSampleVars struct below.  This will prevent any

		problems with naming and scope.

##

##

*/

struct {

	/* Your fields of choice go here */

	LbEightByte numModDecays1Calls;		/* number of calls to the first ModDecays function */

	LbEightByte numModDecays2Calls;		/* number of calls to the second ModDecays function */

} addrandUsrSampleVars;





/*********************************************************************************

*

*		Name:			addrandUsrSampleInitialize

*

*		Summary:	Initialize the counter variables.

*

*			This function is assigned to the function pointer

*			addrandUsrInitializeFPtr, which is called near the

*			end of adrandInitialize.  This happens only once at the

*			beginning of a run.

*			Other functions one might want to point to at this point

*			could include:  a function to read in parameters from

*			a user-defined parameter file (see user function example 3);

*			a function to set the addrandom user function pointers

*			dynamically (e.g., there could be several options for one

*			of the function pointers, with the one to be used

*			decided on the basis of the detector parameters); or

*			initialization of user variables, as is done in this

*			sample.

*			

* 

*

*		Arguments:

*			DetRunTimeParamsTy	*detParams

*						This is a ptr to the detector parameters. It

*						is assumed that you will treat this as 

*						a read-only structure. Modifications should 

*						only be made after a careful study of 

*						the parameters involved.

*

*		Function return: None.

*

*********************************************************************************/

void addrandUsrSampleInitialize(DetRunTimeParamsTy	*detParams)

{

	/* Do your initialization here */

	

	DetRunTimeParamsTy		*myDetParams;

	

	myDetParams = detParams;	/* Removes unused parameter compiler warning */

	

	/* initialize counters for the number of calls to the ModDecay functions */

	addrandUsrSampleVars.numModDecays1Calls = 0;

	addrandUsrSampleVars.numModDecays1Calls = 0;

	

}



/*********************************************************************************

*

*		Name:		addrandUsrSampleModDecays1

*

*		Summary:	Increment a counter - keeps track of how many times

*				this function is called.

*

*				This function is assigned to the function pointer

*				addrandUsrModDecays1FPtr, which is called before the 

*				photons/decays are time-windowed.  All the decays within 

*				the coincidence resolution time of each other (or 

*				actually longer: the window is extended each time another

*				photon arrives less than the res. time after the previous

*				photon) are passed in the 

*				timeWindowDetectionsTy structure.  The photons and decays

*				may be altered/deleted as desired, e.g. to account for

*				deadtime, apply different triples handling than SimSET's

*				(currently all triples are deleted, even if they occur

*				in separate rings of the tomograph).  Make sure you

*				understand the data structure first!

*					

*

*		Arguments:

*			timeWindowDetectionsTy	*timeWindowDetections - 

*					All the decays in the current time window

*

*		Function return: None

*

*********************************************************************************/

void addrandUsrSampleModDecays1(timeWindowDetectionsTy	*timeWindowDetectionsPtr,

							PhoHFileHkTy			*adrandRandomsFileHk, 

							LbUsFourByte			*numDecaysWritten, 

							LbUsFourByte			*numDecaysUnchanged, 

							LbUsFourByte			*numDecaysRandom, 

							LbUsFourByte			*numLostCorrectWindow, 

							LbUsFourByte			*numDecaysLostTriples)

{

	timeWindowDetectionsTy		*myTimeWindowDetectionsPtr;

	

	myTimeWindowDetectionsPtr = timeWindowDetectionsPtr;	/* Removes unused parameter compiler warning */

	

	

	/* Make your alterations here */

	do {

		

		/* increment the number of calls to this function */

		addrandUsrSampleVars.numModDecays1Calls += 1;

		

	} while (false);

	

	return;

}



/*********************************************************************************

*

*			Name:			addrandUsrSampleModDecays2

*

*		Summary:	Increment a counter - keeps track of how many times

*				this function is called.

*

*				This function is assigned to the function pointer

*				addrandUsrModDecays2FPtr, which is called  after the 

*				photons/decays are time-windowed.  Only the decays that are 

*				going to be written out will be passed to this routine - many

*				of the decays/time windows passed to addrandUsrSampleModDecays1

*				will not generate a call to this function.  This 

*				function gives the user a chance to change/reject the 

*				coincidences that are getting written out.

*					

*		Arguments:

*			timeWindowDetectionsTy	*timeWindowDetections - 

*					All the decays in the current time window

*

*		Function return: Boolean

*			True to accept the coincidence, False to reject it.

*

*********************************************************************************/

Boolean addrandUsrSampleModDecays2(timeWindowDetectionsTy	*timeWindowDetectionsPtr,

							PhoHFileHkTy			*adrandRandomsFileHk, 

							LbUsFourByte			*numDecaysWritten, 

							LbUsFourByte			*numDecaysUnchanged, 

							LbUsFourByte			*numDecaysRandom, 

							LbUsFourByte			*numLostCorrectWindow, 

							LbUsFourByte			*numDecaysLostTriples)



{	

	Boolean	acceptCoincidence = false;	/* Should coincidence be added to history list? */

	timeWindowDetectionsTy			*myTimeWindowDetectionsPtr;

	

	myTimeWindowDetectionsPtr = timeWindowDetectionsPtr;	/* Removes unused parameter compiler warning */

	

	

	/* Make your alterations here */

	do {

		

		/* increment the number of calls to this function */

		addrandUsrSampleVars.numModDecays2Calls += 1;

		

		/* Set acceptance to true if we made it to here */

		acceptCoincidence = true;

	

	} while (false);

	

	return (acceptCoincidence);

}



/*********************************************************************************

*

*		Name:		addrandUsrSampleTerminate

*

*		Summary:	Report the counter variables.

*

*			This function is assigned to the function pointer

*			addrandUsrTerminateFPtr, which is called at the

*			beginning of adrandTerminate in addrandoms.c.  This 

*			function pointer is the one to use to do any termination 

*			required by the user functions (e.g. report statistics,

*			write to/close any user files).

*

*		Arguments:

*			DetRunTimeParamsTy	*detParams

*						This is a ptr to the detector parameters. It

*						is assumed that you will treat this as 

*						a read-only structure. Modifications should 

*						only be made after a careful study of 

*						the parameters involved.

*

*		Function return: None.

*

*********************************************************************************/

void addrandUsrSampleTerminate(DetRunTimeParamsTy	*detParams)

{

	DetRunTimeParamsTy		*myDetParams;

	

	myDetParams = detParams;	/* Removes unused parameter compiler warning */

	

	/* Do your finishing here */

	/* we just report the number calls to the ModDecay functions */

	LbInPrintf("\n\n\n\tSAMPLE ADDRAND USER MODULE REPORT");

	LbInPrintf("\n\tNumber of calls to addrandUsrSampleModDecays1 is %lld",

				addrandUsrSampleVars.numModDecays1Calls);

	LbInPrintf("\n\tNumber of calls to addrandUsrSampleModDecays2 is %lld\n\n",

				addrandUsrSampleVars.numModDecays2Calls);

				

}



