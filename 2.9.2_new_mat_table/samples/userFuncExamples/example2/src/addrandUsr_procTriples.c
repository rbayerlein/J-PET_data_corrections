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
*			Module Name:		addrandUsr_procTriples.c
*			Revision Number:	0.1
*			Date last revised:	24 May 2013
*			Programmer:			Robert Harrison
*			Date Originated:	19 June 2012
*
*			Module Overview:	This is an include file for addrandUsr.c that,
*								when activated, changes SimSET's
*								handling of cases where three or more photons
*								are detected within the coincidence time window
*								(triples).  SimSET's default action is to discard
*								all photons in such windows;  the main function 
*								in this module, addrandUsr_procTriplesProcess,
*								instead accepts all combinations of two photons
*								that fall within the user-specified coincidence
*								resolution time of each other.
*
*								TO ACTIVATE THIS USER MODULE: in the local
*								prototypes section of addrandUsr.c, reassign  
*								the pointers addrandUsrInitializeFPtr,
*								addrandUsrModDecays1FPtr, and 
*								addrandUsrTerminateFPtr from NULL to point to
*								the appropriate functions in 
*								addrandUsr_procTriples.c.  SIMSET MUST THEN
*								BE RECOMPILED TO USE THE CODE.
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
**********************************************************************************/


/* global variables for procTriples */
/*
##
##
	Local globals declared within the procTriplesVars struct below.  This will prevent any
		problems with naming and scope.
##
##
*/
struct {
	/* Your fields of choice go here */
	double			adrandTimingWindow;		/* the user-supplied timing window converted to seconds */
	LbEightByte		num_single_decay_from_triples;	/* number of single decay events recovered from triples processing */
	LbEightByte		num_randoms_from_triples;	/* number of randoms added from triples processing */
} procTriplesVars;


/* Local Prototypes for the functions declared here */
void procTriplesInitialize(DetRunTimeParamsTy	*detParams);
void procTriplesProcess(	timeWindowDetectionsTy	*timeWindowDetectionsPtr,
							PhoHFileHkTy			*adrandRandomsFileHk, 
							LbUsFourByte			*numDecaysWritten, 
							LbUsFourByte			*numDecaysUnchanged, 
							LbUsFourByte			*numDecaysRandom, 
							LbUsFourByte			*numLostCorrectWindow, 
							LbUsFourByte			*numDecaysLostTriples);
void procTriplesTerminate(DetRunTimeParamsTy	*detParams);


/*********************************************************************************
*
*		Name:			procTriplesInitialize
*
*		Summary:	We initialize two counter variables and copy the coincidence
*				timing window to a local variable.
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
void procTriplesInitialize(DetRunTimeParamsTy	*detParams)
{
	/* Do your initialization here */
	
	/* initialize counters for the number of events added with this triples processing */
	procTriplesVars.adrandTimingWindow = detParams->CoincidenceTimingWindowNS * 1.0E-9;
	procTriplesVars.num_single_decay_from_triples = 0;
	procTriplesVars.num_randoms_from_triples = 0;
	
}

/*********************************************************************************
*
*		Name:		procTriplesProcess
*
*		Summary:	This routine is called before the photons/decays are
*				time-windowed.  All the decays within a time-window of each
*				other (or actually longer: the window is extended each
*				time another photon arrives) are passed in the 
*				timeWindowDetectionsTy structure.  The photons and decays
*				may be altered/deleted as desired. 
*				
*				In this example we go through any window with more than
*				three photons in it and create an event for every possible
*				pair of photons.  For two photons from the same decay the
*				event will be labelled as normal PET events, for two photons 
*				from different decays it will be labeled random.  These events
*				are written to the addrandoms history file.  The time window
*				will be returned with no decays to avoid further processing in 
*				adrandProcessTimeWindow.
*
*				We ignore windows with exactly one or exactly two photons:
*				these windows are correctly handled by the calling function.
*				 
*				Note that any operations specified by a function pointed to
*				by addrandUsrModDecays2FPtr will not be carried out on decays
*				created in this function.  If the user creates such a function
*				they may wish to call it here in addition to the call in
*				adrandProcessTimeWindow.
*
*		Arguments:
*			timeWindowDetectionsTy	*timeWindowDetections - 
*					All the decays in the current time window
*			PhoHFileHkTy			*adrandRandomsFileHk -
*					pointer to addrandoms output history file
*			LbUsFourByte			*numDecaysWritten -
*					addrandoms counter for the number of coincidences written
*					to the history file
*			LbUsFourByte			*numDecaysUnchanged -
*					addrandoms counter for the number of single decay 
*					coincidences written to the history file
*			LbUsFourByte			*numDecaysRandom -
*					addrandoms counter for the number of random 
*					coincidences written to the history file
*			LbUsFourByte			*numLostCorrectWindow -
*					addrandoms counter for the number of coincidences 
*					from a time window not written to the history file
*					because the time difference between photons was >
*					adrandTimingWindow
*			LbUsFourByte			*numDecaysLostTriples -
*					addrandoms counter giving the number of decays 
*					that were not written out because there were >2 photons
*					in the time window (as happens in the default processing,
*					but not when using this user module)
*
*		Function return: None
*
*********************************************************************************/
void procTriplesProcess(	timeWindowDetectionsTy	*timeWindowDetectionsPtr,
							PhoHFileHkTy			*adrandRandomsFileHk, 
							LbUsFourByte			*numDecaysWritten, 
							LbUsFourByte			*numDecaysUnchanged, 
							LbUsFourByte			*numDecaysRandom, 
							LbUsFourByte			*numLostCorrectWindow, 
							LbUsFourByte			*numDecaysLostTriples)
{
	LbUsFourByte		decayNum;				/* loop counter */
	LbUsFourByte		totalPhotons = 0;		/* total photons in this time window */
	LbUsFourByte		photonNum, photonNum2;	/* loop counters */
	PHG_Decay			triple_decay;			/* a decay created from a triple */
	
	typedef struct {
		LbUsFourByte		decayNum;			/* decay number for this photon */
		double				decayTime;			/* decayTime for this photon's decay */
		PHG_Position		decayLocation;		/* Origination of decay */
		PHG_DetectedPhoton	*photonPtr;			/* pointer to photon */
	}triple_photonTy;
	
	triple_photonTy	*triple_photon_list;		/* list of the photons to match up as triples */
	
	/* Make your alterations here */
	do {
		
		/* count up the total photons in this time window */
		for ( decayNum = 0; decayNum < timeWindowDetectionsPtr->numDecays; decayNum++ ) {
			
			/* there should be at most one blue and one pink photon per decay - 
			  randoms processing requires that importance sampling be turned off */
			if (	(timeWindowDetectionsPtr->decays[decayNum].numBluePhotons > 1) ||
					(timeWindowDetectionsPtr->decays[decayNum].numPinkPhotons > 1) ) {
				PhgAbort("procTriplesProcess: too many photons in decay - IS should be false.", false);
			}
			
			totalPhotons += timeWindowDetectionsPtr->decays[decayNum].numBluePhotons;
			totalPhotons += timeWindowDetectionsPtr->decays[decayNum].numPinkPhotons;
			
		}
		
		if ( totalPhotons < 3 ) { break; }; /* these windows handled by mainline code */
		
		/* allocate a list of all the photons in the time window */
		if ((triple_photon_list = (triple_photonTy *)LbMmAlloc(
					sizeof(triple_photonTy) * totalPhotons)) == 0) {
			PhgAbort("procTriplesProcess: unable to allocate triple_photon_list.", false);
		}
		
		/* fill a list of all the photons in the time window */
		photonNum = 0;
		for ( decayNum = 0; decayNum < timeWindowDetectionsPtr->numDecays; decayNum++ ) {
			
			if (timeWindowDetectionsPtr->decays[decayNum].numBluePhotons == 1) {
				triple_photon_list[photonNum].decayNum = decayNum;
				triple_photon_list[photonNum].decayTime = timeWindowDetectionsPtr->decays[decayNum].decay.decayTime;
				triple_photon_list[photonNum].decayLocation = timeWindowDetectionsPtr->decays[decayNum].decay.location;
				triple_photon_list[photonNum].photonPtr = &(timeWindowDetectionsPtr->decays[decayNum].bluePhotons[0]);
				photonNum++;
			};
			
			if (timeWindowDetectionsPtr->decays[decayNum].numPinkPhotons == 1) {
				triple_photon_list[photonNum].decayNum = decayNum;
				triple_photon_list[photonNum].decayTime = timeWindowDetectionsPtr->decays[decayNum].decay.decayTime;
				triple_photon_list[photonNum].decayLocation = timeWindowDetectionsPtr->decays[decayNum].decay.location;
				triple_photon_list[photonNum].photonPtr = &(timeWindowDetectionsPtr->decays[decayNum].pinkPhotons[0]);
				photonNum++;
			};
			
		}
		
		/* the number of photons put in the triple_photon_list should be the same as totalPhotons */
		if (photonNum != totalPhotons) {
			PhgAbort("procTriplesProcess: wrong number of photons in triple_photon_list.", false);
		}

		/* step through the photons in the time window, matching them up whenever they
		 fall within the coincidence window (note that it is possible for photons within
		 a single timeWindowDetections to be detected with a time difference greater than
		 the coindidence resolving time, see addrandoms.c for details). Events are
		 labeled random unless they come from the same decay */
		
		/* setup the decay fields that won't change as we step thru the photons... */
		triple_decay.startWeight = 1.0;
		triple_decay.decayTime = triple_photon_list[0].decayTime;
		
		/* adjust  photons' time_since_creation to account for the decayTime 
		  difference between the decays for the 0th photon and  other photons */
		for ( photonNum = 1; photonNum < totalPhotons; photonNum++ ) {
			
			triple_photon_list[photonNum].photonPtr->time_since_creation +=
				triple_photon_list[photonNum].decayTime - triple_photon_list[0].decayTime;
			
		}
		
		/* now through all photon pairs */
		for ( photonNum = 0; photonNum < (totalPhotons - 1); photonNum++ ) {
			
			for ( photonNum2 = (photonNum + 1); photonNum2 < totalPhotons; photonNum2++ ) {
				
				/* check to see if the two photons  pass within time window */
				if ( fabs( triple_photon_list[photonNum2].photonPtr->time_since_creation -
							triple_photon_list[photonNum].photonPtr->time_since_creation )
							<= procTriplesVars.adrandTimingWindow ) {
					
					if ( triple_photon_list[photonNum].decayNum == triple_photon_list[photonNum2].decayNum ) {
						
						/* this is a single decay event that would have been skipped with the normal triples processing */
						triple_decay.decayType = PhgEn_Positron;
						triple_decay.location = triple_photon_list[photonNum].decayLocation;
						
						procTriplesVars.num_single_decay_from_triples++;
						*numDecaysUnchanged = *numDecaysUnchanged + 1;
						
					} else {
						
						/* this is a random event that would have been skipped with the normal triples processing */
						triple_decay.decayType = PhgEn_PETRandom;
						triple_decay.location = triple_photon_list[photonNum].decayLocation;	/* meaningless! */
						
						/* assign the first photon as the blue, second as the pink */
						LbFgSet(triple_photon_list[photonNum].photonPtr->flags, PHGFg_PhotonBlue);
						LbFgClear(triple_photon_list[photonNum2].photonPtr->flags, PHGFg_PhotonBlue);
						
						procTriplesVars.num_randoms_from_triples++;
						*numDecaysRandom = *numDecaysRandom + 1;
						
					}
						
					/* write this decay to output list mode file */
					PhoHFileRewriteDetections( adrandRandomsFileHk, &triple_decay,
						triple_photon_list[photonNum].photonPtr, 1,
						triple_photon_list[photonNum2].photonPtr, 1);
					
					/* increment the number of decays written out and the number of randoms created */
					*numDecaysWritten = *numDecaysWritten + 1;
											
				} else {
					
					/* keep track of how many randoms are lost to this second time windowing */
					*numLostCorrectWindow = *numLostCorrectWindow + 1;
					
				}
			}
			
		}
		
		/* set the number of decays in this time window to 0 so that the decays won't be reprocessed -
		  note we only get here when totalPhotons is >= 3 and such windows have just been processed */
		timeWindowDetectionsPtr->numDecays = 0;
		
		/* free up memory */
		if (triple_photon_list != 0) {
			LbMmFree((void **)&triple_photon_list);
		}
		
	} while (false);
	
	return;
}

/*********************************************************************************
*
*		Name:		procTriplesTerminate
*
*		Summary:	Reports the final values of the local counter variables.
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
void procTriplesTerminate(DetRunTimeParamsTy	*detParams)
{
	DetRunTimeParamsTy		*myDetParams;
	
	myDetParams = detParams;	/* Removes unused parameter compiler warning */
	
	/* Do your finishing here */
	/* we just report the number of single decay and random events were added from triples processing */
	LbInPrintf("\n\n\n\tTRIPLES PROCESSING USER MODULE REPORT");
	LbInPrintf("\n\tNumber of single decay events created from triples is %lld",
				procTriplesVars.num_single_decay_from_triples);
	LbInPrintf("\n\tNumber of random events created from triples is %lld\n\n",
				procTriplesVars.num_randoms_from_triples);

}


