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
*			Module Name:		addrandUsr.h
*			Revision Number:	2.2
*			Date last revised:	25 June 2012
*			Programmer:			R Harrison
*			Date Originated:	17 February 2006
*
*			Module Overview:	Definitions for addrandUsr.c.
*
*			References:			
*
**********************************************************************************
*
*			Global functions defined:	
*
*			Global variables defined:
*					addrandUsrParamsFPtr
*					addrandUsrModDecays1FPtr
*					addrandUsrModDecays2FPtr
*					addrandUsrTerminateFPtr
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
*			Revision date:		25 June 2012
*
*			Revision description:	Added parameters to the two detection functions
*
**********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Steven Gillispie
*
*			Revision date:		10 January 2012
*
*			Revision description:	Changed form of user functions to pointers
*
*********************************************************************************/

#ifndef ADDRAND_USR_HDR
#define ADDRAND_USR_HDR

#ifdef ADDRAND_USR
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* DEFINITIONS */
typedef 	void addrandUsrParamsFType(DetRunTimeParamsTy	*detParams);
typedef 	void addrandUsrDets1FType(	timeWindowDetectionsTy	*timeWindowDetectionsPtr,
										PhoHFileHkTy			*adrandRandomsFileHk, 
										LbUsFourByte			*numDecaysWritten, 
										LbUsFourByte			*numDecaysUnchanged, 
										LbUsFourByte			*numDecaysRandom, 
										LbUsFourByte			*numLostCorrectWindow, 
										LbUsFourByte			*numDecaysLostTriples);
typedef 	Boolean addrandUsrDets2FType(	timeWindowDetectionsTy	*timeWindowDetectionsPtr,
											PhoHFileHkTy			*adrandRandomsFileHk, 
											LbUsFourByte			*numDecaysWritten, 
											LbUsFourByte			*numDecaysUnchanged, 
											LbUsFourByte			*numDecaysRandom, 
											LbUsFourByte			*numLostCorrectWindow, 
											LbUsFourByte			*numDecaysLostTriples);

/* GLOBALS */
extern addrandUsrParamsFType	*addrandUsrInitializeFPtr;
extern addrandUsrDets1FType		*addrandUsrModDecays1FPtr;
extern addrandUsrDets2FType		*addrandUsrModDecays2FPtr;
extern addrandUsrParamsFType	*addrandUsrTerminateFPtr;

/* PROTOTYPES */

#undef LOCALE
#endif /* ADDRAND_USR_HDR */
