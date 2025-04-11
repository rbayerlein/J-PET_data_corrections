/*********************************************************************************
*                                                                                *
*                       Source code developed by the                             *
*           Imaging Research Laboratory - University of Washington               *
*               (C) Copyright 1992-2000 Department of Radiology                  *
*                           University of Washington                             *
*                              All Rights Reserved                               *
*                                                                                *
*********************************************************************************/
 
/*********************************************************************************
*
*			Module Name:		ProdTbl.h
*			Revision Number:	1.0
*			Date last revised:	Wednesday, September 9, 1992
*			Programmer:			Steven Vannoy
*			Date Originated:	Wednesday, September 9, 1992
*
*			Module Overview:	Definitions for ProdTbl.c.
*
*			References:			'Productivity Table Processes' PHG design.
*
**********************************************************************************
*
*			Global functions defined:	
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
*********************************************************************************/

#ifndef PROD_TABLE
#define PROD_TABLE

#ifdef PRODUCTIVITY_TABLE
	#define	LOCALE
#else
	#define LOCALE	extern
#endif

/* CONSTANTS */
#define PRODTBL_TOT_STRAT_CELLS		60					/* Count of stratification cells */
#define PRODTBL_ACC_STRAT_CELLS		48					/* Count of stratification cells within acceptance angle */
#define PRODTBL_NOTACC_STRAT_CELLS	12					/* Count of stratification cells outside acceptance angle */


#define PRODTBLFg_Scatter			LBFlag0				/* Reference scatter productivity */
#define PRODTBLFg_Primary			LBFlag1				/* Reference primary productivity */

/* TYPES */
	/* Initial Productivity Array */
typedef struct {
	double	startOfBoundary;						/* Start of cells boundary */
	double	endOfBoundary;							/* End of cells boundary */
	double	cellProductivity;						/* Productivity of this cell */
} ProdTblProdCellTy;
typedef ProdTblProdCellTy		ProdTblProdArrayTy[PRODTBL_TOT_STRAT_CELLS];		/* Productivity cells for stratified angles */

typedef struct {
	double zMin;									/* Minimum z for this slice */
	double zMax;									/* Maximum z for this slice */
	ProdTblProdArrayTy	productivity;				/* Productivities for this slice */
}ProdTblProdTblElemTy, *ProdTblProdTblTy;


	/* Productivity table information */
typedef struct {
	char			*inputFileName;				/* Productivity input file */
	char			*outputFileName;			/* Productivity output file */
	LbUsFourByte	numSlices;					/* Number of slices in object/productivity table */
	double			acceptanceAngle;			/* The acceptance angle */
} ProdTblProdTblInfoTy;


/* GLOBALS */
LOCALE	LbUsFourByte				ProdTblNumAngleCells;			/* Number of stratification cells */
LOCALE	LbUsFourByte				ProdTblNumAccAngleCells;		/* Number of acceptable stratification cells */
LOCALE	LbUsFourByte				ProdTblNumUnAccAngleCells;		/* Number of unacceptable stratification cells */
LOCALE	LbUsFourByte				ProdTblNumSlices;				/* Number of slices in table */
LOCALE	ProdTblProdTblTy			ProdTblPrimProdTbl;				/* Our productivity table for primary photons */
LOCALE	ProdTblProdTblTy			ProdTblScatProdTbl;				/* Our productivity table for scattered photons */
LOCALE	ProdTblProdTblTy			ProdTblMaxProdTbl;				/* Our productivity table for scattered photons */


/* MACROS */
/*********************************************************************************
*
*			Name:			PRODTBLGetProdTblAngleEnd
*
*			Summary:		Return the endOfBoundary for a given slice/angle.
*
*			Arguments:
*				LbUsFourByte	sliceIndex	- Current slice.
*				LbUsFourByte	angleIndex	- Current angle.
*
*			Function return: Boundary end of given angle cell.
*
*********************************************************************************/
#define PRODTBLGetProdTblAngleEnd(sliceIndex, angleIndex) \
	(ProdTblPrimProdTbl[(sliceIndex)].productivity[(angleIndex)].endOfBoundary)
	
/*********************************************************************************
*
*			Name:			PRODTBLGetProdTblAngleSize
*
*			Summary:		Return the size of the given angle.
*
*			Arguments:
*				LbUsFourByte	sliceIndex	- Current slice.
*				LbUsFourByte	angleIndex	- Current angle.
*
*			Function return: Size of given angle cell.
*
*********************************************************************************/
#define PRODTBLGetProdTblAngleSize(sliceIndex, angleIndex) \
	(fabs(ProdTblPrimProdTbl[(sliceIndex)].productivity[(angleIndex)].startOfBoundary - \
	ProdTblPrimProdTbl[(sliceIndex)].productivity[(angleIndex)].endOfBoundary))

/*********************************************************************************
*
*			Name:			PRODTBLGetProdTblAngleStart
*
*			Summary:		Return the startOfBoundary for a given slice/angle.
*
*			Arguments:
*				LbUsFourByte	sliceIndex	- Current slice.
*				LbUsFourByte	angleIndex	- Current angle.
*
*			Function return: Boundary start of given angle cell.
*
*********************************************************************************/
#define PRODTBLGetProdTblAngleStart(sliceIndex, angleIndex) \
	(ProdTblPrimProdTbl[(sliceIndex)].productivity[(angleIndex)].startOfBoundary)

/*********************************************************************************
*
*			Name:			PRODTBLGetPrimCellProductivity
*
*			Summary:		Return the productivity for a given cell.
*
*			Arguments:
*				LbUsFourByte	sliceIndex	- Current slice.
*				LbUsFourByte	angleIndex	- Current angle.
*
*			Function return: Productivity of given angle cell.
*
*********************************************************************************/
#define PRODTBLGetPrimCellProductivity(sliceIndex, angleIndex) \
	ProdTblPrimProdTbl[(sliceIndex)].productivity[(angleIndex)].cellProductivity

/*********************************************************************************
*
*			Name:			PRODTBLGetScatCellProductivity
*
*			Summary:		Return the productivity for a given cell.
*
*			Arguments:
*				LbUsFourByte	sliceIndex	- Current slice.
*				LbUsFourByte	angleIndex	- Current angle.
*
*			Function return: Productivity of given angle cell.
*
*********************************************************************************/
#define PRODTBLGetScatCellProductivity(sliceIndex, angleIndex) \
	ProdTblScatProdTbl[(sliceIndex)].productivity[(angleIndex)].cellProductivity

/*********************************************************************************
*
*			Name:			PRODTBLGetMaxCellProductivity
*
*			Summary:		Return the productivity for a given cell.
*
*			Arguments:
*				LbUsFourByte	sliceIndex	- Current slice.
*				LbUsFourByte	angleIndex	- Current angle.
*
*			Function return: Productivity of given angle cell.
*
*********************************************************************************/
#define PRODTBLGetMaxCellProductivity(sliceIndex, angleIndex) \
	ProdTblMaxProdTbl[(sliceIndex)].productivity[(angleIndex)].cellProductivity

/*********************************************************************************
*
*			Name:			PRODTBLGetNumAngleCells
*			Summary:		Return the number of stratification cells.
*
*			Arguments:
*
*			Function return: Number of stratification cells
*
*********************************************************************************/
#define PRODTBLGetNumAngleCells()  (ProdTblNumAngleCells)


/* PROTOTYPES */
void 			ProdTblAddDetectedProductivity(PHG_TrackingPhoton *trackingPhotonPtr,
					LbUsFourByte whichOne);
void 			ProdTblAddStartingProductivity(PHG_TrackingPhoton *trackingPhotonPtr,
					LbUsFourByte whichOne);	
Boolean 		ProdTblCloseTable(ProdTblProdTblInfoTy *prodTableInfoPtr);
Boolean 		ProdTblCreateTable(ProdTblProdTblInfoTy *prodTablInfoPtr);	
LbFourByte		ProdTblGetAngleIndex(LbFourByte sliceIndex, double z_cosine);
void			ProdTblGetZmaxZmin(LbUsFourByte sliceIndex, double *zMax, double *zMin);
Boolean			ProdTblInitialize(void);
void			ProdTblTerminate(ProdTblProdTblInfoTy *prodTableInfoPtr);

#ifdef PHG_DEBUG
LOCALE	void	ProdTblDumpObjects(void);
#endif
#undef LOCALE
#endif /* PROD_TABLE */
