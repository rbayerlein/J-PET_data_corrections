/********************************************************************************
*																				*
*                       Source code developed by the							*
*            Imaging Research Laboratory - University of Washington				*
*                (C) Copyright 1992-2013 Department of Radiology				*
*                         University of Washington								*
*                            All Rights Reserved								*
*																				*
*********************************************************************************/

/********************************************************************************
*
*     Module Name:        LbTypes.h
*     Revision Number:    2.6
*     Date last revised:  24 July 2013
*     Programmer:         Steven Vannoy
*     Date Originated:    Tuesday, July 21, 1992
*
*     Module Overview:	This module provides the core data definitions for all
*						library code.
*						MACHINE DEPENDENT						
*
*     References:         
*
**********************************************************************************
*
*     Global functions defined:
*     	LbEnGetOptions
*
*     Global variables defined:   
*
**********************************************************************************
*
*     Revision Section (Also update version number, if relevant)
*
*     Programmer(s):		
*
*     Revision date:		
*
*     Revision description:
*
**********************************************************************************
*
*     Revision Section (Also update version number, if relevant)
*
*     Programmer(s):		Steven Gillispie
*
*     Revision date:		10 January 2013
*
*     Revision description:	Changed LINUX definitions to GEN_UNIX definitions
*
**********************************************************************************
*
*     Revision Section (Also update version number, if relevant)
*
*     Programmer(s):		Steven Gillispie
*
*     Revision date:		23 June 2005
*
*     Revision description:	Modified to use sys/types.h, if available
*
*********************************************************************************
*
*			Revision Section (Also update version number, if relevant)
*
*			Programmer(s):		Robert Harrison
*
*			Revision date:		19 Nov 2004
*
*			Revision description:	Added 8 byte integer definitions
*
**********************************************************************************/

#ifndef LB_TYPE_HDR
#define LB_TYPE_HDR


/* Define this value if your system has the file /usr/include/sys/types.h
    (otherwise comment out the following line).  */
#define LB_TYPE_USE_SYS_INTS 


#include <float.h>
#ifdef LB_TYPE_USE_SYS_INTS
#include <sys/types.h>
#endif


/* Constants */
#define	LBFlag0		0x1
#define	LBFlag1		0x2
#define	LBFlag2		0x4
#define	LBFlag3		0x8
#define	LBFlag4		0x10
#define	LBFlag5		0x20
#define	LBFlag6		0x40
#define	LBFlag7		0x80
#define LBFlag8		0x100
#define LBFlag9		0x200
#define LBFlag10	0x400
#define	LBFlag11	0x800
#define LBFlag12	0x1000
#define LBFlag13	0x2000
#define LBFlag14	0x4000
#define LBFlag15	0x8000
#define LBFlag16	0x10000
#define LBFlag17	0x20000
#define LBFlag18	0x40000
#define LBFlag19	0x80000
#define LBFlag20	0x100000
#define LBFlag21	0x200000
#define LBFlag22	0x400000
#define LBFlag23	0x800000
#define LBFlag24	0x1000000
#define	LBFlag25	0x2000000
#define LBFlag26	0x4000000
#define LBFlag27	0x8000000
#define LBFlag28	0x10000000
#define LBFlag29	0x20000000
#define LBFlag30	0x40000000
#define LBFlag31	0x80000000

#define	LBEnMxArgLen	32

/* Define LIMITS */
#define LBONEBYTE_MAX		127
#define LBTWOBYTE_MAX		32767
#define LBFOURBYTE_MAX		2147483647
#define LBFOURBYTE_MIN		-LBFOURBYTE_MAX
#define LBEIGHTBYTE_MAX		9223372036854775807
#define LBEIGHTBYTE_MIN		-LBEIGHTBYTE_MAX

#define LBUSONEBYTE_MAX		255U
#define LBUSTWOBYTE_MAX		65535U
#define LBUSFOURBYTE_MAX	4294967295U
#define LBUSEIGHTBYTE_MAX	18446744073709551615

/* The floating point range comes from float.h */
#define LBFLOAT_MAX			FLT_MAX
#define LBFLOAT_MIN			FLT_MIN
#define LBDOUBLE_MAX		DBL_MAX
#define LBDOUBLE_MIN		DBL_MIN

/* Flags for types of GUI systems */
#ifdef __MWERKS__
	#define GUIOS
	#define MACGUIOS
#endif
#ifdef MPW
	#define GUIOS
	#define MACGUIOS
#endif
#ifdef OSX
	#define GUIOS
	#define MACGUIOS
#endif
#ifdef COCOA
	#define GUIOS
	#define MACGUIOS
#endif

#ifdef LB_TYPE_USE_SYS_INTS
	/* Use the system-provided types in /usr/include/sys/types.h */
	/*   (see the #define for LB_TYPE_USE_SYS_INTS above).  */
	
	typedef int8_t			LbOneByte;
	typedef	int16_t			LbTwoByte;
	typedef	int32_t			LbFourByte;
	typedef	int64_t			LbEightByte;
	
	typedef	u_int8_t		LbUsOneByte;
	typedef	u_int16_t		LbUsTwoByte;
	typedef	u_int32_t		LbUsFourByte;
	typedef	u_int64_t		LbUsEightByte;
	
	typedef void			*LbUserData;
	
	/* Comment out the following line declaring type Boolean if your system complains --
		it must already have the type Boolean defined */
	typedef LbTwoByte		Boolean;
	#define true	1
	#define	false	0
	
	/* Set our type flag */
	#define OS_SUPPORTED
#else
	#ifdef GEN_UNIX
		typedef signed char		LbOneByte;
		typedef	short			LbTwoByte;
		typedef	long			LbFourByte;
		typedef	long long		LbEightByte;
		
		typedef	unsigned char	LbUsOneByte;
		typedef	unsigned short	LbUsTwoByte;
		typedef	unsigned long	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void			*LbUserData;

		/* Set our type flag */
		#define OS_SUPPORTED

	#endif /* GEN_UNIX */

	#ifdef DARWIN
		typedef signed char	LbOneByte;
		typedef	short	LbTwoByte;
		typedef	long	LbFourByte;
		typedef	long long		LbEightByte;
		
		typedef	unsigned char	LbUsOneByte;
		typedef	unsigned short	LbUsTwoByte;
		typedef	unsigned long	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void	*LbUserData;
		/* Set our type flag */
		#define OS_SUPPORTED
	#endif

	#ifdef OSX
		typedef signed char			LbOneByte;
		typedef	short				LbTwoByte;
		typedef	long				LbFourByte;
		typedef	long long			LbEightByte;
		
		typedef	unsigned char		LbUsOneByte;
		typedef	unsigned short		LbUsTwoByte;
		typedef	unsigned long		LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void				*LbUserData;
		
		/* Set our type flag */
		#define OS_SUPPORTED
	#endif

	#ifdef COCOA
		typedef signed char			LbOneByte;
		typedef	short				LbTwoByte;
		typedef	long				LbFourByte;
		typedef	long long			LbEightByte;
		
		typedef	unsigned char		LbUsOneByte;
		typedef	unsigned short		LbUsTwoByte;
		typedef	unsigned long		LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void				*LbUserData;
		
		/* Set our type flag */
		#define OS_SUPPORTED
	#endif

	/* Cross Platform Types NOTE: as of 5/19/98 we use Metrowerks version of Boolean */
	#ifndef __MWERKS__
		#ifndef MPW
			#ifndef OSX
				typedef int Boolean;
				#define true	1
				#define	false	0
			#endif
		#endif
	#endif
	
	#ifdef FREE_BSD
		typedef signed char	LbOneByte;
		typedef short		LbTwoByte;
		typedef long		LbFourByte;
		typedef	long long		LbEightByte;

		typedef unsigned char	LbUsOneByte;
		typedef unsigned short	LbUsTwoByte;
		typedef unsigned long	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;

		typedef void		*LbUserData;

		/* Set our type flag */
		#define OS_SUPPORTED
	#endif /* FREE_BSD */

	/* Platform Dependant Types */
	#ifdef SGI_OS
		typedef signed char	LbOneByte;
		typedef	short	LbTwoByte;
		typedef	long	LbFourByte;
		typedef	long long		LbEightByte;
		
		typedef	unsigned char	LbUsOneByte;
		typedef	unsigned short	LbUsTwoByte;
		typedef	unsigned long	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void	*LbUserData;

		/* Set our type flag */
		#define OS_SUPPORTED
	#endif /* SGI_OS */

	#ifdef SGI64_OS
		typedef signed char	LbOneByte;
		typedef	short	LbTwoByte;
		typedef	int		LbFourByte;
		typedef	long long		LbEightByte;
		
		typedef	unsigned char	LbUsOneByte;
		typedef	unsigned short	LbUsTwoByte;
		typedef	unsigned int	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void	*LbUserData;

		/* Set our type flag */
		#define OS_SUPPORTED
	#endif /* SGI64_OS */

	#ifdef ALPHA_OS
		typedef signed char	LbOneByte;
		typedef	short	LbTwoByte;
		typedef	int		LbFourByte;
		typedef	long long		LbEightByte;
		
		typedef	unsigned char	LbUsOneByte;
		typedef	unsigned short	LbUsTwoByte;
		typedef	unsigned int	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void 			*LbUserData;
		
		/* Set our type flag */
		#define OS_SUPPORTED
	#endif /* ALPHA_OS */

	#ifdef DEC_ULTRIX
		typedef signed char	LbOneByte;
		typedef	short	LbTwoByte;
		typedef	long	LbFourByte;
		typedef	long long		LbEightByte;
		
		typedef	unsigned char	LbUsOneByte;
		typedef	unsigned short	LbUsTwoByte;
		typedef	unsigned long	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void	*LbUserData;
		
		/* Set our type flag */
		#define OS_SUPPORTED
	#endif /* DEC_ULTRIX */


	#ifdef DGUX
		typedef signed char	LbOneByte;
		typedef	short	LbTwoByte;
		typedef	long	LbFourByte;
		typedef	long long		LbEightByte;
		
		typedef	unsigned char	LbUsOneByte;
		typedef	unsigned short	LbUsTwoByte;
		typedef	unsigned long	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void	*LbUserData;

		/* Set our type flag */
		#define OS_SUPPORTED
		
	#endif /*	DGUX */

	#ifdef SUN_OS
		typedef signed char	LbOneByte;
		typedef	short	LbTwoByte;
		typedef	long	LbFourByte;
		typedef	long long		LbEightByte;
		
		typedef	unsigned char	LbUsOneByte;
		typedef	unsigned short	LbUsTwoByte;
		typedef	unsigned long	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void	*LbUserData;

		/* Set our type flag */
		#define OS_SUPPORTED
		
	#endif	/* SUN_OS */

	#ifdef MPW
		typedef signed char	LbOneByte;
		typedef	short	LbTwoByte;
		typedef	long	LbFourByte;
		typedef	long long		LbEightByte;
		
		typedef	unsigned char	LbUsOneByte;
		typedef	unsigned short	LbUsTwoByte;
		typedef	unsigned long	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void	*LbUserData;

		/* Set our type flag */
		#define OS_SUPPORTED
		
	#endif	/* MPW */

	#ifdef __MWERKS__
		typedef signed char	LbOneByte;
		typedef	short	LbTwoByte;
		typedef	long	LbFourByte;
		typedef	long long		LbEightByte;
		
		typedef	unsigned char	LbUsOneByte;
		typedef	unsigned short	LbUsTwoByte;
		typedef	unsigned long	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void	*LbUserData;

		/* This should be defined in the Mac headers, but is not */
		#define M_PI	PI
		/* Set our type flag */
		#define OS_SUPPORTED
		
	#endif	/* __MWERKS__ */

	#ifdef RS6000
		typedef signed char	LbOneByte;
		typedef	short	LbTwoByte;
		typedef	long	LbFourByte;
		typedef	long long		LbEightByte;
		
		typedef	unsigned char	LbUsOneByte;
		typedef	unsigned short	LbUsTwoByte;
		typedef	unsigned long	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void	*LbUserData;

		/* Set our type flag */
		#define OS_SUPPORTED
		
	#endif	/* RS6000 */

	#ifdef HP_OS
		typedef signed char	LbOneByte;
		typedef	short	LbTwoByte;
		typedef	long	LbFourByte;
		typedef	long long		LbEightByte;
		
		typedef	unsigned char	LbUsOneByte;
		typedef	unsigned short	LbUsTwoByte;
		typedef	unsigned long	LbUsFourByte;
		typedef	unsigned long long	LbUsEightByte;
		
		typedef void	*LbUserData;

		/* Set our type flag */
		#define OS_SUPPORTED
		
	#endif	/* HP_OS */

	#ifdef WINNT
		typedef __int8				LbOneByte;
		typedef	__int16				LbTwoByte;
		typedef	__int32				LbFourByte;
		typedef	__int64				LbEightByte;
		
		typedef	unsigned __int8		LbUsOneByte;
		typedef	unsigned __int16	LbUsTwoByte;
		typedef	unsigned __int32	LbUsFourByte;
		typedef	unsigned __int64	LbUsEightByte;
		
		typedef void				*LbUserData;

		typedef LbTwoByte			Boolean;
		#define true	1
		#define	false	0
		
		/* Set our type flag */
		#define OS_SUPPORTED

	#endif	/* WINNT */

#endif


/* Verify they have defined a supported OS */
#ifndef OS_SUPPORTED
You have not declared an OS-dependent clause in LbTypes.h
#endif

#endif /* LB_TYPE_HDR */

