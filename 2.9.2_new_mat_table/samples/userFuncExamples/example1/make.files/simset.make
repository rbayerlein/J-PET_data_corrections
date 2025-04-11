#******************************************************************************
#
#		     Source code developed by the			
#	 Imaging Research Laboratory - University of Washington		
#	     (C) Copyright 1996-2013 Department of Radiology		
#			 University of Washington			
#			    All Rights Reserved				
#									
#*****************************************************************************}

#******************************************************************************
#
#	     Module Name:       simset.make
#	     Revision Number:   1.4
#	     Date last revised: 10 January 2013
#	     Programmer:        Steven Vannoy, Steven Gillispie
#	     Date Originated:   March 20, 1996
#
#	     Module Overview:   Make file for SimSET
#
#	     References:        None
#
#******************************************************************************
#
#	     Global functions defined:		     None
#
#	     Global variables defined:		     None
#
#******************************************************************************
#
#	     Revision Section (Also update version number, if relevant)
#
#	     Programmer(s):
#
#	     Revision date:
#
#	     Revision description:
#
#*****************************************************************************}


# Setup paths 

# Change SIMSET_PATH to the directory where you installed SimSET
# (In the current version, SIMSET_PATH cannot contain any spaces.)
SIMSET_PATH = /Users/useruser/Desktop/2.9.2

# If you are building a separated user function version of SimSET, 
#   change SIMSET_PATH_USR to the top directory of your separate version.
# Note:  SIMSET_PATH_USR may contain spaces, but they must be 
#   preceded by a backslash and no quote marks should be used.
#   Example:  /this/path\ name/contains/a/space
# You will also need to change the user function file instructions at the end of this file.
# The default is to have everything in the same directory path.

# !! so for this user function example we change SIMSET_PATH_USR to the example directory
SIMSET_PATH_USR = $(SIMSET_PATH)/samples/userFuncExamples/example1

MKFILES = $(SIMSET_PATH_USR)/make.files
PHG_SRC = $(SIMSET_PATH)/src
PHG_SRC_USR = $(SIMSET_PATH_USR)/src
OBJ_DIR = $(SIMSET_PATH_USR)/obj
PROGRAMS = $(SIMSET_PATH_USR)/bin


# Set OS Flags; all Unix systems should set GEN_UNIX.
# Choose one of the following for CFLAGS:
# (Below we provide sample flags for Alphas, Macintoshes using OSX, and Linux.
# These flags will depend on the system environment--even if you have a listed 
# OS you may need to add/change flags.  If you can provide flags for another 
# operating system, we would be most grateful.)

# SUN_OS
# SGI_OS
# SGI64_OS

# ALPHA_OS
# Suggested compiler flags for ALPHA
# OS_CFLAGS = -DGEN_UNIX -DALPHA_OS -g -std1 -fptm su

# DEC_ULTRIX
# DGUX
# MPW
# _MWERKS_
# RS6000 (IBM)
# HP_OS
# WINNT

# Suggested compiler flags for Macintosh Unix (Darwin)
OS_CFLAGS = -DDARWIN -g

# LINUX
# Suggested compiler flags for LINUX
# OS_CFLAGS = -DGEN_UNIX -DLINUX


COMPILER = cc

# Select the debug or nodebug CFLAGS option:  debug has added data checking
#  and is recommended
# Choose between the last two options instead if your compiler does not yet support -iquote
CFLAGS = ${OS_CFLAGS} -c -DPHG_DEBUG -iquote"${PHG_SRC}" -DkSimSET_Path='"${SIMSET_PATH}"' 
# CFLAGS = ${OS_CFLAGS} -c -iquote"${PHG_SRC}" -DkSimSET_Path='"${SIMSET_PATH}"' 
# 
# CFLAGS = ${OS_CFLAGS} -c -DPHG_DEBUG -I"${PHG_SRC}" -DkSimSET_Path='"${SIMSET_PATH}"' 
# CFLAGS = ${OS_CFLAGS} -c -I"${PHG_SRC}" -DkSimSET_Path='"${SIMSET_PATH}"' 


PROGRAM = ${PROGRAMS}/simset
MKFILE = ${MKFILES}/simset.make


# This is a rather lazy, and not entirely accurate way
# of specifying the header dependencies, but that's the
# way it is for now.
PHG_LB_HDRS = \
	${PHG_SRC}/SystemDependent.h\
	${PHG_SRC}/LbTypes.h\
	${PHG_SRC}/LbError.h\
	${PHG_SRC}/LbDebug.h\
	${PHG_SRC}/LbEnvironment.h\
	${PHG_SRC}/LbFile.h\
	${PHG_SRC}/LbMath.h\
	${PHG_SRC}/LbMemory.h\
	${PHG_SRC}/LbParamFile.h\
	${PHG_SRC}/LbInterface.h\
	${PHG_SRC}/LbConvert.h\
	${PHG_SRC}/LbSort.h\
	${PHG_SRC}/LbTiming.h\
	${PHG_SRC}/Lb2DGeometry.h\
	${PHG_SRC}/MT19937.h

PHG_HDRS = \
	$(PHG_LB_HDRS) \
	${PHG_SRC}/Photon.h\
	${PHG_SRC}/PhgIsotopes.h\
	${PHG_SRC}/PhgParams.h\
	${PHG_SRC}/PhgMath.h\
	${PHG_SRC}/CylPos.h\
	${PHG_SRC}/ProdTbl.h\
	${PHG_SRC}/SubObj.h\
	${PHG_SRC}/EmisList.h\
	${PHG_SRC}/PhoTrk.h\
	${PHG_SRC}/PhgBin.h\
	${PHG_SRC}/ColTypes.h\
	${PHG_SRC}/Collimator.h\
	${PHG_SRC}/DetTypes.h\
	${PHG_SRC}/Detector.h\
	${PHG_SRC}/DetGeometric.h\
	${PHG_SRC}/DetPlanar.h\
	${PHG_SRC}/DetCylinder.h\
	${PHG_SRC}/DetBlock.h\
	${PHG_SRC}/ColSlat.h\
	${PHG_SRC}/addrandoms.h\
	${PHG_SRC}/addrandUsr.h\
	${PHG_SRC}/phg.h

OBJECTS = \
	${OBJ_DIR}/simset.o \
	${OBJ_DIR}/bin.phg.o \
	${OBJ_DIR}/convert.o \
	${OBJ_DIR}/collapse.o \
	${OBJ_DIR}/build_att.o \
	${OBJ_DIR}/build_coh.o \
	${OBJ_DIR}/combine.bin.o \
	${OBJ_DIR}/combine.hist.o \
	${OBJ_DIR}/reverse.bytes.o \
	${OBJ_DIR}/display.header.o \
	${OBJ_DIR}/print.header.o \
	${OBJ_DIR}/convert.coh.o \
	${OBJ_DIR}/convert.header.o \
	${OBJ_DIR}/makeindexfile.o \
	${OBJ_DIR}/phg.swap.o \
	${OBJ_DIR}/ttest.o \
	${OBJ_DIR}/strip.header.o \
	${OBJ_DIR}/migrate.o \
	${OBJ_DIR}/calc.attenuation.o \
	${OBJ_DIR}/atten.correct.o \
	${OBJ_DIR}/collapse3d.o \
	${OBJ_DIR}/extract.o \
	${OBJ_DIR}/extract.lines.o \
	${OBJ_DIR}/bcomp.o \
	${OBJ_DIR}/breakpoint.swap.o \
	${OBJ_DIR}/line3d.o \
	${OBJ_DIR}/scale.o \
	${OBJ_DIR}/reorder.o \
	${OBJ_DIR}/timesort.o \
	${OBJ_DIR}/addrandoms.o \
	${OBJ_DIR}/addrandUsr.o \
	${OBJ_DIR}/resampledecaytime.o \
	${OBJ_DIR}/Collimator.o \
	${OBJ_DIR}/ColParams.o \
	${OBJ_DIR}/ColUsr.o \
	${OBJ_DIR}/ColSlat.o \
	${OBJ_DIR}/Detector.o \
	${OBJ_DIR}/DetGeometric.o \
	${OBJ_DIR}/DetPlanar.o \
	${OBJ_DIR}/DetCylinder.o \
	${OBJ_DIR}/DetBlock.o \
	${OBJ_DIR}/DetParams.o \
	${OBJ_DIR}/DetUsr.o \
	${OBJ_DIR}/EmisList.o \
	${OBJ_DIR}/PhgBin.o \
	${OBJ_DIR}/PhgIsotopes.o \
	${OBJ_DIR}/PhgMath.o \
	${OBJ_DIR}/PhgParams.o \
	${OBJ_DIR}/PhgUsrBin.o \
	${OBJ_DIR}/PhoHFile.o \
	${OBJ_DIR}/PhoHStat.o \
	${OBJ_DIR}/ProdTbl.o \
	${OBJ_DIR}/SubObj.o \
	${OBJ_DIR}/PhgHdr.o \
	${OBJ_DIR}/PhoTrk.o \
	${OBJ_DIR}/CylPos.o \
	${OBJ_DIR}/phg.o \
	${OBJ_DIR}/UNCCollimator.o \
	${OBJ_DIR}/LbDebug.o \
	${OBJ_DIR}/LbEnvironment.o \
	${OBJ_DIR}/LbFile.o \
	${OBJ_DIR}/LbError.o \
	${OBJ_DIR}/LbInterface.o \
	${OBJ_DIR}/LbMemory.o \
	${OBJ_DIR}/LbParamFile.o \
	${OBJ_DIR}/LbHeader.o \
	${OBJ_DIR}/LbConvert.o \
	${OBJ_DIR}/LbSort.o \
	${OBJ_DIR}/LbTiming.o \
	${OBJ_DIR}/Lb2DGeometry.o \
	${OBJ_DIR}/MT19937.o



# Linking instructions

${PROGRAM}:: ${MKFILE} $(OBJECTS)
	echo linking
	${COMPILER}  -o ${PROGRAM} $(OBJECTS)   -lm 



# Compiling instructions

${OBJ_DIR}/simset.o: ${MKFILE} ${PHG_SRC}/simset.c \
				$(PHG_HDRS)
	${COMPILER}  ${CFLAGS} -o ${OBJ_DIR}/simset.o  ${CFLAGS} ${PHG_SRC}/simset.c

${OBJ_DIR}/bin.phg.o: ${MKFILE} ${PHG_SRC}/bin.phg.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/bin.phg.o ${PHG_SRC}/bin.phg.c

${OBJ_DIR}/convert.o: ${MKFILE} ${PHG_SRC}/convert.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/convert.o ${PHG_SRC}/convert.c

${OBJ_DIR}/build_att.o: ${MKFILE} ${PHG_SRC}/build_att.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/build_att.o ${PHG_SRC}/build_att.c

${OBJ_DIR}/build_coh.o: ${MKFILE} ${PHG_SRC}/build_coh.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/build_coh.o ${PHG_SRC}/build_coh.c

${OBJ_DIR}/combine.bin.o: ${MKFILE} ${PHG_SRC}/combine.bin.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/combine.bin.o ${PHG_SRC}/combine.bin.c

${OBJ_DIR}/combine.hist.o: ${MKFILE} ${PHG_SRC}/combine.hist.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/combine.hist.o ${PHG_SRC}/combine.hist.c

${OBJ_DIR}/reverse.bytes.o: ${MKFILE} ${PHG_SRC}/reverse.bytes.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/reverse.bytes.o ${PHG_SRC}/reverse.bytes.c

${OBJ_DIR}/display.header.o: ${MKFILE} ${PHG_SRC}/display.header.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/display.header.o ${PHG_SRC}/display.header.c

${OBJ_DIR}/makeindexfile.o: ${MKFILE} ${PHG_SRC}/makeindexfile.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/makeindexfile.o ${PHG_SRC}/makeindexfile.c

${OBJ_DIR}/phg.swap.o: ${MKFILE} ${PHG_SRC}/phg.swap.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/phg.swap.o ${PHG_SRC}/phg.swap.c

${OBJ_DIR}/ttest.o: ${MKFILE} ${PHG_SRC}/ttest.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/ttest.o ${PHG_SRC}/ttest.c

${OBJ_DIR}/convert.coh.o: ${MKFILE} ${PHG_SRC}/convert.coh.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/convert.coh.o ${PHG_SRC}/convert.coh.c

${OBJ_DIR}/convert.header.o: ${MKFILE} ${PHG_SRC}/convert.header.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/convert.header.o ${PHG_SRC}/convert.header.c

${OBJ_DIR}/strip.header.o: ${MKFILE} ${PHG_SRC}/strip.header.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/strip.header.o ${PHG_SRC}/strip.header.c

${OBJ_DIR}/print.header.o: ${MKFILE} ${PHG_SRC}/print.header.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/print.header.o ${PHG_SRC}/print.header.c

${OBJ_DIR}/migrate.o: ${MKFILE} ${PHG_SRC}/migrate.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/migrate.o ${PHG_SRC}/migrate.c

${OBJ_DIR}/calc.attenuation.o: ${MKFILE} ${PHG_SRC}/calc.attenuation.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/calc.attenuation.o ${PHG_SRC}/calc.attenuation.c

${OBJ_DIR}/atten.correct.o: ${MKFILE} ${PHG_SRC}/atten.correct.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/atten.correct.o ${PHG_SRC}/atten.correct.c

${OBJ_DIR}/collapse.o: ${MKFILE} ${PHG_SRC}/collapse.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/collapse.o ${PHG_SRC}/collapse.c

${OBJ_DIR}/collapse3d.o: ${MKFILE} ${PHG_SRC}/collapse3d.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/collapse3d.o ${PHG_SRC}/collapse3d.c

${OBJ_DIR}/extract.o: ${MKFILE} ${PHG_SRC}/extract.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/extract.o ${PHG_SRC}/extract.c

${OBJ_DIR}/extract.lines.o: ${MKFILE} ${PHG_SRC}/extract.lines.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/extract.lines.o ${PHG_SRC}/extract.lines.c

${OBJ_DIR}/bcomp.o: ${MKFILE} ${PHG_SRC}/bcomp.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/bcomp.o ${PHG_SRC}/bcomp.c

${OBJ_DIR}/breakpoint.swap.o: ${MKFILE} ${PHG_SRC}/breakpoint.swap.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/breakpoint.swap.o ${PHG_SRC}/breakpoint.swap.c

${OBJ_DIR}/line3d.o: ${MKFILE} ${PHG_SRC}/line3d.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/line3d.o ${PHG_SRC}/line3d.c

${OBJ_DIR}/scale.o: ${MKFILE} ${PHG_SRC}/scale.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/scale.o ${PHG_SRC}/scale.c

${OBJ_DIR}/reorder.o: ${MKFILE} ${PHG_SRC}/reorder.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/reorder.o ${PHG_SRC}/reorder.c

${OBJ_DIR}/timesort.o: ${MKFILE} ${PHG_SRC}/timesort.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/timesort.o ${PHG_SRC}/timesort.c

${OBJ_DIR}/addrandoms.o: ${MKFILE} ${PHG_SRC}/addrandoms.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/addrandoms.o ${PHG_SRC}/addrandoms.c

${OBJ_DIR}/resampledecaytime.o: ${MKFILE} ${PHG_SRC}/resampledecaytime.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/resampledecaytime.o ${PHG_SRC}/resampledecaytime.c

${OBJ_DIR}/Collimator.o: ${MKFILE} ${PHG_SRC}/Collimator.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/Collimator.o ${PHG_SRC}/Collimator.c

${OBJ_DIR}/ColParams.o: ${MKFILE} ${PHG_SRC}/ColParams.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/ColParams.o ${PHG_SRC}/ColParams.c

${OBJ_DIR}/ColSlat.o: ${MKFILE} ${PHG_SRC}/ColSlat.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/ColSlat.o ${PHG_SRC}/ColSlat.c

${OBJ_DIR}/Detector.o: ${MKFILE} ${PHG_SRC}/Detector.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/Detector.o ${PHG_SRC}/Detector.c

${OBJ_DIR}/DetGeometric.o: ${MKFILE} ${PHG_SRC}/DetGeometric.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/DetGeometric.o ${PHG_SRC}/DetGeometric.c

${OBJ_DIR}/DetPlanar.o: ${MKFILE} ${PHG_SRC}/DetPlanar.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/DetPlanar.o ${PHG_SRC}/DetPlanar.c

${OBJ_DIR}/DetCylinder.o: ${MKFILE} ${PHG_SRC}/DetCylinder.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/DetCylinder.o ${PHG_SRC}/DetCylinder.c

${OBJ_DIR}/DetBlock.o: ${MKFILE} ${PHG_SRC}/DetBlock.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/DetBlock.o ${PHG_SRC}/DetBlock.c

${OBJ_DIR}/DetParams.o: ${MKFILE} ${PHG_SRC}/DetParams.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/DetParams.o ${PHG_SRC}/DetParams.c

${OBJ_DIR}/PhgBin.o: ${MKFILE} ${PHG_SRC}/PhgBin.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/PhgBin.o ${PHG_SRC}/PhgBin.c

${OBJ_DIR}/PhgIsotopes.o: ${MKFILE} ${PHG_SRC}/PhgIsotopes.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/PhgIsotopes.o ${PHG_SRC}/PhgIsotopes.c

${OBJ_DIR}/PhgMath.o: ${MKFILE} ${PHG_SRC}/PhgMath.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/PhgMath.o ${PHG_SRC}/PhgMath.c

${OBJ_DIR}/PhgParams.o: ${MKFILE} ${PHG_SRC}/PhgParams.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/PhgParams.o ${PHG_SRC}/PhgParams.c

${OBJ_DIR}/PhoHFile.o: ${MKFILE} ${PHG_SRC}/PhoHFile.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/PhoHFile.o ${PHG_SRC}/PhoHFile.c

${OBJ_DIR}/PhoHStat.o: ${MKFILE} ${PHG_SRC}/PhoHStat.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/PhoHStat.o ${PHG_SRC}/PhoHStat.c

${OBJ_DIR}/ProdTbl.o: ${MKFILE} ${PHG_SRC}/ProdTbl.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/ProdTbl.o ${PHG_SRC}/ProdTbl.c

${OBJ_DIR}/SubObj.o: ${MKFILE} ${PHG_SRC}/SubObj.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/SubObj.o ${PHG_SRC}/SubObj.c

${OBJ_DIR}/UNCCollimator.o: ${MKFILE} ${PHG_SRC}/UNCCollimator.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/UNCCollimator.o ${PHG_SRC}/UNCCollimator.c

${OBJ_DIR}/PhgHdr.o: ${MKFILE} ${PHG_SRC}/PhgHdr.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/PhgHdr.o ${PHG_SRC}/PhgHdr.c

${OBJ_DIR}/PhoTrk.o: ${MKFILE} ${PHG_SRC}/PhoTrk.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/PhoTrk.o ${PHG_SRC}/PhoTrk.c

${OBJ_DIR}/EmisList.o: ${MKFILE} ${PHG_SRC}/EmisList.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/EmisList.o ${PHG_SRC}/EmisList.c

${OBJ_DIR}/CylPos.o: ${MKFILE} ${PHG_SRC}/CylPos.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/CylPos.o ${PHG_SRC}/CylPos.c

${OBJ_DIR}/phg.o: ${MKFILE} ${PHG_SRC}/phg.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/phg.o ${PHG_SRC}/phg.c



# Compile all library source files.

${OBJ_DIR}/LbDebug.o: ${MKFILE} ${PHG_SRC}/LbDebug.c \
				 $(PHG_LB_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/LbDebug.o ${PHG_SRC}/LbDebug.c

${OBJ_DIR}/LbEnvironment.o: ${MKFILE} ${PHG_SRC}/LbEnvironment.c \
				 $(PHG_LB_HDRS)
	${COMPILER}  ${CFLAGS} -o ${OBJ_DIR}/LbEnvironment.o ${PHG_SRC}/LbEnvironment.c

${OBJ_DIR}/LbFile.o: ${MKFILE} ${PHG_SRC}/LbFile.c \
				 $(PHG_LB_HDRS)
	${COMPILER}  ${CFLAGS} -o ${OBJ_DIR}/LbFile.o ${PHG_SRC}/LbFile.c

${OBJ_DIR}/LbError.o: ${MKFILE} ${PHG_SRC}/LbError.c \
				 $(PHG_LB_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/LbError.o  ${PHG_SRC}/LbError.c

${OBJ_DIR}/LbInterface.o: ${MKFILE} ${PHG_SRC}/LbInterface.c \
				 $(PHG_LB_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/LbInterface.o  ${PHG_SRC}/LbInterface.c

${OBJ_DIR}/LbMemory.o: ${MKFILE} ${PHG_SRC}/LbMemory.c \
				 $(PHG_LB_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/LbMemory.o ${PHG_SRC}/LbMemory.c

${OBJ_DIR}/LbParamFile.o: ${MKFILE} ${PHG_SRC}/LbParamFile.c \
				 $(PHG_LB_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/LbParamFile.o ${PHG_SRC}/LbParamFile.c

${OBJ_DIR}/LbHeader.o: ${MKFILE} ${PHG_SRC}/LbHeader.c \
				 $(PHG_LB_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/LbHeader.o ${PHG_SRC}/LbHeader.c

${OBJ_DIR}/LbConvert.o: ${MKFILE} ${PHG_SRC}/LbConvert.c \
				 $(PHG_LB_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/LbConvert.o ${PHG_SRC}/LbConvert.c

${OBJ_DIR}/LbSort.o: ${MKFILE} ${PHG_SRC}/LbSort.c \
				 $(PHG_LB_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/LbSort.o ${PHG_SRC}/LbSort.c

${OBJ_DIR}/LbTiming.o: ${MKFILE} ${PHG_SRC}/LbTiming.c \
				 $(PHG_LB_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/LbTiming.o ${PHG_SRC}/LbTiming.c

${OBJ_DIR}/Lb2DGeometry.o: ${MKFILE} ${PHG_SRC}/Lb2DGeometry.c \
				 $(PHG_LB_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/Lb2DGeometry.o ${PHG_SRC}/Lb2DGeometry.c

${OBJ_DIR}/MT19937.o: ${MKFILE} ${PHG_SRC}/MT19937.c \
				 $(PHG_LB_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/MT19937.o ${PHG_SRC}/MT19937.c



# User function compiling instructions:
# Change PHG_SRC to PHG_SRC_USR for the particular files for which you have your own versions.
# (There are two changes to make for each user file.)

# !! so for this user function example we change PHG_SRC to PHG_SRC_USR for addrandUsr.c
${OBJ_DIR}/addrandUsr.o: ${MKFILE} ${PHG_SRC_USR}/addrandUsr.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/addrandUsr.o ${PHG_SRC_USR}/addrandUsr.c

${OBJ_DIR}/ColUsr.o: ${MKFILE} ${PHG_SRC}/ColUsr.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/ColUsr.o ${PHG_SRC}/ColUsr.c

${OBJ_DIR}/DetUsr.o: ${MKFILE} ${PHG_SRC}/DetUsr.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/DetUsr.o ${PHG_SRC}/DetUsr.c

${OBJ_DIR}/PhgUsrBin.o: ${MKFILE} ${PHG_SRC}/PhgUsrBin.c \
				 $(PHG_HDRS)
	${COMPILER} ${CFLAGS} -o ${OBJ_DIR}/PhgUsrBin.o ${PHG_SRC}/PhgUsrBin.c

