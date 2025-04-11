#!/bin/sh

# This file can be used as a shell script for building SimSET.
# Note that it may not work conveniently if the user has
#   modified their environment.
#
# THIS SHELL SCRIPT SHOULD BE EXECUTED FROM WITHIN THE SimSET DIRECTORY

# Check for spaces in the SimSET directory path
dirspc=`pwd | grep ' '` 
# echo $dirspc
if test -n "$dirspc"; then
    echo "The directory path for SimSET cannot contain any spaces, but this one does:"
    echo "  "`pwd`
    exit;
fi

# Remove any current object files
rm -f obj/*.o

# Change working directory to make.files
cd make.files

# Make SimSET in full
make -f simset.make

# Go into bin and create links
cd ../bin

# First create a link to the phg only
ln -s -f simset phg

# Now create a link for each utility
ln -s -f simset bin
ln -s -f simset phgbin
ln -s -f simset convert
ln -s -f simset buildatt
ln -s -f simset combinebin
ln -s -f simset combinehist
ln -s -f simset reversebytes
ln -s -f simset displayheader
ln -s -f simset makeindexfile
ln -s -f simset ttest
ln -s -f simset stripheader
ln -s -f simset migrate
ln -s -f simset phgswap
ln -s -f simset printheader
ln -s -f simset convertheader
ln -s -f simset calcattenuation
ln -s -f simset attncorrect
ln -s -f simset collapse
ln -s -f simset extract
ln -s -f simset collapse3d
ln -s -f simset bcomp
ln -s -f simset extractlines
ln -s -f simset buildcoh
ln -s -f simset line3d
ln -s -f simset scale
ln -s -f simset convertcoh
ln -s -f simset timesort
ln -s -f simset addrandoms
ln -s -f simset resampledecaytimes


# Finally, change back to the original directory
cd ..


