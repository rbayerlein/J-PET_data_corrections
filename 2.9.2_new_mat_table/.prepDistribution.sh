##################################################
#
#
#       This script clears the bin and obj
#		directories.  It also clears results
#		out of the samples directory.  This
#       reduces the size of the distribution
#       package and insures that users will
#       not accidentally run the results script
#       and compare the base results with
#       UW results instead of their own.
#
#
##################################################

cd bin
rm *
cd ..

cd obj
rm *
cd ..

cd samples/fastTest
.clearResults.sh
cd ../..

cd samples/userFuncExamples
.clearBinaries.sh
.clearResults.sh
cd ../..

