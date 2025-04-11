#!/bin/bash

# make file for all executables required by JPET image reconstruction framework
#
# REQUIRED: 
# - MATLAB installed version 2020b or higher
# - Downloaded the material tables for SimSET and saved them in <path/to/jpet>/2.9.2_new_mat_table/phg.data/coh.tables/

# ATTENTION:
# Castor is NOT installed when running this script.
# please install it youself and choose the configuration you need

################################################################
# rbayerlein@ucdavis.edu
# April 2025
################################################################

# lm2sino5d_castor
cd lm2sino5d_castor
make clean
make
cd ..

# delayed2randoms4d
cd delayed2randoms4d
make clean
make
cd ..

# add_randoms_factors
cd add_randoms_factors
make clean
make
cd ..

# castor
# Please install Castor using the README inside the castor folder

# SimSET
cd 2.9.2_new_mat_table
chmod u+x make_all.sh	# make executable for user
./make_all.sh
cd ..

# hist2lm
cd data_processing/hist2lm/src
make clean
make
cd ../../../

# lm2sino5d
cd lm2sino5d
make clean
make
cd ..

# add_scatter_factors
cd add_scatter_factors

make clean
make
cd ..

echo "================================================="
echo "--> All executables successfully installed. Done."

