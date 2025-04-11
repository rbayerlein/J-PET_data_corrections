#!/bin/sh

# This file compiles the programs the materials directory.

cc -o bin/mat src/mat.c
cc -o bin/calc_ad src/calc_ad.c -lm
cc -o bin/calc_comp_ad src/calc_comp_ad.c -lm

# create links to programs
ln -s -f bin/mat mat
ln -s -f bin/calc_ad calc_ad
ln -s -f bin/calc_comp_ad calc_comp_ad
