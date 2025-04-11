#!/bin/bash

# user-configurable parameters

num_threads=14
coh_tables=/home/rbayerlein/Code/Recon/scatter_correction/simset/2.9.2_original/phg.data/coh.tables

# main code

echo "./mat start time: $(date)"
parallel -j ${num_threads} './mat < {}' ::: $coh_tables/water*.dat
wait
echo "./mat end time: $(date)"
echo "./calc_ad start time: $(date)"
parallel -j ${num_threads} './calc_ad < {} > {.}.ad' ::: $coh_tables/water*.dat
wait
echo "./calc_ad end time: $(date)"
echo "./calc_comp_ad start time: $(date)"
parallel -j ${num_threads} './calc_comp_ad < {} > {.}.comp_ad' ::: $coh_tables/water*.dat
wait
echo "./calc_comp_ad end time: $(date)"

