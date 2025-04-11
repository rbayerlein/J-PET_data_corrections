#!/bin/bash

fname_in=/home/rbayerlein/Code/Recon/scatter_correction/simset/2.9.2_original/phg.data/phg_ad_files

n_mat=70

for ((i = 0; i < n_mat; i++)); do
  echo "@simset/phg.data/coh.tables/temp${i}.ad" >> ${fname_in}
done
