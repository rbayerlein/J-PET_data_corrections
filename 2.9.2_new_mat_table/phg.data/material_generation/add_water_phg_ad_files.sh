#!/bin/bash

fname_in=/home/rbayerlein/Code/Recon/scatter_correction/simset/2.9.2_original/phg.data/phg_ad_files
fname_orig=${fname_in}_orig
cp $fname_in $fname_orig

n_mat=22001
strt=1

if echo "86572881a9e8aab9370262b4a44596ad ${fname_in}" | md5sum -c ; then
  echo "Modifying file..."
  for ((i=$strt; i<=$n_mat; i+=10)); do
    echo "@simset/phg.data/coh.tables/water$(printf %05d $i).ad" >> ${fname_in}
  done
  echo "Done."
else
  echo "File modified previously. Skipping..."
fi
