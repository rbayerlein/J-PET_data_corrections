#!/bin/bash

fname_in=/home/rbayerlein/Code/Recon/scatter_correction/simset/2.9.2_new_mat_table/phg.data/materials_index.txt

fname_orig=${fname_in}_orig
cp $fname_in $fname_orig

n_mat=2201
strt=1

if echo "529ab32030a9589fc316b506f0323b3a ${fname_in}" | md5sum -c ; then
  echo "Modifying file..."
  for i in $(seq 1 $n_mat); do
	  let inc=$strt+($i-1)*10
    echo "$((i+99))		water$(printf %05d $inc)" >> ${fname_in}
  done
  echo "Done."
else
  echo "File modified previously. Skipping..."
fi
