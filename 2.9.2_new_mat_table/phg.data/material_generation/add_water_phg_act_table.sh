#!/bin/bash

dir_in=/home/rbayerlein/Code/Recon/scatter_correction/simset/2.9.2_new_mat_table/phg.data
fname_in=${dir_in}/phg_act_table

fname_orig=${fname_in}_orig
cp $fname_in $fname_orig

n_mat=2301

if echo "67fb683097805b5055807fb22fd8542c ${fname_in}" | md5sum -c ; then
  echo "Modifying file..."
  echo ${n_mat} > ${fname_in} # change the number of materials in the table
  for ((i=0; i<=${n_mat}-1; i+=1)); do
    echo $(printf %0.6f $(bc <<< "scale=7; $i/1000000")) >> ${fname_in}
  done
  echo "Done."
else
  echo "File modified previously. Skipping..."
fi
