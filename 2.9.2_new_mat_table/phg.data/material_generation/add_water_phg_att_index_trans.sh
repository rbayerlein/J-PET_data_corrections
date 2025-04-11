#!/bin/bash

fname_in=/home/rbayerlein/Code/Recon/scatter_correction/simset/2.9.2_original/phg.data/phg_att_index_trans

fname_in_orig=${fname_in}_orig
cp $fname_in $fname_in_orig

n_mat=2201
start=1

if echo "baa7130bdd559150f8565ca35dd24003 ${fname_in}" | md5sum -c ; then
  echo "Modifying file..."
  # remove materials #100 to #255
  head -n 110 ${fname_in} > ${fname_in}.tmp # create tmp file for materials #0 to #99
  mv ${fname_in}.tmp ${fname_in} # overwrite original file
  for ((i = $start; i<=$n_mat; i+=1)); do
    echo "$((i+99))	$((i+99))" >> ${fname_in}
  done
  echo "Done."
else
  echo "File modified previously. Skipping..."
fi
