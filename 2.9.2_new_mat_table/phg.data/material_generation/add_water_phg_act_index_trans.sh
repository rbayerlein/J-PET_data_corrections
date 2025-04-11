#!/bin/bash

fname_in=/home/rbayerlein/Code/Recon/scatter_correction/simset/2.9.2_new_mat_table/phg.data/phg_act_index_trans
fname_orig=${fname_in}_orig
cp $fname_in $fname_orig

n_mat=2201

if echo "897d185a763fc7bdfa965f131bc9e9ab ${fname_in}" | md5sum -c ; then
  echo "Modifying file..."
  # remove materials #100 to #255
  head -n 110 ${fname_in} > ${fname_in}.tmp # create tmp file for materials #0 to #99
  mv ${fname_in}.tmp ${fname_in} # overwrite original file
  for i in $(seq 1 $n_mat); do
    echo "$((i+99))	$((i+99))" >> ${fname_in}
  done
  echo "Done."
else
  echo "File modified previously. Skipping..."
fi
