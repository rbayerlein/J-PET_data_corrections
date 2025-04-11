#!/bin/bash

dir_in=/home/rbayerlein/Code/Recon/scatter_correction/simset/2.9.2_original/phg.data
fname_in=${dir_in}/phg_att_table

fname_in_orig=${fname_in}_orig
cp $fname_in $fname_in_orig

n_mat=22001
start=00001

if echo "ed6dc01d256e35c0c0706815a8e3acec ${fname_in}" | md5sum -c ; then
  echo "Modifying file..."
  sed -i "1s/.*/$((n_mat+100))/" ${fname_in} # change the number of materials in the table
  for ((i = $start; i <= $n_mat; i+=10)); do
    cat ${dir_in}/coh.tables/water$(printf %05d $i) >> ${fname_in}
  done
  echo "Done."
else
  echo "File modified previously. Skipping..."
fi
