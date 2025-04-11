#!/bin/bash

fname_in=phg_act_table

n_mat=22100

echo "Creating file..."
echo ${n_mat} > ${fname_in} # add the number of materials in the table
for i in $(seq 0 $((n_mat-1))); do
  echo "0.000000" >> ${fname_in}
done
echo "Done."

