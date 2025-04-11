#!/bin/bash

dir_in=../coh.tables
n_mat=70

for ((i = 0; i < n_mat; i++)); do
  cp -v ${dir_in}/air ${dir_in}/temp${i}
  cp -v ${dir_in}/air.ad ${dir_in}/temp${i}.ad
  cp -v ${dir_in}/air.comp_ad ${dir_in}/temp${i}.comp_ad
  cp -v ${dir_in}/air.dat ${dir_in}/temp${i}.dat
  sed -i "s/air/temp${i}/g" ${dir_in}/temp${i}*
done
