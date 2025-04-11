#!/bin/bash

# user-configurable parameters
# start = starting density (g/cm3) / 10000
# end = starting density (g/cm3) / 10000
# 10000 is arbitrary

start=00001 # g/cm3 / 10000
end=22001 # g/cm3 / 10000
outfolder=/home/rbayerlein/Code/Recon/scatter_correction/simset/2.9.2_original/phg.data/coh.tables

# main codei
ct=0
for ((i = start; i <= end; i+=10)); do
  echo $(printf %0.4f $(bc <<< "scale=5; $i/10000")),water$(printf %05d $i) > $outfolder/water$(printf %05d $i).dat
  echo "1,0.11189" >> $outfolder/water$(printf %05d $i).dat
  echo "8,0.88811" >> $outfolder/water$(printf %05d $i).dat
  echo "-1,-1" >> $outfolder/water$(printf %05d $i).dat
  let ct=ct+1
done
echo "number of materials produced: $ct"
echo "written to $outfolder/water..."
