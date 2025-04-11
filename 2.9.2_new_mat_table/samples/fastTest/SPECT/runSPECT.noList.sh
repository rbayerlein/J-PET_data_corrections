echo "Starting SPECT simulation (with no list mode) at:"
date

echo
echo "For current status of the simulation see SPECT/SPECT.noList.phgout"

nice ../../../bin/phg SPECT.noList.phgin > SPECT.noList.phgout

echo
echo "Ending SPECT simulation at:"
date
echo "---------------------------------"
echo
