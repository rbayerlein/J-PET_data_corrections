echo "Starting SPECT fanbeam (SPECTfan) simulation at:"
date

echo
echo "For current status of the simulation see SPECTfan/SPECTfan.phgout"

nice ../../../bin/phg SPECTfan.phgin > SPECTfan.phgout

echo
echo "Ending SPECT fanbeam simulation at:"
date
echo "---------------------------------"
echo
