echo "Starting SPECT fanbeam (SPECTfan) base simulation at:"
date

echo
echo "For current status of the simulation see SPECTfan/.baseSim/SPECTfan.base.phgout"

nice ../../../../../2.6.2.6/bin/phg SPECTfan.base.phgin > SPECTfan.base.phgout

echo
echo "Ending SPECT fanbeam simulation at:"
date
echo "---------------------------------"
echo
