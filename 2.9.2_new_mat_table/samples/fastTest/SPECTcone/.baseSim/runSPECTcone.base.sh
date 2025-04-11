echo "Starting SPECT conebeam (SPECTcone) base simulation at:"
date

echo
echo "For current status of the simulation see SPECTcone/.baseSim/SPECTcone.base.phgout"

nice ../../../../../2.6.2.6/bin/phg SPECTcone.base.phgin > SPECTcone.base.phgout

echo
echo "Ending SPECT conebeam simulation at:"
date
echo "---------------------------------"
echo
