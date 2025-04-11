echo "Starting SPECT conebeam (SPECTcone) simulation at:"
date

echo
echo "For current status of the simulation see SPECTcone/SPECTcone.phgout"

nice ../../../bin/phg SPECTcone.phgin > SPECTcone.phgout

echo
echo "Ending SPECT conebeam simulation at:"
date
echo "---------------------------------"
echo
