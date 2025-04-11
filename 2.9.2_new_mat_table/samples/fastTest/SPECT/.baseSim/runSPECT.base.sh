echo "Starting SPECT base simulation at:"
date

echo
echo "For current status of the simulation see SPECT/.baseSim/SPECT.base.phgout"

nice ../../../../../2.6.2.6/bin/phg SPECT.base.phgin > SPECT.base.phgout

echo
echo "Ending SPECT simulation at:"
date
echo "---------------------------------"
echo
