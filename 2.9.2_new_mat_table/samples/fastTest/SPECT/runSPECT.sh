echo "Starting SPECT simulation at:"
date

echo
echo "For current status of the simulation see SPECT/SPECT.phgout"

nice ../../../bin/phg SPECT.phgin > SPECT.phgout

echo
echo "Starting SPECT history file processing at:"
date

echo
echo "For current status of the simulation see SPECT/SPECT.histout"

nice ../../../bin/phgbin -p SPECT.histin > SPECT.histout

rm SPECT.phghist

echo
echo "Ending SPECT simulation at:"
date
echo "---------------------------------"
echo
