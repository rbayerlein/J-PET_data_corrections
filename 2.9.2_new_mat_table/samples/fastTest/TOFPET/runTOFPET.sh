echo "Starting TOFPET full simulation at:"
date

echo
echo "For current status of the simulation see TOFPET/TOFPET.phgout"

nice ../../../bin/phg TOFPET.phgin > TOFPET.phgout

echo
echo "Ending TOFPET simulation at:"
date
echo "---------------------------------"
echo
