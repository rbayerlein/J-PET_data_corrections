echo "Starting TOFPET base simulation at:"
date

echo
echo "For current status of the simulation see TOFPET/.baseSim/TOFPET.base.phgout"

nice ../../../../../2.9/bin/phg TOFPET.base.phgin > TOFPET.base.phgout

echo
echo "Ending TOFPET simulation at:"
date
echo "---------------------------------"
echo
