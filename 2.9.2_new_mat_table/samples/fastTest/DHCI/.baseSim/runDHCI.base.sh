echo "Starting DHCI base simulation at:"
date

echo
echo "For current status of the simulation see DHCI/.baseSim/DHCI.base.phgout"

nice ../../../../../2.6.2.6/bin/phg DHCI.base.phgin > DHCI.base.phgout

echo
echo "Ending DHCI simulation at:"
date
echo "---------------------------------"
echo

