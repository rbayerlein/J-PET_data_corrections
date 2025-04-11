echo "Starting DHCI simulation (no list mode) at:"
date

echo
echo "For current status of the simulation see DHCI/DHCI.noList.phgout"

nice ../../../bin/phg DHCI.noList.phgin > DHCI.noList.phgout

echo
echo "Ending DHCI simulation at:"
date
echo "---------------------------------"
echo

