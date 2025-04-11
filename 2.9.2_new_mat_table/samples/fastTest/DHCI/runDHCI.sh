echo "Starting DHCI simulation at:"
date

echo
echo "For current status of the simulation see DHCI/DHCI.phgout"

nice ../../../bin/phg DHCI.phgin > DHCI.phgout

echo
echo "Starting DHCI history file processing at:"
date

echo
echo "For current status of the simulation see DHCI/DHCI.histout"

nice ../../../bin/bin -c DHCI.histin > DHCI.histout

rm DHCI.colhist

echo
echo "Ending DHCI simulation at:"
date
echo "---------------------------------"
echo

