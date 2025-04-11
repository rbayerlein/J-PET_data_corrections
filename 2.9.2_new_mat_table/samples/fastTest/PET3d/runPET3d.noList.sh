echo "Starting PET3d simulation (no list mode data) at:"
date

echo
echo "For current status of the simulation see PET3d/PET3d.noList.phgout"

nice ../../../bin/phg PET3d.noList.phgin > PET3d.noList.phgout

echo
echo "Ending PET3d simulation at:"
date
echo "---------------------------------"
echo

