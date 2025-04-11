echo "Starting PET2d full simulation at:"
date

echo
echo "For current status of the simulation see PET2d/PET2d.phgout"

nice ../../../bin/phg PET2d.phgin > PET2d.phgout

echo
echo "Ending PET2d simulation at:"
date
echo "---------------------------------"
echo
