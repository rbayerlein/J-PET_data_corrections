echo "Starting PET3d base simulation at:"
date

echo
echo "For current status of the simulation see PET3d/.baseSim/PET3d.base.phgout"

nice ../../../../../2.6.2.6/bin/phg PET3d.base.phgin > PET3d.base.phgout

echo
echo "Ending PET3d simulation at:"
date
echo "---------------------------------"
echo

