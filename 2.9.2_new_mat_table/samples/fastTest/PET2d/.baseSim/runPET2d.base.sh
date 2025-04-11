echo "Starting PET2d base simulation at:"
date

echo
echo "For current status of the simulation see PET2d/.baseSim/PET2d.base.phgout"

nice ../../../../../2.6.2.6/bin/phg PET2d.base.phgin > PET2d.base.phgout

echo
echo "Ending PET2d simulation at:"
date
echo "---------------------------------"
echo
