echo "Starting PETblock base simulation at:"
date

echo
echo "For current status of the simulation see PETblock/.baseSim/PETblock.base.phgout"

nice ../../../../../2.9.1/bin/phg PETblock.base.phgin > PETblock.base.phgout

echo
echo "Ending PETblock simulation at:"
date
echo "---------------------------------"
echo

