echo "Starting PETblock full simulation at:"
date

echo
echo "For current status of the simulation see PETblock/PETblock.phgout"

nice ../../../bin/phg PETblock.phgin > PETblock.phgout

echo
echo "Ending PETblock simulation at:"
date
echo "---------------------------------"
echo

