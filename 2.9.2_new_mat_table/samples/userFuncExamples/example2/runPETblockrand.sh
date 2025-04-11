echo "Starting PETblockrand full simulation at:"
date

echo
echo "For current status of the simulation see PETblockrand.phgout"

nice bin/phg PETblockrand.phgin > PETblockrand.phgout

echo
echo "Time-sorting list-mode data at:"
date

echo
echo "For current status of the simulation see PETblockrand.tsout"

nice bin/timesort -d PETblockrand.tsparms > PETblockrand.tsout

rm PETblockrand.dethist

echo
echo "Adding randoms to list-mode data at"
date

echo
echo "For current status of the simulation see PETblockrand.arout"

nice bin/addrandoms PETblockrand.arparms > PETblockrand.arout

rm PETblockrand.dethist*

echo
echo "Binning list-mode data at:"
date

echo
echo "For current status of the simulation see PETblockrand.binout"

nice bin/bin -d PETblockrand.binin > PETblockrand.binout

rm PETblockrand.arhist

echo
echo "Ending PETblockrand simulation at:"
date
echo "---------------------------------"
echo

