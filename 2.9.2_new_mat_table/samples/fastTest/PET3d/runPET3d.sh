echo "Starting PET3d full simulation at:"
date

echo
echo "For current status of the simulation see PET3d/PET3d.phgout"

nice ../../../bin/phg PET3d.phgin > PET3d.phgout

echo
echo "Starting PET3d history file processing at:"
date

echo
echo "For current status of the simulation see PET3d/PET3d.histout"

nice ../../../bin/bin -d PET3d.histin > PET3d.histout

rm PET3d.dethist

echo
echo "Ending PET3d simulation at:"
date
echo "---------------------------------"
echo

