echo "Starting randPET full simulation at:"
date

echo
echo "For current status of the simulation see randPET/randPET.phgout"

nice ../../../bin/phg randPET.phgin > randPET.phgout

echo
echo "Time-sorting list-mode data at:"
date

echo
echo "For current status of the simulation see randPET/randPET.tsout"

nice ../../../bin/timesort -d randPET.tsparms > randPET.tsout

rm randPET.dethist

echo
echo "Adding randoms to list-mode data at"
date

echo
echo "For current status of the simulation see randPET/randPET.arout"

nice ../../../bin/addrandoms randPET.arparms > randPET.arout

rm randPET.dethist*

echo
echo "Binning list-mode data at:"
date

echo
echo "For current status of the simulation see randPET/randPET.binout"

nice ../../../bin/bin -d randPET.binin > randPET.binout

rm randPET.arhist

echo
echo "Ending randPET simulation at:"
date
echo "---------------------------------"
echo

