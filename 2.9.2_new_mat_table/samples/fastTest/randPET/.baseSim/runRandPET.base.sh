echo "Starting randPET base simulation at:"
date

echo
echo "For current status of the simulation see randPET/.baseSim/randPET.base.phgout"

nice ../../../../../2.9/bin/phg randPET.base.phgin > randPET.base.phgout

echo
echo "Time-sorting list-mode data at:"
date

echo
echo "For current status of the simulation see randPET/.baseSim/randPET.base.tsout"

nice ../../../../../2.9/bin/timesort -d randPET.base.tsparms > randPET.base.tsout

echo
echo "Adding randoms to list-mode data at"
date

echo
echo "For current status of the simulation see randPET/.baseSim/randPET.base.arout"

nice ../../../../../2.9/bin/addrandoms randPET.base.arparms > randPET.base.arout

echo
echo "Binning list-mode data at:"
date

echo
echo "For current status of the simulation see randPET/.baseSim/randPET.base.binout"

nice ../../../../../2.9/bin/bin -d randPET.base.binin > randPET.base.binout

echo
echo "Ending randPET base simulation at:"
date
echo "---------------------------------"
echo

