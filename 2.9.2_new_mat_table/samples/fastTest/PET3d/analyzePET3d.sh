echo "COMPARISON OF 3D PET (PET3d) DATA TO BASELINE (base) DATA"
echo "____________________________________________"
echo
echo "t-test comparing PET3d run to baseline data:"
nice ../../../bin/ttest < PET3d.ttestparms > PET3d.ttest.results
grep "Minimum T-Test" PET3d.ttest.results
grep "Maximum T-Test" PET3d.ttest.results
echo
grep "0.0 <= |T-Test|" PET3d.ttest.results
echo "____________________________________________"
echo
echo "t-test comparing PET3d history file run to baseline data:"
nice ../../../bin/ttest < PET3d.hist.ttestparms > PET3d.hist.ttest.results
grep "Minimum T-Test" PET3d.hist.ttest.results
grep "Maximum T-Test" PET3d.hist.ttest.results
echo
grep "0.0 <= |T-Test|" PET3d.hist.ttest.results
echo "____________________________________________"
echo "____________________________________________"
