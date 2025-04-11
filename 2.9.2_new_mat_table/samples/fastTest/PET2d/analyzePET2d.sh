echo "COMPARISON OF 2D PET (PET2d) DATA TO BASELINE (base) DATA"
echo "____________________________________________"
echo
echo "t-test comparing PET2d run to baseline data:"
nice ../../../bin/ttest < PET2d.ttestparms > PET2d.ttest.results
grep "Minimum T-Test" PET2d.ttest.results
grep "Maximum T-Test" PET2d.ttest.results
echo
grep "0.0 <= |T-Test|" PET2d.ttest.results
echo "____________________________________________"
echo "____________________________________________"
echo
echo
