echo "COMPARISON OF SPECT DATA TO BASELINE (base) DATA"
echo "____________________________________________"
echo
echo "t-test comparing SPECT run to baseline data:"
nice ../../../bin/ttest < SPECT.ttestparms > SPECT.ttest.results
grep "Minimum T-Test" SPECT.ttest.results
grep "Maximum T-Test" SPECT.ttest.results
echo
grep "0.0 <= |T-Test|" SPECT.ttest.results
echo "____________________________________________"
echo "____________________________________________"
echo
echo
