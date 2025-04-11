echo "COMPARISON OF DHCI DATA TO BASELINE (base) DATA"
echo "____________________________________________"
echo
echo "t-test comparing DHCI run to baseline data:"
nice ../../../bin/ttest < DHCI.ttestparms > DHCI.ttest.results
grep "Minimum T-Test" DHCI.ttest.results
grep "Maximum T-Test" DHCI.ttest.results
echo
grep "0.0 <= |T-Test|" DHCI.ttest.results
echo "____________________________________________"
echo
echo "t-test comparing DHCI history file run to baseline data:"
nice ../../../bin/ttest < DHCI.hist.ttestparms > DHCI.hist.ttest.results
grep "Minimum T-Test" DHCI.hist.ttest.results
grep "Maximum T-Test" DHCI.hist.ttest.results
echo
grep "0.0 <= |T-Test|" DHCI.hist.ttest.results
echo "____________________________________________"
echo "____________________________________________"
