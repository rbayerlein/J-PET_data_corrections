echo "COMPARISON OF TOFPET (TOFPET) DATA TO BASELINE (base) DATA"
echo "____________________________________________"
echo
echo "t-test comparing TOFPET run to baseline data:"
nice ../../../bin/ttest < TOFPET.ttestparms > TOFPET.ttest.results
grep "Minimum T-Test" TOFPET.ttest.results
grep "Maximum T-Test" TOFPET.ttest.results
echo
grep "0.0 <= |T-Test|" TOFPET.ttest.results
echo "____________________________________________"
echo "____________________________________________"
echo
echo
