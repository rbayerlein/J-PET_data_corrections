echo "COMPARISON OF CONE-BEAM SPECT (SPECTcone) DATA TO BASELINE (base) DATA"
echo "____________________________________________"
echo
echo "t-test comparing SPECTcone run to baseline data:"
nice ../../../bin/ttest < SPECTcone.ttestparms > SPECTcone.ttest.results
grep "Minimum T-Test" SPECTcone.ttest.results
grep "Maximum T-Test" SPECTcone.ttest.results
echo
grep "0.0 <= |T-Test|" SPECTcone.ttest.results
echo "____________________________________________"
echo "____________________________________________"
echo
echo
