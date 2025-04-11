echo "COMPARISON OF FAN-BEAM SPECT (SPECTfan) DATA TO BASELINE (base) DATA"
echo "____________________________________________"
echo
echo "t-test comparing SPECTfan run to baseline data:"
nice ../../../bin/ttest < SPECTfan.ttestparms > SPECTfan.ttest.results
grep "Minimum T-Test" SPECTfan.ttest.results
grep "Maximum T-Test" SPECTfan.ttest.results
echo
grep "0.0 <= |T-Test|" SPECTfan.ttest.results
echo "____________________________________________"
echo "____________________________________________"
echo
echo
