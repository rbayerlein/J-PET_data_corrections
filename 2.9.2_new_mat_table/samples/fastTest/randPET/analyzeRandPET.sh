echo "COMPARISON OF PET DATA WITH RANDOMS (randPET) TO BASELINE (base) DATA"
echo "____________________________________________"
echo
echo "t-test comparing randPET run to baseline data:"
nice ../../../bin/ttest < randPET.ttestparms > randPET.ttest.results
grep "Minimum T-Test" randPET.ttest.results
grep "Maximum T-Test" randPET.ttest.results
echo
grep "0.0 <= |T-Test|" randPET.ttest.results
echo "____________________________________________"
echo "____________________________________________"
