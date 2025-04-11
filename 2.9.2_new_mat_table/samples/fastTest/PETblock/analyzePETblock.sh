echo "COMPARISON OF PET BLOCK (PETblock) DATA TO BASELINE (base) DATA"
echo "____________________________________________"
echo
echo "t-test comparing PETblock run to baseline data:"
nice ../../../bin/ttest < PETblock.ttestparms > PETblock.ttest.results
grep "Minimum T-Test" PETblock.ttest.results
grep "Maximum T-Test" PETblock.ttest.results
echo
grep "0.0 <= |T-Test|" PETblock.ttest.results
echo "____________________________________________"
echo "____________________________________________"
