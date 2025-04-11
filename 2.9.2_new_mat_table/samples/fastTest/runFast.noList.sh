echo
echo "This macro runs eight simulations designed to provide a fast"
echo "test of many major SimSET features.  The macro is a space-saving"
echo "version of runFast.sh that does not test list-mode features."
echo "The simulations are 3d PET (PET3d), 2d PET (PET2d), block-detector"
echo "PET (PETblock), time-of-flight PET (TOFPET), dual-head"
echo "coincidence imaging (DHCI), SPECT (SPECT), cone-beam SPECT"
echo "(SPECTcone), and fan-beam SPECT (SPECTfan).  Note that the"
echo "randPET simulation is skipped as the randoms simulation requires"
echo "the use of a list-mode file.  Further details about each of"
echo "the simulations are available in the %%ReadMe%% files in each"
echo "subdirectory."
echo

cd PET3d
./runPET3d.noList.sh

cd ../PET2d
./runPET2d.sh

cd ../PETblock
./runPETblock.sh

cd ../TOFPET
./runTOFPET.sh

cd ../DHCI
./runDHCI.noList.sh

cd ../SPECT
./runSPECT.noList.sh

cd ../SPECTfan
./runSPECTfan.sh

cd ../SPECTcone
./runSPECTcone.sh

cd ..

echo
echo "runFast.noList is complete.  An overview of the test results"
echo "can be viewed using resultsFast.noList.sh."
echo
