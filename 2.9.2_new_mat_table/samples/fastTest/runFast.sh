echo
echo "This macro runs nine simulations designed to provide a fast"
echo "test of many major SimSET features.  The simulations are"
echo "3d PET (PET3d), 2d PET (PET2d), block-detector PET (PETblock),"
echo "time-of-flight PET (TOFPET), a PET simulation with randoms"
echo "simulation incorporated (randPET), dual-head coincidence imaging"
echo "(DHCI), SPECT (SPECT), cone-beam SPECT (SPECTcone), and fan-"
echo "beam SPECT (SPECTfan).  Further details about each of these"
echo "simulations are available in the %%ReadMe%% files in each"
echo "subdirectory."
echo

cd PET3d
./runPET3d.sh

cd ../PET2d
./runPET2d.sh

cd ../PETblock
./runPETblock.sh

cd ../TOFPET
./runTOFPET.sh

cd ../randPET
./runRandPET.sh

cd ../DHCI
./runDHCI.sh

cd ../SPECT
./runSPECT.sh

cd ../SPECTfan
./runSPECTfan.sh

cd ../SPECTcone
./runSPECTcone.sh

cd ..

echo
echo "runFast is complete.  An overview of the test results can be"
echo "viewed using resultsFast.sh."
echo
