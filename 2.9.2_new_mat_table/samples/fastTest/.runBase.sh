echo
echo "This macro runs nine simulations to create the baseline data"
echo "which the fastTest simulations use for comparison.  For more"
echo "information on the fastTest suite, see the %%ReadMe%% file in"
echo "this directory."
echo

cd PET3d/.baseSim
./runPET3d.base.sh

cd ../../PET2d/.baseSim
./runPET2d.base.sh

cd ../../PETblock/.baseSim
./runPETblock.base.sh

cd ../../TOFPET/.baseSim
./runTOFPET.base.sh

cd ../../randPET/.baseSim
./runRandPET.base.sh

cd ../../DHCI/.baseSim
./runDHCI.base.sh

cd ../../SPECT/.baseSim
./runSPECT.base.sh

cd ../../SPECTfan/.baseSim
./runSPECTfan.base.sh

cd ../../SPECTcone/.baseSim
./runSPECTcone.base.sh

cd ../..

echo
echo "runBase is complete."
echo
