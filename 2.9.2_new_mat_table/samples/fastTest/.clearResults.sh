##################################################
#
#
#       This script clears the results from
#       previous fastTest runs.  This reduces
#       the size of the distribution package
#       and insures that users will not
#       accidentally run the results script
#       and compare the base results with
#       UW results instead of their own.
#
#
##################################################

cd DHCI
rm -f DHCI.ct
rm -f DHCI.wt
rm -f DHCI.wtsq
rm -f DHCI.hist.ct
rm -f DHCI.hist.wt
rm -f DHCI.hist.wtsq
rm -f DHCI.colhist
rm -f DHCI.histout DHCI.noList.phgout
rm -f DHCI.phgout runDHCI.out
cd ..

cd PET2d
rm -f PET2d.ct
rm -f PET2d.wt
rm -f PET2d.wtsq
rm -f PET2d.histout PET2d.phgout runPET2d.out
cd ..

cd PET3d
rm -f PET3d.ct
rm -f PET3d.wt
rm -f PET3d.wtsq
rm -f PET3d.hist.ct
rm -f PET3d.hist.wt
rm -f PET3d.hist.wtsq
rm -f PET3d.dethist
rm -f PET3d.histout PET3d.phgout
rm -f PET3d.noList.phgout runPET3d.out
cd ..

cd PETblock
rm -f PETblock.ct
rm -f PETblock.wt
rm -f PETblock.wtsq
rm -f PETblock.phgout runPETblock.out
cd ..

cd randPET
rm -f randPET.ct
rm -f randPET.wt
rm -f randPET.wtsq
rm -f randPET.dethist
rm -f randPET.dethist.srtd
rm -f randPET.arhist
rm -f randPET.arout randPET.binout
rm -f randPET.phgout randPET.tsout
cd .baseSim
rm -f randPET.base.dethist
rm -f randPET.base.dethist.srtd
rm -f randPET.base.arhist
cd ../..

cd SPECT
rm -f SPECT.ct
rm -f SPECT.wt
rm -f SPECT.wtsq
rm -f SPECT.hist.ct
rm -f SPECT.hist.wt
rm -f SPECT.hist.wtsq
rm -f SPECT.phghist
rm -f SPECT.histout SPECT.noList.phgout
rm -f SPECT.phgout
cd ..

cd SPECTcone
rm -f SPECTcone.ct
rm -f SPECTcone.wt
rm -f SPECTcone.wtsq
rm -f SPECTcone.phgout
cd ..

cd SPECTfan
rm -f SPECTfan.ct
rm -f SPECTfan.wt
rm -f SPECTfan.wtsq
rm -f SPECTfan.phgout
cd ..

cd TOFPET
rm -f TOFPET.ct
rm -f TOFPET.wt
rm -f TOFPET.wtsq
rm -f TOFPET.phgout
cd ..

rm -f runFast.noList.out runFast.out
