##########################################################
#
#	This script sets up links so that fastTest
#	uses the bigendian data for the object files
#	(in the object subdirectory) and for the baseline
#	data (in DHCI, PET2d, PET3d, SPECT, SPECTcone,
#	SPECTfan).
#
##########################################################
cd DHCI
rm -f DHCI.base.ct
rm -f DHCI.base.wt
rm -f DHCI.base.wtsq
ln -s -f DHCI.base.ct.bigend DHCI.base.ct
ln -s -f DHCI.base.wt.bigend DHCI.base.wt
ln -s -f DHCI.base.wtsq.bigend DHCI.base.wtsq
cd ..

cd PET2d
rm -f PET2d.base.ct
rm -f PET2d.base.wt
rm -f PET2d.base.wtsq
ln -s -f PET2d.base.ct.bigend PET2d.base.ct
ln -s -f PET2d.base.wt.bigend PET2d.base.wt
ln -s -f PET2d.base.wtsq.bigend PET2d.base.wtsq
cd ..

cd PET3d
rm -f PET3d.base.ct
rm -f PET3d.base.wt
rm -f PET3d.base.wtsq
ln -s -f PET3d.base.ct.bigend PET3d.base.ct
ln -s -f PET3d.base.wt.bigend PET3d.base.wt
ln -s -f PET3d.base.wtsq.bigend PET3d.base.wtsq
cd ..

cd PETblock
rm -f PETblock.base.ct
rm -f PETblock.base.wt
rm -f PETblock.base.wtsq
ln -s -f PETblock.base.ct.bigend PETblock.base.ct
ln -s -f PETblock.base.wt.bigend PETblock.base.wt
ln -s -f PETblock.base.wtsq.bigend PETblock.base.wtsq
cd ..

cd randPET
rm -f randPET.base.ct
rm -f randPET.base.wt
rm -f randPET.base.wtsq
ln -s -f randPET.base.ct.bigend randPET.base.ct
ln -s -f randPET.base.wt.bigend randPET.base.wt
ln -s -f randPET.base.wtsq.bigend randPET.base.wtsq
cd ..

cd SPECT
rm -f SPECT.base.ct
rm -f SPECT.base.wt
rm -f SPECT.base.wtsq
ln -s -f SPECT.base.ct.bigend SPECT.base.ct
ln -s -f SPECT.base.wt.bigend SPECT.base.wt
ln -s -f SPECT.base.wtsq.bigend SPECT.base.wtsq
cd ..

cd SPECTcone
rm -f SPECTcone.base.ct
rm -f SPECTcone.base.wt
rm -f SPECTcone.base.wtsq
ln -s -f SPECTcone.base.ct.bigend SPECTcone.base.ct
ln -s -f SPECTcone.base.wt.bigend SPECTcone.base.wt
ln -s -f SPECTcone.base.wtsq.bigend SPECTcone.base.wtsq
cd ..

cd SPECTfan
rm -f SPECTfan.base.ct
rm -f SPECTfan.base.wt
rm -f SPECTfan.base.wtsq
ln -s -f SPECTfan.base.ct.bigend SPECTfan.base.ct
ln -s -f SPECTfan.base.wt.bigend SPECTfan.base.wt
ln -s -f SPECTfan.base.wtsq.bigend SPECTfan.base.wtsq
cd ..

cd TOFPET
rm -f TOFPET.base.ct
rm -f TOFPET.base.wt
rm -f TOFPET.base.wtsq
ln -s -f TOFPET.base.ct.bigend TOFPET.base.ct
ln -s -f TOFPET.base.wt.bigend TOFPET.base.wt
ln -s -f TOFPET.base.wtsq.bigend TOFPET.base.wtsq
cd ..

cd object
rm -f fastTest.act_indexes
rm -f fastTest20.act_indexes
rm -f fastTest.att_indexes
rm -f PETblock.act_indexes
rm -f PETblock.att_indexes
ln -s -f fastTest.act_indexes.bigend fastTest.act_indexes
ln -s -f fastTest20.act_indexes.bigend fastTest20.act_indexes
ln -s -f fastTest.att_indexes.bigend fastTest.att_indexes
ln -s -f PETblock.act_indexes.bigend PETblock.act_indexes
ln -s -f PETblock.att_indexes.bigend PETblock.att_indexes
cd ..

