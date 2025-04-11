##########################################################
#
#	This script sets up links so that fastTest
#	uses the littleendian data for the object files
#	(in the object subdirectory) and for the baseline
#	data (in DHCI, PET2d, PET3d, SPECT, SPECTcone,
#	SPECTfan).
#
##########################################################
cd DHCI
rm -f DHCI.base.ct
rm -f DHCI.base.wt
rm -f DHCI.base.wtsq
ln -s -f DHCI.base.ct.littleend DHCI.base.ct
ln -s -f DHCI.base.wt.littleend DHCI.base.wt
ln -s -f DHCI.base.wtsq.littleend DHCI.base.wtsq
cd ..

cd PET2d
rm -f PET2d.base.ct
rm -f PET2d.base.wt
rm -f PET2d.base.wtsq
ln -s -f PET2d.base.ct.littleend PET2d.base.ct
ln -s -f PET2d.base.wt.littleend PET2d.base.wt
ln -s -f PET2d.base.wtsq.littleend PET2d.base.wtsq
cd ..

cd PET3d
rm -f PET3d.base.ct
rm -f PET3d.base.wt
rm -f PET3d.base.wtsq
ln -s -f PET3d.base.ct.littleend PET3d.base.ct
ln -s -f PET3d.base.wt.littleend PET3d.base.wt
ln -s -f PET3d.base.wtsq.littleend PET3d.base.wtsq
cd ..

cd PETblock
rm -f PETblock.base.ct
rm -f PETblock.base.wt
rm -f PETblock.base.wtsq
ln -s -f PETblock.base.ct.littleend PETblock.base.ct
ln -s -f PETblock.base.wt.littleend PETblock.base.wt
ln -s -f PETblock.base.wtsq.littleend PETblock.base.wtsq
cd ..

cd randPET
rm -f randPET.base.ct
rm -f randPET.base.wt
rm -f randPET.base.wtsq
ln -s -f randPET.base.ct.littleend randPET.base.ct
ln -s -f randPET.base.wt.littleend randPET.base.wt
ln -s -f randPET.base.wtsq.littleend randPET.base.wtsq
cd ..

cd SPECT
rm -f SPECT.base.ct
rm -f SPECT.base.wt
rm -f SPECT.base.wtsq
ln -s -f SPECT.base.ct.littleend SPECT.base.ct
ln -s -f SPECT.base.wt.littleend SPECT.base.wt
ln -s -f SPECT.base.wtsq.littleend SPECT.base.wtsq
cd ..

cd SPECTcone
rm -f SPECTcone.base.ct
rm -f SPECTcone.base.wt
rm -f SPECTcone.base.wtsq
ln -s -f SPECTcone.base.ct.littleend SPECTcone.base.ct
ln -s -f SPECTcone.base.wt.littleend SPECTcone.base.wt
ln -s -f SPECTcone.base.wtsq.littleend SPECTcone.base.wtsq
cd ..

cd SPECTfan
rm -f SPECTfan.base.ct
rm -f SPECTfan.base.wt
rm -f SPECTfan.base.wtsq
ln -s -f SPECTfan.base.ct.littleend SPECTfan.base.ct
ln -s -f SPECTfan.base.wt.littleend SPECTfan.base.wt
ln -s -f SPECTfan.base.wtsq.littleend SPECTfan.base.wtsq
cd ..

cd TOFPET
rm -f TOFPET.base.ct
rm -f TOFPET.base.wt
rm -f TOFPET.base.wtsq
ln -s -f TOFPET.base.ct.littleend TOFPET.base.ct
ln -s -f TOFPET.base.wt.littleend TOFPET.base.wt
ln -s -f TOFPET.base.wtsq.littleend TOFPET.base.wtsq
cd ..

cd object
rm -f fastTest.act_indexes
rm -f fastTest20.act_indexes
rm -f fastTest.att_indexes
rm -f PETblock.act_indexes
rm -f PETblock.att_indexes
ln -s -f fastTest.act_indexes.littleend fastTest.act_indexes
ln -s -f fastTest20.act_indexes.littleend fastTest20.act_indexes
ln -s -f fastTest.att_indexes.littleend fastTest.att_indexes
ln -s -f PETblock.act_indexes.littleend PETblock.act_indexes
ln -s -f PETblock.att_indexes.littleend PETblock.att_indexes
cd ..

