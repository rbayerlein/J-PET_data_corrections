##########################################################
#
#	This script sets up links so that the user
#	function examples use little endian
#	data for their activity and attenuation objects.
#
##########################################################

cd example1
rm -f PETblockrand.act_indexes
rm -f PETblockrand.att_indexes
ln -s -f PETblockrand.act_indexes.littleend PETblockrand.act_indexes
ln -s -f PETblockrand.att_indexes.littleend PETblockrand.att_indexes
cd ..

cd example2
rm -f PETblockrand.act_indexes
rm -f PETblockrand.att_indexes
ln -s -f PETblockrand.act_indexes.littleend PETblockrand.act_indexes
ln -s -f PETblockrand.att_indexes.littleend PETblockrand.att_indexes
cd ..

