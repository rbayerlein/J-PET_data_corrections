##########################################################
#
#       This script sets up links so that the user
#       function examples use big endian
#       data for their activity and attenuation objects.
#
##########################################################
rm -f example1/PETblockrand.act_indexes
rm -f example1/PETblockrand.att_indexes
rm -f example2/PETblockrand.act_indexes
rm -f example2/PETblockrand.att_indexes

cd example1
ln -s -f PETblockrand.act_indexes.bigend PETblockrand.act_indexes
ln -s -f PETblockrand.att_indexes.bigend PETblockrand.att_indexes
cd ..

cd example2
ln -s -f PETblockrand.act_indexes.bigend PETblockrand.act_indexes
ln -s -f PETblockrand.att_indexes.bigend PETblockrand.att_indexes
cd ..

