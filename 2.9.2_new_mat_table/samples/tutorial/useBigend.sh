##########################################################
#
#       This script sets up links so that the user
#       function examples use big endian
#       data for their activity and attenuation objects.
#
##########################################################


rm -f tutor1.act_index
rm -f tutor1.att_index
rm -f tutor2.att_index
ln -s -f tutor1.act_index.bigend tutor1.act_index
ln -s -f tutor1.att_index.bigend tutor1.att_index
ln -s -f tutor2.att_index.bigend tutor2.att_index

