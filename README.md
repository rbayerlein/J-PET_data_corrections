# J-PET_data_corrections
J-PET image reconstruction framework. First version.
Required scripts for creating sinograms and adding scatter and random factors to the castor list mode files.
installation of SimSET and Castor software in one location for streamlined code package and compact execution

# REQUIRED: 
- MATLAB installed version 2020b or higher
- Downloaded the material tables for SimSET and saved them in <path/to/jpet>/2.9.2_new_mat_table/phg.data/coh.tables/
- link to the material tables: [link_to_BOX_folder](https://ucdavis.box.com/s/0namcmjj14yeqlxbxsuzl1jmo2h4jzqm)
- Downloaded the attenuation tables for SimSET and saved them in <path/to/jpet>/2.9.2_new_mat_table/phg.data/
- [link_to_BOX_folder](https://ucdavis.box.com/s/yhvu7xpbdn25dakdnonoyhqoy5f53ote)
# Author/owner:
rbayerlein@ucdavis.edu
manishdasind@gmail.com

# Date
April 2025

---

# Details
## install all executbles using the following commands:
`chmod u+x make_all_executable.sh`
`./u+x make_all_executable.sh`

## install Castor software:
go to the castor folder and open Readme. 
Follow the instructions :)

# ATTENTION:
- Castor is NOT installed when running the make_all_executables script.
- please install it youself and choose the configuration you need
