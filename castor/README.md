CASToR repository

installation isntructions:

create folder called build

cd build

cmake ../

# if not installed yet
sudo apt install cmake-curses-gui

ccmake ../

enable option CASTOR_OMP
choose the install directory in the last line
configure: c
generate and install: g

make install

export the path to install/bin/castor-recon in the bashrc file
