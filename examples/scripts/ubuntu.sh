# !/bin/sh


#Created on July 10, 2016
#author: Andrew Abi-Mansour

# -*- coding: utf8 -*-
# -------------------------------------------------------------------------
#
# Shell script for compilation/installation of LIGGGHTS + PyGran from source 
# code on Ubuntu 18.04.
#
# --------------------------------------------------------------------------
#
#   This program is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 2 of the License, or
#   (at your option) any later version.

#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.

#   You should have received a copy of the GNU General Public License
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.

# --------


echo "Updating system and installing dependencies ..."
sudo apt-get update && sudo apt-get install gcc libopenmpi-dev python3-setuptools python3-pip ipython3 git python3-matplotlib libvtk6-dev -y

echo "Installing PyGran dependencies ..."
pip3 install numpy cython --user

echo "Cloning LIGGGHTS ..."
git clone https://github.com/CFDEMproject/LIGGGHTS-PUBLIC.git

echo "Compiling LIGGGHTS ..."
cd LIGGGHTS-PUBLIC/src && make clean-all
make auto
make -f Makefile.shlib auto

echo "Downloading PyGran examples"
cd ../..
git clone https://github.com/Andrew-AbiMansour/PyGran.git
python3 setup.py install --user
cd examples/scripts/DEM/compaction/

echo "Testing PyGran: running compaction problem"
mpirun python3 compaction.py

echo "Testing done! PyGran was successfully tested on your system."
