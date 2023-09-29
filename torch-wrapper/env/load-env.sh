#!/bin/bash

module load intel-compiler/2021.5.0

export SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
export PATH_TO_LIBTORCH=$SCRIPT_DIR/libtorch


mkdir $SCRIPT_DIR/../../build
#cd ../../build
#cmake -DCMAKE_PREFIX_PATH=$PATH_TO_LIBTORCH 


#module load python3/3.9.2
#module load netcdf/4.8.0p
#module load cdo/1.9.10
#module load gcc/11.1.0
#module load openmpi/4.1.1
#module load ncl/6.6.2
#module load nco/4.9.2
#module load ncview/2.1.7
#module load proj/8.1.1
#module load geos/3.8.0
#module load gdal/3.0.2

#python3 -m venv install env
#source env/bin/activate
#python3 -m pip install -r requirements.txt


#source env/bin/activate


