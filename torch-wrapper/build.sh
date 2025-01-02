#!/bin/bash

echo "Make sure to run the install script first (see envs/nci-gadi/)"
source ./env/load-env.sh

#export PATH_TO_LIBTORCH="/g/data/up6/daf561/repos/torchclim/torch-wrapper/env/libtorch/"
cd build
cmake --build . --target clean

#cmake -DCMAKE_PREFIX_PATH=$PATH_TO_LIBTORCH ..
cmake -DONNXRUNTIME_ROOTDIR=$PATH_TO_ONNX -DCMAKE_BUILD_TYPE=DEBUG -DCFLAGS=-trace ..

#cmake --build . --config Release
cmake --build . --config Debug


