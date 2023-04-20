#!/bin/bash

echo "Make sure to run the install script first (see envs/nci-gadi/)"
source ./envs/nci-gadi/load-env.sh

cd build
cmake -DCMAKE_PREFIX_PATH=$PATH_TO_LIBTORCH ..
cmake --build . --config Release


