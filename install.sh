#!/usr/bin/env bash
mkdir build
cd build
cmake .. -DCMAKE_CXX_COMPILER=${1:-g++} -DCMAKE_CC_COMPILER=${2:-gcc} -DCMAKE_BUILD_TYPE=Debug
cmake --build .
ctest -VV
cd ..
