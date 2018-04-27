#!/usr/bin/env bash
mkdir build
cd build
cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CC_COMPILER=gcc -DCMAKE_BUILD_TYPE=Debug
cmake --build .
ctest -VV
cd ..
