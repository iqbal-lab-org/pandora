#!/bin/bash
set -evu

export CMAKE_OPTIONS="-DCMAKE_BUILD_TYPE=$BUILD_TYPE "

echo "$CMAKE_OPTIONS"
mkdir build
cd build || exit 1
{ cmake "$CMAKE_OPTIONS" ..  && make -j4; } || exit 1
export PATH="${PWD}:${PATH}"
./pandora --help
ctest -V || exit 1
