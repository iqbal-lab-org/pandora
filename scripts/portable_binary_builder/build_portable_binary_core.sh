#!/bin/bash

set -e

# Activate Holy Build Box environment.
source /hbb_exe/activate

set -eux

# install packages needed to compile
yum install wget git binutils-devel -y

# compile pandora
cd io
mkdir build_portable_executable
cd build_portable_executable
cmake -DPRINT_STACKTRACE=True -DCMAKE_BUILD_TYPE=Release ..
make VERBOSE=1 -j 4
ctest -VV

# verify if the binary is portable
set +e
libcheck pandora
set -e

# print ldd output for us to check the dependencies also
ldd pandora

# copy binary to host filesystem
cp pandora /io/pandora-linux-precompiled
