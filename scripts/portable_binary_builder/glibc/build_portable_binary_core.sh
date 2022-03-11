#!/bin/bash

set -eux

# compile pandora
cd /pandora
mkdir build_portable_executable_glibc
cd build_portable_executable_glibc
cmake \
  -DHUNTER_JOBS_NUMBER=4 \
  -DCMAKE_EXE_LINKER_FLAGS="-static" \
  -DCMAKE_BUILD_TYPE=Release \
  -Dspoa_use_simde=ON \
  -Dspoa_optimize_for_native=OFF \
  ..
make VERBOSE=1 -j 4
ctest -VV

# copy binary to host filesystem
cp pandora /pandora/pandora-linux-precompiled-glibc