mkdir build
cd build
cmake .. -DCMAKE_CXX_COMPILER=g++ -DCMAKE_CC_COMPILER=gcc
cmake --build .
ctest -VV
