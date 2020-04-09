#!/bin/bash
cp CMakeLists.txt build/
#cmake -S . -B .
cmake .
#make -C ./build/ -j8
make -j8
mv build/test0 test0
./test0