#!/bin/bash
mkdir -p build
cd build
make clean
cmake ..
cd ..
make -C ./build/ -j8
