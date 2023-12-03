#!/bin/bash

make gen_svd
mkdir -p build

cd build
rm -rf *
cmake ..
make -j6 LIBRARY_PATH=/opt/homebrew/Cellar/icu4c/73.2/lib
cd ..

