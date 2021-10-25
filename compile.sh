#!/bin/bash

if [[ -d build ]] 
then
  echo "Removing dir"
  rm -rf build
fi

mkdir build 
cp test/build_clothoid.txt build/
cp test/dubinsTest.txt build/
cd build
cmake .. -D$1 -DTEST:STRING=$2
make 
#ctest $2
