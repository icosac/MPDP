#!/bin/bash

if [[ -d build ]] 
then
  echo "Removing dir"
  rm -rf build
fi

mkdir build 
cp test/build_clothoid.txt build/
cp test/dubinsTest.txt build/
cmake . --preset=$1 -DTEST:STRING=$2
cmake --build --preset=$1
#ctest $2
