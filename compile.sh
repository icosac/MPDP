#!/bin/bash

#if [[ -d build ]]
#then
#  echo "Removing dir"
#  rm -rf build
#fi
#
#mkdir build
#cp test/build_clothoid.txt build/
#cp test/dubinsTest.txt build/
#cmake . --preset=$1
#cmake --build --preset=$1$2

#!/usr/bin/env bash

if [ -d build ]; then
  rm -rf build
fi

# case "$OSTYPE" in
#   solaris*) echo "SOLARIS" ;;
#   darwin*)  echo "OSX" ;;
#   linux*)   echo "LINUX" ;;
#   bsd*)     echo "BSD" ;;
#   msys*)    echo "WINDOWS" ;;
#   cygwin*)  echo "ALSO WINDOWS" ;;
#   *)        echo "nknown: $OSTYPE" ;;
# esac

if [[ $OSTYPE =~ darwin* ]]; then
  echo "Changing compilers for Mac"
  export CC=/usr/local/bin/gcc-11
  export CXX=/usr/local/bin/g++-11
fi

mkdir build
cmake -B build -S .
cmake --build build --target all -- -j
#ctest --test-dir build

./build/MPMDCC_exec
