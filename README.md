# Framework to compute multipoint using dynamic programming, even GPU enhanced.

## Compile

### CMake 
The code can be compiled with cmake. 

If not already present, create the directory which will contain the compiled code:

```shell
mkdir build
```

Check the version of cmake:

```shell
cmake --version
```

If the version is less than 3.20:

- <3.20 then use these commands to generate the code and compile it:
```shell
cmake -DCMAKE_TOOLCHAIN_FILE=CmakeUnix.cmake -B build -S .
cmake --build build -- -j4
```
- \>=3.20 then use these commands to generate the code and compile it:
```shell
cmake --presets default -B build -S .
cmake --build build -- -j4
```

From my understanding, you do not have a CUDA capable device, in case you were, please use the following commands to 
generate and compile the code for CUDA:

- <3.20 then use these commands to generate the code and compile it:
```shell
cmake -DCMAKE_TOOLCHAIN_FILE=CmakeUnixCUDA.cmake -B build -S .
cmake --build build -- -j4
```
- \>=3.20 then use these commands to generate the code and compile it:
```shell
cmake --presets cuda -B build -S .
cmake --build build -- -j4
```

Finally, you can find the executable insider the `build` directory. If you compiled:

- for C++, then it will be named `MPMDCC_exec`,
- for CUDA, then it will be named `MPMDCU_exec`
