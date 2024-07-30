<h1 style="text-align: center;">Multi-Point Dynamic Programming Framework (MPDP)</h1>

----------
- [Introduction](#introduction)
- [Download](#download)
- [Dependencies](#dependencies)
- [Compile](#compile)
  - [CMake](#cmake)
    - [CUDA support](#cuda-support)
  - [Makefile](#makefile)
- [Usage](#usage)
  - [Point-to-point](#PP)
  - [Multi-point](#MP)
  - [Splitting](#splitting)
- [ROS2 Integration](#ros2-integration)
- [Documentation](#documentation)
- [References](#references)
- [Contributors](#contributors)
- [Acknowledgements](#acknowledgements)
- [Licence](#licence)

----------

# Introduction

This framework uses dynamic programming to solve sampling-based multi-point interpolation. It's based on the papers
present in the References section.

The framework contains general code to interpolate a series of points by finding the angles through a precise sampling.

The two types of paths currently supported and integrated are [Markov-Dubins](#dubins) paths and [Reeds Shepp](#rs) (RS), but the general
structure of the code allows for easily implementing any new path by extending the class `Curve`.

# Download

To download the code, you can clone the repository:

```shell
git clone --recurse-submodules https://www.github.com/icosac/mpdp.git
```
The external module that is downloaded is [Google Test](https://github.com/google/googletest) for testing purposes.

Neither of the two modules are essential for the code to work. 

# Dependencies

The code is written in C++ and requires the following libraries:
- [OpenMP](https://www.openmp.org/) for parallelization;
- g++ version to support at least C++11;
- If you want to compile the code with CUDA support, you need to have a CUDA-enabled GPU and the CUDA toolkit installed.
- [asymptote](https://asymptote.sourceforge.io/) will be needed for compiling the .asy file when plotting the paths.

# Compile

## CMake 
The code can be compiled with cmake. 

If not already present, create the directory which will contain the compiled code:

```shell
mkdir build
```

Check the version of cmake:

```shell
cmake --version
```

If the version is:

- <3.20 then use these commands to generate the code and compile it:
```shell
cmake -DCMAKE_TOOLCHAIN_FILE=CmakeUnix.cmake -B build -S .
cmake --build build --parallel
``` 

- \>=3.20 then use these commands to generate the code and compile it:
```shell
cmake --presets default -B build -S .
cmake --build build --parallel
```

Other possibilities are:
- <3.20: `CMakeMSVC.cmake` for Windows compiling;
- \>=3.20 other available presets are: 

  - `Debug` for compiling with extra printings and debugger capabilities;
  - `GTEST` to compile a series of tests;
  - `ninja` to compile on Windows.

Finally, you can find the executable and the library inside the `build` directory. If you compiled:

- for C++, then they will be named `MPMDCC_exec` and `libMPDPCC.a`, respectively;
- for CUDA, then they will be named `MPMDCU_exec`, `libMPDPCU.a`, respectively.

### CUDA support

To compile with CUDA support, run the previous commands adding `Cuda` to either the cmake file or to the preset name.

- <3.20 then use these commands to generate the code and compile it:
```shell
cmake -DCMAKE_TOOLCHAIN_FILE=CMakeUnixCuda.cmake -B build -S .
cmake --build build --parallel
```
- \>=3.20 then use these commands to generate the code and compile it:
```shell
cmake --presets ReleaseCuda -B build -S .
cmake --build build --parallel
```

Notice that:
- CUDA support is available only for Unix systems.
- It's possible that you need to change the CUDA_ARCHITECTURES in the CMake file to match your GPU architecture.

## Makefile

Notice that a Makefile is available for compiling the code, but it has been deprecated and **is no longer maintained**. 
If you require it for any reason, you can start from that and make the due changes to compile the code.

# Usage

The best way to understand how to use the library is to look at the examples in the `examples` directory.

There are 3 main usage scenarios:

- Single point-to-point (PP) interpolation;
- Multi-point (MP) interpolation;
- Splitting a curve (mainly used for ROS2 applications).

In all cases, both the Dubins and the RS paths are supported.

## PP 

Simply create an object of type `Dubins` or `RS` providing the start and end configurarations (x, y, theta) and the curvature of the path. 
Then call the `solve` method to find the optimal path.

## MP

The class `DP` can be used to interpolate a series of points. You should provide:

- The configurations to interpolate through (the angle is calculated automatically and can be set to 0 initially)'.
- A vector of `bool` to specify if the angle should be calculated or kept fixed with the provided value.
- The number of discretizations for the angle.
- The number of refinement steps for the dynamic programming.
- Parameters for the path, e.g., maximum curvature.

Please refer to the documentation for more details on the arguments to pass. 

The number of discretizations `k` defines the number of samples to consider for each angle. During the first iteration, the whole space ($2\pi$) is divided into `k` values which are then tested. 
If the number of refinements is greater than 0, then a smaller space around the best angle found before is sampled by considering `k` new values. 

Please refer to (Frego et al., 2020) for more details on the algorithm.

## Splitting

Not implemented yet 

# ROS2 Integration

The code has been integrated with ROS2 offering the interpolation of the paths as a service. The code is available in the `mpdp_ros2` repository.

# Documentation

For the documentation, please visit this [link](https://icosac.github.io/mpdp/).

If the page is not available, you can generate the documentation by running Doxygen with the following command within the main directory:

```shell
doxygen
```

# References

M. Frego, P. Bevilacqua, E. Saccon, L. Palopoli and D. Fontanelli, "An Iterative Dynamic Programming Approach to the Multipoint Markov-Dubins Problem," in _IEEE Robotics and Automation Letters_, vol. 5, no. 2, pp. 2483-2490, April 2020, doi: <a href="https://www.doi.org/10.1109/LRA.2020.2972787">10.1109/LRA.2020.2972787</a>.

E. Saccon, P. Bevilacqua, D. Fontanelli, M. Frego, L. Palopoli and R. Passerone, "Robot Motion Planning: can GPUs be a Game Changer?," _2021 IEEE 45th Annual Computers, Software, and Applications Conference (COMPSAC)_, Madrid, Spain, 2021, pp. 21-30, doi: <a href="https://www.doi.org/10.1109/COMPSAC51774.2021.00015">10.1109/COMPSAC51774.2021.00015</a>.

<p id="dubins">Dubins, Lester E.. "On Curves of Minimal Length with a Constraint on Average Curvature, and with Prescribed Initial and Terminal Positions and Tangents." _American Journal of Mathematics_ 79 (1957): 497.</p>

<p id="rs">Reeds, James A. and Larry A. Shepp. "Optimal paths for a car that goes both forwards and backwards" _Pacific Journal of Mathematics_ 145 (1990): 367-393, <a href="http://dx.doi.org/10.2140/pjm.1990.145.367">10.2140/pjm.1990.145.367</a>.</p>


# Contributors

- [Enrico Saccon](enricosaccon96@gmail.com)
- [Marco Frego](marco.frego@unibz.it)

# Acknowledgements

The work is based on previous work from Marco Frego and Paolo Bevilacqua. Also, the work was largely inspired by Enrico Bertolazzi and his work on the [Clothoids path](https://github.com/ebertolazzi/Clothoids), from which the `asyplot.{cc,hh}` and the `clothoidLib.asyplot` files for plotting have been taken.

# Licence

The code is released under the AGPL3 licence. Please refer to the LICENCE file for more details.
