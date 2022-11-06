# CMAKE generated file: DO NOT EDIT!
# Generated by "Unix Makefiles" Generator, CMake Version 3.20

# Delete rule output on recipe failure.
.DELETE_ON_ERROR:

#=============================================================================
# Special targets provided by cmake.

# Disable implicit rules so canonical targets will work.
.SUFFIXES:

# Disable VCS-based implicit rules.
% : %,v

# Disable VCS-based implicit rules.
% : RCS/%

# Disable VCS-based implicit rules.
% : RCS/%,v

# Disable VCS-based implicit rules.
% : SCCS/s.%

# Disable VCS-based implicit rules.
% : s.%

.SUFFIXES: .hpux_make_needs_suffix_list

# Command-line flag to silence nested $(MAKE).
$(VERBOSE)MAKESILENT = -s

#Suppress display of executed commands.
$(VERBOSE).SILENT:

# A target that is always out of date.
cmake_force:
.PHONY : cmake_force

#=============================================================================
# Set environment variables for the build.

# The shell in which to execute make rules.
SHELL = /bin/sh

# The CMake executable.
CMAKE_COMMAND = /opt/clion/bin/cmake/linux/bin/cmake

# The command to remove a file.
RM = /opt/clion/bin/cmake/linux/bin/cmake -E rm -f

# Escaping for special characters.
EQUALS = =

# The top-level source directory on which CMake was run.
CMAKE_SOURCE_DIR = /home/enrico/Projects/researchProject

# The top-level build directory on which CMake was run.
CMAKE_BINARY_DIR = /home/enrico/Projects/researchProject/cmake-cu-release

# Include any dependencies generated for this target.
include CMakeFiles/MPMDCU_exec.dir/depend.make
# Include any dependencies generated by the compiler for this target.
include CMakeFiles/MPMDCU_exec.dir/compiler_depend.make

# Include the progress variables for this target.
include CMakeFiles/MPMDCU_exec.dir/progress.make

# Include the compile flags for this target's objects.
include CMakeFiles/MPMDCU_exec.dir/flags.make

CMakeFiles/MPMDCU_exec.dir/exec/main.cu.o: CMakeFiles/MPMDCU_exec.dir/flags.make
CMakeFiles/MPMDCU_exec.dir/exec/main.cu.o: ../exec/main.cu
CMakeFiles/MPMDCU_exec.dir/exec/main.cu.o: CMakeFiles/MPMDCU_exec.dir/compiler_depend.ts
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --progress-dir=/home/enrico/Projects/researchProject/cmake-cu-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_1) "Building CUDA object CMakeFiles/MPMDCU_exec.dir/exec/main.cu.o"
	/opt/cuda/bin/nvcc -forward-unknown-to-host-compiler $(CUDA_DEFINES) $(CUDA_INCLUDES) $(CUDA_FLAGS) -MD -MT CMakeFiles/MPMDCU_exec.dir/exec/main.cu.o -MF CMakeFiles/MPMDCU_exec.dir/exec/main.cu.o.d -x cu -dc /home/enrico/Projects/researchProject/exec/main.cu -o CMakeFiles/MPMDCU_exec.dir/exec/main.cu.o

CMakeFiles/MPMDCU_exec.dir/exec/main.cu.i: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Preprocessing CUDA source to CMakeFiles/MPMDCU_exec.dir/exec/main.cu.i"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_PREPROCESSED_SOURCE

CMakeFiles/MPMDCU_exec.dir/exec/main.cu.s: cmake_force
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green "Compiling CUDA source to assembly CMakeFiles/MPMDCU_exec.dir/exec/main.cu.s"
	$(CMAKE_COMMAND) -E cmake_unimplemented_variable CMAKE_CUDA_CREATE_ASSEMBLY_SOURCE

# Object files for target MPMDCU_exec
MPMDCU_exec_OBJECTS = \
"CMakeFiles/MPMDCU_exec.dir/exec/main.cu.o"

# External object files for target MPMDCU_exec
MPMDCU_exec_EXTERNAL_OBJECTS =

CMakeFiles/MPMDCU_exec.dir/cmake_device_link.o: CMakeFiles/MPMDCU_exec.dir/exec/main.cu.o
CMakeFiles/MPMDCU_exec.dir/cmake_device_link.o: CMakeFiles/MPMDCU_exec.dir/build.make
CMakeFiles/MPMDCU_exec.dir/cmake_device_link.o: libMPMDCU.a
CMakeFiles/MPMDCU_exec.dir/cmake_device_link.o: CMakeFiles/MPMDCU_exec.dir/dlink.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/enrico/Projects/researchProject/cmake-cu-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_2) "Linking CUDA device code CMakeFiles/MPMDCU_exec.dir/cmake_device_link.o"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MPMDCU_exec.dir/dlink.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MPMDCU_exec.dir/build: CMakeFiles/MPMDCU_exec.dir/cmake_device_link.o
.PHONY : CMakeFiles/MPMDCU_exec.dir/build

# Object files for target MPMDCU_exec
MPMDCU_exec_OBJECTS = \
"CMakeFiles/MPMDCU_exec.dir/exec/main.cu.o"

# External object files for target MPMDCU_exec
MPMDCU_exec_EXTERNAL_OBJECTS =

MPMDCU_exec: CMakeFiles/MPMDCU_exec.dir/exec/main.cu.o
MPMDCU_exec: CMakeFiles/MPMDCU_exec.dir/build.make
MPMDCU_exec: libMPMDCU.a
MPMDCU_exec: CMakeFiles/MPMDCU_exec.dir/cmake_device_link.o
MPMDCU_exec: CMakeFiles/MPMDCU_exec.dir/link.txt
	@$(CMAKE_COMMAND) -E cmake_echo_color --switch=$(COLOR) --green --bold --progress-dir=/home/enrico/Projects/researchProject/cmake-cu-release/CMakeFiles --progress-num=$(CMAKE_PROGRESS_3) "Linking CUDA executable MPMDCU_exec"
	$(CMAKE_COMMAND) -E cmake_link_script CMakeFiles/MPMDCU_exec.dir/link.txt --verbose=$(VERBOSE)

# Rule to build all files generated by this target.
CMakeFiles/MPMDCU_exec.dir/build: MPMDCU_exec
.PHONY : CMakeFiles/MPMDCU_exec.dir/build

CMakeFiles/MPMDCU_exec.dir/clean:
	$(CMAKE_COMMAND) -P CMakeFiles/MPMDCU_exec.dir/cmake_clean.cmake
.PHONY : CMakeFiles/MPMDCU_exec.dir/clean

CMakeFiles/MPMDCU_exec.dir/depend:
	cd /home/enrico/Projects/researchProject/cmake-cu-release && $(CMAKE_COMMAND) -E cmake_depends "Unix Makefiles" /home/enrico/Projects/researchProject /home/enrico/Projects/researchProject /home/enrico/Projects/researchProject/cmake-cu-release /home/enrico/Projects/researchProject/cmake-cu-release /home/enrico/Projects/researchProject/cmake-cu-release/CMakeFiles/MPMDCU_exec.dir/DependInfo.cmake --color=$(COLOR)
.PHONY : CMakeFiles/MPMDCU_exec.dir/depend

