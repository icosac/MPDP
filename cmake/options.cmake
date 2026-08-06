# List of available options

option(COMPILE_CXX_EXEC "Compile executable with C++" ON)
option(CUDA_ON "Compile with CUDA" OFF)
option(COMPILE_CUDA_EXEC "Compile executable with CUDA (only available if CUDA_ON is set)" OFF)
option(COMPILE_EXAMPLES "Compile examples" ON)
option(COMPILE_TESTS "Compile tests" OFF)
option(MPDP_OPENMP "Link the C++ library against OpenMP (the DP has a parallel for)" OFF)
set(TEST "GTEST" CACHE STRING "Choose the test framework [GTEST/BOOST]")
