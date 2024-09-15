# List of available options

option(COMPILE_CXX_EXEC "Compile executable with C++" ON)
option(CUDA_ON "Compile with CUDA" OFF)
option(COMPILE_CUDA_EXEC "Compile executable with CUDA (only available if CUDA_ON is set)" OFF)
option(COMPILE_EXAMPLES "Compile examples" ON)
option(COMPILE_TEST "Compile tests" OFF)
option(TEST "Choose the test framework [GTEST/BOOST]" GTEST)