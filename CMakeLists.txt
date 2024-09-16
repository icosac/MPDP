cmake_minimum_required(VERSION 3.22 FATAL_ERROR)
include(GoogleTest)

project(MPDP LANGUAGES CXX)

# Write a compile_commands.json file with information about the compiler calls
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(MPDP_DRAW OFF)

set(THIS_IS_A_TEST ON)

include(cmake/options.cmake)

# Print options
message(STATUS "COMPILE_CXX_EXEC: ${COMPILE_CXX_EXEC}")
message(STATUS "CUDA_ON: ${CUDA_ON}")
message(STATUS "COMPILE_CUDA_EXEC: ${COMPILE_CUDA_EXEC}")
message(STATUS "COMPILE_EXAMPLES: ${COMPILE_EXAMPLES}")
message(STATUS "COMPILE_TESTS: ${COMPILE_TESTS}")

include(cmake/dependencies.cmake)
include(cmake/build_cxx.cmake)

if(CUDA_ON)
  include(cmake/build_cuda.cmake)
endif()

if(COMPILE_EXAMPLES)
  include(cmake/build_examples.cmake)
endif()

if(COMPILE_TESTS)
  include(cmake/build_tests.cmake)
endif()


# include(ExternalLibs/RSPredict/CMakeLists.txt)
# include(ExternalLibs/SVMQuadraticPredict/CMakeLists.txt)
# include(ExternalLibs/NNWideCompactPredict/CMakeLists.txt)
# include_directories(${Matlab_INCLUDE_DIRS})
# include_directories(ExternalLibs/RSPredict)
# include_directories(ExternalLibs/SVMQuadraticPredict)
# include_directories(ExternalLibs/NNWideCompactPredict)
# target_link_libraries(${APP_EXEC_CXX} PRIVATE ${LIB}) # RSPredict SVMQuadraticPredict NNWideCompactPredict)
