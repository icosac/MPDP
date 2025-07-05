cmake_minimum_required(VERSION 3.14)
project(HNSWlibEigenExample)

# Use C++17 or higher
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# FetchContent to download HNSWlib automatically
include(FetchContent)
FetchContent_Declare(
    hnswlib
    GIT_REPOSITORY https://github.com/nmslib/hnswlib.git
    GIT_TAG master  # Use the latest master branch
)
FetchContent_MakeAvailable(hnswlib)

# FetchContent to download Eigen automatically
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG master  
)
FetchContent_MakeAvailable(eigen)

# Create an executable
add_executable(knn_train train.cpp)
add_executable(knn_load load.cpp)

message(STATUS "hnswlib_SOURCE_DIR: ${hnswlib_SOURCE_DIR}")
message(STATUS "eigen_SOURCE_DIR: ${eigen_SOURCE_DIR}")

include_directories(${eigen_SOURCE_DIR})
include_directories(${hnswlib_SOURCE_DIR})

