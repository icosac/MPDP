cmake_minimum_required(VERSION 3.14)
project(TorchInference)

set(CMAKE_CXX_STANDARD 17)

include(FetchContent)

# Fetch Eigen
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4
)
FetchContent_MakeAvailable(eigen)

FetchContent_Declare(
    Torch
    # URL https://download.pytorch.org/libtorch/cpu/libtorch-win-shared-with-deps-2.0.0%2Bcpu.zip
    URL https://download.pytorch.org/libtorch/cpu/libtorch-macos-arm64-2.6.0.zip
    # URL https://download.pytorch.org/libtorch/cpu/libtorch-win-shared-with-deps-debug-2.0.0%2Bcpu.zip
    # OVERRIDE_FIND_PACKAGE
)
FetchContent_MakeAvailable(Torch)
message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
list(APPEND CMAKE_PREFIX_PATH "${torch_SOURCE_DIR}")
message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")
set(CMAKE_PREFIX_PATH "${torch_SOURCE_DIR}")
message(STATUS "CMAKE_PREFIX_PATH: ${CMAKE_PREFIX_PATH}")

find_package(Torch REQUIRED)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${TORCH_CXX_FLAGS}")

add_executable(inference inference.cpp)
target_link_libraries(inference "${TORCH_LIBRARIES}" Eigen3::Eigen)
set_property(TARGET inference PROPERTY CXX_STANDARD 17)
