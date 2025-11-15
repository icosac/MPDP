cmake_minimum_required(VERSION 3.16)
project(PolynomialTools LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

# Fetch Eigen so the target builds even if Eigen is not preinstalled.
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4
)
FetchContent_MakeAvailable(eigen)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)

add_library(polynomials_lib STATIC
    polynomials/man01.cc
    polynomials/man07.cc
    polynomials/man10.cc
    polynomials/man11.cc
    polynomials/man12.cc
    polynomials/man13.cc
)
target_link_libraries(polynomials_lib PRIVATE Eigen3::Eigen)

add_executable(polynomials polynomials)
target_link_libraries(polynomials PRIVATE Eigen3::Eigen)
target_link_libraries(polynomials PRIVATE polynomials_lib)
