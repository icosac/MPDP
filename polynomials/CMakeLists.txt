cmake_minimum_required(VERSION 3.16)
project(PolynomialTools LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

include(FetchContent)

message(STATUS "Getting Eigen library...")

# Fetch Eigen so the target builds even if Eigen is not preinstalled.
FetchContent_Declare(
    eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG 3.4
)
FetchContent_MakeAvailable(eigen)

find_package(Eigen3 3.4 REQUIRED NO_MODULE)
if(NOT TARGET MPDPCC AND NOT TARGET MPDPCC::MPDPCC)
    if(DEFINED MPDPCC_ROOT)
        find_package(MPDPCC REQUIRED CONFIG PATHS ${MPDPCC_ROOT} NO_DEFAULT_PATH)
    else()
        find_package(MPDPCC REQUIRED CONFIG)
    endif()
endif()


if(TARGET MPDPCC::MPDPCC)
    set(MPDPCC_LIB_TARGET MPDPCC::MPDPCC)
else()
    set(MPDPCC_LIB_TARGET MPDPCC)
endif()


add_library(polynomials_lib STATIC
    polynomials/man01.cc
    polynomials/man07.cc
    polynomials/man10.cc
    polynomials/man11.cc
    polynomials/man12.cc
    polynomials/man13.cc
    polynomials/man14.cc
    polynomials/man_unimplemented.cc
)
target_include_directories(polynomials_lib PRIVATE ${CMAKE_SOURCE_DIR}/polynomials)
target_link_libraries(polynomials_lib PRIVATE Eigen3::Eigen)
target_link_libraries(polynomials_lib PRIVATE ${MPDPCC_LIB_TARGET})

add_executable(polynomials polynomials.cc)
target_link_libraries(polynomials PRIVATE Eigen3::Eigen)
target_link_libraries(polynomials PRIVATE polynomials_lib)
target_link_libraries(polynomials PRIVATE ${MPDPCC_LIB_TARGET})
