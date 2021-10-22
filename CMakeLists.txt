cmake_minimum_required(VERSION 3.10)
#set(CMAKE_CUDA_COMPILER "/usr/local/cuda-10.0/bin/nvcc") #Sometimes one needs to set the path to the compiler
include(GoogleTest)

project(DynamicCurve LANGUAGES CXX)
set(LIB MPMD)
set(APP_EXEC MPMD_exec)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} "-Wall -O3 -fopenmp")
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

IF(${CUDA_ON} MATCHES ON)
  enable_language(CUDA)
  file(GLOB MPMDExec "exec/main.cu")
  add_definitions(-DCUDA_ON)
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} -std=c++11 -O3 -arch=sm_75 -DCUDA_ON --default-stream per-thread --compiler-options -Wall")
  set(CUDA_SEPARABLE_COMPILATION ON)
  set(SRC "srcCU")
  set(SRCF "${SRC}/*.cu")
ELSE()
  file(GLOB MPMDExec "exec/main.cc")
  set(SRC "srcCC")
  set(SRCF "${SRC}/*.cc")
ENDIF()

message(STATUS "Compiling C++ with: " -std=c++${CMAKE_CXX_STANDARD} ${CMAKE_CXX_FLAGS})
message(STATUS "Compiling CUDA with: " ${CMAKE_CUDA_FLAGS})

include_directories(${SRC}/include include)

file(GLOB MPMDSrc ${SRCF})

message(STATUS "${MPMDSrc}")

add_library(${LIB} STATIC ${MPMDSrc})

#add_executable(${APP_SRC} ${MPMDSrc})
add_executable(${APP_EXEC} ${MPMDExec})

target_link_libraries(${APP_EXEC} ${LIB})
target_include_directories(${APP_EXEC} PRIVATE src/include/)

#add_subdirectory(test)

IF(${TEST} MATCHES "GTEST")
    #GTEST_TEST
    add_definitions(-DGTEST)
    add_subdirectory(ExternalLibs/googletest)
    enable_testing()

    add_executable(ConfigurationTest test/configurationTest.cc)
    add_executable(DubinsTest test/dubinsTest.cc)
    add_executable(ClothoidG1Test test/clothoidG1Test.cc)

    target_link_libraries(ConfigurationTest gtest gtest_main)
    target_link_libraries(DubinsTest gtest gtest_main ${LIB})
    target_link_libraries(ClothoidG1Test gtest gtest_main ${LIB})
    target_include_directories(ConfigurationTest PRIVATE ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
    target_include_directories(DubinsTest PRIVATE ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
    target_include_directories(ClothoidG1Test PRIVATE ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})

    gtest_discover_tests(ConfigurationTest)
    gtest_discover_tests(DubinsTest)
    gtest_discover_tests(ClothoidG1Test)

ELSEIF(${TEST} MATCHES "BOOST")
    #BOOST_TEST
    add_definitions(-DBOOST)
    add_definitions(-DBOOST_TEST_DYN_LINK)
    
    enable_testing()
    find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
    set(Boost_USE_MULTITHREADED TRUE)
    
    add_executable(ConfigurationTest test/configurationTest.cc)
    add_executable(ClothoidG1Test test/clothoidG1Test.cc)
    add_executable(DubinsTest test/dubinsTest.cc)
    target_link_libraries(ConfigurationTest ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    target_link_libraries(ClothoidG1Test ${LIB} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
    target_link_libraries(DubinsTest ${LIB} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

    add_test(NAME ConfigurationTest COMMAND ConfigurationTest)
    add_test(NAME ClothoidG1Test COMMAND ClothoidG1Test)
    add_test(NAME DubinsTest COMMAND DubinsTest)
ENDIF()
