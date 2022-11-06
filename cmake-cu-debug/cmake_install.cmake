# Install script for directory: /home/enrico/Projects/researchProject

# Set the install prefix
if(NOT DEFINED CMAKE_INSTALL_PREFIX)
  set(CMAKE_INSTALL_PREFIX "/usr/local")
endif()
string(REGEX REPLACE "/$" "" CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}")

# Set the install configuration name.
if(NOT DEFINED CMAKE_INSTALL_CONFIG_NAME)
  if(BUILD_TYPE)
    string(REGEX REPLACE "^[^A-Za-z0-9_]+" ""
           CMAKE_INSTALL_CONFIG_NAME "${BUILD_TYPE}")
  else()
    set(CMAKE_INSTALL_CONFIG_NAME "")
  endif()
  message(STATUS "Install configuration: \"${CMAKE_INSTALL_CONFIG_NAME}\"")
endif()

# Set the component getting installed.
if(NOT CMAKE_INSTALL_COMPONENT)
  if(COMPONENT)
    message(STATUS "Install component: \"${COMPONENT}\"")
    set(CMAKE_INSTALL_COMPONENT "${COMPONENT}")
  else()
    set(CMAKE_INSTALL_COMPONENT)
  endif()
endif()

# Install shared libraries without execute permission?
if(NOT DEFINED CMAKE_INSTALL_SO_NO_EXE)
  set(CMAKE_INSTALL_SO_NO_EXE "0")
endif()

# Is this installation the result of a crosscompile?
if(NOT DEFINED CMAKE_CROSSCOMPILING)
  set(CMAKE_CROSSCOMPILING "FALSE")
endif()

# Set default install directory permissions.
if(NOT DEFINED CMAKE_OBJDUMP)
  set(CMAKE_OBJDUMP "/usr/bin/objdump")
endif()

if("x${CMAKE_INSTALL_COMPONENT}x" STREQUAL "xUnspecifiedx" OR NOT CMAKE_INSTALL_COMPONENT)
  list(APPEND CMAKE_ABSOLUTE_DESTINATION_FILES
   "/home/enrico/Projects/researchProject/cmake-cu-debug/inc/IOSettings.hh;/home/enrico/Projects/researchProject/cmake-cu-debug/inc/IOUtils.hh;/home/enrico/Projects/researchProject/cmake-cu-debug/inc/asyplot.hh;/home/enrico/Projects/researchProject/cmake-cu-debug/inc/tests.hh;/home/enrico/Projects/researchProject/cmake-cu-debug/inc/timeperf.hh;/home/enrico/Projects/researchProject/cmake-cu-debug/inc/typedefs.hh;/home/enrico/Projects/researchProject/cmake-cu-debug/inc/utilities.hh;/home/enrico/Projects/researchProject/cmake-cu-debug/inc/configuration.cuh;/home/enrico/Projects/researchProject/cmake-cu-debug/inc/constants.cuh;/home/enrico/Projects/researchProject/cmake-cu-debug/inc/curve.cuh;/home/enrico/Projects/researchProject/cmake-cu-debug/inc/dp.cuh;/home/enrico/Projects/researchProject/cmake-cu-debug/inc/dubins.cuh;/home/enrico/Projects/researchProject/cmake-cu-debug/inc/utils.cuh")
  if(CMAKE_WARN_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(WARNING "ABSOLUTE path INSTALL DESTINATION : ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
  if(CMAKE_ERROR_ON_ABSOLUTE_INSTALL_DESTINATION)
    message(FATAL_ERROR "ABSOLUTE path INSTALL DESTINATION forbidden (by caller): ${CMAKE_ABSOLUTE_DESTINATION_FILES}")
  endif()
file(INSTALL DESTINATION "/home/enrico/Projects/researchProject/cmake-cu-debug/inc" TYPE FILE FILES
    "/home/enrico/Projects/researchProject/./include/IOSettings.hh"
    "/home/enrico/Projects/researchProject/./include/IOUtils.hh"
    "/home/enrico/Projects/researchProject/./include/asyplot.hh"
    "/home/enrico/Projects/researchProject/./include/tests.hh"
    "/home/enrico/Projects/researchProject/./include/timeperf.hh"
    "/home/enrico/Projects/researchProject/./include/typedefs.hh"
    "/home/enrico/Projects/researchProject/./include/utilities.hh"
    "/home/enrico/Projects/researchProject/srcCU/include/configuration.cuh"
    "/home/enrico/Projects/researchProject/srcCU/include/constants.cuh"
    "/home/enrico/Projects/researchProject/srcCU/include/curve.cuh"
    "/home/enrico/Projects/researchProject/srcCU/include/dp.cuh"
    "/home/enrico/Projects/researchProject/srcCU/include/dubins.cuh"
    "/home/enrico/Projects/researchProject/srcCU/include/utils.cuh"
    )
endif()

if(CMAKE_INSTALL_COMPONENT)
  set(CMAKE_INSTALL_MANIFEST "install_manifest_${CMAKE_INSTALL_COMPONENT}.txt")
else()
  set(CMAKE_INSTALL_MANIFEST "install_manifest.txt")
endif()

string(REPLACE ";" "\n" CMAKE_INSTALL_MANIFEST_CONTENT
       "${CMAKE_INSTALL_MANIFEST_FILES}")
file(WRITE "/home/enrico/Projects/researchProject/cmake-cu-debug/${CMAKE_INSTALL_MANIFEST}"
     "${CMAKE_INSTALL_MANIFEST_CONTENT}")
