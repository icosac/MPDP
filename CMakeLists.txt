cmake_minimum_required(VERSION 3.12 FATAL_ERROR)
include(GoogleTest)

project(MPDP LANGUAGES CXX)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

if(${CUDA_ON} MATCHES ON)
#  include(CMakeMSCuda.cmake)
  set(LIB "MPDPCU")
  set(APP_EXEC "MPDPCU_exec")
  #Define files to be compiled
  enable_language(CUDA)

  set(SRC "srcCU")
  set(SRCF "${SRC}/*.cu")
  file(GLOB MPDPSrc ${SRCF})
  file(GLOB MPDPExec "exec/main.cu")
  file(GLOB includeFiles "${SRC}/include/*.cuh" "./include/*.hh")

  #What should be compiled
  include_directories(${SRC}/include include)
  add_library(${LIB} ${MPDPSrc})
  add_executable(${APP_EXEC} ${MPDPExec})
  target_link_libraries(${APP_EXEC} PUBLIC ${LIB})

  #What should be installed
  install(FILES ${includeFiles} DESTINATION ${CMAKE_BINARY_DIR}/inc)

  #Set compiler options
  set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} ${CUDA_FLAGS}")
  set_target_properties(${LIB} PROPERTIES CUDA_SEPARABLE_COMPILATION "${CUDA_SEP_COMP}")
  set_target_properties(${LIB} PROPERTIES CUDA_ARCHITECTURES "${CUDA_ARCH}")
  set_target_properties(${APP_EXEC} PROPERTIES CUDA_SEPARABLE_COMPILATION "${CUDA_SEP_COMP}")
  set_target_properties(${APP_EXEC} PROPERTIES CUDA_ARCHITECTURES "${CUDA_ARCH}")

  message(FATAL "Compiling CUDA with: " ${CUDA_FLAGS})
  message(STATUS "${MPDPSrc}")
else()
  set(LIB "MPDPCC")
  set(APP_EXEC "MPDPCC_exec")
  #Define files to be compiled
  set(SRC "srcCC")
  set(SRCF "${SRC}/*.cc")
  file(GLOB MPDPSrc ${SRCF})
  file(GLOB MPDPExec "exec/main.cc")
  file(GLOB includeFiles "${SRC}/include/*.hh" "./include/*.hh")  
  
  message(STATUS "Compiling C++ with: " ${CXX_FLAGS})
  message(STATUS "${MPDPSrc}")

  #Compile 
  include_directories(${SRC}/include include)
  add_library(${LIB} STATIC ${MPDPSrc})
  add_executable(${APP_EXEC} ${MPDPExec})
  target_link_libraries(${APP_EXEC} ${LIB})
  
  #What should be installed
  install(FILES ${includeFiles} DESTINATION ${CMAKE_BINARY_DIR}/inc)
  
  find_package(OpenMP REQUIRED)
  if (OpenMP_CXX_FOUND)
    target_link_libraries(${LIB} PUBLIC OpenMP::OpenMP_CXX)
  endif()
endif()


if(DEFINED TEST)
  set(TEST_NAMES
      configurationTest
      dubinsTest
      clothoidG1Test
      RSTest)

  if(${TEST} MATCHES "GTEST")
    #GTEST_TEST
    add_definitions(-DGTEST)
    add_subdirectory(ExternalLibs/googletest)
    enable_testing()

    foreach( S ${TEST_NAMES} )
      message(STATUS "Adding test ${S}")
      add_executable(${S} test/${S}.cc)
      target_link_libraries(${S} gtest gtest_main ${LIB})
      target_include_directories(${S} PRIVATE ${gtest_SOURCE_DIR}/include ${gtest_SOURCE_DIR})
      gtest_discover_tests(${S})
    endforeach()

  elseif(${TEST} MATCHES "BOOST")
    #BOOST_TEST
    add_definitions(-DBOOST)
    add_definitions(-DBOOST_TEST_DYN_LINK)

    enable_testing()
    find_package(Boost COMPONENTS system filesystem unit_test_framework REQUIRED)
    set(Boost_USE_MULTITHREADED TRUE)

    foreach( S ${TEST_NAMES} )
      message(STATUS "Adding test ${S}")
      add_executable(${S} test/${S}.cc)
      target_link_libraries(${S} ${LIB} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
      add_test(NAME ${S} COMMAND ${S})
    endforeach()

#    add_executable(ConfigurationTest test/configurationTest.cc)
#    add_executable(ClothoidG1Test test/clothoidG1Test.cc)
#    add_executable(DubinsTest test/dubinsTest.cc)
#    target_link_libraries(ConfigurationTest ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
#    target_link_libraries(ClothoidG1Test ${LIB} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
#    target_link_libraries(DubinsTest ${LIB} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
#
#    add_test(NAME ConfigurationTest COMMAND ConfigurationTest)
#    add_test(NAME ClothoidG1Test COMMAND ClothoidG1Test)
#    add_test(NAME DubinsTest COMMAND DubinsTest)
  endif()

  add_custom_command(TARGET dubinsTest POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_SOURCE_DIR}/test/dubinsTest.txt ${CMAKE_CURRENT_BINARY_DIR}/test/dubinsTest.txt)

  add_custom_command(TARGET RSTest POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_SOURCE_DIR}/test/RSDataset.txt ${CMAKE_CURRENT_BINARY_DIR}/test/RSDataset.txt)

endif()
