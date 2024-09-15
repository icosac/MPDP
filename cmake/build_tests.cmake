# Compile tests

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
      add_executable(${S} test/${S}.cc
              examples/3PMD/3PMD.cc
              examples/3PMD/3PMD.hh)
      target_link_libraries(${S} ${LIB} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
      add_test(NAME ${S} COMMAND ${S})
  endforeach()

endif()

add_custom_command(TARGET dubinsTest POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_SOURCE_DIR}/test/dubinsTest.txt ${CMAKE_CURRENT_BINARY_DIR}/test/dubinsTest.txt)

add_custom_command(TARGET RSTest POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy
      ${CMAKE_SOURCE_DIR}/test/RSDataset.txt ${CMAKE_CURRENT_BINARY_DIR}/test/RSDataset.txt)
