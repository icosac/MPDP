# Compile examples

# 3PMD
file(GLOB_RECURSE 3PMD_SOURCES examples/3PMD/*.cc)
add_executable(3PMDDemo ${3PMD_SOURCES})
target_link_libraries(3PMDDemo ${LIB_CXX})
target_include_directories(3PMDDemo PUBLIC
        ${SRC_CXX}/include
        ${CMAKE_SOURCE_DIR}/examples/3PMD
        ${CMAKE_SOURCE_DIR}/include)


