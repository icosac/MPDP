# Compile examples

set(DEMOS_LIST 3PMD MPMD P2PDubinsDataset P2PRSDataset)

message(STATUS "Compiling the following demos: ${DEMOS_LIST}")

foreach(DEMO ${DEMOS_LIST})
        file(GLOB_RECURSE ${DEMO}_SOURCES examples/${DEMO}/*.cc)
        add_executable(${DEMO}Demo ${${DEMO}_SOURCES})
        target_link_libraries(${DEMO}Demo ${LIB_CXX})
        target_include_directories(${DEMO}Demo PUBLIC
                ${SRC_CXX}/include
                ${CMAKE_SOURCE_DIR}/include
                ${CMAKE_SOURCE_DIR}/examples/${DEMO})
endforeach(DEMO)


