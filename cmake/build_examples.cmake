# Compile examples

set(DEMOS_LIST DPS 3PMD MPMD P2PDubinsDataset P2PRSDataset)

message(STATUS "Compiling the following demos: ${DEMOS_LIST} # remember to add Demo at the end of the name")

foreach(DEMO ${DEMOS_LIST})
        # file(GLOB_RECURSE ${DEMO}_SOURCES examples/${DEMO}/*.cc)
        file(GLOB ${DEMO}_SOURCES examples/${DEMO}/*.cc)
        add_executable(${DEMO}Demo ${${DEMO}_SOURCES})
        target_link_libraries(${DEMO}Demo ${LIB_CXX})
        target_include_directories(${DEMO}Demo PUBLIC
                ${SRC_CXX}/include
                ${CMAKE_SOURCE_DIR}/include
                ${CMAKE_SOURCE_DIR}/examples/${DEMO})
endforeach(DEMO)

##########################################################################################
# MPMDBenchmark - CPU vs GPU comparison
##########################################################################################
# This demo does not fit the loop above, which builds one executable from every .cc in a
# directory. It needs several: the problem generator, the CPU driver and the GPU driver
# stay separate processes exchanging a problem file, which keeps CUDA context creation out
# of the CPU timings and guarantees both solvers see byte-identical inputs.

set(BENCH_DIR "${CMAKE_SOURCE_DIR}/examples/MPMDBenchmark")

# Writes the shared problem set. Standard library only, no MPDP dependency.
add_executable(MPMDBenchmarkGen "${BENCH_DIR}/gen_problems.cc")
target_include_directories(MPMDBenchmarkGen PRIVATE "${BENCH_DIR}")

# The CPU reference, linked against the C++ library only.
add_executable(MPMDBenchmarkCPU "${BENCH_DIR}/bench_cpu.cc")
target_include_directories(MPMDBenchmarkCPU PRIVATE "${BENCH_DIR}")
target_link_libraries(MPMDBenchmarkCPU PRIVATE ${LIB_CXX})

set(MPMD_BENCHMARK_TARGETS MPMDBenchmarkGen MPMDBenchmarkCPU)

# The GPU solver, linked against the CUDA library only.
if(CUDA_ON)
        add_executable(MPMDBenchmarkGPU "${BENCH_DIR}/bench_gpu.cu")
        target_include_directories(MPMDBenchmarkGPU PRIVATE "${BENCH_DIR}")
        target_link_libraries(MPMDBenchmarkGPU PRIVATE ${LIB_CU})
        set_target_properties(MPMDBenchmarkGPU PROPERTIES
                CUDA_SEPARABLE_COMPILATION ${CUDA_SEP_COMP}
                CUDA_ARCHITECTURES "${CUDA_ARCH}"
                CUDA_STANDARD 17
                CUDA_STANDARD_REQUIRED ON
        )
        list(APPEND MPMD_BENCHMARK_TARGETS MPMDBenchmarkGPU)

        add_executable(MPMDBenchmarkRSCheck "${BENCH_DIR}/rs_check.cu")
        target_include_directories(MPMDBenchmarkRSCheck PRIVATE
                "${BENCH_DIR}"
                "${SRC_CXX}/include"
                "${SRC_CU}/include"
                "${CMAKE_SOURCE_DIR}/include")
        target_link_libraries(MPMDBenchmarkRSCheck PRIVATE ${LIB_CXX})
        set_target_properties(MPMDBenchmarkRSCheck PROPERTIES
                CUDA_ARCHITECTURES "${CUDA_ARCH}"
                CUDA_STANDARD 17
                CUDA_STANDARD_REQUIRED ON
        )
        list(APPEND MPMD_BENCHMARK_TARGETS MPMDBenchmarkRSCheck)
else()
        message(STATUS "CUDA_ON is off: MPMDBenchmark builds the CPU baseline only")
endif()

# Convenience target so the whole demo builds with one name:
#   cmake --build <dir> --target MPMDBenchmark
add_custom_target(MPMDBenchmark DEPENDS ${MPMD_BENCHMARK_TARGETS})

message(STATUS "MPMDBenchmark executables: ${MPMD_BENCHMARK_TARGETS}")
message(STATUS "  run it with examples/MPMDBenchmark/run_benchmark.sh")
