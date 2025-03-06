cmake_minimum_required(VERSION 3.14)
project(BaggedTrees)

# Enable C++17 or higher
set(CMAKE_CXX_STANDARD 17)
# set(CMAKE_CXX_STANDARD_REQUIRED True)

# Download and install ONNX Runtime
include(FetchContent)
FetchContent_Declare(
    onnxruntime
    GIT_REPOSITORY https://github.com/microsoft/onnxruntime.git
    GIT_TAG v1.20.1  # Set to the latest stable version
)

# Extract and set ONNX runtime paths
FetchContent_MakeAvailable(onnxruntime)
set(ONNXRUNTIME_INCLUDE_DIR ${onnxruntime_SOURCE_DIR}/include)
set(ONNXRUNTIME_LIB_DIR ${onnxruntime_SOURCE_DIR}/lib)

# Add executable
add_executable(bagged_trees bagged_trees.cc)

# Include ONNX Runtime headers
target_include_directories(bagged_trees PRIVATE ${ONNXRUNTIME_INCLUDE_DIR})

# Link ONNX Runtime library
target_link_directories(bagged_trees PRIVATE ${ONNXRUNTIME_LIB_DIR})
target_link_libraries(bagged_trees PRIVATE onnxruntime)
