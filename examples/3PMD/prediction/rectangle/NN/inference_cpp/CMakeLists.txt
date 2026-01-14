cmake_minimum_required(VERSION 3.11)
project(MyOnnxProject)

include(FetchContent)

set(CMAKE_CXX_STANDARD 14)

option(MY_INFERENCE_USE_CUDA "Download the CUDA-enabled ONNX Runtime and enable CUDA execution provider" ON)

# Define the URL for pre-built ONNX Runtime libraries
if (UNIX AND NOT APPLE)
  if (MY_INFERENCE_USE_CUDA)
    message(STATUS "Downloading onnxruntime-linux-x64-gpu")
    set(ONNXRUNTIME_URL "https://github.com/microsoft/onnxruntime/releases/download/v1.21.0/onnxruntime-linux-x64-gpu-1.21.0.tgz")
  else()
    message(STATUS "Downloading onnxruntime-linux-x64")
    set(ONNXRUNTIME_URL "https://github.com/microsoft/onnxruntime/releases/download/v1.21.0/onnxruntime-linux-x64-1.21.0.tgz")
  endif()
else()
  message(STATUS "Downloading onnxruntime-osx-arm64")
  set(ONNXRUNTIME_URL "https://github.com/microsoft/onnxruntime/releases/download/v1.21.0/onnxruntime-osx-arm64-1.21.0.tgz")
endif()

# Add -O3 optimization flag
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# Download and extract the pre-built ONNX Runtime
FetchContent_Declare(
  onnxruntime
  URL ${ONNXRUNTIME_URL}
)
FetchContent_MakeAvailable(onnxruntime)

set(CMAKE_OSX_ARCHITECTURES arm64)

# Print the source directory to see what was downloaded
message(STATUS "ONNX Runtime source dir: ${onnxruntime_SOURCE_DIR}")

# Add your executable
add_executable(my_inference_app src/inference.cpp src/model_classification.cpp src/model.cpp src/model_multi_task.cpp)

# Also add the root directory and some common patterns
target_include_directories(my_inference_app PRIVATE
  ${onnxruntime_SOURCE_DIR}/include
  ${CMAKE_SOURCE_DIR}/include
)

# Find the potential library directories
execute_process(
  COMMAND find ${onnxruntime_SOURCE_DIR} -name "libonnxruntime.so" -o -name "libonnxruntime.a" -o -name "libonnxruntime.dylib" 
  OUTPUT_VARIABLE ONNX_LIBS
)
message(STATUS "ONNX libraries found: ${ONNX_LIBS}")

# Try to locate the lib directory
execute_process(
  COMMAND find ${onnxruntime_SOURCE_DIR} -type d -name "lib"
  OUTPUT_VARIABLE ONNX_LIB_DIRS
)
message(STATUS "ONNX lib directories: ${ONNX_LIB_DIRS}")


target_link_directories(my_inference_app PRIVATE ${onnxruntime_SOURCE_DIR}/lib)

# Link against ONNX Runtime
target_link_libraries(my_inference_app PRIVATE onnxruntime)

if (MY_INFERENCE_USE_CUDA AND UNIX AND NOT APPLE)
  target_compile_definitions(my_inference_app PRIVATE USE_CUDA)
elseif(MY_INFERENCE_USE_CUDA)
  message(WARNING "CUDA execution provider requested, but this platform does not provide a GPU build. Falling back to CPU.")
endif()
