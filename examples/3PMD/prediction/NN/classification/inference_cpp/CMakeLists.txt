cmake_minimum_required(VERSION 3.11)
project(MyOnnxProject)

include(FetchContent)

# Define the URL for pre-built ONNX Runtime binaries (for Linux x64)
set(ONNXRUNTIME_URL "https://github.com/microsoft/onnxruntime/releases/download/v1.21.0/onnxruntime-linux-x64-1.21.0.tgz")

# Download and extract the pre-built ONNX Runtime
FetchContent_Declare(
  onnxruntime
  URL ${ONNXRUNTIME_URL}
  URL_HASH ""  # You can add a hash for security if needed
)
FetchContent_MakeAvailable(onnxruntime)

# Print the source directory to see what was downloaded
message(STATUS "ONNX Runtime source dir: ${onnxruntime_SOURCE_DIR}")

# Execute a command to list the directory structure
execute_process(
  COMMAND find ${onnxruntime_SOURCE_DIR} -type f -name "*.h" | grep -i onnx
  OUTPUT_VARIABLE ONNX_HEADERS
)
message(STATUS "ONNX headers found: ${ONNX_HEADERS}")

# Find the include directory structure
execute_process(
  COMMAND find ${onnxruntime_SOURCE_DIR} -type d -name "include"
  OUTPUT_VARIABLE ONNX_INCLUDE_DIRS
)
message(STATUS "ONNX include directories: ${ONNX_INCLUDE_DIRS}")

# Add your executable
add_executable(my_inference_app src/inference.cpp)

# Add all potential include directories
foreach(dir ${ONNX_INCLUDE_DIRS})
  string(STRIP ${dir} dir_stripped)
  target_include_directories(my_inference_app PRIVATE ${dir_stripped})
endforeach()

# Also add the root directory and some common patterns
target_include_directories(my_inference_app PRIVATE
  ${onnxruntime_SOURCE_DIR}
  ${onnxruntime_SOURCE_DIR}/include
  ${onnxruntime_SOURCE_DIR}/onnxruntime-linux-x64-1.16.3/include
  ${onnxruntime_SOURCE_DIR}/onnxruntime-linux-x64/include
)

# Find the potential library directories
execute_process(
  COMMAND find ${onnxruntime_SOURCE_DIR} -name "*.so" -o -name "*.a" | grep -i onnx
  OUTPUT_VARIABLE ONNX_LIBS
)
message(STATUS "ONNX libraries found: ${ONNX_LIBS}")

# Try to locate the lib directory
execute_process(
  COMMAND find ${onnxruntime_SOURCE_DIR} -type d -name "lib"
  OUTPUT_VARIABLE ONNX_LIB_DIRS
)
message(STATUS "ONNX lib directories: ${ONNX_LIB_DIRS}")

# Add all potential library directories
foreach(dir ${ONNX_LIB_DIRS})
  string(STRIP ${dir} dir_stripped)
  target_link_directories(my_inference_app PRIVATE ${dir_stripped})
endforeach()

# Link against ONNX Runtime
target_link_libraries(my_inference_app PRIVATE onnxruntime)