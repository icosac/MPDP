# Compile the CUDA library and, optionally, the CUDA executable.

enable_language(CUDA)
find_package(CUDAToolkit REQUIRED)

set(LIB_CU "MPDPCU")

set(SRC_CU "${PROJECT_SOURCE_DIR}/srcCU")
file(GLOB MPDPSrcCU "${SRC_CU}/*.cu")
file(GLOB includeFilesCU "${SRC_CU}/include/*.cuh" "${PROJECT_SOURCE_DIR}/include/*.hh")


if(NOT DEFINED CUDA_ARCH OR CUDA_ARCH STREQUAL "")
  set(CUDA_ARCH "native")
  message(STATUS "CUDA_ARCH not set, defaulting to: ${CUDA_ARCH}")
endif()


if(CUDA_ARCH STREQUAL "native" AND CMAKE_VERSION VERSION_LESS 3.24)
  set(_detected_arch "")
  find_program(NVIDIA_SMI_EXECUTABLE nvidia-smi)
  if(NVIDIA_SMI_EXECUTABLE)
    execute_process(
            COMMAND ${NVIDIA_SMI_EXECUTABLE} --query-gpu=compute_cap --format=csv,noheader
            OUTPUT_VARIABLE _smi_out
            ERROR_QUIET
            RESULT_VARIABLE _smi_rc
            OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(_smi_rc EQUAL 0 AND NOT _smi_out STREQUAL "")
      # "8.9\n8.9" (one line per GPU) -> "89"
      string(REPLACE "\n" ";" _smi_list "${_smi_out}")
      foreach(_cap IN LISTS _smi_list)
        string(STRIP "${_cap}" _cap)
        string(REPLACE "." "" _cap "${_cap}")
        if(_cap MATCHES "^[0-9]+$")
          list(APPEND _detected_arch "${_cap}")
        endif()
      endforeach()
      list(REMOVE_DUPLICATES _detected_arch)
    endif()
  endif()

  if(_detected_arch)
    set(CUDA_ARCH "${_detected_arch}")
    message(STATUS "CMake ${CMAKE_VERSION} cannot pass 'native' through; detected: ${CUDA_ARCH}")
  else()
    set(CUDA_ARCH "70;75;80;86;89")
    message(WARNING "Could not detect the local CUDA architecture, building for: ${CUDA_ARCH}")
  endif()
endif()

if(NOT DEFINED CUDA_SEP_COMP OR CUDA_SEP_COMP STREQUAL "")
  set(CUDA_SEP_COMP OFF)
endif()

add_library(${LIB_CU} STATIC ${MPDPSrcCU})

target_include_directories(${LIB_CU} PUBLIC
        $<BUILD_INTERFACE:${SRC_CU}/include>
        $<BUILD_INTERFACE:${PROJECT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_compile_definitions(${LIB_CU} PUBLIC CUDA_ON)

target_compile_options(${LIB_CU} PRIVATE
        $<$<COMPILE_LANGUAGE:CUDA>:--expt-relaxed-constexpr>
)

set_target_properties(${LIB_CU} PROPERTIES
        CUDA_SEPARABLE_COMPILATION ${CUDA_SEP_COMP}
        CUDA_ARCHITECTURES "${CUDA_ARCH}"
        CUDA_STANDARD 17
        CUDA_STANDARD_REQUIRED ON
        POSITION_INDEPENDENT_CODE ON
)

if(DEFINED CUDA_FLAGS AND NOT CUDA_FLAGS STREQUAL "")
  separate_arguments(CUDA_FLAGS_LIST NATIVE_COMMAND "${CUDA_FLAGS}")
  target_compile_options(${LIB_CU} PRIVATE
          $<$<COMPILE_LANGUAGE:CUDA>:${CUDA_FLAGS_LIST}>)
endif()

target_link_libraries(${LIB_CU} PUBLIC CUDA::cudart)

message(STATUS "Compiling CUDA library ${LIB_CU} for architectures: ${CUDA_ARCH}")
message(STATUS "CUDA sources: ${MPDPSrcCU}")

install(FILES ${includeFilesCU} DESTINATION ${CMAKE_BINARY_DIR}/inc)

if(COMPILE_CUDA_EXEC)
  set(APP_EXEC_CU "MPDPCU_exec")
  message(STATUS "Compiling CUDA executable ${APP_EXEC_CU}")

  add_executable(${APP_EXEC_CU} exec/main.cu)
  target_link_libraries(${APP_EXEC_CU} PUBLIC ${LIB_CU})
  set_target_properties(${APP_EXEC_CU} PROPERTIES
          CUDA_SEPARABLE_COMPILATION ${CUDA_SEP_COMP}
          CUDA_ARCHITECTURES "${CUDA_ARCH}"
          CUDA_STANDARD 17
          CUDA_STANDARD_REQUIRED ON
  )
endif()
