# Compile CUDA library and executable

enable_language(CUDA)

set(LIB_CU "MPDPCU")

#Define files to be compiled
set(SRC "srcCU")
set(SRCF "${SRC}/*.cu")
file(GLOB MPDPSrc ${SRCF})
file(GLOB includeFiles "${SRC}/include/*.cuh" "./include/*.hh")

#What should be compiled
include_directories(${SRC}/include include)
add_library(${LIB_CU} ${MPDPSrc})

#What should be installed
install(FILES ${includeFiles} DESTINATION ${CMAKE_BINARY_DIR}/inc)

#Set compiler options
set(CMAKE_CUDA_FLAGS ${CMAKE_CUDA_FLAGS} ${CUDA_FLAGS})
set_target_properties(${LIB_CU} PROPERTIES CUDA_SEPARABLE_COMPILATION ${CUDA_SEP_COMP})
set_target_properties(${LIB_CU} PROPERTIES CUDA_ARCHITECTURES ${CUDA_ARCH})

message(STATUS "Compiling CUDA with: " ${CUDA_FLAGS})
message(STATUS "Compiling CUDA library: " ${MPDPSrc})

if (COMPILE_CUDA_EXEC)
  # CUDA executable
  set(APP_EXEC_CU "MPDPCU_exec")
  file(GLOB MPDPExec "exec/main.cu")

  message(STATUS "Compiling CUDA executable ${APP_EXEC_CU}")

  add_executable(${APP_EXEC_CU} ${MPDPExec})
  target_link_libraries(${APP_EXEC_CU} PUBLIC ${LIB})
  set_target_properties(${APP_EXEC_CU} PROPERTIES CUDA_SEPARABLE_COMPILATION ${CUDA_SEP_COMP})
  set_target_properties(${APP_EXEC_CU} PROPERTIES CUDA_ARCHITECTURES ${CUDA_ARCH})
endif()