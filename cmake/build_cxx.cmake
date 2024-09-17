# Compile C++ library and executable

set(LIB_CXX "MPDPCC")
set(SRC_CXX "${CMAKE_SOURCE_DIR}/srcCC")
file(GLOB MPDPSrc "${SRC_CXX}/*.cc")

message(STATUS "Compiling C++ library with source files: ${MPDPSrc}")
message(STATUS "Compiling with options " ${CMAKE_CXX_FLAGS})

add_library(${LIB_CXX} STATIC ${MPDPSrc})

target_include_directories(${LIB_CXX} PUBLIC
        $<BUILD_INTERFACE:${SRC_CXX}/include>
        $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

# Optional external library if MPDP_DRAW is enabled
if (${MPDP_DRAW})
    include_directories(ExternalLibs/clothoids)
    add_compile_definitions(MPDP_DRAW)
endif()

# Compile the executable if the flag is set
if(COMPILE_CXX_EXEC)
    set(APP_EXEC_CXX "MPDPCC_exec")
    file(GLOB MPDPExec "exec/main.cc")  # Use ${SRC_CXX} for the correct path

    message(STATUS "Compiling C++ executable ${APP_EXEC_CXX} with source: ${MPDPExec}")

    # Define the executable
    add_executable(${APP_EXEC_CXX} exec/main.cc)

    # Set include directories for the executable (build and install)
    target_include_directories(${APP_EXEC_CXX} PUBLIC
            $<BUILD_INTERFACE:${SRC_CXX}/include>
            $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
            $<INSTALL_INTERFACE:include>
    )

    # Link the executable with the library
    target_link_libraries(${APP_EXEC_CXX} PUBLIC ${LIB_CXX})

    # Link OpenMP if found
    if (OpenMP_CXX_FOUND)
        target_link_libraries(${LIB_CXX} PUBLIC OpenMP::OpenMP_CXX)
    endif()
endif()

# Define installation directories using standard CMake variables
set(CMAKE_INSTALL_LIBDIR "${CMAKE_BINARY_DIR}/lib")
set(CMAKE_INSTALL_BINDIR "${CMAKE_BINARY_DIR}/bin")
set(CMAKE_INSTALL_INCLUDEDIR "${CMAKE_BINARY_DIR}/include")
set(CMAKE_INSTALL_CMAKEDIR "${CMAKE_INSTALL_LIBDIR}/cmake")

# Log the install destinations
message(STATUS "Installing C++ library to ${CMAKE_INSTALL_LIBDIR}")
message(STATUS "Installing C++ headers to ${CMAKE_INSTALL_INCLUDEDIR}")
message(STATUS "Installing C++ cmake files to ${CMAKE_INSTALL_CMAKEDIR}")

# Install the library
install(TARGETS ${LIB_CXX}
        EXPORT ${LIB_CXX}
        ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
        LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

# Install the include directories
install(DIRECTORY ${SRC_CXX}/include
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)
install(DIRECTORY ${CMAKE_SOURCE_DIR}/include
        DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

# Install the CMake export file for the library
install(EXPORT ${LIB_CXX}
        FILE ${LIB_CXX}.cmake
        NAMESPACE MPDPCC::
        DESTINATION ${CMAKE_INSTALL_CMAKEDIR}
)

# Create a Config file for find_package
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
        "${CMAKE_CURRENT_BINARY_DIR}/MPDPCCConfigVersion.cmake"
        VERSION 1.0
        COMPATIBILITY AnyNewerVersion
)

# Install the Config file
install(FILES
        "${CMAKE_CURRENT_BINARY_DIR}/MPDPCCConfig.cmake"       # The main config file
        "${CMAKE_CURRENT_BINARY_DIR}/MPDPCCConfigVersion.cmake" # Version file
        DESTINATION lib/cmake/MPDPCC                           # Destination for config files
)

# Create the MPDPCCConfig.cmake file
configure_file(MPDPCCConfig.cmake.in
        "${CMAKE_CURRENT_BINARY_DIR}/MPDPCCConfig.cmake" @ONLY)