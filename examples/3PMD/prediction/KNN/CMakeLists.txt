cmake_minimum_required(VERSION 3.14)
project(KNN_Annoy)

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# FetchContent for Annoy
include(FetchContent)
FetchContent_Declare(
    annoy
    GIT_REPOSITORY https://github.com/spotify/annoy.git
    GIT_TAG main  # You can also specify a commit hash or release tag
)
FetchContent_MakeAvailable(annoy)

# Add executable
add_executable(knn KNN.cpp)

# Include Annoy headers
target_include_directories(knn PRIVATE ${annoy_SOURCE_DIR}/src)
