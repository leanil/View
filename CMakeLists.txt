cmake_minimum_required(VERSION 3.2.2)

project(View)

list(APPEND CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/Modules)
include(FindComputeCpp)

# Set the host compiler C++ standard to C++ 14
# This should be set as a global property when we update the CMake version
set(CMAKE_CXX_STANDARD 14)

# Set include ComputeCpp directories
include_directories(SYSTEM ${COMPUTECPP_INCLUDE_DIRECTORY})

add_subdirectory(View)

