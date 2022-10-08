# Setting for the version of CMake
cmake_minimum_required(VERSION 3.16.1)
# Setting for clang++
set(CMAKE_CXX_COMPILER /usr/bin/clang++)
# Setting for the project name and coding language
project(D3Q19_opencl CXX)
# Making an execution file named ns.out from ns.cpp
add_executable(ns.opencl ns.cpp)
# Setting for the option of g++
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -std=c++11")

# For Boost
find_package(Boost REQUIRED)

# For OpenCL
#set(OpenCL_LIBRARY /opt/amdgpu-pro/lib/x86_64-linux-gnu/libOpenCL.so) # If necessary
find_package(OpenCL REQUIRED)

include_directories(${OpenCL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
target_link_libraries(testSTL ${OpenCL_LIBRARIES})
