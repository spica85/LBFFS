# Setting for the version of CMake
cmake_minimum_required(VERSION 3.16.1)
# Setting for clang++
set(CMAKE_CXX_COMPILER /bin/clang++)
# Setting for the project name and coding language
project(D3Q19_opencl CXX)
# Making an execution file named ns.out from ns.cpp
add_executable(ns.opencl ns.cpp)
# Setting for the option of g++
set(CMAKE_CXX_FLAGS "-O3 -fopenmp -std=c++11")

# For OpenCL
set(OpenCL_LIBRARY /opt/amdgpu-pro/lib/x86_64-linux-gnu/libOpenCL.so)
set(OpenCL_INCLUDE_DIR /usr/include/CL)
find_package(OpenCL REQUIRED)
include_directories(${OpenCL_INCLUDE_DIRS})
target_link_libraries(ns.opencl ${OpenCL_LIBRARIES})
configure_file(kernel.cl ../build/kernel.cl COPYONLY)

