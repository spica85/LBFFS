# Setting for the version of CMake
cmake_minimum_required(VERSION 3.16.1)
# Setting for the project name and coding language
project(D3Q19_cmake CXX)
# Making an execution file named ns.out from ns.cpp
add_executable(ns.out ns.cpp LKS.cpp input.cpp)
# Setting for the option of g++
set(CMAKE_CXX_FLAGS "-O3 -march=native -fopenmp -std=c++11")

# For ArrayFire
find_package(ArrayFire)
target_link_libraries(ns.out ArrayFire::afopencl) # afopencl,afcpu
