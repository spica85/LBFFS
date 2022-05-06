# Setting for the version of CMake
cmake_minimum_required(VERSION 3.16.1)
# Setting for clang++
set(CMAKE_CXX_COMPILER /bin/clang++)
# Setting for the project name and coding language
project(D3Q19_cmake CXX)
# Making an execution file named ns.out from ns.cpp
add_executable(ns.out ns.cpp)
# Setting for the option of g++
set(CMAKE_CXX_FLAGS "-O3 -Xpreprocessor  -fopenmp -lomp -std=c++11 -L/opt/homebrew/lib -I/opt/homebrew/include")
