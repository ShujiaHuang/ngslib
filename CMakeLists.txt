# The first tutorial about CMakeLists.txt:
# <https://progcz.com/posts/cmake-tutorial-write-cmakelists-txt-from-the-very-beginning/>

cmake_minimum_required(VERSION 2.8.12)
set(CMAKE_C_COMPILER "gcc")
set(CMAKE_CXX_COMPILER "g++")
set(CMAKE_CXX_FLAGS "-std=c++0x -O3 -pthread")

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Set the name of project and version info.
project(ngslib VERSION 1.0.0)

# Configure a header file to pass some of the CMake settings
# to the source codes.

# Add the binary tree to the search path for include files
