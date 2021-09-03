# The first tutorial about CMakeLists.txt:
# <https://progcz.com/posts/cmake-tutorial-write-cmakelists-txt-from-the-very-beginning/>

cmake_minimum_required(VERSION 3.15)

# Set the name of project and version info.
project(ngslib VERSION 1.0.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Configure a header file to pass some of the CMake settings
# to the source codes.

# Add the binary tree to the search path for include files
