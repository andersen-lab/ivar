cmake_minimum_required(VERSION 3.10)

# Project name and version
project(ivar VERSION 2.0)

# Enable testing
enable_testing()

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -Wall -O2")

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)


# Installation rules
install(TARGETS ivar DESTINATION bin)