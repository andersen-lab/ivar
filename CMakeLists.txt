cmake_minimum_required(VERSION 3.10)

# Project name and version
project(ivar VERSION 2.0)

# Use ccache to speed up rebuilds when available (no-op if not installed)
find_program(CCACHE_PROGRAM ccache)
if(CCACHE_PROGRAM)
    set(CMAKE_CXX_COMPILER_LAUNCHER "${CCACHE_PROGRAM}")
endif()

# Enable testing
enable_testing()

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

add_compile_options(-Wall)
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()
set(CMAKE_CXX_FLAGS_RELEASE "-O2 -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g3")

# Add subdirectories
add_subdirectory(src)
add_subdirectory(tests)


# Installation rules
install(TARGETS ivar DESTINATION bin)