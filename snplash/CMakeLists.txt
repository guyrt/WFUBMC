cmake_minimum_required(VERSION 2.0)

project(SNPLASH)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -g -c -Wall -Wunused $(VERSION_INFO) $(RUN_IN_PARALLEL)")

# The version number.
set (SNPLASH_VERSION_MAJOR 1)
set (SNPLASH_VERSION_MINOR 0)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/src/snplashConfig.h.in"
  "${PROJECT_BINARY_DIR}/snplashConfig.h"
  )

add_subdirectory (src) 

