cmake_minimum_required(VERSION 2.4)

project(SNPLASH)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -g -Wall -Wunused")

# Edit for the Intel compiler
if (ICPC)
    set(CMAKE_CXX_FLAGS "-O3 -no-prec-div  -openmp -g -Wall -Wunused -wd981 -wd383 ")
endif (ICPC) 

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/..)

# The version number.
set (SNPLASH_VERSION_MAJOR 1)
set (SNPLASH_VERSION_MINOR 0)

# configure a header file to pass some of the CMake settings
# to the source code
configure_file (
  "${PROJECT_SOURCE_DIR}/src/snplashConfig.h"
  "${PROJECT_BINARY_DIR}/src/snplashConfig.h"
  )

add_subdirectory (src) 

