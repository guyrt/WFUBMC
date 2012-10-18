# CMakeLists.txt
#
# Authored by: Richard T. Guy <gu@cs.toronto.edu>
# Modified by: David R. McWillams <dmcwilli@wakehealth.edu>

# Code for modules system created by David W. Chen, Wake Forest
# University

cmake_minimum_required(VERSION 2.4)

project(SNPLASH)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -g -Wall -Wunused")

# Edit appropriately for the Intel compiler.  The options below were
# chosen using the Intel 'Link Line Advisor v2.2' at
# http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/

### As of 2012-10-17
# Intel Product: Intel(R) Composer XE 2011
# OS: Linux
# Compiler: Intel(R) C/C++
# Architecture: Intel(R) 64
# Dynamic or static linking: Dynamic
# Interface layer: LP64 (32-bit integer)
# Sequential or multi-threaded layer: Multi-threaded
# OpenMP library: Intel(R) (libiomp5)
# Cluster library: [NONE]
# MPI library: [NONE]
# Fortran 95 interfaces: [NONE]
# Link with Intel MKL libraries explicitly: [UNCHECKED]

if (ICPC)
    # set(CMAKE_CXX_FLAGS "-O3 -no-prec-div  -openmp -g -Wall -Wunused -wd981 -wd383 ")
    set(COMP_OPTS "-I${MKLROOT}/include")
    set(CMAKE_CXX_FLAGS "-xSSE4.2 -O3 -no-prec-div -openmp -mkl=parallel -g -Wall -Wunused -wd981 -wd383 ${COMP_OPTS}")
endif (ICPC) 

# SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/..)

# The version number.
set (SNPLASH_VERSION_MAJOR 1)
set (SNPLASH_VERSION_MINOR 0)

# Configure a header file to pass some of the CMake settings
# to the source code.
configure_file (
  "${PROJECT_SOURCE_DIR}/src/snplashConfig.h"
  "${PROJECT_BINARY_DIR}/src/snplashConfig.h"
  )

add_subdirectory (src) 

# This is the directory that holds the doxygen doxyfile template (doxy.conf.in)
SET( DOC_INPUT_DIR ${PROJECT_SOURCE_DIR}/docs)

# If doxygen exists, use the doc/CMakeLists.txt file for further instructions. 
FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND)
	ADD_SUBDIRECTORY(docs)
	SET( DOC_OUTPUT_DIR ${PROJECT_SOURCE_DIR}/docs )
	MESSAGE(STATUS "Doxy out: ${DOC_OUTPUT_DIR}")
	MESSAGE(STATUS "Doxy in: ${DOC_INPUT_DIR}")
ELSE (DOXYGEN_FOUND)
	MESSAGE(STATUS "WARNING: Doxygen not found - doc won't be created")
ENDIF (DOXYGEN_FOUND)

