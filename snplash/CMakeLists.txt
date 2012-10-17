# CMakeLists.txt
#
# Created by Richard T. Guy <gu@cs.toronto.edu>
# Modifications by David R. McWillams <dmcwilli@wakehealth.edu>
#

cmake_minimum_required(VERSION 2.4)

project(SNPLASH)

SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR})

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -g -Wall -Wunused")

# Edit appropriately for the Intel compiler.  The options below were
# chosen using the Intel 'Link Line Advisor v2.2' at
# http://software.intel.com/en-us/articles/intel-mkl-link-line-advisor/

if (ICPC)
    # set(CMAKE_CXX_FLAGS "-O3 -no-prec-div  -openmp -g -Wall -Wunused -wd981 -wd383 ")
    set(LINKER_OPTS "-L$(MKLROOT)/lib/intel64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm")
    set(COMP_OPTS "-I$(MKLROOT)/include")
    set(CMAKE_CXX_FLAGS "-O3 -no-prec-div -g Wall $COMP_OPTS $LINKER_OPTS")
endif (ICPC) 

#SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/..)

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

