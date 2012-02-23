cmake_minimum_required(VERSION 2.8)

project(SNPLASH)

set(CMAKE_CXX_FLAGS "-O3 -fopenmp -g -Wall -Wunused")

# Edit for the Intel compiler
if (ICPC)
    set(CMAKE_CXX_FLAGS "-O3 -no-prec-div  -openmp -g -Wall -Wunused -wd981 -wd383 ")
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
SET( DOC_INPUT_DIR ${CYCLUS_SRC_DIR}/doc)

# If doxygen exists, use the doc/CMakeLists.txt file for further instructions. 
FIND_PACKAGE(Doxygen)
IF (DOXYGEN_FOUND)
	ADD_SUBDIRECTORY(docs)
	SET( DOC_OUTPUT_DIR ${CMAKE_BINARY_DIR}/docs )
ELSE (DOXYGEN_FOUND)
	MESSAGE(STATUS "WARNING: Doxygen not found - doc won't be created")
ENDIF (DOXYGEN_FOUND)

