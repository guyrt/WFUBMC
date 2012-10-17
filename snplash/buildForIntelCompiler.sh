#!bin/sh

# buildForIntelCompiler
#
# Create the environment and directory structure to compile SNPLASH
# with the Intel compiler suite.

mkdir build
cd build

module load compilers/intel-2012-lp64

# ICPC=`which icpc`

# Intel compiler command
MPICC=`which MPICC`

cmake -D CMAKE_CXX_COMPILER=$MPICC -D ICPC=TRUE ..

echo "Now cd into build and type make. The executable will be in this directory."


