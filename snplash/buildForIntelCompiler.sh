mkdir build
cd build

ICPC=`which icpc`

cmake -D CMAKE_CXX_COMPILER=$ICPC -D ICPC=TRUE ..

echo "Now cd into build and type make. The executable will be in this directory."


