
FILE(GLOB SNPLASH_CORE ${CMAKE_CURRENT_SOURCE_DIR}/engine/*.cpp) 
FILE(GLOB SNPLASH_CORE_T ${CMAKE_CURRENT_SOURCE_DIR}/param/*.cpp)
SET(SNPLASH_CORE ${SNPLASH_CORE} ${SNPLASH_CORE_T})
FILE(GLOB SNPLASH_CORE_T ${CMAKE_CURRENT_SOURCE_DIR}/reader/*.cpp)
SET(SNPLASH_CORE ${SNPLASH_CORE} ${SNPLASH_CORE_T})
FILE(GLOB SNPLASH_CORE_T ${CMAKE_CURRENT_SOURCE_DIR}/logger/*.cpp)
SET(SNPLASH_CORE ${SNPLASH_CORE} ${SNPLASH_CORE_T})
FILE(GLOB SNPLASH_CORE_T ${CMAKE_CURRENT_SOURCE_DIR}/engine/*/*.cpp)
SET(SNPLASH_CORE ${SNPLASH_CORE} ${SNPLASH_CORE_T})

#MESSAGE(${SNPLASH_CORE})


add_subdirectory (reader)
add_subdirectory (param) 
add_subdirectory (engine)
add_subdirectory (logger)

add_executable (snplash snpadt.cpp) 

target_link_libraries (snplash machlearn) 
target_link_libraries (snplash intertwolog) 
target_link_libraries (snplash logger) 
target_link_libraries (snplash ld) 
target_link_libraries (snplash dandelion)
target_link_libraries (snplash snpgwa)  
target_link_libraries (snplash qsnpgwa)  
target_link_libraries (snplash expectmax)  

target_link_libraries (snplash engineout)

target_link_libraries (snplash engineutils) 
target_link_libraries (snplash linalg) 
target_link_libraries (snplash coreengine) 

target_link_libraries (snplash reader) 
target_link_libraries (snplash param) 


#-----------
# TESTS
#-----------

# Be sure to clear these each time
set ( GENERATED_TEST_LIST "" CACHE INTERNAL "" )
set ( GENERATED_FILTER_LIST "" CACHE INTERNAL "" )
set ( GENERATED_TEST_SOURCE_LIST "" CACHE INTERNAL "" )

# The gtest include file directories.
SET(GTEST_INCLUDE_DIRS "${CMAKE_CURRENT_SOURCE_DIR}/Testing/GoogleTest")

# The gtest library directories.
SET(GTEST_LIBRARY_DIRS "${CMAKE_CURRENT_BINARY_DIR}/Testing")

# Alerts CMake to any tests we need to worry about
INCLUDE_DIRECTORIES(${GTEST_INCLUDE_DIRS})

# Create an executable for all the gtests
# Options for testing
option ( USE_TESTING "Build testing" ON )
if ( USE_TESTING )
  enable_testing()
  include(CTest)
endif()

ADD_SUBDIRECTORY(Testing)
SET(TestSource 
  ${TestSource} 
  # Need to change this next line.
  ${SNPLASH_CORE}
  ${SNPLASH_TEST_CORE}
  )

# Build Google Test
find_package(Threads)
if (CMAKE_USE_PTHREADS_INIT)  # The pthreads library is available.
  set(cxx_base_flags "${cxx_base_flags} -DGTEST_HAS_PTHREAD=1")
endif()

# Build Google Testing
set ( GTest
  Testing/GoogleTest/gtest/gtest-all.cc
)
include_directories ( Testing ) # necessary?

add_library(gtest ${GTest})
if (CMAKE_USE_PTHREADS_INIT)
  target_link_libraries(gtest ${CMAKE_THREAD_LIBS_INIT})
endif()

INCLUDE_DIRECTORIES(Testing/Googletest) # necessary?
add_executable(snplashUnitTestDriver
  Testing/snplashUnitTestDriver.cpp ${TestSource}
  ${GENERATED_TEST_SOURCE}
)

target_link_libraries ( snplashUnitTestDriver dl ${SNPLASH_LIBRARIES} 
  ${INTERFACE_LIBRARIES} dl ${LIBS} gtest)

# Here we add all the tests
add_google_tests ( ${EXECUTABLE_OUTPUT_PATH}/snplashUnitTestDriver ${TestSource} 
  ${INTERFACE_LIBRARIES})





