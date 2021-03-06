#Find tests generated in the binary directory

# To add a new file, just add it to this list.  Any GoogleTests inside will be automatically
# added to ctest.
set ( SNPLASH_TEST_CORE 
  ${CMAKE_CURRENT_SOURCE_DIR}/LR_Engine_Test.cpp 
  ${CMAKE_CURRENT_SOURCE_DIR}/StringUtils_Test.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/Statistics_Test.cpp
  PARENT_SCOPE)

SET(INTERFACE_LIBRARIES ${SNPLASH_TEST_CORE} PARENT_SCOPE )

# Add all the tests by parsing the source code
# This macro searches for GoogleTest macros and adds them as tests automatically
macro(ADD_GOOGLE_TESTS executable)
  # Add the generated tests
  set ( AllTestsHits ${GENERATED_TEST_LIST} )

  foreach ( source ${ARGN} )
    file(READ "${source}" contents)

    # Find all test lists
    string(REGEX MATCHALL "TEST_?F?\\(([A-Za-z_0-9 ,_]+)\\)" AllTests ${contents})

    # Convert the C++ code into a short test name
    foreach(hit ${AllTests})
      string(REGEX REPLACE ".*\\(( *[A-Za-z_0-9]+)[, ]*([A-Za-z_0-9]+) *\\).*" "\\1.\\2" test_name ${hit})
      set ( AllTestsHits ${AllTestsHits} ${test_name} )
    endforeach(hit)
  endforeach(source)

  MESSAGE("--tests: |${AllTestsHits}|")
  if(AllTestsHits STREQUAL "") 
    list ( SORT AllTestsHits)
  endif(AllTestsHits STREQUAL "")

  # Add all the remaining tests to CTest's test list
  foreach(hit ${AllTestsHits})
    add_test(${hit} ${executable} --gtest_filter=${hit})
  endforeach(hit)
endmacro(ADD_GOOGLE_TESTS)

