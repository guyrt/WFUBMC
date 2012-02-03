
add_subdirectory (reader)
add_subdirectory (param) 
add_subdirectory (engine)
add_subdirectory (logger)

add_executable (snplash snpadt.cpp) 

# Link the executable to the libraries. 
link_directories (${HELLO_BINARY_DIR}/reader) 
target_link_libraries (snplash reader) 
target_link_libraries (snplash params) 
target_link_libraries (snplash linalg) 
target_link_libraries (snplash machlearn) 
target_link_libraries (snplash coreengine) 
target_link_libraries (snplash intertwolog) 
target_link_libraries (snplash engineutils) 
target_link_libraries (snplash logger) 
target_link_libraries (snplash ld) 
target_link_libraries (snplash dandelion)
target_link_libraries (snplash snpgwa)  
target_link_libraries (snplash engineout)
target_link_libraries (snplash qsnpgwa)  
target_link_libraries (snplash expectmax)  



