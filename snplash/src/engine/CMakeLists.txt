
# Add each subdir.
add_subdirectory(dandelion)
add_subdirectory(em)
add_subdirectory(intertwolog)
add_subdirectory(ld)
add_subdirectory(linalg)
add_subdirectory(machineLearning)
add_subdirectory(output)
add_subdirectory(qsnpgwa)
add_subdirectory(snpgwa)
add_subdirectory(utils)

# This library contains core items.
add_library ( coreengine data_plugin.cpp randwh.cpp snp_data.cpp )

