include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# List the final name and cpp for each exe
HydraAddExecutable(SimpleFit SimpleFit.cpp)
if(ROOT_FOUND)
    # All the targets, cuda, omp, and tbb, use this interface
    target_link_libraries(SimpleFit INTERFACE Root)
endif()

