# List the final name and cpp for each exe

HydraAddExecutable(PhSp PhSp.cpp)
if(ROOT_FOUND)
    # All the targets, cuda, omp, and tbb, use this interface
    target_link_libraries(PhSp INTERFACE Root)
endif()
