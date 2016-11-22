include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# List the final name and cpp for each exe
HydraAddExecutable(NewFit NewFit.cpp)
if(ROOT_FOUND)
    # All the targets, cuda, omp, and tbb, use this interface
    target_link_libraries(NewFit INTERFACE ${ROOT_LIBRARIES} Threads::Threads)
    target_include_directories(NewFit INTERFACE ${ROOT_INCLUDE_DIR})
    target_compile_definitions(NewFit INTERFACE FOUND_ROOT)
endif()

