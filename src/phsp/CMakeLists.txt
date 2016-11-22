# List the final name and cpp for each exe

HydraAddExecutable(PhSp PhSp.cpp)
if(ROOT_FOUND)
    # All the targets, cuda, omp, and tbb, use this interface
    target_link_libraries(PhSp INTERFACE ${ROOT_LIBRARIES} Threads::Threads)
    target_include_directories(PhSp INTERFACE ${ROOT_INCLUDE_DIR})
    target_compile_definitions(PhSp INTERFACE FOUND_ROOT)
endif()
