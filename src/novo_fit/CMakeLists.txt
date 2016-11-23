include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# List the final name and cpp for each exe
HydraAddExecutable(NovoFit NovoFit.cpp)
if(ROOT_FOUND)
    target_link_libraries(NovoFit INTERFACE Root)
endif()

