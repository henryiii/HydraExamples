include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# List the final name and cpp for each exe
# Doesn't build on CUDA, doesn't work anywhere else
if(NOT CUDA_FOUND)
    HydraAddExecutable(NovoFit NovoFit.cpp)
    if(ROOT_FOUND)
        target_link_libraries(NovoFit INTERFACE Root)
    endif()
endif()

# List the final name and cpp for each exe
HydraAddExecutable(NovoMultiFit NovoMultiFit.cpp)
if(ROOT_FOUND)
    target_link_libraries(NovoMultiFit INTERFACE Root)
endif()
