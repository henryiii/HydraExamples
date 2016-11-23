
# List the final name and cpp for each exe

# This was buggy until a recent version of CMake
find_package(Threads)

find_package(ROOT)
if(ROOT_FOUND)
    add_library(Root INTERFACE)
    set(root_libs
        ${ROOT_Hist_LIBRARY}
        ${ROOT_Minuit2_LIBRARY}
        ${ROOT_Core_LIBRARY}
        ${ROOT_MathCore_LIBRARY}
        ${ROOT_Thread_LIBRARY}
        ${ROOT_Physics_LIBRARY}
        ${ROOT_Gpad_LIBRARY})
    message(STATUS "Setting root libraries: ${root_libs}")
    target_link_libraries(Root INTERFACE ${root_libs} Threads::Threads)
    target_include_directories(Root INTERFACE ${ROOT_INCLUDE_DIR})
    target_compile_definitions(Root INTERFACE FOUND_ROOT)
endif()

add_subdirectory(phsp)
add_subdirectory(simple_fit)
add_subdirectory(novo_fit)
