
# List the final name and cpp for each exe

# This was buggy until a recent version of CMake
find_package(Threads)

add_subdirectory(phsp)
add_subdirectory(simple_fit)
add_subdirectory(novo_fit)
