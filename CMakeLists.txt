cmake_minimum_required(VERSION 3.3)

# project name
project(HydraTemplate CXX)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

SET(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(HYDRA_ROOT "../Hydra")

find_package(Hydra REQUIRED)

HydraSetup()

find_package(ROOT)

add_subdirectory(src)



