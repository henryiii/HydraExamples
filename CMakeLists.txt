cmake_minimum_required(VERSION 3.3)

# project name
project(HydraExamples CXX C)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR})

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING
      "Choose the type of build, options are: Debug Release
RelWithDebInfo MinSizeRel."
      FORCE)
endif(NOT CMAKE_BUILD_TYPE)


set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})
set(HYDRA_ROOT "../Hydra")

find_package(Hydra REQUIRED)
find_package(ROOT)

include_directories(include)

add_subdirectory(src)
option(PYTHON_BIND OFF)
if(PYTHON_BIND)
    add_subdirectory(python)
endif()

