include_directories(${HYDRA_INCLUDE_DIR})

add_subdirectory(pybind11)
pybind11_add_module(hydra_example main.cpp)

set_target_properties(hydra_example PROPERTIES 
        COMPILE_FLAGS "-DTHRUST_HOST_SYSTEM=THRUST_HOST_SYSTEM_CPP -DTHRUST_DEVICE_SYSTEM=THRUST_DEVICE_SYSTEM_CPP ${HYDRA_CXX_FLAGS}")
