#pragma once

#include <string>

namespace hydra {

const std::string CURRENT_DEVICE = 
#if THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_CUDA
"cuda"
#elif THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_OMP
"omp"
#elif THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_TBB
"tbb"
#elif THRUST_DEVICE_SYSTEM==THRUST_DEVICE_SYSTEM_CPP
"cpp"
#else
"unknown"
#endif
;
}
