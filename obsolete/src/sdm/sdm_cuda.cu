// including it first not to require pthread option to nvcc
#include <blitz/array.h>

#define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_CUDA
#define ICICLE_THRUST_DEVICE_SYSTEM sdm::cuda
#include "sdm.cpp"
