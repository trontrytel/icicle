#include "../cfg/cfg_thrust.hpp"
#if defined(_OPENMP) && defined(USE_THRUST)
#  define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP
#  define ICICLE_THRUST_DEVICE_SYSTEM sdm::openmp
#  include "sdm.cpp"
#endif
