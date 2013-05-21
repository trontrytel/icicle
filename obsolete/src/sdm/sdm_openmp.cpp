#include "../cfg/cfg_thrust.hpp"
#include "../cfg/cfg_boost_odeint.hpp"
#if defined(_OPENMP) && defined(USE_THRUST) && defined(USE_BOOST_ODEINT)
#  define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP
#  define ICICLE_THRUST_DEVICE_SYSTEM sdm::openmp
#  include "sdm.cpp"
#endif
