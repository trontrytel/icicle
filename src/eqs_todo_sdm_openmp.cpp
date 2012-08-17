#if defined(USE_THRUST)
#  define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP
#else
#  define THRUST_DEVICE_SYSTEM 0 // TODO!!!!
#endif
#include "eqs_todo_sdm.inc"

