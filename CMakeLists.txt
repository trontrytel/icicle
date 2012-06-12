## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date November 2011 - February 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

# http://www.cmake.org/Wiki/CMake_Version_Compatibility_Matrix/StandardCMakeModules
# cmake_minimum_required(VERSION 2.6) # for simplified if-constructs
cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState

# TODO: -fprofile-use, -fprofile-generate
set(TMP "-std=c++0x") # C++11 (name understandable for pre-2011 compilers)
set(CMAKE_CXX_FLAGS_RELEASE "${TMP} -DNDEBUG -Wfatal-errors" CACHE STRING "Release options." FORCE)
if (CMAKE_COMPILER_IS_GNUCXX)
  set(CMAKE_CXX_FLAGS_RELEASE "${TMP} -Ofast" CACHE STRING "Release options." FORCE)
else()
  set(CMAKE_CXX_FLAGS_RELEASE "${TMP} -O3" CACHE STRING "Release options." FORCE)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "${TMP} -DBZ_DEBUG" CACHE STRING "Debug options." FORCE)
unset(TMP)

if(NOT (
  "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR 
  "${CMAKE_BUILD_TYPE}" STREQUAL "Release" 
))
  message(FATAL_ERROR "Please specify the build type, i.e.:
    \$ cmake . -DCMAKE_BUILD_TYPE=Debug   \# (CXX_FLAGS=\"\$CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_DEBUG}\")
    or
    \$ cmake . -DCMAKE_BUILD_TYPE=Release \# (CXX_FLAGS=\"\$CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS_RELEASE}\")")
endif()

project(icicle CXX)

set(TMP "src/cfg.hpp")
foreach(d .;doc;tests;src)
  set(TMP "${TMP};${d}/Makefile;${d}/CMakeCache.txt;${d}/cmake_install.cmake;${d}/CMakeFiles;${d}/Testing;${d}/CTestTestfile.cmake")
endforeach()
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${TMP}")

# CMake includes
include(FeatureSummary)
include(CheckCXXSourceCompiles) 
include(CheckCXXSourceRuns)
include(CMakePushCheckState)

# describing available options
function(icicle_option var desc onoff) # functions use local scope (macros use global)
  option(${var} ${desc} ${onoff})
  add_feature_info(${var} ${var} ${desc})
endfunction()

# checking if a given piece of C++ code works
function(icicle_check stat what bhvr msg code) # functions use local scope (macros use global)
  # we use CheckCXXSource* macros not to require a C compiler
  check_cxx_source_compiles(${code} "icicle_${what}") # ${what} is printed
  if(NOT icicle_${what}) 
    if(NOT "${bhvr}" STREQUAL "HIDDEN")
      set(bar "**************************************************************************")
      set(msg "${bar} Test ${what} failed. ${msg}")
      set(msg "${msg} Consult CMakeFiles/CMakeError.log for details.")
      set(msg "${msg} After fixing the problem, delete CMakeCache.txt and rerun cmake)")
      set(msg "${msg} ${bar}")
      message(${bhvr} ${msg})
    endif()
    set(${stat} FALSE PARENT_SCOPE)
  else()
    set(${stat} TRUE PARENT_SCOPE)
  endif() 
endfunction()

# checking for a Boost library
function(icicle_check_boost_lib stat libname bhvr msg code)
  cmake_push_check_state()
  foreach(sfx "-mt" "") # tagged installation layout issue
    set(libs "-l${libname}${sfx}")
    if(ARGC EQUAL 6) 
      set(libs "${libs} -l${ARGV5}${sfx}")
    endif()
    set(CMAKE_REQUIRED_LIBRARIES "${libs}")
    if("${sfx}" STREQUAL "-mt")
      icicle_check(status "${libname}" "HIDDEN" "${msg}" "${code}")
      if(status)
        break()
      endif()
    else()
      message(STATUS "Retrying without the -mt suffix...")
      unset("icicle_{$libname}" CACHE)
      icicle_check(status "${libname}" ${bhvr} "${msg}" "${code}")
    endif()
  endforeach()
  if(status)
    target_link_libraries(icicle ${CMAKE_REQUIRED_LIBRARIES})
  endif()
  cmake_pop_check_state()
  set(${stat} ${status} PARENT_SCOPE)
endfunction()

# FLOAT, DOUBLE, LDOUBLE and FLOAT128
icicle_option(FLOAT "'float' type support (i.e. single precision, 32 bits)" ON)
icicle_option(DOUBLE "'double' type support (i.e. double precision, 64 bits)" ON)
icicle_option(LDOUBLE "'long double' type support (compiler-dependent, 80 bits for GCC on x86_64)" ON)
#icicle_option(FLOAT128 "'__float128' type support (GNU extension, 128 bits)" ON) #TODO: check if __float128 is actually supported by the compiler...

set(tmp "")
set(tmp "${tmp};src/adv_leapfrog.cpp")
set(tmp "${tmp};src/adv_upstream.cpp")
set(tmp "${tmp};src/adv_mpdata.cpp")
set(tmp "${tmp};src/adv_mpdata_fct.cpp")
set(tmp "${tmp};src/slv_serial.cpp")
set(tmp "${tmp};src/slv_parallel_serial.cpp")
set(tmp "${tmp};src/slv_parallel_openmp.cpp")
set(tmp "${tmp};src/slv_parallel_threads.cpp")
set(tmp "${tmp};src/slv_parallel_distmem_fork.cpp")
set(tmp "${tmp};src/slv_parallel_distmem_mpi.cpp")
set(tmp "${tmp};src/eqs_harmonic_oscillator.cpp")
set(tmp "${tmp};src/eqs_shallow_water.cpp")
set(tmp "${tmp};src/eqs_scalar_advection.cpp")
set(tmp "${tmp};src/eqs_isentropic.cpp")
set(tmp "${tmp};src/eqs_todo_sdm.cpp")
set(tmp "${tmp};src/eqs_todo_bulk_ode.cpp")
set(tmp "${tmp};src/vel_momeq_extrapol.cpp")
set(tmp "${tmp};src/vel_func_test.cpp")
set(tmp "${tmp};src/vel_func_stream_rasinski.cpp")
set(tmp "${tmp};src/inf.cpp")
set(tmp "${tmp};src/ini_netcdf.cpp")
set(tmp "${tmp};src/ini_func_cone.cpp")
set(tmp "${tmp};src/ini_func_gauss.cpp")
set(tmp "${tmp};src/ini_func_boxcar.cpp")
set(tmp "${tmp};src/ini_func_pks_wwg_1989.cpp")
set(tmp "${tmp};src/out_netcdf.cpp")
set(tmp "${tmp};src/out_gnuplot.cpp")
set(tmp "${tmp};src/out_debug.cpp")
set(tmp "${tmp};src/icicle.cpp")

if(NOT FLOAT)
  message(STATUS "Building without float support (-DFLOAT=NO specified)")
else()
  set(USE_FLOAT TRUE)
endif()
if(NOT DOUBLE)
  message(STATUS "Building without double support (-DDOUBLE=NO specified)")
else()
  set(USE_DOUBLE TRUE)
endif()    
if(NOT LDOUBLE)
  message(STATUS "Building without long double support (-DLDOUBLE=NO specified)")
else()
  set(USE_LDOUBLE TRUE)
endif()
#      if(NOT "${FLOAT128}" STREQUAL "NO")
#        set(USE_FLOAT128 TRUE)
#      endif()
# TODO __fp16 (http://gcc.gnu.org/onlinedocs/gcc/Half_002dPrecision.html)
if("${tmp}" STREQUAL "src/icicle.cpp")
  message(FATAL_ERROR "-DFLOAD=NO -DDOUBLE=NO -DLDBOULBE=NO and -DFLOAT128=NO may not be specified all at once!")
endif() 
add_executable(icicle ${tmp})

# Blitz++ (mandatory)
target_link_libraries(icicle debug "-lblitz") # TODO: check if debugging actually works, if not suggest recompiling blitz with --enable-shared
icicle_check(status "blitz" FATAL_ERROR ""
  "#include <blitz/array.h>\n int main() {}"
)
check_cxx_source_runs("
  #undef NDEBUG
  #include <blitz/array.h>
  int main() {
    assert(std::string(BZ__compiler_name) == std::string(\"${CMAKE_CXX_COMPILER}\"));
  }
" icicle_blitz_compiler)
if (NOT icicle_blitz_compiler)
  message(WARNING "Blitz++ seem not to have been configured with the same compiler!
    Consider reconfiguring Blitz++ with CXX=${CMAKE_CXX_COMPILER}"
  )
endif()

# Boost.lexical_cast (mandatory)
icicle_check(status "boost_lexical_cast" FATAL_ERROR ""
  "#include <boost/lexical_cast.hpp>\n int main() {}"
)

# Boost.Units (mandatory)
icicle_check(status "boost_units" FATAL_ERROR ""
  "#include <boost/units/systems/si.hpp>\n int main() {}"
)

# Boost.program_options (mandatory)
icicle_check_boost_lib(
  status "boost_program_options" FATAL_ERROR ""
  "#include <boost/program_options.hpp>\n int main() {boost::program_options::options_description d\;}"
)

# Boost.Thread (optional)
if (DEFINED THREADS AND NOT THREADS)
  message(STATUS "Building without Boost.Thread support (-DTHREADS=NO specified)")
else()
  icicle_check_boost_lib(
    status "boost_thread" WARNING "Threads parallelisation will not be available."
    "#include <boost/thread.hpp>\n int main() { boost::thread th\; }"
  )
  set(USE_BOOST_THREAD ${status})
  if(status)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  endif()
endif()
icicle_option(THREADS "threads paralellisation using Boost.Thread" "${status}")

# Boost.Timer (optional) 
if (DEFINED TIMER AND NOT TIMER)
  message(STATUS "Building without Boost.Timer support (-DTIMER=NO specified)")
else()
  icicle_check_boost_lib(
    status "boost_timer" WARNING "Simulation timing will not be available."
    "#include <boost/timer/timer.hpp>\n int main() { boost::timer::auto_cpu_timer t\; }" 
    "boost_system"
  )
  set(USE_BOOST_TIMER ${status}) # TODO: the same name as for -DTIMER=NO ?
endif()                                                                
icicle_option(TIMER "Simulation timing using Boost.Timer" "${status}")

# Boost.MPI (optional)
if (DEFINED MPI AND NOT MPI)
  message(STATUS "Building without Boost.MPI support (-DMPI=NO specified)")
else()
  icicle_check_boost_lib(
    status "boost_mpi" WARNING 
    "MPI parallelisation will not be available. Have you used a proper MPI compiler wrapper? (e.g. -DCMAKE_CXX_COMPILER=mpic++)"
    "#include <boost/mpi.hpp>\n int main() {boost::mpi::communicator().send(1,0,0.)\;}"
  )
  set(USE_BOOST_MPI ${status})
endif()
icicle_option(MPI "MPI paralellisation using Boost.MPI" "${status}")

# Boost.odeint (optional)
if (DEFINED ODEINT AND NOT ODEINT) 
  message(STATUS "Building without Boost.odeint support (-DODEINT=NO specified)")
else()
  icicle_check(status "boost_odeint" WARNING 
    "Features relying on ODE integration will not be available."
    "#include <boost/numeric/odeint.hpp>\n int main() {}" 
  )
  set(USE_BOOST_ODEINT ${status})
endif()
icicle_option(ODEINT "ODE integration with Boost.odeint" "${status}")

# OpenMP (optional, not using FIND_PACKAGE(OpenMP) as it relies on a C compiler)
if (DEFINED OPENMP AND NOT OPENMP)   
  message(STATUS "Building without OpenMP support (-DOPENMP=NO specified)")
else()
  set(TMP ${CMAKE_CXX_FLAGS}) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  icicle_check(status "openmp" WARNING 
    "OpenMP parallelisation will not be available." 
    "#include <omp.h>\n int main() { omp_set_num_threads(2)\;}"
  )
  if(NOT status) 
    set(CMAKE_CXX_FLAGS ${TMP}) 
  endif()   
endif()
icicle_option(OPENMP "OpenMP paralellisation" "${status}")

# NetCDF (optional)
if(DEFINED NETCDF AND NOT NETCDF) 
  message(STATUS "Building without netCDF support (-DNETCDF=NO specified)")
else()
  cmake_push_check_state()
  set(CMAKE_REQUIRED_LIBRARIES "-lnetcdf_c++4")
  icicle_check(status "netcdf_cxx4" WARNING 
    "NetCDF input/output will not be available."
    "#include <netcdf>\n int main() { netCDF::NcFile *f\; }" 
  )
  set(USE_NETCDF ${status})
  if(status) 
    target_link_libraries(icicle ${CMAKE_REQUIRED_LIBRARIES})
  endif()
  cmake_pop_check_state()
endif() 
icicle_option(NETCDF "netCDF input/output" "${status}")

# Thrust (optional)
if(DEFINED THRUST AND NOT THRUST) 
  message(STATUS "Building without Thrust support (-DTHRUST=NO specified)")
else()
  icicle_check(status "thrust" WARNING "Features relying on Thrust will not be available." "
    #define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP
    #include <thrust/device_vector.h>
    int main() { 
      thrust::device_vector<int> v(1)\; 
    }
  ")
  set(USE_THRUST ${status})
  set(THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP) #TODO: use somehow the information on OpenMP availability!
endif() 
icicle_option(THRUST "Thrust (OpenMP/CUDA parallelisation)" "${status}")

# TODO: document
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug") 
  set(INFO_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  set(INFO_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
else()
  set(INFO_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()
execute_process(COMMAND "whoami" OUTPUT_VARIABLE INFO_WHOAMI OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND "hostname" OUTPUT_VARIABLE INFO_HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND "${CMAKE_CXX_COMPILER}" "--version" COMMAND "head" "-1" OUTPUT_VARIABLE INFO_CXX_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/cfg.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/src/cfg.hpp)

#find_package(Doxygen)
#if (DOXYGEN_FOUND)
#  #add_subdirectory(doc) (FIXME - cmake-controlled doc build)
#else()
#  message(WARNING "Doxygent not found - 'make doc' will not be available.")
#endif()

enable_testing()

add_subdirectory(tests)
#add_subdirectory(src) 

message(STATUS "-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --")
feature_summary(WHAT ENABLED_FEATURES DESCRIPTION "Enabled Features (to disable use -DFEATURE=NO):")
feature_summary(WHAT DISABLED_FEATURES DESCRIPTION "Disabled Features (to enable use -DFEATURE=ON):")
message(STATUS "Compilation set-up:
compiler: ${CMAKE_CXX_COMPILER}
version: ${INFO_CXX_VERSION}
flags: ${INFO_CXX_FLAGS}")
message(STATUS "-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --")
