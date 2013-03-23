## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date November 2011 - February 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

# http://www.cmake.org/Wiki/CMake_Version_Compatibility_Matrix/StandardCMakeModules
# cmake_minimum_required(VERSION 2.6) # for simplified if-constructs
# cmake_minimum_required(VERSION 2.8) # for FindCuda

cmake_minimum_required(VERSION 2.8.6) # for CMakePushCheckState
project(icicle CXX)

# CMake includes
include(FeatureSummary)
include(CheckCXXSourceCompiles) 
include(CheckCXXSourceRuns)
include(CMakePushCheckState)

find_package(OpenMP)

set(CXX11_FLAG "-std=c++0x ${OpenMP_CXX_FLAGS}") # C++11 (name understandable for pre-2011 compilers)
set(CMAKE_CXX_FLAGS "${OpenMP_CXX_FLAGS}")
set(CMAKE_CXX_FLAGS_RELEASE "${CXX11_FLAG} -DNDEBUG -Wfatal-errors" CACHE STRING "Release options." FORCE)
if(CMAKE_COMPILER_IS_GNUCXX)
#  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -Ofast" CACHE STRING "Release options." FORCE)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2" CACHE STRING "Release options." FORCE)
#  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3" CACHE STRING "Release options." FORCE)
else()
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3" CACHE STRING "Release options." FORCE)
endif()
set(CMAKE_CXX_FLAGS_DEBUG "${CXX11_FLAG} -DBZ_DEBUG" CACHE STRING "Debug options." FORCE)
unset(CXX11_FLAG)

if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
  # TODO: use the actual flags for tests
elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  # TODO: use the actual flags for tests
else()
  message(FATAL_ERROR "Please specify the build type, i.e.:
    \$ cmake . -DCMAKE_BUILD_TYPE=Debug   \# (CXX_FLAGS: \"${CMAKE_CXX_FLAGS_DEBUG} \$CMAKE_CXX_FLAGS\")
    or
    \$ cmake . -DCMAKE_BUILD_TYPE=Release \# (CXX_FLAGS: \"${CMAKE_CXX_FLAGS_RELEASE} \$CMAKE_CXX_FLAGS\")")
endif()

set(TMP "src/cfg_types.hpp") # // TODO: other cfg_*
foreach(d .;doc;tests;src)
  set(TMP "${TMP};${d}/Makefile;${d}/CMakeCache.txt;${d}/cmake_install.cmake;${d}/CMakeFiles;${d}/Testing;${d}/CTestTestfile.cmake")
endforeach()
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${TMP}")

# describing available options
function(icicle_option var desc) # functions use local scope, macros use global
  option(${var} ${desc} ${${var}})
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

set(tmp "")
set(tmp "${tmp};src/slv_serial.cpp")
set(tmp "${tmp};src/slv_parallel_serial.cpp")
set(tmp "${tmp};src/slv_parallel_distmem_fork.cpp")
set(tmp "${tmp};src/slv_parallel_distmem_mpi.cpp")
set(tmp "${tmp};src/eqs/eqs_scalar_advection.cpp")
set(tmp "${tmp};src/eqs/eqs_isentropic.cpp")
set(tmp "${tmp};src/eqs/eqs_todo_bulk_ode.cpp")
set(tmp "${tmp};src/eqs/eqs_todo_mm.cpp")
set(tmp "${tmp};src/eqs/eqs_todo_sdm.cpp")
set(tmp "${tmp};src/sdm/sdm_openmp.cpp")
set(tmp "${tmp};src/adv_leapfrog.cpp")
set(tmp "${tmp};src/adv_upstream.cpp")
set(tmp "${tmp};src/adv_mpdata.cpp")
set(tmp "${tmp};src/adv_mpdata_fct.cpp")
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
set(tmp "${tmp};src/icicle.cpp")

add_executable(icicle ${tmp})

# FLOAT, DOUBLE, LDOUBLE and FLOAT128
set(USE_FLOAT FALSE)
set(USE_DOUBLE TRUE)
set(USE_LDOUBLE FALSE)
set(USE_FLOAT128 FALSE)      

# checking if fork() requires librt
set(code "#include <boost/interprocess/ipc/message_queue.hpp>\n int main() {  boost::interprocess::message_queue q(boost::interprocess::open_or_create, \"n\", 1, sizeof(bool))\; }")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
icicle_check(status "ipc" WARNING "" "${code}")
if(NOT status)
  cmake_push_check_state()
  set(CMAKE_REQUIRED_LIBRARIES "-lrt")
  icicle_check(status "ipc_librt" FATAL_ERROR "ipc not usable" "${code}")
  target_link_libraries(icicle ${CMAKE_REQUIRED_LIBRARIES})
  cmake_pop_check_state()
endif()
unset(code)

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

find_package(Boost COMPONENTS iostreams system program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})
target_link_libraries(icicle ${Boost_LIBRARIES})

# Boost.lexical_cast (mandatory)
icicle_check(status "boost_lexical_cast" FATAL_ERROR ""
  "#include <boost/lexical_cast.hpp>\n int main() {}"
)

# Boost.Units (mandatory)
icicle_check(status "boost_units" FATAL_ERROR ""
  "#include <boost/units/systems/si.hpp>\n int main() {}"
)

## Boost.program_options (mandatory)
#icicle_check_boost_lib(
#  status "boost_program_options" FATAL_ERROR ""
#  "#include <boost/program_options.hpp>\n int main() {boost::program_options::options_description d\;}"
#)

# Boost.Timer (optional) 
set(USE_BOOST_TIMER FALSE) 
if (DEFINED BOOST_TIMER AND NOT BOOST_TIMER)
  message(STATUS "Building without Boost.Timer support (-DBOOST_TIMER=NO specified or CMake have not found it)")
else()
  icicle_check_boost_lib(
    USE_BOOST_TIMER "boost_timer" WARNING "Simulation timing will not be available."
    "#include <boost/timer/timer.hpp>\n int main() { boost::timer::auto_cpu_timer t\; }" 
    "boost_system"
  )
endif()                                                                
icicle_option(USE_BOOST_TIMER "Simulation timing using Boost.Timer")

# Boost.MPI (optional)
set(USE_BOOST_MPI FALSE)
if (DEFINED BOOST_MPI AND NOT BOOST_MPI)
  message(STATUS "Building without Boost.MPI support (-DBOOST_MPI=NO specified or CMake have not found it)")
else()
  icicle_check_boost_lib(
    USE_BOOST_MPI "boost_mpi" WARNING 
    "MPI parallelisation will not be available. Have you used a proper MPI compiler wrapper? (e.g. -DCMAKE_CXX_COMPILER=mpic++)"
    "#include <boost/mpi.hpp>\n int main() {boost::mpi::communicator().send(1,0,0.)\;}"
  )
endif()
icicle_option(USE_BOOST_MPI "MPI paralellisation using Boost.MPI")

# Boost.odeint (optional)
set(USE_BOOST_ODEINT FALSE)
if (DEFINED BOOST_ODEINT AND NOT BOOST_ODEINT) 
  message(STATUS "Building without Boost.odeint support (-DODEINT=NO specified or CMake have not found it)")
else()
  icicle_check(USE_BOOST_ODEINT "boost_odeint" WARNING 
    "Features relying on ODE integration will not be available."
    "#include <boost/numeric/odeint.hpp>\n int main() {}" 
  )
  if (USE_BOOST_ODEINT)
    set(TMP ${CMAKE_CXX_FLAGS}) 
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    icicle_check(USE_BOOST_ODEINT "boost_fusion" WARNING 
      "Boost.fusion seems not functional (buggy Boost 1.50?)"
      "#include <boost/fusion/container.hpp>\n int main() {}" 
    )
    set(CMAKE_CXX_FLAGS ${TMP}) 
  endif()
endif()
icicle_option(USE_BOOST_ODEINT "ODE integration with Boost.odeint")

# NetCDF (optional)
set(USE_NETCDF FALSE)
if(DEFINED NETCDF AND NOT NETCDF) 
  message(STATUS "Building without netCDF support (-DNETCDF=NO specified or CMake did not found it)")
else()
  cmake_push_check_state()
  set(CMAKE_REQUIRED_LIBRARIES "-lnetcdf_c++4")
  icicle_check(USE_NETCDF "netcdf_cxx4" WARNING 
    "NetCDF input/output will not be available."
    "#include <netcdf>\n int main() { netCDF::NcFile *f\; }" 
  )
  if(USE_NETCDF) 
    target_link_libraries(icicle ${CMAKE_REQUIRED_LIBRARIES})
  endif()
  cmake_pop_check_state()
endif() 
icicle_option(USE_NETCDF "netCDF input/output")

# Thrust (optional)
if(DEFINED THRUST AND NOT THRUST) 
  message(STATUS "Building without Thrust support (-DTHRUST=NO specified or CMake have not found it)")
else()
  icicle_check(USE_THRUST "thrust" WARNING "Features relying on Thrust will not be available." "
    #define THRUST_DEVICE_SYSTEM THRUST_DEVICE_SYSTEM_OMP
    #include <thrust/device_vector.h>
    int main() { 
      thrust::device_vector<int> v(1)\; 
    }
  ")
endif() 
icicle_option(USE_THRUST "Thrust (OpenMP parallelisation)")


# TODO: document
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug") 
  set(INFO_CXX_FLAGS "${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_CXX_FLAGS}")
elseif("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
  set(INFO_CXX_FLAGS "${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_CXX_FLAGS}")
else()
  #TODO: assert false...
endif()
execute_process(COMMAND "whoami" OUTPUT_VARIABLE INFO_WHOAMI OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND "hostname" OUTPUT_VARIABLE INFO_HOSTNAME OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND "${CMAKE_CXX_COMPILER}" "--version" COMMAND "head" "-1" OUTPUT_VARIABLE INFO_CXX_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/cfg/cfg_types.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/src/cfg/cfg_types.hpp)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/cfg/cfg_info.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/src/cfg/cfg_info.hpp)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/cfg/cfg_netcdf.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/src/cfg/cfg_netcdf.hpp)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/cfg/cfg_thrust.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/src/cfg/cfg_thrust.hpp)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/cfg/cfg_boost_mpi.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/src/cfg/cfg_boost_mpi.hpp)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/cfg/cfg_boost_odeint.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/src/cfg/cfg_boost_odeint.hpp)
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/src/cfg/cfg_boost_timer.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/src/cfg/cfg_boost_timer.hpp)

enable_testing()

add_subdirectory(tests)
add_subdirectory(examples)

message(STATUS "-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --")
message(STATUS "-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --")
feature_summary(WHAT ENABLED_FEATURES DESCRIPTION "Enabled Features (to disable use -DFEATURE=NO where USE_FEATURE=ON):")
feature_summary(WHAT DISABLED_FEATURES DESCRIPTION "Disabled Features (to enable use -DFEATURE=ON where USE_FEATURE=OFF):")
message(STATUS "Compilation set-up:
compiler: ${CMAKE_CXX_COMPILER}
version: ${INFO_CXX_VERSION}
flags: ${INFO_CXX_FLAGS}")
message(STATUS "-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --")
message(STATUS "-- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --")
