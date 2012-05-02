## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date November 2011 - February 2012
#  @section LICENSE
#    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)

cmake_minimum_required(VERSION 2.6) # for simplified if-constructs
#TODO: check which version is needed for FeatureSummary

# TODO: -fprofile-use, -fprofile-generate
set(TMP "-std=c++0x -DTHRUST_DEVICE_BACKEND=THRUST_DEVICE_BACKEND_OMP") # C++11 (name understandable for pre-2011 compilers)
#set(TMP "${TMP} -Wfatal-errors") # stop on first error
set(CMAKE_CXX_FLAGS_RELEASE "${TMP} -Ofast -DNDEBUG" CACHE STRING "Release options." FORCE)
set(CMAKE_CXX_FLAGS_DEBUG "${TMP} -DBZ_DEBUG" CACHE STRING "Debug options." FORCE)
unset(TMP)

if(NOT (
  "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR 
  "${CMAKE_BUILD_TYPE}" STREQUAL "Release" 
))
  message(FATAL_ERROR "Please specify the build type, i.e.:
    \$ cmake . -DCMAKE_BUILD_TYPE=Debug   \# (CXX_FLAGS=\"\$CXXFLAGS ${CMAKE_CXX_FLAGS_DEBUG}\")
    or
    \$ cmake . -DCMAKE_BUILD_TYPE=Release \# (CXX_FLAGS=\"\$CXXFLAGS ${CMAKE_CXX_FLAGS_RELEASE}\")")
endif()

project(icicle CXX)

set(TMP "src/cfg.hpp")
foreach(d .;doc;tests;src)
  set(TMP "${TMP};${d}/Makefile;${d}/CMakeCache.txt;${d}/cmake_install.cmake;${d}/CMakeFiles;${d}/Testing;${d}/CTestTestfile.cmake")
endforeach()
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${TMP}")

# we use CheckCXXSource* macros only not to require a C compiler
include(CheckCXXSourceCompiles) 
include(CheckCXXSourceRuns)

# describing available options
include(FeatureSummary)
macro(icicle_option var desc onoff)
  option(${var} ${desc} ${onoff})
  add_feature_info(${var} ${var} ${desc})
endmacro()

#### <-DFLOAT=NO, -DDOUBLE=NO, -DLDOUBLE=NO> ###########################
#TODO: check if __float128 is actually supported by the compiler...
icicle_option(FLOAT "'float' type support (i.e. single precision, 32 bits)" ON)
icicle_option(DOUBLE "'double' type support (i.e. double precision, 64 bits)" ON)
icicle_option(LDOUBLE "'long double' type support (compiler-dependent, 80 bits for GCC on x86_64)" ON)
#icicle_option(FLOAT128 "'__float128' type support (GNU extension, 128 bits)" ON)
if(NOT "${FLOAT}" STREQUAL "NO")                                       
  set(USE_FLOAT TRUE)
  if(NOT "${DOUBLE}" STREQUAL "NO")
    set(USE_DOUBLE TRUE)
    if(NOT "${LDOUBLE}" STREQUAL "NO")
      set(USE_LDOUBLE TRUE)
#      if(NOT "${FLOAT128}" STREQUAL "NO")
#        set(USE_FLOAT128 TRUE)
#        add_executable(icicle src/icicle.cpp src/mdl_flt.cpp src/mdl_dbl.cpp src/mdl_ldb.cpp src/mdl_128.cpp)
#      else()
        add_executable(icicle src/icicle.cpp src/mdl_flt.cpp src/mdl_dbl.cpp src/mdl_ldb.cpp)
#      endif()
    else()
      add_executable(icicle src/icicle.cpp src/mdl_flt.cpp src/mdl_dbl.cpp)
    endif()
  else()
    add_executable(icicle src/icicle.cpp src/mdl_flt.cpp)
  endif()
elseif(NOT "${DOUBLE}" STREQUAL "NO")
  set(USE_DOUBLE TRUE)
  if(NOT "${LDOUBLE}" STREQUAL "NO")
    set(USE_LDOUBLE TRUE)
#    if(NOT "${FLOAT128}" STREQUAL "NO")
#      add_executable(icicle src/icicle.cpp src/mdl_dbl.cpp src/mdl_ldb.cpp src/mdl_128.cpp)
#    else()
      add_executable(icicle src/icicle.cpp src/mdl_dbl.cpp src/mdl_ldb.cpp)
#    endif()
  else()
    add_executable(icicle src/icicle.cpp src/mdl_dbl.cpp)
  endif()
elseif(NOT "${LDOUBLE}" STREQUAL "NO")
  set(USE_LDOUBLE TRUE)
#  if(NOT "${FLOAT128}" STREQUAL "NO")
#    add_executable(icicle src/mdl_ldb.cpp src/mdl_128.cpp)
#  else()
    add_executable(icicle src/mdl_ldb.cpp)
#  endif()
#elseif(NOT "${FLOAT128}" STREQUAL "NO")
#  add_executable(icicle src/mdl_128.cpp)
else()
  message(FATAL_ERROR "-DFLOAD=NO -DDOUBLE=NO -DLDBOULBE=NO and -DFLOAT128=NO may not be specified all at once!")
endif() 
if("${CMAKE_BUILD_TYPE}" STREQUAL "Debug") 
  target_link_libraries(icicle "-lblitz") # TODO: could be done with the "debug" keyword to target_link_libraries()
  # TODO: check if debugging actually works, if not suggest recompiling blitz with --enable-shared
endif()
########################################################################

#### <Boost.units> (mandatory) #########################################
check_cxx_source_compiles(
  "#include <boost/units/systems/si.hpp>\n int main() {}" 
  Find_Boost_Units
)
if(NOT Find_Boost_Units) 
  message(FATAL_ERROR "Boost.Units not found! 
    (cf. CMakeFiles/CMakeError.log, fix, rm CMakeCache.txt and rerun)"
  )
endif() 
#### </Boost.units> ####################################################

#### <Boost.lexical_cast> (mandatory) ##################################
check_cxx_source_compiles(
  "#include <boost/lexical_cast.hpp>\n int main() {}" 
  Find_Boost_Lexical_Cast
)
if(NOT Find_Boost_Lexical_Cast) 
  message(FATAL_ERROR "Boost.lexical_cast not found! 
    (cf. CMakeFiles/CMakeError.log, fix, rm CMakeCache.txt and rerun)"
  )
endif() 
#### </Boost.lexical_cast> #############################################

#### <Blitz> (mandatory) ###############################################
check_cxx_source_compiles( 
  "#include <blitz/array.h>\n int main() {}"
  Find_Blitz
)
if(NOT Find_Blitz) 
  message(FATAL_ERROR "Blitz++ not found!
    (cf. CMakeFiles/CMakeError.log, fix, rm CMakeCache.txt and rerun)"
  )
endif()                                                                
check_cxx_source_runs("
  #include <blitz/array.h>
  int main() { 
    assert(std::string(BZ__compiler_name) == std::string(\"${CMAKE_CXX_COMPILER}\"));
  }
  " Blitz_Configured_With_The_Same_Compiler
)
if (NOT Blitz_Configured_With_The_Same_Compiler)
  message(WARNING "Blitz++ seem not to have been configured with the same compiler!
    Consider reconfiguring Blitz++ with CXX=${CMAKE_CXX_COMPILER} 
    (cf. CMakeFiles/CMakeError.log, fix, rm CMakeCache.txt and rerun)"
  )
endif() 
#### </Blitz> ##########################################################

#### <Boost.program_options> (mandatory) ###############################
# TODO: could it just be set(CMAKE_REQUIRED_LIBRARIES "-laqq-mt -laqq") ???
foreach(sfx "-mt" "") # tagged installation layout issue               
  set(CMAKE_REQUIRED_LIBRARIES "-lboost_program_options${sfx}")  # TODO? should we use CMAKE_PUSH_CHECK_STATE() and CMAKE_POP_CHECK_STATE() ??
  check_cxx_source_compiles(
    "#include <boost/program_options.hpp>\n int main() { 
      boost::program_options::options_description d; 
    }" 
    Find_Boost_Program_Options
  )
  if (Find_Boost_Program_Options) 
    break()
  else()
    message(STATUS "Retrying without the -mt suffix...")
    unset(Find_Boost_Program_Options CACHE)
  endif()
endforeach()
if(NOT Find_Boost_Program_Options) 
  message(FATAL_ERROR "Boost.Program_Options not found! 
    (cf. CMakeFiles/CMakeError.log, fix, rm CMakeCache.txt and rerun)"
  )
endif()
target_link_libraries(icicle ${CMAKE_REQUIRED_LIBRARIES})
#### </Boost.program_options> ##########################################

#### <OpenMP> (optional) ###############################################
# Not using FIND_PACKAGE(OpenMP) as it relies on C compiler rather than C++ one (TODO: support other options for other compilers?)
if ("${OPENMP}" STREQUAL "NO")   
  message(STATUS "Building without OpenMP (-DOPENMP=NO specified)")
else()
  set(TMP ${CMAKE_CXX_FLAGS}) 
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
  check_cxx_source_compiles(
    "#include <omp.h>\n int main() { omp_set_num_threads(2); }" 
    Find_OpenMP
  )
  if(NOT Find_OpenMP)
    message(WARNING "OpenMP not found (cf. CMakeFiles/CMakeError.log) 
      OpenMP parallelisation will not be available.")
    set(CMAKE_CXX_FLAGS ${TMP}) 
  endif()   
endif()
icicle_option(OPENMP "OpenMP paralellisation" "${Find_OpenMP}")
#### </OpenMP> #########################################################

#### <NetCDF> (optional) ###############################################
if("${NETCDF}" STREQUAL "NO") 
  message(STATUS "Building without netCDF (-DNETCDF=NO specified)")
else()
  set(CMAKE_REQUIRED_LIBRARIES "-lnetcdf_c++4")
  check_cxx_source_compiles(
    "#include <netcdf>\n int main() { netCDF::NcFile *f; }" 
    Find_NetCDF
  )
  if(NOT Find_NetCDF) 
    message(WARNING "
      NetCDF >=4.1.1 not found! (cf. CMakeFiles/CMakeError.log)
      Without it the NetCDF output will not be available.
      After installing it erase the CMakeCache.txt before re-running cmake
    ")
  else()
    set(USE_NETCDF ${Find_NetCDF})
    target_link_libraries(icicle ${CMAKE_REQUIRED_LIBRARIES})
  endif()
endif() 
icicle_option(NETCDF "netCDF output" "${Find_NetCDF}")
#### </NetCDF> #########################################################

#### <Boost.Thread> (optional) #########################################
if ("${THREADS}" STREQUAL "NO") 
  message(STATUS "Building without Boost.Thread (-DTHREADS=NO specified)")
else()
  foreach(sfx "-mt" "") # tagged installation layout issue  
    set(CMAKE_REQUIRED_LIBRARIES "-lboost_thread${sfx}")
    check_cxx_source_compiles(
      "#include <boost/thread.hpp>\n int main() { boost::thread th; }" 
      Find_Boost_Thread
    )
    if (Find_Boost_Thread) 
      break()
    else()
      message(STATUS "Retrying without the -mt suffix...")
      unset(Find_Boost_Thread CACHE)
    endif()
  endforeach()
  if(NOT Find_Boost_Thread) 
    message(WARNING "
      Boost.Thread not found! (cf. CMakeFiles/CMakeError.log) 
      Threads parallelisation will not be available.
      After installing it erase the CMakeCache.txt before re-running cmake
    ")
  else()
    set(USE_BOOST_THREAD ${Find_Boost_Thread})
    target_link_libraries(icicle ${CMAKE_REQUIRED_LIBRARIES})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  endif()
endif()                                                                
icicle_option(THREADS "threads paralellisation using Boost.Thread" "${Find_Boost_Thread}")
#### </Boost.Thread> ###################################################

#### <Boost.MPI> (optional) ############################################
if ("${MPI}" STREQUAL "NO")
  message(STATUS "Building without Boost.MPI (-DMPI=NO specified)")
else()
  foreach(sfx "-mt" "") # tagged installation layout issue 
    set(CMAKE_REQUIRED_LIBRARIES "-lboost_mpi${sfx}")
    check_cxx_source_compiles(
      "#include <boost/mpi.hpp>
       int main() {boost::mpi::communicator().send(1,0,0.);}"
      Find_Boost_MPI
    )
    if (Find_Boost_MPI) 
      break()
    else()
      message(STATUS "Retrying without the -mt suffix...")
      unset(Find_Boost_MPI CACHE)
    endif()
  endforeach()
  if(NOT Find_Boost_MPI)
    message(WARNING "
      Boost.MPI not found! (cf. CMakeFiles/CMakeError.log)
      MPI parallelisation will not be available.
      Possible reasons:
      - a wrapper like mpic++ needs to be used (e.g. export CXX=mpic++)
      - Boost.MPI needs to be installed:
        (after installation delete CMakeCache.txt before running cmake)
    ")
  else()
    set(USE_BOOST_MPI ${Find_Boost_MPI})
    target_link_libraries(icicle ${CMAKE_REQUIRED_LIBRARIES})
  endif()
endif()
icicle_option(MPI "MPI paralellisation using Boost.MPI" "${Find_Boost_MPI}")
#### </Boost.MPI> ######################################################

#### <Boost.Timer> (optional) ##########################################
if ("${TIMER}" STREQUAL "NO")
  message(STATUS "Building without Boost.Timer (-DTIMER=NO specified)")
else()
  foreach(sfx "-mt" "") # tagged installation layout issue     
    set(CMAKE_REQUIRED_LIBRARIES "-lboost_timer${sfx} -lboost_system${sfx}")
    check_cxx_source_compiles(
      "#include <boost/timer/timer.hpp>
       int main() { boost::timer::auto_cpu_timer t; }" 
      Find_Boost_Timer
    )
    if (Find_Boost_Timer) 
      break()
    else()
      message(STATUS "Retrying without the -mt suffix...")
      unset(Find_Boost_Timer CACHE)
    endif()
  endforeach()
  if(NOT Find_Boost_Timer) 
    message(WARNING "
      Boost.Timer (or Boost.System) not found! (cf. CMakeFiles/CMakeError.log) 
      Simulation timing will not be available.
      (after installation delete CMakeCache.txt before running cmake)
    ")
  else()
    set(USE_BOOST_TIMER ${Find_Boost_Timer})
    target_link_libraries(icicle ${CMAKE_REQUIRED_LIBRARIES})
  endif()
endif()                                                                
icicle_option(TIMER "Simulation timing using Boost.Timer" "${Find_Boost_Timer}")
#### </Boost.Timer> ####################################################

#### <Boost.odeint> ####################################################
if ("${ODEINT}" STREQUAL "NO") 
  message(STATUS "Building without Boost.odeint (-DODEINT=NO specified)")
else()
  check_cxx_source_compiles(
    "#include <boost/numeric/odeint.hpp>\n int main() {}" 
    Find_Boost_odeint
  )
  if(NOT Find_Boost_odeint) 
    message(WARNING "
      Boost.odeint not found! (cf. CMakeFiles/CMakeError.log) 
      ODE integration will not be available.
      (after installation delete CMakeCache.txt before running cmake)
    ")
  else()
    set(USE_BOOST_ODEINT ${Find_Boost_odeint})
  endif() 
endif()
icicle_option(ODEINT "ODE integration with Boost.odeint" "${Find_Boost_odeint}")
#### </Boost.odeint> ###################################################

#### icicle options ####################################################
# TODO!
# icicle_option(EQS_ISENTROPIC "..." "${???}")
########################################################################

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
