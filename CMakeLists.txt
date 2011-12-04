## @file
#  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
#  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
#  @copyright University of Warsaw
#  @date November 2011
#  @section LICENSE
#    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)

cmake_minimum_required(VERSION 2.6) # for simplified if-constructs

set(CMAKE_CXX_FLAGS_RELEASE "-O3 -ffast-math -DNDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g") # TODO: -DBZ_DEBUG (needs linking with libblitz?

if(NOT (
  "${CMAKE_BUILD_TYPE}" STREQUAL "Debug" OR 
  "${CMAKE_BUILD_TYPE}" STREQUAL "Release" OR
  "${CMAKE_BUILD_TYPE}" STREQUAL "None"
))
  message(FATAL_ERROR "Please specify the build type, i.e.:
    \$ cmake . -DCMAKE_BUILD_TYPE=Debug   \# (CXX_FLAGS=\"\$CXXFLAGS ${CMAKE_CXX_FLAGS_DEBUG}\")
    or
    \$ cmake . -DCMAKE_BUILD_TYPE=Release \# (CXX_FLAGS=\"\$CXXFLAGS ${CMAKE_CXX_FLAGS_RELEASE}\")
    or
    \$ cmake . -DCMAKE_BUILD_TYPE=None    \# (CXX_FLAGS=\$CXXFLAGS)")
endif()

project(icicles CXX)

set(TMP "src/config.hpp")
foreach(d .;doc;tests;src)
  set(TMP "${TMP};${d}/Makefile;${d}/CMakeCache.txt;${d}/cmake_install.cmake;${d}/CMakeFiles;${d}/Testing;${d}/CTestTestfile.cmake")
endforeach()
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${TMP}")

include(CheckCXXSourceCompiles) 

#### <Boost.units> (mandatory) #########################################
check_cxx_source_compiles(                                             #
  "#include <boost/units/systems/si.hpp>\n int main() {}" 
  Find_Boost_Units
)
if(NOT Find_Boost_Units) 
  message(FATAL_ERROR "Boost.Units not found! 
    (cf. CMakeFiles/CMakeError.log, fix, rm CMakeCache.txt and rerun)"
  )
endif()                                                                #
#### </Boost.units> ####################################################

#### <Boost.lexical_cast> (mandatory) ##################################
check_cxx_source_compiles(                                             #
  "#include <boost/lexical_cast.hpp>\n int main() {}" 
  Find_Boost_Lexical_Cast
)
if(NOT Find_Boost_Lexical_Cast) 
  message(FATAL_ERROR "Boost.lexical_cast not found! 
    (cf. CMakeFiles/CMakeError.log, fix, rm CMakeCache.txt and rerun)"
  )
endif()                                                                #
#### </Boost.lexical_cast> #############################################

#### <Blitz> (mandatory) ###############################################
check_cxx_source_compiles(                                             #
  "#include <blitz/array.h>\n int main() {}" 
  Find_Blitz
)
if(NOT Find_Blitz) 
  message(FATAL_ERROR "Blitz++ not found!
    (cf. CMakeFiles/CMakeError.log, fix, rm CMakeCache.txt and rerun)"
  )
endif()                                                                #
#### </Blitz> ##########################################################

#### <Boost.program_options> (mandatory) ###############################
set(CMAKE_REQUIRED_LIBRARIES "-lboost_program_options")                #
check_cxx_source_compiles(
  "#include <boost/program_options.hpp>\n int main() { 
    boost::program_options::options_description d; 
  }" 
  Find_Boost_Program_Options
)
if(NOT Find_Boost_Program_Options) 
  message(FATAL_ERROR "Boost.Program_Options not found! 
    (cf. CMakeFiles/CMakeError.log, fix, rm CMakeCache.txt and rerun)"
  )
endif()
set(LINK_FLAGS "${LINK_FLAGS} ${CMAKE_REQUIRED_LIBRARIES}")            #
#### </Boost.program_options> ##########################################

#### <OpenMP> (optional) ###############################################
if ("${OPENMP}" STREQUAL "NO")                                         #   
  message("-- Building without OpenMP (-DOPENMP=NO specified)")
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
endif()                                                                #
#### </OpenMP> #########################################################

#### <NetCDF> (optional) ###############################################
if("${NETCDF}" STREQUAL "NO") # handles -DNETCDF=NO                    #
  message("-- Building without netCDF (-DNETCDF=NO specified)")
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
        suitable Debian/Ubuntu package: libnetcdf-dev
        suitable Macports port: netcdf
        erase the CMakeCache.txt before re-running cmake
    ")
  else()
    set(USE_NETCDF ${Find_NetCDF})
    set(LINK_FLAGS "${LINK_FLAGS} ${CMAKE_REQUIRED_LIBRARIES}") 
  endif()
endif()                                                                #
#### </NetCDF> #########################################################

#### <Boost.Thread> (optional) #########################################
if ("${THREADS}" STREQUAL "NO") # handles -DTHREADS=NO                 #
  message("-- Building without Boost.Thread (-DTHREADS=NO specified)")
else()
  set(CMAKE_REQUIRED_LIBRARIES "-lboost_thread-mt")
  check_cxx_source_compiles(
    "#include <boost/thread.hpp>\n int main() { boost::thread th; }" 
    Find_Boost_Thread
  )
  if(NOT Find_Boost_Thread) 
    message(WARNING "
      Boost.Thread not found! (cf. CMakeFiles/CMakeError.log) 
      Threads parallelisation will not be available.
        suitable Debian/Ubuntu package: libboost-dev
        suitable Macports port: boost / boost+openmpi
        (after installation delete CMakeCache.txt before running cmake)
    ")
  else()
    set(USE_BOOST_THREAD ${Find_Boost_Thread})
    set(LINK_FLAGS "${LINK_FLAGS} ${CMAKE_REQUIRED_LIBRARIES}") 
  endif()
endif()                                                                #
#### </Boost.Thread> ###################################################

#### <Boost.MPI> (optional) ############################################
if ("${MPI}" STREQUAL "NO") # handles -DMPI=NO
  message("-- Building without Boost.MPI (-DMPI=NO specified)")
else()
  set(CMAKE_REQUIRED_LIBRARIES "-lboost_mpi-mt")
  check_cxx_source_compiles(
    "#include <boost/mpi.hpp>
     int main() {boost::mpi::communicator().send(1,0,0.);}"
    Find_Boost_MPI
  )
  if(NOT Find_Boost_MPI)
    message(WARNING "
      Boost.MPI not found! (cf. CMakeFiles/CMakeError.log)
      MPI parallelisation will not be available.
      Possible reasons:
      - a wrapper like mpic++ need to be used (e.d. export CXX=mpic++)
      - Boost.MPI needs to be installed:
        suitable Debian/Ubuntu package: ... TODO ...
        suitable Macports port: ... TODO ...
        (after installation delete CMakeCache.txt before running cmake)
    ")
  else()
    set(USE_BOOST_MPI ${Find_Boost_MPI})
    set(LINK_FLAGS "${LINK_FLAGS} ${CMAKE_REQUIRED_LIBRARIES}")
  endif()
endif()
#### </Boost.MPI> ######################################################

add_executable(icicles src/icicles.cpp)
set_target_properties(icicles PROPERTIES LINK_FLAGS ${LINK_FLAGS})
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/src/config.hpp)

find_package(Doxygen)
if (DOXYGEN_FOUND)
  #add_subdirectory(doc) (FIXME - cmake-controlled doc build)
else()
  message(WARNING "Doxygent not found - 'make doc' will not be available.")
endif()

enable_testing()

add_subdirectory(tests)
#add_subdirectory(src) # has to be done after the config.hpp file is created! (...?)
