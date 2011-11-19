## @file
#  @author 
#    Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
#    Anna Jaruga (University of Warsaw) 
#  @date November 2011
#  @section LICENSE
#    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)

cmake_minimum_required(VERSION 2.6) # for simplified if-constructs
project(icicle CXX)
add_executable(icicle src/icicle.cpp)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES 
  "common.hpp;Makefile;CMakeCache.txt;cmake_install.cmake;CMakeFiles;Testing;CTestTestfile.cmake;tests/Makefile;tests/CMakeCache.txt;tests/cmake_install.cmake;tests/CMakeFiles;tests/Testing;tests/CTestTestfile.cmake;doc/Makefile;doc/CMakeCache.txt;doc/cmake_install.cmake;doc/CMakeFiles;doc/Testing;doc/CTestTestfile.cmake") 
include(CheckCXXSourceCompiles) 


#### Boost units (mandatory) 
check_cxx_source_compiles("#include <boost/units/systems/si.hpp>\n main() {}" Find_Boost_Units)
if(NOT Find_Boost_Units) 
  message(FATAL_ERROR "Boost.Units not found! (cf. CMakeFiles/CMakeError.log)")
endif()


#### Boost lexical cast (mandatory) 
check_cxx_source_compiles("#include <boost/lexical_cast.hpp>\n main() {}" Find_Boost_Lexical_Cast)
if(NOT Find_Boost_Lexical_Cast) 
  message(FATAL_ERROR "Boost.lexical_cast not found! (cf. CMakeFiles/CMakeError.log)")
endif()


#### Blitz (mandatory) 
check_cxx_source_compiles("#include <blitz/array.h>\n main() {}" Find_Blitz)
if(NOT Find_Blitz) 
  message(FATAL_ERROR "Blitz++ not found! (cf. CMakeFiles/CMakeError.log)")
endif()


#### Boost program options (mandatory) 
set(CMAKE_REQUIRED_LIBRARIES "-lboost_program_options")
check_cxx_source_compiles(
  "#include <boost/program_options.hpp>\n main() { boost::program_options::options_description d; }" 
  Find_Boost_Program_Options
)
if(NOT Find_Boost_Program_Options) 
  message(FATAL_ERROR "Boost.Program_Options not found! (cf. CMakeFiles/CMakeError.log)")
endif()
set(LINK_FLAGS "${LINK_FLAGS} -lboost_program_options") 


#### OpenMP (optional) 
set(TMP ${CMAKE_CXX_FLAGS})
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
check_cxx_source_compiles("#include <omp.h>\n main() { omp_set_num_threads(2); }" Find_OpenMP)
if(NOT Find_OpenMP)
  message(WARNING "OpenMP not found (cf. CMakeFiles/CMakeError.log) - OpenMP parallelisation will not be available.")
  set(CMAKE_CXX_FLAGS ${TMP})
endif()


#### NetCDF (optional)
#TODO handle -DNETCDF=no
set(TMP ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES "-lnetcdf_c++")
check_cxx_source_compiles(
  "#include <netcdfcpp.h>\n main() { NcFile *f; }" 
  Find_NetCDF
)
if(NOT Find_NetCDF) 
  message(WARNING "
    NetCDF not found! (cf. CMakeFiles/CMakeError.log) 
    NetCDF output will not be available.
      suitable Debian/Ubuntu package: libnetcdf-dev
      suitable Macports port: netcdf
      erase the CMakeCache.txt before re-running cmake
  ")
else()
  set(USE_NETCDF ${Find_NetCDF})
  set(LINK_FLAGS "${LINK_FLAGS} ${CMAKE_REQUIRED_LIBRARIES}") 
endif()
set(CMAKE_REQUIRED_LIBRARIES ${TMP})


#### Boost.Thread (optional)
#TODO handle -D...THREADS=no
set(TMP ${CMAKE_REQUIRED_LIBRARIES})
set(CMAKE_REQUIRED_LIBRARIES "-lboost_thread-mt")
check_cxx_source_compiles(
  "#include <boost/thread.hpp>\n main() { boost::thread th; }" 
  Find_Boost_Thread
)
if(NOT Find_Boost_Thread) 
  message(WARNING "
    Boost.Thread not found! (cf. CMakeFiles/CMakeError.log) 
    Threads parallelisation will not be available.
      suitable Debian/Ubuntu package: libboost-dev
      suitable Macports port: boost / boost+openmpi
      erase the CMakeCache.txt before re-running cmake
  ")
else()
  set(USE_BOOST_THREAD ${Find_Boost_Thread})
  set(LINK_FLAGS "${LINK_FLAGS} ${CMAKE_REQUIRED_LIBRARIES}") 
endif()
set(CMAKE_REQUIRED_LIBRARIES ${TMP})

set_target_properties(icicle PROPERTIES LINK_FLAGS ${LINK_FLAGS}) 
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.hpp.in ${CMAKE_CURRENT_BINARY_DIR}/src/config.hpp)

enable_testing()
add_subdirectory(tests)

string(COMPARE EQUAL ${CMAKE_BUILD_TYPE} "Debug" dbg)
if(dbg)
  find_package(Doxygen)
  if (DOXYGEN_FOUND)
    #add_subdirectory(doc) (FIXME - cmake-controlled doc build)
  else()
    message(WARNING "Doxygent not found - 'make doc' will not be available.")
  endif()

  add_subdirectory(src)
endif()
