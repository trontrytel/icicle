/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef CONFIG_H
#  define CONFIG_H

// compiler information
#define INFO_COMPILER "/opt/local/bin/clang++"
#define INFO_COMPILER_VERSION "clang version 2.9 (tags/RELEASE_29/final)"
#define INFO_COMPILER_CXX_FLAGS "-O0 -g -DBZ_DEBUG -I/opt/local/include  -pthread -pthread"
#define INFO_COMPILER_LINK_FLAGS " -lblitz -lboost_program_options-mt -lnetcdf_c++4 -lboost_thread-mt -lboost_timer-mt -lboost_system-mt"
#define INFO_COMPILER_USER "slayoo@eyrie.prac.igf"
#define INFO_COMPILER_SYSTEM "Darwin-10.6.0 i386"
#define INFO_COMPILER_TARGET "Darwin-10.6.0 i386"

// compile-time options
#define USE_NETCDF 1
#define USE_BOOST_THREAD 1
/* #undef USE_BOOST_MPI */
#define USE_BOOST_TIMER 1

#endif
