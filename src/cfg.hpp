/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief compile-time configuration (CMake-generated)
 */
#ifndef CONFIG_H
#  define CONFIG_H

// compiler information
#define INFO_COMPILER "/usr/lib/gcc-snapshot/bin/g++"
#define INFO_COMPILER_VERSION "g++ (Ubuntu/Linaro 20120625-0ubuntu1) 4.8.0 20120625 (experimental) [trunk revision 188931]"
#define INFO_COMPILER_CXX_FLAGS "-std=c++0x -DNDEBUG -Wfatal-errors -Ofast  -fopenmp"
#define INFO_COMPILER_LINK_FLAGS ""
#define INFO_COMPILER_USER "ania@laptok"
#define INFO_COMPILER_SYSTEM "Linux-3.2.0-27-generic x86_64"
#define INFO_COMPILER_TARGET "Linux-3.2.0-27-generic x86_64"

// TODO: why not to use the same syntax above and below???

// compile-time options
#define USE_FLOAT 1
#define USE_DOUBLE 1
#define USE_LDOUBLE 1
/* #undef USE_FLOAT128 */

#define USE_NETCDF 1
/* #undef USE_THRUST */
/* #undef USE_CUDA */
/* #undef USE_BOOST_THREAD */
/* #undef USE_BOOST_MPI */
#define USE_BOOST_TIMER 1
#define USE_BOOST_ODEINT 1

#endif
