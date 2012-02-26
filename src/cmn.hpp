/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief common includes, using statements, macro definitions etc
 */

#ifndef CMN_HPP
#  define CMN_HPP

// STL includes
#  include <cmath>
#  include <iostream>
#  include <string>
#  include <vector>
#  include <map>
#  include <memory>
using std::string;
using std::vector;
using std::ostringstream;
using std::map;
using std::unique_ptr;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::cos;
using std::sin;
using std::isfinite;
using std::pair;

// the Lynton Appel's netCDF-4 C++ API (since netCDF 4.1.1)
#    ifdef USE_BOOST_MPI 
// fixes preprocessor macro redefinition conflict with MPI
// cf. http://www.unidata.ucar.edu/mailing_lists/archives/netcdfgroup/2009/msg00350.html
#      include "mpi/mpi.h"
#      define MPI_INCLUDED
#    endif    
#    include <netcdf>
using netCDF::NcFile;
using netCDF::NcVar;
using netCDF::NcDim;
using netCDF::ncFloat;
using netCDF::exceptions::NcException;

// overloading the default d-tor with a virtual one (enforces execution of child d-tors)
class root { public: virtual ~root() {} };

// error reporting
#  define error_macro(msg) \
{ \
  cerr << "-- error: " << msg << endl; \
  throw exception(); \
}
#  define warning_macro(msg) \
{ \
  cerr << "-- warning: " << msg << endl; \
}

// Boost.Units
#  include <boost/units/systems/si.hpp>
#  include <boost/units/cmath.hpp> 
#  include <boost/units/io.hpp>
using boost::units::quantity; 
namespace si = boost::units::si;
using boost::units::pow;

typedef boost::units::multiply_typeof_helper<
    si::velocity,
    si::length
  >::type velocity_times_length;

// Boost.ptr_vector // TODO: should it be checked at 'configure' step?
#  include <boost/ptr_container/ptr_vector.hpp>
using boost::ptr_vector;
using boost::nullable;

#endif
