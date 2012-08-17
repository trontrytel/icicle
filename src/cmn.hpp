/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief common includes, using statements, macro definitions etc
 */

#pragma once

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

// error reporting
#  include <exception>
using std::exception;
#  define error_macro(msg) \
{ \
  std::cerr << "-- error: " << msg << std::endl; \
  throw exception(); \
}
#  define warning_macro(msg) \
{ \
  std::cerr << "-- warning: " << msg << std::endl; \
}
#  define notice_macro(msg) \
{ \
  std::cerr << "-- notice: " << msg << std::endl; \
}

// Boost.Units
#  include <boost/units/systems/si.hpp>
#  include <boost/units/cmath.hpp> 
#  include <boost/units/io.hpp>
namespace si = boost::units::si;
using boost::units::quantity; 
using boost::units::pow;
using boost::units::multiply_typeof_helper;
using boost::units::divide_typeof_helper;
using boost::units::power_typeof_helper;
using boost::units::static_rational;

typedef multiply_typeof_helper< // TODO: get rid of it
    si::velocity,
    si::pressure
  >::type velocity_times_pressure;

typedef multiply_typeof_helper< // TODO: get rid of it
    si::acceleration,
    si::length
  >::type specific_energy;

// Boost.ptr_vector // TODO: should it be checked at 'configure' step?
#  include <boost/ptr_container/ptr_vector.hpp>
using boost::ptr_vector;
using boost::nullable;
#  include <boost/ptr_container/ptr_unordered_map.hpp>
using boost::ptr_unordered_map;
#  include <boost/assign/ptr_map_inserter.hpp>
using boost::assign::ptr_map_insert;
