/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
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
using std::auto_ptr;
using std::cout;
using std::cerr;
using std::endl;
using std::exception;
using std::cos;
using std::sin;

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

typedef boost::units::multiply_typeof_helper<
    si::velocity,
    si::length
  >::type velocity_times_length;

// Boost.ptr_vector // TODO: should it be checked at 'configure' step?
#  include <boost/ptr_container/ptr_vector.hpp>
using boost::ptr_vector;

#endif
