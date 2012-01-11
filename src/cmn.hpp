/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */

#ifndef CMN_HPP
#  define CMN_HPP

// most common includes
#  include <boost/units/systems/si.hpp>
#  include <boost/units/io.hpp>
using namespace boost::units;

// TODO: ...
using namespace std;

// overloading the default d-tor with a virtual one (enforces execution of child d-tors)
class root { public: virtual ~root() {} };

// error reporting
#  define error_macro(msg) \
{ \
  std::cerr << "-- error: " << msg << std::endl; \
  throw std::exception(); \
}
#  define warning_macro(msg) \
{ \
  std::cerr << "-- warning: " << msg << std::endl; \
}

// some non-standard units
typedef multiply_typeof_helper<
    si::velocity,
    si::length
  >::type velocity_times_length;

#endif
