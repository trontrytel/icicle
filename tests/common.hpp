#pragma once

#include <iostream>
#include <set>
#include <string>
#include <sstream>
#include <vector>

using std::set;
using std::string;
using std::ostringstream;
using std::cerr;
using std::endl;
using std::exception;
using std::vector;

#include <boost/units/systems/si.hpp>
namespace si = boost::units::si;
using boost::units::quantity;


// error reporting
#define error_macro(msg) \
{ \
  cerr << "error: " << msg << endl; \
  throw exception(); \
}

#define notice_macro(msg) \
{ \
  cerr << " info: " << msg << endl; \
}
