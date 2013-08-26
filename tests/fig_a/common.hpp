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


// bin sizes for calc and plot
vector<quantity<si::length>> bins_dry()
{
  vector<quantity<si::length>> ret;
  // dry radius bins: .001 ... .01 ... .1 ... 1 (30 bins in total)
  for (int i = 0; i < 30; ++i)
    ret.push_back(1e-6 * pow(10, -3 + i * .1) * si::metres);
  return ret;
}

vector<quantity<si::length>> bins_wet()
{
    vector<quantity<si::length>> ret;
    for (int i = 0; i < 50; ++i)
      ret.push_back(1e-6 * pow(10, -3 + i * .2) * si::metres);
    return ret;
  }
