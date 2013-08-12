#pragma once

#include <iostream>
#include <set>
#include <string>
#include <sstream>
using std::set;
using std::string;
using std::ostringstream;
using std::cerr;
using std::endl;
using std::exception;

#define error_macro(msg) \
{ \
  cerr << "error: " << msg << endl; \
  throw exception(); \
}

#define notice_macro(msg) \
{ \
  cerr << " info: " << msg << endl; \
}
