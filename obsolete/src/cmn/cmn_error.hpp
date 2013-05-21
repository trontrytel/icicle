/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - February 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief error reporting macros
 */

#pragma once

#include <exception>
#include <iostream>

#define error_macro(msg) \
{ \
  std::cerr << "-- error: " << msg << std::endl; \
  throw std::exception(); \
}

#define warning_macro(msg) \
{ \
  std::cerr << "-- warning: " << msg << std::endl; \
}

#define notice_macro(msg) \
{ \
  std::cerr << "-- notice: " << msg << std::endl; \
}
