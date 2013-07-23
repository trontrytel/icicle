/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date December 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the @ref inf class - a facility for gathering simulation info
 */
#pragma once

#include "mtx.hpp"

#include <string>
using std::string;

#include <map>
using std::map;

/// @brief a facility for gathering simulation info - timing, compiler configurration, etc
class inf 
{
  private: string options;
  public: inf(const string &options)
    : options(options)
  { }

  public: map<string,string> get_map();
};
