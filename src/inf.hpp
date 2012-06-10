/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date December 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the @ref inf class - a facility for gathering simulation info
 */
#pragma once

#include "cmn.hpp"
#include "mtx.hpp"

#ifdef USE_BOOST_TIMER
#  include <boost/timer/timer.hpp>
#endif

/// @brief a facility for gathering simulation info - timing, compiler configurration, etc
class inf 
{
#ifdef USE_BOOST_TIMER
  private: boost::timer::cpu_timer tmr;
#endif
  private: string options;
  public: inf(const string &options)
    : options(options)
  { }

  public: map<string,string> get_map();
};
