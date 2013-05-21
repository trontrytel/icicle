/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "mtx.hpp" 

#include <string>
using std::string;

template <typename real_t>
class out 
{
  public: virtual void record( 
    const string &name,
    const mtx::arr<real_t> &psi, // TODO: rename: that handles aux vars as well... 
    const mtx::idx &ijk, 
    const unsigned long t
  ) = 0;
  
  // ensuring derived-class destructor will be called
  public: virtual ~out() {}
};
