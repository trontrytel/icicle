/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "out.hpp"

template <typename real_t>
class out_debug : public out<real_t>
{
  public: void record(
    const string &name,
    const mtx::arr<real_t> &psi,
    const mtx::idx &ijk, 
    const unsigned long t
  );
};
