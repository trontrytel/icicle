/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date Februaryy 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef RHS_HPP
#  define RHS_HPP

#  include "cmn.hpp" 
#  include "mtx.hpp"

template <typename real_t>
class rhs 
{
  public: virtual void operator()(
    mtx::arr<real_t> &R, 
    mtx::arr<real_t> **psi, 
    mtx::idx &ijk
  ) { assert(false); };
};  

#endif
