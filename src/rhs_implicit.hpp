/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef RHS_IMPLICIT_HPP
#  define RHS_IMPLICIT_HPP

#  include "rhs.hpp" 

template <typename real_t>
class rhs_implicit : public rhs<real_t>
{
  public: virtual void explicit_part(
    mtx::arr<real_t> &R, 
    const ptr_map<string, mtx::arr<real_t>> &aux,
    const mtx::arr<real_t> * const * const psi, 
    const quantity<si::time, real_t> t
  )
  {}
};  
#endif
