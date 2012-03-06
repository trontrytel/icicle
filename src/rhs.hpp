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

/// cf. the discussion of implicit (linear) and explicit terms of the rhs in
/// section 2.2 of Prusa, Smolarkiewicz & Wyszogrodzki 2008 (Computers & Fluids)
template <typename real_t>
class rhs : root
{
  // TODO: stencil-extent-like method?
  public: virtual void explicit_part(
    mtx::arr<real_t> &R, 
    mtx::arr<real_t> **psi, 
    mtx::idx &ijk
  ) 
  {}

  public: virtual void implicit_part()
  {}
};  

#endif
