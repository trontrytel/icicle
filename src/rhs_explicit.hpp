/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef RHS_EXPLICIT_HPP
#  define RHS_EXPLICIT_HPP

#  include "rhs.hpp" 

template <typename real_t>
class rhs_explicit : public rhs<real_t>
{
  public: real_t implicit_part(
    const quantity<si::time, real_t> dt 
  ) const
  { 
    return real_t(0); 
  }
};  
#endif
