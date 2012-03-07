/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date March 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef RHS_LINEAR_HPP
#  define RHS_LINEAR_HPP

#  include "rhs.hpp" 

template <typename real_t>
class rhs_linear : rhs<real_t>
{
  public: int eqid;
  public: real_t C;
  public: rhs_linear(int eqid, real_t C)
    : eqid(eqid), C(C)
  {}
};  
#endif
