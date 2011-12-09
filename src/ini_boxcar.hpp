/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef INI_BOXCAR_HPP
#  define INI_BOXCAR_HPP

#  include "ini.hpp" 

template <typename real_t>
class ini_boxcar : public ini<real_t>
{
  private: quantity<si::length, real_t> a, b;
  private: quantity<si::dimensionless, real_t> A, A0;

  public: ini_boxcar(
    const quantity<si::length, real_t> &a,
    const quantity<si::length, real_t> &b,
    const quantity<si::dimensionless, real_t> &A,
    const quantity<si::dimensionless, real_t> &A0
  )
    : a(a), b(b), A(A), A0(A0)
  {}

  public: quantity<si::dimensionless, real_t> psi(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) 
  {
    if (
      x < a || x > b || 
      y < a || y > b || 
      z < a || z > b
    ) return A0;
    return A0 + A;
  }
};
#endif
