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
  private: quantity<si::length, real_t> ax, bx, ay, by, az, bz;
  private: quantity<si::dimensionless, real_t> A, A0;

  public: ini_boxcar(
    const quantity<si::length, real_t> &ax,
    const quantity<si::length, real_t> &bx,
    const quantity<si::length, real_t> &ay,
    const quantity<si::length, real_t> &by,
    const quantity<si::length, real_t> &az,
    const quantity<si::length, real_t> &bz,
    const quantity<si::dimensionless, real_t> &A,
    const quantity<si::dimensionless, real_t> &A0
  )
    : ax(ax), bx(bx), ay(ay), by(by), az(az), bz(bz), A(A), A0(A0)
  {}

  public: quantity<si::dimensionless, real_t> psi(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) 
  {
    if (
      x < ax || x > bx || 
      y < ay || y > by || 
      z < az || z > bz
    ) return A0;
    return A0 + A;
  }
};
#endif
