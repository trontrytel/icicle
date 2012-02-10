/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef INI_FUNC_GAUSS_HPP
#  define INI_FUNC_GAUSS_HPP

#  include "ini_func.hpp" 

template <typename real_t>
class ini_func_gauss : public ini_func<real_t>
{
  private: quantity<si::dimensionless, real_t> A;
  private: quantity<si::length, real_t> x0, y0, z0, sx, sy, sz;

  public: ini_func_gauss(
    const grd<real_t> &grid,
    const quantity<si::dimensionless, real_t> &A,
    const quantity<si::length, real_t> &x0,
    const quantity<si::length, real_t> &y0,
    const quantity<si::length, real_t> &z0,
    const quantity<si::length, real_t> &sx,
    const quantity<si::length, real_t> &sy,
    const quantity<si::length, real_t> &sz
  )
    : ini_func<real_t>(grid), A(A), x0(x0), y0(y0), z0(z0), sx(sx), sy(sy), sz(sz)
  {}

  public: quantity<si::dimensionless, real_t> psi(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) 
  {
    return A * exp(
      real_t(-.5) * (
        pow<2>(x-x0) * pow<-2>(sx) +
        pow<2>(y-y0) * pow<-2>(sy) +
        pow<2>(z-z0) * pow<-2>(sz)
      )
    );
  }
};
#endif
