/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef INI_FUNC_CONE_HPP
#  define INI_FUNC_CONE_HPP

#  include "ini_func.hpp" 

template <typename real_t>
class ini_func_cone : public ini_func<real_t>
{
  private: quantity<si::length, real_t> h, x0, z0, r, h0;

  public: ini_func_cone(
    const grd<real_t> &grid,
    const quantity<si::length, real_t> &h,
    const quantity<si::length, real_t> &x0,
    const quantity<si::length, real_t> &z0,
    const quantity<si::length, real_t> &r,
    const quantity<si::length, real_t> &h0
  )
    : ini_func<real_t>(grid), h(h), x0(x0), z0(z0), r(r), h0(h0)
  {}

  public: quantity<si::dimensionless, real_t> psi(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &,
    const quantity<si::length, real_t> &z
  ) 
  {
    // TODO: units
    real_t tmp = pow((x-x0) / si::metres, 2) + pow((z-z0) / si::metres, 2);
    tmp /= pow(r/h, 2);
    tmp = h/si::meters - sqrt(tmp);
    return (h0 / si::metres) + (tmp > 0 ? tmp : 0.);
  }
};
#endif