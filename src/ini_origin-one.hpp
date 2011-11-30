/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef INI_ORIGIN_ONE_HPP
#  define INI_ORIGIN_ONE_HPP

#  include "ini.hpp" 

template <typename real_t>
class ini_origin_one : public ini<real_t>
{
  private: quantity<si::length, real_t> dxhlf, dyhlf, dzhlf;

  public: ini_origin_one(
    const quantity<si::length, real_t> &dx,
    const quantity<si::length, real_t> &dy,
    const quantity<si::length, real_t> &dz
  )
    : dxhlf(real_t(.5) * dx), dyhlf(real_t(.5) * dy), dzhlf(real_t(.5) * dz)
  {}

  public: quantity<si::dimensionless, real_t> psi(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) 
  {
    return (x == dxhlf && y == dyhlf && z == dzhlf) ? real_t(1) : real_t(0);
  }
};
#endif
