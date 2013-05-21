/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    1D initial condition from Smolarkiewicz & Grabowski 1989
 *    (section 4, eq. 22)
 */
#pragma once
#include "ini_func.hpp" 
#include "grd.hpp"

template <typename real_t>
class ini_func_pks_wwg_1989 : public ini_func<real_t>
{
  private: quantity<si::length, real_t> dx;

  public: ini_func_pks_wwg_1989(const grd<real_t> &grid)
    : ini_func<real_t>(grid), dx(grid.dx())
  {}

  private: real_t psi_0(real_t x) const;

  public: quantity<si::dimensionless, real_t> psi(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &,
    const quantity<si::length, real_t> &
  ) const;
};
