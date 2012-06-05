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
#ifndef INI_FUNC_PKS_WWG_1989_HPP
#  define INI_FUNC_PKS_WWG_1989_HPP

#  include "ini_func.hpp" 
#  include "grd.hpp"

template <typename real_t>
class ini_func_pks_wwg_1989 : public ini_func<real_t>
{
  private: quantity<si::length, real_t> dx;

  public: ini_func_pks_wwg_1989(const grd<real_t> &grid)
    : ini_func<real_t>(grid), dx(grid.dx())
  {
  }

  private: real_t psi_0(real_t x) const
  {
    if (x >= 8 && x <= 28) return -1.;
    if (x > 28 && x <= 39) return 1;
    return 0.;
  }

  public: quantity<si::dimensionless, real_t> psi(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &,
    const quantity<si::length, real_t> &
  ) const
  {
    return 2 + psi_0(x / si::metres) 
      * (1 + .3 * sin(2 * M_PI /  9. * (x/dx))) 
      * (1 + .4 * sin(2 * M_PI / 10. * (x/dx)));
  }
};
#endif
