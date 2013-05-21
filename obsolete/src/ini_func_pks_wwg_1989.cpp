/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include "ini_func_pks_wwg_1989.hpp" 

template <typename real_t>
real_t ini_func_pks_wwg_1989<real_t>::psi_0(real_t x) const
{
  if (x >= 8 && x <= 28) return -1.;
  if (x > 28 && x <= 39) return 1;
  return 0.;
}

template <typename real_t>
quantity<si::dimensionless, real_t> ini_func_pks_wwg_1989<real_t>::psi(
  const quantity<si::length, real_t> &x,
  const quantity<si::length, real_t> &,
  const quantity<si::length, real_t> &
) const
{
  return 2 + psi_0(x / si::metres) 
    * (1 + .3 * sin(2 * M_PI /  9. * (x/dx))) 
    * (1 + .4 * sin(2 * M_PI / 10. * (x/dx)));
}

// explicit instantiations
#define ICICLE_INSTANTIATE_CLASS ini_func_pks_wwg_1989
#include "cmn/cmn_instant.hpp"
