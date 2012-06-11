/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the ini_func_cone class
 */
#include "cfg.hpp"
#include "ini_func_cone.hpp" 

template <typename real_t>
ini_func_cone<real_t>::ini_func_cone(
  const grd<real_t> &grid,
  const quantity<si::length, real_t> &h,
  const quantity<si::length, real_t> &x0,
  const quantity<si::length, real_t> &z0,
  const quantity<si::length, real_t> &r,
  const quantity<si::length, real_t> &h0
)
  : ini_func<real_t>(grid), h(h), x0(x0), z0(z0), r(r), h0(h0)
{}

template <typename real_t>
quantity<si::dimensionless, real_t> ini_func_cone<real_t>::psi(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &,
    const quantity<si::length, real_t> &z
) const
{
    // TODO: units
    real_t tmp = pow((x-x0) / si::metres, 2) + pow((z-z0) / si::metres, 2);
    tmp /= pow(r/h, 2);
    tmp = h/si::meters - sqrt(tmp);
    return (h0 / si::metres) + (tmp > 0 ? tmp : 0.);
}

// explicit instantiations
#if defined(USE_FLOAT)
template class ini_func_cone<float>;
#endif
#if defined(USE_DOUBLE)
template class ini_func_cone<double>;
#endif
#if defined(USE_LDOUBLE)
template class ini_func_cone<long double>;
#endif
