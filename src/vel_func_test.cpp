/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include "vel_func_test.hpp"

template <typename real_t>
vel_func_test<real_t>::vel_func_test(
    const grd<real_t> &grid,
    quantity<si::frequency, real_t> omega,
    quantity<si::length, real_t> x0,
    quantity<si::length, real_t> z0,
    quantity<si::velocity, real_t> vv
)
: vel_func<real_t>(grid), omega(omega), x0(x0), z0(z0), vv(vv)
{ }
 
template <typename real_t>
quantity<si::velocity, real_t> vel_func_test<real_t>::u(
  const quantity<si::length, real_t> &, 
  const quantity<si::length, real_t> &, 
  const quantity<si::length, real_t> &z
) const
{   
  return -omega*(z-z0);
}

template <typename real_t>
quantity<si::velocity, real_t> vel_func_test<real_t>::w(
  const quantity<si::length, real_t> &x, 
  const quantity<si::length, real_t> &, 
  const quantity<si::length, real_t> &
) const
{   
  return omega*(x-x0);
}

template <typename real_t>
quantity<si::velocity, real_t> vel_func_test<real_t>::v(
  const quantity<si::length, real_t> &,
  const quantity<si::length, real_t> &,
  const quantity<si::length, real_t> &
) const
{
  return vv;
}

// explicit instantiations
#include "cfg/cfg_types.hpp"
#if defined(USE_FLOAT)
template class vel_func_test<float>;
#endif
#if defined(USE_DOUBLE)
template class vel_func_test<double>;
#endif
#if defined(USE_LDOUBLE)
template class vel_func_test<long double>;
#endif
