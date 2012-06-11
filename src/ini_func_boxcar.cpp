/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011 - January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the ini_func_boxcar class
 */
#include "cfg.hpp"
#include "ini_func_boxcar.hpp" 

template <typename real_t>
ini_func_boxcar<real_t>::ini_func_boxcar(
  const grd<real_t> &grid,
  const quantity<si::length, real_t> &ax,
  const quantity<si::length, real_t> &bx,
  const quantity<si::length, real_t> &ay,
  const quantity<si::length, real_t> &by,
  const quantity<si::length, real_t> &az,
  const quantity<si::length, real_t> &bz,
  const quantity<si::dimensionless, real_t> &A,
  const quantity<si::dimensionless, real_t> &A0
)
  : ini_func<real_t>(grid), ax(ax), bx(bx), ay(ay), by(by), az(az), bz(bz), A(A), A0(A0)
{}

template <typename real_t>
quantity<si::dimensionless, real_t> ini_func_boxcar<real_t>::psi(
  const quantity<si::length, real_t> &x,
  const quantity<si::length, real_t> &y,
  const quantity<si::length, real_t> &z
) const
{
  if (
    x < ax || x > bx || 
    y < ay || y > by || 
    z < az || z > bz
  ) return A0;
  return A0 + A;
}

// explicit instantiations
#if defined(USE_FLOAT)
template class ini_func_boxcar<float>;
#endif
#if defined(USE_DOUBLE)
template class ini_func_boxcar<double>;
#endif
#if defined(USE_LDOUBLE)
template class ini_func_boxcar<long double>;
#endif