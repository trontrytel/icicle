/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date January 2012
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include "cfg.hpp"
#include "ini_func_gauss.hpp" 

template <typename real_t>
ini_func_gauss<real_t>::ini_func_gauss(
  const grd<real_t> &grid,
  const quantity<si::dimensionless, real_t> &A,
  const quantity<si::dimensionless, real_t> &A0,
  const quantity<si::length, real_t> &x0,
  const quantity<si::length, real_t> &y0,
  const quantity<si::length, real_t> &z0,
  const quantity<si::length, real_t> &sx,
  const quantity<si::length, real_t> &sy,
  const quantity<si::length, real_t> &sz
)
  : ini_func<real_t>(grid), A(A), A0(A0), x0(x0), y0(y0), z0(z0), sx(sx), sy(sy), sz(sz)
{}

template <typename real_t>
quantity<si::dimensionless, real_t> ini_func_gauss<real_t>::psi(
  const quantity<si::length, real_t> &x,
  const quantity<si::length, real_t> &y,
  const quantity<si::length, real_t> &z
) const 
{
  return A0 + A * exp(
    real_t(-.5) * (
      pow<2>(x-x0) * pow<-2>(sx) +
      pow<2>(y-y0) * pow<-2>(sy) +
      pow<2>(z-z0) * pow<-2>(sz)
    )
  );
}

// explicit instantiations
#if defined(USE_FLOAT)
template class ini_func_gauss<float>;
#endif
#if defined(USE_DOUBLE)
template class ini_func_gauss<double>;
#endif
#if defined(USE_LDOUBLE)
template class ini_func_gauss<long double>;
#endif
