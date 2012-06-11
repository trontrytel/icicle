/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#include "cfg.hpp"
#include "vel_func_stream_rasinski.hpp"

template <typename real_t>
vel_func_stream_rasinski<real_t>::vel_func_stream_rasinski(
  const grd<real_t> &grid,
  const string file,
  quantity<divide_typeof_helper<si::momentum, si::area>::type, real_t> A
) :
  vel_func_stream<real_t>(grid, file),
  pi(real_t(4)*atan(real_t(1))), 
  X(real_t(grid.nx()) * grid.dx()), 
  Z(real_t(grid.ny()) * grid.dy()), 
  A(A)
{}
 
template <typename real_t>
quantity<si::velocity, real_t> vel_func_stream_rasinski<real_t>::u(
  const quantity<si::length, real_t> &x, 
  const quantity<si::length, real_t> &z, 
  const quantity<si::length, real_t> &
) const 
{   
  // rho * u = - \frac{\partial \psi}{\partial z}
  return pi * (A / Z / this->rho(z)) * real_t(cos(2*pi * x / X) * cos(pi * z / Z));
}

template <typename real_t>
quantity<si::velocity, real_t> vel_func_stream_rasinski<real_t>::v(
  const quantity<si::length, real_t> &x, 
  const quantity<si::length, real_t> &z, 
  const quantity<si::length, real_t> &
) const
{   
  // rho * w = \frac{\partial \psi}{\partial x}
  return 2 * pi * (A / X / this->rho(z)) * real_t(sin(2*pi * x / X) * sin(pi * z / Z)); 
}

// explicit instantiations
#if defined(USE_FLOAT)
template class vel_func_stream_rasinski<float>;
#endif
#if defined(USE_DOUBLE)
template class vel_func_stream_rasinski<double>;
#endif
#if defined(USE_LDOUBLE)
template class vel_func_stream_rasinski<long double>;
#endif
