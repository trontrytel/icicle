/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    velocity field representing an idealized eddy circulation prescribed
 *    using the streamfunction \psi in the form:
 *
 *    \psi(x,z) = -A sin\left(
 *      \pi \frac{z}{Z}
 *    \right) cos\left(
 *      2\pi\frac{x}{X}
 *    \right) 
 *
 *    (similar to eq. 2 in Rasinski et al. 2011, Atmos. Res. 102)
 */
#ifndef VEL_FUNC_STREAM_RASINSKI_HPP
#  define VEL_FUNC_STREAM_RASINSKI_HPP

#  include "vel_func_stream.hpp"

template <typename real_t>
class vel_func_stream_rasinski : public vel_func_stream<real_t>
{
  private: quantity<si::dimensionless, real_t> pi;
  private: quantity<si::length, real_t> X, Z;
  private: quantity<divide_typeof_helper<si::momentum, si::area>::type, real_t> A;

  public: vel_func_stream_rasinski(
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
 
  public: virtual quantity<si::velocity, real_t> u(
    const quantity<si::length, real_t> &x, 
    const quantity<si::length, real_t> &z, 
    const quantity<si::length, real_t> &
  )   
  {   
    // rho * u = - \frac{\partial \psi}{\partial z}
    return pi * (A / Z / this->rho(z)) * real_t(cos(2*pi * x / X) * cos(pi * z / Z));
  }

  public: virtual quantity<si::velocity, real_t> v(
    const quantity<si::length, real_t> &x, 
    const quantity<si::length, real_t> &z, 
    const quantity<si::length, real_t> &
  )   
  {   
    // rho * w = \frac{\partial \psi}{\partial x}
    return 2 * pi * (A / X / this->rho(z)) * real_t(sin(2*pi * x / X) * sin(pi * z / Z)); 
  }
};
#endif
