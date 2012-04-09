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
 *      \frac{\pi}{2}\frac{z}{Z_{clb}}
 *    \right) cos\left(
 *      2\pi\frac{x}{X}
 *    \right)    for z \le Z_{clb}
 *
 *    \psi(x,z) = -A sin\left(
 *      \frac{\pi}{2}\left(\frac{z-Z_{top}}{Z_{top}-Z_{clb}}+1\right)
 *    \right) cos\left(
 *      2\pi \frac{x}{X}
 *    \right)    for z \gt Z_{clb}
 *
 *    (eq. 2 in Rasinski et al. 2011, Atmos. Res. 102)
 */
#ifndef VEL_FUNC_STREAM_RASINSKI_HPP
#  define VEL_FUNC_STREAM_RASINSKI_HPP

#  include "vel_func_stream.hpp"

template <typename real_t>
class vel_func_stream_rasinski : public vel_func_stream<real_t>
{
  private: quantity<si::dimensionless, real_t> pi;
  private: quantity<si::length, real_t> X, Z_clb, Z_top;
  private: quantity<divide_typeof_helper<si::momentum, si::area>::type, real_t> A;

  public: vel_func_stream_rasinski(
    const grd<real_t> &grid,
    const string file,
    quantity<si::length, real_t> Z_clb,
    quantity<si::length, real_t> Z_top, 
    quantity<divide_typeof_helper<si::momentum, si::area>::type, real_t> A
  ) :
    vel_func_stream<real_t>(grid, file),
    pi(real_t(4)*atan(real_t(1))), 
    X(real_t(grid.nx()) * grid.dx()), 
    Z_clb(Z_clb), 
    Z_top(Z_top), 
    A(A)
  {}
 
  public: virtual quantity<si::velocity, real_t> u(
    const quantity<si::length, real_t> &x, 
    const quantity<si::length, real_t> &z, 
    const quantity<si::length, real_t> &
  )   
  {   
    // rho * u = - \frac{\partial \psi}{\partial z}
    return z <= Z_clb
      ?  pi/2 * (A /      Z_clb      / this->rho(z)) * real_t(cos(2*pi * x / X) * cos(pi/2 *    z        /      Z_clb     ))  
      : -pi/2 * (A / (Z_clb - Z_top) / this->rho(z)) * real_t(cos(2*pi * x / X) * sin(pi/2 * (z - Z_clb) / (Z_clb - Z_top)));
  }

  public: virtual quantity<si::velocity, real_t> v(
    const quantity<si::length, real_t> &x, 
    const quantity<si::length, real_t> &z, 
    const quantity<si::length, real_t> &
  )   
  {   
    // rho * w = \frac{\partial \psi}{\partial x}
    return z <= Z_clb 
      ? pi/2 * (A / X / this->rho(z)) * real_t(sin(2*pi * x / X) * sin(pi/2 *    z        /      Z_clb     ))  
      : pi/2 * (A / X / this->rho(z)) * real_t(sin(2*pi * x / X) * cos(pi/2 * (z - Z_clb) / (Z_clb - Z_top)));
  }
};
#endif
