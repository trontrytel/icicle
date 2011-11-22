/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
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
#ifndef VEL_2D_XZ_RASINSKI_HPP
#  define VEL_2D_XZ_RASINSKI_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting
#  include "vel_2d-xz.hpp"

template <typename real_t>
class vel_2d_xz_rasinski : public vel_2d_xz<real_t>
{
  private: quantity<si::dimensionless, real_t> pi;
  private: quantity<si::length, real_t> X, Z_clb, Z_top;
  private: quantity<velocity_times_length, real_t> A;

  public: vel_2d_xz_rasinski(
    quantity<si::length, real_t> X, quantity<si::length, real_t> Z_clb,
    quantity<si::length, real_t> Z_top, quantity<velocity_times_length, real_t> A
  )
    : pi(real_t(4)*atan(real_t(1))), X(X), Z_clb(Z_clb), Z_top(Z_top), A(A)
  {}
 
  public: virtual quantity<si::velocity> u(
    const quantity<si::length, real_t> &x, 
    const quantity<si::length, real_t> &, 
    const quantity<si::length, real_t> &z
  )   
  {   
    // u = - \frac{\partial \psi}{\partial z}
    // TODO: in the paper they divide it by the anelastic density
    return z <= Z_clb
      ?  pi/2 * A /      Z_clb      * cos(2*pi * x / X) * cos(pi/2 *    z        /      Z_clb     )
      : -pi/2 * A / (Z_clb - Z_top) * cos(2*pi * x / X) * sin(pi/2 * (z - Z_clb) / (Z_clb - Z_top)); 
  }

  public: virtual quantity<si::velocity> w(
    const quantity<si::length, real_t> &x, 
    const quantity<si::length, real_t> &, 
    const quantity<si::length, real_t> &z
  )   
  {   
    // w = \frac{\partial \psi}{\partial x}
    // TODO: in the paper they divide it by the anelastic density
    return z <= Z_clb 
      ? pi/2 * A / X * sin(2*pi * x / X) * sin(pi/2 *    z        /      Z_clb     )
      : pi/2 * A / X * sin(2*pi * x / X) * cos(pi/2 * (z - Z_clb) / (Z_clb - Z_top));
  }
};
#endif
