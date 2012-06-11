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
#pragma once
#include "vel_func_stream.hpp"

template <typename real_t>
class vel_func_stream_rasinski : public vel_func_stream<real_t>
{
  private: const quantity<si::dimensionless, real_t> pi;
  private: const quantity<si::length, real_t> X, Z;
  private: const quantity<divide_typeof_helper<si::momentum, si::area>::type, real_t> A;

  public: vel_func_stream_rasinski(
    const grd<real_t> &grid,
    const string file,
    quantity<divide_typeof_helper<si::momentum, si::area>::type, real_t> A
  );
 
  public: quantity<si::velocity, real_t> u(
    const quantity<si::length, real_t> &x, 
    const quantity<si::length, real_t> &z, 
    const quantity<si::length, real_t> &
  ) const;

  public: quantity<si::velocity, real_t> v(
    const quantity<si::length, real_t> &x, 
    const quantity<si::length, real_t> &z, 
    const quantity<si::length, real_t> &
  ) const;
};
