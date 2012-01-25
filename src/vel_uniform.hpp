/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    uniform velocity field
 */
#ifndef VEL_UNIFORM_HPP
#  define VEL_UNIFORM_HPP

#  include "cmn.hpp" // root class, error reporting
#  include "vel.hpp"

template <typename real_t>
class vel_uniform : public vel<real_t>
{
  public: bool is_constant() { return true; }

  private: quantity<si::velocity, real_t> uu, vv, ww;

  public: vel_uniform(
    quantity<si::velocity, real_t> uu,
    quantity<si::velocity, real_t> vv,
    quantity<si::velocity, real_t> ww
  )
    : uu(uu), vv(vv), ww(ww)
  {}

  public: virtual quantity<si::velocity, real_t> u(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) { return uu; }
  public: virtual quantity<si::velocity, real_t> v(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) { return vv; }
  public: virtual quantity<si::velocity, real_t> w(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) { return ww; }
};
#endif
