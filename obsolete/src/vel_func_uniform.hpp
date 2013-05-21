/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    uniform velocity field
 */
#pragma once
#include "vel_func.hpp"

template <typename real_t>
class vel_func_uniform : public vel_func<real_t>
{
  private: const quantity<si::velocity, real_t> uu, vv, ww;

  public: vel_func_uniform(
    const grd<real_t> &grid,
    const quantity<si::velocity, real_t> uu,
    const quantity<si::velocity, real_t> vv,
    const quantity<si::velocity, real_t> ww
  )
    : vel_func<real_t>(grid), uu(uu), vv(vv), ww(ww)
  {}

  public: quantity<si::velocity, real_t> u(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) const { return uu; }
  public: quantity<si::velocity, real_t> v(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) const { return vv; }
  public: quantity<si::velocity, real_t> w(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) const { return ww; }
};
