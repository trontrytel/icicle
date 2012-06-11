/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#pragma once
#include "vel_func.hpp"

template <typename real_t>
class vel_func_test : public vel_func<real_t>
{
  private: const quantity<si::frequency, real_t> omega;
  private: const quantity<si::length, real_t> x0,z0;
  private: const quantity<si::velocity, real_t> vv;

  public: vel_func_test(
    const grd<real_t> &grid,
    quantity<si::frequency, real_t> omega,
    quantity<si::length, real_t> x0,
    quantity<si::length, real_t> z0,
    quantity<si::velocity, real_t> vv
  );
 
  public: quantity<si::velocity, real_t> u(
    const quantity<si::length, real_t> &, 
    const quantity<si::length, real_t> &, 
    const quantity<si::length, real_t> &z
  ) const;

  public: quantity<si::velocity, real_t> w(
    const quantity<si::length, real_t> &x, 
    const quantity<si::length, real_t> &, 
    const quantity<si::length, real_t> &
  ) const;

  public: quantity<si::velocity, real_t> v(
    const quantity<si::length, real_t> &,
    const quantity<si::length, real_t> &,
    const quantity<si::length, real_t> &
  ) const;
};
