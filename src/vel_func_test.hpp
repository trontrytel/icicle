/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef VEL_FUNC_TEST_HPP
#  define VEL_FUNC_TEST_HPP

#  include "vel_func.hpp"

template <typename real_t>
class vel_func_test : public vel_func<real_t>
{
  private: quantity<si::frequency, real_t> omega;
  private: quantity<si::length, real_t> x0,z0;
  private: quantity<si::velocity, real_t> vv;

  public: vel_func_test(
    grd<real_t> *grid,
    quantity<si::frequency, real_t> omega,
    quantity<si::length, real_t> x0,
    quantity<si::length, real_t> z0,
    quantity<si::velocity, real_t> vv
  )
    : vel_func<real_t>(grid), omega(omega), x0(x0), z0(z0), vv(vv)
  { }
 
  public: virtual quantity<si::velocity, real_t> u(
    const quantity<si::length, real_t> &, 
    const quantity<si::length, real_t> &, 
    const quantity<si::length, real_t> &z
  )   
  {   
    return -omega*(z-z0);
  }

  public: virtual quantity<si::velocity, real_t> w(
    const quantity<si::length, real_t> &x, 
    const quantity<si::length, real_t> &, 
    const quantity<si::length, real_t> &
  )   
  {   
    return omega*(x-x0);
  }

  public: virtual quantity<si::velocity, real_t> v(
    const quantity<si::length, real_t> &,
    const quantity<si::length, real_t> &,
    const quantity<si::length, real_t> &
  )
  {
    return vv;
  }
};
#endif
