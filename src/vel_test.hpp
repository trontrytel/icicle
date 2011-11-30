/** @file
 *  @author Anna Jaruga <ajaruga@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef VEL_TEST_HPP
#  define VEL_TEST_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting
#  include "vel.hpp"

template <typename real_t>
class vel_test : public vel<real_t>
{
  private: quantity<si::frequency, real_t> omega;
  private: quantity<si::length, real_t> x0,z0;

  public: vel_test(
    quantity<si::frequency, real_t> omega,
    quantity<si::length, real_t> x0,
    quantity<si::length, real_t> z0
  )
    : omega(omega), x0(x0), z0(z0)
  {}
 
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
    return real_t(0.) * si::metres / si::seconds;
  }
};
#endif
