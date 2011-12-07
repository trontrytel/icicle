/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    base class for velocity fields
 */
#ifndef VEL_HPP
#  define VEL_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting

template <typename real_t>
class vel : root
{
  public: virtual quantity<si::velocity, real_t> u(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) = 0;
  public: virtual quantity<si::velocity, real_t> v(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) = 0;
  public: virtual quantity<si::velocity, real_t> w(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) = 0;
};
#endif