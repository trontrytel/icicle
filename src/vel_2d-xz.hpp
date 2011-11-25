/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    2D velocity field
 */
#ifndef VEL_2D_XZ_HPP
#  define VEL_2D_XZ_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting
#  include "vel.hpp"

template <typename real_t>
class vel_2d_xz : public vel<real_t>
{
  public: virtual quantity<si::velocity> v(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) 
  {  
    error_macro("Y component requested with a 2D (XZ) velcotiy field!");
  }
};
#endif
