/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    uniform velocity field
 */
#ifndef VEL_STREAMFUNC_HPP
#  define VEL_STREAMFUNC_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting
#  include "vel.hpp"

template <typename real_t>
class vel_streamfunc : public vel<real_t>
{
  public: virtual quantity<si::velocity> v(
    const quantity<si::length, real_t> &x,
    const quantity<si::length, real_t> &y,
    const quantity<si::length, real_t> &z
  ) 
  {  
    error_macro("specifying velocity field using a stream function is valid for 2D flows only!");
  }
};
#endif
