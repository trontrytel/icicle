/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef GRD_2D_XZ_HPP
#  define GRD_2D_XZ_HPP

#  include "grd.hpp"

template<typename real_t>
class grd_2d_xz : public grd<real_t>
{
  public: quantity<si::length, real_t> j2y(int) 
  { 
    error_macro("j2y() requested for a 2D (XZ) grid")
  }
  public: quantity<si::length, real_t> dy() 
  { 
    error_macro("dy() requested for a 2D (XZ) grid")
  }

  public: virtual quantity<si::length, real_t> u_y(int i, int j, int k) 
  {
    error_macro("u_y() requested for a 2D (XZ) grid")
  }

  public: virtual quantity<si::length, real_t> v_y(int i, int j, int k) 
  {
    error_macro("v_y() requested for a 2D (XZ) grid")
  }

  public: virtual quantity<si::length, real_t> w_y(int i, int j, int k) 
  {
    error_macro("w_y() requested for a 2D (XZ) grid")
  }

  public: virtual quantity<si::length, real_t> v_x(int i, int j, int k) 
  {
    error_macro("v_x() requested for a 2D (XZ) grid")
  }

  public: virtual quantity<si::length, real_t> v_z(int i, int j, int k) 
  {
    error_macro("v_z() requested for a 2D (XZ) grid")
  }
};
#endif
