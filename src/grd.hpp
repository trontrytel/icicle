/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef GRD_HPP
#  define GRD_HPP

#  include "config.hpp" // USE_* defines
#  include "common.hpp" // root class, error reporting

template<typename real_t>
class grd : root
{
  public: static const int m_half = 0;
  public: static const int p_half = 1;
  public: static inline quantity<si::length, real_t> i2x(int i, quantity<si::length, real_t> dx) 
  { 
    // i(-1/2, dx) = 0 m
    // i(  0 , dx) = dx/2 
    // i( 1/2, dx) = dx
    // ...
    return (real_t(i) + real_t(.5)) * dx;
  }

  public: static inline quantity<si::length, real_t> j2y(int j, quantity<si::length, real_t> dy) 
  { 
    return i2x(j, dy); 
  }

  public: static inline quantity<si::length, real_t> k2z(int k, quantity<si::length, real_t> dz) 
  { 
    return i2x(k, dz); 
  }
};
#endif
