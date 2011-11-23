/** @file
 *  @author Sylwester Arabas (University of Warsaw) <slayoo@igf.fuw.edu.pl>
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    the Arakawa-C staggered grid
 */
#ifndef GRD_2D_XZ_ARAKAWA_C_HPP
#  define GRD_2D_XZ_ARAKAWA_C_HPP

#  include "grd_2d-xz.hpp"

template<typename real_t>
class grd_2d_xz_arakawa_c : public grd_2d_xz<real_t>
{
  public: static const int m_half = 0;
  public: static const int p_half = 1;

  private: quantity<si::length, real_t> dx_, dz_;

  public: grd_2d_xz_arakawa_c(
    quantity<si::length, real_t> dx,
    quantity<si::length, real_t> dz
  )
    : dx_(dx), dz_(dz)
  {}

  public: quantity<si::length, real_t> i2x(int i) 
  { 
    // i(-1/2, dx) = 0 m
    // i(  0 , dx) = dx/2 
    // i( 1/2, dx) = dx
    // ...
    return (real_t(i) + real_t(.5)) * dx_;
  }

  public: quantity<si::length, real_t> k2z(int k)
  { 
    // ditto
    return (real_t(k) + real_t(.5)) * dz_;
  }

  public: quantity<si::length, real_t> dx()
  {
    return dx_;
  }

  public: quantity<si::length, real_t> dz()
  {
    return dz_;
  }

  //public: Range rng_sclr(int first, int last) 
  //{ 
  //  return Range(first, last); 
  //}

  public: Range rng_vctr(int first, int last) 
  { 
    return Range(first - m_half, last + p_half); 
  }

  public: quantity<si::length, real_t> u_x(int i, int j, int k) { return real_t(i) * dx_; }
  public: quantity<si::length, real_t> u_z(int i, int j, int k) { return (k + real_t(.5)) * dz_; }
  public: quantity<si::length, real_t> w_x(int i, int j, int k) { return (i + real_t(.5)) * dx_; }
  public: quantity<si::length, real_t> w_z(int i, int j, int k) { return real_t(k) * dz_; }

};
#endif
