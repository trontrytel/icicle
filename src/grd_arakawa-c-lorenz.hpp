/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 *  @section DESCRIPTION
 *    the Arakawa-C/Lorenz staggered grid
 */
#ifndef GRD_ARAKAWA_C_LORENZ_HPP
#  define GRD_ARAKAWA_C_LORENZ_HPP

#  include "grd.hpp"

template<typename real_t>
class grd_arakawa_c_lorenz : public grd<real_t>
{
  public: static const int m_half = 0;
  public: static const int p_half = 1;

  private: quantity<si::length, real_t> dx_, dy_, dz_;
  private: int nx_, ny_, nz_;

  public: grd_arakawa_c_lorenz(
    quantity<si::length, real_t> dx,
    quantity<si::length, real_t> dy,
    quantity<si::length, real_t> dz,
    int nx, int ny, int nz
  )
    : dx_(dx), dy_(dy), dz_(dz), nx_(nx), ny_(ny), nz_(nz)
  {}

  public: quantity<si::length, real_t> i2x(int i) 
  { 
    // i(-1/2, dx) = 0 m
    // i(  0 , dx) = dx/2 
    // i( 1/2, dx) = dx
    // ...
    return (real_t(i) + real_t(.5)) * dx_;
  }

  // ditto
  public: quantity<si::length, real_t> j2y(int j) { return (real_t(j) + real_t(.5)) * dy_; }
  public: quantity<si::length, real_t> k2z(int k) { return (real_t(k) + real_t(.5)) * dz_; }

  public: quantity<si::length, real_t> dx() { return dx_; }
  public: quantity<si::length, real_t> dy() { return dy_; }
  public: quantity<si::length, real_t> dz() { return dz_; }

  public: int nx() { return nx_; }
  public: int ny() { return ny_; }
  public: int nz() { return nz_; }

  public: mtx::rng rng_sclr(int first, int last, int halo) 
  { 
    return mtx::rng(first - halo, last + halo); 
  }
  public: mtx::rng rng_vctr(int first, int last, int halo) 
  { 
    assert(halo > 0);
    return mtx::rng(first - m_half - halo, last + p_half + halo); 
  }

  // coordinates at which u(x,y,z), v(x,y,z), w(x,y,z) are evaluated
  public: quantity<si::length, real_t> u_x(int i, int  , int  ) { return real_t(i) * dx_; }
  public: quantity<si::length, real_t> u_y(int  , int j, int  ) { return (j + real_t(.5)) * dy_; }
  public: quantity<si::length, real_t> u_z(int  , int  , int k) { return (k + real_t(.5)) * dz_; }

  public: quantity<si::length, real_t> v_x(int i, int  , int  ) { return (i + real_t(.5)) * dx_; }
  public: quantity<si::length, real_t> v_y(int  , int j, int  ) { return real_t(j) * dy_; }
  public: quantity<si::length, real_t> v_z(int  , int  , int k) { return (k + real_t(.5)) * dz_; }

  public: quantity<si::length, real_t> w_x(int i, int  , int  ) { return (i + real_t(.5)) * dx_; }
  public: quantity<si::length, real_t> w_y(int  , int j, int  ) { return (j + real_t(.5)) * dy_; }
  public: quantity<si::length, real_t> w_z(int  , int  , int k) { return real_t(k) * dz_; }

  public: quantity<si::length, real_t> x(int i, int  , int  ) { return (i + real_t(.5)) * dx_; };
  public: quantity<si::length, real_t> y(int  , int j, int  ) { return (j + real_t(.5)) * dy_; };
  public: quantity<si::length, real_t> z(int  , int  , int k) { return (k + real_t(.5)) * dz_;};

};
#endif
