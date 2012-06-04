/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the grd_arakawa_c_lorenz class representing the Arakawa-C/Lorenz staggered grid
 */
#ifndef GRD_ARAKAWA_C_LORENZ_HPP
#  define GRD_ARAKAWA_C_LORENZ_HPP

#  include "grd.hpp"

/// @brief the Arakawa-C/Lorenz staggered grid
template<typename real_t>
class grd_arakawa_c_lorenz : public grd<real_t> // TODO: rmerge it with grid... KISS!
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

  public: quantity<si::length, real_t> i2x(int i) const
  { 
    // i(-1/2, dx) = 0 m
    // i(  0 , dx) = dx/2 
    // i( 1/2, dx) = dx
    // ...
    return (real_t(i) + real_t(.5)) * dx_;
  }

  // ditto
  public: quantity<si::length, real_t> j2y(int j) const { return (real_t(j) + real_t(.5)) * dy_; }
  public: quantity<si::length, real_t> k2z(int k) const { return (real_t(k) + real_t(.5)) * dz_; }

  public: const quantity<si::length, real_t> dx() const { return dx_; }
  public: const quantity<si::length, real_t> dy() const { return dy_; }
  public: const quantity<si::length, real_t> dz() const { return dz_; }

  public: const int nx() const { return nx_; }
  public: const int ny() const { return ny_; }
  public: const int nz() const { return nz_; }

  private: mtx::rng rng_sclr(int first, int last, int halo) const
  { 
    return mtx::rng(first - halo, last + halo); 
  }

  public: mtx::rng rng_vctr(int first, int last, int halo) const
  { 
    assert(halo > 0);
    return mtx::rng(first - m_half - halo, last + p_half + halo); 
  }

  public: mtx::idx rng_sclr(int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int halo) const
  {
    return mtx::idx_ijk(
      rng_sclr(i_min, i_max, halo),
      rng_sclr(j_min, j_max, halo),
      rng_sclr(k_min, k_max, halo)
    );
  }

  public: mtx::idx rng_vctr_x(const mtx::idx &ijk, int halo) const
  {
    return mtx::idx_ijk(
      rng_vctr(ijk.lbound(0), ijk.ubound(0), halo),
      rng_sclr(ijk.lbound(1), ijk.ubound(1), halo),
      rng_sclr(ijk.lbound(2), ijk.ubound(2), halo)
    );
  }

  public: mtx::idx rng_vctr_y(const mtx::idx &ijk, int halo) const
  {
    return mtx::idx_ijk(
      rng_sclr(ijk.lbound(0), ijk.ubound(0), halo),
      rng_vctr(ijk.lbound(1), ijk.ubound(1), halo),
      rng_sclr(ijk.lbound(2), ijk.ubound(2), halo)
    );
  }

  public: mtx::idx rng_vctr_z(const mtx::idx &ijk, int halo) const
  {
    return mtx::idx_ijk(
      rng_sclr(ijk.lbound(0), ijk.ubound(0), halo),
      rng_sclr(ijk.lbound(1), ijk.ubound(1), halo),
      rng_vctr(ijk.lbound(2), ijk.ubound(2), halo)
    );
  }

  // coordinates at which u(x,y,z), v(x,y,z), w(x,y,z) are evaluated
  public: quantity<si::length, real_t> u_x(int i, int  , int  ) const { return real_t(i) * dx_; }
  public: quantity<si::length, real_t> u_y(int  , int j, int  ) const { return (j + real_t(.5)) * dy_; }
  public: quantity<si::length, real_t> u_z(int  , int  , int k) const { return (k + real_t(.5)) * dz_; }

  public: quantity<si::length, real_t> v_x(int i, int  , int  ) const { return (i + real_t(.5)) * dx_; }
  public: quantity<si::length, real_t> v_y(int  , int j, int  ) const { return real_t(j) * dy_; }
  public: quantity<si::length, real_t> v_z(int  , int  , int k) const { return (k + real_t(.5)) * dz_; }

  public: quantity<si::length, real_t> w_x(int i, int  , int  ) const { return (i + real_t(.5)) * dx_; }
  public: quantity<si::length, real_t> w_y(int  , int j, int  ) const { return (j + real_t(.5)) * dy_; }
  public: quantity<si::length, real_t> w_z(int  , int  , int k) const { return real_t(k) * dz_; }

  public: quantity<si::length, real_t> x(int i, int  , int  ) const { return (i + real_t(.5)) * dx_; };
  public: quantity<si::length, real_t> y(int  , int j, int  ) const { return (j + real_t(.5)) * dy_; };
  public: quantity<si::length, real_t> z(int  , int  , int k) const { return (k + real_t(.5)) * dz_;};

};
#endif
