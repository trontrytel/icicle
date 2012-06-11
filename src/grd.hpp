/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPLv3+ (see the COPYING file or http://www.gnu.org/licenses/)
 *  @brief contains definition of the @ref grd class - a base class for all grids
 */
#pragma once
#include "cmn.hpp"
#include "mtx.hpp"

/// @brief base class for grids
template<typename real_t>
class grd
{
  public: virtual const quantity<si::length, real_t> dx() const = 0;
  public: virtual const quantity<si::length, real_t> dy() const = 0;
  public: virtual const quantity<si::length, real_t> dz() const = 0;

  public: virtual const int nx() const = 0;
  public: virtual const int ny() const = 0;
  public: virtual const int nz() const = 0;

  // returns ranges to be passed as constructors to arr
  // first and last reflect scalar indices
  public: virtual mtx::idx rng_sclr(int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int halo) const = 0;
  public: virtual mtx::rng rng_vctr(int i_min, int i_max, int halo) const = 0; // TODO: make private...
  public: virtual mtx::idx rng_vctr_x(const mtx::idx &ijk, int halo) const = 0;
  public: virtual mtx::idx rng_vctr_y(const mtx::idx &ijk, int halo) const = 0;
  public: virtual mtx::idx rng_vctr_z(const mtx::idx &ijk, int halo) const = 0;

  public: virtual quantity<si::length, real_t> x(int i, int j, int k) const = 0;
  public: virtual quantity<si::length, real_t> y(int i, int j, int k) const = 0;
  public: virtual quantity<si::length, real_t> z(int i, int j, int k) const = 0;

  public: virtual quantity<si::length, real_t> u_x(int i, int j, int k) const = 0;
  public: virtual quantity<si::length, real_t> u_y(int i, int j, int k) const = 0;
  public: virtual quantity<si::length, real_t> u_z(int i, int j, int k) const = 0;
  public: virtual quantity<si::length, real_t> v_x(int i, int j, int k) const = 0;
  public: virtual quantity<si::length, real_t> v_y(int i, int j, int k) const = 0;
  public: virtual quantity<si::length, real_t> v_z(int i, int j, int k) const = 0;
  public: virtual quantity<si::length, real_t> w_x(int i, int j, int k) const = 0;
  public: virtual quantity<si::length, real_t> w_y(int i, int j, int k) const = 0;
  public: virtual quantity<si::length, real_t> w_z(int i, int j, int k) const = 0;

  // TODO: what about non-carthesian coordinates?
  public: virtual quantity<si::length, real_t> i2x(int i) const = 0;
  public: virtual quantity<si::length, real_t> j2y(int j) const = 0;
  public: virtual quantity<si::length, real_t> k2z(int k) const = 0;
};
