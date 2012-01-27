/** @file
 *  @author Sylwester Arabas <slayoo@igf.fuw.edu.pl>
 *  @copyright University of Warsaw
 *  @date November 2011
 *  @section LICENSE
 *    GPL v3 (see the COPYING file or http://www.gnu.org/licenses/)
 */
#ifndef GRD_HPP
#  define GRD_HPP

#  include "cmn.hpp" // root class, error reporting
#  include "vel.hpp"
#  include "ini.hpp"
#  include "mtx.hpp"

template<typename real_t>
class grd : root
{
  public: virtual quantity<si::length, real_t> dx() = 0;
  public: virtual quantity<si::length, real_t> dy() = 0;
  public: virtual quantity<si::length, real_t> dz() = 0;

  public: virtual int nx() = 0;
  public: virtual int ny() = 0;
  public: virtual int nz() = 0;

  // returns ranges to be passed as constructors to arr
  // first and last reflect scalar indices
  public: virtual mtx::idx rng_sclr(int i_min, int i_max, int j_min, int j_max, int k_min, int k_max, int halo) = 0;
  public: virtual mtx::rng rng_vctr(int i_min, int i_max, int halo) = 0; // TODO: make private...
  public: virtual mtx::idx rng_vctr_x(const mtx::idx &ijk, int halo) = 0;
  public: virtual mtx::idx rng_vctr_y(const mtx::idx &ijk, int halo) = 0;
  public: virtual mtx::idx rng_vctr_z(const mtx::idx &ijk, int halo) = 0;

  public: virtual quantity<si::length, real_t> x(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> y(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> z(int i, int j, int k) = 0;

  public: virtual quantity<si::length, real_t> u_x(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> u_y(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> u_z(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> v_x(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> v_y(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> v_z(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> w_x(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> w_y(int i, int j, int k) = 0;
  public: virtual quantity<si::length, real_t> w_z(int i, int j, int k) = 0;

  // ... Nothing is real and nothing to get hung about.
  // courant fields forever ...
  public: virtual void populate_courant_fields(
    mtx::arr<real_t> *Cx, 
    mtx::arr<real_t> *Cy, 
    mtx::arr<real_t> *Cz, 
    vel<real_t> *velocity,
    quantity<si::time, real_t> dt
  ) 
  {
    quantity<si::length, real_t> 
      mx = real_t(nx()) * dx(),
      my = real_t(ny()) * dy(),
      mz = real_t(nz()) * dz();
    for (int i = Cx->lbound(mtx::i); i <= Cx->ubound(mtx::i); ++i)
      for (int j = Cx->lbound(mtx::j); j <= Cx->ubound(mtx::j); ++j)
        for (int k = Cx->lbound(mtx::k); k <= Cx->ubound(mtx::k); ++k)
          (*Cx)(i, j, k) = dt / dx() *
            velocity->u(fmod(u_x(i, j, k), mx), fmod(u_y(i, j, k), my), fmod(u_z(i, j, k), mz));
    for (int i = Cy->lbound(mtx::i); i <= Cy->ubound(mtx::i); ++i)
      for (int j = Cy->lbound(mtx::j); j <= Cy->ubound(mtx::j); ++j)
        for (int k = Cy->lbound(mtx::k); k <= Cy->ubound(mtx::k); ++k)
          (*Cy)(i, j, k) = dt / dy() *
            velocity->v(fmod(v_x(i, j, k), mx), fmod(v_y(i, j, k), my), fmod(v_z(i, j, k), mz));
    for (int i = Cz->lbound(mtx::i); i <= Cz->ubound(mtx::i); ++i)
      for (int j = Cz->lbound(mtx::j); j <= Cz->ubound(mtx::j); ++j)
        for (int k = Cz->lbound(mtx::k); k <= Cz->ubound(mtx::k); ++k)
          (*Cz)(i, j, k) = dt / dz() *
            velocity->w(fmod(w_x(i, j, k), mx), fmod(w_y(i, j, k), my), fmod(w_z(i, j, k), mz));
  }

  public: virtual void populate_scalar_field(
    const mtx::idx &ijk,
    mtx::arr<real_t> *psi, 
    ini<real_t> *intcond
  )
  {
    for (int i = ijk.lbound(mtx::i); i <= ijk.ubound(mtx::i); ++i)
      for (int j = ijk.lbound(mtx::j); j <= ijk.ubound(mtx::j); ++j)
        for (int k = ijk.lbound(mtx::k); k <= ijk.ubound(mtx::k); ++k)
          (*psi)(i,j,k) = intcond->psi(x(i,j,k), y(i,j,k), z(i,j,k));
  }
};
#endif
